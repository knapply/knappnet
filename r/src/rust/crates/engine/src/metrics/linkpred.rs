//! Link prediction scores for pairs of nodes, from the neighbourhoods of
//! the two nodes.
//!
//! For nodes u and v with neighbour sets N(u) and N(v), common neighbours
//! CN(u, v) = (N(u) ∩ N(v)) \ {u, v} and simple degrees k:
//!
//! | Method | Score | Source |
//! |---|---|---|
//! | common neighbours | \|CN\| | `liben2007link` |
//! | Jaccard | \|CN\| / \|N(u) ∪ N(v)\| (0 for an empty union) | `liben2007link` |
//! | Adamic-Adar | Σ_{w ∈ CN} 1 / ln k_w | `adamic2003friends` |
//! | resource allocation | Σ_{w ∈ CN} 1 / k_w | `zhou2009predicting` |
//! | preferential attachment | k_u k_v | `liben2007link`, `barabasi1999emergence` |
//! | common-neighbour centrality | α \|CN\| + (1 − α) n / d(u, v) | `ahmad2020missing` |
//! | CN Soundarajan-Hopcroft | \|CN\| + \|{w ∈ CN : c_w = c_u}\| if c_u = c_v, else \|CN\| | `soundarajan2012using` |
//! | RA Soundarajan-Hopcroft | Σ_{w ∈ CN, c_w = c_u} 1 / k_w if c_u = c_v, else 0 | `soundarajan2012using` |
//! | within-inter-cluster | \|W\| / (\|CN \ W\| + δ), W = {w ∈ CN : c_w = c_u}, if c_u = c_v, else 0 | `valverde2012link` |
//!
//! c is a community label per node. The formulas, and every convention
//! below, are those of NetworkX 3.7's `algorithms.link_prediction`
//! (`hagberg2008exploring`), which the scores reproduce exactly.
//!
//! Conventions:
//! - Self-loops: a looped node is in its own neighbour set (it counts in
//!   unions) and its loop adds 2 to its degree, as in a NetworkX `Graph`;
//!   u and v are never their own common neighbours.
//! - Parallel edges are merged: scores depend only on which pairs are
//!   adjacent (NetworkX does not accept multigraphs here).
//! - Direction is ignored (the underlying undirected graph; NetworkX
//!   accepts only undirected graphs).
//! - Weights are ignored.
//! - Distances d(u, v) for common-neighbour centrality are hop counts; an
//!   unreachable pair has n / d = 0. The pair (u, u) is an error there, as
//!   in NetworkX. For other methods (u, u) is scored by the same formulas.
//! - Adamic-Adar of (u, u) with a common neighbour of degree 1 divides by
//!   ln 1 = 0 and gives +inf (NetworkX raises `ZeroDivisionError`).
//! - Sums run over common neighbours in increasing index order, so a score
//!   is identical for every thread count; they can differ from NetworkX's
//!   (which sums in set order) in the last bits.
//!
//! All non-edges are the unordered pairs u < v of distinct non-adjacent
//! nodes, in lexicographic order (NetworkX's `non_edges` yields the same
//! pairs in set order). For given pairs, common neighbours are population
//! counts of ANDed bit rows ([`crate::bitset::NeighborSets`]) on dense
//! networks and merges of sorted rows otherwise; for all non-edges, each
//! node's two-step neighbourhood is tallied once, so a row of n scores
//! costs the sum of its neighbours' degrees plus n.

use crate::bitset::NeighborSets;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;

/// A link prediction score.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinkMethod {
    /// Number of common neighbours.
    CommonNeighbors,
    /// Jaccard coefficient.
    Jaccard,
    /// Adamic-Adar index.
    AdamicAdar,
    /// Resource allocation index.
    ResourceAllocation,
    /// Preferential attachment (product of degrees).
    PreferentialAttachment,
    /// Common-neighbour centrality with weight `alpha`.
    CommonNeighborCentrality,
    /// Common neighbours with the Soundarajan-Hopcroft community bonus.
    CnSoundarajanHopcroft,
    /// Resource allocation within the pair's community.
    RaSoundarajanHopcroft,
    /// Within-inter-cluster ratio with smoothing `delta`.
    WithinInterCluster,
}

impl LinkMethod {
    /// Decodes the C interface codes: 0 = common neighbours, 1 = Jaccard,
    /// 2 = Adamic-Adar, 3 = resource allocation, 4 = preferential
    /// attachment, 5 = common-neighbour centrality, 6 = CN
    /// Soundarajan-Hopcroft, 7 = RA Soundarajan-Hopcroft, 8 =
    /// within-inter-cluster.
    pub fn from_code(code: i32) -> Result<LinkMethod> {
        use LinkMethod::*;
        Ok(match code {
            0 => CommonNeighbors,
            1 => Jaccard,
            2 => AdamicAdar,
            3 => ResourceAllocation,
            4 => PreferentialAttachment,
            5 => CommonNeighborCentrality,
            6 => CnSoundarajanHopcroft,
            7 => RaSoundarajanHopcroft,
            8 => WithinInterCluster,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown link prediction method code {code}"
                )))
            }
        })
    }

    /// Whether the method reads community labels.
    pub fn needs_communities(self) -> bool {
        matches!(
            self,
            LinkMethod::CnSoundarajanHopcroft
                | LinkMethod::RaSoundarajanHopcroft
                | LinkMethod::WithinInterCluster
        )
    }
}

/// Options for [`link_scores`] and [`link_scores_non_edges`].
#[derive(Debug, Clone, Copy)]
pub struct LinkOptions {
    /// The score.
    pub method: LinkMethod,
    /// Weight of the common-neighbour count in common-neighbour
    /// centrality (NetworkX's default is 0.8).
    pub alpha: f64,
    /// Smoothing term of the within-inter-cluster ratio, > 0 (NetworkX's
    /// default is 0.001).
    pub delta: f64,
    /// Threads (0 = all cores).
    pub threads: usize,
}

impl Default for LinkOptions {
    fn default() -> Self {
        LinkOptions {
            method: LinkMethod::Jaccard,
            alpha: 0.8,
            delta: 0.001,
            threads: 0,
        }
    }
}

struct Scorer<'a> {
    sets: NeighborSets,
    communities: Option<&'a [i64]>,
    opts: LinkOptions,
    n: usize,
}

impl<'a> Scorer<'a> {
    fn new(g: &Graph, communities: Option<&'a [i64]>, opts: LinkOptions) -> Result<Scorer<'a>> {
        let n = g.node_count();
        if opts.method.needs_communities() {
            match communities {
                None => {
                    return Err(Error::InvalidArgument(
                        "this link prediction method needs a community label for every node".into(),
                    ))
                }
                Some(c) if c.len() != n => {
                    return Err(Error::LengthMismatch {
                        what: "communities",
                        expected: n,
                        found: c.len(),
                    })
                }
                _ => {}
            }
        }
        if opts.method == LinkMethod::WithinInterCluster
            && (opts.delta.is_nan() || opts.delta <= 0.0)
        {
            return Err(Error::InvalidArgument(
                "delta must be greater than zero".into(),
            ));
        }
        if opts.method == LinkMethod::CommonNeighborCentrality && !opts.alpha.is_finite() {
            return Err(Error::InvalidArgument("alpha must be finite".into()));
        }
        Ok(Scorer {
            sets: NeighborSets::new(g, true),
            communities,
            opts,
            n,
        })
    }

    /// The score of (u, v); `dist` is the hop distance (common-neighbour
    /// centrality only; `u32::MAX` when unreachable).
    fn score(&self, u: usize, v: usize, dist: u32) -> f64 {
        let s = &self.sets;
        let cn_count = || {
            let mut c = s.common_count(u, v);
            if s.adjacent(u, u) && s.adjacent(v, u) {
                c -= 1;
            }
            if u != v && s.adjacent(u, v) && s.adjacent(v, v) {
                c -= 1;
            }
            c
        };
        let each_cn = |f: &mut dyn FnMut(usize)| {
            s.for_each_common(u, v, |w| {
                if w != u && w != v {
                    f(w)
                }
            })
        };
        let comm = |w: usize| self.communities.map_or(0, |c| c[w]);
        match self.opts.method {
            LinkMethod::CommonNeighbors => cn_count() as f64,
            LinkMethod::Jaccard => {
                let both = s.common_count(u, v);
                let union = s.len(u) + s.len(v) - both;
                if union == 0 {
                    0.0
                } else {
                    cn_count() as f64 / union as f64
                }
            }
            LinkMethod::AdamicAdar => {
                let mut t = 0.0;
                each_cn(&mut |w| t += 1.0 / (s.simple_degree(w) as f64).ln());
                t
            }
            LinkMethod::ResourceAllocation => {
                let mut t = 0.0;
                each_cn(&mut |w| t += 1.0 / s.simple_degree(w) as f64);
                t
            }
            LinkMethod::PreferentialAttachment => {
                (s.simple_degree(u) as f64) * (s.simple_degree(v) as f64)
            }
            LinkMethod::CommonNeighborCentrality => {
                let cn = cn_count() as f64;
                let alpha = self.opts.alpha;
                if alpha == 1.0 {
                    return cn;
                }
                // Evaluated as NetworkX writes it, so the result is
                // bitwise the same: alpha * cn + ((1 - alpha) * n) / d.
                let d = if dist == u32::MAX {
                    f64::INFINITY
                } else {
                    dist as f64
                };
                alpha * cn + (1.0 - alpha) * self.n as f64 / d
            }
            LinkMethod::CnSoundarajanHopcroft => {
                let (cu, cv) = (comm(u), comm(v));
                let mut total = 0usize;
                let mut same = 0usize;
                each_cn(&mut |w| {
                    total += 1;
                    if comm(w) == cu {
                        same += 1;
                    }
                });
                if cu == cv {
                    (total + same) as f64
                } else {
                    total as f64
                }
            }
            LinkMethod::RaSoundarajanHopcroft => {
                let (cu, cv) = (comm(u), comm(v));
                if cu != cv {
                    return 0.0;
                }
                let mut t = 0.0;
                each_cn(&mut |w| {
                    if comm(w) == cu {
                        t += 1.0 / s.simple_degree(w) as f64
                    }
                });
                t
            }
            LinkMethod::WithinInterCluster => {
                let (cu, cv) = (comm(u), comm(v));
                if cu != cv {
                    return 0.0;
                }
                let (mut within, mut inter) = (0usize, 0usize);
                each_cn(&mut |w| {
                    if comm(w) == cu {
                        within += 1
                    } else {
                        inter += 1
                    }
                });
                within as f64 / (inter as f64 + self.opts.delta)
            }
        }
    }

    fn needs_distances(&self) -> bool {
        self.opts.method == LinkMethod::CommonNeighborCentrality && self.opts.alpha != 1.0
    }
}

/// Hop distances from `s` in the undirected graph of `sets`.
fn bfs(sets: &NeighborSets, s: usize, dist: &mut [u32], queue: &mut Vec<u32>) {
    dist.iter_mut().for_each(|d| *d = u32::MAX);
    queue.clear();
    dist[s] = 0;
    queue.push(s as u32);
    let mut head = 0;
    while head < queue.len() {
        let v = queue[head] as usize;
        head += 1;
        for &w in sets.row(v) {
            if dist[w as usize] == u32::MAX {
                dist[w as usize] = dist[v] + 1;
                queue.push(w);
            }
        }
    }
}

/// Scores of the pairs `(us[i], vs[i])`. `communities` holds one label per
/// node for the methods that use them (ignored otherwise).
pub fn link_scores(
    g: &Graph,
    us: &[u32],
    vs: &[u32],
    communities: Option<&[i64]>,
    opts: LinkOptions,
) -> Result<Vec<f64>> {
    if us.len() != vs.len() {
        return Err(Error::LengthMismatch {
            what: "second nodes of the pairs",
            expected: us.len(),
            found: vs.len(),
        });
    }
    let n = g.node_count();
    for (&u, &v) in us.iter().zip(vs) {
        if u as usize >= n || v as usize >= n {
            return Err(Error::InvalidArgument(format!(
                "pair ({u}, {v}) names a node outside 0..{n}"
            )));
        }
        if u == v && opts.method == LinkMethod::CommonNeighborCentrality {
            return Err(Error::InvalidArgument(
                "common-neighbour centrality is not defined for a node paired with itself".into(),
            ));
        }
    }
    let sc = Scorer::new(g, communities, opts)?;
    let mut out = vec![0.0; us.len()];
    if sc.needs_distances() {
        // Group the pairs by first node: one search per distinct node.
        let mut start = vec![0usize; n + 1];
        for &u in us {
            start[u as usize + 1] += 1;
        }
        for i in 0..n {
            start[i + 1] += start[i];
        }
        let mut order = vec![0usize; us.len()];
        let mut cursor = start.clone();
        for (i, &u) in us.iter().enumerate() {
            order[cursor[u as usize]] = i;
            cursor[u as usize] += 1;
        }
        let sources: Vec<usize> = (0..n).filter(|&u| start[u + 1] > start[u]).collect();
        let work = parallel::search_work(sources.len(), n, sc.sets.entries());
        let threads = parallel::threads_for(opts.threads, sources.len(), work);
        let rows = parallel::map_items(
            sources.len(),
            threads,
            || (vec![0u32; n], Vec::new()),
            |(dist, queue), k| {
                let u = sources[k];
                bfs(&sc.sets, u, dist, queue);
                order[start[u]..start[u + 1]]
                    .iter()
                    .map(|&i| sc.score(u, vs[i] as usize, dist[vs[i] as usize]))
                    .collect::<Vec<f64>>()
            },
        );
        for (k, row) in rows.into_iter().enumerate() {
            let u = sources[k];
            for (&i, x) in order[start[u]..start[u + 1]].iter().zip(row) {
                out[i] = x;
            }
        }
        return Ok(out);
    }
    let avg = 1 + sc.sets.entries() / n.max(1);
    let work = (us.len() as u64).saturating_mul(avg as u64);
    let threads = parallel::threads_for(opts.threads, us.len(), work);
    parallel::fill_chunks(
        &mut out,
        1024,
        threads,
        || (),
        |_, c, part| {
            for (k, x) in part.iter_mut().enumerate() {
                let i = c * 1024 + k;
                *x = sc.score(us[i] as usize, vs[i] as usize, 0);
            }
        },
    );
    Ok(out)
}

/// Number of non-edges: unordered pairs of distinct nodes that are not
/// adjacent in either direction.
pub fn non_edge_count(g: &Graph) -> usize {
    let sets = NeighborSets::new(g, false);
    non_edge_count_of(&sets)
}

fn non_edge_count_of(sets: &NeighborSets) -> usize {
    let n = sets.node_count();
    let loops = (0..n).filter(|&v| sets.has_loop(v)).count();
    let edges = (sets.entries() - loops) / 2;
    n * n.saturating_sub(1) / 2 - edges
}

/// Scores of every non-edge u < v, in lexicographic order, written to
/// `us`, `vs`, and `scores`, which must each hold [`non_edge_count`]
/// values.
pub fn link_scores_non_edges(
    g: &Graph,
    communities: Option<&[i64]>,
    opts: LinkOptions,
    us: &mut [u32],
    vs: &mut [u32],
    scores: &mut [f64],
) -> Result<()> {
    let sc = Scorer::new(g, communities, opts)?;
    let n = sc.n;
    let total = non_edge_count_of(&sc.sets);
    for (what, len) in [("us", us.len()), ("vs", vs.len()), ("scores", scores.len())] {
        if len != total {
            return Err(Error::LengthMismatch {
                what,
                expected: total,
                found: len,
            });
        }
    }
    // Row u holds the non-adjacent v > u; split the outputs by row.
    let mut slots = Vec::with_capacity(n);
    let (mut ur, mut vr, mut sr) = (us, vs, scores);
    for u in 0..n {
        let row = sc.sets.row(u);
        let later_adjacent = row.len() - row.partition_point(|&w| w as usize <= u);
        let len = n - 1 - u - later_adjacent;
        let (a, rest_u) = std::mem::take(&mut ur).split_at_mut(len);
        let (b, rest_v) = std::mem::take(&mut vr).split_at_mut(len);
        let (c, rest_s) = std::mem::take(&mut sr).split_at_mut(len);
        ur = rest_u;
        vr = rest_v;
        sr = rest_s;
        slots.push(std::sync::Mutex::new(Some((a, b, c))));
    }
    let work = parallel::search_work(n, n, sc.sets.entries())
        .saturating_add((n as u64).saturating_mul(n as u64));
    let threads = parallel::threads_for(opts.threads, n, work);
    let needs_dist = sc.needs_distances();
    parallel::for_each_claimed(
        n,
        threads,
        || {
            (
                vec![u32::MAX; if needs_dist { n } else { 0 }],
                Vec::new(),
                Accumulator::new(n),
            )
        },
        |(dist, queue, acc), u| {
            let taken = slots[u].lock().unwrap_or_else(|e| e.into_inner()).take();
            let (a, b, c) = taken.expect("each row is filled once");
            if needs_dist {
                bfs(&sc.sets, u, dist, queue);
            }
            acc.gather(&sc, u);
            let row = sc.sets.row(u);
            let mut r = row.partition_point(|&w| w as usize <= u);
            let mut k = 0;
            for v in u + 1..n {
                if r < row.len() && row[r] as usize == v {
                    r += 1;
                    continue;
                }
                a[k] = u as u32;
                b[k] = v as u32;
                c[k] = acc.score(&sc, u, v, dist.get(v).copied().unwrap_or(0));
                k += 1;
            }
            acc.clear();
        },
    );
    Ok(())
}

/// Common-neighbour tallies of one node u with every node reached in two
/// steps: for x, the number of common neighbours w ≠ u, how many of them
/// share u's community, and the method's sum over them (1 / k_w, 1 /
/// ln k_w, or 1 / k_w over the same community), each added in increasing
/// order of w as [`Scorer::score`] adds them, so the scores are bitwise
/// the same. One pass over the two-step neighbourhood replaces an
/// intersection per pair: O(Σ_{w ∈ N(u)} k_w + n) per node.
struct Accumulator {
    count: Vec<u32>,
    same: Vec<u32>,
    sum: Vec<f64>,
    touched: Vec<u32>,
}

impl Accumulator {
    fn new(n: usize) -> Accumulator {
        Accumulator {
            count: vec![0; n],
            same: vec![0; n],
            sum: vec![0.0; n],
            touched: Vec::new(),
        }
    }

    fn gather(&mut self, sc: &Scorer, u: usize) {
        let s = &sc.sets;
        let comm = |w: usize| sc.communities.map_or(0, |c| c[w]);
        let cu = comm(u);
        for &w in s.row(u) {
            let w = w as usize;
            if w == u {
                continue;
            }
            let same_w = comm(w) == cu;
            let term = match sc.opts.method {
                LinkMethod::AdamicAdar => 1.0 / (s.simple_degree(w) as f64).ln(),
                LinkMethod::ResourceAllocation => 1.0 / s.simple_degree(w) as f64,
                LinkMethod::RaSoundarajanHopcroft if same_w => 1.0 / s.simple_degree(w) as f64,
                _ => 0.0,
            };
            for &x in s.row(w) {
                let x = x as usize;
                if self.count[x] == 0 {
                    self.touched.push(x as u32);
                }
                self.count[x] += 1;
                if same_w {
                    self.same[x] += 1;
                }
                if term != 0.0 {
                    self.sum[x] += term;
                }
            }
        }
    }

    fn clear(&mut self) {
        for &x in &self.touched {
            let x = x as usize;
            self.count[x] = 0;
            self.same[x] = 0;
            self.sum[x] = 0.0;
        }
        self.touched.clear();
    }

    /// The score of the non-edge (u, v) from the tallies.
    fn score(&self, sc: &Scorer, u: usize, v: usize, dist: u32) -> f64 {
        let s = &sc.sets;
        let comm = |w: usize| sc.communities.map_or(0, |c| c[w]);
        let cn = self.count[v] as usize;
        let same = self.same[v] as usize;
        match sc.opts.method {
            LinkMethod::CommonNeighbors => cn as f64,
            LinkMethod::Jaccard => {
                let union = s.len(u) + s.len(v) - cn;
                if union == 0 {
                    0.0
                } else {
                    cn as f64 / union as f64
                }
            }
            LinkMethod::AdamicAdar | LinkMethod::ResourceAllocation => self.sum[v],
            LinkMethod::PreferentialAttachment => {
                (s.simple_degree(u) as f64) * (s.simple_degree(v) as f64)
            }
            LinkMethod::CommonNeighborCentrality => {
                let alpha = sc.opts.alpha;
                if alpha == 1.0 {
                    return cn as f64;
                }
                let d = if dist == u32::MAX {
                    f64::INFINITY
                } else {
                    dist as f64
                };
                alpha * cn as f64 + (1.0 - alpha) * sc.n as f64 / d
            }
            LinkMethod::CnSoundarajanHopcroft => {
                if comm(u) == comm(v) {
                    (cn + same) as f64
                } else {
                    cn as f64
                }
            }
            LinkMethod::RaSoundarajanHopcroft => {
                if comm(u) == comm(v) {
                    self.sum[v]
                } else {
                    0.0
                }
            }
            LinkMethod::WithinInterCluster => {
                if comm(u) != comm(v) {
                    return 0.0;
                }
                same as f64 / ((cn - same) as f64 + sc.opts.delta)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(method: LinkMethod) -> LinkOptions {
        LinkOptions {
            method,
            threads: 1,
            ..Default::default()
        }
    }

    #[test]
    fn hand_checked_scores_on_a_kite() {
        // 0-1, 0-2, 1-2, 1-3, 2-3, 3-4; pair (0, 3): CN = {1, 2}.
        let g =
            Graph::from_edges(5, &[0, 0, 1, 1, 2, 3], &[1, 2, 2, 3, 3, 4], None, false).unwrap();
        let s = |m| link_scores(&g, &[0], &[3], Some(&[0, 0, 0, 0, 1]), opts(m)).unwrap()[0];
        assert_eq!(s(LinkMethod::CommonNeighbors), 2.0);
        // N(0) ∪ N(3) = {1, 2, 4}.
        assert_eq!(s(LinkMethod::Jaccard), 2.0 / 3.0);
        assert_eq!(s(LinkMethod::ResourceAllocation), 2.0 / 3.0);
        assert!((s(LinkMethod::AdamicAdar) - 2.0 / 3f64.ln()).abs() < 1e-15);
        assert_eq!(s(LinkMethod::PreferentialAttachment), 6.0);
        // 0.8 * 2 + 0.2 * 5 / 2.
        assert!((s(LinkMethod::CommonNeighborCentrality) - 2.1).abs() < 1e-15);
        assert_eq!(s(LinkMethod::CnSoundarajanHopcroft), 4.0);
        assert_eq!(s(LinkMethod::RaSoundarajanHopcroft), 2.0 / 3.0);
        assert!((s(LinkMethod::WithinInterCluster) - 2000.0).abs() < 1e-9);
    }

    #[test]
    fn non_edges_are_listed_in_order() {
        let g = Graph::from_edges(4, &[0, 1, 2], &[1, 2, 2], None, false).unwrap();
        let k = non_edge_count(&g);
        assert_eq!(k, 4);
        let (mut a, mut b, mut c) = (vec![0; k], vec![0; k], vec![0.0; k]);
        link_scores_non_edges(&g, None, opts(LinkMethod::Jaccard), &mut a, &mut b, &mut c).unwrap();
        assert_eq!(a, [0, 0, 1, 2]);
        assert_eq!(b, [2, 3, 3, 3]);
        // N(0) = {1}, N(2) = {1, 2}: union {1, 2}, CN {1}.
        assert_eq!(c[0], 0.5);
        assert_eq!(c[1], 0.0);
    }
}
