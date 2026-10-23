//! Voronoi communities.
//!
//! Deritei, Lazar, Papp, Jarai-Szabo, Sumi, Varga, Ravasz, and
//! Ercsey-Ravasz (2014, `deritei2014community`) partition a network into
//! the graph Voronoi cells of "generator" nodes of high local density;
//! Molnar, Marton, Horvat, and Ercsey-Ravasz (2024, `molnar2024community`)
//! extend the method to weighted and directed networks, choosing the
//! radius that maximizes modularity. This is the method of igraph's
//! `cluster_voronoi()`. In their notation:
//!
//! - Edge lengths: l_ij = f(w_ij) / C_ij, where f(w_ij) is the given
//!   length of the edge (1 when none is given) and C_ij = (z_ij + 1) /
//!   min(k_i - 1, k_j - 1) is the edge clustering coefficient of Radicchi,
//!   Castellano, Cecconi, Loreto, and Parisi (2004, `radicchi2004defining`),
//!   z_ij the number of common neighbors and k the degree in the underlying
//!   simple undirected network. Edges in triangles become shorter.
//! - Local relative density: rho_i = s_i m_i / (m_i + k_i), where s_i is
//!   the strength (weighted degree, in plus out), m_i the number of edges
//!   among i and its neighbors, and k_i the number of edges leaving that
//!   set.
//! - Generators for radius r: repeatedly take the densest node not yet
//!   excluded, make it a generator, and exclude every node within distance
//!   less than r of it, until every node is a generator or excluded.
//! - Each node joins the cell of its nearest generator.
//! - Unless r is given, it is chosen to maximize the modularity of the
//!   resulting partition.
//!
//! Conventions:
//! - Directed and undirected networks. `mode` sets the direction of
//!   distances: `Out` measures them from a generator along the arcs (the
//!   generator reaches the node), `In` against them (the node reaches the
//!   generator), and `All` ignores direction. The same direction is used
//!   for excluding nodes and for assigning them to cells. Undirected
//!   networks ignore `mode`.
//! - Lengths (if given) are one finite, non-negative value per edge; the
//!   ECC, degrees, and neighborhood counts ignore direction, parallel
//!   edges, and self-loops. When min(k_i - 1, k_j - 1) = 0 (an edge to a
//!   node of degree one) C_ij is infinite and the length is 0, the limit
//!   of the formula, so such a node is always in its neighbor's cell.
//!   Among parallel arcs the shortest is used.
//! - Weights (if used) are added as given, must be non-negative, and
//!   enter the strengths and the modularity; lengths and weights are
//!   independent, as in igraph. Self-loops are left out of strengths.
//! - Ties: generators are taken in decreasing density, ties going to the
//!   lower node index; a node at equal distance from several generators
//!   joins the one taken first (the paper chooses at random). Distances
//!   are sums of floating-point lengths, so two distances (and two
//!   modularities) within a relative 1e-9 of each other are taken as
//!   equal: a node at distance r from a generator (to that precision) is
//!   not excluded, and ties to that precision follow the rules above. A
//!   node with no edges has density 0 (0/0 is taken as 0).
//! - Choice of radius: the generator set changes only where r crosses a
//!   shortest-path distance, so the candidates are every distinct finite
//!   distance between two nodes (and infinity, which leaves one generator
//!   per set of nodes it reaches); when there are more than
//!   1,000 distinct distances, 1,000 of them evenly spaced by rank are
//!   tried. The candidate of highest modularity is used, ties going to the
//!   smallest radius. The paper (and igraph) instead search with Brent's
//!   method, a local optimizer, so where they differ this choice has the
//!   higher modularity. Finding the distances takes one shortest-path
//!   search per node (in parallel over `threads`) and memory for the
//!   distinct distance values.
//! - Modularity is [`crate::metrics::community::modularity`] (directed for
//!   directed networks, resolution 1), NaN without edges. Groups are
//!   labelled 0, 1, ... in order of their first node; generators are
//!   reported in the order taken.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use super::relabel;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::community::modularity;

/// Direction of distances in [`voronoi`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VoronoiMode {
    /// From the generator along the arcs.
    Out,
    /// From the node to the generator.
    In,
    /// Direction ignored.
    All,
}

impl VoronoiMode {
    /// From the C ABI code: 1 = out, 2 = in, 3 = all.
    pub fn from_code(code: i32) -> Result<VoronoiMode> {
        match code {
            1 => Ok(VoronoiMode::Out),
            2 => Ok(VoronoiMode::In),
            3 => Ok(VoronoiMode::All),
            _ => Err(Error::InvalidArgument(format!(
                "unknown mode code {code} (1 = out, 2 = in, 3 = all)"
            ))),
        }
    }
}

/// Result of [`voronoi`].
#[derive(Debug, Clone, PartialEq)]
pub struct Voronoi {
    /// Group (cell) of each node.
    pub membership: Vec<usize>,
    /// Generators, in the order taken (densest first).
    pub generators: Vec<usize>,
    /// Modularity of the partition.
    pub modularity: f64,
    /// Radius used.
    pub radius: f64,
}

/// Arcs in the direction distances follow, as a CSR.
struct Arcs {
    offsets: Vec<usize>,
    to: Vec<usize>,
    len: Vec<f64>,
}

#[derive(PartialEq)]
struct Item(f64, usize, usize);
impl Eq for Item {}
impl PartialOrd for Item {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Item {
    // Min-heap on (distance, rank, node).
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .0
            .total_cmp(&self.0)
            .then(other.1.cmp(&self.1))
            .then(other.2.cmp(&self.2))
    }
}

impl Arcs {
    fn row(&self, v: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let r = self.offsets[v]..self.offsets[v + 1];
        self.to[r.clone()]
            .iter()
            .copied()
            .zip(self.len[r].iter().copied())
    }

    /// Distances from `s`, visiting only nodes closer than `limit`; calls
    /// `f(node, distance)` for each settled node.
    fn search(&self, s: usize, limit: f64, dist: &mut [f64], f: &mut impl FnMut(usize, f64)) {
        let mut heap = BinaryHeap::new();
        let mut settled = Vec::new();
        dist[s] = 0.0;
        heap.push(Item(0.0, 0, s));
        while let Some(Item(d, _, v)) = heap.pop() {
            if d > dist[v] {
                continue;
            }
            settled.push(v);
            f(v, d);
            for (u, l) in self.row(v) {
                let nd = d + l;
                if nd < dist[u] && nd < limit {
                    dist[u] = nd;
                    heap.push(Item(nd, 0, u));
                }
            }
        }
        for v in settled {
            dist[v] = f64::INFINITY;
        }
        // Nodes reached but never settled cannot exist: every pushed node is
        // popped. Reset is complete.
    }
}

/// Relative tolerance under which two distances (or modularities) are
/// taken as equal.
const TOL: f64 = 1e-9;

fn same(a: f64, b: f64) -> bool {
    a == b || (a.is_finite() && b.is_finite() && (a - b).abs() <= TOL * a.abs().max(b.abs()))
}

fn evaluate(
    g: &Graph,
    arcs: &Arcs,
    rev: &Arcs,
    order: &[usize],
    radius: f64,
    weighted: bool,
) -> Result<Cells> {
    let n = g.node_count();
    let mut excluded = vec![false; n];
    let mut dist = vec![f64::INFINITY; n];
    let mut generators = Vec::new();
    let reach = if radius.is_finite() {
        radius * (1.0 + 2.0 * TOL)
    } else {
        radius
    };
    for &v in order {
        if excluded[v] {
            continue;
        }
        generators.push(v);
        excluded[v] = true;
        arcs.search(v, reach, &mut dist, &mut |u, d| {
            if d < radius && !same(d, radius) {
                excluded[u] = true;
            }
        });
    }
    // Distance to the nearest generator, then the first generator (in the
    // order taken) among those at that distance, through the shortest-path
    // arcs.
    let mut near = vec![f64::INFINITY; n];
    let mut heap = BinaryHeap::new();
    for &s in &generators {
        near[s] = 0.0;
        heap.push(Item(0.0, 0, s));
    }
    while let Some(Item(d, _, v)) = heap.pop() {
        if d > near[v] {
            continue;
        }
        for (u, l) in arcs.row(v) {
            if d + l < near[u] {
                near[u] = d + l;
                heap.push(Item(d + l, 0, u));
            }
        }
    }
    let mut rank = vec![usize::MAX; n];
    for (r, &s) in generators.iter().enumerate() {
        rank[s] = r;
    }
    let mut by_distance: Vec<usize> = (0..n).collect();
    by_distance.sort_by(|&a, &b| near[a].total_cmp(&near[b]).then(a.cmp(&b)));
    loop {
        let mut changed = false;
        for &v in &by_distance {
            for (u, l) in rev.row(v) {
                if rank[u] < rank[v] && same(near[u] + l, near[v]) {
                    rank[v] = rank[u];
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let membership = relabel(&rank);
    let memb: Vec<i64> = membership.iter().map(|&c| c as i64).collect();
    let q = modularity(g, &memb, 1.0, weighted)?;
    Ok((membership, generators, q))
}

const MAX_CANDIDATES: usize = 1000;

/// Membership, generators, and modularity for one radius.
type Cells = (Vec<usize>, Vec<usize>, f64);

/// Voronoi communities. `lengths` (one per edge, or `None` for 1) give
/// edge lengths; `radius` is used when given, otherwise chosen to maximize
/// modularity. See the module documentation.
pub fn voronoi(
    g: &Graph,
    lengths: Option<&[f64]>,
    weighted: bool,
    mode: VoronoiMode,
    radius: Option<f64>,
    threads: usize,
) -> Result<Voronoi> {
    let n = g.node_count();
    let m = g.edge_count();
    let w = crate::weights::additive_weights(g, weighted, "Voronoi communities")?;
    if let Some(l) = lengths {
        if l.len() != m {
            return Err(Error::LengthMismatch {
                what: "lengths",
                expected: m,
                found: l.len(),
            });
        }
        if let Some(e) = l.iter().position(|&x| !(x.is_finite() && x >= 0.0)) {
            return Err(Error::InvalidArgument(format!(
                "edge {e} has length {}; lengths must be finite and non-negative",
                l[e]
            )));
        }
    }
    if let Some(r) = radius {
        if r.is_nan() || r < 0.0 {
            return Err(Error::InvalidArgument(format!(
                "the radius must be non-negative, not {r}"
            )));
        }
    }
    // Underlying simple undirected network: degrees, common neighbors,
    // neighborhood edge counts.
    let simple = g.simple_undirected(false);
    let deg: Vec<usize> = (0..n).map(|v| simple.neighbors(v).len()).collect();
    let mut mark = vec![usize::MAX; n];
    // Triangles at each node, for the density.
    let mut tri = vec![0usize; n];
    for (v, t) in tri.iter_mut().enumerate() {
        for &u in simple.neighbors(v) {
            mark[u as usize] = v;
        }
        for &u in simple.neighbors(v) {
            for &x in simple.neighbors(u as usize) {
                if (x as usize) > (u as usize) && mark[x as usize] == v {
                    *t += 1;
                }
            }
        }
    }
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut strength = vec![0.0f64; n];
    for e in 0..m {
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        if u != v {
            let x = w.map_or(1.0, |w| w[e]);
            strength[u] += x;
            strength[v] += x;
        }
    }
    let density: Vec<f64> = (0..n)
        .map(|v| {
            let inside = deg[v] + tri[v];
            let boundary: usize = deg[v]
                + simple
                    .neighbors(v)
                    .iter()
                    .map(|&u| deg[u as usize])
                    .sum::<usize>()
                - 2 * inside;
            if inside + boundary == 0 {
                0.0
            } else {
                strength[v] * inside as f64 / (inside + boundary) as f64
            }
        })
        .collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| density[b].total_cmp(&density[a]).then(a.cmp(&b)));
    // Edge lengths divided by the ECC.
    let mut seen = vec![false; n];
    let mut ecc_len = |e: usize| -> f64 {
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        let base = lengths.map_or(1.0, |l| l[e]);
        let lo = deg[u].min(deg[v]);
        if lo <= 1 {
            return 0.0;
        }
        for &x in simple.neighbors(u) {
            seen[x as usize] = true;
        }
        let z = simple
            .neighbors(v)
            .iter()
            .filter(|&&x| seen[x as usize])
            .count();
        for &x in simple.neighbors(u) {
            seen[x as usize] = false;
        }
        base * (lo - 1) as f64 / (z + 1) as f64
    };
    let directed = g.is_directed();
    let mut arcs_list: Vec<(usize, usize, f64)> = Vec::with_capacity(2 * m);
    for e in 0..m {
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        if u == v {
            continue;
        }
        let l = ecc_len(e);
        match (directed, mode) {
            (true, VoronoiMode::Out) => arcs_list.push((u, v, l)),
            (true, VoronoiMode::In) => arcs_list.push((v, u, l)),
            _ => {
                arcs_list.push((u, v, l));
                arcs_list.push((v, u, l));
            }
        }
    }
    let mut offsets = vec![0usize; n + 1];
    for &(u, _, _) in &arcs_list {
        offsets[u + 1] += 1;
    }
    for i in 0..n {
        offsets[i + 1] += offsets[i];
    }
    let mut fill = offsets.clone();
    let mut to = vec![0usize; arcs_list.len()];
    let mut len = vec![0.0f64; arcs_list.len()];
    for &(u, v, l) in &arcs_list {
        to[fill[u]] = v;
        len[fill[u]] = l;
        fill[u] += 1;
    }
    let arcs = Arcs { offsets, to, len };
    let rev = {
        let mut offsets = vec![0usize; n + 1];
        for &(_, v, _) in &arcs_list {
            offsets[v + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let mut fill = offsets.clone();
        let mut to = vec![0usize; arcs_list.len()];
        let mut len = vec![0.0f64; arcs_list.len()];
        for &(u, v, l) in &arcs_list {
            to[fill[v]] = u;
            len[fill[v]] = l;
            fill[v] += 1;
        }
        Arcs { offsets, to, len }
    };
    if let Some(r) = radius {
        let (membership, generators, q) = evaluate(g, &arcs, &rev, &order, r, weighted)?;
        return Ok(Voronoi {
            membership,
            generators,
            modularity: q,
            radius: r,
        });
    }
    // Distinct finite distances between distinct nodes.
    let parts = crate::parallel::for_each_item(
        n,
        threads,
        || (vec![f64::INFINITY; n], Vec::<f64>::new()),
        |(dist, vals), s| {
            arcs.search(s, f64::INFINITY, dist, &mut |v, d| {
                if v != s {
                    vals.push(d);
                }
            });
            if vals.len() > 4 * n + 1024 {
                vals.sort_by(f64::total_cmp);
                vals.dedup();
            }
        },
    );
    let mut values: Vec<f64> = parts.into_iter().flat_map(|(_, v)| v).collect();
    values.sort_by(f64::total_cmp);
    values.dedup_by(|b, a| same(*a, *b));
    if values.len() > MAX_CANDIDATES {
        let k = values.len();
        values = (0..MAX_CANDIDATES)
            .map(|i| values[i * (k - 1) / (MAX_CANDIDATES - 1)])
            .collect();
    }
    values.push(f64::INFINITY);
    let results = crate::parallel::for_each_item(
        values.len(),
        threads,
        Vec::new,
        |out: &mut Vec<(usize, Result<Cells>)>, i| {
            out.push((i, evaluate(g, &arcs, &rev, &order, values[i], weighted)));
        },
    );
    let mut all: Vec<(usize, Result<Cells>)> = results.into_iter().flatten().collect();
    all.sort_by_key(|x| x.0);
    let mut best: Option<(f64, Vec<usize>, Vec<usize>, f64)> = None;
    for (i, r) in all {
        let (membership, generators, q) = r?;
        let better = match &best {
            None => true,
            // NaN (no edges) never beats the first candidate.
            Some(b) => q > b.0 && !same(q, b.0),
        };
        if better {
            best = Some((q, membership, generators, values[i]));
        }
    }
    let (q, membership, generators, r) = best.expect("at least one candidate radius");
    Ok(Voronoi {
        membership,
        generators,
        modularity: q,
        radius: r,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_cliques() -> Graph {
        // Two 4-cliques joined by one edge (3-4).
        let mut s = Vec::new();
        let mut t = Vec::new();
        for base in [0i64, 4] {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    s.push(base + i);
                    t.push(base + j);
                }
            }
        }
        s.push(3);
        t.push(4);
        Graph::from_edges(8, &s, &t, None, false).unwrap()
    }

    #[test]
    fn two_cliques_split() {
        let g = two_cliques();
        let r = voronoi(&g, None, false, VoronoiMode::All, None, 1).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 0, 1, 1, 1, 1]);
        assert_eq!(r.generators.len(), 2);
        // Threads do not change the result.
        assert_eq!(
            r,
            voronoi(&g, None, false, VoronoiMode::All, None, 2).unwrap()
        );
        // Radius 0 makes every node a generator.
        let r0 = voronoi(&g, None, false, VoronoiMode::All, Some(0.0), 1).unwrap();
        assert_eq!(r0.membership, (0..8).collect::<Vec<_>>());
        let rinf = voronoi(&g, None, false, VoronoiMode::All, Some(f64::INFINITY), 1).unwrap();
        assert_eq!(rinf.membership, vec![0; 8]);
    }

    #[test]
    fn edge_cases() {
        let g = Graph::from_edges(3, &[], &[], None, false).unwrap();
        let r = voronoi(&g, None, false, VoronoiMode::All, None, 1).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.modularity.is_nan());
        let g = two_cliques();
        assert!(voronoi(&g, Some(&[1.0]), false, VoronoiMode::All, None, 1).is_err());
        assert!(voronoi(&g, None, false, VoronoiMode::All, Some(-1.0), 1).is_err());
    }
}
