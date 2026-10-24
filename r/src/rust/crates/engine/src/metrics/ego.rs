//! Ego networks: neighborhoods and egocentric measures.
//!
//! An ego network is a focal node (ego), the nodes it is tied to (its
//! alters), and the ties among them (Wasserman and Faust 1994,
//! `wasserman1994social`, section 2.3; Perry, Pescosolido, and Borgatti
//! 2018, `perry2018egocentric`). With `order` k the alters are every node
//! within k steps of ego, which generalizes the first-order neighborhood.
//!
//! Measures computed for every ego at once, in parallel:
//! - size: the number of alters;
//! - ties: the number of ties among the alters (ego's own ties excluded);
//! - pairs: the number of pairs of alters, k(k - 1)/2, or k(k - 1) ordered
//!   pairs in a directed network, for k alters;
//! - density: ties / pairs (ego network density as in Marsden 2002,
//!   `marsden2002egocentric`, and UCINET's ego network measures);
//! - components: the number of weak components among the alters once ego
//!   is removed (isolated alters are components of size 1);
//! - ego betweenness: the betweenness of ego within its own ego network
//!   (Everett and Borgatti 2005, `everett2005ego`; Marsden 2002), computed
//!   with Brandes' accumulation (Brandes 2001, `brandes2001faster`) on the
//!   ego network alone. For a first-order ego network in an undirected
//!   network, where all alters are adjacent to ego, it is the sum over
//!   non-adjacent pairs of alters of 1 / (number of their common
//!   neighbours in the ego network, ego included), Everett and Borgatti's
//!   (2005) formula; that form is used when it is cheaper (dense ego
//!   networks). Other dense ego networks (directed, or of order above 1)
//!   use the distances and shortest-path counts to and from ego: ego is on
//!   a shortest path from i to j exactly when d(i, ego) + d(ego, j) =
//!   d(i, j), with sigma(i, ego) sigma(ego, j) of the sigma_ij paths
//!   through it (Brandes 2001, section 3), so one search from each alter,
//!   stopped at depth d(i, ego) + max_j d(ego, j) and run on bit rows,
//!   replaces the full accumulation. All forms give the same value up to
//!   rounding.
//!
//! Burt's structural-hole measures, which are also egocentric, live in
//! [`crate::metrics::brokerage`].
//!
//! Conventions:
//! - Alters in directed networks are chosen by `mode`: `Out` takes the
//!   nodes ego reaches along outgoing arcs (successors), `In` the nodes that
//!   reach ego (predecessors), and `All` ignores direction. This is igraph's
//!   `ego()`/`make_ego_graph(mode = ...)`; NetworkX's `ego_graph()` uses
//!   successors, or both directions with `undirected = True`. Undirected
//!   networks ignore the mode.
//! - Distances count steps, so weights are ignored by every function here.
//! - The ego network is the subgraph induced by ego and its alters: every
//!   tie among them counts, whatever its direction and whatever `mode` chose
//!   the alters. Ties keep their direction, so in a directed network an
//!   alter pair joined both ways contributes two ties.
//! - Ties among alters are distinct ties between distinct nodes: parallel
//!   edges count once and self-loops are ignored, the adjacency-matrix
//!   reading that UCINET uses. Ego betweenness is computed on the same
//!   binary ties, so parallel edges are not distinct shortest paths.
//! - Ego betweenness counts each unordered pair of alters once in an
//!   undirected network (as betweenness does) and each ordered pair in a
//!   directed one. It is not normalized; its maximum is the number of pairs.
//! - An ego without alters has size 0, ties 0, pairs 0, components 0, ego
//!   betweenness 0, and density NaN. One alter gives density NaN (no
//!   pairs).
//! - `mindist` in [`neighborhood`] drops nodes closer than `mindist` steps,
//!   as igraph's `mindist`; `mindist = 1` drops ego itself.

use crate::bitset::{batches_pay, for_each_bit, BitMatrix, MultiBfs, RowPaths, Step, View, BATCH};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;
use crate::parallel;
use crate::union_find::UnionFind;

/// Which nodes form a neighborhood.
#[derive(Debug, Clone, Copy)]
pub struct NeighborhoodOptions {
    /// Largest number of steps from ego.
    pub order: usize,
    /// Direction of the steps in directed graphs.
    pub mode: DegreeMode,
    /// Smallest number of steps from ego; 0 keeps ego, 1 drops it.
    pub mindist: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for NeighborhoodOptions {
    fn default() -> Self {
        NeighborhoodOptions {
            order: 1,
            mode: DegreeMode::All,
            mindist: 0,
            threads: 0,
        }
    }
}

/// Options for [`ego_summary`].
#[derive(Debug, Clone, Copy)]
pub struct EgoOptions {
    /// Alters are the nodes within this many steps of ego (at least 1).
    pub order: usize,
    /// Direction of the steps in directed graphs.
    pub mode: DegreeMode,
    /// Whether to compute ego betweenness, the costliest measure.
    pub betweenness: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for EgoOptions {
    fn default() -> Self {
        EgoOptions {
            order: 1,
            mode: DegreeMode::All,
            betweenness: true,
            threads: 0,
        }
    }
}

/// Egocentric measures of every node, in node order. See the module
/// documentation for definitions.
#[derive(Debug, Clone, PartialEq)]
pub struct EgoSummary {
    /// Number of alters.
    pub size: Vec<f64>,
    /// Ties among alters.
    pub ties: Vec<f64>,
    /// Pairs (ordered, if directed) of alters.
    pub pairs: Vec<f64>,
    /// `ties / pairs`, NaN without pairs.
    pub density: Vec<f64>,
    /// Weak components among the alters once ego is removed.
    pub components: Vec<f64>,
    /// Betweenness of ego within its ego network, if requested.
    pub betweenness: Option<Vec<f64>>,
}

fn check_ego(g: &Graph, ego: usize) -> Result<()> {
    if ego >= g.node_count() {
        return Err(Error::InvalidArgument(format!(
            "node index {ego} is out of range for a network with {} nodes",
            g.node_count()
        )));
    }
    Ok(())
}

/// Per-thread breadth-first search state, reused across egos.
struct Bfs {
    /// Steps from the current ego; `u32::MAX` if not reached.
    dist: Vec<u32>,
    /// Reached nodes in order of discovery; doubles as the queue.
    order: Vec<u32>,
}

impl Bfs {
    fn new(n: usize) -> Bfs {
        Bfs {
            dist: vec![u32::MAX; n],
            order: Vec::new(),
        }
    }

    /// Fills `order` with ego and every node within `order_max` steps, in
    /// order of discovery (so in non-decreasing distance).
    fn run(&mut self, g: &Graph, ego: usize, order_max: usize, mode: DegreeMode) {
        self.order.clear();
        self.dist[ego] = 0;
        self.order.push(ego as u32);
        let directed = g.is_directed();
        let limit = order_max.min(u32::MAX as usize - 1) as u32;
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let d = self.dist[v];
            if d >= limit {
                continue;
            }
            let follow_out = !directed || mode != DegreeMode::In;
            let follow_in = directed && mode != DegreeMode::Out;
            if follow_out {
                for &w in g.out().neighbors(v) {
                    if self.dist[w as usize] == u32::MAX {
                        self.dist[w as usize] = d + 1;
                        self.order.push(w);
                    }
                }
            }
            if follow_in {
                for &w in g.inc().neighbors(v) {
                    if self.dist[w as usize] == u32::MAX {
                        self.dist[w as usize] = d + 1;
                        self.order.push(w);
                    }
                }
            }
        }
    }

    fn reset(&mut self) {
        for &v in &self.order {
            self.dist[v as usize] = u32::MAX;
        }
    }
}

/// The nodes within `opts.order` steps of `ego` and at least
/// `opts.mindist` steps away, in ascending node order.
pub fn neighborhood(g: &Graph, ego: usize, opts: NeighborhoodOptions) -> Result<Vec<u32>> {
    check_ego(g, ego)?;
    let mut bfs = Bfs::new(g.node_count());
    bfs.run(g, ego, opts.order, opts.mode);
    let mut out: Vec<u32> = bfs
        .order
        .iter()
        .copied()
        .filter(|&v| bfs.dist[v as usize] as usize >= opts.mindist)
        .collect();
    out.sort_unstable();
    Ok(out)
}

/// [`neighborhood`] for several egos, in parallel, in the order given.
pub fn neighborhoods(
    g: &Graph,
    egos: &[usize],
    opts: NeighborhoodOptions,
) -> Result<Vec<Vec<u32>>> {
    for &e in egos {
        check_ego(g, e)?;
    }
    let n = g.node_count();
    let ball = parallel::ball_steps(n, g.out().arc_count(), opts.order);
    let threads = parallel::threads_for(
        opts.threads,
        egos.len(),
        (egos.len() as u64).saturating_mul(ball),
    );
    let partials = parallel::for_each_item(
        egos.len(),
        threads,
        || (Bfs::new(n), Vec::new()),
        |(bfs, found): &mut (Bfs, Vec<(usize, Vec<u32>)>), i| {
            bfs.run(g, egos[i], opts.order, opts.mode);
            let mut members: Vec<u32> = bfs
                .order
                .iter()
                .copied()
                .filter(|&v| bfs.dist[v as usize] as usize >= opts.mindist)
                .collect();
            members.sort_unstable();
            bfs.reset();
            found.push((i, members));
        },
    );
    let mut out = vec![Vec::new(); egos.len()];
    for (_, found) in partials {
        for (i, members) in found {
            out[i] = members;
        }
    }
    Ok(out)
}

/// Neighborhood sizes for the given egos (every node if `None`), in
/// parallel. This is igraph's `ego_size()`: with `mindist = 0` ego itself
/// is counted.
pub fn neighborhood_sizes(
    g: &Graph,
    egos: Option<&[usize]>,
    opts: NeighborhoodOptions,
) -> Result<Vec<usize>> {
    let n = g.node_count();
    if let Some(list) = egos {
        for &e in list {
            check_ego(g, e)?;
        }
    }
    let count = egos.map_or(n, |l| l.len());
    let sources: Vec<u32> = match egos {
        Some(l) => l.iter().map(|&e| e as u32).collect(),
        None => (0..n as u32).collect(),
    };
    // Breadth-first search from 64 egos at a time
    // (`crate::bitset::MultiBfs`) over the arcs the mode follows; one
    // batch is one work item, so the result does not depend on the thread
    // count.
    let view = View::new(g, opts.mode);
    let max_level = opts.order.min(u32::MAX as usize - 1) as u32;
    let ball = parallel::ball_steps(n, view.arc_count(), opts.order);
    let threads = parallel::threads_for(opts.threads, count, (count as u64).saturating_mul(ball));
    if !batches_pay(view, max_level) {
        // Small neighborhoods that rarely overlap: one search per ego.
        let partials = parallel::for_each_item(
            count,
            threads,
            || (Bfs::new(n), Vec::new()),
            |(bfs, found): &mut (Bfs, Vec<(usize, usize)>), i| {
                bfs.run(g, sources[i] as usize, opts.order, opts.mode);
                let size = bfs
                    .order
                    .iter()
                    .filter(|&&v| bfs.dist[v as usize] as usize >= opts.mindist)
                    .count();
                bfs.reset();
                found.push((i, size));
            },
        );
        let mut out = vec![0usize; count];
        for (_, found) in partials {
            for (i, size) in found {
                out[i] = size;
            }
        }
        return Ok(out);
    }
    let partials = parallel::for_each_item(
        (count + BATCH - 1) / BATCH,
        threads,
        || (MultiBfs::new(n), Vec::new()),
        |(ms, found): &mut (MultiBfs, Vec<(usize, [usize; BATCH])>), b| {
            let lo = b * BATCH;
            let hi = (lo + BATCH).min(count);
            let mut sizes = [usize::from(opts.mindist == 0); BATCH];
            ms.run(view, &sources[lo..hi], max_level, |step| {
                if let Step::Reached { level, sources, .. } = step {
                    if level as usize >= opts.mindist {
                        for_each_bit(sources, |i| sizes[i] += 1);
                    }
                }
            });
            found.push((lo, sizes));
        },
    );
    let mut out = vec![0usize; count];
    for (_, found) in partials {
        for (lo, sizes) in found {
            let hi = (lo + BATCH).min(count);
            out[lo..hi].copy_from_slice(&sizes[..hi - lo]);
        }
    }
    Ok(out)
}

/// Largest first-order ego network whose ego betweenness is computed from
/// pairs of alters (its bit rows then take at most 8 MB).
const PAIR_MAX: usize = 8192;

/// Per-thread state for [`ego_summary`].
struct EgoWork {
    bfs: Bfs,
    /// Position of each node in the current ego network; `u32::MAX` if
    /// outside it.
    local: Vec<u32>,
    offsets: Vec<usize>,
    adj: Vec<u32>,
    sets: UnionFind,
    hops: Vec<u32>,
    sigma: Vec<f64>,
    delta: Vec<f64>,
    queue: Vec<u32>,
    /// Bit rows of the ego network's adjacency, for
    /// [`EgoWork::ego_betweenness_pairs`].
    rows: Vec<u64>,
    /// The ego network's arcs and their transpose as bit matrices, and
    /// search scratch space, for [`EgoWork::ego_betweenness_rows`].
    out_bits: BitMatrix,
    in_bits: BitMatrix,
    paths: RowPaths,
    /// Distance and number of shortest paths from ego, and to ego.
    from_ego: Vec<(u32, f64)>,
    to_ego: Vec<(u32, f64)>,
    /// The nodes at each distance from ego, as bits.
    ego_levels: Vec<u64>,
    /// (ego, [size, ties, pairs, density, components, betweenness]).
    found: Vec<(u32, [f64; 6])>,
}

impl EgoWork {
    fn new(n: usize) -> EgoWork {
        EgoWork {
            bfs: Bfs::new(n),
            local: vec![u32::MAX; n],
            offsets: Vec::new(),
            adj: Vec::new(),
            sets: UnionFind::new(0),
            hops: Vec::new(),
            sigma: Vec::new(),
            delta: Vec::new(),
            queue: Vec::new(),
            rows: Vec::new(),
            out_bits: BitMatrix::new(0),
            in_bits: BitMatrix::new(0),
            paths: RowPaths::new(),
            from_ego: Vec::new(),
            to_ego: Vec::new(),
            ego_levels: Vec::new(),
            found: Vec::new(),
        }
    }

    fn measure(&mut self, g: &Graph, ego: usize, opts: &EgoOptions) {
        self.bfs.run(g, ego, opts.order, opts.mode);
        // Ego is local node 0; alters follow in order of discovery.
        let members = &self.bfs.order;
        let s = members.len();
        for (i, &v) in members.iter().enumerate() {
            self.local[v as usize] = i as u32;
        }

        // Binary adjacency of the induced ego network: outgoing arcs (both
        // directions for undirected graphs), distinct, without loops.
        self.offsets.clear();
        self.adj.clear();
        self.offsets.push(0);
        for (a, &v) in members.iter().enumerate() {
            let start = self.adj.len();
            for &w in g.out().neighbors(v as usize) {
                let b = self.local[w as usize];
                if b != u32::MAX && b as usize != a {
                    self.adj.push(b);
                }
            }
            self.adj[start..].sort_unstable();
            let mut keep = start;
            for k in start..self.adj.len() {
                if k == start || self.adj[k] != self.adj[keep - 1] {
                    self.adj[keep] = self.adj[k];
                    keep += 1;
                }
            }
            self.adj.truncate(keep);
            self.offsets.push(self.adj.len());
        }

        let k = (s - 1) as f64;
        let directed = g.is_directed();
        let mut arcs = 0usize;
        for a in 1..s {
            for &b in &self.adj[self.offsets[a]..self.offsets[a + 1]] {
                if b != 0 {
                    arcs += 1;
                }
            }
        }
        let ties = if directed { arcs } else { arcs / 2 } as f64;
        let pairs = if directed {
            k * (k - 1.0)
        } else {
            k * (k - 1.0) / 2.0
        };
        let density = if pairs > 0.0 { ties / pairs } else { f64::NAN };

        // Weak components among alters by union-find (ego excluded).
        self.sets.reset(s);
        let mut components = s.saturating_sub(1);
        for a in 1..s {
            for idx in self.offsets[a]..self.offsets[a + 1] {
                let b = self.adj[idx];
                if b == 0 {
                    continue;
                }
                if self.sets.union(a, b as usize).is_some() {
                    components -= 1;
                }
            }
        }

        let betweenness = if !opts.betweenness || s <= 2 {
            0.0
        } else if !directed && opts.order == 1 && s <= PAIR_MAX && s * s < 256 * self.adj.len() {
            self.ego_betweenness_pairs(s)
        } else if s <= PAIR_MAX && s * ((s + 63) / 64) + 32 < self.adj.len() {
            // A search on bit rows costs about s / 64 words per node plus a
            // fixed overhead, against one step per arc for a list-based
            // search; measured break-even on G(n, p) ego networks.
            self.ego_betweenness_rows(s, directed)
        } else {
            self.ego_betweenness(s, directed)
        };

        for &v in &self.bfs.order {
            self.local[v as usize] = u32::MAX;
        }
        self.bfs.reset();
        self.found.push((
            ego as u32,
            [k, ties, pairs, density, components as f64, betweenness],
        ));
    }

    /// Ego betweenness of a first-order undirected ego network from its
    /// pairs of alters (Everett and Borgatti 2005, `everett2005ego`):
    /// every alter is adjacent to ego, so two non-adjacent alters i and j
    /// are at distance 2, their geodesics are their common neighbours in
    /// the ego network (ego among them), and ego's share is 1 / (A^2)_ij;
    /// adjacent alters contribute nothing. Rows of the ego network's
    /// adjacency are held as bits, so each pair costs s / 64 word ANDs:
    /// O(s^3 / 64) in all, against O(s m) for Brandes' accumulation, which
    /// makes it the cheaper choice on dense ego networks.
    fn ego_betweenness_pairs(&mut self, s: usize) -> f64 {
        let words = (s + 63) / 64;
        self.rows.clear();
        self.rows.resize(s * words, 0);
        for a in 0..s {
            let row = &mut self.rows[a * words..(a + 1) * words];
            for &b in &self.adj[self.offsets[a]..self.offsets[a + 1]] {
                row[b as usize / 64] |= 1 << (b % 64);
            }
        }
        pair_shares(&self.rows, words, s)
    }

    /// Ego betweenness of any ego network, directed or not, from the
    /// distances and shortest-path counts through ego.
    ///
    /// Ego (local node 0) lies on a shortest path from i to j exactly when
    /// d(i, ego) + d(ego, j) = d(i, j), and then
    /// sigma_ij(ego) = sigma(i, ego) sigma(ego, j) (Brandes 2001,
    /// `brandes2001faster`, section 3), so ego's share of the pair is
    /// sigma(i, ego) sigma(ego, j) / sigma_ij. One search from ego along
    /// the arcs and one against them give d(ego, .), sigma(ego, .),
    /// d(., ego), and sigma(., ego). Since d(i, j) <= d(i, ego) + d(ego, j)
    /// always, the search from each source i only needs to run to depth
    /// d(i, ego) + max_j d(ego, j), and j counts exactly when it is found
    /// at depth d(i, ego) + d(ego, j); sources that cannot reach ego are
    /// skipped. In a first-order ego network with `mode = "out"` every
    /// d(ego, j) is 1, so a search stops one level past ego. Searches run
    /// on bit rows ([`RowPaths`]), which costs O(s^2 / 64) word operations
    /// per source instead of O(m) for a list-based search with Brandes'
    /// backward pass. The sum is the same as [`EgoWork::ego_betweenness`]
    /// up to rounding.
    fn ego_betweenness_rows(&mut self, s: usize, directed: bool) -> f64 {
        self.out_bits.reset(s);
        self.in_bits.reset(s);
        for a in 0..s {
            for &b in &self.adj[self.offsets[a]..self.offsets[a + 1]] {
                self.out_bits.set(a, b as usize);
                self.in_bits.set(b as usize, a);
            }
        }
        let words = self.out_bits.words_per_row();
        self.paths
            .run(&self.out_bits, &self.in_bits, s, 0, u32::MAX);
        let reach = self.paths.depth() - 1;
        self.from_ego.clear();
        self.from_ego
            .extend((0..s).map(|v| (self.paths.level(v), self.paths.sigma(v))));
        self.ego_levels.clear();
        for l in 0..=reach {
            self.ego_levels.extend_from_slice(self.paths.level_bits(l));
        }
        self.paths
            .run(&self.in_bits, &self.out_bits, s, 0, u32::MAX);
        self.to_ego.clear();
        self.to_ego
            .extend((0..s).map(|v| (self.paths.level(v), self.paths.sigma(v))));
        let mut total = 0.0;
        for src in 1..s {
            let (h, through) = self.to_ego[src];
            if h == u32::MAX || reach == 0 {
                continue;
            }
            let h = h as usize;
            // The nodes ego can be on the way to: those at depth
            // d(src, ego) + d(ego, j) for some j, which are the only
            // counts needed.
            let ego_levels = &self.ego_levels;
            self.paths.run_for(
                &self.out_bits,
                &self.in_bits,
                s,
                src,
                (h + reach) as u32,
                |l, w| {
                    if l > h && l - h <= reach {
                        ego_levels[(l - h) * words + w]
                    } else {
                        0
                    }
                },
            );
            for l in h + 1..=h + reach {
                for &j in self.paths.nodes_at(l) {
                    let (d, onward) = self.from_ego[j as usize];
                    if d as usize == l - h {
                        total += through * onward / self.paths.sigma(j as usize);
                    }
                }
            }
        }
        if directed {
            total
        } else {
            total / 2.0
        }
    }

    /// Dependency of local node 0 (ego) summed over every other source,
    /// by breadth-first Brandes accumulation over successors.
    fn ego_betweenness(&mut self, s: usize, directed: bool) -> f64 {
        self.hops.clear();
        self.hops.resize(s, u32::MAX);
        self.sigma.clear();
        self.sigma.resize(s, 0.0);
        self.delta.clear();
        self.delta.resize(s, 0.0);
        let mut total = 0.0;
        for src in 1..s {
            self.queue.clear();
            self.queue.push(src as u32);
            self.hops[src] = 0;
            self.sigma[src] = 1.0;
            let mut head = 0;
            while head < self.queue.len() {
                let v = self.queue[head] as usize;
                head += 1;
                let next = self.hops[v] + 1;
                for &w in &self.adj[self.offsets[v]..self.offsets[v + 1]] {
                    let w = w as usize;
                    if self.hops[w] == u32::MAX {
                        self.hops[w] = next;
                        self.queue.push(w as u32);
                    }
                    if self.hops[w] == next {
                        self.sigma[w] += self.sigma[v];
                    }
                }
            }
            for i in (0..self.queue.len()).rev() {
                let v = self.queue[i] as usize;
                let next = self.hops[v] + 1;
                let mut dv = 0.0;
                for &w in &self.adj[self.offsets[v]..self.offsets[v + 1]] {
                    let w = w as usize;
                    if self.hops[w] == next {
                        dv += (1.0 + self.delta[w]) / self.sigma[w];
                    }
                }
                self.delta[v] = self.sigma[v] * dv;
            }
            total += self.delta[0];
            for &v in &self.queue {
                let v = v as usize;
                self.hops[v] = u32::MAX;
                self.sigma[v] = 0.0;
                self.delta[v] = 0.0;
            }
        }
        if directed {
            total
        } else {
            total / 2.0
        }
    }
}

/// Sum over non-adjacent pairs of alters (local nodes 1..s) of
/// 1 / (number of common neighbours), from bit rows of `words` words.
/// Population counts dominate, so on x86-64 processors with the POPCNT
/// instruction this runs through a copy compiled to use it; both copies
/// perform the same operations on the same integers.
fn pair_shares(rows: &[u64], words: usize, s: usize) -> f64 {
    #[inline(always)]
    fn body(rows: &[u64], words: usize, s: usize) -> f64 {
        let mut total = 0.0;
        for i in 1..s {
            let ri = &rows[i * words..(i + 1) * words];
            for j in i + 1..s {
                if ri[j / 64] >> (j % 64) & 1 == 1 {
                    continue;
                }
                let rj = &rows[j * words..(j + 1) * words];
                let common: u32 = ri.iter().zip(rj).map(|(x, y)| (x & y).count_ones()).sum();
                total += 1.0 / common as f64;
            }
        }
        total
    }
    #[cfg(target_arch = "x86_64")]
    {
        #[target_feature(enable = "popcnt")]
        unsafe fn fast(rows: &[u64], words: usize, s: usize) -> f64 {
            body(rows, words, s)
        }
        if std::is_x86_feature_detected!("popcnt") {
            // SAFETY: the processor supports POPCNT.
            return unsafe { fast(rows, words, s) };
        }
    }
    body(rows, words, s)
}

/// Egocentric measures of every node. See the module documentation.
pub fn ego_summary(g: &Graph, opts: EgoOptions) -> Result<EgoSummary> {
    if opts.order == 0 {
        return Err(Error::InvalidArgument(
            "the order of an ego network must be at least 1".into(),
        ));
    }
    let n = g.node_count();
    // Building an ego network scans the arcs of its members; ego
    // betweenness then costs about one more step per member and arc for
    // every member.
    let ball = parallel::ball_steps(n, g.out().arc_count(), opts.order);
    let size = (ball as f64).sqrt() as u64 + 1;
    let per_ego = if opts.betweenness {
        ball.saturating_mul(size)
    } else {
        ball
    };
    let threads = parallel::threads_for(opts.threads, n, (n as u64).saturating_mul(per_ego));
    // Ego networks differ widely in size (hubs), so threads claim egos;
    // every result is stored by ego.
    let partials = parallel::for_each_claimed(
        n,
        threads,
        || EgoWork::new(n),
        |ws, ego| ws.measure(g, ego, &opts),
    );
    let mut out = EgoSummary {
        size: vec![0.0; n],
        ties: vec![0.0; n],
        pairs: vec![0.0; n],
        density: vec![f64::NAN; n],
        components: vec![0.0; n],
        betweenness: opts.betweenness.then(|| vec![0.0; n]),
    };
    for ws in partials {
        for (ego, m) in ws.found {
            let e = ego as usize;
            out.size[e] = m[0];
            out.ties[e] = m[1];
            out.pairs[e] = m[2];
            out.density[e] = m[3];
            out.components[e] = m[4];
            if let Some(b) = out.betweenness.as_mut() {
                b[e] = m[5];
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairwise_ego_betweenness_equals_brandes() {
        // Undirected networks with parallel edges and loops, sparse to
        // dense, so that both formulas see every kind of ego network.
        let mut rng = crate::rng::Stream::new(7, 0);
        for &(n, m) in &[(12usize, 20usize), (30, 200), (40, 700), (25, 400)] {
            let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let g = Graph::from_edges(n, &s, &t, None, false).unwrap();
            let opts = EgoOptions::default();
            let mut w = EgoWork::new(n);
            for ego in 0..n {
                // `measure` leaves the ego network's adjacency in place.
                w.measure(&g, ego, &opts);
                let size = w.offsets.len() - 1;
                if size <= 2 {
                    continue;
                }
                let (a, b) = (
                    w.ego_betweenness_pairs(size),
                    w.ego_betweenness(size, false),
                );
                assert!(
                    (a - b).abs() <= 1e-12 * b.abs().max(1.0),
                    "ego {ego}: {a} vs {b}"
                );
            }
        }
    }

    #[test]
    fn row_ego_betweenness_equals_brandes() {
        // Directed and undirected networks with parallel edges and loops,
        // every mode, first and second order, sparse to dense.
        let mut rng = crate::rng::Stream::new(11, 0);
        for &(n, m) in &[
            (12usize, 20usize),
            (30, 120),
            (40, 700),
            (25, 400),
            (70, 2500),
        ] {
            for directed in [true, false] {
                let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
                let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
                let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
                for mode in [DegreeMode::Out, DegreeMode::In, DegreeMode::All] {
                    for order in [1, 2] {
                        let opts = EgoOptions {
                            order,
                            mode,
                            ..EgoOptions::default()
                        };
                        let mut w = EgoWork::new(n);
                        for ego in 0..n {
                            w.measure(&g, ego, &opts);
                            let size = w.offsets.len() - 1;
                            if size <= 2 {
                                continue;
                            }
                            let (a, b) = (
                                w.ego_betweenness_rows(size, directed),
                                w.ego_betweenness(size, directed),
                            );
                            assert!(
                                (a - b).abs() <= 1e-12 * b.abs().max(1.0),
                                "n {n} m {m} {mode:?} order {order} ego {ego}: {a} vs {b}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn row_paths_count_shortest_paths() {
        // A directed 3 x 3 grid from corner 0 has binomial path counts, and
        // a long path forces many levels across word boundaries.
        let mut out = BitMatrix::new(9);
        let mut inc = BitMatrix::new(9);
        for r in 0..3 {
            for c in 0..3 {
                let v = 3 * r + c;
                for w in [(c < 2).then(|| v + 1), (r < 2).then(|| v + 3)]
                    .into_iter()
                    .flatten()
                {
                    out.set(v, w);
                    inc.set(w, v);
                }
            }
        }
        let mut p = RowPaths::new();
        p.run(&out, &inc, 9, 0, u32::MAX);
        assert_eq!(p.depth(), 5);
        assert_eq!(p.sigma(8), 6.0);
        assert_eq!(p.sigma(4), 2.0);
        assert_eq!(p.nodes_at(2), &[2, 4, 6]);
        p.run(&out, &inc, 9, 4, 1);
        assert_eq!(
            (p.level(5), p.level(8), p.level(0)),
            (1, u32::MAX, u32::MAX)
        );
        let n = 150;
        let mut out = BitMatrix::new(n);
        let mut inc = BitMatrix::new(n);
        for v in 0..n - 1 {
            out.set(v, v + 1);
            inc.set(v + 1, v);
        }
        p.run(&out, &inc, n, 0, u32::MAX);
        assert_eq!(
            (p.depth(), p.level(n - 1), p.sigma(n - 1)),
            (n, n as u32 - 1, 1.0)
        );
        p.run(&out, &inc, n, n - 1, u32::MAX);
        assert_eq!((p.depth(), p.level(0)), (1, u32::MAX));
    }

    fn g(n: usize, s: &[i64], t: &[i64], directed: bool) -> Graph {
        Graph::from_edges(n, s, t, None, directed).unwrap()
    }

    #[test]
    fn neighborhoods_follow_mode_and_order() {
        // 0 -> 1 -> 2 -> 3, 4 -> 0.
        let d = g(5, &[0, 1, 2, 4], &[1, 2, 3, 0], true);
        let o = |order, mode, mindist| NeighborhoodOptions {
            order,
            mode,
            mindist,
            threads: 1,
        };
        assert_eq!(
            neighborhood(&d, 0, o(1, DegreeMode::Out, 0)).unwrap(),
            vec![0, 1]
        );
        assert_eq!(
            neighborhood(&d, 0, o(1, DegreeMode::In, 0)).unwrap(),
            vec![0, 4]
        );
        assert_eq!(
            neighborhood(&d, 0, o(1, DegreeMode::All, 1)).unwrap(),
            vec![1, 4]
        );
        assert_eq!(
            neighborhood(&d, 0, o(2, DegreeMode::Out, 2)).unwrap(),
            vec![2]
        );
        assert_eq!(
            neighborhood_sizes(&d, None, o(2, DegreeMode::All, 0)).unwrap(),
            vec![4, 5, 4, 3, 3]
        );
        let many = neighborhoods(&d, &[3, 0], o(1, DegreeMode::In, 0)).unwrap();
        assert_eq!(many, vec![vec![2, 3], vec![0, 4]]);
        assert!(neighborhood(&d, 5, o(1, DegreeMode::All, 0)).is_err());
    }

    #[test]
    fn kite_like_ego_measures() {
        // Ego 0 with alters 1..=4; ties 1-2 (twice, parallel) and 3-4,
        // plus a loop at 1 and an outside node 5 tied to 4.
        let u = g(
            6,
            &[0, 0, 0, 0, 1, 1, 3, 1, 4],
            &[1, 2, 3, 4, 2, 2, 4, 1, 5],
            false,
        );
        let s = ego_summary(&u, EgoOptions::default()).unwrap();
        assert_eq!(s.size[0], 4.0);
        assert_eq!(s.ties[0], 2.0);
        assert_eq!(s.pairs[0], 6.0);
        assert_eq!(s.density[0], 2.0 / 6.0);
        assert_eq!(s.components[0], 2.0);
        // Ego brokers the four pairs across the two dyads.
        assert_eq!(s.betweenness.as_ref().unwrap()[0], 4.0);
        // Node 5 has one alter: no pairs.
        assert_eq!(s.size[5], 1.0);
        assert!(s.density[5].is_nan());
        assert_eq!(s.components[5], 1.0);
    }

    #[test]
    fn star_ego_betweenness_is_the_number_of_pairs() {
        let star = g(5, &[0, 0, 0, 0], &[1, 2, 3, 4], false);
        let s = ego_summary(&star, EgoOptions::default()).unwrap();
        assert_eq!(s.betweenness.unwrap(), vec![6.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(s.components[0], 4.0);
        let d = g(3, &[1, 0], &[0, 2], true);
        let s = ego_summary(&d, EgoOptions::default()).unwrap();
        // Only 1 -> 0 -> 2 passes through 0.
        assert_eq!(s.betweenness.as_ref().unwrap()[0], 1.0);
        assert_eq!(s.pairs[0], 2.0);
        let out = EgoOptions {
            mode: DegreeMode::Out,
            ..Default::default()
        };
        let s = ego_summary(&d, out).unwrap();
        assert_eq!(s.size, vec![1.0, 1.0, 0.0]);
        assert!(s.density[2].is_nan());
        assert_eq!(s.components[2], 0.0);
    }

    #[test]
    fn rejects_order_zero() {
        let u = g(2, &[0], &[1], false);
        let bad = EgoOptions {
            order: 0,
            ..Default::default()
        };
        assert!(ego_summary(&u, bad).is_err());
    }
}
