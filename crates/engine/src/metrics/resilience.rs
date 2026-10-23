//! Key players, group centrality, and flow betweenness: the measures of
//! node importance used to plan the disruption of covert networks (Everton
//! 2012, `everton2012disrupting`, chapter 7).
//!
//! Sources:
//! - Key players (Borgatti 2006, `borgatti2006identifying`). KPP-Neg asks
//!   which k nodes, removed, most fragment the network; KPP-Pos asks which
//!   k nodes reach most of it. With S the set and n the number of nodes of
//!   the network being measured:
//!   - fragmentation of G - S, eq. 3: F = 1 - sum_{i != j} r_ij / (n(n - 1)),
//!     with r_ij = 1 if i reaches j;
//!   - distance-weighted fragmentation of G - S, eq. 9:
//!     DF = 1 - sum_{i != j} (1 / d_ij) / (n(n - 1)), with 1/inf = 0;
//!   - m-reach, eq. 12: the number of nodes outside S within distance m of
//!     some member of S (a count, not normalized, as in the paper);
//!   - distance-weighted reach, eq. 14: DR = sum_j (1 / d_Sj) / n, where
//!     d_Sj is the distance from the nearest member of S and members count
//!     as distance 1 ("regard all distances within the KP set to be
//!     unity", p. 29).
//!     The paper writes F and DF over unordered pairs of an undirected graph
//!     with a factor 2; the ordered-pair forms above are identical there and
//!     extend them to directed paths. Optimization is Borgatti's greedy swap
//!     heuristic (Fig. 6): from a random set, make the single swap of a
//!     member for a non-member that most improves the fit, until no swap
//!     improves it; restarts run in parallel from independent random streams
//!     and the best result wins. Exact enumeration of every set is available
//!     for small problems.
//! - Group degree, closeness, and betweenness (Everett and Borgatti 1999,
//!   `everett1999centrality`). Group betweenness counts, for every pair of
//!   non-members, the fraction of shortest paths that pass through at least
//!   one member, from path counts with and without the members (one
//!   search per source, as in Brandes 2008, `brandes2008variants`).
//! - Flow betweenness (Freeman, Borgatti and White 1991,
//!   `freeman1991centrality`): for each node v, the sum over pairs of other
//!   nodes i, j of f(i, j) - f_v(i, j), where f is the maximum flow from i
//!   to j and f_v the maximum flow with v removed; that is, the part of
//!   the maximum flow that cannot avoid v. "Normalized" divides by the
//!   total maximum flow between pairs not involving v (the standardization
//!   UCINET reports; sna `normflow`), and "fractional" sums
//!   (f - f_v) / f over pairs with f > 0 (Koschutzki et al. 2005,
//!   `koschutzki2005centrality`; sna `fracflow`). Maximum flows use
//!   Dinic's blocking-flow algorithm (Dinic 1970, `dinic1970algorithm`);
//!   by the max-flow min-cut theorem (Ford and Fulkerson 1956,
//!   `ford1956maximal`) f(i, j) is also the capacity of a minimum i-j cut.
//!
//! Conventions:
//! - Self-loops never lie on a path or carry flow and are ignored
//!   throughout. Parallel edges do not change distances; they count as
//!   distinct shortest paths in group betweenness (as in [`betweenness`])
//!   and add their capacities in flows.
//! - Key players and group closeness and betweenness are path-based:
//!   weights, when used, are distances, or strengths converted to lengths
//!   1/w ([`crate::weights`]), and must be positive. Breadth-first search is
//!   used for binary networks and uniform weights. Group degree is defined
//!   on binary ties and ignores weights.
//! - Directed networks: paths follow arc directions (`mode = out`), or the
//!   reverse (`in`), or ignore them (`all`).
//! - Flows: weights, when used, are capacities and must be non-negative;
//!   without weights every edge has capacity 1. An undirected edge carries
//!   flow either way up to its capacity. Vertices have unlimited capacity.
//!   Pairs are ordered in directed networks and unordered in undirected
//!   ones (sna's `gmode = "graph"`).
//! - Ties: greedy swaps prefer the lowest node indices; key player sets
//!   tie-break to the lexicographically smallest set; scores within a
//!   relative 1e-9 of each other are ties.
//! - Undefined values are NaN: fragmentation of a network with fewer than
//!   two nodes, normalized group betweenness with fewer than two
//!   non-members, normalized flow betweenness with no flow between other
//!   pairs.
//! - Cost: KPP-Neg evaluates all pairs of the remaining network for every
//!   candidate swap (O(k (n - k) n m) per greedy step); flow betweenness
//!   needs a maximum flow for every pair and, for every node carrying flow
//!   between that pair, another with the node removed: up to O(n^3)
//!   maximum flows, run in parallel over sources.
//!
//! [`betweenness`]: crate::metrics::betweenness::betweenness

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;
use crate::parallel;
use crate::rng::Stream;
use crate::union_find::UnionFind;
use crate::weights::{additive_weights, path_lengths, Lengths, PathAlgorithm, WeightKind};

const EPSILON: f64 = 1e-10;
const TIE: f64 = 1e-9;

#[inline]
fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPSILON * a.abs().max(b.abs())
}

/// Whether `a` beats `b` by more than the tie tolerance.
#[inline]
fn better(a: f64, b: f64) -> bool {
    a > b && (a - b) > TIE * a.abs().max(b.abs())
}

fn check_nodes(n: usize, nodes: &[u32], what: &str) -> Result<()> {
    if let Some(&bad) = nodes.iter().find(|&&v| v as usize >= n) {
        return Err(Error::InvalidArgument(format!(
            "{what} node index {bad} is out of range for a network with {n} nodes"
        )));
    }
    let mut seen = vec![false; n];
    for &v in nodes {
        if seen[v as usize] {
            return Err(Error::InvalidArgument(format!(
                "{what} node index {v} appears more than once"
            )));
        }
        seen[v as usize] = true;
    }
    Ok(())
}

// ---- Shortest paths on a network with removed nodes -------------------------

/// A graph oriented for a path mode, with arc lengths. `unit` is the common
/// length when every weight is equal, since breadth-first search counts hops.
struct PathGraph {
    g: Graph,
    lengths: Lengths,
    unit: f64,
}

fn path_graph(g: &Graph, mode: DegreeMode, weighted: bool, kind: WeightKind) -> Result<PathGraph> {
    let h = if !g.is_directed() || mode == DegreeMode::Out {
        g.clone()
    } else {
        let w = g.edge_weights();
        let s: Vec<i64> = g.edge_sources().iter().map(|&x| x as i64).collect();
        let t: Vec<i64> = g.edge_targets().iter().map(|&x| x as i64).collect();
        match mode {
            DegreeMode::In => Graph::from_edges(g.node_count(), &t, &s, w, true)?,
            _ => Graph::from_edges(g.node_count(), &s, &t, w, false)?,
        }
    };
    let lengths = path_lengths(&h, weighted, kind, PathAlgorithm::Auto)?;
    let unit = match (&lengths, weighted) {
        (Lengths::Uniform, true) => {
            h.out()
                .all_weights()
                .and_then(|w| w.first())
                .map_or(1.0, |&w| match kind {
                    WeightKind::Distance => w,
                    WeightKind::Strength => 1.0 / w,
                })
        }
        _ => 1.0,
    };
    Ok(PathGraph {
        g: h,
        lengths,
        unit,
    })
}

#[derive(Clone, Copy, PartialEq)]
struct HeapItem {
    dist: f64,
    node: u32,
}

impl Eq for HeapItem {}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap on distance, then node index.
        other
            .dist
            .partial_cmp(&self.dist)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Scratch space for repeated searches.
struct Search {
    dist: Vec<f64>,
    order: Vec<u32>,
    settled: Vec<bool>,
    heap: BinaryHeap<HeapItem>,
}

impl Search {
    fn new(n: usize) -> Search {
        Search {
            dist: vec![f64::INFINITY; n],
            order: Vec::with_capacity(n),
            settled: vec![false; n],
            heap: BinaryHeap::new(),
        }
    }

    fn reset(&mut self) {
        for &v in &self.order {
            self.dist[v as usize] = f64::INFINITY;
            self.settled[v as usize] = false;
        }
        self.order.clear();
        self.heap.clear();
    }

    /// Distances from the nearest of `sources`, skipping removed nodes.
    /// `order` lists reached nodes by non-decreasing distance. Distances are
    /// in units of the common length for uniform lengths (multiply by
    /// `PathGraph::unit`).
    fn run(&mut self, p: &PathGraph, sources: &[u32], removed: &[bool]) {
        self.reset();
        let adj = p.g.out();
        match &p.lengths {
            Lengths::Uniform => {
                for &s in sources {
                    let s = s as usize;
                    if !removed[s] && self.dist[s].is_infinite() {
                        self.dist[s] = 0.0;
                        self.order.push(s as u32);
                    }
                }
                let mut head = 0;
                while head < self.order.len() {
                    let v = self.order[head] as usize;
                    head += 1;
                    let next = self.dist[v] + 1.0;
                    for &w in adj.neighbors(v) {
                        let w = w as usize;
                        if !removed[w] && self.dist[w].is_infinite() {
                            self.dist[w] = next;
                            self.order.push(w as u32);
                        }
                    }
                }
            }
            lens => {
                let len = |k: usize| match lens {
                    Lengths::Integer(l, _) => l[k] as f64,
                    Lengths::Real(l) => l[k],
                    Lengths::Uniform => 1.0,
                };
                let mut touched: Vec<u32> = Vec::new();
                for &s in sources {
                    let s = s as usize;
                    if !removed[s] && self.dist[s].is_infinite() {
                        self.dist[s] = 0.0;
                        touched.push(s as u32);
                        self.heap.push(HeapItem {
                            dist: 0.0,
                            node: s as u32,
                        });
                    }
                }
                while let Some(HeapItem { dist: d, node }) = self.heap.pop() {
                    let v = node as usize;
                    if self.settled[v] || d > self.dist[v] {
                        continue;
                    }
                    self.settled[v] = true;
                    self.order.push(v as u32);
                    for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                        let w = w as usize;
                        if removed[w] || self.settled[w] {
                            continue;
                        }
                        let alt = d + len(k);
                        if alt < self.dist[w] {
                            if self.dist[w].is_infinite() {
                                touched.push(w as u32);
                            }
                            self.dist[w] = alt;
                            self.heap.push(HeapItem {
                                dist: alt,
                                node: w as u32,
                            });
                        }
                    }
                }
                // Nodes reached but never settled cannot exist (every pushed
                // node is eventually popped), so `order` covers `touched`.
                debug_assert_eq!(touched.len(), self.order.len());
            }
        }
    }
}

// ---- Key players -------------------------------------------------------------

/// Which key player criterion to optimize.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KppMeasure {
    /// KPP-Neg: fragmentation F of the network without the set.
    Fragmentation,
    /// KPP-Neg: distance-weighted fragmentation DF of the network without
    /// the set.
    DistanceFragmentation,
    /// KPP-Pos: number of non-members within distance m of the set.
    MReach,
    /// KPP-Pos: distance-weighted reach DR.
    DistanceReach,
}

impl KppMeasure {
    /// Decodes the C interface codes: 0 = fragmentation, 1 = distance-
    /// weighted fragmentation, 2 = m-reach, 3 = distance-weighted reach.
    pub fn from_code(code: i32) -> Result<KppMeasure> {
        match code {
            0 => Ok(KppMeasure::Fragmentation),
            1 => Ok(KppMeasure::DistanceFragmentation),
            2 => Ok(KppMeasure::MReach),
            3 => Ok(KppMeasure::DistanceReach),
            _ => Err(Error::InvalidArgument(format!(
                "unknown key player measure code {code}"
            ))),
        }
    }
}

/// How to search for the best set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum KppMethod {
    /// Borgatti's greedy swap heuristic with random restarts.
    #[default]
    Greedy,
    /// Every set of size k (feasible only for small problems).
    Exact,
}

impl KppMethod {
    /// Decodes the C interface codes: 0 = greedy, 1 = exact.
    pub fn from_code(code: i32) -> Result<KppMethod> {
        match code {
            0 => Ok(KppMethod::Greedy),
            1 => Ok(KppMethod::Exact),
            _ => Err(Error::InvalidArgument(format!(
                "unknown key player method code {code}"
            ))),
        }
    }
}

/// Largest number of sets exact enumeration will evaluate.
pub const EXACT_MAX_SETS: u128 = 50_000_000;

/// Options for [`kpp_score`] and [`key_players`].
#[derive(Debug, Clone, Copy)]
pub struct KppOptions {
    /// The criterion.
    pub measure: KppMeasure,
    /// Distance limit for m-reach (hops when unweighted).
    pub m: f64,
    /// Direction of paths in directed networks.
    pub mode: DegreeMode,
    /// Use edge weights as path lengths.
    pub weighted: bool,
    /// Whether weights are distances or strengths.
    pub weights_as: WeightKind,
    /// Greedy search or exact enumeration.
    pub method: KppMethod,
    /// Number of greedy restarts.
    pub restarts: usize,
    /// Seed of the random starting sets.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for KppOptions {
    fn default() -> Self {
        KppOptions {
            measure: KppMeasure::Fragmentation,
            m: 1.0,
            mode: DegreeMode::Out,
            weighted: false,
            weights_as: WeightKind::Distance,
            method: KppMethod::Greedy,
            restarts: 10,
            seed: 1,
            threads: 0,
        }
    }
}

/// The best set found.
#[derive(Debug, Clone, PartialEq)]
pub struct KeyPlayers {
    /// Members, in increasing node order.
    pub set: Vec<u32>,
    /// Value of the criterion for the set.
    pub score: f64,
    /// Number of sets evaluated.
    pub evaluations: u64,
}

/// Evaluates key player criteria, reusing scratch space.
struct Scorer<'a> {
    p: &'a PathGraph,
    opts: KppOptions,
    search: Search,
    removed: Vec<bool>,
    sets: UnionFind,
    size: Vec<u32>,
}

impl<'a> Scorer<'a> {
    fn new(p: &'a PathGraph, opts: KppOptions) -> Scorer<'a> {
        let n = p.g.node_count();
        Scorer {
            p,
            opts,
            search: Search::new(n),
            removed: vec![false; n],
            sets: UnionFind::new(n),
            size: vec![0; n],
        }
    }

    fn score(&mut self, set: &[u32]) -> f64 {
        match self.opts.measure {
            KppMeasure::Fragmentation | KppMeasure::DistanceFragmentation => self.negative(set),
            KppMeasure::MReach | KppMeasure::DistanceReach => self.positive(set),
        }
    }

    fn negative(&mut self, set: &[u32]) -> f64 {
        let n = self.p.g.node_count();
        for &v in set {
            self.removed[v as usize] = true;
        }
        let rest = n - set.len();
        let out = if rest < 2 {
            f64::NAN
        } else if self.opts.measure == KppMeasure::Fragmentation && !self.p.g.is_directed() {
            self.undirected_fragmentation(rest)
        } else {
            let distance = self.opts.measure == KppMeasure::DistanceFragmentation;
            let mut total = 0.0;
            for s in 0..n {
                if self.removed[s] {
                    continue;
                }
                self.search.run(self.p, &[s as u32], &self.removed);
                if distance {
                    for &v in &self.search.order[1..] {
                        total += 1.0 / (self.search.dist[v as usize] * self.p.unit);
                    }
                } else {
                    total += (self.search.order.len() - 1) as f64;
                }
            }
            let pairs = rest as f64 * (rest as f64 - 1.0);
            1.0 - total / pairs
        };
        for &v in set {
            self.removed[v as usize] = false;
        }
        out
    }

    /// F from component sizes, exact in integers (Borgatti 2006, eq. 4).
    fn undirected_fragmentation(&mut self, rest: usize) -> f64 {
        let g = &self.p.g;
        let n = g.node_count();
        self.sets.reset(n);
        for v in self.size.iter_mut() {
            *v = 1;
        }
        for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
            if self.removed[a as usize] || self.removed[b as usize] {
                continue;
            }
            if let Some((lo, hi)) = self.sets.union(a as usize, b as usize) {
                self.size[lo] += self.size[hi];
            }
        }
        let mut within: u128 = 0;
        for v in 0..n {
            if !self.removed[v] && self.sets.is_root(v) {
                let s = self.size[v] as u128;
                within += s * (s - 1);
            }
        }
        let pairs = rest as u128 * (rest as u128 - 1);
        (pairs - within) as f64 / pairs as f64
    }

    fn positive(&mut self, set: &[u32]) -> f64 {
        let n = self.p.g.node_count();
        if n == 0 {
            return f64::NAN;
        }
        self.search.run(self.p, set, &self.removed);
        let unit = self.p.unit;
        match self.opts.measure {
            KppMeasure::MReach => {
                // A small relative slack keeps real-valued distances that
                // equal m up to rounding inside the limit.
                let limit = self.opts.m * (1.0 + EPSILON);
                self.search.order[set.len()..]
                    .iter()
                    .filter(|&&v| self.search.dist[v as usize] * unit <= limit)
                    .count() as f64
            }
            _ => {
                let members = set.len() as f64;
                let others: f64 = self.search.order[set.len()..]
                    .iter()
                    .map(|&v| 1.0 / (self.search.dist[v as usize] * unit))
                    .sum();
                (members + others) / n as f64
            }
        }
    }
}

fn kpp_prepare(g: &Graph, opts: &KppOptions) -> Result<PathGraph> {
    if opts.measure == KppMeasure::MReach && (opts.m.is_nan() || opts.m < 0.0) {
        return Err(Error::InvalidArgument(
            "the m-reach distance limit must be a non-negative number".into(),
        ));
    }
    let mode = if g.is_directed() {
        opts.mode
    } else {
        DegreeMode::Out
    };
    path_graph(g, mode, opts.weighted, opts.weights_as)
}

/// Value of a key player criterion for the set `set`. See the module
/// documentation.
pub fn kpp_score(g: &Graph, set: &[u32], opts: KppOptions) -> Result<f64> {
    check_nodes(g.node_count(), set, "key player")?;
    let p = kpp_prepare(g, &opts)?;
    let mut sorted = set.to_vec();
    sorted.sort_unstable();
    Ok(Scorer::new(&p, opts).score(&sorted))
}

/// KPP-Neg and distance-weighted reach are maximized; so is m-reach.
/// Every criterion here is "larger is better".
fn kpp_better(a: f64, b: f64) -> bool {
    if b.is_nan() {
        return !a.is_nan();
    }
    better(a, b)
}

/// Finds a set of `k` key players. See the module documentation.
pub fn key_players(g: &Graph, k: usize, opts: KppOptions) -> Result<KeyPlayers> {
    let n = g.node_count();
    if k == 0 || k > n {
        return Err(Error::InvalidArgument(format!(
            "the number of key players must be between 1 and the number of nodes ({n})"
        )));
    }
    let p = kpp_prepare(g, &opts)?;
    match opts.method {
        KppMethod::Exact => exact_key_players(&p, n, k, opts),
        KppMethod::Greedy => greedy_key_players(&p, n, k, opts),
    }
}

fn binomial(n: usize, k: usize) -> u128 {
    let k = k.min(n - k);
    let mut c: u128 = 1;
    for i in 0..k {
        c = c * (n - i) as u128 / (i + 1) as u128;
        if c > u64::MAX as u128 {
            return c;
        }
    }
    c
}

/// The `rank`-th k-subset of 0..n in lexicographic order.
fn unrank(mut rank: u128, n: usize, k: usize, out: &mut [u32]) {
    let mut start = 0;
    for (slot, item) in out.iter_mut().enumerate().take(k) {
        let left = k - slot - 1;
        let mut v = start;
        loop {
            let c = binomial(n - v - 1, left);
            if rank < c {
                break;
            }
            rank -= c;
            v += 1;
        }
        *item = v as u32;
        start = v + 1;
    }
}

fn next_combination(set: &mut [u32], n: usize) -> bool {
    let k = set.len();
    let mut i = k;
    while i > 0 {
        i -= 1;
        if (set[i] as usize) < n - k + i {
            set[i] += 1;
            for j in i + 1..k {
                set[j] = set[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

fn exact_key_players(p: &PathGraph, n: usize, k: usize, opts: KppOptions) -> Result<KeyPlayers> {
    let total = binomial(n, k);
    if total > EXACT_MAX_SETS {
        return Err(Error::InvalidArgument(format!(
            "exact search would evaluate {total} sets (more than {EXACT_MAX_SETS}); \
             use the greedy method"
        )));
    }
    let threads = parallel::resolve_threads(opts.threads, total as usize);
    let chunks = (threads * 8).min(total as usize).max(1);
    let per = total.div_euclid(chunks as u128);
    let extra = total % chunks as u128;
    let parts = parallel::for_each_item(
        chunks,
        threads,
        || (Scorer::new(p, opts), Vec::<(usize, Vec<u32>, f64)>::new()),
        |(scorer, found), c| {
            let c128 = c as u128;
            let begin = c128 * per + c128.min(extra);
            let len = per + u128::from(c128 < extra);
            if len == 0 {
                return;
            }
            let mut set = vec![0u32; k];
            unrank(begin, n, k, &mut set);
            let mut best: Option<(Vec<u32>, f64)> = None;
            for i in 0..len {
                let s = scorer.score(&set);
                if best.as_ref().map_or(true, |(_, b)| kpp_better(s, *b)) {
                    best = Some((set.clone(), s));
                }
                if i + 1 < len {
                    next_combination(&mut set, n);
                }
            }
            if let Some((set, s)) = best {
                found.push((c, set, s));
            }
        },
    );
    let mut all: Vec<(usize, Vec<u32>, f64)> = parts.into_iter().flat_map(|p| p.1).collect();
    all.sort_by_key(|x| x.0);
    // Chunks are in lexicographic order, so keeping the first of tied
    // scores keeps the lexicographically smallest set.
    let mut best: Option<(Vec<u32>, f64)> = None;
    for (_, set, s) in all {
        if best.as_ref().map_or(true, |(_, b)| kpp_better(s, *b)) {
            best = Some((set, s));
        }
    }
    let (set, score) = best.expect("at least one set");
    Ok(KeyPlayers {
        set,
        score,
        evaluations: total as u64,
    })
}

fn greedy_key_players(p: &PathGraph, n: usize, k: usize, opts: KppOptions) -> Result<KeyPlayers> {
    let restarts = opts.restarts.max(1);
    let parts = parallel::for_each_item(
        restarts,
        opts.threads,
        || {
            (
                Scorer::new(p, opts),
                Vec::<(usize, Vec<u32>, f64, u64)>::new(),
            )
        },
        |(scorer, found), r| {
            let mut rng = Stream::new(opts.seed, r as u64);
            let mut perm = vec![0usize; n];
            rng.permutation(&mut perm);
            let mut set: Vec<u32> = perm[..k].iter().map(|&v| v as u32).collect();
            set.sort_unstable();
            let mut member = vec![false; n];
            for &v in &set {
                member[v as usize] = true;
            }
            let mut fit = scorer.score(&set);
            let mut evaluations = 1u64;
            let mut trial = set.clone();
            loop {
                // Best single swap; ties keep the first found, scanning
                // members and then non-members in increasing order.
                let mut best: Option<(usize, u32, f64)> = None;
                for i in 0..k {
                    for v in 0..n as u32 {
                        if member[v as usize] {
                            continue;
                        }
                        trial.copy_from_slice(&set);
                        trial[i] = v;
                        trial.sort_unstable();
                        let s = scorer.score(&trial);
                        evaluations += 1;
                        let reference = best.map_or(fit, |b| b.2);
                        if kpp_better(s, reference) {
                            best = Some((i, v, s));
                        }
                    }
                }
                match best {
                    Some((i, v, s)) => {
                        member[set[i] as usize] = false;
                        member[v as usize] = true;
                        set[i] = v;
                        set.sort_unstable();
                        fit = s;
                    }
                    None => break,
                }
            }
            found.push((r, set, fit, evaluations));
        },
    );
    let mut all: Vec<(usize, Vec<u32>, f64, u64)> = parts.into_iter().flat_map(|p| p.1).collect();
    all.sort_by_key(|x| x.0);
    let evaluations = all.iter().map(|x| x.3).sum();
    let mut best: Option<(Vec<u32>, f64)> = None;
    for (_, set, s, _) in all {
        let replace = match &best {
            None => true,
            Some((bset, b)) => kpp_better(s, *b) || (!kpp_better(*b, s) && set < *bset),
        };
        if replace {
            best = Some((set, s));
        }
    }
    let (set, score) = best.expect("at least one restart");
    Ok(KeyPlayers {
        set,
        score,
        evaluations,
    })
}

// ---- Group centrality -----------------------------------------------------------

/// Which group centrality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupMeasure {
    /// Non-members adjacent to at least one member.
    Degree,
    /// Distances from the group to non-members.
    Closeness,
    /// Shortest paths between non-members passing through the group.
    Betweenness,
}

impl GroupMeasure {
    /// Decodes the C interface codes: 0 = degree, 1 = closeness,
    /// 2 = betweenness.
    pub fn from_code(code: i32) -> Result<GroupMeasure> {
        match code {
            0 => Ok(GroupMeasure::Degree),
            1 => Ok(GroupMeasure::Closeness),
            2 => Ok(GroupMeasure::Betweenness),
            _ => Err(Error::InvalidArgument(format!(
                "unknown group measure code {code}"
            ))),
        }
    }
}

/// Options for [`group_centrality`].
#[derive(Debug, Clone, Copy)]
pub struct GroupOptions {
    /// Degree: which neighbors count; closeness: direction of paths from
    /// the group (`Out`), to it (`In`), or ignoring direction (`All`).
    /// Ignored by betweenness and in undirected networks.
    pub mode: DegreeMode,
    /// Use edge weights as path lengths (closeness and betweenness).
    pub weighted: bool,
    /// Whether weights are distances or strengths.
    pub weights_as: WeightKind,
    /// Degree: divide by n - k. Closeness: (n - k) / sum instead of
    /// 1 / sum. Betweenness: divide by the number of pairs of non-members.
    pub normalized: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for GroupOptions {
    fn default() -> Self {
        GroupOptions {
            mode: DegreeMode::Out,
            weighted: false,
            weights_as: WeightKind::Distance,
            normalized: true,
            threads: 0,
        }
    }
}

/// Group degree, closeness, or betweenness of `group`. See the module
/// documentation.
///
/// Closeness ignores non-members the group cannot reach, as NetworkX does,
/// and is 0 when it reaches none.
pub fn group_centrality(
    g: &Graph,
    group: &[u32],
    measure: GroupMeasure,
    opts: GroupOptions,
) -> Result<f64> {
    let n = g.node_count();
    check_nodes(n, group, "group")?;
    if group.is_empty() {
        return Err(Error::InvalidArgument("the group has no members".into()));
    }
    let k = group.len();
    let mut member = vec![false; n];
    for &v in group {
        member[v as usize] = true;
    }
    let rest = (n - k) as f64;
    let mode = if g.is_directed() {
        opts.mode
    } else {
        DegreeMode::Out
    };
    match measure {
        GroupMeasure::Degree => {
            let mut hit = vec![false; n];
            for &v in group {
                let v = v as usize;
                if mode != DegreeMode::In {
                    for &w in g.out().neighbors(v) {
                        hit[w as usize] = true;
                    }
                }
                if g.is_directed() && mode != DegreeMode::Out {
                    for &w in g.inc().neighbors(v) {
                        hit[w as usize] = true;
                    }
                }
            }
            let count = (0..n).filter(|&w| hit[w] && !member[w]).count() as f64;
            Ok(if opts.normalized { count / rest } else { count })
        }
        GroupMeasure::Closeness => {
            let p = path_graph(g, mode, opts.weighted, opts.weights_as)?;
            let mut search = Search::new(n);
            search.run(&p, group, &vec![false; n]);
            let total: f64 = search.order[k..]
                .iter()
                .map(|&v| search.dist[v as usize] * p.unit)
                .sum();
            Ok(if total == 0.0 {
                0.0
            } else if opts.normalized {
                rest / total
            } else {
                1.0 / total
            })
        }
        GroupMeasure::Betweenness => {
            let raw = group_betweenness(g, &member, opts)?;
            let nk = n - k;
            if !opts.normalized {
                return Ok(raw);
            }
            if nk < 2 {
                return Ok(f64::NAN);
            }
            let pairs = nk as f64 * (nk as f64 - 1.0);
            Ok(if g.is_directed() {
                raw / pairs
            } else {
                2.0 * raw / pairs
            })
        }
    }
}

/// Sum over pairs of non-members of the fraction of shortest paths through
/// at least one member. For each non-member source, sigma counts all
/// shortest paths and `avoid` those using no member; the difference is the
/// number through the group. Pairs are ordered in directed networks and
/// unordered in undirected ones.
fn group_betweenness(g: &Graph, member: &[bool], opts: GroupOptions) -> Result<f64> {
    let n = g.node_count();
    let lengths = path_lengths(g, opts.weighted, opts.weights_as, PathAlgorithm::Auto)?;
    let p = PathGraph {
        g: g.clone(),
        lengths,
        unit: 1.0,
    };
    let sources: Vec<u32> = (0..n as u32).filter(|&v| !member[v as usize]).collect();
    let parts = parallel::for_each_item(
        sources.len(),
        opts.threads,
        || (Search::new(n), vec![0.0f64; n], vec![0.0f64; n], 0.0f64),
        |(search, sigma, avoid, total), i| {
            let s = sources[i] as usize;
            search.run(&p, &[s as u32], &vec![false; n]);
            sigma[s] = 1.0;
            avoid[s] = 1.0;
            let adj = p.g.out();
            for idx in 0..search.order.len() {
                let v = search.order[idx] as usize;
                let dv = search.dist[v];
                for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                    let w = w as usize;
                    if w == v {
                        continue;
                    }
                    let len = match &p.lengths {
                        Lengths::Uniform => 1.0,
                        Lengths::Integer(l, _) => l[k] as f64,
                        Lengths::Real(l) => l[k],
                    };
                    let on_path = match &p.lengths {
                        Lengths::Real(_) => {
                            search.settled[w]
                                && search.dist[w] > dv
                                && approx_eq(dv + len, search.dist[w])
                        }
                        _ => dv + len == search.dist[w],
                    };
                    if on_path {
                        sigma[w] += sigma[v];
                        if !member[v] {
                            avoid[w] += avoid[v];
                        }
                    }
                }
            }
            for &t in &search.order {
                let t = t as usize;
                if t != s && !member[t] {
                    *total += 1.0 - avoid[t] / sigma[t];
                }
            }
            for &v in &search.order {
                sigma[v as usize] = 0.0;
                avoid[v as usize] = 0.0;
            }
        },
    );
    let total: f64 = parts.iter().map(|p| p.3).sum();
    Ok(if g.is_directed() { total } else { total / 2.0 })
}

// ---- Flow betweenness ------------------------------------------------------------------

/// A residual network for Dinic's algorithm. Arc 2e is edge e in its
/// stated direction and arc 2e + 1 its reverse; an undirected edge gives
/// both arcs its capacity.
struct FlowNet {
    nodes: usize,
    head: Vec<u32>,
    tail: Vec<u32>,
    cap: Vec<f64>,
    residual: Vec<f64>,
    offsets: Vec<usize>,
    arcs: Vec<u32>,
    level: Vec<u32>,
    next: Vec<usize>,
    queue: Vec<u32>,
    eps: f64,
}

impl FlowNet {
    fn new(g: &Graph, capacity: Option<&[f64]>) -> FlowNet {
        let n = g.node_count();
        let m = g.edge_count();
        let mut head = Vec::with_capacity(2 * m);
        let mut tail = Vec::with_capacity(head.capacity());
        let mut cap = Vec::with_capacity(head.capacity());
        let mut max_cap = 0.0f64;
        for e in 0..m {
            let (a, b) = (g.edge_sources()[e], g.edge_targets()[e]);
            let c = if a == b {
                0.0
            } else {
                capacity.map_or(1.0, |w| w[e])
            };
            max_cap = max_cap.max(c);
            tail.push(a);
            head.push(b);
            cap.push(c);
            tail.push(b);
            head.push(a);
            cap.push(if g.is_directed() { 0.0 } else { c });
        }
        let mut offsets = vec![0usize; n + 1];
        for &t in &tail {
            offsets[t as usize + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let mut cursor = offsets.clone();
        let mut arcs = vec![0u32; tail.len()];
        for (k, &t) in tail.iter().enumerate() {
            arcs[cursor[t as usize]] = k as u32;
            cursor[t as usize] += 1;
        }
        FlowNet {
            nodes: n,
            residual: cap.clone(),
            head,
            tail,
            cap,
            offsets,
            arcs,
            level: vec![u32::MAX; n],
            next: vec![0; n],
            queue: Vec::with_capacity(n),
            eps: 1e-12 * max_cap.max(1.0),
        }
    }

    fn reset(&mut self) {
        self.residual.copy_from_slice(&self.cap);
    }

    /// Breadth-first levels from `s` over arcs with residual capacity.
    fn levels(&mut self, s: usize, t: usize, blocked: Option<usize>) -> bool {
        self.level.iter_mut().for_each(|l| *l = u32::MAX);
        self.queue.clear();
        self.level[s] = 0;
        self.queue.push(s as u32);
        let mut i = 0;
        while i < self.queue.len() {
            let v = self.queue[i] as usize;
            i += 1;
            for &a in &self.arcs[self.offsets[v]..self.offsets[v + 1]] {
                let a = a as usize;
                let w = self.head[a] as usize;
                if self.residual[a] > self.eps && self.level[w] == u32::MAX && Some(w) != blocked {
                    self.level[w] = self.level[v] + 1;
                    self.queue.push(w as u32);
                }
            }
        }
        self.level[t] != u32::MAX
    }

    /// Maximum flow from `s` to `t`, optionally with node `blocked` removed.
    /// Starts from the current residual capacities.
    fn run(&mut self, s: usize, t: usize, blocked: Option<usize>) -> f64 {
        let mut flow = 0.0;
        let mut path: Vec<u32> = Vec::new();
        while self.levels(s, t, blocked) {
            for v in 0..self.nodes {
                self.next[v] = self.offsets[v];
            }
            // Iterative depth-first search for blocking flow, advancing
            // through each node's current arc.
            path.clear();
            let mut v = s;
            loop {
                if v == t {
                    let mut bottleneck = f64::INFINITY;
                    for &a in &path {
                        bottleneck = bottleneck.min(self.residual[a as usize]);
                    }
                    for &a in &path {
                        let a = a as usize;
                        self.residual[a] -= bottleneck;
                        self.residual[a ^ 1] += bottleneck;
                    }
                    flow += bottleneck;
                    // Retreat to the tail of the first saturated arc.
                    let cut = path
                        .iter()
                        .position(|&a| self.residual[a as usize] <= self.eps)
                        .unwrap_or(0);
                    path.truncate(cut);
                    v = path.last().map_or(s, |&a| self.head[a as usize] as usize);
                    continue;
                }
                let mut advanced = false;
                while self.next[v] < self.offsets[v + 1] {
                    let a = self.arcs[self.next[v]] as usize;
                    let w = self.head[a] as usize;
                    if self.residual[a] > self.eps && self.level[w] == self.level[v] + 1 {
                        path.push(a as u32);
                        v = w;
                        advanced = true;
                        break;
                    }
                    self.next[v] += 1;
                }
                if advanced {
                    continue;
                }
                // Dead end: no augmenting path continues through v.
                self.level[v] = u32::MAX;
                match path.pop() {
                    None => break,
                    Some(a) => {
                        v = self.tail[a as usize] as usize;
                        self.next[v] += 1;
                    }
                }
            }
        }
        flow
    }
}

fn flow_capacities(g: &Graph, use_capacity: bool) -> Result<Option<&[f64]>> {
    additive_weights(g, use_capacity, "a flow capacity")
}

/// How flow betweenness is scaled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FlowScale {
    /// Freeman, Borgatti and White's total.
    #[default]
    Raw,
    /// Divided by the total maximum flow between pairs not involving the
    /// node (UCINET; sna `normflow`).
    Normalized,
    /// Sum over pairs of the fraction of their maximum flow that must pass
    /// through the node (sna `fracflow`).
    Fractional,
}

impl FlowScale {
    /// Decodes the C interface codes: 0 = raw, 1 = normalized,
    /// 2 = fractional.
    pub fn from_code(code: i32) -> Result<FlowScale> {
        match code {
            0 => Ok(FlowScale::Raw),
            1 => Ok(FlowScale::Normalized),
            2 => Ok(FlowScale::Fractional),
            _ => Err(Error::InvalidArgument(format!(
                "unknown flow betweenness scale code {code}"
            ))),
        }
    }
}

/// Flow betweenness of every node. See the module documentation.
pub fn flow_betweenness(
    g: &Graph,
    use_capacity: bool,
    scale: FlowScale,
    threads: usize,
) -> Result<Vec<f64>> {
    let n = g.node_count();
    let cap = flow_capacities(g, use_capacity)?;
    let directed = g.is_directed();
    let base = FlowNet::new(g, cap);
    // (raw, fractional, flow on pairs involving the node), and the total
    // flow over all pairs.
    let parts = parallel::for_each_item(
        n,
        threads,
        || {
            (
                base_clone(&base),
                vec![0.0f64; n],
                vec![0.0f64; n],
                vec![0.0f64; n],
                0.0f64,
                vec![0.0f64; n],
            )
        },
        |(net, raw, frac, involved, total, through), i| {
            let targets: Box<dyn Iterator<Item = usize>> = if directed {
                Box::new((0..n).filter(move |&j| j != i))
            } else {
                Box::new(i + 1..n)
            };
            for j in targets {
                net.reset();
                let f = net.run(i, j, None);
                if f <= net.eps {
                    continue;
                }
                *total += f;
                involved[i] += f;
                involved[j] += f;
                // Nodes carrying no flow in this maximum flow lose nothing
                // when removed, since the same flow avoids them.
                through.iter_mut().for_each(|x| *x = 0.0);
                for (a, (&h, &c)) in net.head.iter().zip(&net.cap).enumerate() {
                    let used = c - net.residual[a];
                    if used > net.eps {
                        through[h as usize] += used;
                    }
                }
                let carrying: Vec<usize> = (0..n)
                    .filter(|&v| v != i && v != j && through[v] > net.eps)
                    .collect();
                for v in carrying {
                    net.reset();
                    let fv = net.run(i, j, Some(v));
                    let lost = (f - fv).max(0.0);
                    if lost > net.eps * f.max(1.0) {
                        raw[v] += lost;
                        frac[v] += lost / f;
                    }
                }
            }
        },
    );
    let mut raw = vec![0.0; n];
    let mut frac = vec![0.0; n];
    let mut involved = vec![0.0; n];
    let mut total = 0.0;
    for p in parts {
        for v in 0..n {
            raw[v] += p.1[v];
            frac[v] += p.2[v];
            involved[v] += p.3[v];
        }
        total += p.4;
    }
    Ok(match scale {
        FlowScale::Raw => raw,
        FlowScale::Fractional => frac,
        FlowScale::Normalized => (0..n)
            .map(|v| {
                let denom = total - involved[v];
                if denom > 0.0 {
                    raw[v] / denom
                } else {
                    f64::NAN
                }
            })
            .collect(),
    })
}

fn base_clone(b: &FlowNet) -> FlowNet {
    FlowNet {
        nodes: b.nodes,
        head: b.head.clone(),
        tail: b.tail.clone(),
        cap: b.cap.clone(),
        residual: b.residual.clone(),
        offsets: b.offsets.clone(),
        arcs: b.arcs.clone(),
        level: b.level.clone(),
        next: b.next.clone(),
        queue: Vec::with_capacity(b.nodes),
        eps: b.eps,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(n: usize, s: &[i64], t: &[i64], directed: bool) -> Graph {
        Graph::from_edges(n, s, t, None, directed).unwrap()
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn kpp_scores_by_hand() {
        // Path 0-1-2-3-4.
        let p = g(5, &[0, 1, 2, 3], &[1, 2, 3, 4], false);
        let opts = KppOptions::default();
        // Removing 2 leaves two pairs {0,1}, {3,4}: 4 of 12 ordered pairs
        // connected.
        assert!(close(kpp_score(&p, &[2], opts).unwrap(), 8.0 / 12.0));
        let df = KppOptions {
            measure: KppMeasure::DistanceFragmentation,
            ..opts
        };
        // Same: each connected pair is at distance 1.
        assert!(close(kpp_score(&p, &[2], df).unwrap(), 8.0 / 12.0));
        // Removing 0: path of 4; sum of 1/d over ordered pairs is
        // 2 (3 + 2/2 + 1/3) = 26/3, over 12 pairs.
        assert!(close(kpp_score(&p, &[0], df).unwrap(), 1.0 - 26.0 / 36.0));
        let dr = KppOptions {
            measure: KppMeasure::DistanceReach,
            ..opts
        };
        // From 2: itself 1, two at distance 1, two at 2: (1 + 2 + 1) / 5.
        assert!(close(kpp_score(&p, &[2], dr).unwrap(), 4.0 / 5.0));
        let mr = KppOptions {
            measure: KppMeasure::MReach,
            m: 1.0,
            ..opts
        };
        assert_eq!(kpp_score(&p, &[0, 4], mr).unwrap(), 2.0);
        let mr2 = KppOptions { m: 2.0, ..mr };
        assert_eq!(kpp_score(&p, &[0], mr2).unwrap(), 2.0);
        assert!(kpp_score(&p, &[0, 0], opts).is_err());
        assert!(kpp_score(&g(2, &[0], &[1], false), &[0], opts)
            .unwrap()
            .is_nan());
    }

    #[test]
    fn directed_fragmentation_counts_reachable_ordered_pairs() {
        // 0 -> 1 -> 2: 3 of 6 ordered pairs reachable.
        let d = g(3, &[0, 1], &[1, 2], true);
        let opts = KppOptions::default();
        let empty: [u32; 0] = [];
        assert!(close(kpp_score(&d, &empty, opts).unwrap(), 0.5));
        let all = KppOptions {
            mode: DegreeMode::All,
            ..opts
        };
        assert!(close(kpp_score(&d, &empty, all).unwrap(), 0.0));
    }

    #[test]
    fn greedy_and_exact_agree_on_a_barbell() {
        // Two triangles joined through node 3: 0-1-2 and 4-5-6.
        let b = g(
            7,
            &[0, 1, 2, 2, 3, 4, 5, 6],
            &[1, 2, 0, 3, 4, 5, 6, 4],
            false,
        );
        // Removing 2 or 3 leaves DF = 0.6 either way, so the tie goes to the
        // smaller index. F singles out the cut point 3; DR prefers 2, which
        // reaches 0, 1, and 3 in one step (31/42 against 5/7).
        for (measure, best) in [
            (KppMeasure::Fragmentation, 3),
            (KppMeasure::DistanceFragmentation, 2),
            (KppMeasure::DistanceReach, 2),
        ] {
            let opts = KppOptions {
                measure,
                restarts: 5,
                threads: 2,
                ..Default::default()
            };
            let exact = key_players(
                &b,
                1,
                KppOptions {
                    method: KppMethod::Exact,
                    ..opts
                },
            )
            .unwrap();
            let greedy = key_players(&b, 1, opts).unwrap();
            assert_eq!(exact.set, vec![best], "{measure:?}");
            assert_eq!(greedy.set, exact.set);
            assert!(close(greedy.score, exact.score));
        }
    }

    #[test]
    fn greedy_is_reproducible_across_threads() {
        let s: Vec<i64> = (0..30).map(|i| i % 12).collect();
        let t: Vec<i64> = (0..30).map(|i| (i * 7 + 3) % 12).collect();
        let h = g(12, &s, &t, false);
        let opts = KppOptions {
            measure: KppMeasure::DistanceFragmentation,
            restarts: 6,
            seed: 9,
            ..Default::default()
        };
        let a = key_players(&h, 3, KppOptions { threads: 1, ..opts }).unwrap();
        let b = key_players(&h, 3, KppOptions { threads: 3, ..opts }).unwrap();
        assert_eq!(a, b);
        let e = key_players(
            &h,
            3,
            KppOptions {
                method: KppMethod::Exact,
                ..opts
            },
        )
        .unwrap();
        // The exact optimum is never worse than the greedy result.
        assert!(e.score >= a.score - 1e-12);
    }

    #[test]
    fn combinations_unrank_in_order() {
        let mut set = vec![0u32; 3];
        let mut walk = vec![0, 1, 2];
        for r in 0..binomial(6, 3) {
            unrank(r, 6, 3, &mut set);
            assert_eq!(set, walk);
            next_combination(&mut walk, 6);
        }
    }

    #[test]
    fn group_centrality_by_hand() {
        // Path 0-1-2-3-4, group {1, 3}.
        let p = g(5, &[0, 1, 2, 3], &[1, 2, 3, 4], false);
        let opts = GroupOptions::default();
        let deg = group_centrality(&p, &[1, 3], GroupMeasure::Degree, opts).unwrap();
        assert!(close(deg, 1.0));
        let clo = group_centrality(&p, &[1, 3], GroupMeasure::Closeness, opts).unwrap();
        assert!(close(clo, 1.0));
        // Pairs of non-members {0,2}, {0,4}, {2,4}: every path crosses the
        // group.
        let bet = group_centrality(&p, &[1, 3], GroupMeasure::Betweenness, opts).unwrap();
        assert!(close(bet, 1.0));
        let raw = group_centrality(
            &p,
            &[1, 3],
            GroupMeasure::Betweenness,
            GroupOptions {
                normalized: false,
                ..opts
            },
        )
        .unwrap();
        assert!(close(raw, 3.0));
        // A 4-cycle 0-1-2-3: two shortest paths 0 -> 2, one through 1.
        let c4 = g(4, &[0, 1, 2, 3], &[1, 2, 3, 0], false);
        let raw = group_centrality(
            &c4,
            &[1],
            GroupMeasure::Betweenness,
            GroupOptions {
                normalized: false,
                ..opts
            },
        )
        .unwrap();
        assert!(close(raw, 0.5));
    }

    #[test]
    fn flow_betweenness_by_hand() {
        // Path 0-1-2: all flow between 0 and 2 passes through 1.
        let p = g(3, &[0, 1], &[1, 2], false);
        let raw = flow_betweenness(&p, false, FlowScale::Raw, 1).unwrap();
        assert_eq!(raw, vec![0.0, 1.0, 0.0]);
        let norm = flow_betweenness(&p, false, FlowScale::Normalized, 1).unwrap();
        assert!(close(norm[1], 1.0));
        // A 4-cycle: flow 2 between opposite nodes; removing a middle node
        // leaves 1, so each loses 1 of 2.
        let c4 = g(4, &[0, 1, 2, 3], &[1, 2, 3, 0], false);
        let raw = flow_betweenness(&c4, false, FlowScale::Raw, 2).unwrap();
        // Node 1 lies between 0 and 2 (loses 1), and between pairs of
        // adjacent nodes {0, 3}? flow 2 there, 1 without node 1: loses 1;
        // {2, 3}: likewise 1. Total 3.
        assert!(raw.iter().all(|&x| close(x, 3.0)), "{raw:?}");
        let frac = flow_betweenness(&c4, false, FlowScale::Fractional, 1).unwrap();
        assert!(frac.iter().all(|&x| close(x, 1.5)));
    }
}
