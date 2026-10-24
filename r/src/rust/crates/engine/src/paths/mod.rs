//! Shortest paths: the paths themselves, not only their lengths.
//!
//! The measures in [`crate::metrics`] need path lengths and path counts;
//! this module returns the paths, predecessor structures, and the
//! algorithms for lengths that the rest of the engine does not allow
//! (zero and negative lengths, arbitrary lengths supplied per edge, a
//! heuristic for goal-directed search). Everything works on a
//! [`PathGraph`]: the arcs a search may follow, oriented for a mode, with
//! one length per arc. A caller builds it once and runs any number of
//! searches on it.
//!
//! | Problem | Function | Algorithm |
//! |---|---|---|
//! | one source, or several at once | [`search::single_source`], [`search::multi_source`] | breadth-first search for unit lengths; otherwise Dijkstra's algorithm with an indexed 4-ary heap (`dijkstra1959note`, `johnson1975priority`) |
//! | every source | [`search::all_pairs`] | one search per source, in parallel |
//! | all predecessors on shortest paths | [`search::predecessors`] | one scan of the arcs after a search |
//! | one pair | [`bidirectional::bidirectional`] | bidirectional breadth-first search or Dijkstra (`nicholson1966finding`) |
//! | one pair with a heuristic | [`astar::astar`] | A* (`hart1968formal`) |
//! | negative lengths | [`negative::bellman_ford`], [`negative::goldberg_radzik`] | Bellman-Ford-Moore with a FIFO queue (`bellman1958routing`, `cherkassky1996shortest`); Goldberg and Radzik's topological-scan algorithm (`goldberg1993heuristic`) |
//! | negative cycles | [`negative::find_negative_cycle`], [`negative::has_negative_cycle`] | the same searches |
//! | every pair, negative lengths | [`negative::johnson`] | Johnson's reweighting (`johnson1977efficient`) |
//! | every pair, dense | [`dense::floyd_warshall`] | Floyd-Warshall (`floyd1962algorithm`, `warshall1962theorem`) |
//! | every shortest path of a pair | [`enumerate::all_shortest_paths`] | enumeration over the shortest-path DAG |
//! | k shortest simple paths | [`simple::k_shortest_simple_paths`] | Yen's algorithm (`yen1971finding`) with Lawler's restriction of spur nodes (`lawler1972procedure`) |
//!
//! Conventions shared by every function here:
//! - Nodes are `0..n`. Paths are node sequences that start at the source.
//! - Direction: a [`PathGraph`] built with `Out` follows arcs forward, with
//!   `In` backward (paths *to* the source in the original direction), and
//!   with `All` in both directions. Undirected networks ignore the mode.
//! - Lengths: unit lengths (unweighted), the edge weights read as distances
//!   or as strengths (length 1/w), or one length per edge supplied by the
//!   caller (for example -log p for edges that work with probability p, so
//!   that a shortest path is a most likely path). An infinite length hides
//!   the edge, the way a NetworkX weight function hides an edge by
//!   returning `None`. Lengths may be zero; searches that need
//!   non-negative lengths say so and refuse negative ones, and the
//!   negative-length algorithms accept them.
//! - Parallel edges: a path uses the shortest of them; with
//!   [`PathGraphOptions::collapse_parallel`] only that one is kept, which
//!   is how NetworkX reads a multigraph. Self-loops never shorten a path
//!   unless their length is negative, which is a negative cycle.
//! - Ties: two path lengths are equal when they agree to a relative
//!   tolerance of 1e-10 ([`Tie::Relative`], the engine's default and
//!   igraph's policy), or only when they are equal as floating-point
//!   numbers ([`Tie::Exact`], NetworkX). Integer lengths are summed
//!   exactly either way. Which of several shortest paths a search returns
//!   is deterministic: among nodes at the same distance, the one with the
//!   smaller index is settled first, and a node keeps the first
//!   predecessor that reached its final distance. NetworkX settles ties
//!   in the order nodes were inserted into its adjacency dictionaries, so
//!   on networks with tied paths it can return another, equally short,
//!   path.
//! - Searches are designed as building blocks for further path problems:
//!   the label-setting core takes the arrival label of an arc as a
//!   function ([`search::Searcher::run_with`]), so a time-dependent arrival
//!   function with the FIFO property (`dreyfus1969appraisal`) uses the same
//!   code, and every search can start from many sources at once and
//!   report which source each node is closest to.
//!
//! NetworkX (`hagberg2008exploring`) documents the semantics these
//! functions reproduce; its source was read to pin down edge cases, and
//! each module says where it differs.

pub mod astar;
pub mod bidirectional;
pub mod dense;
pub mod enumerate;
pub mod negative;
pub mod search;
pub mod simple;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;
use crate::weights::WeightKind;

/// Marker for "no node" in parent and predecessor arrays.
pub const NO_NODE: u32 = u32::MAX;

/// Relative tolerance of [`Tie::Relative`].
pub const TIE_TOLERANCE: f64 = 1e-10;

/// When two path lengths count as equal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Tie {
    /// Equal to a relative tolerance of [`TIE_TOLERANCE`] (igraph, and the
    /// rest of the engine).
    #[default]
    Relative,
    /// Equal as floating-point numbers (NetworkX).
    Exact,
}

impl Tie {
    /// Decodes the C interface codes: 0 = relative, 1 = exact.
    pub fn from_code(code: i32) -> Result<Tie> {
        match code {
            0 => Ok(Tie::Relative),
            1 => Ok(Tie::Exact),
            _ => Err(Error::InvalidArgument(format!("unknown tie code {code}"))),
        }
    }

    /// Whether `a` and `b` are equal under this rule.
    #[inline]
    pub fn equal(self, a: f64, b: f64) -> bool {
        a == b
            || (self == Tie::Relative
                && a.is_finite()
                && b.is_finite()
                && (a - b).abs() <= TIE_TOLERANCE * a.abs().max(b.abs()))
    }
}

/// Where arc lengths come from.
#[derive(Debug, Clone, Copy)]
pub enum ArcLengths<'a> {
    /// Every arc has length 1 (weights are ignored).
    Unit,
    /// The graph's edge weights, as distances or as strengths (length 1/w,
    /// which requires positive weights).
    Weights(WeightKind),
    /// One length per edge, in input order. Infinite lengths hide the edge.
    PerEdge(&'a [f64]),
}

/// How to build a [`PathGraph`].
#[derive(Debug, Clone, Copy)]
pub struct PathGraphOptions<'a> {
    /// Direction of the paths in directed graphs.
    pub mode: DegreeMode,
    /// Arc lengths.
    pub lengths: ArcLengths<'a>,
    /// Keep only the shortest of parallel arcs (the first such edge in input
    /// order), as NetworkX reads a multigraph. Counts of shortest paths then
    /// count node sequences rather than edge sequences.
    pub collapse_parallel: bool,
}

impl Default for PathGraphOptions<'_> {
    fn default() -> Self {
        PathGraphOptions {
            mode: DegreeMode::Out,
            lengths: ArcLengths::Unit,
            collapse_parallel: false,
        }
    }
}

/// The arcs a path may follow, with their lengths, in compressed sparse
/// row form: the arcs leaving `v` are `out_arcs(v)`, sorted by head and
/// then by edge index, and the arcs entering `v` are listed by
/// `in_arcs(v)`. Each arc remembers the edge it comes from.
#[derive(Debug, Clone)]
pub struct PathGraph {
    n: usize,
    offsets: Vec<usize>,
    heads: Vec<u32>,
    tails: Vec<u32>,
    lens: Vec<f64>,
    edges: Vec<u32>,
    rev_offsets: Vec<usize>,
    rev_arcs: Vec<u32>,
    unit: bool,
    min_len: f64,
    integral: bool,
    edge_count: usize,
}

impl PathGraph {
    /// Builds the arcs of `g` for `opts`. See the module documentation.
    pub fn new(g: &Graph, opts: PathGraphOptions) -> Result<PathGraph> {
        let n = g.node_count();
        let m = g.edge_count();
        let (src, tgt) = (g.edge_sources(), g.edge_targets());
        let edge_len: Vec<f64> = match opts.lengths {
            ArcLengths::Unit => vec![1.0; m],
            ArcLengths::Weights(kind) => {
                let w = g.edge_weights().ok_or(Error::MissingWeights)?;
                match kind {
                    WeightKind::Distance => w.to_vec(),
                    WeightKind::Strength => {
                        if let Some(edge) = w.iter().position(|&x| x <= 0.0) {
                            return Err(Error::InvalidArgument(format!(
                                "edge {edge} has weight {}; weights read as strengths become \
                                 lengths 1/w and must be strictly positive",
                                w[edge]
                            )));
                        }
                        w.iter().map(|x| 1.0 / x).collect()
                    }
                }
            }
            ArcLengths::PerEdge(l) => {
                if l.len() != m {
                    return Err(Error::LengthMismatch {
                        what: "edge lengths",
                        expected: m,
                        found: l.len(),
                    });
                }
                if let Some(edge) = l.iter().position(|x| x.is_nan() || *x == f64::NEG_INFINITY) {
                    return Err(Error::InvalidArgument(format!(
                        "edge {edge} has length {}; lengths must be numbers, and +Inf hides an edge",
                        l[edge]
                    )));
                }
                l.to_vec()
            }
        };
        let forward = !g.is_directed() || opts.mode != DegreeMode::In;
        let backward = !g.is_directed() || opts.mode != DegreeMode::Out;

        // Arcs in edge order (forward before backward for each edge), then
        // two stable counting sorts: by head, then by tail.
        let mut arcs: Vec<(u32, u32, u32)> = Vec::with_capacity(2 * m);
        for e in 0..m {
            if edge_len[e] == f64::INFINITY {
                continue;
            }
            let (s, t) = (src[e], tgt[e]);
            if forward {
                arcs.push((s, t, e as u32));
            }
            if backward && (s != t || !forward) {
                arcs.push((t, s, e as u32));
            }
        }
        let arcs = counting_sort(counting_sort(arcs, n, |a| a.1), n, |a| a.0);

        let mut heads = Vec::with_capacity(arcs.len());
        let mut tails = Vec::with_capacity(arcs.len());
        let mut lens = Vec::with_capacity(arcs.len());
        let mut edges = Vec::with_capacity(arcs.len());
        for &(t, h, e) in &arcs {
            let len = edge_len[e as usize];
            if opts.collapse_parallel {
                let k = heads.len();
                if k > 0 && tails[k - 1] == t && heads[k - 1] == h {
                    if len < lens[k - 1] {
                        lens[k - 1] = len;
                        edges[k - 1] = e;
                    }
                    continue;
                }
            }
            heads.push(h);
            tails.push(t);
            lens.push(len);
            edges.push(e);
        }
        let mut offsets = vec![0usize; n + 1];
        for &t in &tails {
            offsets[t as usize + 1] += 1;
        }
        for v in 0..n {
            offsets[v + 1] += offsets[v];
        }
        // Arcs entering each node, ordered by tail (arcs are already sorted
        // by tail, so a stable counting sort by head keeps that order).
        let mut rev_offsets = vec![0usize; n + 1];
        for &h in &heads {
            rev_offsets[h as usize + 1] += 1;
        }
        for v in 0..n {
            rev_offsets[v + 1] += rev_offsets[v];
        }
        let mut cursor = rev_offsets.clone();
        let mut rev_arcs = vec![0u32; heads.len()];
        for (a, &h) in heads.iter().enumerate() {
            rev_arcs[cursor[h as usize]] = a as u32;
            cursor[h as usize] += 1;
        }
        let unit = lens.iter().all(|&x| x == 1.0);
        let min_len = lens.iter().copied().fold(f64::INFINITY, f64::min);
        let integral = lens
            .iter()
            .all(|&x| x.fract() == 0.0 && x.abs() < 2f64.powi(52));
        Ok(PathGraph {
            n,
            offsets,
            heads,
            tails,
            lens,
            edges,
            rev_offsets,
            rev_arcs,
            unit,
            min_len,
            integral,
            edge_count: m,
        })
    }

    /// The same arcs with other lengths (one per arc, in arc order).
    pub(crate) fn with_lengths(&self, lens: Vec<f64>) -> PathGraph {
        debug_assert_eq!(lens.len(), self.arc_count());
        let mut pg = self.clone();
        pg.unit = lens.iter().all(|&x| x == 1.0);
        pg.min_len = lens.iter().copied().fold(f64::INFINITY, f64::min);
        pg.integral = lens
            .iter()
            .all(|&x| x.fract() == 0.0 && x.abs() < 2f64.powi(52));
        pg.lens = lens;
        pg
    }

    /// Number of nodes.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.n
    }

    /// Number of edges of the graph this was built from.
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.edge_count
    }

    /// Number of arcs.
    #[inline]
    pub fn arc_count(&self) -> usize {
        self.heads.len()
    }

    /// The arcs leaving `v`, as a range of arc indices.
    #[inline]
    pub fn out_arcs(&self, v: usize) -> std::ops::Range<usize> {
        self.offsets[v]..self.offsets[v + 1]
    }

    /// The arcs entering `v`, as arc indices ordered by tail.
    #[inline]
    pub fn in_arcs(&self, v: usize) -> &[u32] {
        &self.rev_arcs[self.rev_offsets[v]..self.rev_offsets[v + 1]]
    }

    /// Head of arc `a`.
    #[inline]
    pub fn head(&self, a: usize) -> u32 {
        self.heads[a]
    }

    /// Tail of arc `a`.
    #[inline]
    pub fn tail(&self, a: usize) -> u32 {
        self.tails[a]
    }

    /// Length of arc `a`.
    #[inline]
    pub fn len(&self, a: usize) -> f64 {
        self.lens[a]
    }

    /// The edge (input index) arc `a` comes from.
    #[inline]
    pub fn edge(&self, a: usize) -> u32 {
        self.edges[a]
    }

    /// Whether the graph has no arcs.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.heads.is_empty()
    }

    /// Whether every arc has length 1, so breadth-first search applies.
    #[inline]
    pub fn is_unit(&self) -> bool {
        self.unit
    }

    /// Whether every length is an integer small enough that sums of up to
    /// n of them are exact in floating point.
    #[inline]
    pub fn is_integral(&self) -> bool {
        self.integral
    }

    /// The smallest arc length (infinity without arcs).
    #[inline]
    pub fn min_len(&self) -> f64 {
        self.min_len
    }

    /// An error unless every length is non-negative, naming `what` needs it.
    pub fn require_nonnegative(&self, what: &str) -> Result<()> {
        if self.min_len < 0.0 {
            let a = self.lens.iter().position(|&x| x < 0.0).unwrap_or(0);
            return Err(Error::InvalidArgument(format!(
                "edge {} has negative length {}; {what} needs non-negative lengths \
                 (use Bellman-Ford, Goldberg-Radzik, or Johnson's algorithm)",
                self.edges[a], self.lens[a]
            )));
        }
        Ok(())
    }

    /// An error unless every length is strictly positive.
    pub fn require_positive(&self, what: &str) -> Result<()> {
        if self.min_len <= 0.0 {
            let a = self.lens.iter().position(|&x| x <= 0.0).unwrap_or(0);
            return Err(Error::InvalidArgument(format!(
                "edge {} has length {}; {what} needs strictly positive lengths",
                self.edges[a], self.lens[a]
            )));
        }
        Ok(())
    }

    /// An error unless `v` is a node.
    pub fn check_node(&self, v: usize, what: &str) -> Result<()> {
        if v >= self.n {
            return Err(Error::InvalidArgument(format!(
                "{what} node index {v} is out of range for a network with {} nodes",
                self.n
            )));
        }
        Ok(())
    }

    /// The first arc from `u` to `v` of smallest length, if any.
    pub fn arc_between(&self, u: usize, v: usize) -> Option<usize> {
        let r = self.out_arcs(u);
        let row = &self.heads[r.clone()];
        let start = row.partition_point(|&h| (h as usize) < v);
        let mut best: Option<usize> = None;
        for (i, &h) in row[start..].iter().enumerate() {
            if h as usize != v {
                break;
            }
            let a = r.start + start + i;
            if best.map_or(true, |b| self.lens[a] < self.lens[b]) {
                best = Some(a);
            }
        }
        best
    }

    /// The length of a path given as a node sequence, using the shortest
    /// arc between consecutive nodes, or `None` if some step is not an arc.
    pub fn path_length(&self, path: &[u32]) -> Option<f64> {
        let mut total = 0.0;
        for p in path.windows(2) {
            total += self.lens[self.arc_between(p[0] as usize, p[1] as usize)?];
        }
        Some(total)
    }
}

/// A stable counting sort of `items` by `key(item)` in `0..n`.
fn counting_sort<T: Copy>(items: Vec<T>, n: usize, key: impl Fn(&T) -> u32) -> Vec<T> {
    if items.is_empty() {
        return items;
    }
    let mut start = vec![0usize; n + 1];
    for it in &items {
        start[key(it) as usize + 1] += 1;
    }
    for i in 0..n {
        start[i + 1] += start[i];
    }
    let mut out = items.clone();
    for it in items {
        let k = key(&it) as usize;
        out[start[k]] = it;
        start[k] += 1;
    }
    out
}

/// A collection of node lists (paths, predecessor lists, cycles), stored
/// flat: list `i` is `nodes[offsets[i]..offsets[i + 1]]`, with an optional
/// length per list.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NodeLists {
    /// Start of each list in `nodes`, with a final entry equal to
    /// `nodes.len()`.
    pub offsets: Vec<usize>,
    /// The members of every list, concatenated.
    pub nodes: Vec<u32>,
    /// One value per list (a path length), or empty.
    pub values: Vec<f64>,
}

impl NodeLists {
    /// An empty collection.
    pub fn new() -> NodeLists {
        NodeLists {
            offsets: vec![0],
            nodes: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Appends a list.
    pub fn push(&mut self, list: &[u32]) {
        self.nodes.extend_from_slice(list);
        self.offsets.push(self.nodes.len());
    }

    /// Appends a list with its value.
    pub fn push_with(&mut self, list: &[u32], value: f64) {
        self.push(list);
        self.values.push(value);
    }

    /// Number of lists.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Whether there are no lists.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// List `i`.
    pub fn get(&self, i: usize) -> &[u32] {
        &self.nodes[self.offsets[i]..self.offsets[i + 1]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_arcs_for_each_mode() {
        // 0 -> 1 (2), 0 -> 1 (1), 1 -> 2 (4), 2 -> 2 (3).
        let g = Graph::from_edges(
            3,
            &[0, 0, 1, 2],
            &[1, 1, 2, 2],
            Some(&[2.0, 1.0, 4.0, 3.0]),
            true,
        )
        .unwrap();
        let w = ArcLengths::Weights(WeightKind::Distance);
        let out = PathGraph::new(
            &g,
            PathGraphOptions {
                lengths: w,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.arc_count(), 4);
        assert_eq!(out.out_arcs(0), 0..2);
        assert_eq!((out.edge(0), out.edge(1)), (0, 1));
        assert_eq!(out.arc_between(0, 1), Some(1));
        assert_eq!(out.path_length(&[0, 1, 2]), Some(5.0));
        assert_eq!(out.path_length(&[2, 1]), None);
        assert_eq!(out.in_arcs(1), &[0, 1]);
        let inn = PathGraph::new(
            &g,
            PathGraphOptions {
                mode: DegreeMode::In,
                lengths: w,
                collapse_parallel: true,
            },
        )
        .unwrap();
        assert_eq!(inn.arc_count(), 3);
        assert_eq!(inn.path_length(&[2, 1, 0]), Some(5.0));
        assert_eq!(inn.edge(inn.arc_between(1, 0).unwrap()), 1);
        let all = PathGraph::new(
            &g,
            PathGraphOptions {
                mode: DegreeMode::All,
                ..Default::default()
            },
        )
        .unwrap();
        // Two arcs per non-loop edge, one for the loop.
        assert_eq!(all.arc_count(), 7);
        assert!(all.is_unit());
        let hidden = [1.0, f64::INFINITY, 1.0, 1.0];
        let h = PathGraph::new(
            &g,
            PathGraphOptions {
                lengths: ArcLengths::PerEdge(&hidden),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(h.arc_count(), 3);
        assert!(PathGraph::new(
            &g,
            PathGraphOptions {
                lengths: ArcLengths::PerEdge(&[1.0, f64::NAN, 1.0, 1.0]),
                ..Default::default()
            },
        )
        .is_err());
    }

    #[test]
    fn tie_rules() {
        assert!(Tie::Relative.equal(0.1 + 0.2, 0.3));
        assert!(!Tie::Exact.equal(0.1 + 0.2, 0.3));
        assert!(Tie::Exact.equal(f64::INFINITY, f64::INFINITY));
        assert!(!Tie::Relative.equal(1.0, 1.0 + 1e-6));
    }
}
