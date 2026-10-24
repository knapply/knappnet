//! Variants of shortest-path betweenness: betweenness and edge
//! betweenness restricted to sets of sources and targets (which includes
//! betweenness estimated from a sample of sources), with or without the
//! endpoints of each path; percolation centrality; and load centrality
//! with NetworkX's edge load.
//!
//! Sources:
//! - **Betweenness for sources S and targets T** (Brandes 2008,
//!   `brandes2008variants`, sections 3.3 and 4): the sum over s in S and
//!   t in T, s != t, of sigma_st(v) / sigma_st, by Brandes' accumulation
//!   (Brandes 2001, `brandes2001faster`) with the dependency of t counted
//!   only when t is a target:
//!   delta_s(v) = sum over w with v a predecessor of w of
//!   (sigma_sv / sigma_sw)(1[w in T] + delta_s(w)). Endpoints (Brandes
//!   2008, section 3.2) add, for each pair, 1 to s and to t. Edges receive
//!   (sigma_sv / sigma_sw)(1[w in T] + delta_s(w)) per arc (v, w). With
//!   S = every node and T = every node this is ordinary betweenness; with S
//!   a random sample it is the estimator of Brandes and Pich (2007,
//!   `brandes2007centrality`).
//! - **Percolation centrality** (Piraveenan, Prokopenko, and Hossain 2013,
//!   `piraveenan2013percolation`): with percolation states x,
//!   PC(v) = 1/(n - 2) sum over s != v != t of
//!   (sigma_st(v) / sigma_st) x_s / (sum of x - x_v), computed with the same
//!   accumulation, each source's dependencies weighted.
//! - **Load centrality** (Goh, Kahng, and Kim 2001, `goh2001universal`;
//!   Brandes 2008, section 3.6): a unit of flow from each source to each
//!   other node is split equally among the predecessors of every node it
//!   passes, instead of in proportion to path counts:
//!   delta_s(v) = sum over w with v a predecessor of w of
//!   (1 + delta_s(w)) / |P_s(w)|, summed over sources. **Edge load** is
//!   NetworkX's `edge_load_centrality`, which NetworkX documents as a
//!   concept "not analysed or discussed outside of NetworkX"; it is
//!   reproduced exactly (see [`edge_load`]).
//!
//! Conventions:
//! - Direction: paths follow arcs in directed networks. Self-loops never
//!   lie on a shortest path.
//! - Parallel edges ([`Parallel`]): `Distinct` counts each parallel edge as
//!   a separate shortest path and gives each edge its own value, as igraph
//!   and the engine's [`crate::metrics::betweenness`]; `Collapse` counts a
//!   pair of adjacent nodes once (node sequences), as NetworkX, and edge
//!   values of a pair are shared equally by its parallel edges of smallest
//!   length (NetworkX's handling of multigraph keys); other parallel edges
//!   and self-loops get 0.
//! - Weights are lengths (distances, or strengths read as 1/w) and must be
//!   strictly positive. Tied path lengths are equal to a relative tolerance
//!   of 1e-10 ([`Tie::Relative`]) or only when exactly equal
//!   ([`Tie::Exact`], NetworkX).
//! - Scaling ([`Rescale`]): `Unordered` (the default) halves the sums over
//!   ordered pairs in undirected networks, so S = T = every node gives
//!   [`crate::metrics::betweenness::betweenness`]; normalized, node values
//!   divide by (n - 1)(n - 2) ordered pairs (n(n - 1) with endpoints) and
//!   edge values by n(n - 1), before halving. `NetworkX` is the rescaling of
//!   NetworkX's `betweenness_centrality` and `edge_betweenness_centrality`,
//!   including their correction for a sample of k sources (given as
//!   `sources`), and `NetworkXSubset` that of
//!   `betweenness_centrality_subset` and
//!   `edge_betweenness_centrality_subset`.
//! - NetworkX's `edge_betweenness_centrality_subset` divides the dependency
//!   of a node that is not a target equally among its predecessors,
//!   instead of in proportion to path counts as its documented formula
//!   (and Brandes 2008) does; [`EdgeSplit::NetworkX`] reproduces that, and
//!   the default [`EdgeSplit::Proportional`] follows the formula. The two
//!   agree whenever T contains every node or shortest paths do not branch.
//! - Percolation centrality is undefined (NaN) for n <= 2 (NetworkX divides
//!   by n - 2) and, for a node v reached from some source, when the states
//!   of all nodes other than v sum to 0.
//! - Load centrality sums over ordered pairs, also in undirected networks,
//!   as NetworkX's `load_centrality` does; normalized, it divides by
//!   (n - 1)(n - 2) (when n > 2). Predecessors are nodes, so parallel edges
//!   count once. `cutoff` limits paths to length at most `cutoff`. With
//!   `networkx_ties`, a node one of whose predecessors is the source passes
//!   nothing to its other predecessors, which is what NetworkX's
//!   `load_centrality` does (it stops at the source in each predecessor
//!   list, and the source always comes first); this differs from the
//!   definition only on weighted networks where the direct arc from the
//!   source ties with a longer path.
//! - Sources are summed in fixed groups merged in group order
//!   ([`crate::parallel::sum_in_order`]), so results are the same for
//!   every thread count.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;
use crate::parallel;
use crate::paths::{ArcLengths, PathGraph, PathGraphOptions, Tie};
use crate::weights::{IndexedHeap, WeightKind};

/// How parallel edges enter path counts. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Parallel {
    /// Each parallel edge is a distinct path (igraph).
    #[default]
    Distinct,
    /// Adjacent nodes are joined once (NetworkX).
    Collapse,
}

impl Parallel {
    /// Decodes the C interface codes: 0 = distinct, 1 = collapse.
    pub fn from_code(code: i32) -> Result<Parallel> {
        match code {
            0 => Ok(Parallel::Distinct),
            1 => Ok(Parallel::Collapse),
            _ => Err(Error::InvalidArgument(format!(
                "unknown code {code} for parallel edges"
            ))),
        }
    }
}

/// How sums over pairs are scaled. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Rescale {
    /// Unordered pairs in undirected networks (igraph).
    #[default]
    Unordered,
    /// NetworkX's `betweenness_centrality` and `edge_betweenness_centrality`.
    NetworkX,
    /// NetworkX's `*_betweenness_centrality_subset`.
    NetworkXSubset,
}

impl Rescale {
    /// Decodes the C interface codes: 0 = unordered, 1 = NetworkX,
    /// 2 = NetworkX subset.
    pub fn from_code(code: i32) -> Result<Rescale> {
        match code {
            0 => Ok(Rescale::Unordered),
            1 => Ok(Rescale::NetworkX),
            2 => Ok(Rescale::NetworkXSubset),
            _ => Err(Error::InvalidArgument(format!(
                "unknown rescaling code {code}"
            ))),
        }
    }
}

/// How an edge's share of a non-target dependency is computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EdgeSplit {
    /// In proportion to path counts (the definition).
    #[default]
    Proportional,
    /// Equally among predecessors (NetworkX's subset edge betweenness).
    NetworkX,
}

impl EdgeSplit {
    /// Decodes the C interface codes: 0 = proportional, 1 = NetworkX.
    pub fn from_code(code: i32) -> Result<EdgeSplit> {
        match code {
            0 => Ok(EdgeSplit::Proportional),
            1 => Ok(EdgeSplit::NetworkX),
            _ => Err(Error::InvalidArgument(format!(
                "unknown edge split code {code}"
            ))),
        }
    }
}

/// Options for [`betweenness_subset`].
#[derive(Debug, Clone, Copy, Default)]
pub struct SubsetOptions<'a> {
    /// Sources (every node when `None`).
    pub sources: Option<&'a [u32]>,
    /// Targets (every node when `None`).
    pub targets: Option<&'a [u32]>,
    /// Use edge weights as lengths.
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Count the endpoints of each path (node values only).
    pub endpoints: bool,
    /// Normalize (see [`Rescale`]).
    pub normalized: bool,
    /// Scaling convention.
    pub rescale: Rescale,
    /// Treatment of parallel edges.
    pub parallel: Parallel,
    /// Edge shares of non-target dependencies.
    pub edge_split: EdgeSplit,
    /// When path lengths tie.
    pub tie: Tie,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Node and edge values of [`betweenness_subset`].
#[derive(Debug, Clone, PartialEq)]
pub struct SubsetBetweenness {
    /// One value per node.
    pub nodes: Vec<f64>,
    /// One value per edge, in input order.
    pub edges: Vec<f64>,
}

fn path_graph(g: &Graph, weighted: bool, kind: WeightKind, collapse: bool) -> Result<PathGraph> {
    let lengths = if weighted {
        ArcLengths::Weights(kind)
    } else {
        ArcLengths::Unit
    };
    let pg = PathGraph::new(
        g,
        PathGraphOptions {
            mode: DegreeMode::Out,
            lengths,
            collapse_parallel: collapse,
        },
    )?;
    if weighted {
        pg.require_positive("shortest-path betweenness")?;
    }
    Ok(pg)
}

fn node_set(n: usize, nodes: Option<&[u32]>, what: &str) -> Result<Vec<u32>> {
    match nodes {
        None => Ok((0..n as u32).collect()),
        Some(list) => {
            if let Some(&v) = list.iter().find(|&&v| v as usize >= n) {
                return Err(Error::InvalidArgument(format!(
                    "{what} node index {v} is out of range for a network with {n} nodes"
                )));
            }
            Ok(list.to_vec())
        }
    }
}

/// The shortest-path DAG from one source: settling order, distances, path
/// counts, and the arcs from predecessors.
struct Dag {
    dist: Vec<f64>,
    sigma: Vec<f64>,
    rank: Vec<u32>,
    order: Vec<u32>,
    preds: Vec<Vec<u32>>,
    touched: Vec<u32>,
    heap: IndexedHeap,
}

impl Dag {
    fn new(n: usize) -> Dag {
        Dag {
            dist: vec![f64::INFINITY; n],
            sigma: vec![0.0; n],
            rank: vec![u32::MAX; n],
            order: Vec::with_capacity(n),
            preds: vec![Vec::new(); n],
            touched: Vec::new(),
            heap: IndexedHeap::new(n),
        }
    }

    fn reset(&mut self) {
        for &v in &self.touched {
            let v = v as usize;
            self.dist[v] = f64::INFINITY;
            self.sigma[v] = 0.0;
            self.rank[v] = u32::MAX;
            self.preds[v].clear();
        }
        self.touched.clear();
        self.order.clear();
        self.heap.clear();
    }

    /// Searches from `s` up to distance `cutoff`.
    fn run(&mut self, pg: &PathGraph, s: usize, tie: Tie, cutoff: f64) {
        self.reset();
        self.dist[s] = 0.0;
        self.sigma[s] = 1.0;
        self.touched.push(s as u32);
        if pg.is_unit() {
            self.rank[s] = 0;
            self.order.push(s as u32);
            let mut head = 0;
            while head < self.order.len() {
                let v = self.order[head] as usize;
                head += 1;
                if v != s {
                    self.sigma[v] = self.preds[v]
                        .iter()
                        .map(|&a| self.sigma[pg.tail(a as usize) as usize])
                        .sum();
                }
                let next = self.dist[v] + 1.0;
                if next > cutoff {
                    continue;
                }
                for a in pg.out_arcs(v) {
                    let w = pg.head(a) as usize;
                    if self.rank[w] == u32::MAX {
                        self.dist[w] = next;
                        self.rank[w] = self.order.len() as u32;
                        self.order.push(w as u32);
                        self.touched.push(w as u32);
                    }
                    if self.dist[w] == next {
                        self.preds[w].push(a as u32);
                    }
                }
            }
            return;
        }
        self.heap.push_or_decrease(s as u32, 0.0);
        while let Some((d, v)) = self.heap.pop() {
            let v = v as usize;
            self.rank[v] = self.order.len() as u32;
            self.order.push(v as u32);
            if v != s {
                self.sigma[v] = self.preds[v]
                    .iter()
                    .map(|&a| self.sigma[pg.tail(a as usize) as usize])
                    .sum();
            }
            for a in pg.out_arcs(v) {
                let w = pg.head(a) as usize;
                if self.rank[w] != u32::MAX {
                    continue;
                }
                let alt = d + pg.len(a);
                if alt > cutoff {
                    continue;
                }
                let cur = self.dist[w];
                if cur.is_infinite() {
                    self.touched.push(w as u32);
                    self.dist[w] = alt;
                    self.preds[w].push(a as u32);
                    self.heap.push_or_decrease(w as u32, alt);
                } else if tie.equal(alt, cur) {
                    self.preds[w].push(a as u32);
                } else if alt < cur {
                    self.dist[w] = alt;
                    self.preds[w].clear();
                    self.preds[w].push(a as u32);
                    self.heap.push_or_decrease(w as u32, alt);
                }
            }
        }
    }
}

/// Per-thread state of [`betweenness_subset`].
struct SubsetWork {
    dag: Dag,
    delta: Vec<f64>,
    delta_edge: Vec<f64>,
}

/// Betweenness and edge betweenness for sources and targets. See the
/// module documentation.
pub fn betweenness_subset(g: &Graph, opts: SubsetOptions) -> Result<SubsetBetweenness> {
    let n = g.node_count();
    let m = g.edge_count();
    let collapse = opts.parallel == Parallel::Collapse;
    let pg = path_graph(g, opts.weighted, opts.weights_as, collapse)?;
    let sources = node_set(n, opts.sources, "source")?;
    let targets = node_set(n, opts.targets, "target")?;
    let mut is_target = vec![false; n];
    for &t in &targets {
        is_target[t as usize] = true;
    }
    let quirk = opts.edge_split == EdgeSplit::NetworkX;
    let mut acc_nodes = vec![0.0; n];
    let mut acc_edges = vec![0.0; m];
    let threads = parallel::threads_for(
        opts.threads,
        sources.len(),
        parallel::search_work(sources.len(), n, pg.arc_count()),
    );
    parallel::sum_in_order(
        sources.len(),
        threads,
        || SubsetWork {
            dag: Dag::new(n),
            delta: vec![0.0; n],
            delta_edge: vec![0.0; n],
        },
        || (vec![0.0; n], vec![0.0; m]),
        |ws, acc: &mut (Vec<f64>, Vec<f64>), k| {
            let s = sources[k] as usize;
            ws.dag.run(&pg, s, opts.tie, f64::INFINITY);
            let dag = &ws.dag;
            let mut reached_targets = 0usize;
            for &w in dag.order.iter().rev() {
                let w = w as usize;
                let tw = if is_target[w] && w != s { 1.0 } else { 0.0 };
                let coeff = (tw + ws.delta[w]) / dag.sigma[w];
                let np = dag.preds[w].len() as f64;
                for &a in &dag.preds[w] {
                    let v = pg.tail(a as usize) as usize;
                    let c = dag.sigma[v] * coeff;
                    ws.delta[v] += c;
                    let ce = if quirk {
                        let c = if is_target[w] {
                            dag.sigma[v] / dag.sigma[w] * (1.0 + ws.delta_edge[w])
                        } else {
                            ws.delta_edge[w] / np
                        };
                        ws.delta_edge[v] += c;
                        c
                    } else {
                        c
                    };
                    acc.1[pg.edge(a as usize) as usize] += ce;
                }
                if w != s {
                    acc.0[w] += ws.delta[w];
                    if opts.endpoints && tw > 0.0 {
                        acc.0[w] += 1.0;
                        reached_targets += 1;
                    }
                }
            }
            if opts.endpoints {
                acc.0[s] += reached_targets as f64;
            }
            for &w in &dag.order {
                ws.delta[w as usize] = 0.0;
                ws.delta_edge[w as usize] = 0.0;
            }
        },
        |acc| {
            for (x, y) in acc_nodes.iter_mut().zip(acc.0.iter_mut()) {
                *x += *y;
                *y = 0.0;
            }
            for (x, y) in acc_edges.iter_mut().zip(acc.1.iter_mut()) {
                *x += *y;
                *y = 0.0;
            }
        },
    );

    // Scaling.
    let directed = g.is_directed();
    let mut source_set = vec![false; n];
    for &s in &sources {
        source_set[s as usize] = true;
    }
    let sampled = opts.sources.map(|_| sources.len());
    match opts.rescale {
        Rescale::Unordered => {
            let half = if directed { 1.0 } else { 0.5 };
            let (node_pairs, edge_pairs) = if opts.normalized {
                let nf = n as f64;
                let np = if opts.endpoints {
                    nf * (nf - 1.0)
                } else {
                    (nf - 1.0) * (nf - 2.0)
                };
                (np, nf * (nf - 1.0))
            } else {
                (1.0, 1.0)
            };
            let fn_ = if node_pairs > 0.0 {
                half / node_pairs
            } else {
                half
            };
            let fe = if edge_pairs > 0.0 {
                half / edge_pairs
            } else {
                half
            };
            acc_nodes.iter_mut().for_each(|x| *x *= fn_);
            acc_edges.iter_mut().for_each(|x| *x *= fe);
        }
        Rescale::NetworkX => {
            networkx_rescale(
                &mut acc_nodes,
                n,
                opts.normalized,
                directed,
                opts.endpoints,
                sampled,
                |v| source_set[v],
            );
            networkx_rescale(
                &mut acc_edges,
                n,
                opts.normalized,
                directed,
                true,
                sampled,
                |_| false,
            );
        }
        Rescale::NetworkXSubset => {
            networkx_rescale(
                &mut acc_nodes,
                n,
                opts.normalized,
                directed,
                false,
                None,
                |_| false,
            );
            networkx_rescale(
                &mut acc_edges,
                n,
                opts.normalized,
                directed,
                true,
                None,
                |_| false,
            );
        }
    }
    if collapse {
        let lengths: Vec<f64> = match (opts.weighted, g.edge_weights()) {
            (true, Some(w)) => match opts.weights_as {
                WeightKind::Distance => w.to_vec(),
                WeightKind::Strength => w.iter().map(|x| 1.0 / x).collect(),
            },
            _ => vec![1.0; m],
        };
        share_among_parallel(g, &lengths, &mut acc_edges);
    }
    Ok(SubsetBetweenness {
        nodes: acc_nodes,
        edges: acc_edges,
    })
}

/// NetworkX's `_rescale` for betweenness values (`networkx/algorithms/
/// centrality/betweenness.py`): `sampled` is the number of sampled sources
/// (None without sampling) and `is_source(v)` whether item v is one.
fn networkx_rescale(
    values: &mut [f64],
    n: usize,
    normalized: bool,
    directed: bool,
    endpoints: bool,
    sampled: Option<usize>,
    is_source: impl Fn(usize) -> bool,
) {
    let big_n = if endpoints { n as f64 } else { n as f64 - 1.0 };
    if big_n < 2.0 {
        return;
    }
    let correction = if directed { 1.0 } else { 2.0 };
    let k_source = sampled.map_or(big_n, |k| k as f64);
    if sampled.is_none() || endpoints {
        let scale = if normalized {
            1.0 / (k_source * (big_n - 1.0))
        } else {
            big_n / (k_source * correction)
        };
        if scale != 1.0 {
            values.iter_mut().for_each(|x| *x *= scale);
        }
        return;
    }
    let (scale_source, scale_other) = if normalized {
        (
            if k_source > 1.0 {
                1.0 / ((k_source - 1.0) * (big_n - 1.0))
            } else {
                f64::NAN
            },
            1.0 / (k_source * (big_n - 1.0)),
        )
    } else {
        (
            if k_source > 1.0 {
                big_n / ((k_source - 1.0) * correction)
            } else {
                f64::NAN
            },
            big_n / (k_source * correction),
        )
    };
    for (v, x) in values.iter_mut().enumerate() {
        *x *= if is_source(v) {
            scale_source
        } else {
            scale_other
        };
    }
}

/// Shares each adjacent pair's edge value equally among its parallel edges
/// of smallest length (NetworkX's `_add_edge_keys`); other parallel edges
/// and self-loops get 0. `values` holds the pair's value on the edge the
/// collapsed path graph kept; `lengths` are the edge lengths.
fn share_among_parallel(g: &Graph, lengths: &[f64], values: &mut [f64]) {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let m = g.edge_count();
    let mut groups: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for e in 0..m {
        let (a, b) = (src[e], tgt[e]);
        let k = if g.is_directed() || a <= b {
            (a, b)
        } else {
            (b, a)
        };
        groups.entry(k).or_default().push(e);
    }
    let mut out = vec![0.0; m];
    for (k, edges) in groups {
        if k.0 == k.1 {
            continue;
        }
        let total: f64 = edges.iter().map(|&e| values[e]).sum();
        let shortest = edges
            .iter()
            .map(|&e| lengths[e])
            .fold(f64::INFINITY, f64::min);
        let ties: Vec<usize> = edges
            .iter()
            .copied()
            .filter(|&e| lengths[e] == shortest)
            .collect();
        for &e in &ties {
            out[e] = total / ties.len() as f64;
        }
    }
    values.copy_from_slice(&out);
}

/// Options for [`percolation_centrality`].
#[derive(Debug, Clone, Copy, Default)]
pub struct PercolationOptions {
    /// Use edge weights as lengths.
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Treatment of parallel edges.
    pub parallel: Parallel,
    /// When path lengths tie.
    pub tie: Tie,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Percolation centrality with percolation states `states` (one per node,
/// in [0, 1] in the original definition). See the module documentation.
pub fn percolation_centrality(
    g: &Graph,
    states: &[f64],
    opts: PercolationOptions,
) -> Result<Vec<f64>> {
    let n = g.node_count();
    if states.len() != n {
        return Err(Error::LengthMismatch {
            what: "percolation states",
            expected: n,
            found: states.len(),
        });
    }
    if let Some(v) = states.iter().position(|x| !x.is_finite()) {
        return Err(Error::InvalidArgument(format!(
            "the percolation state of node {v} is not a finite number"
        )));
    }
    let pg = path_graph(
        g,
        opts.weighted,
        opts.weights_as,
        opts.parallel == Parallel::Collapse,
    )?;
    let total: f64 = states.iter().sum();
    let mut pc = vec![0.0; n];
    let threads =
        parallel::threads_for(opts.threads, n, parallel::search_work(n, n, pg.arc_count()));
    parallel::sum_in_order(
        n,
        threads,
        || (Dag::new(n), vec![0.0; n]),
        || vec![0.0; n],
        |(dag, delta), acc: &mut Vec<f64>, s| {
            dag.run(&pg, s, opts.tie, f64::INFINITY);
            for &w in dag.order.iter().rev() {
                let w = w as usize;
                let coeff = (1.0 + delta[w]) / dag.sigma[w];
                for &a in &dag.preds[w] {
                    let v = pg.tail(a as usize) as usize;
                    delta[v] += dag.sigma[v] * coeff;
                }
                if w != s {
                    let denom = total - states[w];
                    acc[w] += if denom == 0.0 {
                        f64::NAN
                    } else {
                        delta[w] * states[s] / denom
                    };
                }
            }
            for &w in &dag.order {
                delta[w as usize] = 0.0;
            }
        },
        |acc| {
            for (x, y) in pc.iter_mut().zip(acc.iter_mut()) {
                *x += *y;
                *y = 0.0;
            }
        },
    );
    let scale = if n > 2 {
        1.0 / (n as f64 - 2.0)
    } else {
        f64::NAN
    };
    Ok(pc.into_iter().map(|x| x * scale).collect())
}

/// Options for [`load_centrality`].
#[derive(Debug, Clone, Copy)]
pub struct LoadOptions {
    /// Use edge weights as lengths.
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Divide by (n - 1)(n - 2).
    pub normalized: bool,
    /// Longest path length considered (infinity for no limit).
    pub cutoff: f64,
    /// Reproduce NetworkX's handling of the source among several
    /// predecessors (see the module documentation).
    pub networkx_ties: bool,
    /// When path lengths tie.
    pub tie: Tie,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for LoadOptions {
    fn default() -> Self {
        LoadOptions {
            weighted: false,
            weights_as: WeightKind::Distance,
            normalized: false,
            cutoff: f64::INFINITY,
            networkx_ties: false,
            tie: Tie::Relative,
            threads: 0,
        }
    }
}

/// Load centrality of every node. See the module documentation.
pub fn load_centrality(g: &Graph, opts: LoadOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    if opts.cutoff.is_nan() {
        return Err(Error::InvalidArgument("the cutoff must not be NaN".into()));
    }
    let pg = path_graph(g, opts.weighted, opts.weights_as, true)?;
    let mut load = vec![0.0; n];
    let threads =
        parallel::threads_for(opts.threads, n, parallel::search_work(n, n, pg.arc_count()));
    parallel::sum_in_order(
        n,
        threads,
        || (Dag::new(n), vec![0.0; n]),
        || vec![0.0; n],
        |(dag, between), acc: &mut Vec<f64>, s| {
            dag.run(&pg, s, opts.tie, opts.cutoff);
            for &w in &dag.order {
                between[w as usize] = 1.0;
            }
            for &v in dag.order.iter().rev() {
                let v = v as usize;
                if v == s {
                    continue;
                }
                let preds = &dag.preds[v];
                let share = between[v] / preds.len() as f64;
                let from_source = preds.iter().any(|&a| pg.tail(a as usize) as usize == s);
                if !(opts.networkx_ties && from_source) {
                    for &a in preds {
                        let x = pg.tail(a as usize) as usize;
                        if x != s {
                            between[x] += share;
                        }
                    }
                }
            }
            for &w in &dag.order {
                let w = w as usize;
                if w != s {
                    acc[w] += between[w] - 1.0;
                }
                between[w] = 0.0;
            }
        },
        |acc| {
            for (x, y) in load.iter_mut().zip(acc.iter_mut()) {
                *x += *y;
                *y = 0.0;
            }
        },
    );
    if opts.normalized && n > 2 {
        let scale = 1.0 / ((n as f64 - 1.0) * (n as f64 - 2.0));
        load.iter_mut().for_each(|x| *x *= scale);
    }
    Ok(load)
}

/// NetworkX's edge load (`edge_load_centrality`), reproduced exactly: two
/// values per edge, in input order, for the ordered pairs (source, target)
/// and (target, source) of the edge (NetworkX keys its result by both).
/// Weights are ignored; `cutoff` limits the breadth-first searches to that
/// many steps (infinity for no limit).
///
/// For each source, every ordered pair of adjacent nodes starts at 1, and
/// nodes are taken in order of decreasing distance; for a node v and each
/// predecessor w of v, every predecessor x of w receives
/// b(v, w) / |P(w)| on the pair (w, x) and b(w, v) / |P(w)| on (x, w).
/// The values are summed over sources, so every pair has at least n.
pub fn edge_load(g: &Graph, cutoff: f64, threads: usize) -> Result<Vec<f64>> {
    let n = g.node_count();
    let m = g.edge_count();
    if cutoff.is_nan() {
        return Err(Error::InvalidArgument("the cutoff must not be NaN".into()));
    }
    let pg = path_graph(g, false, WeightKind::Distance, true)?;
    // Ordered pairs of adjacent nodes, each with an index.
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut index: HashMap<(u32, u32), usize> = HashMap::new();
    let mut slot = Vec::with_capacity(2 * m);
    for e in 0..m {
        for p in [(src[e], tgt[e]), (tgt[e], src[e])] {
            let next = index.len();
            slot.push(*index.entry(p).or_insert(next));
        }
    }
    let pairs = index.len();
    let mut total = vec![0.0; pairs];
    let threads = parallel::threads_for(threads, n, parallel::search_work(n, n, pg.arc_count()));
    parallel::sum_in_order(
        n,
        threads,
        || (Dag::new(n), vec![1.0; pairs]),
        || vec![0.0; pairs],
        |(dag, between), acc: &mut Vec<f64>, s| {
            dag.run(&pg, s, Tie::Exact, cutoff);
            let pair = |a: u32, b: u32| index[&(a, b)];
            for &v in dag.order.iter().rev() {
                let v = v as usize;
                for &a in &dag.preds[v] {
                    let w = pg.tail(a as usize) as usize;
                    let pw = &dag.preds[w];
                    if pw.is_empty() {
                        continue;
                    }
                    let np = pw.len() as f64;
                    let b_vw = between[pair(v as u32, w as u32)];
                    let b_wv = between[pair(w as u32, v as u32)];
                    for &a2 in pw {
                        let x = pg.tail(a2 as usize);
                        between[pair(w as u32, x)] += b_vw / np;
                        between[pair(x, w as u32)] += b_wv / np;
                    }
                }
            }
            for (t, b) in acc.iter_mut().zip(between.iter_mut()) {
                *t += *b;
                *b = 1.0;
            }
        },
        |acc| {
            for (x, y) in total.iter_mut().zip(acc.iter_mut()) {
                *x += *y;
                *y = 0.0;
            }
        },
    );
    Ok(slot.into_iter().map(|i| total[i]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::betweenness::{betweenness, BetweennessOptions};
    use crate::metrics::reach::{edge_betweenness, EdgeBetweennessOptions};

    fn g(n: usize, s: &[i64], t: &[i64], w: Option<&[f64]>, directed: bool) -> Graph {
        Graph::from_edges(n, s, t, w, directed).unwrap()
    }

    #[test]
    fn full_subset_is_ordinary_betweenness() {
        let h = g(
            6,
            &[0, 0, 1, 2, 3, 4, 1, 1],
            &[1, 2, 3, 3, 4, 5, 3, 1],
            Some(&[1.0, 2.0, 2.0, 1.0, 1.0, 3.0, 2.0, 1.0]),
            false,
        );
        for weighted in [false, true] {
            let sub = betweenness_subset(
                &h,
                SubsetOptions {
                    weighted,
                    threads: 2,
                    ..Default::default()
                },
            )
            .unwrap();
            let bc = betweenness(
                &h,
                BetweennessOptions {
                    weighted,
                    ..Default::default()
                },
            )
            .unwrap();
            let eb = edge_betweenness(
                &h,
                EdgeBetweennessOptions {
                    weighted,
                    ..Default::default()
                },
            )
            .unwrap();
            for (a, b) in sub.nodes.iter().zip(&bc) {
                assert!((a - b).abs() < 1e-12, "{:?} {:?}", sub.nodes, bc);
            }
            for (a, b) in sub.edges.iter().zip(&eb) {
                assert!((a - b).abs() < 1e-12, "{:?} {:?}", sub.edges, eb);
            }
        }
    }

    #[test]
    fn load_on_a_path_and_a_diamond() {
        // Path 0-1-2: node 1 carries load between 0 and 2 in both
        // directions.
        let p = g(3, &[0, 1], &[1, 2], None, false);
        assert_eq!(
            load_centrality(&p, LoadOptions::default()).unwrap(),
            vec![0.0, 2.0, 0.0]
        );
        // Diamond 0-1-3, 0-2-3: equal split.
        let d = g(4, &[0, 0, 1, 2], &[1, 2, 3, 3], None, false);
        assert_eq!(
            load_centrality(&d, LoadOptions::default()).unwrap(),
            vec![1.0, 1.0, 1.0, 1.0]
        );
    }

    #[test]
    fn percolation_with_equal_states_is_scaled_betweenness() {
        let h = g(5, &[0, 1, 2, 3, 1], &[1, 2, 3, 4, 3], None, false);
        let pc = percolation_centrality(&h, &[1.0; 5], PercolationOptions::default()).unwrap();
        let bc = betweenness(&h, BetweennessOptions::default()).unwrap();
        // x_s / (sum - x_v) = 1/4 and ordered pairs double the unordered
        // betweenness: PC = 2 bc / 4 / (n - 2).
        for (a, b) in pc.iter().zip(&bc) {
            assert!((a - 2.0 * b / 4.0 / 3.0).abs() < 1e-12);
        }
    }
}
