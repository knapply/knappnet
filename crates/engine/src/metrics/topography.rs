//! Network topography: density, components, fragmentation, and Freeman
//! centralization.
//!
//! Sources:
//! - Density is the share of possible ties that are present (Wasserman and
//!   Faust 1994, `wasserman1994social`, section 4.2).
//! - Strong components use Tarjan's depth-first algorithm (Tarjan 1972,
//!   `tarjan1972depth`), written with an explicit stack so deep graphs cannot
//!   overflow the call stack. Weak components use breadth-first search over
//!   both arc directions.
//! - Fragmentation is the proportion of ordered pairs of distinct nodes that
//!   lie in different components (Borgatti 2006, `borgatti2006identifying`):
//!   F = 1 - sum_k s_k (s_k - 1) / (n (n - 1)), where s_k are component
//!   sizes.
//! - Centralization is Freeman's index (Freeman 1978, `freeman1978centrality`):
//!   sum_i (max_j c_j - c_i), divided by the largest value that sum can take
//!   in a network with the same number of nodes. The maxima of each reference
//!   tool are recorded in [`crate::conventions`].
//!
//! Conventions:
//! - Density counts every edge, including parallel edges and self-loops, and
//!   divides by n(n - 1) (directed) or n(n - 1)/2 (undirected); with
//!   `loops = true` the denominator also counts the n possible self-loops,
//!   n^2 or n(n + 1)/2. This is igraph's `edge_density()`, and NetworkX's
//!   `density()` for n > 1. With parallel edges or loops the value can
//!   exceed 1. sna's `gden()` works on the adjacency matrix, which is
//!   [`simplify`] followed by this function. A network with no nodes, or
//!   one node and `loops = false`, has density NaN (igraph and sna; NetworkX
//!   returns 0).
//! - Weights are ignored by every measure here.
//! - Components: in undirected networks weak and strong components coincide.
//!   Isolated nodes are components of size 1. Component labels are numbered
//!   from 0 in order of each component's lowest-numbered node.
//! - Fragmentation of a network with fewer than two nodes is NaN (no pairs).

use crate::conventions::Convention;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::betweenness::{betweenness, BetweennessOptions};
use crate::metrics::closeness::{closeness, ClosenessOptions};
use crate::metrics::degree::{degree, DegreeMode};

/// Share of possible edges present. See the module documentation.
pub fn density(g: &Graph, loops: bool) -> f64 {
    let n = g.node_count() as f64;
    let m = g.edge_count() as f64;
    if n == 0.0 || (n == 1.0 && !loops) {
        return f64::NAN;
    }
    let possible = match (g.is_directed(), loops) {
        (true, false) => n * (n - 1.0),
        (true, true) => n * n,
        (false, false) => n * (n - 1.0) / 2.0,
        (false, true) => n * (n + 1.0) / 2.0,
    };
    m / possible
}

/// Which weight a set of parallel edges keeps when [`simplify`] merges them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ParallelWeight {
    /// The largest (UCINET's symmetrizing by maximum).
    #[default]
    Max,
    /// The smallest: the shortest of several distances, which is what
    /// NetworkX's shortest-path searches use on a multigraph.
    Min,
}

impl ParallelWeight {
    /// Decodes the C interface codes: 0 = max, 1 = min.
    pub fn from_code(code: i32) -> Result<ParallelWeight> {
        match code {
            0 => Ok(ParallelWeight::Max),
            1 => Ok(ParallelWeight::Min),
            _ => Err(Error::InvalidArgument(format!(
                "unknown parallel-edge weight code {code}"
            ))),
        }
    }
}

/// The adjacency matrix of `g` as a graph: parallel edges become one edge,
/// keeping the largest or smallest of their weights (if any), and
/// self-loops are dropped. With `undirected = true` a directed graph is
/// symmetrized first (i and j are joined if either arc is present), which
/// with [`ParallelWeight::Max`] is UCINET's symmetrizing by maximum. This
/// is how sna and UCINET, which work on adjacency matrices, read a network.
pub fn simplify(g: &Graph, undirected: bool, keep: ParallelWeight) -> Result<Graph> {
    let directed = g.is_directed() && !undirected;
    let w = g.edge_weights();
    let mut arcs: Vec<(u32, u32, f64)> = g
        .edge_sources()
        .iter()
        .zip(g.edge_targets())
        .enumerate()
        .filter(|(_, (s, t))| s != t)
        .map(|(e, (&s, &t))| {
            let x = w.map_or(1.0, |w| w[e]);
            if directed || s < t {
                (s, t, x)
            } else {
                (t, s, x)
            }
        })
        .collect();
    arcs.sort_unstable_by_key(|a| (a.0, a.1));
    let mut merged: Vec<(u32, u32, f64)> = Vec::with_capacity(arcs.len());
    for (s, t, x) in arcs {
        match merged.last_mut() {
            Some(last) if last.0 == s && last.1 == t => {
                last.2 = match keep {
                    ParallelWeight::Max => last.2.max(x),
                    ParallelWeight::Min => last.2.min(x),
                }
            }
            _ => merged.push((s, t, x)),
        }
    }
    let src: Vec<i64> = merged.iter().map(|p| p.0 as i64).collect();
    let tgt: Vec<i64> = merged.iter().map(|p| p.1 as i64).collect();
    let weights: Option<Vec<f64>> = w.map(|_| merged.iter().map(|p| p.2).collect());
    Graph::from_edges(g.node_count(), &src, &tgt, weights.as_deref(), directed)
}

/// Which kind of component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ComponentKind {
    /// Connected when arc directions are ignored.
    #[default]
    Weak,
    /// Every node reaches every other along directed paths.
    Strong,
}

impl ComponentKind {
    /// Decodes the C interface codes: 0 = weak, 1 = strong.
    pub fn from_code(code: i32) -> Result<ComponentKind> {
        match code {
            0 => Ok(ComponentKind::Weak),
            1 => Ok(ComponentKind::Strong),
            _ => Err(Error::InvalidArgument(format!(
                "unknown component kind code {code}"
            ))),
        }
    }
}

/// Component membership and sizes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Components {
    /// Component label of every node, numbered from 0 in order of each
    /// component's lowest-numbered node.
    pub membership: Vec<u32>,
    /// Number of nodes in each component, indexed by label.
    pub sizes: Vec<usize>,
}

impl Components {
    /// Number of components.
    pub fn count(&self) -> usize {
        self.sizes.len()
    }
}

/// Weak or strong components. See the module documentation.
pub fn components(g: &Graph, kind: ComponentKind) -> Components {
    let raw = if kind == ComponentKind::Strong && g.is_directed() {
        tarjan(g)
    } else {
        weak(g)
    };
    // Relabel in order of first appearance, which is the order of each
    // component's lowest-numbered node.
    let mut relabel = vec![u32::MAX; g.node_count()];
    let mut sizes = Vec::new();
    let membership = raw
        .into_iter()
        .map(|c| {
            let c = c as usize;
            if relabel[c] == u32::MAX {
                relabel[c] = sizes.len() as u32;
                sizes.push(0);
            }
            let label = relabel[c];
            sizes[label as usize] += 1;
            label
        })
        .collect();
    Components { membership, sizes }
}

fn weak(g: &Graph) -> Vec<u32> {
    let n = g.node_count();
    let mut label = vec![u32::MAX; n];
    let mut queue: Vec<u32> = Vec::new();
    let mut next = 0u32;
    for s in 0..n {
        if label[s] != u32::MAX {
            continue;
        }
        label[s] = next;
        queue.clear();
        queue.push(s as u32);
        let mut head = 0;
        while head < queue.len() {
            let v = queue[head] as usize;
            head += 1;
            let both = g.out().neighbors(v).iter().chain(if g.is_directed() {
                g.inc().neighbors(v)
            } else {
                &[]
            });
            for &w in both {
                if label[w as usize] == u32::MAX {
                    label[w as usize] = next;
                    queue.push(w);
                }
            }
        }
        next += 1;
    }
    label
}

/// Tarjan's algorithm with an explicit stack. Returns a component number
/// per node (in the order Tarjan completes them, relabeled by the caller).
fn tarjan(g: &Graph) -> Vec<u32> {
    let n = g.node_count();
    let adj = g.out();
    const UNSEEN: u32 = u32::MAX;
    let mut index = vec![UNSEEN; n];
    let mut low = vec![0u32; n];
    let mut on_stack = vec![false; n];
    let mut comp = vec![UNSEEN; n];
    let mut stack: Vec<u32> = Vec::new();
    // Each frame is a node and the position of the next arc to explore.
    let mut frames: Vec<(u32, usize)> = Vec::new();
    let mut counter = 0u32;
    let mut ncomp = 0u32;

    for root in 0..n {
        if index[root] != UNSEEN {
            continue;
        }
        frames.push((root as u32, 0));
        index[root] = counter;
        low[root] = counter;
        counter += 1;
        stack.push(root as u32);
        on_stack[root] = true;

        while let Some(&mut (v, ref mut pos)) = frames.last_mut() {
            let v = v as usize;
            let nbrs = adj.neighbors(v);
            if *pos < nbrs.len() {
                let w = nbrs[*pos] as usize;
                *pos += 1;
                if index[w] == UNSEEN {
                    index[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push(w as u32);
                    on_stack[w] = true;
                    frames.push((w as u32, 0));
                } else if on_stack[w] {
                    low[v] = low[v].min(index[w]);
                }
                continue;
            }
            // All arcs of v explored.
            frames.pop();
            if let Some(&(parent, _)) = frames.last() {
                let p = parent as usize;
                low[p] = low[p].min(low[v]);
            }
            if low[v] == index[v] {
                loop {
                    let w = stack.pop().expect("Tarjan stack underflow") as usize;
                    on_stack[w] = false;
                    comp[w] = ncomp;
                    if w == v {
                        break;
                    }
                }
                ncomp += 1;
            }
        }
    }
    comp
}

/// Proportion of ordered pairs of distinct nodes in different components
/// (Borgatti 2006). See the module documentation.
pub fn fragmentation(g: &Graph, kind: ComponentKind) -> f64 {
    let n = g.node_count() as u128;
    if n < 2 {
        return f64::NAN;
    }
    let c = components(g, kind);
    let within: u128 = c
        .sizes
        .iter()
        .map(|&s| (s as u128) * (s as u128).saturating_sub(1))
        .sum();
    let pairs = n * (n - 1);
    // Both counts are exact integers; divide once.
    (pairs - within) as f64 / pairs as f64
}

/// Which centrality a centralization index summarizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CentralizationMeasure {
    /// Degree (counts of edge endpoints).
    Degree,
    /// Shortest-path betweenness, not normalized.
    Betweenness,
    /// Normalized closeness.
    Closeness,
}

impl CentralizationMeasure {
    /// Decodes the C interface codes: 0 = degree, 1 = betweenness,
    /// 2 = closeness.
    pub fn from_code(code: i32) -> Result<CentralizationMeasure> {
        match code {
            0 => Ok(CentralizationMeasure::Degree),
            1 => Ok(CentralizationMeasure::Betweenness),
            2 => Ok(CentralizationMeasure::Closeness),
            _ => Err(Error::InvalidArgument(format!(
                "unknown centralization measure code {code}"
            ))),
        }
    }
}

/// Options for [`centralization`].
#[derive(Debug, Clone, Copy)]
pub struct CentralizationOptions {
    /// The centrality being summarized.
    pub measure: CentralizationMeasure,
    /// For directed graphs: which degree, or which direction of paths for
    /// closeness. Ignored for betweenness and for undirected graphs.
    pub mode: DegreeMode,
    /// Whether self-loops count toward degree and the degree maximum
    /// (igraph only; sna and UCINET never count them).
    pub loops: bool,
    /// Divide by the theoretical maximum.
    pub normalized: bool,
    /// Whose maxima and closeness rule to use. NetworkX has no
    /// centralization and is rejected.
    pub convention: Convention,
    /// Worker threads for betweenness and closeness; 0 uses all cores.
    pub threads: usize,
}

impl Default for CentralizationOptions {
    fn default() -> Self {
        CentralizationOptions {
            measure: CentralizationMeasure::Degree,
            mode: DegreeMode::All,
            loops: true,
            normalized: true,
            convention: Convention::Igraph,
            threads: 0,
        }
    }
}

/// Freeman centralization of a network. See the module documentation and
/// [`crate::conventions`].
///
/// The convention selects the theoretical maximum and, for closeness, how
/// unreachable nodes are treated. Under the sna and UCINET conventions the
/// network is read as an adjacency matrix first ([`simplify`]), so parallel
/// edges collapse and self-loops are ignored, as those tools do.
pub fn centralization(g: &Graph, opts: CentralizationOptions) -> Result<f64> {
    use CentralizationMeasure as M;
    let conv = opts.convention;
    if conv == Convention::Networkx {
        return Err(Error::InvalidArgument(
            "NetworkX has no centralization index; use the igraph, sna, or UCINET convention"
                .into(),
        ));
    }
    let simple;
    let g = if conv.uses_adjacency_matrix() {
        simple = simplify(g, false, ParallelWeight::Max)?;
        &simple
    } else {
        g
    };
    let directed = g.is_directed();
    let mode = if directed { opts.mode } else { DegreeMode::All };
    let loops = opts.loops && conv == Convention::Igraph;
    let scores = match opts.measure {
        M::Degree => {
            if conv == Convention::Ucinet && directed && mode == DegreeMode::All {
                return Err(Error::InvalidArgument(
                    "UCINET reports degree centralization for in- and out-degree only".into(),
                ));
            }
            degree(g, mode, loops)
        }
        M::Betweenness => betweenness(
            g,
            BetweennessOptions {
                threads: opts.threads,
                ..Default::default()
            },
        )?,
        M::Closeness => {
            let copts = ClosenessOptions {
                mode,
                normalized: true,
                unreachable: match conv {
                    Convention::Sna => crate::metrics::closeness::Unreachable::Zero,
                    _ => crate::metrics::closeness::Unreachable::Ignore,
                },
                threads: opts.threads,
                ..Default::default()
            };
            let c = closeness(g, copts)?;
            // UCINET reports closeness centralization only when every node
            // reaches every other.
            if conv == Convention::Ucinet {
                let n = g.node_count();
                let full = crate::metrics::closeness::distance_profiles(g, &copts)?
                    .iter()
                    .all(|p| p.reached + 1 == n);
                if !full {
                    return Ok(f64::NAN);
                }
            }
            c
        }
    };
    let n = g.node_count();
    if n == 0 {
        return Ok(f64::NAN);
    }
    let max = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = scores.iter().sum();
    // NaN scores (closeness of a node that reaches no other) propagate.
    let max = if scores.iter().any(|x| x.is_nan()) {
        f64::NAN
    } else {
        max
    };
    let raw = n as f64 * max - total;
    if !opts.normalized {
        return Ok(raw);
    }
    let tmax = centralization_max(n, directed, opts.measure, mode, opts.loops, conv)?;
    Ok(raw / tmax)
}

/// Theoretical maximum of the centralization sum for a network with `n`
/// nodes, as each tool defines it (see [`crate::conventions`]). `loops`
/// matters only for degree under the igraph convention.
pub fn centralization_max(
    n: usize,
    directed: bool,
    measure: CentralizationMeasure,
    mode: DegreeMode,
    loops: bool,
    convention: Convention,
) -> Result<f64> {
    use CentralizationMeasure as M;
    if n == 0 {
        return Ok(f64::NAN);
    }
    let n = n as f64;
    let mode = if directed { mode } else { DegreeMode::All };
    let loops = loops && convention == Convention::Igraph;
    Ok(match (measure, convention) {
        (_, Convention::Networkx) => {
            return Err(Error::InvalidArgument(
                "NetworkX has no centralization index".into(),
            ))
        }
        (M::Degree, _) => match (directed, mode, loops) {
            (false, _, false) => (n - 1.0) * (n - 2.0),
            (false, _, true) => (n - 1.0) * n,
            (true, DegreeMode::All, false) => 2.0 * (n - 1.0) * (n - 2.0),
            (true, DegreeMode::All, true) => 2.0 * (n - 1.0) * (n - 1.0),
            (true, _, false) => (n - 1.0) * (n - 1.0),
            (true, _, true) => (n - 1.0) * n,
        },
        (M::Betweenness, _) => {
            let t = (n - 1.0) * (n - 1.0) * (n - 2.0);
            if directed {
                t
            } else {
                t / 2.0
            }
        }
        (M::Closeness, Convention::Ucinet) => (n - 1.0) * (n - 2.0) / (2.0 * n - 3.0),
        (M::Closeness, _) => {
            if mode == DegreeMode::All {
                (n - 1.0) * (n - 2.0) / (2.0 * n - 3.0)
            } else {
                (n - 1.0) * (1.0 - 1.0 / n)
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(n: usize, s: &[i64], t: &[i64], directed: bool) -> Graph {
        Graph::from_edges(n, s, t, None, directed).unwrap()
    }

    #[test]
    fn density_formulas() {
        let u = g(4, &[0, 1, 2], &[1, 2, 3], false);
        assert_eq!(density(&u, false), 0.5);
        assert_eq!(density(&u, true), 0.3);
        let d = g(3, &[0, 1, 1], &[1, 2, 1], true);
        assert_eq!(density(&d, false), 0.5);
        assert_eq!(density(&d, true), 3.0 / 9.0);
        assert!(density(&g(1, &[], &[], false), false).is_nan());
        assert_eq!(density(&g(1, &[], &[], false), true), 0.0);
        assert!(density(&g(0, &[], &[], true), true).is_nan());
    }

    #[test]
    fn simplify_collapses_and_symmetrizes() {
        let w = [1.0, 4.0, 2.0, 9.0, 3.0];
        let d = Graph::from_edges(3, &[0, 0, 1, 1, 2], &[1, 1, 0, 1, 0], Some(&w), true).unwrap();
        let s = simplify(&d, false, ParallelWeight::Max).unwrap();
        assert_eq!(s.edge_sources(), &[0, 1, 2]);
        assert_eq!(s.edge_targets(), &[1, 0, 0]);
        assert_eq!(s.edge_weights().unwrap(), &[4.0, 2.0, 3.0]);
        assert!(s.is_directed());
        let u = simplify(&d, true, ParallelWeight::Max).unwrap();
        assert_eq!(u.edge_count(), 2);
        assert_eq!(u.edge_weights().unwrap(), &[4.0, 3.0]);
        assert!(!u.is_directed());
        let m = simplify(&d, true, ParallelWeight::Min).unwrap();
        assert_eq!(m.edge_weights().unwrap(), &[1.0, 3.0]);
    }

    #[test]
    fn weak_and_strong_components() {
        // 0 -> 1 -> 2 -> 0, 2 -> 3, 4 isolated, 5 <-> 6.
        let d = g(7, &[0, 1, 2, 2, 5, 6], &[1, 2, 0, 3, 6, 5], true);
        let w = components(&d, ComponentKind::Weak);
        assert_eq!(w.membership, vec![0, 0, 0, 0, 1, 2, 2]);
        assert_eq!(w.sizes, vec![4, 1, 2]);
        let s = components(&d, ComponentKind::Strong);
        assert_eq!(s.membership, vec![0, 0, 0, 1, 2, 3, 3]);
        assert_eq!(s.sizes, vec![3, 1, 1, 2]);
        assert_eq!(s.count(), 4);
    }

    #[test]
    fn deep_path_does_not_overflow() {
        let n = 200_000;
        let s: Vec<i64> = (0..n as i64 - 1).collect();
        let t: Vec<i64> = (1..n as i64).collect();
        let d = g(n, &s, &t, true);
        assert_eq!(components(&d, ComponentKind::Strong).count(), n);
        assert_eq!(components(&d, ComponentKind::Weak).count(), 1);
    }

    #[test]
    fn fragmentation_counts_pairs() {
        // Components of sizes 3, 3, 1: 12 of 42 ordered pairs connected.
        let u = g(7, &[0, 1, 2, 3, 4, 5], &[1, 2, 0, 4, 5, 3], false);
        assert!((fragmentation(&u, ComponentKind::Weak) - 30.0 / 42.0).abs() < 1e-15);
        assert!(fragmentation(&g(1, &[], &[], false), ComponentKind::Weak).is_nan());
        assert_eq!(
            fragmentation(&g(3, &[0, 1], &[1, 2], false), ComponentKind::Weak),
            0.0
        );
    }

    #[test]
    fn star_is_maximally_centralized() {
        let star = g(6, &[0, 0, 0, 0, 0], &[1, 2, 3, 4, 5], false);
        for measure in [
            CentralizationMeasure::Degree,
            CentralizationMeasure::Betweenness,
            CentralizationMeasure::Closeness,
        ] {
            for convention in [Convention::Igraph, Convention::Sna, Convention::Ucinet] {
                let opts = CentralizationOptions {
                    measure,
                    convention,
                    loops: false,
                    ..Default::default()
                };
                let c = centralization(&star, opts).unwrap();
                assert!((c - 1.0).abs() < 1e-12, "{measure:?} {convention:?}: {c}");
            }
        }
        let cycle = g(5, &[0, 1, 2, 3, 4], &[1, 2, 3, 4, 0], false);
        let c = centralization(&cycle, CentralizationOptions::default()).unwrap();
        assert_eq!(c, 0.0);
        let nx = CentralizationOptions {
            convention: Convention::Networkx,
            ..Default::default()
        };
        assert!(centralization(&cycle, nx).is_err());
    }
}
