//! Minimum node and edge cuts: smallest sets of nodes or edges whose
//! removal separates two nodes, or disconnects the network.
//!
//! Sources: Menger (1927, `menger1927allgemeinen`) and Ford and Fulkerson
//! (1956, `ford1956maximal`) for the duality with disjoint paths; the
//! cuts are read off maximum flows in the edge network and the split
//! network (Even and Tarjan 1975, `even1975network`), and the global cuts
//! are searched over the same pairs as the global connectivities in
//! [`super::values`] (Even 1975, `even1975algorithm`; Esfahanian 2013,
//! `esfahanian2013connectivity`).
//!
//! Minimum cuts are generally not unique. Two conventions are offered
//! ([`CutConvention`]):
//! - `Engine`: the s-t edge cut is the set of arcs leaving the largest
//!   source side of a minimum cut (canonical, independent of the flow
//!   found). The s-t node cut is read off the split network with arcs of
//!   infinite capacity between nodes, so it consists of nodes only: the
//!   nodes whose in-copy is reachable from s in the residual network and
//!   whose out-copy is not (the cut nearest s). Adjacent s and t have no
//!   node cut, which is an error.
//! - `Networkx`: NetworkX's `minimum_st_edge_cut` and
//!   `minimum_st_node_cut` (`hagberg2008exploring`). The edge cut is the
//!   same as the engine's. The node cut is computed on NetworkX's split
//!   network, where the arcs between nodes also have capacity 1: every
//!   arc leaving the largest source side of a minimum cut is taken, and
//!   the cut is the set of their endpoints other than s and t; for
//!   adjacent s and t it is empty. Because NetworkX's minimum cut side is
//!   canonical, this reproduces its result exactly.
//!
//! Global cuts: an empty set for a network that is not connected (weakly,
//! if directed); NetworkX raises an error instead. For a complete graph,
//! which no node set disconnects, the global node cut is the set of
//! neighbours of the first node of least degree, as in NetworkX. Among
//! several minimum cuts the first one found in the order of the pairs in
//! [`super::values`] is returned (NetworkX's choice among equally small
//! cuts depends on the iteration order of Python sets and can differ).
//!
//! Edge cuts are reported as node pairs (u, v) of the simple graph, u on
//! the source side (for an undirected network, one pair per edge), sorted.

use super::values::even_pairs;
use super::{check_pair, Simple};
use crate::error::{Error, Result};
use crate::flow::FlowNetwork;
use crate::graph::Graph;
use crate::parallel;

/// Per-thread state of a search over pairs: a network and the cuts found,
/// tagged with the index of their pair.
type Found<T> = (FlowNetwork, Vec<(usize, Vec<T>)>);

/// Which minimum cut to report. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CutConvention {
    /// The engine's convention.
    #[default]
    Engine,
    /// NetworkX's.
    Networkx,
}

impl CutConvention {
    /// Decodes the C interface codes: 0 = engine, 1 = NetworkX.
    pub fn from_code(code: i32) -> Result<CutConvention> {
        match code {
            0 => Ok(CutConvention::Engine),
            1 => Ok(CutConvention::Networkx),
            _ => Err(Error::InvalidArgument(format!(
                "unknown cut convention code {code}"
            ))),
        }
    }
}

/// Arcs of the simple graph leaving `side`.
fn arcs_leaving(simple: &Simple, side: &[bool]) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for u in 0..simple.n {
        if !side[u] {
            continue;
        }
        for &v in &simple.out[u] {
            if !side[v as usize] {
                out.push((u as u32, v));
            }
        }
    }
    out
}

fn st_edge_cut(simple: &Simple, net: &mut FlowNetwork, s: usize, t: usize) -> Vec<(u32, u32)> {
    net.reset();
    net.max_flow(s, t, None).expect("unit capacities");
    let side: Vec<bool> = net.reaching(&[t]).into_iter().map(|x| !x).collect();
    arcs_leaving(simple, &side)
}

/// Node cut on a prepared split network (`infinite` = arcs between nodes
/// of infinite capacity, the engine convention).
fn st_node_cut(
    simple: &Simple,
    net: &mut FlowNetwork,
    arcs: &[(u32, u32)],
    s: usize,
    t: usize,
    convention: CutConvention,
) -> Option<Vec<u32>> {
    let n = simple.n;
    if simple.has_arc(s, t) {
        return match convention {
            CutConvention::Engine => None,
            CutConvention::Networkx => Some(Vec::new()),
        };
    }
    net.reset();
    net.max_flow(n + s, t, None).expect("no adjacent terminals");
    let mut cut = match convention {
        CutConvention::Engine => {
            let side = net.reachable_from(&[n + s]);
            (0..n)
                .filter(|&v| side[v] && !side[n + v])
                .map(|v| v as u32)
                .collect::<Vec<u32>>()
        }
        CutConvention::Networkx => {
            let side: Vec<bool> = net.reaching(&[t]).into_iter().map(|x| !x).collect();
            let mut nodes = vec![false; n];
            for v in 0..n {
                if side[v] && !side[n + v] {
                    nodes[v] = true;
                }
            }
            for &(u, v) in arcs {
                let (u, v) = (u as usize, v as usize);
                if side[n + u] && !side[v] {
                    nodes[u] = true;
                    nodes[v] = true;
                }
            }
            nodes[s] = false;
            nodes[t] = false;
            (0..n).filter(|&v| nodes[v]).map(|v| v as u32).collect()
        }
    };
    cut.sort_unstable();
    Some(cut)
}

fn node_network(simple: &Simple, convention: CutConvention) -> (FlowNetwork, Vec<(u32, u32)>) {
    simple.node_network(match convention {
        CutConvention::Engine => f64::INFINITY,
        CutConvention::Networkx => 1.0,
    })
}

/// Minimum s-t edge cut. See the module documentation.
pub fn minimum_st_edge_cut(g: &Graph, s: usize, t: usize) -> Result<Vec<(u32, u32)>> {
    check_pair(s, t, g.node_count())?;
    let simple = Simple::new(g);
    let mut net = simple.edge_network();
    Ok(st_edge_cut(&simple, &mut net, s, t))
}

/// Minimum s-t node cut. See the module documentation.
pub fn minimum_st_node_cut(
    g: &Graph,
    s: usize,
    t: usize,
    convention: CutConvention,
) -> Result<Vec<u32>> {
    check_pair(s, t, g.node_count())?;
    let simple = Simple::new(g);
    let (mut net, arcs) = node_network(&simple, convention);
    st_node_cut(&simple, &mut net, &arcs, s, t, convention).ok_or_else(|| {
        Error::InvalidArgument(format!(
            "nodes {s} and {t} are adjacent, so no set of other nodes separates them"
        ))
    })
}

fn nonempty(g: &Graph) -> Result<()> {
    if g.node_count() == 0 {
        Err(Error::InvalidArgument(
            "cuts are not defined for a network with no nodes".into(),
        ))
    } else {
        Ok(())
    }
}

/// A smallest set of nodes whose removal disconnects the network. See the
/// module documentation.
pub fn minimum_node_cut(g: &Graph, convention: CutConvention, threads: usize) -> Result<Vec<u32>> {
    nonempty(g)?;
    let simple = Simple::new(g);
    if !simple.is_connected() {
        return Ok(Vec::new());
    }
    let (v, _, pairs) = even_pairs(&simple);
    // The isolating cut of v: its in- or out-neighbours, whichever are
    // fewer (in-neighbours on a tie, as NetworkX), or its neighbours.
    let mut best: Vec<u32> = if simple.directed && simple.out[v].len() < simple.inc[v].len() {
        simple.out[v].clone()
    } else {
        simple.inc[v].clone()
    };
    let (base, arcs) = node_network(&simple, convention);
    let parts = parallel::for_each_item(
        pairs.len(),
        threads,
        || (base.clone(), Vec::new()),
        |(net, found): &mut Found<u32>, i| {
            let (s, t) = pairs[i];
            if let Some(c) = st_node_cut(&simple, net, &arcs, s, t, convention) {
                found.push((i, c));
            }
        },
    );
    let mut all: Vec<(usize, Vec<u32>)> = parts.into_iter().flat_map(|p| p.1).collect();
    all.sort_by_key(|p| p.0);
    for (_, c) in all {
        if c.len() < best.len() {
            best = c;
        }
    }
    Ok(best)
}

/// A smallest set of edges whose removal disconnects the network. See the
/// module documentation.
pub fn minimum_edge_cut(g: &Graph, threads: usize) -> Result<Vec<(u32, u32)>> {
    nonempty(g)?;
    let simple = Simple::new(g);
    let n = simple.n;
    if !simple.is_connected() {
        return Ok(Vec::new());
    }
    // Start from the edges at a node of least degree.
    let degree = |v: usize| {
        if simple.directed {
            simple.out[v].len().min(simple.inc[v].len())
        } else {
            simple.out[v].len()
        }
    };
    let v = (0..n).min_by_key(|&v| degree(v)).unwrap_or(0);
    let mut best: Vec<(u32, u32)> = if simple.directed && simple.inc[v].len() < simple.out[v].len()
    {
        simple.inc[v].iter().map(|&u| (u, v as u32)).collect()
    } else {
        simple.out[v].iter().map(|&w| (v as u32, w)).collect()
    };
    if n == 1 {
        return Ok(best);
    }
    let pairs: Vec<(usize, usize)> = if simple.directed {
        (0..n).map(|i| (i, (i + 1) % n)).collect()
    } else {
        (1..n).map(|w| (0, w)).collect()
    };
    let base = simple.edge_network();
    let parts = parallel::for_each_item(
        pairs.len(),
        threads,
        || (base.clone(), Vec::new()),
        |(net, found): &mut Found<(u32, u32)>, i| {
            let (s, t) = pairs[i];
            found.push((i, st_edge_cut(&simple, net, s, t)));
        },
    );
    let mut all: Vec<(usize, Vec<(u32, u32)>)> = parts.into_iter().flat_map(|p| p.1).collect();
    all.sort_by_key(|p| p.0);
    for (_, c) in all {
        if c.len() < best.len() {
            best = c;
        }
    }
    best.sort_unstable();
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn cuts_of_two_triangles() {
        // Triangles 0-1-2 and 3-4-5 joined by 2-3.
        let g = graph(
            6,
            &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)],
            false,
        );
        assert_eq!(minimum_st_edge_cut(&g, 0, 5).unwrap(), vec![(2, 3)]);
        assert_eq!(minimum_edge_cut(&g, 2).unwrap(), vec![(2, 3)]);
        // The engine takes the cut nearest the source; NetworkX's cut
        // ({3}, checked with NetworkX 3.7) the one nearest the target.
        assert_eq!(
            minimum_st_node_cut(&g, 0, 5, CutConvention::Engine).unwrap(),
            vec![2]
        );
        assert_eq!(
            minimum_st_node_cut(&g, 0, 5, CutConvention::Networkx).unwrap(),
            vec![3]
        );
        for c in [CutConvention::Engine, CutConvention::Networkx] {
            assert_eq!(minimum_node_cut(&g, c, 1).unwrap().len(), 1);
        }
        assert!(minimum_st_node_cut(&g, 0, 1, CutConvention::Engine).is_err());
        assert!(minimum_st_node_cut(&g, 0, 1, CutConvention::Networkx)
            .unwrap()
            .is_empty());
        let d = graph(3, &[(0, 1)], false);
        assert!(minimum_node_cut(&d, CutConvention::Engine, 1)
            .unwrap()
            .is_empty());
    }
}
