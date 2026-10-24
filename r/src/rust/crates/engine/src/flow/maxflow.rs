//! Maximum flows and minimum cuts of a network, with edge capacities,
//! several sources and sinks, and node capacities.
//!
//! Sources: Ford and Fulkerson (1956, `ford1956maximal`) for the problem
//! and the max-flow min-cut theorem; Dinic (1970, `dinic1970algorithm`)
//! for the algorithm ([`crate::flow::network`]). Several sources or sinks
//! are joined to a super-source or super-sink by arcs of infinite
//! capacity, and node capacities are modelled by splitting every node v
//! into an in-copy and an out-copy joined by an arc of v's capacity (Ford
//! and Fulkerson 1962, `ford1962flows`, section I.11), so one algorithm
//! serves every variant.
//!
//! Conventions:
//! - Capacities: one per edge, in input order, non-negative and possibly
//!   infinite; without capacities every edge has capacity 1. NetworkX
//!   (`hagberg2008exploring`) reads a missing `capacity` attribute as
//!   infinite; a binding reproduces that by passing `inf`.
//! - Self-loops carry no flow and never belong to a cut (NetworkX drops
//!   them from its residual network).
//! - Parallel edges are separate edges whose capacities add up. NetworkX
//!   does not accept multigraphs in its flow functions.
//! - Directed edges carry flow from source to target; an undirected edge
//!   carries flow either way up to its capacity, and its flow is reported
//!   as the net flow from its first to its second endpoint (negative when
//!   the flow runs the other way). NetworkX's flow dictionary for an
//!   undirected graph holds the same net flow under the direction it
//!   runs.
//! - A maximum flow is generally not unique; the flow on each edge is the
//!   one Dinic's algorithm finds, which can differ from NetworkX's while
//!   having the same value.
//! - Minimum cuts are generally not unique either, but two are canonical
//!   and do not depend on which maximum flow was found: the nodes
//!   reachable from the sources in the residual network form the smallest
//!   source side ([`CutSide::Smallest`]), and the complement of the nodes
//!   that reach a sink in it forms the largest ([`CutSide::Largest`]).
//!   NetworkX's `minimum_cut` returns the largest source side.
//! - A maximum flow that can be increased without bound (a path of
//!   infinite capacity from a source to a sink) is an error, as is a node
//!   that is both a source and a sink (NetworkX raises
//!   `NetworkXUnbounded` and `NetworkXError` for these).

use super::network::FlowNetwork;
use crate::error::{Error, Result};
use crate::graph::Graph;

/// Which canonical minimum cut to report. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CutSide {
    /// The nodes reachable from the sources in the residual network: the
    /// smallest source side of any minimum cut.
    #[default]
    Smallest,
    /// Every node that cannot reach a sink in the residual network: the
    /// largest source side of any minimum cut (NetworkX's
    /// `minimum_cut`).
    Largest,
}

impl CutSide {
    /// Decodes the C interface codes: 0 = smallest source side,
    /// 1 = largest source side (NetworkX).
    pub fn from_code(code: i32) -> Result<CutSide> {
        match code {
            0 => Ok(CutSide::Smallest),
            1 => Ok(CutSide::Largest),
            _ => Err(Error::InvalidArgument(format!(
                "unknown cut side code {code}"
            ))),
        }
    }
}

/// Checks a capacity vector of `len` values: non-negative, possibly
/// infinite, never NaN.
pub(crate) fn check_capacities(cap: &[f64], len: usize, what: &'static str) -> Result<()> {
    if cap.len() != len {
        return Err(Error::LengthMismatch {
            what,
            expected: len,
            found: cap.len(),
        });
    }
    if let Some(i) = cap.iter().position(|c| c.is_nan() || *c < 0.0) {
        return Err(Error::InvalidArgument(format!(
            "{what}: value {i} is {}; capacities must be non-negative numbers (inf allowed)",
            cap[i]
        )));
    }
    Ok(())
}

/// The residual network of `g`: pair `e` is edge `e` (a self-loop gets
/// capacity 0). Capacities as in the module documentation.
pub fn edge_network(g: &Graph, capacity: Option<&[f64]>) -> Result<FlowNetwork> {
    let m = g.edge_count();
    if let Some(c) = capacity {
        check_capacities(c, m, "capacities")?;
    }
    let mut net = FlowNetwork::new(g.node_count());
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    for e in 0..m {
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        let c = if u == v {
            0.0
        } else {
            capacity.map_or(1.0, |c| c[e])
        };
        if g.is_directed() {
            net.add_arc(u, v, c);
        } else {
            net.add_edge(u, v, c);
        }
    }
    Ok(net)
}

/// Joins several terminals to a new node by infinite arcs (out of it for
/// sources, into it for sinks); one terminal is used as it is.
fn terminal(net: &mut FlowNetwork, nodes: &[usize], sources: bool) -> usize {
    if nodes.len() == 1 {
        return nodes[0];
    }
    let x = net.add_node();
    for &v in nodes {
        if sources {
            net.add_arc(x, v, f64::INFINITY);
        } else {
            net.add_arc(v, x, f64::INFINITY);
        }
    }
    x
}

fn check_nodes(nodes: &[usize], n: usize, what: &str) -> Result<()> {
    if nodes.is_empty() {
        return Err(Error::InvalidArgument(format!("no {what} given")));
    }
    if let Some(&v) = nodes.iter().find(|&&v| v >= n) {
        return Err(Error::InvalidArgument(format!(
            "{what} {v} is not a node (n = {n})"
        )));
    }
    let mut seen = vec![false; n];
    for &v in nodes {
        if seen[v] {
            return Err(Error::InvalidArgument(format!("{what} {v} is repeated")));
        }
        seen[v] = true;
    }
    Ok(())
}

/// A maximum flow.
#[derive(Debug, Clone, PartialEq)]
pub struct MaxFlow {
    /// The flow value: the net flow out of the sources.
    pub value: f64,
    /// Flow on every edge in input order (net flow from first to second
    /// endpoint for undirected edges; 0 for self-loops).
    pub edge_flow: Vec<f64>,
}

/// A minimum cut.
#[derive(Debug, Clone, PartialEq)]
pub struct MinCut {
    /// Its capacity, equal to the maximum flow value.
    pub value: f64,
    /// Whether each node is on the source side.
    pub source_side: Vec<bool>,
    /// Edges of positive capacity from the source side to the other side
    /// (either way for undirected edges), in increasing order.
    pub edges: Vec<usize>,
}

/// Solves a maximum flow on the edge network of `g`; returns the network
/// and its terminal nodes.
fn solve_edges(
    g: &Graph,
    sources: &[usize],
    sinks: &[usize],
    capacity: Option<&[f64]>,
) -> Result<(FlowNetwork, usize, usize, f64)> {
    let n = g.node_count();
    check_nodes(sources, n, "source")?;
    check_nodes(sinks, n, "sink")?;
    let mut net = edge_network(g, capacity)?;
    let s = terminal(&mut net, sources, true);
    let t = terminal(&mut net, sinks, false);
    if let Some(&v) = sources.iter().find(|v| sinks.contains(v)) {
        return Err(Error::InvalidArgument(format!(
            "node {v} is both a source and a sink"
        )));
    }
    let value = net.max_flow(s, t, None)?;
    Ok((net, s, t, value))
}

/// Maximum flow from `sources` to `sinks`. See the module documentation.
pub fn maximum_flow(
    g: &Graph,
    sources: &[usize],
    sinks: &[usize],
    capacity: Option<&[f64]>,
) -> Result<MaxFlow> {
    let (net, _, _, value) = solve_edges(g, sources, sinks, capacity)?;
    let edge_flow = (0..g.edge_count()).map(|e| net.flow(e)).collect();
    Ok(MaxFlow { value, edge_flow })
}

/// Minimum cut between `sources` and `sinks`. See the module
/// documentation.
pub fn minimum_cut(
    g: &Graph,
    sources: &[usize],
    sinks: &[usize],
    capacity: Option<&[f64]>,
    side: CutSide,
) -> Result<MinCut> {
    let n = g.node_count();
    let (net, s, t, value) = solve_edges(g, sources, sinks, capacity)?;
    let mut source_side = match side {
        CutSide::Smallest => net.reachable_from(&[s]),
        CutSide::Largest => net.reaching(&[t]).into_iter().map(|x| !x).collect(),
    };
    source_side.truncate(n);
    let edges = crossing_edges(g, &source_side, |e| net.capacity(e) > 0.0);
    Ok(MinCut {
        value,
        source_side,
        edges,
    })
}

/// Edges from the source side to the other (either way if undirected)
/// that `keep` accepts, in increasing order.
fn crossing_edges(g: &Graph, side: &[bool], keep: impl Fn(usize) -> bool) -> Vec<usize> {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    (0..g.edge_count())
        .filter(|&e| {
            let (a, b) = (side[src[e] as usize], side[tgt[e] as usize]);
            (if g.is_directed() { a && !b } else { a != b }) && keep(e)
        })
        .collect()
}

/// A minimum cut of a network with node capacities.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeCut {
    /// Its capacity, equal to the maximum flow value.
    pub value: f64,
    /// Nodes in the cut, in increasing order.
    pub nodes: Vec<usize>,
    /// Edges in the cut, in increasing order.
    pub edges: Vec<usize>,
}

/// Split network of `g`: node v's in-copy is v, its out-copy n + v, and
/// pair v is the arc between them; pair n + e is edge e from its first
/// endpoint's out-copy to its second endpoint's in-copy, and for an
/// undirected graph pair n + m + e is the opposite arc.
pub fn split_network(
    g: &Graph,
    node_capacity: Option<&[f64]>,
    edge_capacity: Option<&[f64]>,
) -> Result<FlowNetwork> {
    let n = g.node_count();
    let m = g.edge_count();
    if let Some(c) = node_capacity {
        check_capacities(c, n, "node capacities")?;
    }
    if let Some(c) = edge_capacity {
        check_capacities(c, m, "edge capacities")?;
    }
    let mut net = FlowNetwork::new(2 * n);
    for v in 0..n {
        net.add_arc(v, n + v, node_capacity.map_or(1.0, |c| c[v]));
    }
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let cap = |e: usize| {
        if src[e] == tgt[e] {
            0.0
        } else {
            edge_capacity.map_or(f64::INFINITY, |c| c[e])
        }
    };
    for e in 0..m {
        net.add_arc(n + src[e] as usize, tgt[e] as usize, cap(e));
    }
    if !g.is_directed() {
        for e in 0..m {
            net.add_arc(n + tgt[e] as usize, src[e] as usize, cap(e));
        }
    }
    Ok(net)
}

/// Minimum cut separating `sources` from `sinks` when nodes have
/// capacities (`node_capacity`, default 1) as well as edges
/// (`edge_capacity`, default infinite, so that only nodes are cut). The
/// terminals themselves cannot be cut. `side` picks the canonical cut as
/// for [`minimum_cut`], on the split network.
pub fn minimum_node_cut_capacitated(
    g: &Graph,
    sources: &[usize],
    sinks: &[usize],
    node_capacity: Option<&[f64]>,
    edge_capacity: Option<&[f64]>,
    side: CutSide,
) -> Result<NodeCut> {
    let n = g.node_count();
    let m = g.edge_count();
    check_nodes(sources, n, "source")?;
    check_nodes(sinks, n, "sink")?;
    if let Some(&v) = sources.iter().find(|v| sinks.contains(v)) {
        return Err(Error::InvalidArgument(format!(
            "node {v} is both a source and a sink"
        )));
    }
    let mut net = split_network(g, node_capacity, edge_capacity)?;
    for &v in sources.iter().chain(sinks) {
        net.set_capacity(v, f64::INFINITY, false);
    }
    let outs: Vec<usize> = sources.iter().map(|&s| n + s).collect();
    let s = terminal(&mut net, &outs, true);
    let t = terminal(&mut net, sinks, false);
    let value = net.max_flow(s, t, None).map_err(|e| match e {
        Error::InvalidArgument(msg) if msg.contains("unbounded") => Error::InvalidArgument(
            "no node cut separates the sources from the sinks: a path of infinite capacity \
             (for example an edge of infinite capacity between a source and a sink) joins them"
                .into(),
        ),
        other => other,
    })?;
    let ss = match side {
        CutSide::Smallest => net.reachable_from(&[s]),
        CutSide::Largest => net.reaching(&[t]).into_iter().map(|x| !x).collect(),
    };
    let nodes = (0..n)
        .filter(|&v| ss[v] && !ss[n + v] && net.capacity(v) > 0.0)
        .collect();
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let edges = (0..m)
        .filter(|&e| {
            let (a, b) = (src[e] as usize, tgt[e] as usize);
            let fwd = ss[n + a] && !ss[b];
            let back = !g.is_directed() && ss[n + b] && !ss[a];
            (fwd || back) && net.capacity(n + e) > 0.0
        })
        .collect();
    Ok(NodeCut {
        value,
        nodes,
        edges,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn flows_and_cuts() {
        // Two parallel paths 0-1-3 and 0-2-3 plus a bridge 1-2.
        let g = graph(4, &[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)], true);
        let cap = [3.0, 2.0, 1.0, 2.0, 3.0];
        let f = maximum_flow(&g, &[0], &[3], Some(&cap)).unwrap();
        assert_eq!(f.value, 5.0);
        // Conservation at the inner nodes.
        assert_eq!(f.edge_flow[0], f.edge_flow[2] + f.edge_flow[3]);
        let small = minimum_cut(&g, &[0], &[3], Some(&cap), CutSide::Smallest).unwrap();
        let large = minimum_cut(&g, &[0], &[3], Some(&cap), CutSide::Largest).unwrap();
        assert_eq!(small.value, 5.0);
        assert_eq!(small.source_side, vec![true, false, false, false]);
        assert_eq!(small.edges, vec![0, 1]);
        assert_eq!(large.source_side, vec![true, true, true, false]);
        assert_eq!(large.edges, vec![3, 4]);
        // Several sources and sinks.
        let f = maximum_flow(&g, &[0, 1], &[3], Some(&cap)).unwrap();
        assert_eq!(f.value, 5.0);
        assert!(maximum_flow(&g, &[0], &[0], None).is_err());
    }

    #[test]
    fn node_cuts() {
        // A 4-cycle 0-1-2-3-0 with a chord 1-3: separating 0 from 2 needs
        // nodes 1 and 3.
        let g = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0), (1, 3)], false);
        let c =
            minimum_node_cut_capacitated(&g, &[0], &[2], None, None, CutSide::Smallest).unwrap();
        assert_eq!(c.value, 2.0);
        assert_eq!(c.nodes, vec![1, 3]);
        assert!(c.edges.is_empty());
        // Node 1 is expensive: cut edge 0-1 (capacity 1) instead.
        let c = minimum_node_cut_capacitated(
            &g,
            &[0],
            &[2],
            Some(&[1.0, 5.0, 1.0, 1.0]),
            Some(&[1.0, 9.0, 9.0, 9.0, 9.0]),
            CutSide::Smallest,
        )
        .unwrap();
        assert_eq!(c.value, 2.0);
        assert_eq!(c.nodes, vec![3]);
        assert_eq!(c.edges, vec![0]);
        // Adjacent terminals cannot be separated by nodes alone.
        assert!(
            minimum_node_cut_capacitated(&g, &[0], &[1], None, None, CutSide::Smallest).is_err()
        );
    }
}
