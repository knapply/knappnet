//! Connectivity: local and global node and edge connectivity, minimum
//! node and edge cuts, disjoint paths, k-components, k-edge-components,
//! local bridges, and the component structure of directed networks
//! (condensation, attracting components, semiconnectedness).
//!
//! Each submodule cites its sources and states its conventions. They share
//! these:
//! - Connectivity is a property of the simple graph: parallel edges count
//!   once and self-loops are ignored, as in NetworkX
//!   (`hagberg2008exploring`), whose `Graph` and `DiGraph` cannot hold
//!   parallel edges and whose connectivity functions build auxiliary
//!   networks that ignore loops. Weights are not used.
//! - Directed networks are read with their directions (NetworkX's
//!   behaviour for a `DiGraph`); k-components and local bridges are
//!   defined for undirected networks only.
//! - Every flow is computed by the engine's one maximum-flow
//!   implementation ([`crate::flow`]); node connectivity uses the split
//!   network of Even and Tarjan (1975, `even1975network`), in which node v
//!   becomes an arc v_in -> v_out of capacity 1 and each arc u -> v the arc
//!   u_out -> v_in, also of capacity 1 (NetworkX's auxiliary digraph).
//! - Results never depend on the thread count: parallel loops compute
//!   independent flows and combine integer results by minimum or in a
//!   fixed order.
//!
//! Bridges, cut points, and biconnected components (blocks) were already
//! provided by [`crate::metrics::descriptive::cut_structure`] and
//! [`crate::metrics::factions::bicomponents()`]; this module reuses them.

pub mod components;
pub mod cuts;
pub mod disjoint;
pub mod kcomponents;
pub mod kedge;
pub mod values;

pub use components::{
    attracting_components, condensation, is_semiconnected, local_bridges, strong_components,
    LocalBridges,
};
pub use cuts::{
    minimum_edge_cut, minimum_node_cut, minimum_st_edge_cut, minimum_st_node_cut, CutConvention,
};
pub use disjoint::{edge_disjoint_paths, node_disjoint_paths};
pub use kcomponents::k_components;
pub use kedge::k_edge_components;
pub use values::{
    all_pairs_node_connectivity, average_node_connectivity, edge_connectivity,
    local_edge_connectivity, local_node_connectivity, node_connectivity,
};

use crate::error::{Error, Result};
use crate::flow::FlowNetwork;
use crate::graph::Graph;

/// The simple graph underlying `g` (directions kept): for every node its
/// distinct out-neighbours and in-neighbours other than itself, sorted.
/// For an undirected network both lists are the neighbours.
#[derive(Debug, Clone)]
pub(crate) struct Simple {
    pub n: usize,
    pub directed: bool,
    pub out: Vec<Vec<u32>>,
    pub inc: Vec<Vec<u32>>,
}

impl Simple {
    pub fn new(g: &Graph) -> Simple {
        let n = g.node_count();
        let row = |csr: &crate::graph::Csr, v: usize| -> Vec<u32> {
            let mut r: Vec<u32> = csr
                .neighbors(v)
                .iter()
                .copied()
                .filter(|&w| w as usize != v)
                .collect();
            r.dedup();
            r
        };
        let out: Vec<Vec<u32>> = (0..n).map(|v| row(g.out(), v)).collect();
        let inc = if g.is_directed() {
            (0..n).map(|v| row(g.inc(), v)).collect()
        } else {
            out.clone()
        };
        Simple {
            n,
            directed: g.is_directed(),
            out,
            inc,
        }
    }

    /// The simple graph induced by `nodes` (sorted), relabelled 0..k in
    /// that order.
    pub fn induced(&self, nodes: &[u32]) -> Simple {
        let mut index = vec![u32::MAX; self.n];
        for (i, &v) in nodes.iter().enumerate() {
            index[v as usize] = i as u32;
        }
        let pick = |rows: &Vec<Vec<u32>>| -> Vec<Vec<u32>> {
            nodes
                .iter()
                .map(|&v| {
                    rows[v as usize]
                        .iter()
                        .filter_map(|&w| {
                            let i = index[w as usize];
                            (i != u32::MAX).then_some(i)
                        })
                        .collect()
                })
                .collect()
        };
        let out = pick(&self.out);
        let inc = if self.directed {
            pick(&self.inc)
        } else {
            out.clone()
        };
        Simple {
            n: nodes.len(),
            directed: self.directed,
            out,
            inc,
        }
    }

    /// Whether the arc u -> v (edge {u, v} if undirected) exists.
    pub fn has_arc(&self, u: usize, v: usize) -> bool {
        self.out[u].binary_search(&(v as u32)).is_ok()
    }

    /// Distinct neighbours of `v` in either direction, sorted.
    pub fn neighbors_any(&self, v: usize) -> Vec<u32> {
        if !self.directed {
            return self.out[v].clone();
        }
        let mut r: Vec<u32> = self.out[v].iter().chain(&self.inc[v]).copied().collect();
        r.sort_unstable();
        r.dedup();
        r
    }

    /// Whether the network is connected (weakly, if directed). A network
    /// with no nodes is not.
    pub fn is_connected(&self) -> bool {
        if self.n == 0 {
            return false;
        }
        let mut seen = vec![false; self.n];
        seen[0] = true;
        let mut stack = vec![0usize];
        let mut count = 1;
        while let Some(v) = stack.pop() {
            for &w in self.out[v].iter().chain(&self.inc[v]) {
                let w = w as usize;
                if !seen[w] {
                    seen[w] = true;
                    count += 1;
                    stack.push(w);
                }
            }
        }
        count == self.n
    }

    /// Arcs (u, v) of the simple graph: every ordered pair with an arc if
    /// directed, every pair u < v with an edge if undirected; sorted.
    pub fn arcs(&self) -> Vec<(u32, u32)> {
        let mut a = Vec::new();
        for u in 0..self.n {
            for &v in &self.out[u] {
                if self.directed || (u as u32) < v {
                    a.push((u as u32, v));
                }
            }
        }
        a
    }

    /// The edge network: one arc of capacity 1 per arc, or one undirected
    /// edge of capacity 1 per edge.
    pub fn edge_network(&self) -> FlowNetwork {
        let mut net = FlowNetwork::new(self.n);
        for (u, v) in self.arcs() {
            if self.directed {
                net.add_arc(u as usize, v as usize, 1.0);
            } else {
                net.add_edge(u as usize, v as usize, 1.0);
            }
        }
        net
    }

    /// The split network for node connectivity: v_in = v, v_out = n + v,
    /// pair v the arc v_in -> v_out (capacity 1), then one arc
    /// u_out -> v_in per arc u -> v (both ways if undirected) of capacity
    /// `edge_cap`. Pair n + i is the i-th of those arcs; its endpoints are
    /// in `split_arcs`.
    pub fn node_network(&self, edge_cap: f64) -> (FlowNetwork, Vec<(u32, u32)>) {
        let n = self.n;
        let mut net = FlowNetwork::new(2 * n);
        for v in 0..n {
            net.add_arc(v, n + v, 1.0);
        }
        let mut arcs = Vec::new();
        for u in 0..n {
            for &v in &self.out[u] {
                arcs.push((u as u32, v));
            }
        }
        for &(u, v) in &arcs {
            net.add_arc(n + u as usize, v as usize, edge_cap);
        }
        (net, arcs)
    }
}

/// Checks that `s` and `t` are distinct nodes of a graph with `n` nodes.
pub(crate) fn check_pair(s: usize, t: usize, n: usize) -> Result<()> {
    for (v, what) in [(s, "source"), (t, "target")] {
        if v >= n {
            return Err(Error::InvalidArgument(format!(
                "{what} {v} is not a node (n = {n})"
            )));
        }
    }
    if s == t {
        return Err(Error::InvalidArgument(
            "source and target must be different nodes".into(),
        ));
    }
    Ok(())
}
