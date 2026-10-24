//! Network flows: maximum flows, minimum cuts, minimum-cost flows, and
//! Gomory-Hu trees.
//!
//! One residual network ([`network::FlowNetwork`], Dinic's algorithm)
//! serves every flow and cut computation in the engine's newer modules:
//! the graph-level functions in [`maxflow`] build it with edge
//! capacities, super-sources and super-sinks for several terminals, or
//! split nodes for node capacities; Gomory-Hu trees ([`gomory_hu`]) and
//! connectivity, cuts, disjoint paths, and k-components
//! ([`crate::metrics::connectivity`]) reuse it. Minimum-cost flows use the network simplex
//! method ([`mincost`]).
//!
//! Each module cites its sources and states its conventions; all of them
//! accept infinite capacities and refuse unbounded problems with an error.

pub mod gomory_hu;
pub mod maxflow;
pub mod mincost;
pub mod network;

pub use gomory_hu::{gomory_hu_tree, GomoryHuTree};
pub use maxflow::{
    edge_network, maximum_flow, minimum_cut, minimum_node_cut_capacitated, split_network, CutSide,
    MaxFlow, MinCut, NodeCut,
};
pub use mincost::{max_flow_min_cost, min_cost_flow, CostFlow};
pub use network::FlowNetwork;
