//! Cohesive subsets and factions: bicomponents, n-cliques and n-clans,
//! k-plexes, lambda sets, and UCINET-style factions.
//!
//! Each method lives in its own submodule, which cites its sources and
//! states its conventions. They share these conventions:
//! - Except for factions of directed networks, the network is read as an
//!   undirected simple graph, as [`crate::metrics::cliques`] does: an arc
//!   in either direction joins two nodes, parallel edges count once, and
//!   self-loops are ignored.
//! - Weights are ignored: these are definitions on binary ties.
//! - Node sets list their members in increasing order and are returned in
//!   lexicographic order, so results never depend on the thread count.

pub mod bicomponents;
pub mod kgroups;
pub mod kplex;
pub mod lambda;
pub mod ncliques;

pub use bicomponents::bicomponents;
pub use kgroups::{factions, FactionFit, FactionOptions, Factions};
pub use kplex::{k_plexes, KPlexOptions};
pub use lambda::{edge_connectivity_tree, lambda_connectivity, FlowTree};
pub use ncliques::{n_clans, n_cliques};

use crate::graph::Graph;

/// Rows of the undirected simple graph underlying `g`: for every node, its
/// distinct neighbors other than itself, in increasing order.
pub(crate) fn simple_rows(g: &Graph) -> Vec<Vec<u32>> {
    let s = g.simple_undirected(false);
    (0..g.node_count())
        .map(|v| s.neighbors(v).to_vec())
        .collect()
}

/// Sorts each set and the list of sets, for output.
pub(crate) fn canonical_sets(mut sets: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    for s in sets.iter_mut() {
        s.sort_unstable();
    }
    sets.sort();
    sets
}

/// Whether `rows` (sorted) joins `a` and `b`.
#[inline]
pub(crate) fn adjacent(rows: &[Vec<u32>], a: u32, b: u32) -> bool {
    rows[a as usize].binary_search(&b).is_ok()
}
