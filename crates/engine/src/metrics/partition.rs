//! Community partitions: Girvan-Newman divisive clustering, Newman's
//! leading eigenvector method, and Clauset-Newman-Moore greedy
//! agglomeration ("fast greedy").
//!
//! Each method lives in its own submodule, which cites its sources and
//! states its conventions. Shared conventions:
//! - Undirected networks only, as igraph's `cluster_edge_betweenness()`
//!   (for modularity), `cluster_leading_eigen()`, and
//!   `cluster_fast_greedy()` require; directed networks are rejected rather
//!   than silently symmetrized.
//! - Modularity is that of [`crate::metrics::community::modularity`].
//! - Groups are labelled 0, 1, ... in order of their first node.

pub mod compare;
pub mod dendrogram;
pub mod fast_greedy;
pub mod girvan_newman;
pub mod infomap;
pub mod label_propagation;
pub mod leading_eigen;
pub mod leiden;
pub mod walktrap;

pub use compare::{compare_partitions, split_join_distance, CompareMethod};
pub use dendrogram::{
    check_merges, cut_at_groups, cut_at_max_modularity, cut_dendrogram, dendrogram_modularity,
    merges_from_edge_removals, merges_from_node_pairs,
};
pub use fast_greedy::{fast_greedy, FastGreedy};
pub use girvan_newman::{girvan_newman, GirvanNewman};
pub use infomap::{infomap, map_equation, Infomap};
pub use label_propagation::{label_propagation, LabelPropagationOptions, LpaMode, LpaVariant};
pub use leading_eigen::{leading_eigen, LeadingEigen};
pub use leiden::{leiden, Leiden, LeidenObjective, LeidenOptions};
pub use walktrap::{walktrap, Walktrap};

use crate::error::{Error, Result};
use crate::graph::Graph;

pub(crate) fn require_undirected(g: &Graph, method: &str) -> Result<()> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(format!(
            "{method} is defined here for undirected networks only; build the network \
             as undirected first"
        )));
    }
    Ok(())
}

/// Relabels group labels to 0..k in order of first appearance.
pub(crate) fn relabel(labels: &[usize]) -> Vec<usize> {
    crate::metrics::community::dense_labels(labels).0
}

/// Weights for additive use (modularity), checked to be non-negative.
pub(crate) fn additive_weights(g: &Graph, weighted: bool) -> Result<Option<&[f64]>> {
    crate::weights::additive_weights(g, weighted, "modularity")
}

pub(crate) fn as_labels(m: &[usize]) -> Vec<i64> {
    m.iter().map(|&c| c as i64).collect()
}
