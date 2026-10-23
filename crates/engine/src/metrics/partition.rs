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

mod eigen;
pub mod fast_greedy;
pub mod girvan_newman;
pub mod leading_eigen;

pub(crate) use eigen::symmetric_eigen;
pub use fast_greedy::{fast_greedy, FastGreedy};
pub use girvan_newman::{girvan_newman, GirvanNewman};
pub use leading_eigen::{leading_eigen, LeadingEigen};

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
    let mut map = std::collections::HashMap::new();
    labels
        .iter()
        .map(|&c| {
            let next = map.len();
            *map.entry(c).or_insert(next)
        })
        .collect()
}

/// Weights for additive use (modularity), checked to be non-negative.
pub(crate) fn additive_weights(g: &Graph, weighted: bool) -> Result<Option<&[f64]>> {
    if !weighted {
        return Ok(None);
    }
    let w = g.edge_weights().ok_or(Error::MissingWeights)?;
    if let Some(edge) = w.iter().position(|&x| x < 0.0) {
        return Err(Error::InvalidArgument(format!(
            "edge {edge} has negative weight {}; modularity requires non-negative weights",
            w[edge]
        )));
    }
    Ok(Some(w))
}

pub(crate) fn as_labels(m: &[usize]) -> Vec<i64> {
    m.iter().map(|&c| c as i64).collect()
}
