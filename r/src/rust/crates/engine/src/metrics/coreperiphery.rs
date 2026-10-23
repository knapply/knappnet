//! Network topography beyond components: core-periphery structure,
//! compactness and distance-weighted fragmentation, and the small-world
//! quotient.
//!
//! Each measure lives in its own submodule, which cites its sources and
//! states its conventions.

pub mod compactness;
pub mod continuous;
pub mod discrete;
pub mod small_world;

pub use compactness::compactness;
pub use continuous::{
    continuous_core_periphery, ContinuousCorePeriphery, ContinuousFit, ContinuousOptions,
};
pub use discrete::{core_periphery, CoreLinks, CorePeriphery, CorePeripheryOptions};
pub use small_world::{small_world, SmallWorld, SmallWorldOptions};

use crate::error::Result;
use crate::graph::Graph;
use crate::weights::additive_weights;

/// The adjacency matrix as the core-periphery fits read it: for every node
/// v, its other endpoints u with s1 = the sum of the values of the pairs
/// between v and u and s2 = the sum of their squares. Undirected networks
/// have one pair per unordered pair of nodes (value a_vu); directed
/// networks two, (v, u) and (u, v). A value is 1 for a tie (binary) or the
/// total weight of the edges from one node to the other (weighted).
/// Self-loops are dropped.
pub(crate) struct PairRows {
    pub rows: Vec<Vec<(u32, f64, f64)>>,
    /// 1 for undirected (unordered pairs), 2 for directed (ordered pairs).
    pub per_node_pair: usize,
    /// Number of pairs, sum of values, and sum of squared values.
    pub pairs: f64,
    pub sum: f64,
    pub sum_sq: f64,
}

impl PairRows {
    pub fn new(g: &Graph, weighted: bool) -> Result<PairRows> {
        let n = g.node_count();
        let w = additive_weights(g, weighted, "a core-periphery fit")?;
        let directed = g.is_directed();
        // Value of each ordered (directed) or unordered pair.
        let mut values: Vec<(u32, u32, f64)> = Vec::with_capacity(g.edge_count());
        let (src, tgt) = (g.edge_sources(), g.edge_targets());
        for e in 0..g.edge_count() {
            let (a, b) = (src[e], tgt[e]);
            if a == b {
                continue;
            }
            let x = w.map_or(1.0, |w| w[e]);
            let (a, b) = if directed || a < b { (a, b) } else { (b, a) };
            values.push((a, b, x));
        }
        values.sort_unstable_by_key(|p| (p.0, p.1));
        let mut merged: Vec<(u32, u32, f64)> = Vec::with_capacity(values.len());
        for (a, b, x) in values {
            match merged.last_mut() {
                Some(last) if last.0 == a && last.1 == b => {
                    if weighted {
                        last.2 += x;
                    }
                }
                _ => merged.push((a, b, x)),
            }
        }
        let mut rows: Vec<Vec<(u32, f64, f64)>> = vec![Vec::new(); n];
        let (mut sum, mut sum_sq) = (0.0, 0.0);
        for &(a, b, x) in &merged {
            sum += x;
            sum_sq += x * x;
            rows[a as usize].push((b, x, x * x));
            rows[b as usize].push((a, x, x * x));
        }
        // Directed pairs (a, b) and (b, a) are both listed at each end;
        // combine them into one entry per neighbor.
        for row in rows.iter_mut() {
            row.sort_unstable_by_key(|p| p.0);
            let mut out: Vec<(u32, f64, f64)> = Vec::with_capacity(row.len());
            for &(u, s1, s2) in row.iter() {
                match out.last_mut() {
                    Some(last) if last.0 == u => {
                        last.1 += s1;
                        last.2 += s2;
                    }
                    _ => out.push((u, s1, s2)),
                }
            }
            *row = out;
        }
        let per = if directed { 2 } else { 1 };
        let nf = n as f64;
        Ok(PairRows {
            rows,
            per_node_pair: per,
            pairs: per as f64 * nf * (nf - 1.0) / 2.0,
            sum,
            sum_sq,
        })
    }
}
