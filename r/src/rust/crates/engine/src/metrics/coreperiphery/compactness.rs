//! Compactness and distance-weighted fragmentation.
//!
//! Borgatti (2006, `borgatti2006identifying`) measures how far a network
//! is from being fragmented with the distance-weighted fragmentation
//!
//! ```text
//! DF = 1 - sum_{i != j} (1 / d_ij) / (n (n - 1)),
//! ```
//!
//! where 1 / d_ij = 0 for nodes that cannot reach each other. Its
//! complement, the mean reciprocal distance over ordered pairs of distinct
//! nodes, is the network's compactness (also called global efficiency;
//! Latora and Marchiori 2001, `latora2001efficient`). DF is 0 for a
//! complete network and 1 for a network without ties; unlike plain
//! fragmentation it distinguishes networks that are connected but
//! stretched out.
//!
//! Conventions:
//! - Directed networks use directed distances (d_ij from i to j) over
//!   ordered pairs; `mode` can instead use paths to each node or ignore
//!   direction. Undirected networks count each pair twice, which does not
//!   change the mean. This is NetworkX's `global_efficiency` for
//!   undirected networks.
//! - Distances follow [`crate::metrics::descriptive::distance_rows`]:
//!   without weights they count steps; with weights they are sums of
//!   lengths (weights read as distances, or as strengths converted to
//!   lengths 1/w), which must be strictly positive. Self-loops never
//!   shorten a path and parallel edges count through the shortest.
//! - Networks with fewer than two nodes have no pairs, and the result is
//!   NaN.

use crate::error::Result;
use crate::graph::Graph;
use crate::metrics::descriptive::{distance_rows, PathOptions};

/// Compactness: the mean of 1 / d_ij over ordered pairs of distinct nodes
/// (1 minus Borgatti's distance-weighted fragmentation). See the module
/// documentation.
pub fn compactness(g: &Graph, opts: PathOptions) -> Result<f64> {
    let n = g.node_count();
    if n < 2 {
        return Ok(f64::NAN);
    }
    let mut total = 0.0f64;
    distance_rows(g, None, opts, |k, row| {
        for (j, &d) in row.iter().enumerate() {
            if j != k && d.is_finite() && d > 0.0 {
                total += 1.0 / d;
            }
        }
    })?;
    Ok(total / (n as f64 * (n as f64 - 1.0)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::degree::DegreeMode;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn known_values() {
        // Path a-b-c: distances 1, 1, 2 each way: (1 + 1 + 1/2) * 2 / 6.
        let p = und(3, &[(0, 1), (1, 2)]);
        let c = compactness(&p, PathOptions::default()).unwrap();
        assert!((c - 5.0 / 6.0).abs() < 1e-15);
        let k = und(3, &[(0, 1), (1, 2), (0, 2)]);
        assert_eq!(compactness(&k, PathOptions::default()).unwrap(), 1.0);
        assert_eq!(
            compactness(&und(3, &[]), PathOptions::default()).unwrap(),
            0.0
        );
        assert!(compactness(&und(1, &[]), PathOptions::default())
            .unwrap()
            .is_nan());
        // Directed 0 -> 1 -> 2: reachable pairs (0,1), (1,2), (0,2).
        let d = Graph::from_edges(3, &[0, 1], &[1, 2], None, true).unwrap();
        let c = compactness(&d, PathOptions::default()).unwrap();
        assert!((c - 2.5 / 6.0).abs() < 1e-15);
        let all = PathOptions {
            mode: DegreeMode::All,
            ..Default::default()
        };
        assert!((compactness(&d, all).unwrap() - 5.0 / 6.0).abs() < 1e-15);
        // Weighted distances.
        let w = Graph::from_edges(3, &[0, 1], &[1, 2], Some(&[2.0, 2.0]), false).unwrap();
        let opts = PathOptions {
            weighted: true,
            ..Default::default()
        };
        assert!((compactness(&w, opts).unwrap() - 5.0 / 12.0).abs() < 1e-15);
    }
}
