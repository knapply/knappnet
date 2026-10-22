//! Degree and strength.
//!
//! Conventions (matching igraph and NetworkX defaults):
//! - Parallel edges each count.
//! - In an undirected graph a self-loop adds 2 to its node's degree.
//! - In a directed graph a self-loop adds 1 to out-degree and 1 to in-degree.
//! - `loops = false` ignores self-loops entirely, which reproduces the
//!   `sna` default of ignoring the diagonal (`diag = FALSE`).

use crate::error::{Error, Result};
use crate::graph::Graph;

/// Which edge endpoints to count in a directed graph. Undirected graphs
/// ignore the mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeMode {
    /// Edges leaving the node.
    Out,
    /// Edges entering the node.
    In,
    /// Both.
    All,
}

impl DegreeMode {
    /// Decodes the integer codes used by the C interface: 1 = out, 2 = in,
    /// 3 = all.
    pub fn from_code(code: i32) -> Result<DegreeMode> {
        match code {
            1 => Ok(DegreeMode::Out),
            2 => Ok(DegreeMode::In),
            3 => Ok(DegreeMode::All),
            _ => Err(Error::InvalidArgument(format!(
                "unknown degree mode code {code}"
            ))),
        }
    }
}

/// Number of edge endpoints at each node.
pub fn degree(g: &Graph, mode: DegreeMode, loops: bool) -> Vec<f64> {
    accumulate(g, mode, loops, |_| 1.0)
}

/// Sum of incident edge weights at each node.
pub fn strength(g: &Graph, mode: DegreeMode, loops: bool) -> Result<Vec<f64>> {
    let w = g.edge_weights().ok_or(Error::MissingWeights)?;
    Ok(accumulate(g, mode, loops, |e| w[e]))
}

fn accumulate(g: &Graph, mode: DegreeMode, loops: bool, value: impl Fn(usize) -> f64) -> Vec<f64> {
    let mut out = vec![0.0; g.node_count()];
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let directed = g.is_directed();
    for e in 0..g.edge_count() {
        let (s, t) = (src[e] as usize, tgt[e] as usize);
        if s == t && !loops {
            continue;
        }
        let x = value(e);
        if !directed {
            out[s] += x;
            out[t] += x;
            continue;
        }
        if mode != DegreeMode::In {
            out[s] += x;
        }
        if mode != DegreeMode::Out {
            out[t] += x;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loops_and_multi_edges() {
        // 0->1 twice, 1->1 loop, 2->0
        let g = Graph::from_edges(3, &[0, 0, 1, 2], &[1, 1, 1, 0], None, true).unwrap();
        assert_eq!(degree(&g, DegreeMode::Out, true), vec![2.0, 1.0, 1.0]);
        assert_eq!(degree(&g, DegreeMode::In, true), vec![1.0, 3.0, 0.0]);
        assert_eq!(degree(&g, DegreeMode::All, true), vec![3.0, 4.0, 1.0]);
        assert_eq!(degree(&g, DegreeMode::All, false), vec![3.0, 2.0, 1.0]);

        let u = Graph::from_edges(3, &[0, 0, 1, 2], &[1, 1, 1, 0], None, false).unwrap();
        assert_eq!(degree(&u, DegreeMode::All, true), vec![3.0, 4.0, 1.0]);
    }

    #[test]
    fn handshake_lemma() {
        let g = Graph::from_edges(4, &[0, 1, 2, 3, 3], &[1, 2, 3, 0, 3], None, false).unwrap();
        let total: f64 = degree(&g, DegreeMode::All, true).iter().sum();
        assert_eq!(total, 2.0 * g.edge_count() as f64);
    }
}
