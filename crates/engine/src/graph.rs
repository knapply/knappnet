//! Compressed sparse row (CSR) graph storage.
//!
//! Nodes are the integers `0..n`. The R and Python layers map user-facing
//! node identifiers to these indices, so the engine never sees strings.
//!
//! Outgoing arcs are stored in CSR form: the arcs leaving node `v` are
//! `targets[offsets[v]..offsets[v + 1]]`, sorted by target, with weights (if
//! any) in a parallel array. Directed graphs additionally store incoming arcs
//! in the same layout (equivalently, the CSC form of the adjacency matrix).
//! For undirected graphs every edge is stored once in each endpoint's list,
//! except self-loops, which appear once.
//!
//! Parallel edges are preserved. Algorithms decide how to treat them and
//! document that choice.

use crate::error::{Error, Result};

/// One direction of adjacency in CSR layout.
#[derive(Debug, Clone)]
pub struct Csr {
    offsets: Vec<usize>,
    targets: Vec<u32>,
    weights: Option<Vec<f64>>,
}

impl Csr {
    /// Neighbors of `v`, sorted ascending (parallel edges repeat).
    #[inline]
    pub fn neighbors(&self, v: usize) -> &[u32] {
        &self.targets[self.offsets[v]..self.offsets[v + 1]]
    }

    /// Weights aligned with `neighbors(v)`, if the graph is weighted.
    #[inline]
    pub fn weights(&self, v: usize) -> Option<&[f64]> {
        self.weights
            .as_ref()
            .map(|w| &w[self.offsets[v]..self.offsets[v + 1]])
    }

    #[inline]
    pub fn degree(&self, v: usize) -> usize {
        self.offsets[v + 1] - self.offsets[v]
    }

    /// Builds a CSR from arc endpoints with a counting sort, so each row is
    /// sorted by target in O(n + m) time.
    fn build(n: usize, from: &[u32], to: &[u32], weights: Option<&[f64]>) -> Csr {
        let m = from.len();

        // Pass 1: order arcs by target (stable counting sort).
        let mut by_target = vec![0usize; n + 1];
        for &t in to {
            by_target[t as usize + 1] += 1;
        }
        for i in 0..n {
            by_target[i + 1] += by_target[i];
        }
        let mut order = vec![0usize; m];
        let mut cursor = by_target.clone();
        for (arc, &t) in to.iter().enumerate() {
            order[cursor[t as usize]] = arc;
            cursor[t as usize] += 1;
        }

        // Pass 2: stable counting sort by source. Rows end up sorted by target.
        let mut offsets = vec![0usize; n + 1];
        for &f in from {
            offsets[f as usize + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let mut cursor = offsets.clone();
        let mut targets = vec![0u32; m];
        let mut out_weights = weights.map(|_| vec![0f64; m]);
        for &arc in &order {
            let f = from[arc] as usize;
            let slot = cursor[f];
            targets[slot] = to[arc];
            if let (Some(dst), Some(src)) = (out_weights.as_mut(), weights) {
                dst[slot] = src[arc];
            }
            cursor[f] += 1;
        }

        Csr {
            offsets,
            targets,
            weights: out_weights,
        }
    }
}

/// An immutable graph optimized for analysis.
#[derive(Debug, Clone)]
pub struct Graph {
    n: usize,
    directed: bool,
    sources: Vec<u32>,
    targets: Vec<u32>,
    weights: Option<Vec<f64>>,
    out: Csr,
    /// Incoming arcs. `None` for undirected graphs, where `out` already holds
    /// both directions.
    inc: Option<Csr>,
}

impl Graph {
    /// Builds a graph from an edge list of 0-based node indices.
    ///
    /// `weights`, if given, must have one finite value per edge. Whether
    /// weights may be zero or negative is checked by each algorithm, since
    /// requirements differ (for example, strength allows any finite value,
    /// while shortest paths require positive weights).
    pub fn from_edges(
        n: usize,
        sources: &[i64],
        targets: &[i64],
        weights: Option<&[f64]>,
        directed: bool,
    ) -> Result<Graph> {
        if n > u32::MAX as usize {
            return Err(Error::TooLarge {
                what: "node count",
                value: n,
            });
        }
        if targets.len() != sources.len() {
            return Err(Error::LengthMismatch {
                what: "targets",
                expected: sources.len(),
                found: targets.len(),
            });
        }
        if let Some(w) = weights {
            if w.len() != sources.len() {
                return Err(Error::LengthMismatch {
                    what: "weights",
                    expected: sources.len(),
                    found: w.len(),
                });
            }
            if let Some(edge) = w.iter().position(|x| !x.is_finite()) {
                return Err(Error::NonFiniteWeight { edge });
            }
        }
        let m = sources.len();
        if m > (u32::MAX as usize) / 2 {
            return Err(Error::TooLarge {
                what: "edge count",
                value: m,
            });
        }

        let check = |edge: usize, node: i64| -> Result<u32> {
            if node < 0 || node as u64 >= n as u64 {
                Err(Error::NodeOutOfRange { edge, node, n })
            } else {
                Ok(node as u32)
            }
        };
        let mut src = Vec::with_capacity(m);
        let mut tgt = Vec::with_capacity(m);
        for e in 0..m {
            src.push(check(e, sources[e])?);
            tgt.push(check(e, targets[e])?);
        }
        let weights = weights.map(|w| w.to_vec());

        let (out, inc) = if directed {
            let out = Csr::build(n, &src, &tgt, weights.as_deref());
            let inc = Csr::build(n, &tgt, &src, weights.as_deref());
            (out, Some(inc))
        } else {
            // Each non-loop edge becomes two arcs; a loop stays one arc.
            let loops = src.iter().zip(&tgt).filter(|(a, b)| a == b).count();
            let arcs = 2 * m - loops;
            let mut from = Vec::with_capacity(arcs);
            let mut to = Vec::with_capacity(arcs);
            let mut w = weights.as_ref().map(|_| Vec::with_capacity(arcs));
            for e in 0..m {
                from.push(src[e]);
                to.push(tgt[e]);
                if let (Some(dst), Some(srcw)) = (w.as_mut(), weights.as_ref()) {
                    dst.push(srcw[e]);
                }
                if src[e] != tgt[e] {
                    from.push(tgt[e]);
                    to.push(src[e]);
                    if let (Some(dst), Some(srcw)) = (w.as_mut(), weights.as_ref()) {
                        dst.push(srcw[e]);
                    }
                }
            }
            (Csr::build(n, &from, &to, w.as_deref()), None)
        };

        Ok(Graph {
            n,
            directed,
            sources: src,
            targets: tgt,
            weights,
            out,
            inc,
        })
    }

    #[inline]
    pub fn node_count(&self) -> usize {
        self.n
    }

    #[inline]
    pub fn edge_count(&self) -> usize {
        self.sources.len()
    }

    #[inline]
    pub fn is_directed(&self) -> bool {
        self.directed
    }

    #[inline]
    pub fn is_weighted(&self) -> bool {
        self.weights.is_some()
    }

    /// Edge sources in input order.
    pub fn edge_sources(&self) -> &[u32] {
        &self.sources
    }

    /// Edge targets in input order.
    pub fn edge_targets(&self) -> &[u32] {
        &self.targets
    }

    /// Edge weights in input order.
    pub fn edge_weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// Outgoing adjacency. For undirected graphs, all incident edges.
    #[inline]
    pub fn out(&self) -> &Csr {
        &self.out
    }

    /// Incoming adjacency. For undirected graphs this is the same as `out`.
    #[inline]
    pub fn inc(&self) -> &Csr {
        self.inc.as_ref().unwrap_or(&self.out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directed_rows_are_sorted_and_complete() {
        let g = Graph::from_edges(4, &[2, 0, 0, 1, 0], &[1, 3, 1, 2, 1], None, true).unwrap();
        assert_eq!(g.out().neighbors(0), &[1, 1, 3]);
        assert_eq!(g.out().neighbors(1), &[2]);
        assert_eq!(g.out().neighbors(2), &[1]);
        assert_eq!(g.inc().neighbors(1), &[0, 0, 2]);
        assert_eq!(g.inc().neighbors(3), &[0]);
    }

    #[test]
    fn undirected_stores_both_directions_and_loops_once() {
        let g =
            Graph::from_edges(3, &[0, 1, 2], &[1, 2, 2], Some(&[1.0, 2.0, 3.0]), false).unwrap();
        assert_eq!(g.out().neighbors(1), &[0, 2]);
        assert_eq!(g.out().weights(1).unwrap(), &[1.0, 2.0]);
        assert_eq!(g.out().neighbors(2), &[1, 2]);
        assert_eq!(g.out().weights(2).unwrap(), &[2.0, 3.0]);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(matches!(
            Graph::from_edges(2, &[0], &[2], None, true),
            Err(Error::NodeOutOfRange { .. })
        ));
        assert!(matches!(
            Graph::from_edges(2, &[0], &[1, 1], None, true),
            Err(Error::LengthMismatch { .. })
        ));
        assert!(matches!(
            Graph::from_edges(2, &[0], &[1], Some(&[f64::NAN]), true),
            Err(Error::NonFiniteWeight { .. })
        ));
    }
}
