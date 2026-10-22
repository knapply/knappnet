//! Shortest-path betweenness centrality.
//!
//! Freeman's betweenness (Freeman 1977, `freeman1977set`) computed with
//! Brandes' algorithm (Brandes 2001, `brandes2001faster`): one single-source
//! shortest-path search per node, followed by dependency accumulation in
//! order of non-increasing distance. Total work is O(nm) for unweighted and
//! O(nm + n^2 log n) for weighted graphs.
//!
//! Accumulation walks each node's outgoing arcs and tests whether the arc
//! lies on a shortest path, instead of storing predecessor lists. This is the
//! successor formulation discussed in Brandes (2008, `brandes2008variants`)
//! and saves O(m) memory per worker thread.
//!
//! Conventions:
//! - Endpoints are excluded.
//! - Parallel edges count as distinct shortest paths, as in igraph.
//! - Self-loops never lie on a shortest path and are ignored.
//! - In undirected graphs each unordered pair is counted once (the raw sum is
//!   halved), matching igraph and NetworkX.
//! - Weights are distances and must be strictly positive. Two path lengths
//!   are considered equal when they agree to a relative tolerance of 1e-10,
//!   the same policy igraph uses, so floating-point rounding does not drop
//!   tied shortest paths.
//! - Normalization divides by the number of ordered (directed) or unordered
//!   (undirected) pairs of nodes other than the focal node:
//!   (n - 1)(n - 2) or (n - 1)(n - 2) / 2.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;

const EPSILON: f64 = 1e-10;

/// Options for [`betweenness`].
#[derive(Debug, Clone, Copy, Default)]
pub struct BetweennessOptions {
    /// Use edge weights as distances (they must be strictly positive).
    pub weighted: bool,
    /// Divide by the number of pairs not involving the node.
    pub normalized: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Shortest-path betweenness of every node, in node order.
///
/// See the module documentation for conventions.
pub fn betweenness(g: &Graph, opts: BetweennessOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    if opts.weighted {
        let w = g.edge_weights().ok_or(Error::MissingWeights)?;
        if let Some(edge) = w.iter().position(|&x| x <= 0.0) {
            return Err(Error::NonPositiveWeight {
                edge,
                weight: w[edge],
            });
        }
    }

    let partials = parallel::for_each_item(
        n,
        opts.threads,
        || Workspace::new(n),
        |ws, s| {
            if opts.weighted {
                ws.single_source_weighted(g, s);
            } else {
                ws.single_source_unweighted(g, s);
            }
        },
    );

    let mut bc = vec![0.0; n];
    for ws in partials {
        for (b, x) in bc.iter_mut().zip(ws.bc) {
            *b += x;
        }
    }

    let mut scale = if g.is_directed() { 1.0 } else { 0.5 };
    if opts.normalized && n > 2 {
        let pairs = ((n - 1) * (n - 2)) as f64;
        scale *= if g.is_directed() {
            1.0 / pairs
        } else {
            2.0 / pairs
        };
    }
    if scale != 1.0 {
        for b in &mut bc {
            *b *= scale;
        }
    }
    Ok(bc)
}

#[inline]
fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPSILON * a.abs().max(b.abs())
}

/// Per-thread scratch space, reused across sources.
struct Workspace {
    bc: Vec<f64>,
    sigma: Vec<f64>,
    delta: Vec<f64>,
    hops: Vec<u32>,
    dist: Vec<f64>,
    /// Position at which each node was settled; u32::MAX if not yet settled.
    rank: Vec<u32>,
    /// Nodes in the order they were settled (non-decreasing distance).
    order: Vec<u32>,
    heap: BinaryHeap<HeapItem>,
}

impl Workspace {
    fn new(n: usize) -> Workspace {
        Workspace {
            bc: vec![0.0; n],
            sigma: vec![0.0; n],
            delta: vec![0.0; n],
            hops: vec![u32::MAX; n],
            dist: vec![f64::INFINITY; n],
            rank: vec![u32::MAX; n],
            order: Vec::with_capacity(n),
            heap: BinaryHeap::new(),
        }
    }

    fn single_source_unweighted(&mut self, g: &Graph, s: usize) {
        let adj = g.out();
        self.order.clear();
        self.hops[s] = 0;
        self.sigma[s] = 1.0;
        self.order.push(s as u32);

        // Breadth-first search. `order` doubles as the queue.
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let next = self.hops[v] + 1;
            for &w in adj.neighbors(v) {
                let w = w as usize;
                if self.hops[w] == u32::MAX {
                    self.hops[w] = next;
                    self.order.push(w as u32);
                }
                if self.hops[w] == next {
                    self.sigma[w] += self.sigma[v];
                }
            }
        }

        // Accumulate dependencies over successors, farthest nodes first.
        for i in (0..self.order.len()).rev() {
            let v = self.order[i] as usize;
            let next = self.hops[v] + 1;
            let mut dv = 0.0;
            for &w in adj.neighbors(v) {
                let w = w as usize;
                if self.hops[w] == next {
                    dv += (1.0 + self.delta[w]) / self.sigma[w];
                }
            }
            self.delta[v] = self.sigma[v] * dv;
            if v != s {
                self.bc[v] += self.delta[v];
            }
        }

        for &v in &self.order {
            let v = v as usize;
            self.hops[v] = u32::MAX;
            self.sigma[v] = 0.0;
            self.delta[v] = 0.0;
        }
    }

    fn single_source_weighted(&mut self, g: &Graph, s: usize) {
        let adj = g.out();
        self.order.clear();
        self.heap.clear();
        self.dist[s] = 0.0;
        self.sigma[s] = 1.0;
        self.heap.push(HeapItem {
            dist: 0.0,
            node: s as u32,
        });

        // Dijkstra's algorithm with lazy deletion. Path counts propagate only
        // from settled nodes, so sigma[v] is final when v's arcs are relaxed.
        let mut touched: Vec<u32> = vec![s as u32];
        while let Some(HeapItem { dist: d, node }) = self.heap.pop() {
            let v = node as usize;
            if self.rank[v] != u32::MAX || d > self.dist[v] {
                continue;
            }
            self.rank[v] = self.order.len() as u32;
            self.order.push(v as u32);
            let weights = adj.weights(v).expect("weighted graph");
            for (&w, &len) in adj.neighbors(v).iter().zip(weights) {
                let w = w as usize;
                if self.rank[w] != u32::MAX {
                    continue;
                }
                let alt = d + len;
                let cur = self.dist[w];
                if cur.is_infinite() {
                    touched.push(w as u32);
                    self.dist[w] = alt;
                    self.sigma[w] = self.sigma[v];
                    self.heap.push(HeapItem {
                        dist: alt,
                        node: w as u32,
                    });
                } else if approx_eq(alt, cur) {
                    self.sigma[w] += self.sigma[v];
                } else if alt < cur {
                    self.dist[w] = alt;
                    self.sigma[w] = self.sigma[v];
                    self.heap.push(HeapItem {
                        dist: alt,
                        node: w as u32,
                    });
                }
            }
        }

        // Only nodes settled after v can be its successors on a shortest path.
        // Checking rank keeps accumulation consistent with how sigma was
        // counted, even when an arc is shorter than the comparison tolerance.
        for i in (0..self.order.len()).rev() {
            let v = self.order[i] as usize;
            let rank_v = self.rank[v];
            let dist_v = self.dist[v];
            let weights = adj.weights(v).expect("weighted graph");
            let mut dv = 0.0;
            for (&w, &len) in adj.neighbors(v).iter().zip(weights) {
                let w = w as usize;
                let dist_w = self.dist[w];
                if self.rank[w] > rank_v
                    && self.rank[w] != u32::MAX
                    && approx_eq(dist_v + len, dist_w)
                {
                    dv += (1.0 + self.delta[w]) / self.sigma[w];
                }
            }
            self.delta[v] = self.sigma[v] * dv;
            if v != s {
                self.bc[v] += self.delta[v];
            }
        }

        for &v in &touched {
            let v = v as usize;
            self.dist[v] = f64::INFINITY;
            self.sigma[v] = 0.0;
            self.delta[v] = 0.0;
            self.rank[v] = u32::MAX;
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct HeapItem {
    dist: f64,
    node: u32,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for HeapItem {}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapItem {
    // Reversed so that BinaryHeap (a max-heap) pops the smallest distance.
    // Distances are finite and non-negative, so total_cmp is a true order.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .dist
            .total_cmp(&self.dist)
            .then_with(|| other.node.cmp(&self.node))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bc(
        n: usize,
        s: &[i64],
        t: &[i64],
        w: Option<&[f64]>,
        directed: bool,
        opts: BetweennessOptions,
    ) -> Vec<f64> {
        let g = Graph::from_edges(n, s, t, w, directed).unwrap();
        betweenness(&g, opts).unwrap()
    }

    #[test]
    fn path_graph() {
        // 0 - 1 - 2 - 3: interior nodes lie on 2 pairs each.
        let r = bc(4, &[0, 1, 2], &[1, 2, 3], None, false, Default::default());
        assert_eq!(r, vec![0.0, 2.0, 2.0, 0.0]);
    }

    #[test]
    fn star_normalized_is_one() {
        let r = bc(
            5,
            &[0, 0, 0, 0],
            &[1, 2, 3, 4],
            None,
            false,
            BetweennessOptions {
                normalized: true,
                ..Default::default()
            },
        );
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!(r[1..].iter().all(|&x| x == 0.0));
    }

    #[test]
    fn directed_cycle() {
        // 0->1->2->0: each node is intermediate on exactly one ordered pair.
        let r = bc(3, &[0, 1, 2], &[1, 2, 0], None, true, Default::default());
        assert_eq!(r, vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn split_paths_share_credit() {
        // Diamond 0-1-3, 0-2-3: nodes 1 and 2 each get half of pair (0,3).
        let r = bc(
            4,
            &[0, 0, 1, 2],
            &[1, 2, 3, 3],
            None,
            false,
            Default::default(),
        );
        assert_eq!(r, vec![0.5, 0.5, 0.5, 0.5]);
    }

    #[test]
    fn weights_change_paths() {
        // Triangle where the direct edge 0-2 is longer than going through 1.
        let r = bc(
            3,
            &[0, 1, 0],
            &[1, 2, 2],
            Some(&[1.0, 1.0, 5.0]),
            false,
            BetweennessOptions {
                weighted: true,
                ..Default::default()
            },
        );
        assert_eq!(r, vec![0.0, 1.0, 0.0]);
    }

    #[test]
    fn tied_float_paths_are_detected() {
        // Both routes from 0 to 3 have length 0.3 but differ in the last bit.
        let r = bc(
            4,
            &[0, 1, 0, 2],
            &[1, 3, 2, 3],
            // Route via 1: 0.1 + 0.2 = 0.30000000000000004. Route via 2: 0.15 + 0.15 = 0.3.
            Some(&[0.1, 0.2, 0.15, 0.15]),
            true,
            BetweennessOptions {
                weighted: true,
                ..Default::default()
            },
        );
        assert!(
            (r[1] - 0.5).abs() < 1e-12 && (r[2] - 0.5).abs() < 1e-12,
            "{r:?}"
        );
    }

    #[test]
    fn rejects_non_positive_weights() {
        let g = Graph::from_edges(2, &[0], &[1], Some(&[0.0]), true).unwrap();
        let opts = BetweennessOptions {
            weighted: true,
            ..Default::default()
        };
        assert!(matches!(
            betweenness(&g, opts),
            Err(Error::NonPositiveWeight { .. })
        ));
    }
}
