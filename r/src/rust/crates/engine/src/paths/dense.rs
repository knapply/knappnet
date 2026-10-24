//! All-pairs shortest paths on a dense distance matrix: the Floyd-Warshall
//! algorithm (Floyd 1962, `floyd1962algorithm`; Warshall 1962,
//! `warshall1962theorem`).
//!
//! The matrix starts with the shortest arc between each ordered pair (0 on
//! the diagonal) and, for every intermediate node k in turn, replaces
//! d(u, v) by d(u, k) + d(k, v) when that is strictly smaller, recording
//! the predecessor of v on the path from u. O(n^3) time and O(n^2) memory,
//! independent of the number of arcs, so it suits small dense networks and
//! any sign of length; [`crate::paths::search::all_pairs`] and
//! [`crate::paths::negative::johnson`] are faster on sparse ones.
//!
//! Conventions:
//! - Lengths of any sign. A negative cycle makes some diagonal entry
//!   negative; the result then flags it and its distances are not shortest
//!   path lengths (NetworkX raises `NetworkXUnbounded`). A negative
//!   self-loop is reported the same way.
//! - Parallel arcs count through the shortest; other self-loops are
//!   ignored (the diagonal is 0), as in NetworkX's
//!   `floyd_warshall_predecessor_and_distance`.
//! - `pred[u][v]` is the node before v on the path from u ([`NO_NODE`] when
//!   v is u or unreachable). Ties keep the first path found, as NetworkX
//!   does (an update needs a strictly smaller length), and since row k does
//!   not change while k is the intermediate node, the rows are updated in
//!   parallel with results identical to the sequential order for any
//!   thread count.

use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use crate::parallel;
use crate::paths::{PathGraph, NO_NODE};

/// The result of [`floyd_warshall`].
#[derive(Debug, Clone, PartialEq)]
pub struct DenseResult {
    /// Number of nodes.
    pub n: usize,
    /// Distances, row-major (row = source).
    pub dist: Vec<f64>,
    /// Predecessors, row-major.
    pub pred: Vec<u32>,
    /// Whether a negative cycle was found (distances are then meaningless).
    pub negative_cycle: bool,
}

/// All-pairs distances and predecessors by the Floyd-Warshall algorithm.
/// See the module documentation.
pub fn floyd_warshall(pg: &PathGraph, threads: usize) -> DenseResult {
    let n = pg.node_count();
    let mut dist = vec![f64::INFINITY; n * n];
    let mut pred = vec![NO_NODE; n * n];
    let mut negative_cycle = false;
    for u in 0..n {
        for a in pg.out_arcs(u) {
            let v = pg.head(a) as usize;
            let len = pg.len(a);
            if v == u {
                if len < 0.0 {
                    negative_cycle = true;
                }
                continue;
            }
            if len < dist[u * n + v] {
                dist[u * n + v] = len;
                pred[u * n + v] = u as u32;
            }
        }
        dist[u * n + u] = 0.0;
    }
    if n == 0 {
        return DenseResult {
            n,
            dist,
            pred,
            negative_cycle,
        };
    }
    let threads = parallel::threads_for(threads, n, (n as u64).pow(3));
    if threads <= 1 {
        for k in 0..n {
            let (before, rest) = dist.split_at_mut(k * n);
            let (row_k, after) = rest.split_at_mut(n);
            let (pbefore, prest) = pred.split_at_mut(k * n);
            let (prow_k, pafter) = prest.split_at_mut(n);
            let rows = before.chunks_mut(n).zip(pbefore.chunks_mut(n));
            let rows = rows.chain(after.chunks_mut(n).zip(pafter.chunks_mut(n)));
            for (drow, prow) in rows {
                relax_row(drow, prow, row_k, prow_k, k);
            }
        }
    } else {
        let d: Vec<AtomicU64> = dist.iter().map(|x| AtomicU64::new(x.to_bits())).collect();
        let p: Vec<AtomicU32> = pred.iter().map(|&x| AtomicU32::new(x)).collect();
        parallel::crew(threads, |crew| {
            let mut row_k = vec![0.0f64; n];
            let mut prow_k = vec![0u32; n];
            let mine = crew.share(n);
            for k in 0..n {
                // Row k does not change while k is the intermediate node.
                for v in 0..n {
                    row_k[v] = f64::from_bits(d[k * n + v].load(Ordering::Relaxed));
                    prow_k[v] = p[k * n + v].load(Ordering::Relaxed);
                }
                for u in mine.clone() {
                    if u == k {
                        continue;
                    }
                    let duk = f64::from_bits(d[u * n + k].load(Ordering::Relaxed));
                    if duk == f64::INFINITY {
                        continue;
                    }
                    let (drow, prow) = (&d[u * n..(u + 1) * n], &p[u * n..(u + 1) * n]);
                    for v in 0..n {
                        let alt = duk + row_k[v];
                        if v != k && alt < f64::from_bits(drow[v].load(Ordering::Relaxed)) {
                            drow[v].store(alt.to_bits(), Ordering::Relaxed);
                            prow[v].store(prow_k[v], Ordering::Relaxed);
                        }
                    }
                }
                crew.wait();
            }
        });
        dist = d
            .into_iter()
            .map(|x| f64::from_bits(x.into_inner()))
            .collect();
        pred = p.into_iter().map(|x| x.into_inner()).collect();
    }
    if (0..n).any(|u| dist[u * n + u] < 0.0) {
        negative_cycle = true;
    }
    DenseResult {
        n,
        dist,
        pred,
        negative_cycle,
    }
}

/// Relaxes one row through intermediate node k; returns whether it changed.
#[inline]
fn relax_row(drow: &mut [f64], prow: &mut [u32], row_k: &[f64], prow_k: &[u32], k: usize) -> bool {
    let duk = drow[k];
    if duk == f64::INFINITY {
        return false;
    }
    let mut changed = false;
    for v in 0..drow.len() {
        if v == k {
            continue;
        }
        let alt = duk + row_k[v];
        if alt < drow[v] {
            drow[v] = alt;
            prow[v] = prow_k[v];
            changed = true;
        }
    }
    changed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::paths::negative::johnson;
    use crate::paths::{ArcLengths, PathGraphOptions};
    use crate::rng::Stream;

    #[test]
    fn matches_johnson_for_any_thread_count() {
        for seed in 0..60u64 {
            let mut rng = Stream::new(seed, 9);
            let n = 1 + rng.below(12) as usize;
            let m = rng.below(4 * n as u64) as usize;
            let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let len: Vec<f64> = (0..m).map(|_| rng.below(10) as f64 - 1.0).collect();
            let g = Graph::from_edges(n, &s, &t, None, true).unwrap();
            let pg = PathGraph::new(
                &g,
                PathGraphOptions {
                    lengths: ArcLengths::PerEdge(&len),
                    ..Default::default()
                },
            )
            .unwrap();
            let one = floyd_warshall(&pg, 1);
            let four = floyd_warshall(&pg, 4);
            assert_eq!(one, four);
            match johnson(&pg, None, 1).unwrap() {
                Ok((d, _)) => {
                    assert!(!one.negative_cycle);
                    assert_eq!(one.dist, d, "seed {seed}");
                    // Predecessors trace paths of the right length.
                    for u in 0..n {
                        for v in 0..n {
                            if u == v || one.dist[u * n + v].is_infinite() {
                                continue;
                            }
                            let mut path = vec![v as u32];
                            let mut cur = v;
                            while cur != u {
                                cur = one.pred[u * n + cur] as usize;
                                path.push(cur as u32);
                            }
                            path.reverse();
                            assert_eq!(pg.path_length(&path), Some(one.dist[u * n + v]));
                        }
                    }
                }
                Err(_) => assert!(one.negative_cycle, "seed {seed}"),
            }
        }
    }
}
