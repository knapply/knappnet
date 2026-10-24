//! A* search (Hart, Nilsson, and Raphael 1968, `hart1968formal`): Dijkstra's
//! algorithm ordered by g(v) + h(v), where g is the length of the best path
//! found so far and h(v) estimates the remaining length from v to the
//! target.
//!
//! The heuristic is given as one value per node (the estimate for the
//! target of the query), so any language binding can compute it however it
//! likes. A node is expanded again whenever a shorter path to it is found,
//! so the result is a shortest path whenever the heuristic never
//! overestimates (is admissible), even if it is not consistent; with a
//! consistent heuristic every node is expanded at most once. With an
//! inadmissible heuristic the path may be longer than the shortest, as
//! NetworkX documents for `astar_path`.
//!
//! `cutoff`: a node whose g + h exceeds the cutoff is not queued, as
//! NetworkX's `cutoff` for `astar_path` (where a cutoff of 0 means no
//! cutoff; here an infinite cutoff means none). Lengths must be
//! non-negative. Ties in g + h are settled by node index.

use crate::error::{Error, Result};
use crate::paths::{PathGraph, NO_NODE};
use crate::weights::IndexedHeap;

/// Length and nodes of a path from `s` to `t` found by A* with heuristic
/// values `h` (one per node), or `None` if the target was not reached.
pub fn astar(
    pg: &PathGraph,
    s: u32,
    t: u32,
    h: &[f64],
    cutoff: f64,
) -> Result<Option<(f64, Vec<u32>)>> {
    let n = pg.node_count();
    pg.check_node(s as usize, "source")?;
    pg.check_node(t as usize, "target")?;
    if h.len() != n {
        return Err(Error::LengthMismatch {
            what: "heuristic",
            expected: n,
            found: h.len(),
        });
    }
    if let Some(v) = h.iter().position(|x| x.is_nan()) {
        return Err(Error::InvalidArgument(format!(
            "the heuristic of node {v} is NaN"
        )));
    }
    if cutoff.is_nan() {
        return Err(Error::InvalidArgument("the cutoff must not be NaN".into()));
    }
    pg.require_nonnegative("A* search")?;
    let mut g = vec![f64::INFINITY; n];
    let mut parent = vec![NO_NODE; n];
    let mut heap = IndexedHeap::new(n);
    g[s as usize] = 0.0;
    heap.push_or_decrease(s, h[s as usize]);
    while let Some((_, v)) = heap.pop() {
        let v = v as usize;
        if v == t as usize {
            let mut path = vec![t];
            let mut cur = v;
            while parent[cur] != NO_NODE {
                cur = parent[cur] as usize;
                path.push(cur as u32);
            }
            path.reverse();
            return Ok(Some((g[v], path)));
        }
        let gv = g[v];
        for a in pg.out_arcs(v) {
            let w = pg.head(a) as usize;
            let alt = gv + pg.len(a);
            if alt >= g[w] || alt + h[w] > cutoff {
                continue;
            }
            g[w] = alt;
            parent[w] = v as u32;
            heap.push_or_decrease(w as u32, alt + h[w]);
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::paths::search::{single_source, SearchOptions};
    use crate::paths::{ArcLengths, PathGraphOptions};
    use crate::rng::Stream;
    use crate::weights::WeightKind;

    #[test]
    fn admissible_heuristics_give_shortest_paths() {
        for seed in 0..30u64 {
            let mut rng = Stream::new(seed, 1);
            let (n, m) = (15usize, 35usize);
            let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let w: Vec<f64> = (0..m).map(|_| (1 + rng.below(20)) as f64 / 4.0).collect();
            let g = Graph::from_edges(n, &s, &t, Some(&w), true).unwrap();
            let pg = PathGraph::new(
                &g,
                PathGraphOptions {
                    lengths: ArcLengths::Weights(WeightKind::Distance),
                    ..Default::default()
                },
            )
            .unwrap();
            for target in 0..n as u32 {
                let back = single_source(
                    &pg,
                    target,
                    SearchOptions {
                        reverse: true,
                        ..Default::default()
                    },
                )
                .unwrap();
                // Exact distances (consistent), and an admissible but
                // inconsistent mixture.
                let exact: Vec<f64> = back
                    .dist
                    .iter()
                    .map(|&d| if d.is_finite() { d } else { 0.0 })
                    .collect();
                let mixed: Vec<f64> = exact
                    .iter()
                    .enumerate()
                    .map(|(v, &d)| if v % 3 == 0 { d } else { 0.1 * d })
                    .collect();
                for source in 0..n as u32 {
                    let want = back.dist[source as usize];
                    for h in [&exact, &mixed, &vec![0.0; n]] {
                        match astar(&pg, source, target, h, f64::INFINITY).unwrap() {
                            None => assert!(want.is_infinite()),
                            Some((len, path)) => {
                                assert_eq!(len, want);
                                assert_eq!(pg.path_length(&path), Some(len));
                            }
                        }
                    }
                }
            }
        }
    }
}
