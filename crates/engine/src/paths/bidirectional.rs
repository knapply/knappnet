//! Shortest path between one pair of nodes by bidirectional search
//! (Nicholson 1966, `nicholson1966finding`): one search grows from the
//! source along arcs and one from the target against them, until they
//! meet.
//!
//! - Unit lengths: breadth-first search, expanding the smaller frontier
//!   one whole level at a time; when the frontiers touch, the shortest of
//!   the connections found in that level is returned (NetworkX's
//!   `bidirectional_shortest_path`).
//! - Other non-negative lengths: Dijkstra's algorithm in both directions,
//!   alternating one settled node at a time, keeping the shortest
//!   connection seen, and stopping when some node is settled by both
//!   searches (Nicholson's stopping rule, as in NetworkX's
//!   `bidirectional_dijkstra`).
//!
//! The result is the length and one shortest path, or `None` when the
//! target cannot be reached. Ties are settled by node index, so the path is
//! deterministic; NetworkX may return another path of the same length.

use crate::error::Result;
use crate::paths::{PathGraph, NO_NODE};
use crate::weights::IndexedHeap;

/// Length and nodes of a shortest path from `s` to `t`, or `None` if there
/// is none. Lengths must be non-negative.
pub fn bidirectional(pg: &PathGraph, s: u32, t: u32) -> Result<Option<(f64, Vec<u32>)>> {
    pg.check_node(s as usize, "source")?;
    pg.check_node(t as usize, "target")?;
    if s == t {
        return Ok(Some((0.0, vec![s])));
    }
    if pg.is_unit() {
        return Ok(bfs(pg, s, t));
    }
    pg.require_nonnegative("bidirectional Dijkstra")?;
    Ok(dijkstra(pg, s, t))
}

/// Neighbours of `v` in direction `dir` (0 forward, 1 backward), as
/// (arc, node) pairs.
fn arcs(pg: &PathGraph, v: usize, dir: usize) -> impl Iterator<Item = (usize, usize)> + '_ {
    let (out, inn) = if dir == 0 {
        (pg.out_arcs(v), &[][..])
    } else {
        (0..0, pg.in_arcs(v))
    };
    out.chain(inn.iter().map(|&a| a as usize)).map(move |a| {
        let w = if dir == 0 { pg.head(a) } else { pg.tail(a) };
        (a, w as usize)
    })
}

fn join(parent: &[Vec<u32>; 2], meet: usize) -> Vec<u32> {
    let mut path = vec![meet as u32];
    let mut cur = meet;
    while parent[0][cur] != NO_NODE {
        cur = parent[0][cur] as usize;
        path.push(cur as u32);
    }
    path.reverse();
    let mut cur = meet;
    while parent[1][cur] != NO_NODE {
        cur = parent[1][cur] as usize;
        path.push(cur as u32);
    }
    path
}

fn bfs(pg: &PathGraph, s: u32, t: u32) -> Option<(f64, Vec<u32>)> {
    let n = pg.node_count();
    let mut depth = [vec![u32::MAX; n], vec![u32::MAX; n]];
    let mut parent = [vec![NO_NODE; n], vec![NO_NODE; n]];
    let mut frontier = [vec![s], vec![t]];
    depth[0][s as usize] = 0;
    depth[1][t as usize] = 0;
    while !frontier[0].is_empty() && !frontier[1].is_empty() {
        let dir = if frontier[0].len() <= frontier[1].len() {
            0
        } else {
            1
        };
        let other = 1 - dir;
        let mut next = Vec::new();
        // The best connection found in this level: (total hops, meeting
        // node, the node it is reached from in direction `dir`).
        let mut best: Option<(u32, usize, usize)> = None;
        for &v in &frontier[dir] {
            let v = v as usize;
            let dv = depth[dir][v];
            for (_, w) in arcs(pg, v, dir) {
                if depth[other][w] != u32::MAX {
                    let total = dv + 1 + depth[other][w];
                    if best.map_or(true, |b| total < b.0) {
                        best = Some((total, w, v));
                    }
                }
                if depth[dir][w] == u32::MAX {
                    depth[dir][w] = dv + 1;
                    parent[dir][w] = v as u32;
                    next.push(w as u32);
                }
            }
        }
        if let Some((total, w, v)) = best {
            // The connection is the arc between v (on side `dir`) and w
            // (already reached from the other side). Make the meeting node
            // the one on the forward side and route through the arc.
            let mut parent = parent;
            let meet = if dir == 0 {
                parent[1][v] = w as u32;
                v
            } else {
                parent[1][w] = v as u32;
                w
            };
            return Some((total as f64, join(&parent, meet)));
        }
        frontier[dir] = next;
    }
    None
}

fn dijkstra(pg: &PathGraph, s: u32, t: u32) -> Option<(f64, Vec<u32>)> {
    let n = pg.node_count();
    let mut dist = [vec![f64::INFINITY; n], vec![f64::INFINITY; n]];
    let mut done = [vec![false; n], vec![false; n]];
    let mut parent = [vec![NO_NODE; n], vec![NO_NODE; n]];
    let mut heap = [IndexedHeap::new(n), IndexedHeap::new(n)];
    dist[0][s as usize] = 0.0;
    dist[1][t as usize] = 0.0;
    heap[0].push_or_decrease(s, 0.0);
    heap[1].push_or_decrease(t, 0.0);
    let mut best: Option<(f64, usize)> = None;
    let mut dir = 1;
    loop {
        dir = 1 - dir;
        let Some((d, v)) = heap[dir].pop() else {
            break;
        };
        let v = v as usize;
        done[dir][v] = true;
        if done[1 - dir][v] {
            break;
        }
        for (a, w) in arcs(pg, v, dir) {
            if done[dir][w] {
                continue;
            }
            let alt = d + pg.len(a);
            if alt < dist[dir][w] {
                dist[dir][w] = alt;
                parent[dir][w] = v as u32;
                heap[dir].push_or_decrease(w as u32, alt);
            }
            let other = dist[1 - dir][w];
            if other.is_finite() {
                let total = dist[dir][w] + other;
                if best.map_or(true, |b| total < b.0) {
                    best = Some((total, w));
                }
            }
        }
        if heap[dir].entries().is_empty() {
            break;
        }
    }
    best.map(|(len, meet)| (len, join(&parent, meet)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::paths::search::{single_source, SearchOptions};
    use crate::paths::{ArcLengths, PathGraphOptions};
    use crate::rng::Stream;
    use crate::weights::WeightKind;

    fn random(n: usize, m: usize, seed: u64, weighted: bool, directed: bool) -> PathGraph {
        let mut rng = Stream::new(seed, 0);
        let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
        let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
        let w: Vec<f64> = (0..m).map(|_| (1 + rng.below(9)) as f64 / 2.0).collect();
        let g = Graph::from_edges(n, &s, &t, Some(&w), directed).unwrap();
        let lengths = if weighted {
            ArcLengths::Weights(WeightKind::Distance)
        } else {
            ArcLengths::Unit
        };
        PathGraph::new(
            &g,
            PathGraphOptions {
                lengths,
                ..Default::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn agrees_with_one_directional_search() {
        for seed in 0..40 {
            for &(weighted, directed) in
                &[(false, false), (false, true), (true, false), (true, true)]
            {
                let pg = random(12, 20, seed, weighted, directed);
                for s in 0..12u32 {
                    let tree = single_source(&pg, s, SearchOptions::default()).unwrap();
                    for t in 0..12u32 {
                        let got = bidirectional(&pg, s, t).unwrap();
                        let want = tree.dist[t as usize];
                        match got {
                            None => assert!(want.is_infinite(), "{seed} {s} {t}"),
                            Some((len, path)) => {
                                assert_eq!(len, want, "{seed} {s} {t}");
                                assert_eq!(path[0], s);
                                assert_eq!(*path.last().unwrap(), t);
                                assert_eq!(pg.path_length(&path), Some(len), "{path:?}");
                            }
                        }
                    }
                }
            }
        }
    }
}
