//! Every shortest path between two nodes, and the number of shortest paths
//! from a source to every node.
//!
//! Both work on the predecessor lists of one search from the source
//! ([`crate::paths::search`]; Bellman-Ford when some length is negative):
//! a shortest path to t is t preceded by a shortest path to one of t's
//! predecessors. Paths are node sequences, so parallel edges do not
//! multiply them (NetworkX's `all_shortest_paths` agrees). With zero-length
//! arcs the predecessor lists can contain cycles; only simple paths are
//! listed, as NetworkX does, and counts follow only predecessors whose
//! distance became final earlier.
//!
//! The number of shortest paths can grow exponentially with the size of
//! the network, so [`all_shortest_paths`] takes a limit and fails when it
//! would be exceeded; [`count_shortest_paths`] gives the number first (as a
//! floating-point number, exact up to 2^53). Paths are listed in
//! lexicographic order of their node sequences; NetworkX lists them in an
//! order that depends on its adjacency dictionaries.

use crate::error::{Error, Result};
use crate::paths::negative::bellman_ford_predecessors;
use crate::paths::search::{SearchOptions, Searcher};
use crate::paths::{NodeLists, PathGraph, Tie};

/// Predecessor lists and distances from `s`, choosing Bellman-Ford when a
/// length is negative.
fn dag(pg: &PathGraph, s: u32, tie: Tie) -> Result<(NodeLists, Vec<f64>, Vec<u32>)> {
    let n = pg.node_count();
    if pg.min_len() < 0.0 {
        let (lists, dist) = match bellman_ford_predecessors(pg, &[s], tie)? {
            Ok(x) => x,
            Err(c) => {
                return Err(Error::InvalidArgument(format!(
                    "the network has a negative cycle through node {}, so shortest paths \
                     are unbounded",
                    c[0]
                )))
            }
        };
        let mut order: Vec<u32> = (0..n as u32)
            .filter(|&v| dist[v as usize].is_finite())
            .collect();
        order.sort_by(|&a, &b| {
            dist[a as usize]
                .total_cmp(&dist[b as usize])
                .then(a.cmp(&b))
        });
        let mut rank = vec![u32::MAX; n];
        for (i, &v) in order.iter().enumerate() {
            rank[v as usize] = i as u32;
        }
        return Ok((lists, dist, rank));
    }
    let mut ws = Searcher::new(n);
    ws.run(pg, &[s], SearchOptions::default())?;
    let lists = ws.predecessors(pg, tie, false);
    let dist = (0..n).map(|v| ws.dist(v)).collect();
    let rank = (0..n).map(|v| ws.rank(v)).collect();
    Ok((lists, dist, rank))
}

/// Number of shortest paths (node sequences) from `s` to every node: 1 for
/// `s`, 0 for unreachable nodes. See the module documentation.
pub fn count_shortest_paths(pg: &PathGraph, s: u32, tie: Tie) -> Result<Vec<f64>> {
    pg.check_node(s as usize, "source")?;
    let n = pg.node_count();
    let (lists, _, rank) = dag(pg, s, tie)?;
    let mut order: Vec<u32> = (0..n as u32)
        .filter(|&v| rank[v as usize] != u32::MAX)
        .collect();
    order.sort_by_key(|&v| rank[v as usize]);
    let mut sigma = vec![0.0; n];
    sigma[s as usize] = 1.0;
    for &w in &order {
        let w = w as usize;
        if w == s as usize {
            continue;
        }
        sigma[w] = lists
            .get(w)
            .iter()
            .filter(|&&v| rank[v as usize] < rank[w])
            .map(|&v| sigma[v as usize])
            .sum();
    }
    Ok(sigma)
}

/// Every shortest path from `s` to `t` with their common length, in
/// lexicographic order; empty when `t` is unreachable. Fails when there are
/// more than `max_paths` of them.
pub fn all_shortest_paths(
    pg: &PathGraph,
    s: u32,
    t: u32,
    tie: Tie,
    max_paths: usize,
) -> Result<NodeLists> {
    pg.check_node(s as usize, "source")?;
    pg.check_node(t as usize, "target")?;
    let n = pg.node_count();
    let (lists, dist, _) = dag(pg, s, tie)?;
    let mut out = NodeLists::new();
    if dist[t as usize].is_infinite() {
        return Ok(out);
    }
    // Depth-first search from t back to s over predecessors, keeping the
    // current path simple.
    let mut paths: Vec<Vec<u32>> = Vec::new();
    let mut on_path = vec![false; n];
    let mut stack: Vec<(u32, usize)> = vec![(t, 0)];
    on_path[t as usize] = true;
    while let Some(top) = stack.len().checked_sub(1) {
        let (v, i) = stack[top];
        if v == s {
            let mut p: Vec<u32> = stack.iter().rev().map(|&(x, _)| x).collect();
            p.shrink_to_fit();
            paths.push(p);
            if paths.len() > max_paths {
                return Err(Error::InvalidArgument(format!(
                    "there are more than {max_paths} shortest paths from node {s} to node {t}"
                )));
            }
            on_path[v as usize] = false;
            stack.pop();
            continue;
        }
        let preds = lists.get(v as usize);
        if i < preds.len() {
            stack[top].1 = i + 1;
            let u = preds[i];
            if !on_path[u as usize] {
                on_path[u as usize] = true;
                stack.push((u, 0));
            }
        } else {
            on_path[v as usize] = false;
            stack.pop();
        }
    }
    paths.sort();
    for p in &paths {
        out.push_with(p, dist[t as usize]);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::paths::PathGraphOptions;

    #[test]
    fn grid_paths() {
        // A 3 x 3 grid has C(4, 2) = 6 shortest paths between corners.
        let mut s = Vec::new();
        let mut t = Vec::new();
        for r in 0..3i64 {
            for c in 0..3i64 {
                let v = r * 3 + c;
                if c < 2 {
                    s.push(v);
                    t.push(v + 1);
                }
                if r < 2 {
                    s.push(v);
                    t.push(v + 3);
                }
            }
        }
        // A parallel edge does not add node sequences.
        s.push(0);
        t.push(1);
        let g = Graph::from_edges(9, &s, &t, None, false).unwrap();
        let pg = PathGraph::new(&g, PathGraphOptions::default()).unwrap();
        let p = all_shortest_paths(&pg, 0, 8, Tie::Relative, 100).unwrap();
        assert_eq!(p.len(), 6);
        assert_eq!(p.get(0), &[0, 1, 2, 5, 8]);
        assert_eq!(p.values[0], 4.0);
        let sigma = count_shortest_paths(&pg, 0, Tie::Relative).unwrap();
        assert_eq!(sigma[8], 6.0);
        assert_eq!(sigma[4], 2.0);
        assert!(all_shortest_paths(&pg, 0, 8, Tie::Relative, 5).is_err());
    }
}
