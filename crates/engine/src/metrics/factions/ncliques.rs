//! n-cliques and n-clans.
//!
//! An n-clique is a maximal set of nodes every two of which are joined by
//! a path of length at most n in the whole network (Luce 1950,
//! `luce1950connectivity`; Alba 1973, `alba1973graph`). The paths may pass
//! through nodes outside the set. n-cliques are therefore exactly the
//! maximal cliques of the n-th power of the network, the graph joining two
//! nodes when their distance is at most n, and are listed by running
//! [`crate::metrics::cliques::max_cliques`] on that graph. 1-cliques are
//! cliques.
//!
//! An n-clan is an n-clique whose own subgraph has diameter at most n, so
//! that every pair is joined by a short path inside the group (Mokken
//! 1979, `mokken1979cliques`). n-clans are found by keeping the n-cliques
//! that pass this test.
//!
//! Conventions:
//! - Distances are counted in the undirected simple reading of the network
//!   (see [`crate::metrics::factions`]), without weights, as in UCINET's
//!   n-clique and n-clan routines for symmetric data.
//! - `min_size` drops groups with fewer members; it never reports
//!   non-maximal subsets. With `min_size` at most 1, an isolated node is an
//!   n-clique (and an n-clan) of its own, as for maximal cliques.
//! - Members are listed in increasing order and groups in lexicographic
//!   order. The power graph is built with one breadth-first search per
//!   node, in parallel.

use super::simple_rows;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::cliques::max_cliques;
use crate::parallel;

fn check_order(order: usize) -> Result<()> {
    if order == 0 {
        return Err(Error::InvalidArgument(
            "n (the largest distance within a group) must be at least 1".into(),
        ));
    }
    Ok(())
}

/// Nodes within `order` steps of `source` (excluding it), by breadth-first
/// search restricted to nodes with `inside[v]` (all nodes when `None`).
fn within(
    rows: &[Vec<u32>],
    source: usize,
    order: usize,
    inside: Option<&[bool]>,
    dist: &mut [usize],
    queue: &mut Vec<u32>,
) -> Vec<u32> {
    queue.clear();
    queue.push(source as u32);
    dist[source] = 0;
    let mut head = 0;
    while head < queue.len() {
        let v = queue[head] as usize;
        head += 1;
        if dist[v] == order {
            continue;
        }
        for &w in &rows[v] {
            let w = w as usize;
            if dist[w] == usize::MAX && inside.map_or(true, |s| s[w]) {
                dist[w] = dist[v] + 1;
                queue.push(w as u32);
            }
        }
    }
    let found: Vec<u32> = queue[1..].to_vec();
    for &v in queue.iter() {
        dist[v as usize] = usize::MAX;
    }
    found
}

/// The n-th power of the undirected simple reading of `g`.
fn power_graph(g: &Graph, rows: &[Vec<u32>], order: usize, threads: usize) -> Result<Graph> {
    let n = g.node_count();
    let parts = parallel::for_each_item(
        n,
        threads,
        || (vec![usize::MAX; n], Vec::new(), Vec::new()),
        |(dist, queue, edges): &mut (Vec<usize>, Vec<u32>, Vec<(u32, u32)>), v| {
            for w in within(rows, v, order, None, dist, queue) {
                if (w as usize) > v {
                    edges.push((v as u32, w));
                }
            }
        },
    );
    let mut edges: Vec<(u32, u32)> = parts.into_iter().flat_map(|p| p.2).collect();
    edges.sort_unstable();
    let src: Vec<i64> = edges.iter().map(|e| e.0 as i64).collect();
    let tgt: Vec<i64> = edges.iter().map(|e| e.1 as i64).collect();
    Graph::from_edges(n, &src, &tgt, None, false)
}

/// n-cliques of `g` with at least `min_size` members. See the module
/// documentation.
pub fn n_cliques(
    g: &Graph,
    order: usize,
    min_size: usize,
    threads: usize,
) -> Result<Vec<Vec<u32>>> {
    check_order(order)?;
    let rows = simple_rows(g);
    let power = power_graph(g, &rows, order, threads)?;
    Ok(max_cliques(&power, min_size, threads))
}

/// Whether the subgraph induced by `group` has diameter at most `order`.
fn small_diameter(rows: &[Vec<u32>], group: &[u32], order: usize, n: usize) -> bool {
    let mut inside = vec![false; n];
    for &v in group {
        inside[v as usize] = true;
    }
    let mut dist = vec![usize::MAX; n];
    let mut queue = Vec::new();
    group.iter().all(|&v| {
        within(
            rows,
            v as usize,
            order,
            Some(&inside),
            &mut dist,
            &mut queue,
        )
        .len()
            + 1
            == group.len()
    })
}

/// n-clans of `g` with at least `min_size` members. See the module
/// documentation.
pub fn n_clans(g: &Graph, order: usize, min_size: usize, threads: usize) -> Result<Vec<Vec<u32>>> {
    check_order(order)?;
    let rows = simple_rows(g);
    let power = power_graph(g, &rows, order, threads)?;
    let n = g.node_count();
    let cliques = max_cliques(&power, min_size, threads);
    let keep = parallel::for_each_item(cliques.len(), threads, Vec::new, |acc, i| {
        if small_diameter(&rows, &cliques[i], order, n) {
            acc.push(i);
        }
    });
    let mut keep: Vec<usize> = keep.into_iter().flatten().collect();
    keep.sort_unstable();
    Ok(keep.into_iter().map(|i| cliques[i].clone()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn one_cliques_are_cliques() {
        let g = und(5, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
        assert_eq!(n_cliques(&g, 1, 0, 1).unwrap(), max_cliques(&g, 0, 1));
        assert!(n_cliques(&g, 0, 0, 1).is_err());
    }

    #[test]
    fn path_and_alba_example() {
        // On a path 0-1-2-3-4, the 2-cliques are the windows of three.
        let p = und(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert_eq!(
            n_cliques(&p, 2, 0, 2).unwrap(),
            vec![vec![0, 1, 2], vec![1, 2, 3], vec![2, 3, 4]]
        );
        // A 6-cycle with a pendant: {0, 1, 2, 4} style sets. In a hexagon
        // 0..5, nodes at distance <= 2 of each other form the 2-cliques
        // {i, i+1, i+2} and {0, 2, 4}, {1, 3, 5}.
        let hex = und(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]);
        let two = n_cliques(&hex, 2, 0, 1).unwrap();
        assert!(two.contains(&vec![0, 2, 4]) && two.contains(&vec![1, 3, 5]));
        assert_eq!(two.len(), 8);
        // {0, 2, 4} is a 2-clique but its own subgraph has no edges, so it
        // is not a 2-clan; the windows of three are.
        let clans = n_clans(&hex, 2, 0, 1).unwrap();
        assert_eq!(clans.len(), 6);
        assert!(!clans.contains(&vec![0, 2, 4]));
        assert_eq!(n_clans(&hex, 2, 4, 1).unwrap(), Vec::<Vec<u32>>::new());
    }
}
