//! Edge-disjoint and node-disjoint paths between two nodes.
//!
//! By Menger's theorem (Menger 1927, `menger1927allgemeinen`) the largest
//! number of edge-disjoint s-t paths is lambda(s, t) and the largest
//! number of internally node-disjoint ones is kappa(s, t)
//! ([`super::values`]). A maximum unit-capacity flow in the edge network
//! (or the split network; Even and Tarjan 1975, `even1975network`) is
//! decomposed into that many paths: starting from s, follow unused arcs
//! that carry flow (in the order of the arcs, by increasing neighbour),
//! cutting out any cycle met on the way, until t is reached (Ford and
//! Fulkerson 1962, `ford1962flows`).
//!
//! Conventions: simple graph, directions kept (see [`super`]); an edge from
//! s to t is one of the node-disjoint paths. Paths are node sequences from
//! s to t. When no path exists the result is empty (NetworkX raises
//! `NetworkXNoPath`). The number of paths always equals NetworkX's
//! (`edge_disjoint_paths`, `node_disjoint_paths`), but a maximum flow is
//! not unique, so the paths themselves can differ from NetworkX's, which
//! come from its Edmonds-Karp flow.

use super::{check_pair, Simple};
use crate::error::Result;
use crate::graph::Graph;

/// Decomposes a unit flow given as arcs (tail, head) into paths from `s`
/// to `t`. `arcs` lists every arc carrying one unit, grouped by tail in
/// the order they should be tried.
fn decompose(nodes: usize, arcs: &[(usize, usize)], s: usize, t: usize) -> Vec<Vec<usize>> {
    let mut offsets = vec![0usize; nodes + 1];
    for &(u, _) in arcs {
        offsets[u + 1] += 1;
    }
    for i in 0..nodes {
        offsets[i + 1] += offsets[i];
    }
    let mut heads = vec![0usize; arcs.len()];
    let mut cursor = offsets.clone();
    for &(u, v) in arcs {
        heads[cursor[u]] = v;
        cursor[u] += 1;
    }
    let mut next = offsets.clone();
    let mut on_path = vec![usize::MAX; nodes];
    let mut paths = Vec::new();
    while next[s] < offsets[s + 1] {
        let mut path = vec![s];
        on_path[s] = 0;
        let mut u = s;
        while u != t {
            // Conservation guarantees an unused arc out of every node the
            // walk enters.
            let v = heads[next[u]];
            next[u] += 1;
            if on_path[v] != usize::MAX {
                // Cut out the cycle v ... u.
                for &x in &path[on_path[v] + 1..] {
                    on_path[x] = usize::MAX;
                }
                path.truncate(on_path[v] + 1);
            } else {
                on_path[v] = path.len();
                path.push(v);
            }
            u = v;
        }
        for &x in &path {
            on_path[x] = usize::MAX;
        }
        paths.push(path);
    }
    paths
}

/// Edge-disjoint paths from `s` to `t`, as many as possible. See the
/// module documentation.
pub fn edge_disjoint_paths(g: &Graph, s: usize, t: usize) -> Result<Vec<Vec<usize>>> {
    check_pair(s, t, g.node_count())?;
    let simple = Simple::new(g);
    let mut net = simple.edge_network();
    net.max_flow(s, t, None)?;
    let mut arcs = Vec::new();
    for (k, (u, v)) in simple.arcs().into_iter().enumerate() {
        let f = net.flow(k);
        if f > 0.5 {
            arcs.push((u as usize, v as usize));
        } else if f < -0.5 {
            arcs.push((v as usize, u as usize));
        }
    }
    arcs.sort_unstable();
    Ok(decompose(simple.n, &arcs, s, t))
}

/// Internally node-disjoint paths from `s` to `t`, as many as possible.
/// See the module documentation.
pub fn node_disjoint_paths(g: &Graph, s: usize, t: usize) -> Result<Vec<Vec<usize>>> {
    check_pair(s, t, g.node_count())?;
    let simple = Simple::new(g);
    let n = simple.n;
    let (mut net, split) = simple.node_network(1.0);
    net.max_flow(n + s, t, None)?;
    // Arcs between nodes carrying flow, as (u, v) in the original
    // network; flow through a node is implied.
    let mut arcs: Vec<(usize, usize)> = split
        .iter()
        .enumerate()
        .filter(|&(i, _)| net.flow(n + i) > 0.5)
        .map(|(_, &(u, v))| (u as usize, v as usize))
        .collect();
    arcs.sort_unstable();
    Ok(decompose(n, &arcs, s, t))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn bowtie() {
        // Two triangles sharing node 2: two edge-disjoint paths from 0 to
        // 4 through 2, but only one node-disjoint path.
        let g = graph(5, &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)], false);
        let e = edge_disjoint_paths(&g, 0, 4).unwrap();
        assert_eq!(e.len(), 2);
        assert!(e.iter().all(|p| p[0] == 0 && *p.last().unwrap() == 4));
        let v = node_disjoint_paths(&g, 0, 4).unwrap();
        assert_eq!(v.len(), 1);
        // Adjacent nodes: the edge itself is a path.
        let v = node_disjoint_paths(&g, 0, 1).unwrap();
        assert_eq!(v.len(), 2);
        assert!(v.contains(&vec![0, 1]));
        let d = graph(3, &[(0, 1)], true);
        assert!(edge_disjoint_paths(&d, 1, 0).unwrap().is_empty());
    }
}
