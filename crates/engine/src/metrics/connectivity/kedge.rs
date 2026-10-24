//! k-edge-connected components: the classes of nodes joined by at least k
//! edge-disjoint paths.
//!
//! Definition: nodes u and v are k-edge-connected when lambda(u, v) >= k
//! (in both directions, for a directed network; see [`super::values`]).
//! Since lambda(u, w) >= min(lambda(u, v), lambda(v, w)), this is an
//! equivalence relation, and its classes are the k-edge-connected
//! components (NetworkX's `k_edge_components`, which follows Wang et al.
//! 2015, `wang2015simple`; `hagberg2008exploring`). They partition the
//! nodes, singletons included. Note that a class need not induce a
//! k-edge-connected subgraph: the paths may leave it.
//!
//! Algorithm: for an undirected network, lambda(u, v) is the smallest
//! weight on the path between u and v in a Gomory-Hu tree
//! ([`crate::flow::gomory_hu`]; Gomory and Hu 1961, `gomory1961multi`), so
//! the classes are the components of the tree after dropping its edges of
//! weight below k: n - 1 maximum flows in all. For a directed network,
//! classes are built one at a time: the lowest unassigned node r takes
//! every unassigned node v of its strongly connected component with
//! lambda(r, v) >= k and lambda(v, r) >= k (by transitivity these are
//! exactly r's class), each flow stopping at k.
//!
//! Conventions: simple graph, directions kept, k >= 1 (k = 1 gives the
//! connected or strongly connected components). Labels number the classes
//! from 0 in order of their lowest node.

use super::values::lambda;
use super::Simple;
use crate::error::{Error, Result};
use crate::flow::gomory_hu_tree;
use crate::graph::Graph;
use crate::metrics::connectivity::components::strong_components;
use crate::order::NeighborOrder;
use crate::parallel;

/// Relabels `raw` labels in order of first appearance.
pub(crate) fn first_appearance(raw: &[usize]) -> Vec<u32> {
    let mut map = std::collections::HashMap::new();
    raw.iter()
        .map(|&r| {
            let next = map.len() as u32;
            *map.entry(r).or_insert(next)
        })
        .collect()
}

/// Class label of every node. See the module documentation.
pub fn k_edge_components(g: &Graph, k: u64, threads: usize) -> Result<Vec<u32>> {
    if k == 0 {
        return Err(Error::InvalidArgument("k must be at least 1".into()));
    }
    let simple = Simple::new(g);
    let n = simple.n;
    if n == 0 {
        return Ok(Vec::new());
    }
    if !simple.directed {
        let arcs = simple.arcs();
        let a: Vec<i64> = arcs.iter().map(|x| x.0 as i64).collect();
        let b: Vec<i64> = arcs.iter().map(|x| x.1 as i64).collect();
        let sg = Graph::from_edges(n, &a, &b, None, false)?;
        let tree = gomory_hu_tree(&sg, None)?;
        let mut sets = crate::union_find::UnionFind::new(n);
        for v in 1..n {
            if tree.weight[v] >= k as f64 {
                sets.union(v, tree.parent[v]);
            }
        }
        let roots: Vec<usize> = (0..n).map(|v| sets.find(v)).collect();
        return Ok(first_appearance(&roots));
    }
    let scc = strong_components(g, NeighborOrder::Sorted);
    if k == 1 {
        let raw: Vec<usize> = scc.iter().map(|&c| c as usize).collect();
        return Ok(first_appearance(&raw));
    }
    let base = simple.edge_network();
    let mut class = vec![usize::MAX; n];
    for r in 0..n {
        if class[r] != usize::MAX {
            continue;
        }
        class[r] = r;
        let cands: Vec<usize> = ((r + 1)..n)
            .filter(|&v| class[v] == usize::MAX && scc[v] == scc[r])
            .collect();
        let parts = parallel::for_each_item(
            cands.len(),
            threads,
            || (base.clone(), Vec::new()),
            |(net, ok): &mut (crate::flow::FlowNetwork, Vec<usize>), i| {
                let v = cands[i];
                let limit = Some(k as f64);
                if lambda(net, r, v, limit) >= k && lambda(net, v, r, limit) >= k {
                    ok.push(v);
                }
            },
        );
        for (_, ok) in parts {
            for v in ok {
                class[v] = r;
            }
        }
    }
    Ok(first_appearance(&class))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn classes() {
        // Two triangles joined by a bridge, plus an isolate.
        let g = graph(
            7,
            &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)],
            false,
        );
        assert_eq!(
            k_edge_components(&g, 1, 1).unwrap(),
            vec![0, 0, 0, 0, 0, 0, 1]
        );
        assert_eq!(
            k_edge_components(&g, 2, 1).unwrap(),
            vec![0, 0, 0, 1, 1, 1, 2]
        );
        assert_eq!(
            k_edge_components(&g, 3, 1).unwrap(),
            vec![0, 1, 2, 3, 4, 5, 6]
        );
        // Directed: a complete digraph on 0, 1, 2 and an arc to 3.
        let d = graph(
            4,
            &[(0, 1), (1, 0), (1, 2), (2, 1), (0, 2), (2, 0), (2, 3)],
            true,
        );
        assert_eq!(k_edge_components(&d, 2, 2).unwrap(), vec![0, 0, 0, 1]);
        assert!(k_edge_components(&d, 0, 1).is_err());
    }
}
