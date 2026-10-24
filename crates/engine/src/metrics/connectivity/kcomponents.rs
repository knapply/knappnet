//! k-components: maximal sets of nodes inducing a k-connected subgraph.
//!
//! Definition (Moody and White 2003, `moody2003structural`): a
//! k-component is a maximal set of nodes whose induced subgraph has node
//! connectivity at least k ([`super::values`]), so it has more than k
//! nodes. The k-components of a network form a hierarchy: each
//! (k + 1)-component lies inside a k-component, and two k-components share
//! fewer than k nodes (their union would otherwise be k-connected). The
//! 1-components are the connected components with at least two nodes and
//! the 2-components the biconnected components (blocks) with at least
//! three nodes, so a bridge is a block but not a 2-component.
//!
//! Algorithm (the recursion of Moody and White 2003): every k-connected
//! node set lies inside a (k - 1)-component, and if a node set H has a
//! node cut S with fewer than k nodes, every k-connected subset of H lies
//! in one piece C_i of H - S together with S (removing fewer than k nodes
//! leaves it connected). So the k-components inside H are H itself when
//! kappa(H) >= k, and otherwise the maximal sets among those found
//! recursively in each C_i together with S. Node connectivity and minimum
//! node cuts of induced subgraphs come from [`super::values`] and
//! [`super::cuts`]. Each piece is smaller than H, so the recursion ends;
//! its cost grows with the number of nested cuts and is not polynomial in
//! the worst case.
//!
//! Conventions: undirected networks only (NetworkX's `k_components`
//! refuses directed graphs); the simple graph is used. The result lists,
//! for k = 1, 2, ... up to the largest k with a k-component, each
//! k-component as a sorted node list; components of the same k are in
//! lexicographic order. These are exactly the maximal k-connected node
//! sets, which NetworkX's `k_components` (`hagberg2008exploring`), an
//! implementation of the same method by Torrents and Ferraro (2015,
//! `torrents2015structural`), is meant to return.

use super::cuts::{minimum_node_cut, CutConvention};
use super::values::node_connectivity;
use super::Simple;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::factions::bicomponents;

/// Builds a graph from a simple undirected adjacency.
fn to_graph(s: &Simple) -> Graph {
    let mut a = Vec::new();
    let mut b = Vec::new();
    for u in 0..s.n {
        for &v in &s.out[u] {
            if (u as u32) < v {
                a.push(u as i64);
                b.push(v as i64);
            }
        }
    }
    Graph::from_edges(s.n, &a, &b, None, false).expect("valid induced subgraph")
}

/// Connected components of `s` restricted to the nodes not in `removed`,
/// as sorted lists of indices into `s`.
fn pieces(s: &Simple, removed: &[bool]) -> Vec<Vec<u32>> {
    let mut seen = removed.to_vec();
    let mut out = Vec::new();
    for r in 0..s.n {
        if seen[r] {
            continue;
        }
        seen[r] = true;
        let mut comp = vec![r as u32];
        let mut stack = vec![r];
        while let Some(v) = stack.pop() {
            for &w in &s.out[v] {
                if !seen[w as usize] {
                    seen[w as usize] = true;
                    comp.push(w);
                    stack.push(w as usize);
                }
            }
        }
        comp.sort_unstable();
        out.push(comp);
    }
    out
}

/// Keeps the sets not strictly contained in another (and one copy of
/// each), sorted.
fn maximal(mut sets: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    sets.sort();
    sets.dedup();
    let keep: Vec<bool> = (0..sets.len())
        .map(|i| {
            !(0..sets.len()).any(|j| {
                j != i
                    && sets[j].len() > sets[i].len()
                    && sets[i].iter().all(|x| sets[j].binary_search(x).is_ok())
            })
        })
        .collect();
    sets.into_iter()
        .zip(keep)
        .filter(|(_, k)| *k)
        .map(|(s, _)| s)
        .collect()
}

/// The maximal k-connected subsets of `nodes` (sorted global indices),
/// each with its own node connectivity.
fn within(full: &Simple, nodes: &[u32], k: usize, threads: usize) -> Vec<(Vec<u32>, usize)> {
    if nodes.len() <= k {
        return Vec::new();
    }
    let sub = full.induced(nodes);
    let g = to_graph(&sub);
    let kappa = node_connectivity(&g, threads).expect("non-empty") as usize;
    if kappa >= k {
        return vec![(nodes.to_vec(), kappa)];
    }
    let cut = minimum_node_cut(&g, CutConvention::Engine, threads).expect("non-empty");
    let mut removed = vec![false; sub.n];
    for &c in &cut {
        removed[c as usize] = true;
    }
    let mut found: Vec<(Vec<u32>, usize)> = Vec::new();
    for piece in pieces(&sub, &removed) {
        let mut local: Vec<u32> = piece.iter().chain(&cut).copied().collect();
        local.sort_unstable();
        let global: Vec<u32> = local.iter().map(|&i| nodes[i as usize]).collect();
        found.extend(within(full, &global, k, threads));
    }
    let sets = maximal(found.iter().map(|f| f.0.clone()).collect());
    sets.into_iter()
        .map(|s| {
            let c = found.iter().find(|f| f.0 == s).expect("found").1;
            (s, c)
        })
        .collect()
}

/// Every k-component, by k: element k - 1 holds the k-components. See the
/// module documentation.
pub fn k_components(g: &Graph, threads: usize) -> Result<Vec<Vec<Vec<u32>>>> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(
            "k-components are defined for undirected networks".into(),
        ));
    }
    let full = Simple::new(g);
    let mut levels: Vec<Vec<Vec<u32>>> = Vec::new();
    // k = 1: components with at least two nodes.
    let ones: Vec<Vec<u32>> = pieces(&full, &vec![false; full.n])
        .into_iter()
        .filter(|c| c.len() > 1)
        .collect();
    if ones.is_empty() {
        return Ok(levels);
    }
    let mut ones = ones;
    ones.sort();
    levels.push(ones);
    // k = 2: blocks with at least three nodes.
    let mut current: Vec<(Vec<u32>, usize)> = bicomponents(g)
        .into_iter()
        .filter(|b| b.len() > 2)
        // Their connectivity is at least 2; `within` finds the exact value.
        .map(|b| (b, 2))
        .collect();
    let mut k = 2;
    while !current.is_empty() {
        let mut sets: Vec<Vec<u32>> = current.iter().map(|c| c.0.clone()).collect();
        sets.sort();
        levels.push(sets);
        k += 1;
        let mut next: Vec<(Vec<u32>, usize)> = Vec::new();
        for (nodes, kappa) in &current {
            if *kappa >= k {
                next.push((nodes.clone(), *kappa));
            } else {
                next.extend(within(&full, nodes, k, threads));
            }
        }
        // Components of different parents are distinct maximal sets, but a
        // set found in two parents (sharing it) appears once.
        let mut seen: Vec<Vec<u32>> = Vec::new();
        next.retain(|(s, _)| {
            if seen.contains(s) {
                false
            } else {
                seen.push(s.clone());
                true
            }
        });
        current = next;
    }
    Ok(levels)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn nested_cliques() {
        // K4 on 0..4 sharing node 3 with a triangle 3-4-5, a pendant 5-6,
        // and an isolated node 7.
        let mut e = vec![];
        for a in 0..4 {
            for b in (a + 1)..4 {
                e.push((a, b));
            }
        }
        e.extend([(3, 4), (4, 5), (5, 3), (5, 6)]);
        let g = graph(8, &e);
        let k = k_components(&g, 1).unwrap();
        assert_eq!(k.len(), 3);
        assert_eq!(k[0], vec![vec![0, 1, 2, 3, 4, 5, 6]]);
        assert_eq!(k[1], vec![vec![0, 1, 2, 3], vec![3, 4, 5]]);
        assert_eq!(k[2], vec![vec![0, 1, 2, 3]]);
        assert!(k_components(&graph(3, &[]), 1).unwrap().is_empty());
    }
}
