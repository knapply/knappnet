//! Maximal cliques and the clique number.
//!
//! A clique is a set of nodes that are all pairwise adjacent; it is maximal
//! when no other node is adjacent to all of its members. Maximal cliques are
//! listed with the Bron-Kerbosch algorithm (Bron and Kerbosch 1973,
//! `bron1973algorithm`) with the pivot rule of Tomita, Tanaka, and
//! Takahashi (2006, `tomita2006worst`): at each step the pivot is the
//! candidate or excluded node with the most neighbors among the candidates.
//! The outer level visits nodes in a degeneracy (smallest-last) ordering
//! (Matula and Beck 1983, `matula1983smallest`), as proposed by Eppstein,
//! Loffler, and Strash (2010, `eppstein2010listing`), which bounds the work
//! by O(d n 3^(d/3)) for a graph of degeneracy d. The ordering comes from
//! the core decomposition in [`crate::metrics::cores`].
//!
//! Conventions (matching igraph's `max_cliques()` and NetworkX's
//! `find_cliques()`):
//! - Cliques are found in the undirected simple graph underlying the
//!   network: edge directions are ignored (an arc in either direction makes
//!   two nodes adjacent), parallel edges count once, and self-loops are
//!   ignored. NetworkX's `find_cliques` rejects directed graphs, so its
//!   results are compared on the same projection.
//! - Weights are ignored: cliques are defined on binary ties.
//! - An isolated node is a maximal clique of size 1. A graph with no nodes
//!   has no cliques and clique number 0.
//! - Each clique lists its members in increasing node order, and cliques
//!   are returned in lexicographic order, so results do not depend on the
//!   thread count.
//! - `min_size` drops maximal cliques with fewer members (as igraph's
//!   `min`); it does not report non-maximal subsets of larger cliques.

use crate::graph::{Csr, Graph};
use crate::metrics::cores;
use crate::parallel;

/// The order in which nodes are peeled from the simple graph `adj`: a
/// degeneracy ordering, so each node has few later neighbors.
fn degeneracy_order(adj: &Csr) -> Vec<u32> {
    let n = adj.node_count();
    let deg = (0..n).map(|v| adj.neighbors(v).len()).collect();
    cores::peel(deg, |v, buf| {
        buf.clear();
        buf.extend_from_slice(adj.neighbors(v));
    })
    .1
}

/// Elements of sorted `a` that are also in sorted `b`.
fn intersect(a: &[u32], b: &[u32], out: &mut Vec<u32>) {
    out.clear();
    if a.is_empty() || b.is_empty() {
        return;
    }
    // Binary search when one side is much shorter (for example a small
    // candidate set against a hub's neighbors).
    if a.len() * 8 < b.len() {
        out.extend(a.iter().copied().filter(|x| b.binary_search(x).is_ok()));
        return;
    }
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
}

fn count_common(a: &[u32], b: &[u32]) -> usize {
    if a.len() * 8 < b.len() {
        return a.iter().filter(|x| b.binary_search(x).is_ok()).count();
    }
    let (mut i, mut j, mut c) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                c += 1;
                i += 1;
                j += 1;
            }
        }
    }
    c
}

/// What to do with each maximal clique found.
trait Sink {
    /// Called with every maximal clique (in no particular member order).
    fn report(&mut self, clique: &[u32]);
    /// Branches whose largest possible clique is smaller than this are cut.
    fn bound(&self) -> usize;
}

struct Collect {
    min_size: usize,
    cliques: Vec<Vec<u32>>,
}

impl Sink for Collect {
    fn report(&mut self, clique: &[u32]) {
        if clique.len() >= self.min_size {
            let mut c = clique.to_vec();
            c.sort_unstable();
            self.cliques.push(c);
        }
    }
    fn bound(&self) -> usize {
        self.min_size
    }
}

struct Largest(usize);

impl Sink for Largest {
    fn report(&mut self, clique: &[u32]) {
        self.0 = self.0.max(clique.len());
    }
    fn bound(&self) -> usize {
        // Only strictly larger cliques can change the answer.
        self.0 + 1
    }
}

/// Bron-Kerbosch with Tomita pivoting. `p` and `x` are sorted.
fn expand<S: Sink>(adj: &Csr, r: &mut Vec<u32>, mut p: Vec<u32>, mut x: Vec<u32>, sink: &mut S) {
    if p.is_empty() {
        if x.is_empty() {
            sink.report(r);
        }
        return;
    }
    if r.len() + p.len() < sink.bound() {
        return;
    }
    let pivot = p
        .iter()
        .chain(x.iter())
        .copied()
        .max_by_key(|&u| {
            (
                count_common(&p, adj.neighbors(u as usize)),
                std::cmp::Reverse(u),
            )
        })
        .expect("p is not empty");
    let mut branch = Vec::new();
    intersect(&p, adj.neighbors(pivot as usize), &mut branch);
    // Candidates not adjacent to the pivot: p minus (p ∩ N(pivot)).
    let skip = branch;
    let todo: Vec<u32> = p
        .iter()
        .copied()
        .filter(|w| skip.binary_search(w).is_err())
        .collect();
    let mut np = Vec::new();
    let mut nx = Vec::new();
    for w in todo {
        let nw = adj.neighbors(w as usize);
        intersect(&p, nw, &mut np);
        intersect(&x, nw, &mut nx);
        r.push(w);
        expand(
            adj,
            r,
            std::mem::take(&mut np),
            std::mem::take(&mut nx),
            sink,
        );
        r.pop();
        if r.len() + p.len() < sink.bound() {
            // Removing w only shrinks p further; nothing left can qualify.
            return;
        }
        let at = p.binary_search(&w).expect("w is a candidate");
        p.remove(at);
        let at = x.binary_search(&w).unwrap_err();
        x.insert(at, w);
    }
}

/// Runs the outer, degeneracy-ordered loop on node `order[i]`.
fn from_node<S: Sink>(adj: &Csr, order: &[u32], rank: &[usize], i: usize, sink: &mut S) {
    let v = order[i] as usize;
    let mut p = Vec::new();
    let mut x = Vec::new();
    for &u in adj.neighbors(v) {
        if rank[u as usize] > i {
            p.push(u);
        } else {
            x.push(u);
        }
    }
    let mut r = vec![v as u32];
    expand(adj, &mut r, p, x, sink);
}

fn setup(g: &Graph) -> (Csr, Vec<u32>, Vec<usize>) {
    let adj = g.simple_undirected(false);
    let order = degeneracy_order(&adj);
    let mut rank = vec![0usize; order.len()];
    for (i, &v) in order.iter().enumerate() {
        rank[v as usize] = i;
    }
    (adj, order, rank)
}

/// All maximal cliques with at least `min_size` members, each sorted, in
/// lexicographic order. See the module documentation for conventions.
/// `threads == 0` uses all available cores.
pub fn max_cliques(g: &Graph, min_size: usize, threads: usize) -> Vec<Vec<u32>> {
    let (adj, order, rank) = setup(g);
    let parts = parallel::for_each_item(
        order.len(),
        threads,
        || Collect {
            min_size,
            cliques: Vec::new(),
        },
        |sink, i| from_node(&adj, &order, &rank, i, sink),
    );
    let mut all: Vec<Vec<u32>> = parts.into_iter().flat_map(|s| s.cliques).collect();
    all.sort_unstable();
    all
}

/// Size of the largest clique (0 for a graph with no nodes).
pub fn clique_number(g: &Graph, threads: usize) -> usize {
    let (adj, order, rank) = setup(g);
    // Later nodes in the ordering have smaller forward neighborhoods, so a
    // sequential pass benefits most from the bound; threads each keep their
    // own bound.
    let parts = parallel::for_each_item(
        order.len(),
        threads,
        || Largest(0),
        |sink, i| from_node(&adj, &order, &rank, i, sink),
    );
    parts.into_iter().map(|s| s.0).max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn hand_checked() {
        // Two triangles sharing node 2, a pendant 5, and an isolate 6.
        let g = graph(
            7,
            &[(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (2, 4), (4, 5)],
            false,
        );
        assert_eq!(
            max_cliques(&g, 0, 1),
            vec![vec![0, 1, 2], vec![2, 3, 4], vec![4, 5], vec![6]]
        );
        assert_eq!(max_cliques(&g, 3, 2), vec![vec![0, 1, 2], vec![2, 3, 4]]);
        assert_eq!(clique_number(&g, 1), 3);
    }

    #[test]
    fn projection_ignores_direction_loops_and_multiplicity() {
        let g = graph(3, &[(0, 1), (1, 0), (1, 2), (2, 0), (2, 2), (0, 1)], true);
        assert_eq!(max_cliques(&g, 0, 1), vec![vec![0, 1, 2]]);
        assert_eq!(clique_number(&g, 1), 3);
    }

    #[test]
    fn complete_and_empty() {
        let mut e = Vec::new();
        for i in 0..6 {
            for j in (i + 1)..6 {
                e.push((i, j));
            }
        }
        let k6 = graph(6, &e, false);
        assert_eq!(max_cliques(&k6, 0, 0), vec![vec![0, 1, 2, 3, 4, 5]]);
        assert_eq!(clique_number(&k6, 0), 6);
        let none = graph(0, &[], false);
        assert!(max_cliques(&none, 0, 1).is_empty());
        assert_eq!(clique_number(&none, 1), 0);
        let isolates = graph(3, &[], false);
        assert_eq!(max_cliques(&isolates, 0, 1).len(), 3);
        assert_eq!(clique_number(&isolates, 1), 1);
        assert!(max_cliques(&isolates, 2, 1).is_empty());
    }

    #[test]
    fn complete_tripartite_graph() {
        // K(3,3,3), the complement of three disjoint triangles: a maximal
        // clique takes one node from each part, so there are 3^3 = 27.
        let mut e = Vec::new();
        for i in 0..9i64 {
            for j in (i + 1)..9 {
                if i / 3 != j / 3 {
                    e.push((i, j));
                }
            }
        }
        let g = graph(9, &e, false);
        let c = max_cliques(&g, 0, 3);
        assert_eq!(c.len(), 27);
        assert!(c.iter().all(|k| k.len() == 3));
        assert_eq!(clique_number(&g, 3), 3);
    }
}
