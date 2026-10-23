//! Girvan-Newman divisive clustering.
//!
//! Girvan and Newman (2002, `girvan2002community`) repeatedly remove the
//! edge with the highest edge betweenness, recomputing betweenness after
//! every removal, so that the network falls apart into ever smaller
//! components. Each set of components is a candidate partition; following
//! Newman and Girvan (2004, `newman2004finding`), the one with the highest
//! modularity is chosen. Edge betweenness, the number of shortest paths
//! through an edge (Girvan and Newman 2002; Freeman 1977,
//! `freeman1977set`), is accumulated with Brandes' (2001,
//! `brandes2001faster`) algorithm, and after a removal only the component
//! that contained the edge is recomputed, since no other shortest path
//! changes.
//!
//! Conventions:
//! - Undirected networks only. The network is read as a simple graph:
//!   parallel edges count once and self-loops are ignored, and modularity
//!   is that of this simple graph (NetworkX's `girvan_newman` and
//!   `modularity` on a `Graph`). Weights are not used: betweenness counts
//!   hops.
//! - Ties: the edge removed is the one with the largest betweenness;
//!   values within a relative 1e-9 of the largest are treated as equal, and
//!   among them the edge whose endpoints (smaller node first) come first in
//!   lexicographic order is removed. igraph removes the tied edge with the
//!   lowest edge id and NetworkX the first in its edge iteration order, so
//!   on networks with tied betweenness the three can split differently.
//! - The chosen partition is the first of highest modularity in removal
//!   order (the one with the fewest groups). A network without edges has
//!   modularity NaN and every node is its own group.
//! - Cost O(m^2 n) in the worst case. Shortest-path searches from each
//!   source run in parallel; the tolerance for ties keeps the choice of
//!   edge independent of the thread count.

use std::collections::VecDeque;

use super::{as_labels, relabel, require_undirected};
use crate::error::Result;
use crate::graph::Graph;
use crate::metrics::community::modularity;
use crate::metrics::factions::simple_rows;
use crate::metrics::topography::{components, ComponentKind};
use crate::parallel;

/// Result of [`girvan_newman`].
#[derive(Debug, Clone, PartialEq)]
pub struct GirvanNewman {
    /// Group of each node in the partition of highest modularity.
    pub membership: Vec<usize>,
    /// Its modularity.
    pub modularity: f64,
    /// Edges of the simple graph in the order they were removed, smaller
    /// node first.
    pub removed: Vec<(u32, u32)>,
    /// Modularity before any removal (entry 0) and after each removal.
    pub path: Vec<f64>,
}

/// The remaining simple graph, as adjacency rows of (neighbor, edge id).
struct Remaining {
    rows: Vec<Vec<(u32, usize)>>,
    edges: Vec<(u32, u32)>,
    alive: Vec<bool>,
}

impl Remaining {
    fn component(&self, source: usize) -> Vec<usize> {
        let n = self.rows.len();
        let mut seen = vec![false; n];
        seen[source] = true;
        let mut queue = vec![source];
        let mut head = 0;
        while head < queue.len() {
            let v = queue[head];
            head += 1;
            for &(w, e) in &self.rows[v] {
                if self.alive[e] && !seen[w as usize] {
                    seen[w as usize] = true;
                    queue.push(w as usize);
                }
            }
        }
        queue
    }

    /// Edge betweenness of the edges reached from `sources`, added into
    /// `eb` (which the caller has zeroed for those edges). Each unordered
    /// pair is counted from both ends and the result halved.
    fn betweenness(&self, sources: &[usize], eb: &mut [f64], threads: usize) {
        let n = self.rows.len();
        let m = self.edges.len();
        let parts = parallel::for_each_item(
            sources.len(),
            threads,
            || {
                (
                    vec![0.0f64; m],
                    vec![usize::MAX; n],
                    vec![0.0f64; n],
                    vec![0.0f64; n],
                    Vec::<usize>::new(),
                )
            },
            |(acc, dist, sigma, delta, order), k| {
                let s = sources[k];
                order.clear();
                dist[s] = 0;
                sigma[s] = 1.0;
                let mut queue = VecDeque::from([s]);
                while let Some(v) = queue.pop_front() {
                    order.push(v);
                    for &(w, e) in &self.rows[v] {
                        if !self.alive[e] {
                            continue;
                        }
                        let w = w as usize;
                        if dist[w] == usize::MAX {
                            dist[w] = dist[v] + 1;
                            queue.push_back(w);
                        }
                        if dist[w] == dist[v] + 1 {
                            sigma[w] += sigma[v];
                        }
                    }
                }
                for &w in order.iter().rev() {
                    for &(v, e) in &self.rows[w] {
                        let v = v as usize;
                        if self.alive[e] && dist[v] != usize::MAX && dist[v] + 1 == dist[w] {
                            let c = sigma[v] / sigma[w] * (1.0 + delta[w]);
                            acc[e] += c;
                            delta[v] += c;
                        }
                    }
                }
                for &v in order.iter() {
                    dist[v] = usize::MAX;
                    sigma[v] = 0.0;
                    delta[v] = 0.0;
                }
            },
        );
        for p in parts {
            for (e, x) in p.0.into_iter().enumerate() {
                if x != 0.0 {
                    eb[e] += x / 2.0;
                }
            }
        }
    }
}

fn partition_of(n: usize, edges: &[(u32, u32)], alive: &[bool]) -> Result<Vec<usize>> {
    let (s, t): (Vec<i64>, Vec<i64>) = edges
        .iter()
        .zip(alive)
        .filter(|(_, &a)| a)
        .map(|(&(a, b), _)| (a as i64, b as i64))
        .unzip();
    let h = Graph::from_edges(n, &s, &t, None, false)?;
    Ok(components(&h, ComponentKind::Weak)
        .membership
        .into_iter()
        .map(|c| c as usize)
        .collect())
}

/// Girvan-Newman clustering of an undirected network. See the module
/// documentation.
pub fn girvan_newman(g: &Graph, threads: usize) -> Result<GirvanNewman> {
    require_undirected(g, "Girvan-Newman clustering")?;
    let n = g.node_count();
    let simple = simple_rows(g);
    let mut edges: Vec<(u32, u32)> = Vec::new();
    for (v, row) in simple.iter().enumerate() {
        for &w in row {
            if (w as usize) > v {
                edges.push((v as u32, w));
            }
        }
    }
    let m = edges.len();
    let mut rows: Vec<Vec<(u32, usize)>> = vec![Vec::new(); n];
    for (e, &(a, b)) in edges.iter().enumerate() {
        rows[a as usize].push((b, e));
        rows[b as usize].push((a, e));
    }
    let (s, t): (Vec<i64>, Vec<i64>) = edges.iter().map(|&(a, b)| (a as i64, b as i64)).unzip();
    let simple_graph = Graph::from_edges(n, &s, &t, None, false)?;
    let mut rem = Remaining {
        rows,
        edges: edges.clone(),
        alive: vec![true; m],
    };

    let q_of = |memb: &[usize]| modularity(&simple_graph, &as_labels(memb), 1.0, false);
    let mut current = partition_of(n, &edges, &rem.alive)?;
    let mut best_q = q_of(&current)?;
    let mut best = current.clone();
    let mut path = vec![best_q];
    let mut removed = Vec::with_capacity(m);

    let mut eb = vec![0.0f64; m];
    let all: Vec<usize> = (0..n).collect();
    rem.betweenness(&all, &mut eb, threads);
    for _ in 0..m {
        let max = (0..m)
            .filter(|&e| rem.alive[e])
            .map(|e| eb[e])
            .fold(f64::NEG_INFINITY, f64::max);
        let cut = (0..m)
            .find(|&e| rem.alive[e] && eb[e] >= max - 1e-9 * max.abs())
            .expect("an edge remains");
        rem.alive[cut] = false;
        removed.push(edges[cut]);
        // Recompute betweenness within the component(s) the edge joined.
        let (a, b) = edges[cut];
        let mut affected = rem.component(a as usize);
        let split = !affected.contains(&(b as usize));
        if split {
            affected.extend(rem.component(b as usize));
        }
        for &v in &affected {
            for &(_, e) in &rem.rows[v] {
                eb[e] = 0.0;
            }
        }
        rem.betweenness(&affected, &mut eb, threads);
        if split {
            current = partition_of(n, &edges, &rem.alive)?;
        }
        let q = q_of(&current)?;
        path.push(q);
        if q > best_q || (best_q.is_nan() && !q.is_nan()) {
            best_q = q;
            best.clone_from(&current);
        }
    }
    Ok(GirvanNewman {
        membership: relabel(&best),
        modularity: best_q,
        removed,
        path,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn barbell_splits_at_the_bridge() {
        let g = und(
            6,
            &[
                (0, 1),
                (0, 2),
                (1, 2),
                (3, 4),
                (3, 5),
                (4, 5),
                (2, 3),
                (2, 3),
                (1, 1),
            ],
        );
        let r = girvan_newman(&g, 2).unwrap();
        assert_eq!(r.removed[0], (2, 3));
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        assert!((r.modularity - 5.0 / 14.0).abs() < 1e-15);
        assert_eq!(r.path.len(), 8);
        assert!(r.path[0].abs() < 1e-15);
    }

    #[test]
    fn edge_cases_and_direction() {
        let empty = und(3, &[]);
        let r = girvan_newman(&empty, 1).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.modularity.is_nan());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(girvan_newman(&d, 1).is_err());
    }

    #[test]
    fn thread_count_does_not_matter() {
        let mut e = Vec::new();
        for c in 0..4i64 {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    e.push((4 * c + i, 4 * c + j));
                }
            }
            e.push((4 * c, (4 * c + 5) % 16));
        }
        let g = und(16, &e);
        let a = girvan_newman(&g, 1).unwrap();
        let b = girvan_newman(&g, 4).unwrap();
        assert_eq!(a, b);
        assert_eq!(
            a.membership,
            vec![0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
        );
    }
}
