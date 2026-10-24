//! The k shortest simple paths between two nodes (Yen 1971,
//! `yen1971finding`), and shortest paths that avoid given nodes and edges.
//!
//! Yen's algorithm lists simple paths (no repeated node) in order of
//! length. Each accepted path is the parent of candidates that follow it up
//! to a spur node and then leave it by the shortest path that avoids the
//! nodes before the spur node and the arcs that earlier accepted paths with
//! the same beginning take out of the spur node. Following Lawler (1972,
//! `lawler1972procedure`), spur nodes start at the node where the parent
//! left its own parent, since earlier ones were already explored. Each
//! candidate needs one shortest-path search, so k paths cost O(k n) searches
//! in the worst case.
//!
//! [`shortest_path_avoiding`] is the search Yen's algorithm is built on,
//! offered on its own: a shortest path that avoids a set of nodes and a set
//! of edges, the building block of penalty and exclusion methods for
//! finding several dissimilar paths.
//!
//! Conventions: lengths must be non-negative; paths are node sequences, so
//! parallel edges count once (through the shortest), and an avoided edge
//! removes only that edge; self-loops are never used. Among candidates of
//! equal length the one with the lexicographically smallest node sequence
//! is taken first, which makes the list deterministic; NetworkX's
//! `shortest_simple_paths` takes equal candidates in the order it found
//! them, so paths of equal length can come in another order. Unweighted networks use breadth-first
//! search, so lengths are hop counts.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

use crate::error::Result;
use crate::paths::{NodeLists, PathGraph, NO_NODE};
use crate::weights::IndexedHeap;

/// Scratch space for searches with avoided nodes and arcs.
struct Masked {
    dist: Vec<f64>,
    parent: Vec<u32>,
    done: Vec<bool>,
    touched: Vec<u32>,
    heap: IndexedHeap,
    node_banned: Vec<bool>,
    arc_banned: Vec<bool>,
    queue: Vec<u32>,
}

impl Masked {
    fn new(pg: &PathGraph) -> Masked {
        let n = pg.node_count();
        Masked {
            dist: vec![f64::INFINITY; n],
            parent: vec![NO_NODE; n],
            done: vec![false; n],
            touched: Vec::new(),
            heap: IndexedHeap::new(n),
            node_banned: vec![false; n],
            arc_banned: vec![false; pg.arc_count()],
            queue: Vec::new(),
        }
    }

    fn clear(&mut self) {
        for &v in &self.touched {
            let v = v as usize;
            self.dist[v] = f64::INFINITY;
            self.parent[v] = NO_NODE;
            self.done[v] = false;
        }
        self.touched.clear();
        self.heap.clear();
        self.queue.clear();
    }

    /// Shortest path from s to t avoiding banned nodes and arcs, as a node
    /// sequence with its length.
    fn search(&mut self, pg: &PathGraph, s: usize, t: usize) -> Option<(f64, Vec<u32>)> {
        self.clear();
        if self.node_banned[s] || self.node_banned[t] {
            return None;
        }
        self.dist[s] = 0.0;
        self.touched.push(s as u32);
        let found = if pg.is_unit() {
            self.queue.push(s as u32);
            self.done[s] = true;
            let mut head = 0;
            let mut found = s == t;
            'bfs: while head < self.queue.len() && !found {
                let v = self.queue[head] as usize;
                head += 1;
                for a in pg.out_arcs(v) {
                    let w = pg.head(a) as usize;
                    if self.done[w] || self.arc_banned[a] || self.node_banned[w] {
                        continue;
                    }
                    self.done[w] = true;
                    self.dist[w] = self.dist[v] + 1.0;
                    self.parent[w] = v as u32;
                    self.touched.push(w as u32);
                    self.queue.push(w as u32);
                    if w == t {
                        found = true;
                        break 'bfs;
                    }
                }
            }
            found
        } else {
            self.heap.push_or_decrease(s as u32, 0.0);
            let mut found = false;
            while let Some((d, v)) = self.heap.pop() {
                let v = v as usize;
                self.done[v] = true;
                if v == t {
                    found = true;
                    break;
                }
                for a in pg.out_arcs(v) {
                    let w = pg.head(a) as usize;
                    if self.done[w] || self.arc_banned[a] || self.node_banned[w] {
                        continue;
                    }
                    let alt = d + pg.len(a);
                    if alt < self.dist[w] {
                        if self.dist[w].is_infinite() {
                            self.touched.push(w as u32);
                        }
                        self.dist[w] = alt;
                        self.parent[w] = v as u32;
                        self.heap.push_or_decrease(w as u32, alt);
                    }
                }
            }
            found
        };
        if !found {
            return None;
        }
        let mut path = vec![t as u32];
        let mut cur = t;
        while self.parent[cur] != NO_NODE {
            cur = self.parent[cur] as usize;
            path.push(cur as u32);
        }
        path.reverse();
        Some((self.dist[t], path))
    }
}

/// Length and nodes of a shortest path from `s` to `t` that uses none of
/// `avoid_nodes` (which must not include `s` or `t`, or there is no path)
/// and none of `avoid_edges` (edge indices in input order), or `None`.
pub fn shortest_path_avoiding(
    pg: &PathGraph,
    s: u32,
    t: u32,
    avoid_nodes: &[u32],
    avoid_edges: &[u32],
) -> Result<Option<(f64, Vec<u32>)>> {
    pg.check_node(s as usize, "source")?;
    pg.check_node(t as usize, "target")?;
    for &v in avoid_nodes {
        pg.check_node(v as usize, "avoided")?;
    }
    for &e in avoid_edges {
        if e as usize >= pg.edge_count() {
            return Err(crate::Error::InvalidArgument(format!(
                "edge index {e} is out of range for a network with {} edges",
                pg.edge_count()
            )));
        }
    }
    pg.require_nonnegative("this search")?;
    let mut ws = Masked::new(pg);
    for &v in avoid_nodes {
        ws.node_banned[v as usize] = true;
    }
    if !avoid_edges.is_empty() {
        let banned: HashSet<u32> = avoid_edges.iter().copied().collect();
        for a in 0..pg.arc_count() {
            ws.arc_banned[a] = banned.contains(&pg.edge(a));
        }
    }
    Ok(ws.search(pg, s as usize, t as usize))
}

/// A candidate path, ordered by length and then by node sequence.
#[derive(Debug, Clone, PartialEq)]
struct Candidate {
    len: f64,
    path: Vec<u32>,
    /// Index of the spur node at which it left its parent.
    deviation: usize,
}

impl Eq for Candidate {}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed so that BinaryHeap pops the shortest first.
        other
            .len
            .total_cmp(&self.len)
            .then_with(|| other.path.cmp(&self.path))
    }
}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// The `k` shortest simple paths from `s` to `t` in order of length (fewer
/// if there are fewer), with their lengths. See the module documentation.
pub fn k_shortest_simple_paths(pg: &PathGraph, s: u32, t: u32, k: usize) -> Result<NodeLists> {
    pg.check_node(s as usize, "source")?;
    pg.check_node(t as usize, "target")?;
    pg.require_nonnegative("Yen's algorithm")?;
    let mut out = NodeLists::new();
    if k == 0 {
        return Ok(out);
    }
    let mut ws = Masked::new(pg);
    let Some((len, path)) = ws.search(pg, s as usize, t as usize) else {
        return Ok(out);
    };
    let len = pg.path_length(&path).unwrap_or(len);
    let mut accepted: Vec<Candidate> = vec![Candidate {
        len,
        path,
        deviation: 0,
    }];
    let mut seen: HashSet<Vec<u32>> = HashSet::new();
    seen.insert(accepted[0].path.clone());
    let mut candidates: BinaryHeap<Candidate> = BinaryHeap::new();
    let mut banned_arcs: Vec<usize> = Vec::new();
    while accepted.len() < k {
        let last = accepted.last().expect("one path").clone();
        for i in last.deviation..last.path.len().saturating_sub(1) {
            let root = &last.path[..=i];
            let spur = root[i] as usize;
            // Arcs out of the spur node taken by accepted paths with this root.
            for p in &accepted {
                if p.path.len() > i + 1 && &p.path[..=i] == root {
                    let next = p.path[i + 1] as usize;
                    for a in pg.out_arcs(spur) {
                        if pg.head(a) as usize == next && !ws.arc_banned[a] {
                            ws.arc_banned[a] = true;
                            banned_arcs.push(a);
                        }
                    }
                }
            }
            for &v in &root[..i] {
                ws.node_banned[v as usize] = true;
            }
            if let Some((_, spur_path)) = ws.search(pg, spur, t as usize) {
                let mut path = root[..i].to_vec();
                path.extend_from_slice(&spur_path);
                if !seen.contains(&path) {
                    let len = pg.path_length(&path).expect("a path of arcs");
                    seen.insert(path.clone());
                    candidates.push(Candidate {
                        len,
                        path,
                        deviation: i,
                    });
                }
            }
            for &v in &root[..i] {
                ws.node_banned[v as usize] = false;
            }
            for a in banned_arcs.drain(..) {
                ws.arc_banned[a] = false;
            }
        }
        match candidates.pop() {
            Some(c) => accepted.push(c),
            None => break,
        }
    }
    for c in &accepted {
        out.push_with(&c.path, c.len);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::paths::{ArcLengths, PathGraphOptions};
    use crate::rng::Stream;
    use crate::weights::WeightKind;

    /// Every simple path from s to t by exhaustive search, sorted by
    /// (length, nodes).
    fn all_simple(pg: &PathGraph, s: usize, t: usize) -> Vec<(f64, Vec<u32>)> {
        let mut out = Vec::new();
        let mut path = vec![s as u32];
        fn go(pg: &PathGraph, t: usize, path: &mut Vec<u32>, out: &mut Vec<(f64, Vec<u32>)>) {
            let v = *path.last().unwrap() as usize;
            if v == t {
                out.push((pg.path_length(path).unwrap(), path.clone()));
                return;
            }
            let mut nexts: Vec<u32> = pg.out_arcs(v).map(|a| pg.head(a)).collect();
            nexts.dedup();
            for w in nexts {
                if !path.contains(&w) {
                    path.push(w);
                    go(pg, t, path, out);
                    path.pop();
                }
            }
        }
        go(pg, t, &mut path, &mut out);
        out.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        out
    }

    #[test]
    fn matches_exhaustive_enumeration() {
        for seed in 0..60u64 {
            let mut rng = Stream::new(seed, 5);
            let n = 3 + rng.below(5) as usize;
            let m = rng.below(3 * n as u64) as usize + 2;
            let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let w: Vec<f64> = (0..m).map(|_| (1 + rng.below(4)) as f64).collect();
            let directed = seed % 2 == 0;
            let g = Graph::from_edges(n, &s, &t, Some(&w), directed).unwrap();
            for weighted in [false, true] {
                let lengths = if weighted {
                    ArcLengths::Weights(WeightKind::Distance)
                } else {
                    ArcLengths::Unit
                };
                let pg = PathGraph::new(
                    &g,
                    PathGraphOptions {
                        lengths,
                        ..Default::default()
                    },
                )
                .unwrap();
                let (a, b) = (0usize, n - 1);
                let want = all_simple(&pg, a, b);
                let got = k_shortest_simple_paths(&pg, a as u32, b as u32, want.len() + 3).unwrap();
                assert_eq!(got.len(), want.len(), "seed {seed}");
                for (i, (len, _)) in want.iter().enumerate() {
                    // Lengths are non-decreasing and equal the true k-th
                    // lengths; equal-length paths come in lexicographic order.
                    assert_eq!(got.values[i], *len, "seed {seed} path {i}");
                    assert_eq!(pg.path_length(got.get(i)), Some(*len));
                }
                let mut sorted: Vec<Vec<u32>> =
                    (0..got.len()).map(|i| got.get(i).to_vec()).collect();
                sorted.sort();
                sorted.dedup();
                assert_eq!(sorted.len(), got.len());
            }
        }
    }

    #[test]
    fn avoiding_nodes_and_edges() {
        let g = Graph::from_edges(4, &[0, 1, 0, 2], &[1, 3, 2, 3], None, false).unwrap();
        let pg = PathGraph::new(&g, PathGraphOptions::default()).unwrap();
        let (_, p) = shortest_path_avoiding(&pg, 0, 3, &[1], &[])
            .unwrap()
            .unwrap();
        assert_eq!(p, vec![0, 2, 3]);
        let (_, p) = shortest_path_avoiding(&pg, 0, 3, &[], &[0])
            .unwrap()
            .unwrap();
        assert_eq!(p, vec![0, 2, 3]);
        assert!(shortest_path_avoiding(&pg, 0, 3, &[1, 2], &[])
            .unwrap()
            .is_none());
    }
}
