//! Edge connectivity between pairs of nodes, and lambda sets.
//!
//! The edge connectivity lambda(u, v) of two nodes is the smallest number
//! of edges whose removal disconnects them, which equals the largest
//! number of edge-disjoint paths between them (Menger's theorem; Ford and
//! Fulkerson 1956, `ford1956maximal`). A lambda set is a set of nodes S
//! such that any two members are more edge-connected to each other than
//! either is to any node outside: lambda(a, b) > lambda(a, c) for all a, b
//! in S and c not in S (Borgatti, Everett, and Shirey 1990,
//! `borgatti1990ls`). Because lambda(a, c) >= min(lambda(a, b),
//! lambda(b, c)), the relation "lambda(a, b) >= k" is an equivalence for
//! every k >= 1, and its classes with at least two members are exactly the
//! lambda sets; increasing k gives the nested hierarchy that UCINET's
//! lambda set routine reports.
//!
//! All pairwise connectivities come from n - 1 maximum flow computations
//! with Gusfield's (1990, `gusfield1990very`) method, which builds a tree
//! on the nodes (equivalent to the Gomory-Hu tree for flow values; Gomory
//! and Hu 1961, `gomory1961multi`) such that lambda(u, v) is the smallest
//! edge value on the tree path between u and v. Each flow is found by
//! breadth-first augmenting paths (Edmonds and Karp 1972,
//! `edmonds1972theoretical`), with every edge of capacity 1 in both
//! directions.
//!
//! Conventions:
//! - The undirected simple reading of the network is used (see
//!   [`crate::metrics::factions`]): directions are ignored, parallel edges
//!   count once, loops are ignored, and weights are not used. This is
//!   NetworkX's `edge_connectivity(G, u, v)` on the simple graph.
//! - Nodes in different components have connectivity 0, so the classes
//!   for k >= 1 never join components. The diagonal of the connectivity
//!   matrix is 0.

use std::collections::VecDeque;

use super::simple_rows;
use crate::graph::Graph;
use crate::matrix::Dense;

/// A tree on the nodes in which lambda(u, v) is the smallest `value` on
/// the path from u to v. Node 0 is the root (`parent[0] == 0`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlowTree {
    /// Parent of each node.
    pub parent: Vec<usize>,
    /// Connectivity between each node and its parent (0 for the root).
    pub value: Vec<u64>,
}

/// Residual network with unit capacities in both directions of each edge.
struct Flow {
    head: Vec<usize>,
    to: Vec<u32>,
    cap: Vec<u8>,
}

impl Flow {
    fn new(rows: &[Vec<u32>]) -> Flow {
        let n = rows.len();
        let mut head = vec![0usize; n + 1];
        for (v, r) in rows.iter().enumerate() {
            head[v + 1] = head[v] + r.len();
        }
        let to: Vec<u32> = rows.iter().flatten().copied().collect();
        Flow {
            head,
            to,
            cap: Vec::new(),
        }
    }

    /// Position of the arc w -> v.
    fn reverse(&self, v: usize, w: usize) -> usize {
        let row = &self.to[self.head[w]..self.head[w + 1]];
        self.head[w] + row.binary_search(&(v as u32)).expect("symmetric rows")
    }

    /// Maximum flow from s to t; returns its value and which nodes are on
    /// the source side of a minimum cut.
    fn max_flow(&mut self, s: usize, t: usize) -> (u64, Vec<bool>) {
        let n = self.head.len() - 1;
        // Each undirected edge gives capacity 1 each way; pushing one unit
        // along v -> w frees one unit along w -> v.
        self.cap.clear();
        self.cap.resize(self.to.len(), 1);
        let mut flow = 0u64;
        let mut pred = vec![usize::MAX; n];
        loop {
            for p in pred.iter_mut() {
                *p = usize::MAX;
            }
            pred[s] = s;
            let mut queue = VecDeque::from([s]);
            let mut reached = false;
            while let Some(v) = queue.pop_front() {
                for a in self.head[v]..self.head[v + 1] {
                    let w = self.to[a] as usize;
                    if self.cap[a] > 0 && pred[w] == usize::MAX {
                        pred[w] = a;
                        if w == t {
                            reached = true;
                            break;
                        }
                        queue.push_back(w);
                    }
                }
                if reached {
                    break;
                }
            }
            if !reached {
                let side = pred.iter().map(|&p| p != usize::MAX).collect();
                return (flow, side);
            }
            let mut w = t;
            while w != s {
                let a = pred[w];
                // The arc a runs from its row's node to w.
                let v = self.head.partition_point(|&h| h <= a) - 1;
                self.cap[a] -= 1;
                let back = self.reverse(v, w);
                self.cap[back] += 1;
                w = v;
            }
            flow += 1;
        }
    }
}

/// Gusfield's equivalent flow tree for edge connectivity. See the module
/// documentation.
pub fn edge_connectivity_tree(g: &Graph) -> FlowTree {
    let rows = simple_rows(g);
    let n = rows.len();
    let mut parent = vec![0usize; n];
    let mut value = vec![0u64; n];
    let mut flow = Flow::new(&rows);
    for s in 1..n {
        let t = parent[s];
        let (f, side) = flow.max_flow(s, t);
        value[s] = f;
        for i in (s + 1)..n {
            if side[i] && parent[i] == t {
                parent[i] = s;
            }
        }
    }
    FlowTree { parent, value }
}

/// Edge connectivity between every pair of nodes, as an n-by-n matrix with
/// zeros on the diagonal. See the module documentation.
pub fn lambda_connectivity(g: &Graph) -> Dense {
    let tree = edge_connectivity_tree(g);
    let n = tree.parent.len();
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    for v in 1..n {
        children[tree.parent[v]].push(v);
    }
    let mut out = Dense::zeros(n, n);
    // From each node, walk the tree and carry the path minimum.
    let mut stack: Vec<(usize, usize, u64)> = Vec::new();
    for root in 0..n {
        stack.push((root, usize::MAX, u64::MAX));
        while let Some((v, from, min)) = stack.pop() {
            if v != root {
                out.set(root, v, min as f64);
            }
            let up = tree.parent[v];
            if v != 0 && up != from {
                stack.push((up, v, min.min(tree.value[v])));
            }
            for &c in &children[v] {
                if c != from {
                    stack.push((c, v, min.min(tree.value[c])));
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn small_examples() {
        // Two triangles joined by one edge, plus an isolate.
        let g = und(7, &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)]);
        let l = lambda_connectivity(&g);
        assert_eq!(l.get(0, 1), 2.0);
        assert_eq!(l.get(0, 4), 1.0);
        assert_eq!(l.get(3, 5), 2.0);
        assert_eq!(l.get(0, 6), 0.0);
        assert_eq!(l.get(4, 4), 0.0);
        let k4 = und(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (0, 1)]);
        let l = lambda_connectivity(&k4);
        assert!((0..4).all(|a| (0..4).all(|b| a == b || l.get(a, b) == 3.0)));
    }

    #[test]
    fn tree_matches_pairwise_flows() {
        for seed in 0..20 {
            let mut rng = Stream::new(seed, 0);
            let n = 12;
            let mut e = Vec::new();
            for a in 0..n as i64 {
                for b in (a + 1)..n as i64 {
                    if rng.uniform() < 0.3 {
                        e.push((a, b));
                    }
                }
            }
            let g = und(n, &e);
            let l = lambda_connectivity(&g);
            let rows = simple_rows(&g);
            let mut flow = Flow::new(&rows);
            for a in 0..n {
                for b in 0..n {
                    if a != b {
                        assert_eq!(l.get(a, b), flow.max_flow(a, b).0 as f64);
                    }
                }
            }
        }
    }
}
