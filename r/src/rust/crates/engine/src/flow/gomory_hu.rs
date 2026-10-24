//! Gomory-Hu trees: all pairwise minimum cuts of an undirected network in
//! n - 1 maximum flows.
//!
//! A Gomory-Hu tree (Gomory and Hu 1961, `gomory1961multi`) is a weighted
//! tree on the nodes in which, for every pair u, v, the smallest weight on
//! the tree path between them is the capacity of a minimum u-v cut, and
//! removing that tree edge splits the nodes into the two sides of such a
//! cut. Gusfield's (1990, `gusfield1990very`) method builds one with n - 1
//! maximum flows in the original network, without contracting it: start
//! with a star on node 0; for each other node s in turn, with t its
//! current neighbour in the tree, find a minimum s-t cut (X, V - X), move
//! to s every node of X whose tree neighbour is t, and if t's own tree
//! neighbour lies in X, swap the positions of s and t.
//!
//! NetworkX's `gomory_hu_tree` (`hagberg2008exploring`) runs exactly these
//! steps over its nodes in order, taking for X the source side of the cut
//! its `minimum_cut` returns, which is canonical (the largest source side,
//! [`crate::flow::CutSide::Largest`]); this implementation does the same, so with node
//! order 0..n it returns the same tree, edge for edge.
//!
//! Conventions: undirected networks only (NetworkX refuses directed
//! graphs); edge capacities as in [`crate::flow::maxflow`] (default 1,
//! possibly infinite, parallel edges adding up, self-loops ignored).
//! Nodes in different components are joined by tree edges of weight 0.
//! The steps are sequential, so the thread count plays no part.

use std::collections::HashMap;

use super::maxflow::edge_network;
use crate::error::{Error, Result};
use crate::graph::Graph;

/// A Gomory-Hu tree rooted at node 0: every other node v is joined to
/// `parent[v]` by an edge of weight `weight[v]`.
#[derive(Debug, Clone, PartialEq)]
pub struct GomoryHuTree {
    /// Tree neighbour toward the root (`parent[0] == 0`).
    pub parent: Vec<usize>,
    /// Weight of the edge to the parent (0 for the root).
    pub weight: Vec<f64>,
}

impl GomoryHuTree {
    /// The minimum cut capacity between `u` and `v` (infinite for u == v).
    pub fn min_cut_value(&self, u: usize, v: usize) -> f64 {
        if u == v {
            return f64::INFINITY;
        }
        // Mark u's path to the root with the running minimum.
        let n = self.parent.len();
        let mut best = vec![f64::NAN; n];
        let mut x = u;
        let mut low = f64::INFINITY;
        loop {
            best[x] = low;
            if x == 0 {
                break;
            }
            low = low.min(self.weight[x]);
            x = self.parent[x];
        }
        let mut y = v;
        let mut low_v = f64::INFINITY;
        while best[y].is_nan() {
            low_v = low_v.min(self.weight[y]);
            y = self.parent[y];
        }
        low_v.min(best[y])
    }
}

/// Gomory-Hu tree of an undirected network. See the module documentation.
pub fn gomory_hu_tree(g: &Graph, capacity: Option<&[f64]>) -> Result<GomoryHuTree> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(
            "a Gomory-Hu tree is defined for undirected networks".into(),
        ));
    }
    let n = g.node_count();
    if n == 0 {
        return Err(Error::InvalidArgument(
            "a network with no nodes has no Gomory-Hu tree".into(),
        ));
    }
    let mut net = edge_network(g, capacity)?;
    let root = 0usize;
    let mut tree = vec![root; n];
    let mut labels: HashMap<(usize, usize), f64> = HashMap::new();
    for source in 1..n {
        let target = tree[source];
        net.reset();
        let cut = net.max_flow(source, target, None)?;
        let sink_side = net.reaching(&[target]);
        let in_x = |v: usize| !sink_side[v];
        labels.insert((source, target), cut);
        // Each node's test reads only its own entry, so testing all first
        // is the same as updating in turn.
        let movers: Vec<usize> = (1..n)
            .filter(|&v| v != source && in_x(v) && tree[v] == target)
            .collect();
        for node in movers {
            tree[node] = source;
            let l = labels.get(&(node, target)).copied().unwrap_or(cut);
            labels.insert((node, source), l);
        }
        if target != root && in_x(tree[target]) {
            let tt = tree[target];
            let l = labels.get(&(target, tt)).copied().unwrap_or(f64::NAN);
            labels.insert((source, tt), l);
            labels.insert((target, source), cut);
            tree[source] = tt;
            tree[target] = source;
        }
    }
    let mut weight = vec![0.0; n];
    for v in 1..n {
        weight[v] = labels[&(v, tree[v])];
    }
    tree[root] = root;
    Ok(GomoryHuTree {
        parent: tree,
        weight,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flow::maxflow::maximum_flow;
    use crate::rng::Stream;

    #[test]
    fn every_pair_matches_a_flow() {
        let mut rng = Stream::new(7, 0);
        for _ in 0..20 {
            let n = 8;
            let mut s = Vec::new();
            let mut t = Vec::new();
            let mut c = Vec::new();
            for a in 0..n {
                for b in (a + 1)..n {
                    if rng.uniform() < 0.4 {
                        s.push(a as i64);
                        t.push(b as i64);
                        c.push((1 + rng.below(5)) as f64);
                    }
                }
            }
            let g = Graph::from_edges(n, &s, &t, None, false).unwrap();
            let tree = gomory_hu_tree(&g, Some(&c)).unwrap();
            for u in 0..n {
                for v in (u + 1)..n {
                    let f = maximum_flow(&g, &[u], &[v], Some(&c)).unwrap().value;
                    assert_eq!(tree.min_cut_value(u, v), f, "pair {u} {v}");
                }
            }
        }
    }
}
