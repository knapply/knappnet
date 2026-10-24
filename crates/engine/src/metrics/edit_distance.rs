//! Exact graph edit distance.
//!
//! The graph edit distance between G1 and G2 is the least total cost of a
//! sequence of node and edge substitutions, deletions, and insertions that
//! turns G1 into a graph isomorphic to G2 (Sanfeliu and Fu 1983,
//! `sanfeliu1983distance`). Every edit path corresponds to a mapping that
//! sends each node of G1 to a distinct node of G2 or deletes it, inserts
//! the nodes of G2 left over, and then edits edges as the node mapping
//! dictates: an edge (a, b) of G1 whose endpoints map to an edge of G2 is
//! substituted, any other edge of G1 is deleted, and the edges of G2 not
//! covered are inserted. The distance is the minimum over mappings, found
//! here by depth-first branch and bound over the nodes of G1 (largest
//! degree first), pruned with a lower bound on the cost of the nodes and
//! edges still to be decided, in the manner of the exact depth-first
//! algorithm of Abu-Aisheh et al. (2015, `abuaisheh2015exact`) that
//! NetworkX's `graph_edit_distance` implements (`hagberg2008exploring`).
//! The result is the exact minimum, so it equals NetworkX's (which is
//! exact too when it runs to completion).
//!
//! Costs follow NetworkX's defaults and its `node_match` / `edge_match`
//! convention: deleting or inserting a node or an edge costs 1 (or the
//! given costs); substituting costs 0, or, when labels are given, 0 for
//! equal labels and the substitution cost (1 by default) for different
//! ones, the cost NetworkX assigns when a match function returns false.
//!
//! Conventions:
//! - Both networks must be directed or both undirected; a directed edge
//!   maps only onto an edge in the same direction.
//! - Self-loops are edges like any other: a loop at a maps onto a loop at
//!   the image of a.
//! - Parallel edges are merged (NetworkX's edit distance is defined on
//!   simple graphs); with edge labels they are an error.
//! - Weights are ignored; labels are integers the caller assigns (equal
//!   labels match).
//! - `roots` = (u, v) fixes u ↦ v, as NetworkX's `roots`.
//! - `upper_bound`: no answer (None) when the distance exceeds it.
//! - The search is exponential in the number of nodes; it stops with an
//!   error after `max_steps` mappings rather than return an approximation.
//!   Networks of a dozen nodes are routine.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;

/// Costs of edit operations. Substitution costs apply only between
/// differently labelled elements.
#[derive(Debug, Clone, Copy)]
pub struct EditCosts {
    /// Substituting a node by a node with a different label.
    pub node_subst: f64,
    /// Deleting a node.
    pub node_del: f64,
    /// Inserting a node.
    pub node_ins: f64,
    /// Substituting an edge by an edge with a different label.
    pub edge_subst: f64,
    /// Deleting an edge.
    pub edge_del: f64,
    /// Inserting an edge.
    pub edge_ins: f64,
}

impl Default for EditCosts {
    fn default() -> Self {
        EditCosts {
            node_subst: 1.0,
            node_del: 1.0,
            node_ins: 1.0,
            edge_subst: 1.0,
            edge_del: 1.0,
            edge_ins: 1.0,
        }
    }
}

/// Labels and options for [`graph_edit_distance`].
#[derive(Debug, Clone, Copy, Default)]
pub struct EditOptions<'a> {
    /// Node labels of G1 and G2 (one per node), or none (all nodes match).
    pub node_labels: Option<(&'a [i64], &'a [i64])>,
    /// Edge labels of G1 and G2 (one per edge), or none (all edges match).
    pub edge_labels: Option<(&'a [i64], &'a [i64])>,
    /// Costs.
    pub costs: EditCosts,
    /// A fixed pair (node of G1, node of G2).
    pub roots: Option<(usize, usize)>,
    /// Largest distance of interest.
    pub upper_bound: Option<f64>,
    /// Mappings tried before giving up (0 means 10^8).
    pub max_steps: u64,
}

/// A network as a dense edge table: `edge[a * n + b]` is the label of the
/// edge a -> b (both directions for undirected networks; 0 without
/// labels), or `None` if there is none.
struct Dense {
    n: usize,
    edge: Vec<Option<i64>>,
    label: Vec<i64>,
    degree: Vec<usize>,
    edges: usize,
}

impl Dense {
    fn new(g: &Graph, nodes: Option<&[i64]>, edges: Option<&[i64]>, what: &str) -> Result<Dense> {
        let n = g.node_count();
        if n > 4096 {
            return Err(Error::InvalidArgument(format!(
                "{what} has {n} nodes; exact edit distance is limited to small networks"
            )));
        }
        if let Some(l) = nodes {
            if l.len() != n {
                return Err(Error::LengthMismatch {
                    what: "node labels",
                    expected: n,
                    found: l.len(),
                });
            }
        }
        if let Some(l) = edges {
            if l.len() != g.edge_count() {
                return Err(Error::LengthMismatch {
                    what: "edge labels",
                    expected: g.edge_count(),
                    found: l.len(),
                });
            }
        }
        let mut edge = vec![None; n * n];
        let mut count = 0;
        for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
            let (s, t) = (s as usize, t as usize);
            let code = Some(edges.map_or(0, |l| l[e]));
            if edge[s * n + t].is_some() {
                if edges.is_some() {
                    return Err(Error::InvalidArgument(format!(
                        "{what} has parallel edges, which edge labels cannot tell apart"
                    )));
                }
                continue;
            }
            edge[s * n + t] = code;
            if !g.is_directed() {
                edge[t * n + s] = code;
            }
            count += 1;
        }
        let degree = (0..n)
            .map(|a| {
                (0..n)
                    .filter(|&b| edge[a * n + b].is_some() || edge[b * n + a].is_some())
                    .count()
            })
            .collect();
        Ok(Dense {
            n,
            edge,
            label: nodes.map_or_else(|| vec![0; n], |l| l.to_vec()),
            degree,
            edges: count,
        })
    }

    #[inline]
    fn e(&self, a: usize, b: usize) -> Option<i64> {
        self.edge[a * self.n + b]
    }

    #[inline]
    fn has(&self, a: usize, b: usize) -> bool {
        self.edge[a * self.n + b].is_some()
    }
}

struct Search<'a> {
    g1: &'a Dense,
    g2: &'a Dense,
    directed: bool,
    costs: EditCosts,
    order: Vec<usize>,
    /// Image of each G1 node in `order` position: Some(G2 node) or None
    /// (deleted).
    image: Vec<Option<usize>>,
    used: Vec<bool>,
    best: f64,
    steps: u64,
    max_steps: u64,
    /// Edges of G1 with both endpoints decided, after each depth.
    decided1: Vec<usize>,
}

impl Search<'_> {
    /// Cost of the edge pair (a, b) of G1 against (x, y) of G2, x or y
    /// None for a deleted node or a node not yet placed.
    fn edge_cost(&self, a: usize, b: usize, x: Option<usize>, y: Option<usize>) -> f64 {
        let e1 = self.g1.e(a, b);
        let e2 = match (x, y) {
            (Some(x), Some(y)) => self.g2.e(x, y),
            _ => None,
        };
        match (e1, e2) {
            (Some(p), Some(q)) if p == q => 0.0,
            (Some(_), Some(_)) => self.costs.edge_subst,
            (Some(_), None) => self.costs.edge_del,
            (None, Some(_)) => self.costs.edge_ins,
            (None, None) => 0.0,
        }
    }

    /// Cost added by deciding the node at depth `k` as `x`.
    fn step_cost(&self, k: usize, x: Option<usize>) -> f64 {
        let a = self.order[k];
        let mut c = match x {
            Some(x) if self.g1.label[a] == self.g2.label[x] => 0.0,
            Some(_) => self.costs.node_subst,
            None => self.costs.node_del,
        };
        c += self.edge_cost(a, a, x, x);
        for j in 0..k {
            let b = self.order[j];
            let y = self.image[j];
            c += self.edge_cost(a, b, x, y);
            if self.directed {
                c += self.edge_cost(b, a, y, x);
            }
        }
        c
    }

    /// Cost of inserting the unused nodes of G2 and every edge touching
    /// them.
    fn finish_cost(&self) -> f64 {
        let n2 = self.g2.n;
        let mut c = 0.0;
        for x in 0..n2 {
            if self.used[x] {
                continue;
            }
            c += self.costs.node_ins;
            for y in 0..n2 {
                // Count each edge once: from its unused endpoint, or from
                // the smaller one when both are unused.
                let counts = self.used[y] || y >= x;
                if !counts {
                    continue;
                }
                if self.g2.has(x, y) {
                    c += self.costs.edge_ins;
                }
                if self.directed && y != x && self.g2.has(y, x) {
                    c += self.costs.edge_ins;
                }
            }
        }
        c
    }

    /// A lower bound on the cost still to come at depth `k`.
    fn bound(&self, k: usize, placed: usize) -> f64 {
        let r1 = self.order.len() - k;
        let r2 = self.g2.n - placed;
        let nodes = if r1 > r2 {
            (r1 - r2) as f64 * self.costs.node_del
        } else {
            (r2 - r1) as f64 * self.costs.node_ins
        };
        // Edges not yet decided on either side.
        let e1 = self.g1.edges - self.decided1[k];
        let mut d2 = 0;
        let placed_nodes: Vec<usize> = self.image[..k].iter().flatten().copied().collect();
        for (i, &x) in placed_nodes.iter().enumerate() {
            for &y in &placed_nodes[..=i] {
                if self.g2.has(x, y) {
                    d2 += 1;
                }
                if self.directed && x != y && self.g2.has(y, x) {
                    d2 += 1;
                }
            }
        }
        let e2 = self.g2.edges - d2;
        let edges = if e1 > e2 {
            (e1 - e2) as f64 * self.costs.edge_del
        } else {
            (e2 - e1) as f64 * self.costs.edge_ins
        };
        nodes + edges
    }

    fn run(&mut self, k: usize, placed: usize, cost: f64, fixed: Option<usize>) -> Result<()> {
        self.steps += 1;
        if self.steps > self.max_steps {
            return Err(Error::InvalidArgument(format!(
                "the edit distance search tried more than {} mappings; the networks are too \
                 large for an exact answer",
                self.max_steps
            )));
        }
        if k == self.order.len() {
            let total = cost + self.finish_cost();
            if total < self.best {
                self.best = total;
            }
            return Ok(());
        }
        if cost + self.bound(k, placed) >= self.best {
            return Ok(());
        }
        let a = self.order[k];
        let mut options: Vec<(f64, Option<usize>)> = Vec::new();
        if let (0, Some(x)) = (k, fixed) {
            options.push((self.step_cost(k, Some(x)), Some(x)));
        } else {
            for x in 0..self.g2.n {
                if !self.used[x] {
                    options.push((self.step_cost(k, Some(x)), Some(x)));
                }
            }
            options.push((self.step_cost(k, None), None));
        }
        // Cheapest first, then nodes of similar degree.
        options.sort_by(|p, q| {
            p.0.total_cmp(&q.0).then_with(|| {
                let d = |o: Option<usize>| {
                    o.map_or(usize::MAX, |x| {
                        self.g1.degree[a].abs_diff(self.g2.degree[x])
                    })
                };
                d(p.1).cmp(&d(q.1))
            })
        });
        for (c, x) in options {
            if cost + c >= self.best {
                continue;
            }
            self.image[k] = x;
            if let Some(x) = x {
                self.used[x] = true;
            }
            self.run(k + 1, placed + x.is_some() as usize, cost + c, fixed)?;
            if let Some(x) = x {
                self.used[x] = false;
            }
            self.image[k] = None;
        }
        Ok(())
    }
}

/// The exact graph edit distance between `g1` and `g2`, or `None` when it
/// exceeds `opts.upper_bound`. See the module documentation.
pub fn graph_edit_distance(g1: &Graph, g2: &Graph, opts: EditOptions) -> Result<Option<f64>> {
    if g1.is_directed() != g2.is_directed() {
        return Err(Error::InvalidArgument(
            "both networks must be directed or both undirected".into(),
        ));
    }
    let c = opts.costs;
    for x in [
        c.node_subst,
        c.node_del,
        c.node_ins,
        c.edge_subst,
        c.edge_del,
        c.edge_ins,
    ] {
        if !(x >= 0.0 && x.is_finite()) {
            return Err(Error::InvalidArgument(
                "edit costs must be finite and non-negative".into(),
            ));
        }
    }
    let d1 = Dense::new(
        g1,
        opts.node_labels.map(|p| p.0),
        opts.edge_labels.map(|p| p.0),
        "G1",
    )?;
    let d2 = Dense::new(
        g2,
        opts.node_labels.map(|p| p.1),
        opts.edge_labels.map(|p| p.1),
        "G2",
    )?;
    let mut order: Vec<usize> = (0..d1.n).collect();
    order.sort_by(|&a, &b| d1.degree[b].cmp(&d1.degree[a]).then(a.cmp(&b)));
    let mut fixed = None;
    if let Some((u, v)) = opts.roots {
        if u >= d1.n || v >= d2.n {
            return Err(Error::InvalidArgument(
                "roots must name a node of each network".into(),
            ));
        }
        let pos = order.iter().position(|&a| a == u).expect("u is a node");
        order.remove(pos);
        order.insert(0, u);
        fixed = Some(v);
    }
    // Edges of G1 decided once the first k nodes of the order are.
    let mut pos = vec![0usize; d1.n];
    for (k, &a) in order.iter().enumerate() {
        pos[a] = k;
    }
    let mut decided1 = vec![0usize; d1.n + 1];
    let mut seen: HashMap<(usize, usize), ()> = HashMap::new();
    for a in 0..d1.n {
        for b in 0..d1.n {
            if d1.has(a, b) {
                let key = if g1.is_directed() {
                    (a, b)
                } else {
                    (a.min(b), a.max(b))
                };
                if seen.insert(key, ()).is_none() {
                    decided1[pos[a].max(pos[b]) + 1] += 1;
                }
            }
        }
    }
    for k in 0..d1.n {
        decided1[k + 1] += decided1[k];
    }
    let ub = opts.upper_bound.unwrap_or(f64::INFINITY);
    let mut s = Search {
        g1: &d1,
        g2: &d2,
        directed: g1.is_directed(),
        costs: c,
        order,
        image: vec![None; d1.n],
        used: vec![false; d2.n],
        // Strictly above the bound, so a path costing exactly the bound
        // is found.
        best: if ub.is_finite() {
            ub + ub.abs() * 1e-12 + 1e-12
        } else {
            ub
        },
        steps: 0,
        max_steps: if opts.max_steps == 0 {
            100_000_000
        } else {
            opts.max_steps
        },
        decided1,
    };
    s.run(0, 0, 0.0, fixed)?;
    Ok((s.best <= ub).then_some(s.best))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn small_distances() {
        let tri = graph(3, &[(0, 1), (1, 2), (2, 0)], false);
        let path = graph(3, &[(0, 1), (1, 2)], false);
        let o = EditOptions::default();
        assert_eq!(graph_edit_distance(&tri, &path, o).unwrap(), Some(1.0));
        assert_eq!(graph_edit_distance(&tri, &tri, o).unwrap(), Some(0.0));
        // A path on 3 nodes to one on 4: insert a node and an edge.
        let p4 = graph(4, &[(0, 1), (1, 2), (2, 3)], false);
        assert_eq!(graph_edit_distance(&path, &p4, o).unwrap(), Some(2.0));
        // Reversing an arc: delete it and insert the reverse.
        let a = graph(2, &[(0, 1)], true);
        let b = graph(2, &[(1, 0)], true);
        assert_eq!(graph_edit_distance(&a, &b, o).unwrap(), Some(0.0));
        let c = graph(2, &[(0, 1), (0, 0)], true);
        assert_eq!(graph_edit_distance(&a, &c, o).unwrap(), Some(1.0));
        let ub = EditOptions {
            upper_bound: Some(1.5),
            ..o
        };
        assert_eq!(graph_edit_distance(&path, &p4, ub).unwrap(), None);
    }
}
