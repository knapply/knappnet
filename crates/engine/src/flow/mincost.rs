//! Minimum-cost flows by the primal network simplex method.
//!
//! Problem: in a directed network where node v has demand b(v) (positive
//! for a node that must receive flow, negative for one that supplies it;
//! the demands sum to zero) and edge e has capacity u(e) >= 0 and cost
//! c(e), find a flow 0 <= x(e) <= u(e) whose inflow minus outflow at every
//! node equals its demand and whose cost sum c(e) x(e) is least. The sign
//! of the demand is NetworkX's (`hagberg2008exploring`), so a demand
//! vector can be passed to either unchanged.
//!
//! Method: the network simplex method (Dantzig 1951, `dantzig1951application`)
//! as presented by Ahuja, Magnanti and Orlin (1993, `ahuja1993network`,
//! chapter 11): an artificial root joined to every node by an arc of
//! large cost M starts a feasible spanning tree; each pivot brings in a
//! non-tree arc of negative reduced cost, sends flow around the cycle it
//! closes, and drops a blocking arc. Trees are kept strongly feasible
//! (Cunningham 1976, `cunningham1976network`) by dropping the last
//! blocking arc met when the cycle is traversed in the direction of flow
//! from its apex, which prevents cycling under degeneracy. Entering arcs
//! are chosen by block search (Grigoriadis 1986, `grigoriadis1986efficient`):
//! arcs are scanned cyclically in blocks of about sqrt(m), and the most
//! violating arc of the first block that has one enters. M is one more
//! than n times the largest absolute cost, more than any simple path can
//! cost, so artificial arcs carry flow at the end only if no feasible flow
//! exists.
//!
//! Conventions (those of NetworkX's `network_simplex`, which `min_cost_flow`
//! and `max_flow_min_cost` call):
//! - Directed networks only.
//! - Capacities may be infinite; a missing capacity means infinite in
//!   NetworkX, so a binding passes `inf`. Edges of capacity 0 carry no
//!   flow.
//! - A self-loop carries its full capacity when its cost is negative and
//!   nothing otherwise; a negative-cost loop of infinite capacity makes
//!   the problem unbounded.
//! - Parallel edges are separate edges (NetworkX accepts multigraphs here).
//! - Errors: demands that do not sum to zero or cannot be met
//!   (NetworkX's `NetworkXUnfeasible`), a cycle of negative cost and
//!   infinite capacity (`NetworkXUnbounded`), negative capacities, and
//!   infinite costs or demands.
//! - Optimal flows are generally not unique; the cost is. The flow found
//!   is deterministic but can differ from NetworkX's.
//! - Arithmetic is in `f64`. With integer costs, capacities, and demands
//!   whose sums stay below 2^53 every quantity is an exact integer (the
//!   integrality theorem gives an integer optimal flow); with fractional
//!   data the result is exact up to rounding. NetworkX gives the same
//!   advice.

use super::maxflow::{check_capacities, maximum_flow};
use crate::error::{Error, Result};
use crate::graph::Graph;

/// A minimum-cost flow.
#[derive(Debug, Clone, PartialEq)]
pub struct CostFlow {
    /// Total cost.
    pub cost: f64,
    /// Flow on every edge, in input order.
    pub edge_flow: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Tree,
    Lower,
    Upper,
}

struct Simplex {
    /// Tail and head of every arc (real arcs, then one artificial arc per
    /// node).
    tail: Vec<usize>,
    head: Vec<usize>,
    cap: Vec<f64>,
    cost: Vec<f64>,
    flow: Vec<f64>,
    state: Vec<State>,
    // Tree on nodes 0..=n (n is the root).
    parent: Vec<usize>,
    parent_arc: Vec<usize>,
    depth: Vec<usize>,
    children: Vec<Vec<usize>>,
    pot: Vec<f64>,
    tol: f64,
}

impl Simplex {
    fn reduced(&self, a: usize) -> f64 {
        self.cost[a] + self.pot[self.tail[a]] - self.pot[self.head[a]]
    }

    /// Residual capacity of tree arc `a` in the direction from `from`.
    fn residual(&self, a: usize, from: usize) -> f64 {
        if self.tail[a] == from {
            self.cap[a] - self.flow[a]
        } else {
            self.flow[a]
        }
    }

    fn push(&mut self, a: usize, from: usize, delta: f64) {
        if self.tail[a] == from {
            self.flow[a] += delta;
        } else {
            self.flow[a] -= delta;
        }
    }

    /// Recomputes depth and potentials below `v` (inclusive of its
    /// children, given v's own values).
    fn refresh_below(&mut self, v: usize) {
        let mut stack = vec![v];
        while let Some(u) = stack.pop() {
            for i in 0..self.children[u].len() {
                let c = self.children[u][i];
                let a = self.parent_arc[c];
                self.depth[c] = self.depth[u] + 1;
                self.pot[c] = if self.tail[a] == u {
                    self.pot[u] + self.cost[a]
                } else {
                    self.pot[u] - self.cost[a]
                };
                stack.push(c);
            }
        }
    }

    fn detach(&mut self, child: usize) {
        let p = self.parent[child];
        let list = &mut self.children[p];
        let i = list.iter().position(|&c| c == child).expect("tree child");
        list.swap_remove(i);
    }

    /// One pivot on entering arc `e`. Returns an error when the cycle has
    /// no blocking arc (unbounded).
    fn pivot(&mut self, e: usize) -> Result<()> {
        let (first, second) = if self.state[e] == State::Lower {
            (self.tail[e], self.head[e])
        } else {
            (self.head[e], self.tail[e])
        };
        // Paths up to the apex.
        let (mut a, mut b) = (first, second);
        let mut up_first: Vec<usize> = Vec::new();
        let mut up_second: Vec<usize> = Vec::new();
        while a != b {
            if self.depth[a] >= self.depth[b] {
                up_first.push(a);
                a = self.parent[a];
            } else {
                up_second.push(b);
                b = self.parent[b];
            }
        }
        // Traversal in the direction of flow from the apex: down to
        // `first` (up_first reversed, each arc parent -> v), the entering
        // arc, then up from `second` (each arc v -> parent). The last arc
        // of least residual capacity leaves.
        let mut delta = f64::INFINITY;
        // (arc, node it is traversed from, the tree node below it or
        // usize::MAX for the entering arc)
        let mut leave = (e, first, usize::MAX);
        for &v in up_first.iter().rev() {
            let arc = self.parent_arc[v];
            let r = self.residual(arc, self.parent[v]);
            if r <= delta {
                delta = r;
                leave = (arc, self.parent[v], v);
            }
        }
        let r = if self.state[e] == State::Lower {
            self.cap[e]
        } else {
            self.flow[e]
        };
        if r <= delta {
            delta = r;
            leave = (e, first, usize::MAX);
        }
        for &v in &up_second {
            let arc = self.parent_arc[v];
            let r = self.residual(arc, v);
            if r <= delta {
                delta = r;
                leave = (arc, v, v);
            }
        }
        if delta == f64::INFINITY {
            return Err(Error::InvalidArgument(
                "a cycle of negative cost and infinite capacity makes the minimum cost \
                 unbounded"
                    .into(),
            ));
        }
        // Augment.
        if delta > 0.0 {
            for &v in &up_first {
                let arc = self.parent_arc[v];
                let p = self.parent[v];
                self.push(arc, p, delta);
            }
            self.push(e, first, delta);
            for &v in &up_second {
                let arc = self.parent_arc[v];
                self.push(arc, v, delta);
            }
        }
        let (out_arc, _, below) = leave;
        if out_arc == e {
            self.state[e] = if self.state[e] == State::Lower {
                State::Upper
            } else {
                State::Lower
            };
            return Ok(());
        }
        // The subtree under `below` is cut off; it contains `first` if the
        // leaving arc is on first's side, else `second`. Re-root it at
        // that endpoint and hang it from the other.
        let on_first_side = up_first.contains(&below);
        let (x, y) = if on_first_side {
            (first, second)
        } else {
            (second, first)
        };
        self.state[out_arc] =
            if self.flow[out_arc] <= self.tol.min(0.0) || self.flow[out_arc] == 0.0 {
                State::Lower
            } else if self.cap[out_arc] - self.flow[out_arc] == 0.0 {
                State::Upper
            } else if self.flow[out_arc] < self.cap[out_arc] - self.flow[out_arc] {
                State::Lower
            } else {
                State::Upper
            };
        if self.state[out_arc] == State::Lower {
            self.flow[out_arc] = 0.0;
        } else {
            self.flow[out_arc] = self.cap[out_arc];
        }
        // Reverse the path x -> ... -> below.
        self.detach(below);
        let mut path = vec![x];
        while *path.last().unwrap() != below {
            let v = *path.last().unwrap();
            path.push(self.parent[v]);
        }
        for i in (0..path.len() - 1).rev() {
            // path[i + 1] becomes a child of path[i].
            let (child, newp) = (path[i + 1], path[i]);
            self.detach(newp);
            self.parent[child] = newp;
            self.parent_arc[child] = self.parent_arc[newp];
            self.children[newp].push(child);
        }
        self.parent[x] = y;
        self.parent_arc[x] = e;
        self.children[y].push(x);
        self.state[e] = State::Tree;
        self.depth[x] = self.depth[y] + 1;
        self.pot[x] = if self.tail[e] == y {
            self.pot[y] + self.cost[e]
        } else {
            self.pot[y] - self.cost[e]
        };
        self.refresh_below(x);
        Ok(())
    }

    fn violation(&self, a: usize) -> f64 {
        match self.state[a] {
            State::Tree => 0.0,
            State::Lower => -self.reduced(a),
            State::Upper => self.reduced(a),
        }
    }

    fn solve(&mut self, arcs: usize) -> Result<()> {
        if arcs == 0 {
            return Ok(());
        }
        let block = ((arcs as f64).sqrt().ceil() as usize).max(1);
        let mut start = 0usize;
        loop {
            let mut entering = None;
            let mut scanned = 0usize;
            let mut pos = start;
            while scanned < arcs {
                let end = (scanned + block).min(arcs);
                let mut best = self.tol;
                while scanned < end {
                    let v = self.violation(pos);
                    if v > best {
                        best = v;
                        entering = Some(pos);
                    }
                    pos += 1;
                    if pos == arcs {
                        pos = 0;
                    }
                    scanned += 1;
                }
                if entering.is_some() {
                    break;
                }
            }
            match entering {
                None => return Ok(()),
                Some(e) => {
                    self.pivot(e)?;
                    start = pos;
                }
            }
        }
    }
}

/// Minimum-cost flow meeting `demand` (one value per node; see the module
/// documentation for the sign), with edge `capacity` (default infinite)
/// and `cost` (default 0) in input order.
pub fn min_cost_flow(
    g: &Graph,
    demand: &[f64],
    capacity: Option<&[f64]>,
    cost: Option<&[f64]>,
) -> Result<CostFlow> {
    if !g.is_directed() {
        return Err(Error::InvalidArgument(
            "minimum-cost flow is defined here for directed networks".into(),
        ));
    }
    let n = g.node_count();
    let m = g.edge_count();
    if n == 0 {
        return Err(Error::InvalidArgument("the network has no nodes".into()));
    }
    if demand.len() != n {
        return Err(Error::LengthMismatch {
            what: "demands",
            expected: n,
            found: demand.len(),
        });
    }
    if let Some(i) = demand.iter().position(|d| !d.is_finite()) {
        return Err(Error::InvalidArgument(format!(
            "node {i} has demand {}; demands must be finite",
            demand[i]
        )));
    }
    if let Some(c) = capacity {
        check_capacities(c, m, "capacities")?;
    }
    if let Some(c) = cost {
        if c.len() != m {
            return Err(Error::LengthMismatch {
                what: "costs",
                expected: m,
                found: c.len(),
            });
        }
        if let Some(i) = c.iter().position(|x| !x.is_finite()) {
            return Err(Error::InvalidArgument(format!(
                "edge {i} has cost {}; costs must be finite",
                c[i]
            )));
        }
    }
    let total: f64 = demand.iter().sum();
    let scale: f64 = demand.iter().map(|d| d.abs()).fold(0.0, f64::max);
    if total.abs() > 1e-12 * scale.max(1.0) * n as f64 {
        return Err(Error::InvalidArgument(format!(
            "the demands sum to {total}, not 0, so no flow meets them"
        )));
    }
    let cap_of = |e: usize| capacity.map_or(f64::INFINITY, |c| c[e]);
    let cost_of = |e: usize| cost.map_or(0.0, |c| c[e]);
    let (src, tgt) = (g.edge_sources(), g.edge_targets());

    let mut edge_flow = vec![0.0f64; m];
    let mut loop_cost = 0.0;
    let mut real: Vec<usize> = Vec::new();
    for e in 0..m {
        if src[e] == tgt[e] {
            if cost_of(e) < 0.0 {
                if cap_of(e) == f64::INFINITY {
                    return Err(Error::InvalidArgument(format!(
                        "self-loop {e} has negative cost and infinite capacity, so the minimum \
                         cost is unbounded"
                    )));
                }
                edge_flow[e] = cap_of(e);
                loop_cost += cost_of(e) * cap_of(e);
            }
        } else if cap_of(e) > 0.0 {
            real.push(e);
        }
    }
    let r = n;
    let k = real.len();
    let max_cost = real.iter().map(|&e| cost_of(e).abs()).fold(0.0, f64::max);
    let big = 1.0 + (n as f64) * max_cost.max(1.0);
    let mut s = Simplex {
        tail: Vec::with_capacity(k + n),
        head: Vec::with_capacity(k + n),
        cap: Vec::with_capacity(k + n),
        cost: Vec::with_capacity(k + n),
        flow: vec![0.0; k + n],
        state: vec![State::Lower; k + n],
        parent: vec![r; n + 1],
        parent_arc: vec![usize::MAX; n + 1],
        depth: vec![1; n + 1],
        children: vec![Vec::new(); n + 1],
        pot: vec![0.0; n + 1],
        tol: 1e-12 * big,
    };
    for &e in &real {
        s.tail.push(src[e] as usize);
        s.head.push(tgt[e] as usize);
        s.cap.push(cap_of(e));
        s.cost.push(cost_of(e));
    }
    s.depth[r] = 0;
    for (v, &d) in demand.iter().enumerate() {
        let a = k + v;
        if d > 0.0 {
            s.tail.push(r);
            s.head.push(v);
            s.pot[v] = big;
        } else {
            s.tail.push(v);
            s.head.push(r);
            s.pot[v] = -big;
        }
        s.cap.push(f64::INFINITY);
        s.cost.push(big);
        s.flow[a] = d.abs();
        s.state[a] = State::Tree;
        s.parent_arc[v] = a;
        s.children[r].push(v);
    }
    s.solve(k + n)?;
    let feas_tol = 1e-9 * scale.max(1.0);
    if (0..n).any(|v| s.flow[k + v] > feas_tol) {
        return Err(Error::InvalidArgument(
            "no flow satisfies every node's demand within the capacities".into(),
        ));
    }
    let mut total_cost = loop_cost;
    for (i, &e) in real.iter().enumerate() {
        edge_flow[e] = s.flow[i];
        total_cost += s.cost[i] * s.flow[i];
    }
    Ok(CostFlow {
        cost: total_cost,
        edge_flow,
    })
}

/// A maximum flow from `source` to `sink` of least cost (NetworkX's
/// `max_flow_min_cost`): the maximum flow value F by
/// [`crate::flow::maxflow::maximum_flow`], then [`min_cost_flow`] with
/// demand -F at the source and F at the sink. Returns the cost flow and F.
pub fn max_flow_min_cost(
    g: &Graph,
    source: usize,
    sink: usize,
    capacity: Option<&[f64]>,
    cost: Option<&[f64]>,
) -> Result<(CostFlow, f64)> {
    if !g.is_directed() {
        return Err(Error::InvalidArgument(
            "minimum-cost flow is defined here for directed networks".into(),
        ));
    }
    // NetworkX reads a missing capacity as infinite in both steps.
    let inf;
    let cap = match capacity {
        Some(c) => c,
        None => {
            inf = vec![f64::INFINITY; g.edge_count()];
            &inf[..]
        }
    };
    let value = maximum_flow(g, &[source], &[sink], Some(cap))?.value;
    let mut demand = vec![0.0; g.node_count()];
    demand[source] = -value;
    demand[sink] = value;
    Ok((min_cost_flow(g, &demand, Some(cap), cost)?, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn digraph(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, true).unwrap()
    }

    #[test]
    fn networkx_doc_example() {
        // network_simplex docstring: a -> b (cost 3, cap 4), a -> c (6,
        // 10), b -> d (1, 9), c -> d (2, 5); demand a = -5, d = 5. Cost 24.
        let g = digraph(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        let f = min_cost_flow(
            &g,
            &[-5.0, 0.0, 0.0, 5.0],
            Some(&[4.0, 10.0, 9.0, 5.0]),
            Some(&[3.0, 6.0, 1.0, 2.0]),
        )
        .unwrap();
        assert_eq!(f.cost, 24.0);
        assert_eq!(f.edge_flow, vec![4.0, 1.0, 4.0, 1.0]);
    }

    #[test]
    fn infeasible_unbounded_and_loops() {
        let g = digraph(3, &[(0, 1), (1, 2), (2, 1), (1, 1)]);
        // Demand exceeds capacity.
        assert!(min_cost_flow(&g, &[-3.0, 0.0, 3.0], Some(&[2.0, 2.0, 2.0, 1.0]), None).is_err());
        // Negative cycle 1 -> 2 -> 1 of infinite capacity.
        assert!(min_cost_flow(&g, &[0.0; 3], None, Some(&[0.0, -1.0, 0.0, 0.0])).is_err());
        // Negative loop with finite capacity is saturated.
        let f = min_cost_flow(
            &g,
            &[-1.0, 0.0, 1.0],
            Some(&[5.0, 5.0, 5.0, 2.0]),
            Some(&[1.0, 1.0, 1.0, -3.0]),
        )
        .unwrap();
        assert_eq!(f.cost, 2.0 - 6.0);
        assert_eq!(f.edge_flow, vec![1.0, 1.0, 0.0, 2.0]);
        // Demands not summing to zero.
        assert!(min_cost_flow(&g, &[1.0, 0.0, 0.0], None, None).is_err());
    }

    #[test]
    fn max_flow_of_least_cost() {
        // Two routes of capacity 2 from 0 to 3, costs 1 + 1 and 5 + 5, and
        // a cheap route of capacity 1.
        let g = digraph(4, &[(0, 1), (1, 3), (0, 2), (2, 3), (0, 3)]);
        let (f, v) = max_flow_min_cost(
            &g,
            0,
            3,
            Some(&[2.0, 2.0, 2.0, 2.0, 1.0]),
            Some(&[1.0, 1.0, 5.0, 5.0, 0.0]),
        )
        .unwrap();
        assert_eq!(v, 5.0);
        assert_eq!(f.cost, 4.0 + 20.0);
    }
}
