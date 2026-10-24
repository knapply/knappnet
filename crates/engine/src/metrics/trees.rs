//! Spanning trees of undirected networks, and optimum branchings and
//! spanning arborescences of directed ones.
//!
//! Spanning trees (forests): a minimum spanning forest has least total
//! weight among forests with a tree spanning each component. Three
//! algorithms are offered, which return the same forest: Kruskal's (1956,
//! `kruskal1956shortest`; edges in increasing order, joined unless they
//! close a cycle, with [`crate::union_find`]), Prim's (1957,
//! `prim1957shortest`; grow a tree from each component's lowest node with
//! a heap of crossing edges), and Boruvka's (1926, `boruvka1926jistem`;
//! every component adds its lightest crossing edge in rounds). Edges are
//! compared by (weight, input position), a strict total order, under which
//! the minimum spanning forest is unique, so all three algorithms agree
//! edge for edge whatever the ties in weight. That forest is also the one
//! NetworkX's default `minimum_spanning_edges(algorithm = "kruskal")`
//! (`hagberg2008exploring`) returns when its edges come in the same order,
//! since NetworkX sorts them stably by weight. A maximum spanning forest
//! compares by (minus weight, input position), as NetworkX's stable sort
//! in reverse does.
//!
//! Branchings and arborescences: a branching is a set of arcs with at most
//! one arc into each node and no cycle (a forest of out-trees); a spanning
//! arborescence is a branching with one root that reaches every node.
//! Optimum ones are found by Edmonds' (1967, `edmonds1967optimum`; also
//! Chu and Liu 1965, `chu1965shortest`) algorithm in Tarjan's (1977,
//! `tarjan1977finding`; Camerini, Fratta, and Maffioli 1979,
//! `camerini1979note`) O(m log n) form: every node takes its cheapest
//! incoming arc; a cycle formed this way is contracted into one node
//! whose incoming arcs are repriced by the cost of the cycle arc they
//! would replace (leftist heaps with lazy additions hold each node's
//! incoming arcs, so a contraction merges heaps; Gabow, Galil, Spencer,
//! and Tarjan 1986, `gabow1986efficient`); the contractions are undone in
//! reverse to recover the arcs. Problems are reduced to a minimum
//! arborescence rooted at an artificial root joined to every node:
//! - maximum branching: arcs of positive weight w at cost -w, artificial
//!   arcs at cost 0 (an artificial arc marks a root of the branching);
//! - minimum branching: arcs of negative weight at cost w;
//! - minimum (maximum) spanning arborescence: arcs at cost w (-w), and
//!   artificial arcs at a cost that exceeds any total of real costs (a
//!   lexicographic pair: number of roots, then weight), so the fewest
//!   roots are used first and a spanning arborescence exists exactly when
//!   one artificial arc suffices;
//! - an arborescence rooted at a given node is solved from that node
//!   directly.
//!
//! Conventions:
//! - Weights are edge weights (with `weighted`) or 1. Self-loops never
//!   belong to a tree or branching. Parallel edges are separate
//!   candidates (NetworkX accepts multigraphs here). Spanning trees are
//!   for undirected networks, branchings for directed ones.
//! - Edges are reported by input position in increasing order (spanning
//!   trees: in increasing comparison order, the order Kruskal's algorithm
//!   accepts them), with their total weight.
//! - The optimum weight is unique; the arcs achieving it need not be.
//!   Among equal costs the arc earlier in the input is preferred, so the
//!   result is deterministic, but it can differ from NetworkX's choice.
//! - NetworkX computes spanning arborescences by shifting weights and
//!   calling its maximum branching, which does not always yield a
//!   spanning arborescence even when one exists (it then raises
//!   `NetworkXException`); this module always finds the optimum.
//!   `BranchingKind::Minimal` reproduces NetworkX's `minimal_branching`,
//!   which is a maximum branching under the weights (2 max - min + 1) - w,
//!   max and min taken over all edges.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::union_find::UnionFind;

/// Spanning tree algorithm. All return the same forest (see the module
/// documentation).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SpanningAlgorithm {
    /// Kruskal's algorithm.
    #[default]
    Kruskal,
    /// Prim's algorithm.
    Prim,
    /// Boruvka's algorithm.
    Boruvka,
}

impl SpanningAlgorithm {
    /// Decodes the C interface codes: 0 = Kruskal, 1 = Prim, 2 = Boruvka.
    pub fn from_code(code: i32) -> Result<SpanningAlgorithm> {
        match code {
            0 => Ok(SpanningAlgorithm::Kruskal),
            1 => Ok(SpanningAlgorithm::Prim),
            2 => Ok(SpanningAlgorithm::Boruvka),
            _ => Err(Error::InvalidArgument(format!(
                "unknown spanning tree algorithm code {code}"
            ))),
        }
    }
}

/// A set of edges with its total weight.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeSet {
    /// Input positions of the edges.
    pub edges: Vec<usize>,
    /// Their total weight.
    pub weight: f64,
}

fn edge_weights(g: &Graph, weighted: bool) -> Result<Vec<f64>> {
    if weighted {
        Ok(g.edge_weights().ok_or(Error::MissingWeights)?.to_vec())
    } else {
        Ok(vec![1.0; g.edge_count()])
    }
}

/// Minimum (or, with `maximum`, maximum) spanning forest of an undirected
/// network. See the module documentation.
pub fn spanning_tree(
    g: &Graph,
    weighted: bool,
    maximum: bool,
    algorithm: SpanningAlgorithm,
) -> Result<EdgeSet> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(
            "spanning trees are defined here for undirected networks; see the arborescence \
             functions for directed ones"
                .into(),
        ));
    }
    let n = g.node_count();
    let m = g.edge_count();
    let w = edge_weights(g, weighted)?;
    let sign = if maximum { -1.0 } else { 1.0 };
    let key = |e: usize| (sign * w[e], e);
    let less = |a: usize, b: usize| {
        let (ka, kb) = (key(a), key(b));
        ka.0 < kb.0 || (ka.0 == kb.0 && ka.1 < kb.1)
    };
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut chosen: Vec<usize> = match algorithm {
        SpanningAlgorithm::Kruskal => {
            let mut order: Vec<usize> = (0..m).filter(|&e| src[e] != tgt[e]).collect();
            order.sort_by(|&a, &b| {
                key(a)
                    .0
                    .partial_cmp(&key(b).0)
                    .expect("finite weights")
                    .then(a.cmp(&b))
            });
            let mut sets = UnionFind::new(n);
            order
                .into_iter()
                .filter(|&e| sets.union(src[e] as usize, tgt[e] as usize).is_some())
                .collect()
        }
        SpanningAlgorithm::Prim => {
            // Incidences with the edge behind each, loops left out.
            let inc = crate::order::rows(
                g,
                crate::order::Direction::Out,
                crate::order::NeighborOrder::Sorted,
                true,
                false,
            );
            let mut in_tree = vec![false; n];
            let mut chosen = Vec::new();
            let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<Key>> =
                std::collections::BinaryHeap::new();
            for root in 0..n {
                if in_tree[root] {
                    continue;
                }
                in_tree[root] = true;
                let push = |heap: &mut std::collections::BinaryHeap<std::cmp::Reverse<Key>>,
                            v: usize| {
                    for (&x, &e) in inc.neighbors(v).iter().zip(inc.edges(v)) {
                        let e = e as usize;
                        heap.push(std::cmp::Reverse(Key {
                            w: sign * w[e],
                            e,
                            to: x as usize,
                        }));
                    }
                };
                push(&mut heap, root);
                while let Some(std::cmp::Reverse(k)) = heap.pop() {
                    if in_tree[k.to] {
                        continue;
                    }
                    in_tree[k.to] = true;
                    chosen.push(k.e);
                    push(&mut heap, k.to);
                }
            }
            chosen
        }
        SpanningAlgorithm::Boruvka => {
            let mut sets = UnionFind::new(n);
            let mut chosen = Vec::new();
            loop {
                let mut best = vec![usize::MAX; n];
                for e in 0..m {
                    let (a, b) = (sets.find(src[e] as usize), sets.find(tgt[e] as usize));
                    if a == b {
                        continue;
                    }
                    for c in [a, b] {
                        if best[c] == usize::MAX || less(e, best[c]) {
                            best[c] = e;
                        }
                    }
                }
                let mut added = false;
                for &e in &best {
                    if e != usize::MAX && sets.union(src[e] as usize, tgt[e] as usize).is_some() {
                        chosen.push(e);
                        added = true;
                    }
                }
                if !added {
                    break;
                }
            }
            chosen
        }
    };
    chosen.sort_by(|&a, &b| {
        key(a)
            .0
            .partial_cmp(&key(b).0)
            .expect("finite weights")
            .then(a.cmp(&b))
    });
    let weight = chosen.iter().map(|&e| w[e]).sum();
    Ok(EdgeSet {
        edges: chosen,
        weight,
    })
}

/// Heap entry for Prim's algorithm, ordered by (weight, edge).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Key {
    w: f64,
    e: usize,
    to: usize,
}

impl Eq for Key {}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.w
            .partial_cmp(&other.w)
            .expect("finite weights")
            .then(self.e.cmp(&other.e))
            .then(self.to.cmp(&other.to))
    }
}

/// Which optimum branching or arborescence to find. See the module
/// documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchingKind {
    /// Branching of greatest total weight (only positive arcs help).
    MaximumBranching,
    /// Branching of least total weight (only negative arcs help).
    MinimumBranching,
    /// NetworkX's `minimal_branching`: a maximum branching under the
    /// weights (2 max - min + 1) - w.
    Minimal,
    /// Spanning arborescence of greatest weight.
    MaximumArborescence,
    /// Spanning arborescence of least weight.
    MinimumArborescence,
}

impl BranchingKind {
    /// Decodes the C interface codes: 0 = maximum branching, 1 = minimum
    /// branching, 2 = minimal (NetworkX), 3 = maximum spanning
    /// arborescence, 4 = minimum spanning arborescence.
    pub fn from_code(code: i32) -> Result<BranchingKind> {
        match code {
            0 => Ok(BranchingKind::MaximumBranching),
            1 => Ok(BranchingKind::MinimumBranching),
            2 => Ok(BranchingKind::Minimal),
            3 => Ok(BranchingKind::MaximumArborescence),
            4 => Ok(BranchingKind::MinimumArborescence),
            _ => Err(Error::InvalidArgument(format!(
                "unknown branching kind code {code}"
            ))),
        }
    }
}

/// A cost: a count of artificial arcs (roots), then a weight, compared
/// lexicographically.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
struct Cost {
    roots: i64,
    w: f64,
}

impl std::ops::Add for Cost {
    type Output = Cost;
    fn add(self, o: Cost) -> Cost {
        Cost {
            roots: self.roots + o.roots,
            w: self.w + o.w,
        }
    }
}

impl std::ops::Sub for Cost {
    type Output = Cost;
    fn sub(self, o: Cost) -> Cost {
        Cost {
            roots: self.roots - o.roots,
            w: self.w - o.w,
        }
    }
}

struct Arc {
    from: u32,
    to: u32,
    /// Input position, or `ARTIFICIAL` for an arc from the added root.
    edge: u32,
    cost: Cost,
}

const ARTIFICIAL: u32 = u32::MAX;
const NIL: usize = usize::MAX;
const NIL32: u32 = u32::MAX;

/// Leftist heaps of arcs with lazy cost additions, in one arena; node `i`
/// holds arc `i` (32-bit links keep the arena small).
struct Heaps {
    key: Vec<Cost>,
    arc: Vec<u32>,
    left: Vec<u32>,
    right: Vec<u32>,
    rank: Vec<u32>,
    delta: Vec<Cost>,
}

impl Heaps {
    fn new(cap: usize) -> Heaps {
        Heaps {
            key: Vec::with_capacity(cap),
            arc: Vec::with_capacity(cap),
            left: Vec::with_capacity(cap),
            right: Vec::with_capacity(cap),
            rank: Vec::with_capacity(cap),
            delta: Vec::with_capacity(cap),
        }
    }

    fn single(&mut self, arc: usize, key: Cost) -> usize {
        self.key.push(key);
        self.arc.push(arc as u32);
        self.left.push(NIL32);
        self.right.push(NIL32);
        self.rank.push(1);
        self.delta.push(Cost::default());
        self.key.len() - 1
    }

    fn link(x: u32) -> usize {
        if x == NIL32 {
            NIL
        } else {
            x as usize
        }
    }

    fn settle(&mut self, a: usize) {
        let d = self.delta[a];
        if d == Cost::default() {
            return;
        }
        self.key[a] = self.key[a] + d;
        for c in [self.left[a], self.right[a]] {
            if c != NIL32 {
                let c = c as usize;
                self.delta[c] = self.delta[c] + d;
            }
        }
        self.delta[a] = Cost::default();
    }

    fn rank_of(&self, a: usize) -> u32 {
        if a == NIL {
            0
        } else {
            self.rank[a]
        }
    }

    /// Whether `a` should be above `b`: smaller cost, then earlier arc.
    fn before(&self, a: usize, b: usize, arcs: &[Arc]) -> bool {
        let (ka, kb) = (self.key[a], self.key[b]);
        ka < kb || (ka == kb && arcs[self.arc[a] as usize].edge < arcs[self.arc[b] as usize].edge)
    }

    fn merge(&mut self, a: usize, b: usize, arcs: &[Arc]) -> usize {
        if a == NIL {
            return b;
        }
        if b == NIL {
            return a;
        }
        self.settle(a);
        self.settle(b);
        let (a, b) = if self.before(b, a, arcs) {
            (b, a)
        } else {
            (a, b)
        };
        let r = self.merge(Self::link(self.right[a]), b, arcs);
        self.right[a] = r as u32;
        if self.rank_of(Self::link(self.left[a])) < self.rank_of(r) {
            self.right[a] = self.left[a];
            self.left[a] = r as u32;
        }
        self.rank[a] = self.rank_of(Self::link(self.right[a])) + 1;
        a
    }

    fn pop(&mut self, a: usize, arcs: &[Arc]) -> usize {
        self.settle(a);
        self.merge(Self::link(self.left[a]), Self::link(self.right[a]), arcs)
    }
}

/// Union-find with undo (union by size, no path compression).
struct Rollback {
    parent: Vec<usize>,
    size: Vec<usize>,
    history: Vec<(usize, usize)>,
}

impl Rollback {
    fn new(n: usize) -> Rollback {
        Rollback {
            parent: (0..n).collect(),
            size: vec![1; n],
            history: Vec::new(),
        }
    }

    fn find(&self, mut v: usize) -> usize {
        while self.parent[v] != v {
            v = self.parent[v];
        }
        v
    }

    fn join(&mut self, a: usize, b: usize) -> bool {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return false;
        }
        if self.size[a] < self.size[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b] = a;
        self.size[a] += self.size[b];
        self.history.push((a, b));
        true
    }

    fn time(&self) -> usize {
        self.history.len()
    }

    fn rollback(&mut self, t: usize) {
        while self.history.len() > t {
            let (a, b) = self.history.pop().expect("history");
            self.parent[b] = b;
            self.size[a] -= self.size[b];
        }
    }
}

/// Minimum-cost arborescence rooted at `root` on nodes `0..nodes`: the
/// chosen arc into every other node, or `None` if some node cannot be
/// reached.
fn min_arborescence(nodes: usize, root: usize, arcs: &[Arc]) -> Option<Vec<usize>> {
    let mut heaps = Heaps::new(arcs.len());
    let mut heap = vec![NIL; nodes];
    for (i, a) in arcs.iter().enumerate() {
        if a.from == a.to || a.to as usize == root {
            continue;
        }
        let h = heaps.single(i, a.cost);
        heap[a.to as usize] = heaps.merge(heap[a.to as usize], h, arcs);
    }
    let mut uf = Rollback::new(nodes);
    const NONE: usize = usize::MAX;
    let mut seen = vec![NONE; nodes];
    seen[root] = root;
    let mut in_arc = vec![NONE; nodes];
    // (contracted node, union-find time before, arcs of the cycle)
    let mut cycles: Vec<(usize, usize, Vec<usize>)> = Vec::new();
    for s in 0..nodes {
        let mut u = s;
        // Nodes on the current path with the arc each one chose.
        let mut path: Vec<(usize, usize)> = Vec::new();
        while seen[u] == NONE {
            // The cheapest incoming arc from outside u's own contraction.
            let a = loop {
                let h = heap[u];
                if h == NIL {
                    return None;
                }
                heaps.settle(h);
                let a = heaps.arc[h] as usize;
                let cost = heaps.key[h];
                heap[u] = heaps.pop(h, arcs);
                if uf.find(arcs[a].from as usize) == u {
                    continue;
                }
                // Reprice the remaining incoming arcs relative to this one.
                if heap[u] != NIL {
                    let r = heap[u];
                    heaps.delta[r] = heaps.delta[r] - cost;
                }
                break a;
            };
            path.push((u, a));
            seen[u] = s;
            u = uf.find(arcs[a].from as usize);
            if seen[u] == s {
                // A cycle: contract it into one node.
                let t = uf.time();
                let mut merged = NIL;
                let mut cycle = Vec::new();
                loop {
                    let (w, wa) = path.pop().expect("cycle on the path");
                    merged = heaps.merge(merged, heap[w], arcs);
                    cycle.push(wa);
                    if !uf.join(u, w) {
                        break;
                    }
                }
                u = uf.find(u);
                heap[u] = merged;
                seen[u] = NONE;
                cycles.push((u, t, cycle));
            }
        }
        for &(_, a) in &path {
            let v = uf.find(arcs[a].to as usize);
            in_arc[v] = a;
        }
    }
    for (u, t, cycle) in cycles.into_iter().rev() {
        uf.rollback(t);
        let into = in_arc[u];
        for a in cycle {
            let v = uf.find(arcs[a].to as usize);
            in_arc[v] = a;
        }
        let v = uf.find(arcs[into].to as usize);
        in_arc[v] = into;
    }
    Some(in_arc)
}

/// Optimum branching or spanning arborescence of a directed network. With
/// `root`, the arborescence (for the arborescence kinds) must be rooted
/// at that node. See the module documentation.
pub fn optimum_branching(
    g: &Graph,
    weighted: bool,
    kind: BranchingKind,
    root: Option<usize>,
) -> Result<EdgeSet> {
    if !g.is_directed() {
        return Err(Error::InvalidArgument(
            "branchings and arborescences are defined for directed networks".into(),
        ));
    }
    let n = g.node_count();
    let w = edge_weights(g, weighted)?;
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let arborescence = matches!(
        kind,
        BranchingKind::MaximumArborescence | BranchingKind::MinimumArborescence
    );
    if let Some(r) = root {
        if !arborescence {
            return Err(Error::InvalidArgument(
                "a root can only be given for a spanning arborescence".into(),
            ));
        }
        if r >= n {
            return Err(Error::InvalidArgument(format!(
                "root {r} is not a node (n = {n})"
            )));
        }
    }
    if arborescence && n == 0 {
        return Err(Error::InvalidArgument(
            "a network with no nodes has no spanning arborescence".into(),
        ));
    }
    // Cost of each real arc, or None if it cannot help.
    let shift = if kind == BranchingKind::Minimal {
        let max = w.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min = w.iter().copied().fold(f64::INFINITY, f64::min);
        2.0 * max - min + 1.0
    } else {
        0.0
    };
    let cost_of = |e: usize| -> Option<f64> {
        let x = w[e];
        match kind {
            BranchingKind::MaximumBranching => (x > 0.0).then_some(-x),
            BranchingKind::MinimumBranching => (x < 0.0).then_some(x),
            BranchingKind::Minimal => {
                let y = shift - x;
                (y > 0.0).then_some(-y)
            }
            BranchingKind::MaximumArborescence => Some(-x),
            BranchingKind::MinimumArborescence => Some(x),
        }
    };
    let mut arcs: Vec<Arc> = Vec::new();
    for e in 0..g.edge_count() {
        if src[e] == tgt[e] {
            continue;
        }
        if let Some(c) = cost_of(e) {
            arcs.push(Arc {
                from: src[e],
                to: tgt[e],
                edge: e as u32,
                cost: Cost { roots: 0, w: c },
            });
        }
    }
    let (nodes, top) = match root {
        Some(r) => (n, r),
        None => {
            let extra = Cost {
                roots: if arborescence { 1 } else { 0 },
                w: 0.0,
            };
            for v in 0..n {
                arcs.push(Arc {
                    from: n as u32,
                    to: v as u32,
                    edge: ARTIFICIAL,
                    cost: extra,
                });
            }
            (n + 1, n)
        }
    };
    let chosen = min_arborescence(nodes, top, &arcs).ok_or_else(|| {
        Error::InvalidArgument(match root {
            Some(r) => format!("no spanning arborescence is rooted at node {r}"),
            None => "the network has no spanning arborescence".into(),
        })
    })?;
    let mut edges: Vec<usize> = Vec::new();
    let mut roots = 0;
    for v in 0..n {
        if v == top {
            continue;
        }
        let a = &arcs[chosen[v]];
        if a.edge == ARTIFICIAL {
            roots += 1;
        } else {
            edges.push(a.edge as usize);
        }
    }
    if arborescence && root.is_none() && roots > 1 {
        return Err(Error::InvalidArgument(
            "the network has no spanning arborescence".into(),
        ));
    }
    edges.sort_unstable();
    let weight = edges.iter().map(|&e| w[e]).sum();
    Ok(EdgeSet { edges, weight })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], w: &[f64], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, Some(w), directed).unwrap()
    }

    #[test]
    fn spanning_trees_agree() {
        // A square with a diagonal and ties.
        let g = graph(
            5,
            &[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (3, 3)],
            &[1.0, 2.0, 1.0, 2.0, 1.0, 0.0],
            false,
        );
        for alg in [
            SpanningAlgorithm::Kruskal,
            SpanningAlgorithm::Prim,
            SpanningAlgorithm::Boruvka,
        ] {
            let t = spanning_tree(&g, true, false, alg).unwrap();
            assert_eq!(t.edges, vec![0, 2, 4]);
            assert_eq!(t.weight, 3.0);
            let t = spanning_tree(&g, true, true, alg).unwrap();
            assert_eq!(t.edges, vec![1, 3, 0]);
            assert_eq!(t.weight, 5.0);
        }
    }

    #[test]
    fn branchings() {
        // Cycle 0 -> 1 -> 2 -> 0 with weights 5, 4, 3, and 3 -> 1 (1).
        let g = graph(
            4,
            &[(0, 1), (1, 2), (2, 0), (3, 1), (3, 3)],
            &[5.0, 4.0, 3.0, 1.0, 9.0],
            true,
        );
        let b = optimum_branching(&g, true, BranchingKind::MaximumBranching, None).unwrap();
        assert_eq!(b.edges, vec![0, 1]);
        assert_eq!(b.weight, 9.0);
        let a = optimum_branching(&g, true, BranchingKind::MaximumArborescence, None).unwrap();
        // Root 3: 3 -> 1, 1 -> 2, 2 -> 0.
        assert_eq!(a.edges, vec![1, 2, 3]);
        assert_eq!(a.weight, 8.0);
        let a = optimum_branching(&g, true, BranchingKind::MinimumArborescence, None).unwrap();
        assert_eq!(a.weight, 8.0);
        assert!(optimum_branching(&g, true, BranchingKind::MinimumArborescence, Some(0)).is_err());
        let b = optimum_branching(&g, true, BranchingKind::MinimumBranching, None).unwrap();
        assert!(b.edges.is_empty());
    }
}
