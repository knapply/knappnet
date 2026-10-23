//! Reach, edge betweenness, information centrality, and bridging.
//!
//! Sources:
//! - **k-step reach** is the number of other nodes within k steps of a
//!   node, the measure UCINET calls reach centrality and Borgatti (2006,
//!   `borgatti2006identifying`) uses for sets of key players; its
//!   distance-weighted form (the sum of 1/d over the other nodes, divided
//!   by n - 1 as the average reciprocal distance) is harmonic centrality,
//!   computed by [`crate::metrics::closeness::harmonic_centrality`].
//! - **Edge betweenness** is the number of shortest paths through an edge,
//!   summed over pairs of nodes with each path weighted by the share of the
//!   pair's shortest paths it represents (Girvan and Newman 2002,
//!   `girvan2002community`, after Freeman 1977, `freeman1977set`). It is
//!   computed with Brandes' accumulation for edges (Brandes 2008,
//!   `brandes2008variants`, section 3.4): after a single-source search from
//!   s, an arc (v, w) on a shortest path from s receives
//!   (sigma_sv / sigma_sw)(1 + delta_s(w)).
//! - **Information centrality** (Stephenson and Zelen 1989,
//!   `stephenson1989rethinking`): with A the symmetric (possibly valued)
//!   adjacency matrix and B the matrix with b_ij = 1 - a_ij off the
//!   diagonal and b_ii = 1 + sum_j a_ij on it (B = L + J, the Laplacian plus
//!   a matrix of ones), C = B^-1 and I_i = 1 / (c_ii + (T - 2R)/n), where T
//!   is the trace of C and R any row sum of C (all are equal, to 1/n). This
//!   is the formula documented for sna's `infocent()` (`butts2008sna`).
//!   NetworkX's `information_centrality()` returns I_i / n (it computes
//!   current-flow closeness, which Brandes and Fleischer 2005,
//!   `brandes2005centrality`, show to be information centrality up to that
//!   factor).
//! - **Bridging coefficient** (Hwang, Kim, Ramanathan, and Zhang 2008,
//!   `hwang2008bridging`): BC(v) = (1/d(v)) / sum_{u in N(v)} 1/d(u), where
//!   d is the number of distinct neighbors. Bridging centrality is
//!   betweenness times the bridging coefficient; the product is formed by
//!   the language bindings from [`crate::metrics::betweenness::betweenness`]
//!   and [`bridging_coefficient`].
//!
//! Conventions:
//! - Reach: distances follow `mode` (`Out`: from the node, `In`: to it,
//!   `All`: direction ignored) and the weights as described in
//!   [`crate::weights`]: unweighted networks count steps, weighted ones use
//!   lengths (weights as distances, or strengths converted to 1/w), so k is
//!   then a length. A node is within k when its distance is at most k, up
//!   to a relative tolerance of 1e-10 (the tie tolerance of the weighted
//!   path searches). The node itself is not counted. Normalized reach
//!   divides by n - 1 (UCINET reports reach as a proportion of n - 1);
//!   it is NaN for a single node.
//! - Edge betweenness: one value per edge, in input order. Parallel edges
//!   are distinct edges and distinct shortest paths, and self-loops never
//!   lie on a shortest path and score 0, as in igraph's
//!   `edge_betweenness()`. In undirected networks each unordered pair of
//!   endpoints is counted once (the raw sum is halved), as in igraph and
//!   NetworkX. Weights are lengths, as for node betweenness, with ties
//!   between real path lengths detected to a relative tolerance of 1e-10.
//!   Normalized values divide by the number of ordered pairs of nodes,
//!   n(n - 1), in directed networks and by the number of unordered pairs,
//!   n(n - 1)/2, in undirected ones, as NetworkX's
//!   `edge_betweenness_centrality(normalized = True)` does.
//! - Information centrality is defined for undirected ties: a directed
//!   network is symmetrized with a_ij = max(a_ij, a_ji) (sna's default
//!   `cmode = "weak"` for binary data). Self-loops are ignored. Without
//!   weights a_ij is 1 when the nodes are tied, so parallel edges count
//!   once (the adjacency-matrix reading of sna); with weights a_ij is the
//!   total weight of the edges between i and j, read as a tie strength (a
//!   conductance), used as given, and required to be strictly positive.
//!   Following sna, isolated nodes score 0 and are removed first, so n in
//!   the formula is the number of nodes with at least one tie. The formula
//!   needs the remaining network to be connected (otherwise B is singular);
//!   when it is not, every node with ties gets NaN.
//! - Bridging coefficient: direction, weights, self-loops, and parallel
//!   edges are ignored (d is the number of distinct neighbors other than
//!   the node itself). It is NaN for an isolated node.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;
use crate::metrics::descriptive::{distance_rows, PathOptions};
use crate::parallel;
use crate::weights::{path_lengths, Lengths, PathAlgorithm, WeightKind};

const EPSILON: f64 = 1e-10;

// ---- Reach --------------------------------------------------------------------

/// Options for [`reach`].
#[derive(Debug, Clone, Copy)]
pub struct ReachOptions {
    /// Largest distance counted: a number of steps, or a length when
    /// `weighted` is true. May be infinite.
    pub k: f64,
    /// Direction of the paths in directed graphs.
    pub mode: DegreeMode,
    /// Use edge weights as lengths (they must be strictly positive).
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Divide by n - 1.
    pub normalized: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for ReachOptions {
    fn default() -> Self {
        ReachOptions {
            k: 2.0,
            mode: DegreeMode::Out,
            weighted: false,
            weights_as: WeightKind::Distance,
            normalized: false,
            threads: 0,
        }
    }
}

/// Number (or, normalized, proportion) of other nodes within distance `k`
/// of every node, in node order. See the module documentation.
pub fn reach(g: &Graph, opts: ReachOptions) -> Result<Vec<f64>> {
    if opts.k.is_nan() || opts.k < 0.0 {
        return Err(Error::InvalidArgument(format!(
            "the reach distance k must be a non-negative number, not {}",
            opts.k
        )));
    }
    let n = g.node_count();
    let limit = if opts.k.is_finite() {
        opts.k + EPSILON * opts.k.abs().max(1.0)
    } else {
        f64::INFINITY
    };
    let mut out = vec![0.0; n];
    let popts = PathOptions {
        mode: opts.mode,
        weighted: opts.weighted,
        weights_as: opts.weights_as,
        algorithm: PathAlgorithm::Auto,
        threads: opts.threads,
    };
    distance_rows(g, None, popts, |i, row| {
        let count = row
            .iter()
            .enumerate()
            .filter(|&(j, &d)| j != i && d.is_finite() && d <= limit)
            .count();
        out[i] = count as f64;
    })?;
    if opts.normalized {
        let denom = n as f64 - 1.0;
        for x in &mut out {
            *x = if n > 1 { *x / denom } else { f64::NAN };
        }
    }
    Ok(out)
}

// ---- Edge betweenness -----------------------------------------------------------

/// Options for [`edge_betweenness`].
#[derive(Debug, Clone, Copy, Default)]
pub struct EdgeBetweennessOptions {
    /// Use edge weights as lengths (they must be strictly positive).
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Shortest-path algorithm; `Auto` picks the cheapest exact one.
    pub algorithm: PathAlgorithm,
    /// Divide by the number of (ordered, if directed) pairs of nodes.
    pub normalized: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Shortest-path betweenness of every edge, in input order. See the module
/// documentation.
pub fn edge_betweenness(g: &Graph, opts: EdgeBetweennessOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    let m = g.edge_count();
    let lengths = path_lengths(g, opts.weighted, opts.weights_as, opts.algorithm)?;
    let arc_edge = arc_edges(g);

    let partials = parallel::for_each_item(
        n,
        opts.threads,
        || EdgeWorkspace::new(n, m),
        |ws, s| {
            match &lengths {
                Lengths::Uniform => ws.search_bfs(g, s),
                Lengths::Integer(lens, max) => ws.search_dial(g, lens, *max, s),
                Lengths::Real(lens) => ws.search_dijkstra(g, lens, s),
            }
            ws.accumulate(g, &lengths, &arc_edge);
        },
    );
    let mut eb = vec![0.0; m];
    for ws in partials {
        for (b, x) in eb.iter_mut().zip(ws.eb) {
            *b += x;
        }
    }
    let mut scale = if g.is_directed() { 1.0 } else { 0.5 };
    if opts.normalized {
        let pairs = n as f64 * (n as f64 - 1.0);
        scale = if n < 2 {
            f64::NAN
        } else if g.is_directed() {
            scale / pairs
        } else {
            scale * 2.0 / pairs
        };
    }
    if scale != 1.0 {
        for b in &mut eb {
            *b *= scale;
        }
    }
    Ok(eb)
}

/// The input edge behind every outgoing arc, in the CSR order of
/// `g.out()`. [`Graph::from_edges`] builds rows with two stable counting
/// sorts (by target, then by source), which is a stable sort of the arcs by
/// (source, target) in the order they were generated: each edge in input
/// order, and for undirected non-loop edges its reverse right after it.
fn arc_edges(g: &Graph) -> Vec<u32> {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut arcs: Vec<(u32, u32, u32)> = Vec::with_capacity(g.out().arc_count());
    for e in 0..g.edge_count() {
        arcs.push((src[e], tgt[e], e as u32));
        if !g.is_directed() && src[e] != tgt[e] {
            arcs.push((tgt[e], src[e], e as u32));
        }
    }
    arcs.sort_by_key(|a| (a.0, a.1));
    debug_assert!((0..g.node_count()).all(|v| {
        let range = g.out().arc_range(v);
        arcs[range.clone()]
            .iter()
            .zip(g.out().neighbors(v))
            .all(|(a, &t)| a.0 as usize == v && a.1 == t)
    }));
    arcs.into_iter().map(|a| a.2).collect()
}

#[inline]
fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPSILON * a.abs().max(b.abs())
}

/// Per-thread scratch space for edge betweenness.
struct EdgeWorkspace {
    eb: Vec<f64>,
    sigma: Vec<f64>,
    delta: Vec<f64>,
    hops: Vec<u32>,
    idist: Vec<u64>,
    dist: Vec<f64>,
    rank: Vec<u32>,
    order: Vec<u32>,
    buckets: Vec<Vec<u32>>,
    heap: BinaryHeap<HeapItem>,
}

impl EdgeWorkspace {
    fn new(n: usize, m: usize) -> EdgeWorkspace {
        EdgeWorkspace {
            eb: vec![0.0; m],
            sigma: vec![0.0; n],
            delta: vec![0.0; n],
            hops: vec![u32::MAX; n],
            idist: vec![u64::MAX; n],
            dist: vec![f64::INFINITY; n],
            rank: vec![u32::MAX; n],
            order: Vec::with_capacity(n),
            buckets: Vec::new(),
            heap: BinaryHeap::new(),
        }
    }

    fn search_bfs(&mut self, g: &Graph, s: usize) {
        let adj = g.out();
        self.order.clear();
        self.hops[s] = 0;
        self.sigma[s] = 1.0;
        self.rank[s] = 0;
        self.order.push(s as u32);
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let next = self.hops[v] + 1;
            for &w in adj.neighbors(v) {
                let w = w as usize;
                if self.hops[w] == u32::MAX {
                    self.hops[w] = next;
                    self.rank[w] = self.order.len() as u32;
                    self.order.push(w as u32);
                }
                if self.hops[w] == next {
                    self.sigma[w] += self.sigma[v];
                }
            }
        }
    }

    fn search_dial(&mut self, g: &Graph, lens: &[u32], max: u32, s: usize) {
        let adj = g.out();
        let width = max as usize + 1;
        if self.buckets.len() != width {
            self.buckets = vec![Vec::new(); width];
        }
        self.order.clear();
        self.idist[s] = 0;
        self.sigma[s] = 1.0;
        self.buckets[0].push(s as u32);
        let mut pending = 1usize;
        let mut cur: u64 = 0;
        while pending > 0 {
            let slot = (cur % width as u64) as usize;
            let v = match self.buckets[slot].pop() {
                Some(v) => v as usize,
                None => {
                    cur += 1;
                    continue;
                }
            };
            pending -= 1;
            if self.rank[v] != u32::MAX || self.idist[v] != cur {
                continue;
            }
            self.rank[v] = self.order.len() as u32;
            self.order.push(v as u32);
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                if self.rank[w] != u32::MAX {
                    continue;
                }
                let alt = cur + lens[k] as u64;
                let d = self.idist[w];
                if d == u64::MAX || alt < d {
                    self.idist[w] = alt;
                    self.sigma[w] = self.sigma[v];
                    self.buckets[(alt % width as u64) as usize].push(w as u32);
                    pending += 1;
                } else if alt == d {
                    self.sigma[w] += self.sigma[v];
                }
            }
        }
    }

    fn search_dijkstra(&mut self, g: &Graph, lens: &[f64], s: usize) {
        let adj = g.out();
        self.order.clear();
        self.heap.clear();
        self.dist[s] = 0.0;
        self.sigma[s] = 1.0;
        self.heap.push(HeapItem {
            dist: 0.0,
            node: s as u32,
        });
        while let Some(HeapItem { dist: d, node }) = self.heap.pop() {
            let v = node as usize;
            if self.rank[v] != u32::MAX || d > self.dist[v] {
                continue;
            }
            self.rank[v] = self.order.len() as u32;
            self.order.push(v as u32);
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                if self.rank[w] != u32::MAX {
                    continue;
                }
                let alt = d + lens[k];
                let cur = self.dist[w];
                if cur.is_infinite() {
                    self.dist[w] = alt;
                    self.sigma[w] = self.sigma[v];
                    self.heap.push(HeapItem {
                        dist: alt,
                        node: w as u32,
                    });
                } else if approx_eq(alt, cur) {
                    self.sigma[w] += self.sigma[v];
                } else if alt < cur {
                    self.dist[w] = alt;
                    self.sigma[w] = self.sigma[v];
                    self.heap.push(HeapItem {
                        dist: alt,
                        node: w as u32,
                    });
                }
            }
        }
    }

    /// Brandes' accumulation over the shortest-path arcs, farthest nodes
    /// first, crediting each arc's edge; then clears the per-source state.
    fn accumulate(&mut self, g: &Graph, lengths: &Lengths, arc_edge: &[u32]) {
        let adj = g.out();
        for i in (0..self.order.len()).rev() {
            let v = self.order[i] as usize;
            let rank_v = self.rank[v];
            let mut dv = 0.0;
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                let on_path = match lengths {
                    Lengths::Uniform => {
                        self.hops[w] != u32::MAX && self.hops[w] == self.hops[v] + 1
                    }
                    Lengths::Integer(lens, _) => {
                        self.rank[w] != u32::MAX
                            && self.rank[w] > rank_v
                            && self.idist[v] + lens[k] as u64 == self.idist[w]
                    }
                    Lengths::Real(lens) => {
                        self.rank[w] != u32::MAX
                            && self.rank[w] > rank_v
                            && approx_eq(self.dist[v] + lens[k], self.dist[w])
                    }
                };
                if on_path {
                    let c = self.sigma[v] / self.sigma[w] * (1.0 + self.delta[w]);
                    self.eb[arc_edge[k] as usize] += c;
                    dv += c;
                }
            }
            self.delta[v] = dv;
        }
        for &v in &self.order {
            let v = v as usize;
            self.sigma[v] = 0.0;
            self.delta[v] = 0.0;
            self.hops[v] = u32::MAX;
            self.idist[v] = u64::MAX;
            self.dist[v] = f64::INFINITY;
            self.rank[v] = u32::MAX;
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct HeapItem {
    dist: f64,
    node: u32,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for HeapItem {}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapItem {
    // Reversed so that BinaryHeap pops the smallest distance first.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .dist
            .total_cmp(&self.dist)
            .then_with(|| other.node.cmp(&self.node))
    }
}

// ---- Information centrality ----------------------------------------------------

/// Stephenson and Zelen's information centrality of every node, in node
/// order. With `weighted`, edge weights are tie strengths. See the module
/// documentation for conventions.
pub fn information_centrality(g: &Graph, weighted: bool) -> Result<Vec<f64>> {
    let n = g.node_count();
    let w = if weighted {
        let w = g.edge_weights().ok_or(Error::MissingWeights)?;
        if let Some(edge) = w.iter().position(|&x| x <= 0.0) {
            return Err(Error::InvalidArgument(format!(
                "edge {edge} has weight {}; information centrality reads weights as \
                 tie strengths and requires them to be strictly positive",
                w[edge]
            )));
        }
        Some(w)
    } else {
        None
    };

    // Symmetric tie values among distinct nodes, keyed by ordered pair.
    let mut directed_total: std::collections::HashMap<(u32, u32), f64> =
        std::collections::HashMap::new();
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        if s == t {
            continue;
        }
        let key = if g.is_directed() || s < t {
            (s, t)
        } else {
            (t, s)
        };
        let x = w.map_or(1.0, |w| w[e]);
        let entry = directed_total.entry(key).or_insert(0.0);
        *entry = if weighted { *entry + x } else { 1.0 };
    }
    let mut ties: std::collections::BTreeMap<(u32, u32), f64> = std::collections::BTreeMap::new();
    for (&(s, t), &x) in &directed_total {
        let key = if s < t { (s, t) } else { (t, s) };
        let entry = ties.entry(key).or_insert(0.0);
        *entry = entry.max(x);
    }

    // Nodes with at least one tie, in node order.
    let mut index = vec![usize::MAX; n];
    for &(s, t) in ties.keys() {
        index[s as usize] = 0;
        index[t as usize] = 0;
    }
    let mut kept = Vec::new();
    for (v, slot) in index.iter_mut().enumerate() {
        if *slot == 0 {
            *slot = kept.len();
            kept.push(v);
        }
    }
    let k = kept.len();
    let mut out = vec![0.0; n];
    if k == 0 {
        return Ok(out);
    }

    // Connectivity of the kept nodes, by union-find.
    let mut parent: Vec<usize> = (0..k).collect();
    fn find(p: &mut [usize], mut v: usize) -> usize {
        while p[v] != v {
            p[v] = p[p[v]];
            v = p[v];
        }
        v
    }
    for &(s, t) in ties.keys() {
        let (a, b) = (
            find(&mut parent, index[s as usize]),
            find(&mut parent, index[t as usize]),
        );
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let connected = (0..k).all(|v| find(&mut parent, v) == find(&mut parent, 0));
    if !connected {
        for &v in &kept {
            out[v] = f64::NAN;
        }
        return Ok(out);
    }

    // B = L + J over the kept nodes.
    let mut b = vec![1.0; k * k];
    for (&(s, t), &x) in &ties {
        let (i, j) = (index[s as usize], index[t as usize]);
        b[i * k + j] -= x;
        b[j * k + i] -= x;
        b[i * k + i] += x;
        b[j * k + j] += x;
    }
    let c = match spd_inverse(b, k) {
        Some(c) => c,
        None => {
            for &v in &kept {
                out[v] = f64::NAN;
            }
            return Ok(out);
        }
    };
    let trace: f64 = (0..k).map(|i| c[i * k + i]).sum();
    // Every row of C sums to 1/k in exact arithmetic; the mean of the
    // computed row sums is the most accurate estimate.
    let row: f64 = c.iter().sum::<f64>() / k as f64;
    for (i, &v) in kept.iter().enumerate() {
        out[v] = 1.0 / (c[i * k + i] + (trace - 2.0 * row) / k as f64);
    }
    Ok(out)
}

/// Inverse of a symmetric positive definite matrix (row-major, k by k) by
/// Cholesky factorization (Golub and Van Loan 2013, `golub2013matrix`,
/// algorithm 4.2.1), or `None` if the matrix is not numerically positive
/// definite.
fn spd_inverse(mut a: Vec<f64>, k: usize) -> Option<Vec<f64>> {
    // Factor A = L L^T in place (lower triangle).
    for j in 0..k {
        let mut d = a[j * k + j];
        for p in 0..j {
            d -= a[j * k + p] * a[j * k + p];
        }
        if d.is_nan() || d <= 0.0 {
            return None;
        }
        let d = d.sqrt();
        a[j * k + j] = d;
        for i in j + 1..k {
            let mut s = a[i * k + j];
            for p in 0..j {
                s -= a[i * k + p] * a[j * k + p];
            }
            a[i * k + j] = s / d;
        }
    }
    // Solve L L^T x = e_c for every column c.
    let mut inv = vec![0.0; k * k];
    let mut y = vec![0.0; k];
    for col in 0..k {
        for i in 0..k {
            let mut s = if i == col { 1.0 } else { 0.0 };
            for p in 0..i {
                s -= a[i * k + p] * y[p];
            }
            y[i] = s / a[i * k + i];
        }
        for i in (0..k).rev() {
            let mut s = y[i];
            for p in i + 1..k {
                s -= a[p * k + i] * inv[p * k + col];
            }
            inv[i * k + col] = s / a[i * k + i];
        }
    }
    Some(inv)
}

// ---- Bridging -----------------------------------------------------------------

/// Hwang et al.'s bridging coefficient of every node, in node order: the
/// reciprocal of the node's degree divided by the sum of the reciprocals
/// of its neighbors' degrees, counting distinct neighbors and ignoring
/// direction. NaN for isolated nodes.
pub fn bridging_coefficient(g: &Graph) -> Vec<f64> {
    let n = g.node_count();
    let neighbors = distinct_neighbors(g);
    let deg: Vec<f64> = neighbors.iter().map(|nb| nb.len() as f64).collect();
    (0..n)
        .map(|v| {
            if neighbors[v].is_empty() {
                return f64::NAN;
            }
            let s: f64 = neighbors[v].iter().map(|&u| 1.0 / deg[u as usize]).sum();
            (1.0 / deg[v]) / s
        })
        .collect()
}

/// Distinct neighbors of every node other than itself, direction ignored.
fn distinct_neighbors(g: &Graph) -> Vec<Vec<u32>> {
    let s = g.simple_undirected(false);
    (0..g.node_count())
        .map(|v| s.neighbors(v).to_vec())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::betweenness::{betweenness, BetweennessOptions};

    fn graph(n: usize, e: &[(i64, i64)], w: Option<&[f64]>, directed: bool) -> Graph {
        let s: Vec<i64> = e.iter().map(|x| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        Graph::from_edges(n, &s, &t, w, directed).unwrap()
    }

    fn close(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len(), "{a:?} vs {b:?}");
        for (x, y) in a.iter().zip(b) {
            assert!(
                (x.is_nan() && y.is_nan()) || (x - y).abs() <= 1e-10 * (1.0 + y.abs()),
                "{a:?} vs {b:?}"
            );
        }
    }

    #[test]
    fn reach_counts_nodes_within_k() {
        // Path 0 - 1 - 2 - 3 - 4.
        let g = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4)], None, false);
        let r = |k: f64, normalized| {
            reach(
                &g,
                ReachOptions {
                    k,
                    normalized,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        close(&r(1.0, false), &[1.0, 2.0, 2.0, 2.0, 1.0]);
        close(&r(2.0, false), &[2.0, 3.0, 4.0, 3.0, 2.0]);
        close(&r(f64::INFINITY, true), &[1.0; 5]);
        close(&r(0.0, false), &[0.0; 5]);
        assert!(reach(
            &g,
            ReachOptions {
                k: -1.0,
                ..Default::default()
            }
        )
        .is_err());
        // Directed: 0 -> 1 -> 2; in-reach counts nodes that reach the node.
        let d = graph(3, &[(0, 1), (1, 2)], None, true);
        let out = reach(&d, ReachOptions::default()).unwrap();
        close(&out, &[2.0, 1.0, 0.0]);
        let inn = reach(
            &d,
            ReachOptions {
                mode: DegreeMode::In,
                ..Default::default()
            },
        )
        .unwrap();
        close(&inn, &[0.0, 1.0, 2.0]);
        // Weighted: k is a length.
        let w = graph(3, &[(0, 1), (1, 2)], Some(&[0.5, 2.0]), false);
        let opts = ReachOptions {
            k: 1.0,
            weighted: true,
            ..Default::default()
        };
        close(&reach(&w, opts).unwrap(), &[1.0, 1.0, 0.0]);
    }

    #[test]
    fn edge_betweenness_by_hand() {
        // Path 0 - 1 - 2: each edge carries the pair it joins and (0, 2).
        let g = graph(3, &[(0, 1), (1, 2)], None, false);
        let eb = edge_betweenness(&g, Default::default()).unwrap();
        close(&eb, &[2.0, 2.0]);
        // Directed: 0 -> 1 -> 2 carries (0,1), (0,2) and (1,2), (0,2).
        let d = graph(3, &[(0, 1), (1, 2)], None, true);
        close(
            &edge_betweenness(&d, Default::default()).unwrap(),
            &[2.0, 2.0],
        );
        // Parallel edges split the paths; a loop carries nothing.
        let p = graph(3, &[(0, 1), (0, 1), (1, 2), (2, 2)], None, false);
        close(
            &edge_betweenness(&p, Default::default()).unwrap(),
            &[1.0, 1.0, 2.0, 0.0],
        );
        // Square 0-1-2-3-0: each edge carries its own pair and half of two
        // diagonals.
        let c = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0)], None, false);
        close(
            &edge_betweenness(&c, Default::default()).unwrap(),
            &[2.0; 4],
        );
        let norm = edge_betweenness(
            &c,
            EdgeBetweennessOptions {
                normalized: true,
                ..Default::default()
            },
        )
        .unwrap();
        close(&norm, &[2.0 / 6.0; 4]);
    }

    #[test]
    fn edge_betweenness_sums_to_total_path_length() {
        // On a connected unweighted graph the edge betweenness values sum
        // to the sum of all distances (each shortest path of length d
        // crosses d edges).
        let e = [
            (0, 1),
            (0, 2),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (3, 5),
            (1, 5),
        ];
        let g = graph(6, &e, None, false);
        let eb = edge_betweenness(&g, Default::default()).unwrap();
        let mut total = 0.0;
        distance_rows(&g, None, PathOptions::default(), |i, row| {
            for (j, &d) in row.iter().enumerate() {
                if j > i {
                    total += d;
                }
            }
        })
        .unwrap();
        assert!((eb.iter().sum::<f64>() - total).abs() < 1e-12);
    }

    #[test]
    fn weighted_edge_betweenness_agrees_across_algorithms() {
        let e = [(0, 1), (0, 2), (1, 3), (2, 3), (1, 4), (3, 5), (4, 5)];
        let w = [2.0, 1.0, 1.0, 2.0, 3.0, 2.0, 1.0];
        let g = graph(6, &e, Some(&w), true);
        let run = |algorithm| {
            edge_betweenness(
                &g,
                EdgeBetweennessOptions {
                    weighted: true,
                    algorithm,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let dial = run(PathAlgorithm::Dial);
        close(&dial, &run(PathAlgorithm::Dijkstra));
        // Node betweenness is half the sum over a node's arcs of
        // (in-flow + out-flow) minus the pairs it starts or ends.
        let nb = betweenness(
            &g,
            BetweennessOptions {
                weighted: true,
                ..Default::default()
            },
        )
        .unwrap();
        // Flow conservation: for every node v, sum of incoming arc
        // betweenness = b(v) + (number of nodes reaching v).
        let mut inflow = [0.0; 6];
        for (k, &(_, t)) in e.iter().enumerate() {
            inflow[t as usize] += dial[k];
        }
        let reached = reach(
            &g,
            ReachOptions {
                k: f64::INFINITY,
                mode: DegreeMode::In,
                ..Default::default()
            },
        )
        .unwrap();
        for v in 0..6 {
            assert!((inflow[v] - nb[v] - reached[v]).abs() < 1e-12, "{v}");
        }
    }

    #[test]
    fn information_centrality_matches_the_formula() {
        // Star with centre 0 and 3 leaves. B = L + J; computed by hand:
        // the centre has I = 4/3 * ... checked against the closed form
        // I_i = n / sum_j R_ij with effective resistances R (1 from the
        // centre to each leaf, 2 between leaves).
        let g = graph(4, &[(0, 1), (0, 2), (0, 3)], None, false);
        let ic = information_centrality(&g, false).unwrap();
        close(&ic, &[4.0 / 3.0, 4.0 / 5.0, 4.0 / 5.0, 4.0 / 5.0]);
        // Isolates score 0 and are removed before the formula is applied.
        let g = graph(5, &[(0, 1), (0, 2), (0, 3)], None, false);
        close(
            &information_centrality(&g, false).unwrap(),
            &[4.0 / 3.0, 0.8, 0.8, 0.8, 0.0],
        );
        // Two components with ties: undefined.
        let g = graph(4, &[(0, 1), (2, 3)], None, false);
        assert!(information_centrality(&g, false)
            .unwrap()
            .iter()
            .all(|x| x.is_nan()));
        // Weights are conductances: doubling every weight halves every
        // resistance and doubles every score.
        let e = [(0, 1), (1, 2), (0, 2), (2, 3)];
        let w1 = graph(4, &e, Some(&[1.0, 2.0, 3.0, 1.5]), false);
        let w2 = graph(4, &e, Some(&[2.0, 4.0, 6.0, 3.0]), false);
        let a = information_centrality(&w1, true).unwrap();
        let b = information_centrality(&w2, true).unwrap();
        close(&b, &a.iter().map(|x| 2.0 * x).collect::<Vec<_>>());
        // A directed network is symmetrized, parallel edges collapse.
        let d = graph(4, &[(0, 1), (1, 0), (0, 2), (3, 0), (3, 0)], None, true);
        close(
            &information_centrality(&d, false).unwrap(),
            &[4.0 / 3.0, 0.8, 0.8, 0.8],
        );
    }

    #[test]
    fn bridging_coefficient_by_hand() {
        // Path 0 - 1 - 2: the middle node has degree 2 and neighbors of
        // degree 1: (1/2) / 2 = 1/4; ends: 1 / (1/2) = 2.
        let g = graph(4, &[(0, 1), (1, 2), (1, 2)], None, false);
        let bc = bridging_coefficient(&g);
        close(&bc[..3], &[2.0, 0.25, 2.0]);
        assert!(bc[3].is_nan());
    }
}
