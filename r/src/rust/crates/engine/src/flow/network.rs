//! A general residual network and maximum flow by Dinic's algorithm.
//!
//! Sources:
//! - Maximum flow and the max-flow min-cut theorem: Ford and Fulkerson
//!   (1956, `ford1956maximal`).
//! - Dinic's (1970, `dinic1970algorithm`) blocking-flow algorithm: each
//!   phase builds breadth-first levels from the source over arcs with
//!   residual capacity and saturates every shortest augmenting path; the
//!   source-sink distance grows in every phase, so there are at most n - 1
//!   phases and O(n^2 m) work in all, O(m min(n^(2/3), m^(1/2))) with unit
//!   capacities (Even and Tarjan 1975, `even1975network`).
//!
//! Representation: arcs come in pairs, arc `2k` and arc `2k + 1` being each
//! other's reverse. A directed arc u -> v of capacity c is the pair
//! (u -> v, c) and (v -> u, 0); an undirected edge is (u -> v, c) and
//! (v -> u, c). The flow on an arc is its capacity minus its residual
//! capacity, so the flow on arc `2k + 1` is always minus the flow on arc
//! `2k`, and for an undirected edge the flow on arc `2k` is the net flow
//! from u to v, between -c and c.
//!
//! Capacities may be `+inf`. A source-sink path made only of infinite arcs
//! means the maximum flow is unbounded, which is reported as an error
//! before any flow is pushed; otherwise every augmenting path contains a
//! finite arc, so flows stay finite.
//!
//! Floating point: residual capacities at most `eps` (1e-12 times the
//! largest finite capacity, or 1e-12 when that is below 1) count as zero,
//! so rounding in sums of real-valued capacities cannot leave arcs that
//! look unsaturated. With integer capacities below 2^53 every quantity is
//! an exact integer and `eps` never matters.
//!
//! Every choice (levels by breadth-first search over arcs in the order
//! they were added, depth-first augmentation in the same order) is
//! deterministic, so a network gives the same flow every time.

use crate::error::{Error, Result};

/// A residual network on nodes `0..n` (see the module documentation).
#[derive(Debug, Clone)]
pub struct FlowNetwork {
    n: usize,
    head: Vec<u32>,
    cap: Vec<f64>,
    /// Flow on the forward arc of each pair (the reverse carries minus it).
    pf: Vec<f64>,
    /// Arcs leaving each node, in order of addition (built by `prepare`).
    offsets: Vec<usize>,
    arcs: Vec<u32>,
    prepared: bool,
    eps: f64,
    // Scratch.
    level: Vec<u32>,
    next: Vec<usize>,
    queue: Vec<u32>,
}

impl FlowNetwork {
    /// An empty network on `n` nodes.
    pub fn new(n: usize) -> FlowNetwork {
        FlowNetwork {
            n,
            head: Vec::new(),
            cap: Vec::new(),
            pf: Vec::new(),
            offsets: Vec::new(),
            arcs: Vec::new(),
            prepared: false,
            eps: 1e-12,
            level: Vec::new(),
            next: Vec::new(),
            queue: Vec::new(),
        }
    }

    /// Adds a node and returns its index.
    pub fn add_node(&mut self) -> usize {
        self.n += 1;
        self.prepared = false;
        self.n - 1
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.n
    }

    /// Number of arc pairs added.
    pub fn pair_count(&self) -> usize {
        self.head.len() / 2
    }

    fn push_pair(&mut self, u: usize, v: usize, forward: f64, backward: f64) -> usize {
        debug_assert!(u < self.n && v < self.n);
        debug_assert!(forward >= 0.0 && backward >= 0.0);
        let k = self.head.len() / 2;
        self.head.push(v as u32);
        self.cap.push(forward);
        self.head.push(u as u32);
        self.cap.push(backward);
        self.prepared = false;
        k
    }

    /// Adds an arc u -> v with capacity `cap` (non-negative, possibly
    /// infinite) and returns its pair index.
    pub fn add_arc(&mut self, u: usize, v: usize, cap: f64) -> usize {
        self.push_pair(u, v, cap, 0.0)
    }

    /// Adds an undirected edge {u, v} with capacity `cap` each way and
    /// returns its pair index.
    pub fn add_edge(&mut self, u: usize, v: usize, cap: f64) -> usize {
        self.push_pair(u, v, cap, cap)
    }

    /// Tail of arc `a` (the head of its reverse).
    #[inline]
    fn tail(&self, a: usize) -> usize {
        self.head[a ^ 1] as usize
    }

    fn prepare(&mut self) {
        if self.prepared {
            return;
        }
        let n = self.n;
        let mut offsets = vec![0usize; n + 1];
        for a in 0..self.head.len() {
            offsets[self.tail(a) + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let mut cursor = offsets.clone();
        let mut arcs = vec![0u32; self.head.len()];
        for a in 0..self.head.len() {
            let t = self.tail(a);
            arcs[cursor[t]] = a as u32;
            cursor[t] += 1;
        }
        self.offsets = offsets;
        self.arcs = arcs;
        let max_finite = self
            .cap
            .iter()
            .copied()
            .filter(|c| c.is_finite())
            .fold(0.0f64, f64::max);
        self.eps = 1e-12 * max_finite.max(1.0);
        self.pf = vec![0.0; self.head.len() / 2];
        self.level = vec![u32::MAX; n];
        self.next = vec![0; n];
        self.queue = Vec::with_capacity(n);
        self.prepared = true;
    }

    /// Forgets any flow, so the network can be solved again.
    pub fn reset(&mut self) {
        if self.prepared {
            for f in self.pf.iter_mut() {
                *f = 0.0;
            }
        }
    }

    /// Sets the capacity of the forward arc of pair `k` (the reverse keeps
    /// its own, 0 for an arc, the same for an edge if `both`) and forgets
    /// any flow.
    pub fn set_capacity(&mut self, k: usize, cap: f64, both: bool) {
        self.cap[2 * k] = cap;
        if both {
            self.cap[2 * k + 1] = cap;
        }
        self.prepared = false;
    }

    /// Capacity of the forward arc of pair `k`.
    pub fn capacity(&self, k: usize) -> f64 {
        self.cap[2 * k]
    }

    /// Flow on the forward arc of pair `k` after a solve (for an
    /// undirected edge, the net flow from its first to its second node).
    pub fn flow(&self, k: usize) -> f64 {
        if self.prepared {
            self.pf[k]
        } else {
            0.0
        }
    }

    /// Residual capacity of arc `a`.
    #[inline]
    fn residual(&self, a: usize) -> f64 {
        if a & 1 == 0 {
            self.cap[a] - self.pf[a >> 1]
        } else {
            self.cap[a] + self.pf[a >> 1]
        }
    }

    /// Pushes `b` units along arc `a`; with `saturate`, sets its residual
    /// to exactly zero instead.
    #[inline]
    fn push(&mut self, a: usize, b: f64, saturate: bool) {
        let k = a >> 1;
        self.pf[k] = match (a & 1 == 0, saturate) {
            (true, true) => self.cap[a],
            (true, false) => self.pf[k] + b,
            (false, true) => -self.cap[a],
            (false, false) => self.pf[k] - b,
        };
    }

    /// Whether arc `a` has residual capacity.
    #[inline]
    fn open(&self, a: usize) -> bool {
        self.residual(a) > self.eps
    }

    /// Checks that sources and sinks are valid and disjoint, and that no
    /// path of infinite arcs joins them.
    fn check_terminals(&mut self, sources: &[usize], sinks: &[usize]) -> Result<()> {
        if sources.is_empty() || sinks.is_empty() {
            return Err(Error::InvalidArgument(
                "a flow needs at least one source and one sink".into(),
            ));
        }
        let mut role = vec![0u8; self.n];
        for &s in sources {
            if s >= self.n {
                return Err(Error::InvalidArgument(format!(
                    "source {s} is not a node (n = {})",
                    self.n
                )));
            }
            role[s] = 1;
        }
        for &t in sinks {
            if t >= self.n {
                return Err(Error::InvalidArgument(format!(
                    "sink {t} is not a node (n = {})",
                    self.n
                )));
            }
            if role[t] == 1 {
                return Err(Error::InvalidArgument(format!(
                    "node {t} is both a source and a sink"
                )));
            }
            role[t] = 2;
        }
        // Breadth-first search over infinite arcs.
        let mut seen = vec![false; self.n];
        let mut queue: Vec<usize> = sources.to_vec();
        for &s in sources {
            seen[s] = true;
        }
        let mut i = 0;
        while i < queue.len() {
            let v = queue[i];
            i += 1;
            for &a in &self.arcs[self.offsets[v]..self.offsets[v + 1]] {
                let a = a as usize;
                if self.cap[a] == f64::INFINITY {
                    let w = self.head[a] as usize;
                    if role[w] == 2 {
                        return Err(Error::InvalidArgument(
                            "a path of infinite capacity joins a source and a sink, so the \
                             maximum flow is unbounded"
                                .into(),
                        ));
                    }
                    if !seen[w] {
                        seen[w] = true;
                        queue.push(w);
                    }
                }
            }
        }
        Ok(())
    }

    /// Maximum flow from `s` to `t`. With `limit`, stops once the flow
    /// reaches it (the result is then at least `limit` but need not be
    /// maximum). Flow already in the network is kept and added to, so
    /// call [`FlowNetwork::reset`] first to start from zero.
    pub fn max_flow(&mut self, s: usize, t: usize, limit: Option<f64>) -> Result<f64> {
        self.prepare();
        self.check_terminals(&[s], &[t])?;
        let limit = limit.unwrap_or(f64::INFINITY);
        let mut total = 0.0;
        while total < limit && self.levels(s, t) {
            total += self.blocking(s, t, limit - total);
        }
        Ok(total)
    }

    /// Breadth-first levels from `s` over open arcs; whether `t` is
    /// reached.
    fn levels(&mut self, s: usize, t: usize) -> bool {
        for l in self.level.iter_mut() {
            *l = u32::MAX;
        }
        self.queue.clear();
        self.level[s] = 0;
        self.queue.push(s as u32);
        let mut i = 0;
        while i < self.queue.len() {
            let v = self.queue[i] as usize;
            i += 1;
            if v == t {
                // Nodes at t's level or beyond cannot lie on a shortest
                // path to t.
                break;
            }
            for &a in &self.arcs[self.offsets[v]..self.offsets[v + 1]] {
                let a = a as usize;
                let w = self.head[a] as usize;
                if self.level[w] == u32::MAX && self.open(a) {
                    self.level[w] = self.level[v] + 1;
                    self.queue.push(w as u32);
                }
            }
        }
        self.level[t] != u32::MAX
    }

    /// One blocking flow along level-increasing arcs, by depth-first
    /// search with a current-arc pointer per node; at most `want` units.
    fn blocking(&mut self, s: usize, t: usize, want: f64) -> f64 {
        for v in 0..self.n {
            self.next[v] = self.offsets[v];
        }
        let mut pushed = 0.0;
        let mut path: Vec<u32> = Vec::new();
        let mut v = s;
        loop {
            if v == t {
                // Augment along the path by its bottleneck.
                let mut b = want - pushed;
                for &a in &path {
                    b = b.min(self.residual(a as usize));
                }
                let mut cut_at = path.len();
                for i in 0..path.len() {
                    let a = path[i] as usize;
                    // The bottleneck arc becomes exactly zero.
                    let saturate = self.residual(a) == b;
                    self.push(a, b, saturate);
                    if cut_at == path.len() && !self.open(a) {
                        cut_at = i;
                    }
                }
                pushed += b;
                if pushed >= want || cut_at == path.len() {
                    // Either done, or (only when `want` was the bottleneck)
                    // no arc saturated.
                    break;
                }
                // Retreat to the tail of the first saturated arc.
                v = self.tail(path[cut_at] as usize);
                path.truncate(cut_at);
                continue;
            }
            let mut advanced = false;
            while self.next[v] < self.offsets[v + 1] {
                let a = self.arcs[self.next[v]] as usize;
                let w = self.head[a] as usize;
                if self.open(a) && self.level[w] == self.level[v] + 1 {
                    path.push(a as u32);
                    v = w;
                    advanced = true;
                    break;
                }
                self.next[v] += 1;
            }
            if advanced {
                continue;
            }
            // Dead end: nothing more passes through v in this phase.
            self.level[v] = u32::MAX;
            match path.pop() {
                None => break,
                Some(a) => {
                    v = self.tail(a as usize);
                    self.next[v] += 1;
                }
            }
        }
        pushed
    }

    /// After a solve: the nodes reachable from `s` over open arcs (the
    /// smallest source side of a minimum cut).
    pub fn reachable_from(&self, sources: &[usize]) -> Vec<bool> {
        let mut seen = vec![false; self.n];
        if !self.prepared {
            return seen;
        }
        let mut queue: Vec<usize> = Vec::new();
        for &s in sources {
            if !seen[s] {
                seen[s] = true;
                queue.push(s);
            }
        }
        let mut i = 0;
        while i < queue.len() {
            let v = queue[i];
            i += 1;
            for &a in &self.arcs[self.offsets[v]..self.offsets[v + 1]] {
                let a = a as usize;
                let w = self.head[a] as usize;
                if !seen[w] && self.open(a) {
                    seen[w] = true;
                    queue.push(w);
                }
            }
        }
        seen
    }

    /// After a solve: the nodes from which `t` is reachable over open arcs
    /// (the smallest sink side of a minimum cut).
    pub fn reaching(&self, sinks: &[usize]) -> Vec<bool> {
        let mut seen = vec![false; self.n];
        if !self.prepared {
            return seen;
        }
        let mut queue: Vec<usize> = Vec::new();
        for &t in sinks {
            if !seen[t] {
                seen[t] = true;
                queue.push(t);
            }
        }
        let mut i = 0;
        while i < queue.len() {
            let v = queue[i];
            i += 1;
            // Arc b leaves v; its reverse b ^ 1 enters v from head[b].
            for &b in &self.arcs[self.offsets[v]..self.offsets[v + 1]] {
                let b = b as usize;
                let u = self.head[b] as usize;
                if !seen[u] && self.open(b ^ 1) {
                    seen[u] = true;
                    queue.push(u);
                }
            }
        }
        seen
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn textbook_example() {
        // CLRS figure 26.1: maximum flow 23.
        let mut f = FlowNetwork::new(6);
        for &(u, v, c) in &[
            (0, 1, 16.0),
            (0, 2, 13.0),
            (2, 1, 4.0),
            (1, 3, 12.0),
            (3, 2, 9.0),
            (2, 4, 14.0),
            (4, 3, 7.0),
            (3, 5, 20.0),
            (4, 5, 4.0),
        ] {
            f.add_arc(u, v, c);
        }
        assert_eq!(f.max_flow(0, 5, None).unwrap(), 23.0);
        let side = f.reachable_from(&[0]);
        let cut: f64 = (0..f.pair_count())
            .filter(|&k| side[f.tail(2 * k)] && !side[f.head[2 * k] as usize])
            .map(|k| f.capacity(k))
            .sum();
        assert_eq!(cut, 23.0);
        f.reset();
        assert!(f.max_flow(0, 5, Some(5.0)).unwrap() >= 5.0);
    }

    #[test]
    fn undirected_and_infinite() {
        let mut f = FlowNetwork::new(4);
        f.add_edge(0, 1, f64::INFINITY);
        f.add_edge(1, 2, 2.5);
        f.add_edge(2, 3, f64::INFINITY);
        f.add_edge(0, 3, 1.0);
        assert_eq!(f.max_flow(0, 3, None).unwrap(), 3.5);
        f.reset();
        assert_eq!(f.max_flow(3, 0, None).unwrap(), 3.5);
        assert!((f.flow(1) + 2.5).abs() < 1e-15);
        let mut g = FlowNetwork::new(2);
        g.add_edge(0, 1, f64::INFINITY);
        assert!(g.max_flow(0, 1, None).is_err());
        assert!(g.max_flow(0, 0, None).is_err());
    }
}
