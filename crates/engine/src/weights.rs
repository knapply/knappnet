//! How edge weights are interpreted, and which shortest-path algorithm a
//! set of weights allows.
//!
//! A weight can mean two opposite things. In transportation, computer, and
//! many infrastructure networks it is a *distance* or cost: larger values
//! keep nodes apart. In most social networks it is a *strength* (contact
//! frequency, tie intensity): larger values bring nodes closer. Path-based
//! measures need distances, so strengths are converted with the reciprocal
//! 1/w, following Newman (2001, `newman2001scientific`) and Opsahl,
//! Agneessens, and Skvoretz (2010, `opsahl2010node`). Measures that sum
//! weights, such as strength, use them as given.
//!
//! Once lengths are known, the cheapest exact algorithm is chosen:
//!
//! | Lengths | Algorithm | Cost per source |
//! |---|---|---|
//! | unweighted, or all equal | breadth-first search | O(m) |
//! | positive integers up to [`DIAL_MAX`] | Dial's bucket queue | O(m + D) |
//! | anything else | Dijkstra with a 4-ary indexed heap | O(m log n) |
//!
//! where D is the largest finite distance. Uniform lengths do not change
//! which paths are shortest or how many there are, so they reduce to the
//! unweighted case exactly. Integer lengths are summed in exact integer
//! arithmetic, so ties are never lost to rounding (Dial 1969,
//! `dial1969algorithm`).

use crate::error::{Error, Result};
use crate::graph::Graph;

/// What edge weights mean.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WeightKind {
    /// Larger weights keep nodes farther apart (cost, length, latency).
    #[default]
    Distance,
    /// Larger weights bring nodes closer (tie strength, frequency). Path
    /// lengths use 1/w.
    Strength,
}

impl WeightKind {
    /// Decodes the codes used by the C interface: 0 = distance, 1 = strength.
    pub fn from_code(code: i32) -> Result<WeightKind> {
        match code {
            0 => Ok(WeightKind::Distance),
            1 => Ok(WeightKind::Strength),
            _ => Err(Error::InvalidArgument(format!(
                "unknown weight kind code {code}"
            ))),
        }
    }
}

/// Which shortest-path algorithm to use. `Auto` picks the cheapest exact one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PathAlgorithm {
    /// Choose from the weights (recommended).
    #[default]
    Auto,
    /// Breadth-first search; ignores weights.
    Bfs,
    /// Dial's bucket queue; requires positive integer lengths.
    Dial,
    /// Dijkstra's algorithm with a binary heap.
    Dijkstra,
}

/// Largest integer length for which Dial's algorithm is used automatically.
/// Beyond this, scanning empty buckets can cost more than a heap.
pub const DIAL_MAX: u64 = 1024;

/// Arc lengths prepared for a shortest-path search, aligned with the
/// graph's outgoing CSR arcs.
#[derive(Debug, Clone)]
pub enum Lengths {
    /// Every arc has the same length (or the graph is unweighted).
    Uniform,
    /// Positive integer lengths, with their maximum.
    Integer(Vec<u32>, u32),
    /// Positive real lengths.
    Real(Vec<f64>),
}

/// Converts a graph's weights into arc lengths and chooses an algorithm.
///
/// With `weighted == false`, or a graph without weights, the result is
/// [`Lengths::Uniform`]. Weights must be strictly positive whether they are
/// distances or strengths.
pub fn path_lengths(
    g: &Graph,
    weighted: bool,
    kind: WeightKind,
    algorithm: PathAlgorithm,
) -> Result<Lengths> {
    if !weighted || algorithm == PathAlgorithm::Bfs {
        return Ok(Lengths::Uniform);
    }
    let raw = g.edge_weights().ok_or(Error::MissingWeights)?;
    if let Some(edge) = raw.iter().position(|&x| x <= 0.0) {
        return Err(Error::NonPositiveWeight {
            edge,
            weight: raw[edge],
        });
    }
    let arcs = g.out().all_weights().ok_or(Error::MissingWeights)?;
    let lens: Vec<f64> = match kind {
        WeightKind::Distance => arcs.to_vec(),
        WeightKind::Strength => arcs.iter().map(|w| 1.0 / w).collect(),
    };

    let uniform = lens.windows(2).all(|p| p[0] == p[1]);
    let integer_max = lens.iter().try_fold(0u64, |m, &x| {
        (x.fract() == 0.0 && x >= 1.0 && x <= u32::MAX as f64).then(|| m.max(x as u64))
    });

    match algorithm {
        PathAlgorithm::Auto if uniform => Ok(Lengths::Uniform),
        PathAlgorithm::Auto => match integer_max {
            Some(max) if max <= DIAL_MAX => Ok(integer(lens, max)),
            _ => Ok(Lengths::Real(lens)),
        },
        PathAlgorithm::Dial => match integer_max {
            Some(max) => Ok(integer(lens, max)),
            None => Err(Error::InvalidArgument(
                "Dial's algorithm needs positive integer path lengths".into(),
            )),
        },
        PathAlgorithm::Dijkstra => Ok(Lengths::Real(lens)),
        PathAlgorithm::Bfs => unreachable!("handled above"),
    }
}

/// The edge weights of `g` for additive use by `what` (for example
/// "modularity"), in input order, checked to be non-negative. With
/// `weighted == false` the result is `None`; a weighted request on a graph
/// without weights is an error.
pub fn additive_weights<'a>(g: &'a Graph, weighted: bool, what: &str) -> Result<Option<&'a [f64]>> {
    if !weighted {
        return Ok(None);
    }
    let w = g.edge_weights().ok_or(Error::MissingWeights)?;
    if let Some(edge) = w.iter().position(|&x| x < 0.0) {
        return Err(Error::InvalidArgument(format!(
            "edge {edge} has negative weight {}; {what} requires non-negative weights",
            w[edge]
        )));
    }
    Ok(Some(w))
}

fn integer(lens: Vec<f64>, max: u64) -> Lengths {
    Lengths::Integer(lens.into_iter().map(|x| x as u32).collect(), max as u32)
}

/// The reference for [`RadixHeap`] in tests: a 4-ary min-heap of nodes
/// keyed by tentative distance, with a position
/// index so that a node's key can be lowered in place (Johnson 1975,
/// `johnson1975priority`). Dijkstra's algorithm then holds at most one
/// entry per node instead of one per relaxation, and a 4-ary layout halves
/// the depth of a binary heap, so a sift touches fewer cache lines.
///
/// Entries are ordered by (distance, node index), a total order, so nodes
/// are settled in exactly the order a binary heap with the same ordering
/// and lazy deletion settles them.
///
/// Keys are held as the integers whose order is [`f64::total_cmp`]'s (the
/// IEEE 754 total order), so each comparison is one integer comparison
/// instead of a floating-point total comparison; the order, and so every
/// result, is unchanged.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct IndexedHeap {
    items: Vec<(i64, u32)>,
    /// Position of each node in `items`; `u32::MAX` when absent.
    pos: Vec<u32>,
}

/// The integer whose order among keys is `f64::total_cmp`'s: negative
/// values have their magnitude bits flipped, so they order downward.
#[inline]
fn key(x: f64) -> i64 {
    let b = x.to_bits() as i64;
    b ^ (((b >> 63) as u64) >> 1) as i64
}

/// The inverse of [`key`].
#[inline]
fn unkey(k: i64) -> f64 {
    f64::from_bits((k ^ (((k >> 63) as u64) >> 1) as i64) as u64)
}

#[cfg(test)]
impl IndexedHeap {
    const ARITY: usize = 4;

    /// An empty heap for nodes `0..n`.
    pub(crate) fn new(n: usize) -> IndexedHeap {
        IndexedHeap {
            items: Vec::new(),
            pos: vec![u32::MAX; n],
        }
    }

    #[inline]
    fn less(a: (i64, u32), b: (i64, u32)) -> bool {
        a < b
    }

    /// Inserts `node` with key `dist`, or lowers its key to `dist` if it
    /// is present with a larger key. A larger `dist` is ignored.
    #[inline]
    pub(crate) fn push_or_decrease(&mut self, node: u32, dist: f64) {
        let k = key(dist);
        let p = self.pos[node as usize];
        let i = if p == u32::MAX {
            self.items.push((k, node));
            self.items.len() - 1
        } else {
            let i = p as usize;
            if !Self::less((k, node), self.items[i]) {
                return;
            }
            self.items[i].0 = k;
            i
        };
        self.sift_up(i);
    }

    /// Removes and returns the entry with the smallest key.
    #[inline]
    pub(crate) fn pop(&mut self) -> Option<(f64, u32)> {
        let top = *self.items.first()?;
        self.pos[top.1 as usize] = u32::MAX;
        let last = self.items.pop().expect("not empty");
        if !self.items.is_empty() {
            self.items[0] = last;
            self.sift_down(0);
        }
        Some((unkey(top.0), top.1))
    }

    /// Removes every entry.
    pub(crate) fn clear(&mut self) {
        for &(_, v) in &self.items {
            self.pos[v as usize] = u32::MAX;
        }
        self.items.clear();
    }

    /// The nodes in the heap, in no particular order.
    pub(crate) fn nodes(&self) -> impl Iterator<Item = u32> + '_ {
        self.items.iter().map(|&(_, v)| v)
    }

    fn sift_up(&mut self, mut i: usize) {
        let item = self.items[i];
        while i > 0 {
            let parent = (i - 1) / Self::ARITY;
            let p = self.items[parent];
            if !Self::less(item, p) {
                break;
            }
            self.items[i] = p;
            self.pos[p.1 as usize] = i as u32;
            i = parent;
        }
        self.items[i] = item;
        self.pos[item.1 as usize] = i as u32;
    }

    fn sift_down(&mut self, mut i: usize) {
        let len = self.items.len();
        let item = self.items[i];
        loop {
            let first = i * Self::ARITY + 1;
            if first >= len {
                break;
            }
            let mut best = first;
            for c in first + 1..(first + Self::ARITY).min(len) {
                if Self::less(self.items[c], self.items[best]) {
                    best = c;
                }
            }
            let b = self.items[best];
            if !Self::less(b, item) {
                break;
            }
            self.items[i] = b;
            self.pos[b.1 as usize] = i as u32;
            i = best;
        }
        self.items[i] = item;
        self.pos[item.1 as usize] = i as u32;
    }
}

/// A radix heap (Ahuja, Mehlhorn, Orlin, and Tarjan 1990,
/// `ahuja1990faster`) of nodes keyed by tentative distance, with the same
/// interface and the same order of removal as [`IndexedHeap`]: (distance,
/// node) ascending in the IEEE 754 total order.
///
/// Dijkstra's algorithm removes keys in non-decreasing order, which a radix
/// heap exploits: each (distance, node) pair is read as a 128-bit integer
/// (the order-preserving integer of the distance above the node index),
/// and an entry lives in the bucket numbered by the highest bit in which it
/// differs from the last key removed. Insertion is a bit scan and a push;
/// removal empties the lowest non-empty bucket into lower ones, so an
/// entry moves at most once per bit and most moves touch contiguous memory.
/// Lowering a key adds a new entry and leaves the old one to be skipped
/// when it surfaces. A key below the last one removed (possible only when
/// a floating-point sum rounds back to the last distance and the node index
/// is smaller) goes to bucket 0, from which the smallest is always removed
/// first, so the order stays exact.
#[derive(Debug, Clone)]
pub(crate) struct RadixHeap {
    buckets: Vec<Vec<u128>>,
    /// Current key of each node in the heap; `NONE` when absent.
    current: Vec<u64>,
    last: u128,
    len: usize,
}

impl RadixHeap {
    const NONE: u64 = u64::MAX;

    /// An empty heap for nodes `0..n`.
    pub(crate) fn new(n: usize) -> RadixHeap {
        RadixHeap {
            buckets: vec![Vec::new(); 129],
            current: vec![Self::NONE; n],
            last: 0,
            len: 0,
        }
    }

    /// The unsigned integer whose order is the total order of distances.
    #[inline]
    fn ukey(x: f64) -> u64 {
        (key(x) as u64) ^ (1 << 63)
    }

    #[inline]
    fn bucket(&self, e: u128) -> usize {
        if e <= self.last {
            0
        } else {
            128 - (e ^ self.last).leading_zeros() as usize
        }
    }

    /// Inserts `node` with key `dist`, or lowers its key to `dist` if it
    /// is present with a larger key. A larger `dist` is ignored.
    #[inline]
    pub(crate) fn push_or_decrease(&mut self, node: u32, dist: f64) {
        let k = Self::ukey(dist);
        let cur = self.current[node as usize];
        if cur != Self::NONE {
            if k >= cur {
                return;
            }
        } else {
            self.len += 1;
        }
        self.current[node as usize] = k;
        let e = (k as u128) << 32 | node as u128;
        let b = self.bucket(e);
        self.buckets[b].push(e);
    }

    /// Removes and returns the entry with the smallest key.
    #[inline]
    pub(crate) fn pop(&mut self) -> Option<(f64, u32)> {
        loop {
            if self.len == 0 {
                return None;
            }
            if self.buckets[0].is_empty() {
                let i = (1..129)
                    .find(|&i| !self.buckets[i].is_empty())
                    .expect("a live entry");
                let mut items = std::mem::take(&mut self.buckets[i]);
                self.last = *items.iter().min().expect("not empty");
                for &e in &items {
                    let b = self.bucket(e);
                    self.buckets[b].push(e);
                }
                items.clear();
                // Keep the emptied bucket's storage for later.
                self.buckets[i] = items;
            }
            let zero = &mut self.buckets[0];
            let at = (0..zero.len())
                .min_by_key(|&j| zero[j])
                .expect("bucket 0 is not empty");
            let e = zero.swap_remove(at);
            let node = e as u32;
            let k = (e >> 32) as u64;
            if self.current[node as usize] != k {
                continue; // superseded by a lower key
            }
            self.current[node as usize] = Self::NONE;
            self.len -= 1;
            if self.len == 0 {
                // Only superseded entries remain: drop them, so the next
                // search starts from key 0.
                self.buckets.iter_mut().for_each(Vec::clear);
                self.last = 0;
            }
            return Some((unkey((k ^ (1 << 63)) as i64), node));
        }
    }

    /// Removes every entry.
    pub(crate) fn clear(&mut self) {
        for b in &mut self.buckets {
            for &e in b.iter() {
                self.current[e as u32 as usize] = Self::NONE;
            }
            b.clear();
        }
        self.last = 0;
        self.len = 0;
    }

    /// The nodes in the heap, in no particular order.
    pub(crate) fn nodes(&self) -> impl Iterator<Item = u32> + '_ {
        self.buckets.iter().flatten().filter_map(move |&e| {
            let node = e as u32;
            (self.current[node as usize] == (e >> 32) as u64).then_some(node)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(w: &[f64]) -> Graph {
        let n = w.len() + 1;
        let s: Vec<i64> = (0..w.len() as i64).collect();
        let t: Vec<i64> = (1..=w.len() as i64).collect();
        Graph::from_edges(n, &s, &t, Some(w), true).unwrap()
    }

    #[test]
    fn chooses_the_cheapest_exact_algorithm() {
        let auto = PathAlgorithm::Auto;
        let d = WeightKind::Distance;
        assert!(matches!(
            path_lengths(&g(&[2.5, 2.5]), true, d, auto),
            Ok(Lengths::Uniform)
        ));
        assert!(matches!(
            path_lengths(&g(&[1.0, 3.0]), true, d, auto),
            Ok(Lengths::Integer(_, 3))
        ));
        assert!(matches!(
            path_lengths(&g(&[1.0, 5000.0]), true, d, auto),
            Ok(Lengths::Real(_))
        ));
        assert!(matches!(
            path_lengths(&g(&[0.5, 2.0]), true, d, auto),
            Ok(Lengths::Real(_))
        ));
        assert!(matches!(
            path_lengths(&g(&[0.5, 2.0]), false, d, auto),
            Ok(Lengths::Uniform)
        ));
        // Strengths 1 and 0.5 become lengths 1 and 2.
        assert!(matches!(
            path_lengths(&g(&[1.0, 0.5]), true, WeightKind::Strength, auto),
            Ok(Lengths::Integer(_, 2))
        ));
        assert!(path_lengths(&g(&[0.5]), true, d, PathAlgorithm::Dial).is_err());
        assert!(path_lengths(&g(&[0.0]), true, WeightKind::Strength, auto).is_err());
    }

    #[test]
    fn radix_heap_removes_in_the_order_of_the_indexed_heap() {
        // Dijkstra-like use: keys never fall below the last one removed,
        // except by rounding back onto it; decreases and ties are common.
        let mut rng = crate::rng::Stream::new(3, 0);
        let n = 300;
        let (mut a, mut b) = (IndexedHeap::new(n), RadixHeap::new(n));
        for round in 0..20 {
            let mut last = 0.0f64;
            let mut done = vec![false; n];
            a.push_or_decrease(round as u32, 0.0);
            b.push_or_decrease(round as u32, 0.0);
            loop {
                let (x, y) = (a.pop(), b.pop());
                assert_eq!(
                    x.map(|p| (p.0.to_bits(), p.1)),
                    y.map(|p| (p.0.to_bits(), p.1))
                );
                let Some((d, v)) = x else { break };
                assert!(d >= last);
                last = d;
                done[v as usize] = true;
                for _ in 0..4 {
                    let w = rng.below(n as u64) as u32;
                    if done[w as usize] {
                        continue;
                    }
                    let step = match rng.below(4) {
                        0 => 0.0,
                        1 => 1e-30,
                        2 => rng.below(3) as f64,
                        _ => rng.below(1000) as f64 / 7.0,
                    };
                    a.push_or_decrease(w, d + step);
                    b.push_or_decrease(w, d + step);
                }
            }
        }
    }

    #[test]
    fn heap_keys_follow_the_total_order() {
        let xs = [
            f64::NEG_INFINITY,
            -2.5,
            -0.0,
            0.0,
            1e-300,
            1.0,
            3.25,
            f64::INFINITY,
            f64::NAN,
        ];
        for &a in &xs {
            assert_eq!(unkey(key(a)).to_bits(), a.to_bits());
            for &b in &xs {
                assert_eq!(key(a).cmp(&key(b)), a.total_cmp(&b), "{a} {b}");
            }
        }
    }

    #[test]
    fn indexed_heap_orders_and_decreases() {
        let mut h = IndexedHeap::new(10);
        let mut rng = crate::rng::Stream::new(7, 0);
        let mut best = [f64::INFINITY; 10];
        for _ in 0..200 {
            let v = rng.below(10) as u32;
            let d = (rng.below(20) as f64) / 4.0;
            h.push_or_decrease(v, d);
            best[v as usize] = best[v as usize].min(d);
        }
        let mut expected: Vec<(f64, u32)> = best
            .iter()
            .enumerate()
            .filter(|p| p.1.is_finite())
            .map(|(v, &d)| (d, v as u32))
            .collect();
        expected.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut got = Vec::new();
        while let Some(x) = h.pop() {
            got.push(x);
        }
        assert_eq!(got, expected);
        assert!(h.pos.iter().all(|&p| p == u32::MAX));
        h.push_or_decrease(3, 1.0);
        h.push_or_decrease(4, 2.0);
        assert_eq!(h.nodes().count(), 2);
        h.clear();
        assert!(h.pop().is_none() && h.pos.iter().all(|&p| p == u32::MAX));
    }
}
