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

/// A 4-ary min-heap of nodes keyed by tentative distance, with a position
/// index so that a node's key can be lowered in place (Johnson 1975,
/// `johnson1975priority`). Dijkstra's algorithm then holds at most one
/// entry per node instead of one per relaxation, and a 4-ary layout halves
/// the depth of a binary heap, so a sift touches fewer cache lines.
///
/// Entries are ordered by (distance, node index), a total order, so nodes
/// are settled in exactly the order a binary heap with the same ordering
/// and lazy deletion settles them.
#[derive(Debug, Clone)]
pub(crate) struct IndexedHeap {
    items: Vec<(f64, u32)>,
    /// Position of each node in `items`; `u32::MAX` when absent.
    pos: Vec<u32>,
}

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
    fn less(a: (f64, u32), b: (f64, u32)) -> bool {
        match a.0.total_cmp(&b.0) {
            std::cmp::Ordering::Less => true,
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Equal => a.1 < b.1,
        }
    }

    /// Inserts `node` with key `dist`, or lowers its key to `dist` if it
    /// is present with a larger key. A larger `dist` is ignored.
    #[inline]
    pub(crate) fn push_or_decrease(&mut self, node: u32, dist: f64) {
        let p = self.pos[node as usize];
        let i = if p == u32::MAX {
            self.items.push((dist, node));
            self.items.len() - 1
        } else {
            let i = p as usize;
            if !Self::less((dist, node), self.items[i]) {
                return;
            }
            self.items[i].0 = dist;
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
        Some(top)
    }

    /// Removes every entry.
    pub(crate) fn clear(&mut self) {
        for &(_, v) in &self.items {
            self.pos[v as usize] = u32::MAX;
        }
        self.items.clear();
    }

    /// The entries in no particular order.
    pub(crate) fn entries(&self) -> &[(f64, u32)] {
        &self.items
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
        assert_eq!(h.entries().len(), 2);
        h.clear();
        assert!(h.pop().is_none() && h.pos.iter().all(|&p| p == u32::MAX));
    }
}
