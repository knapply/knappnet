//! A network that can be toggled one dyad at a time.
//!
//! Markov chain Monte Carlo for exponential random graph models proposes a
//! change to a single dyad, scores it, and either applies or rejects it
//! (Snijders 2002, `snijders2002markov`; Hunter et al. 2008,
//! `hunter2008ergm`). The immutable CSR storage of [`crate::graph::Graph`]
//! cannot do that, so this module keeps a mutable representation:
//!
//! - a bit matrix for constant-time tie tests while it fits in 256 MB (the
//!   bound [`crate::bitset`] uses), and a hash map of ties otherwise;
//! - neighbor lists (out- and in-neighbors for directed networks) for
//!   iterating over a node's ties;
//! - the tie list, with each tie's position, so a uniformly random tie can
//!   be drawn in constant time for tie/no-tie proposals;
//! - the set of fixed dyads (structural zeros or ones) that no proposal may
//!   change, and the resulting count of free dyads.
//!
//! Conventions: the network is simple. Self-loops are never allowed, and
//! parallel edges collapse to one tie when a state is built from an edge
//! list (this is also how statnet's `network` objects read an edge list
//! into a binary sociomatrix). Weights are ignored: ERGMs for binary ties
//! are defined on the presence or absence of ties. Undirected ties are
//! stored once, as the dyad (min, max).

use std::collections::HashMap;

use crate::error::{Error, Result};

/// Largest dyad bit matrix kept alongside the hash map, in bytes (the
/// same bound as [`crate::bitset::MAX_BYTES`]).
const MAX_BITS_BYTES: usize = 256 << 20;

/// A square bit matrix that can clear bits as well as set them (the shared
/// [`crate::bitset::BitMatrix`] is write-once, for immutable graphs).
#[derive(Debug, Clone)]
pub(crate) struct DyadBits {
    words_per_row: usize,
    bits: Vec<u64>,
}

impl DyadBits {
    fn new(n: usize) -> DyadBits {
        let words_per_row = n.saturating_add(63) / 64;
        DyadBits {
            words_per_row,
            bits: vec![0; words_per_row * n],
        }
    }

    fn bytes(n: usize) -> usize {
        (n.saturating_add(63) / 64)
            .saturating_mul(8)
            .saturating_mul(n)
    }

    #[inline]
    fn set(&mut self, row: usize, col: usize) {
        self.bits[row * self.words_per_row + col / 64] |= 1u64 << (col % 64);
    }

    #[inline]
    fn clear(&mut self, row: usize, col: usize) {
        self.bits[row * self.words_per_row + col / 64] &= !(1u64 << (col % 64));
    }

    #[inline]
    fn get(&self, row: usize, col: usize) -> bool {
        self.bits[row * self.words_per_row + col / 64] >> (col % 64) & 1 == 1
    }
}

/// Encodes an ordered pair as one 64-bit key.
#[inline]
fn key(i: u32, j: u32) -> u64 {
    (i as u64) << 32 | j as u64
}

/// Which dyads a model may change.
#[derive(Debug, Clone)]
pub struct Constraints {
    fixed: DyadBits,
    n: usize,
    directed: bool,
    /// Every free dyad, as (tail, head), in canonical order.
    free: Vec<(u32, u32)>,
}

impl Constraints {
    /// Marks the dyads where `fixed[i * n + j] != 0` (row-major) as fixed.
    /// For undirected networks the entry for either ordering fixes the dyad.
    pub fn new(n: usize, directed: bool, fixed: &[i32]) -> Result<Constraints> {
        if fixed.len() != n * n {
            return Err(Error::LengthMismatch {
                what: "fixed dyads",
                expected: n * n,
                found: fixed.len(),
            });
        }
        let mut bits = DyadBits::new(n);
        for i in 0..n {
            for j in 0..n {
                if fixed[i * n + j] != 0 {
                    bits.set(i, j);
                    if !directed {
                        bits.set(j, i);
                    }
                }
            }
        }
        let mut free = Vec::new();
        for i in 0..n {
            for j in 0..n {
                if i == j || (!directed && j < i) || bits.get(i, j) {
                    continue;
                }
                free.push((i as u32, j as u32));
            }
        }
        Ok(Constraints {
            fixed: bits,
            n,
            directed,
            free,
        })
    }

    /// Whether dyad (i, j) may not be changed.
    #[inline]
    pub fn is_fixed(&self, i: usize, j: usize) -> bool {
        self.fixed.get(i, j)
    }

    /// The free dyads in canonical order.
    pub fn free_dyads(&self) -> &[(u32, u32)] {
        &self.free
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.n
    }

    /// Whether the constraints were built for a directed network.
    pub fn is_directed(&self) -> bool {
        self.directed
    }
}

/// A binary network that supports toggling one dyad at a time.
#[derive(Debug, Clone)]
pub struct NetState {
    n: usize,
    directed: bool,
    bits: Option<DyadBits>,
    /// Out-neighbors (all neighbors when undirected), in arbitrary order.
    out: Vec<Vec<u32>>,
    /// In-neighbors; empty for undirected networks.
    inc: Vec<Vec<u32>>,
    /// Ties as (tail, head); undirected ties as (min, max).
    ties: Vec<(u32, u32)>,
    /// Position of each tie in `ties`, by canonical key.
    pos: HashMap<u64, usize>,
}

impl NetState {
    /// An empty network.
    pub fn empty(n: usize, directed: bool) -> NetState {
        let bits = if DyadBits::bytes(n) <= MAX_BITS_BYTES {
            Some(DyadBits::new(n))
        } else {
            None
        };
        NetState {
            n,
            directed,
            bits,
            out: vec![Vec::new(); n],
            inc: if directed {
                vec![Vec::new(); n]
            } else {
                Vec::new()
            },
            ties: Vec::new(),
            pos: HashMap::new(),
        }
    }

    /// Builds a state from an edge list. Self-loops are rejected; parallel
    /// edges collapse to one tie.
    pub fn from_edges(n: usize, src: &[i64], tgt: &[i64], directed: bool) -> Result<NetState> {
        if tgt.len() != src.len() {
            return Err(Error::LengthMismatch {
                what: "targets",
                expected: src.len(),
                found: tgt.len(),
            });
        }
        if n > u32::MAX as usize {
            return Err(Error::TooLarge {
                what: "node count",
                value: n,
            });
        }
        let mut s = NetState::empty(n, directed);
        for (e, (&a, &b)) in src.iter().zip(tgt).enumerate() {
            for node in [a, b] {
                if node < 0 || node as u64 >= n as u64 {
                    return Err(Error::NodeOutOfRange { edge: e, node, n });
                }
            }
            if a == b {
                return Err(Error::InvalidArgument(format!(
                    "edge {e} is a self-loop; exponential random graph models are defined \
                     on networks without self-loops"
                )));
            }
            let (i, j) = (a as usize, b as usize);
            if !s.has(i, j) {
                s.toggle(i, j);
            }
        }
        Ok(s)
    }

    /// Number of nodes.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.n
    }

    /// Whether ties are directed.
    #[inline]
    pub fn is_directed(&self) -> bool {
        self.directed
    }

    /// Number of ties.
    #[inline]
    pub fn tie_count(&self) -> usize {
        self.ties.len()
    }

    /// The ties as (tail, head) pairs, in no particular order.
    #[inline]
    pub fn ties(&self) -> &[(u32, u32)] {
        &self.ties
    }

    /// The canonical form of a dyad: as given when directed, (min, max)
    /// otherwise.
    #[inline]
    pub fn canonical(&self, i: usize, j: usize) -> (u32, u32) {
        if self.directed || i < j {
            (i as u32, j as u32)
        } else {
            (j as u32, i as u32)
        }
    }

    /// Whether the tie i -> j (or the undirected tie {i, j}) exists.
    #[inline]
    pub fn has(&self, i: usize, j: usize) -> bool {
        match &self.bits {
            Some(b) => b.get(i, j),
            None => {
                let (a, b) = self.canonical(i, j);
                self.pos.contains_key(&key(a, b))
            }
        }
    }

    /// Out-neighbors of `v` (all neighbors when undirected), unordered.
    #[inline]
    pub fn out_neighbors(&self, v: usize) -> &[u32] {
        &self.out[v]
    }

    /// In-neighbors of `v` (all neighbors when undirected), unordered.
    #[inline]
    pub fn in_neighbors(&self, v: usize) -> &[u32] {
        if self.directed {
            &self.inc[v]
        } else {
            &self.out[v]
        }
    }

    /// Out-degree (degree when undirected).
    #[inline]
    pub fn out_degree(&self, v: usize) -> usize {
        self.out[v].len()
    }

    /// In-degree (degree when undirected).
    #[inline]
    pub fn in_degree(&self, v: usize) -> usize {
        if self.directed {
            self.inc[v].len()
        } else {
            self.out[v].len()
        }
    }

    /// Total degree: in plus out when directed, degree otherwise.
    #[inline]
    pub fn degree(&self, v: usize) -> usize {
        if self.directed {
            self.out[v].len() + self.inc[v].len()
        } else {
            self.out[v].len()
        }
    }

    /// Flips the state of dyad (i, j) and returns whether it is now a tie.
    pub fn toggle(&mut self, i: usize, j: usize) -> bool {
        debug_assert!(i != j);
        let (a, b) = self.canonical(i, j);
        let k = key(a, b);
        if let Some(p) = self.pos.remove(&k) {
            let last = self.ties.len() - 1;
            self.ties.swap(p, last);
            self.ties.pop();
            if p != last {
                let moved = self.ties[p];
                self.pos.insert(key(moved.0, moved.1), p);
            }
            if let Some(bits) = &mut self.bits {
                bits.clear(a as usize, b as usize);
                if !self.directed {
                    bits.clear(b as usize, a as usize);
                }
            }
            remove_from(&mut self.out[a as usize], b);
            if self.directed {
                remove_from(&mut self.inc[b as usize], a);
            } else {
                remove_from(&mut self.out[b as usize], a);
            }
            false
        } else {
            self.pos.insert(k, self.ties.len());
            self.ties.push((a, b));
            if let Some(bits) = &mut self.bits {
                bits.set(a as usize, b as usize);
                if !self.directed {
                    bits.set(b as usize, a as usize);
                }
            }
            self.out[a as usize].push(b);
            if self.directed {
                self.inc[b as usize].push(a);
            } else {
                self.out[b as usize].push(a);
            }
            true
        }
    }

    /// Number of distinct dyads: n(n-1) directed, n(n-1)/2 undirected.
    pub fn dyad_count(&self) -> usize {
        let n = self.n;
        if self.directed {
            n * n.saturating_sub(1)
        } else {
            n * n.saturating_sub(1) / 2
        }
    }

    /// Every dyad in canonical order: (i, j) with i != j directed, i < j
    /// undirected.
    pub fn dyads(&self) -> Vec<(u32, u32)> {
        let n = self.n;
        let mut out = Vec::with_capacity(self.dyad_count());
        for i in 0..n {
            for j in 0..n {
                if i == j || (!self.directed && j < i) {
                    continue;
                }
                out.push((i as u32, j as u32));
            }
        }
        out
    }

    /// Number of nodes `i` and `j` are both tied to (common neighbors,
    /// ignoring direction), excluding the pair itself.
    pub fn shared_partners_undirected(&self, i: usize, j: usize) -> usize {
        let (a, b) = if self.out[i].len() <= self.out[j].len() {
            (i, j)
        } else {
            (j, i)
        };
        self.out[a]
            .iter()
            .filter(|&&u| u as usize != b && self.has(u as usize, b))
            .count()
    }

    /// Number of two-paths i -> k -> j (outgoing two-path shared partners of
    /// the ordered pair).
    pub fn two_paths(&self, i: usize, j: usize) -> usize {
        self.out[i]
            .iter()
            .filter(|&&k| k as usize != j && self.has(k as usize, j))
            .count()
    }
}

#[inline]
fn remove_from(list: &mut Vec<u32>, x: u32) {
    if let Some(p) = list.iter().position(|&y| y == x) {
        list.swap_remove(p);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toggles_maintain_every_view() {
        let mut s = NetState::from_edges(5, &[0, 1, 2, 2], &[1, 2, 3, 3], false).unwrap();
        assert_eq!(s.tie_count(), 3); // the parallel edge collapsed
        assert!(s.has(1, 0) && s.has(2, 3) && !s.has(0, 2));
        assert_eq!(s.degree(2), 2);
        assert!(!s.toggle(3, 2));
        assert_eq!(s.tie_count(), 2);
        assert!(!s.has(2, 3));
        assert_eq!(s.degree(3), 0);
        assert!(s.toggle(4, 0));
        assert_eq!(s.ties().iter().filter(|t| **t == (0, 4)).count(), 1);
        assert_eq!(s.shared_partners_undirected(0, 2), 1);
        assert_eq!(s.dyad_count(), 10);
    }

    #[test]
    fn directed_views() {
        let mut s = NetState::from_edges(4, &[0, 1, 2], &[1, 2, 0], true).unwrap();
        assert!(s.has(0, 1) && !s.has(1, 0));
        assert_eq!(s.in_degree(0), 1);
        assert_eq!(s.out_degree(0), 1);
        assert_eq!(s.two_paths(0, 2), 1);
        s.toggle(1, 0);
        assert_eq!(s.degree(0), 3);
        assert_eq!(s.dyads().len(), 12);
    }

    #[test]
    fn rejects_loops_and_bad_indices() {
        assert!(NetState::from_edges(3, &[0], &[0], true).is_err());
        assert!(NetState::from_edges(3, &[0], &[3], true).is_err());
    }

    #[test]
    fn constraints_list_free_dyads() {
        let mut fixed = vec![0; 9];
        fixed[1] = 1; // (0, 1)
        let c = Constraints::new(3, false, &fixed).unwrap();
        assert!(c.is_fixed(0, 1) && c.is_fixed(1, 0));
        assert_eq!(c.free_dyads(), &[(0, 2), (1, 2)]);
        let c = Constraints::new(3, true, &fixed).unwrap();
        assert!(c.is_fixed(0, 1) && !c.is_fixed(1, 0));
        assert_eq!(c.free_dyads().len(), 5);
    }
}
