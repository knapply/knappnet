//! Greedy agglomerative modularity optimization ("fast greedy").
//!
//! Clauset, Newman, and Moore (2004, `clauset2004finding`), building on
//! Newman (2004, `newman2004fast`), start with every node in its own group
//! and repeatedly merge the two connected groups whose merger increases
//! modularity the most (or decreases it the least), until every component
//! is a single group. With W_IJ the weight between groups I and J, K_I the
//! total degree of group I, and 2m the total degree, merging I and J
//! changes modularity by
//!
//! ```text
//! dQ = 2 (W_IJ / (2m) - gamma K_I K_J / (2m)^2),
//! ```
//!
//! and the partition with the highest modularity along the way is
//! returned. Candidate mergers are kept in a priority queue, and each
//! group's ties to other groups in an ordered map that is merged into the
//! larger of the two maps, as in Clauset, Newman, and Moore (2004).
//! Queue entries are re-evaluated lazily: a merger only lowers the gains of
//! the merged group's other ties (its total degree grew), so a stored gain
//! is an upper bound unless a tie weight changed, and every change of a
//! tie weight pushes a fresh entry. An entry is acted on only when its
//! stored gain equals the current one; otherwise it is pushed again with
//! the current gain. This chooses the same mergers, in the same order, as
//! updating every affected entry at once, without the O(degree) queue
//! updates per merger that dominate on networks with hubs.
//!
//! Conventions:
//! - Undirected networks only (igraph's `cluster_fast_greedy()`; NetworkX's
//!   `greedy_modularity_communities()` also accepts directed networks).
//! - Weights (if used) are added as given and must be non-negative;
//!   parallel edges add up; a self-loop counts toward its node's degree
//!   and its group's internal weight but never toward a merger, as in
//!   [`crate::metrics::community::modularity`]. `gamma` is the resolution
//!   parameter of Reichardt and Bornholdt (2006,
//!   `reichardt2006statistical`).
//! - Ties: among mergers with the same change in modularity, the one whose
//!   two groups have the smallest lowest-numbered members is made first
//!   (each group is identified by its lowest-numbered node). The numerator
//!   2m W_IJ - gamma K_I K_J is compared exactly when weights are whole
//!   numbers and gamma = 1, so ties are detected exactly in unweighted
//!   networks. igraph and NetworkX break ties by their own internal order,
//!   so on networks with tied mergers partitions can differ while having
//!   similar modularity.
//! - Mergers continue past the first that lowers modularity, until no two
//!   groups are connected (as igraph does); the partition returned is the
//!   first of highest modularity in merge order. NetworkX stops at the
//!   first merger that would not increase modularity, which gives the same
//!   partition unless a later merger raises modularity again.
//! - A network without edges has modularity NaN and one group per node.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use super::{additive_weights, as_labels, relabel, require_undirected};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::community::modularity;
use crate::union_find::UnionFind;

/// Result of [`fast_greedy`].
#[derive(Debug, Clone, PartialEq)]
pub struct FastGreedy {
    /// Group of each node in the partition of highest modularity.
    pub membership: Vec<usize>,
    /// Its modularity on the input network.
    pub modularity: f64,
    /// Every merger in order, each group named by its lowest-numbered node.
    pub merges: Vec<(u32, u32)>,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    gain: f64,
    a: u32,
    b: u32,
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Candidate {}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.gain
            .total_cmp(&other.gain)
            .then_with(|| other.a.cmp(&self.a))
            .then_with(|| other.b.cmp(&self.b))
    }
}

/// Greedy modularity communities of an undirected network. See the module
/// documentation.
pub fn fast_greedy(g: &Graph, resolution: f64, weighted: bool) -> Result<FastGreedy> {
    require_undirected(g, "Greedy modularity optimization")?;
    if !(resolution.is_finite() && resolution >= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "resolution must be finite and non-negative, not {resolution}"
        )));
    }
    let w = additive_weights(g, weighted)?;
    let n = g.node_count();
    let mut nbrs: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); n];
    let mut k = vec![0.0f64; n];
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        k[u] += x;
        k[v] += x;
        if u != v {
            *nbrs[u].entry(v).or_insert(0.0) += x;
            *nbrs[v].entry(u).or_insert(0.0) += x;
        }
    }
    let two_m: f64 = k.iter().sum();
    let gain = |wij: f64, ki: f64, kj: f64| two_m * wij - resolution * ki * kj;
    let mut alive = vec![true; n];
    let mut heap = BinaryHeap::new();
    for (a, row) in nbrs.iter().enumerate() {
        for (&b, &wab) in row {
            if b > a {
                heap.push(Candidate {
                    gain: gain(wab, k[a], k[b]),
                    a: a as u32,
                    b: b as u32,
                });
            }
        }
    }
    let _t = crate::linalg::phase::start("fast greedy: mergers");
    let mut merges: Vec<(u32, u32)> = Vec::new();
    // Running sum of the numerators: modularity after the mergers so far is
    // the initial modularity plus 2 total / (2m)^2.
    let mut total = 0.0f64;
    let mut best_total = 0.0f64;
    let mut best_step = 0usize;
    while let Some(c) = heap.pop() {
        let (a, b) = (c.a as usize, c.b as usize);
        if !alive[a] || !alive[b] {
            continue;
        }
        let wab = match nbrs[a].get(&b) {
            Some(&w) => w,
            None => continue,
        };
        // Lazy evaluation: a stored gain is an upper bound on the current
        // one unless the tie weight changed, and every change of a tie
        // weight pushes a fresh entry. An entry is used only when its gain
        // is current; otherwise it is pushed again with the current gain.
        let current = gain(wab, k[a], k[b]);
        if current.to_bits() != c.gain.to_bits() {
            heap.push(Candidate {
                gain: current,
                a: c.a,
                b: c.b,
            });
            continue;
        }
        merges.push((a as u32, b as u32));
        total += c.gain;
        if total > best_total {
            best_total = total;
            best_step = merges.len();
        }
        // Merge b into a (a < b keeps the lowest-numbered node as name).
        let mut from = std::mem::take(&mut nbrs[b]);
        from.remove(&a);
        nbrs[a].remove(&b);
        let b_row: Vec<usize> = from.keys().copied().collect();
        if from.len() > nbrs[a].len() {
            std::mem::swap(&mut from, &mut nbrs[a]);
        }
        for (x, wx) in from {
            *nbrs[a].entry(x).or_insert(0.0) += wx;
        }
        k[a] += k[b];
        alive[b] = false;
        // Only b's former neighbours see a tie weight change (their tie to
        // b becomes, or is added to, their tie to a). The gains of a's
        // other ties can only fall, as K_a grew, and are re-evaluated
        // lazily.
        for x in b_row {
            let wx = nbrs[a][&x];
            nbrs[x].remove(&b);
            nbrs[x].insert(a, wx);
            let (lo, hi) = if x < a { (x, a) } else { (a, x) };
            heap.push(Candidate {
                gain: gain(wx, k[lo], k[hi]),
                a: lo as u32,
                b: hi as u32,
            });
        }
    }
    // Replay the first best_step mergers.
    let mut sets = UnionFind::new(n);
    for &(a, b) in &merges[..best_step] {
        sets.union(a as usize, b as usize);
    }
    let roots: Vec<usize> = (0..n).map(|v| sets.find(v)).collect();
    let membership = relabel(&roots);
    let q = modularity(g, &as_labels(&membership), resolution, weighted)?;
    Ok(FastGreedy {
        membership,
        modularity: q,
        merges,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn barbell_and_cliques() {
        let g = und(6, &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)]);
        let r = fast_greedy(&g, 1.0, false).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        assert!((r.modularity - 5.0 / 14.0).abs() < 1e-14);
        assert_eq!(r.merges.len(), 5);
        let mut e = Vec::new();
        for c in 0..5i64 {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    e.push((4 * c + i, 4 * c + j));
                }
            }
            e.push((4 * c, (4 * c + 4) % 20));
        }
        let r = fast_greedy(&und(20, &e), 1.0, false).unwrap();
        for c in 0..5 {
            assert!((0..4).all(|i| r.membership[4 * c + i] == r.membership[4 * c]));
        }
        assert_eq!(*r.membership.iter().max().unwrap(), 4);
    }

    #[test]
    fn edge_cases() {
        let r = fast_greedy(&und(3, &[]), 1.0, false).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.modularity.is_nan());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(fast_greedy(&d, 1.0, false).is_err());
        // Resolution 0 merges everything connected.
        let g = und(4, &[(0, 1), (2, 3)]);
        assert_eq!(
            fast_greedy(&g, 0.0, false).unwrap().membership,
            vec![0, 0, 1, 1]
        );
    }
}
