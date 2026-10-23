//! Modularity and Louvain community detection.
//!
//! Modularity (Newman and Girvan 2004, `newman2004finding`; Newman 2006,
//! `newman2006modularity`) compares the weight of edges inside groups with
//! its expectation under the configuration model. With edge weights
//! (Newman 2004, `newman2004analysis`) and the resolution parameter gamma of
//! Reichardt and Bornholdt (2006, `reichardt2006statistical`), for an
//! undirected network with total edge weight m,
//!
//! ```text
//! Q = 1/(2m) sum_ij [A_ij - gamma k_i k_j / (2m)] delta(c_i, c_j),
//! ```
//!
//! and for a directed network (Leicht and Newman 2008, `leicht2008community`)
//!
//! ```text
//! Q = 1/m sum_ij [A_ij - gamma k_i^out k_j^in / m] delta(c_i, c_j).
//! ```
//!
//! Conventions for [`modularity`] (matching igraph's `modularity()` and
//! NetworkX's `community.modularity()`):
//! - Parallel edges each count, and their weights add.
//! - In an undirected network a self-loop of weight w adds 2w to the
//!   diagonal of A and to its node's degree (and counts once, as w, toward
//!   m). In a directed network it adds w to A_ii, to out- and in-degree,
//!   and to m.
//! - Weights must be finite and non-negative (igraph rejects negative
//!   weights; NetworkX accepts them). Without weights every edge has
//!   weight 1.
//! - Group labels are arbitrary non-negative integers; only equality
//!   matters.
//! - A network with no edges (or total weight zero) has undefined
//!   modularity and the result is NaN, as in igraph. NetworkX returns 0.
//!
//! [`louvain`] implements the method of Blondel, Guillaume, Lambiotte, and
//! Lefebvre (2008, `blondel2008fast`). Each level moves single nodes to the
//! neighboring group that most increases Q (staying put on ties) until no
//! move improves Q, then merges each group into one node and repeats on the
//! merged network until a level makes no move. Within a level, nodes are
//! processed from a queue rather than in repeated full sweeps: after a node
//! moves, only its neighbors outside its new group are queued again, since
//! only nodes whose neighborhood changed can have a better move (the "fast
//! local move" procedure of Traag, Waltman, and van Eck 2019,
//! `traag2019louvain`). The level still ends only when no single move
//! improves Q, so the result is a local optimum of the same kind, found with
//! far less work on networks where sweeps converge slowly. Conventions:
//! - Undirected networks only, as in igraph's `cluster_louvain()`; directed
//!   networks are rejected so that no silent symmetrization takes place.
//!   (NetworkX optimizes directed modularity instead.)
//! - Weights and loops follow [`modularity`]. Edges of weight zero are kept
//!   but contribute nothing.
//! - Within each level nodes are visited in a random order drawn from
//!   [`crate::rng::Stream`] with the given seed and the level number as the
//!   stream index, so a seed reproduces the result on any platform. The
//!   method runs on one thread.
//! - A move must raise Q by more than a relative tolerance of 1e-12, which
//!   guarantees termination in floating point.
//! - Groups are labelled 0, 1, ... in order of first appearance in node
//!   order. The modularity returned is [`modularity`] of the final
//!   partition on the input network.
//! - Different implementations and orders reach different local optima,
//!   so partitions are not expected to match other tools exactly.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::rng::Stream;
use crate::weights::additive_weights;

fn check_resolution(resolution: f64) -> Result<()> {
    if !(resolution.is_finite() && resolution >= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "resolution must be finite and non-negative, not {resolution}"
        )));
    }
    Ok(())
}

/// Relabels arbitrary labels to 0..k in order of first appearance,
/// returning the new labels and k. Every partition the engine reports
/// is numbered this way.
pub(crate) fn dense_labels<T: Copy + Eq + std::hash::Hash>(labels: &[T]) -> (Vec<usize>, usize) {
    let mut seen: HashMap<T, usize> = HashMap::new();
    let out = labels
        .iter()
        .map(|&l| {
            let k = seen.len();
            *seen.entry(l).or_insert(k)
        })
        .collect();
    (out, seen.len())
}

/// Modularity of the partition `membership` (one label per node). See the
/// module documentation for the definition and conventions.
pub fn modularity(g: &Graph, membership: &[i64], resolution: f64, weighted: bool) -> Result<f64> {
    let n = g.node_count();
    if membership.len() != n {
        return Err(Error::LengthMismatch {
            what: "membership",
            expected: n,
            found: membership.len(),
        });
    }
    if let Some(&bad) = membership.iter().find(|&&c| c < 0) {
        return Err(Error::InvalidArgument(format!(
            "membership labels must be non-negative, found {bad}"
        )));
    }
    check_resolution(resolution)?;
    let w = additive_weights(g, weighted, "modularity")?;
    let (comm, k) = dense_labels(membership);
    Ok(modularity_dense(g, &comm, k, resolution, w))
}

fn modularity_dense(g: &Graph, comm: &[usize], k: usize, gamma: f64, w: Option<&[f64]>) -> f64 {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut inside = vec![0.0; k];
    let mut k_out = vec![0.0; k];
    let mut k_in = vec![0.0; k];
    let mut m = 0.0;
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (cs, ct) = (comm[src[e] as usize], comm[tgt[e] as usize]);
        m += x;
        k_out[cs] += x;
        k_in[ct] += x;
        if cs == ct {
            inside[cs] += x;
        }
    }
    if m <= 0.0 {
        return f64::NAN;
    }
    let mut q = 0.0;
    for c in 0..k {
        if g.is_directed() {
            q += inside[c] / m - gamma * k_out[c] * k_in[c] / (m * m);
        } else {
            let kc = k_out[c] + k_in[c];
            q += inside[c] / m - gamma * (kc / (2.0 * m)) * (kc / (2.0 * m));
        }
    }
    q
}

/// Result of [`louvain`].
#[derive(Debug, Clone, PartialEq)]
pub struct Louvain {
    /// Group of each node, labelled 0, 1, ... by first appearance.
    pub membership: Vec<usize>,
    /// Modularity of `membership` on the input network.
    pub modularity: f64,
    /// Number of aggregation levels performed.
    pub levels: usize,
}

/// A symmetric weighted network in CSR form whose diagonal entries hold
/// A_ii (twice the loop weight). Node strengths are row sums.
struct Level {
    offsets: Vec<usize>,
    targets: Vec<usize>,
    weights: Vec<f64>,
}

impl Level {
    fn n(&self) -> usize {
        self.offsets.len() - 1
    }

    fn row(&self, v: usize) -> (&[usize], &[f64]) {
        let r = self.offsets[v]..self.offsets[v + 1];
        (&self.targets[r.clone()], &self.weights[r])
    }

    fn from_graph(g: &Graph, w: Option<&[f64]>) -> Level {
        let n = g.node_count();
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let (src, tgt) = (g.edge_sources(), g.edge_targets());
        for e in 0..g.edge_count() {
            let x = w.map_or(1.0, |w| w[e]);
            let (a, b) = (src[e] as usize, tgt[e] as usize);
            if a == b {
                rows[a].push((a, 2.0 * x));
            } else {
                rows[a].push((b, x));
                rows[b].push((a, x));
            }
        }
        Level::from_rows(rows)
    }

    /// Builds a level from unsorted rows, merging repeated targets.
    fn from_rows(rows: Vec<Vec<(usize, f64)>>) -> Level {
        let mut offsets = vec![0];
        let mut targets = Vec::new();
        let mut weights = Vec::new();
        for mut row in rows {
            row.sort_by_key(|&(t, _)| t);
            for (t, x) in row {
                if targets.len() > *offsets.last().unwrap() && targets.last() == Some(&t) {
                    *weights.last_mut().unwrap() += x;
                } else {
                    targets.push(t);
                    weights.push(x);
                }
            }
            offsets.push(targets.len());
        }
        Level {
            offsets,
            targets,
            weights,
        }
    }

    /// Merges each group into a single node.
    fn aggregate(&self, comm: &[usize], k: usize) -> Level {
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); k];
        for v in 0..self.n() {
            let (t, x) = self.row(v);
            for (&u, &x) in t.iter().zip(x) {
                rows[comm[v]].push((comm[u], x));
            }
        }
        Level::from_rows(rows)
    }
}

/// One level of local moving. Returns whether any node moved.
fn move_nodes(level: &Level, comm: &mut [usize], gamma: f64, two_m: f64, rng: &mut Stream) -> bool {
    let n = level.n();
    let strength: Vec<f64> = (0..n).map(|v| level.row(v).1.iter().sum()).collect();
    let mut total = strength.clone();
    let mut order = vec![0usize; n];
    rng.permutation(&mut order);
    // Every node starts in the queue, in random order; a node is queued at
    // most once at a time.
    let mut queue: std::collections::VecDeque<usize> = order.into_iter().collect();
    let mut queued = vec![true; n];
    let mut link = vec![0.0; n];
    let mut touched_flag = vec![false; n];
    let mut touched = Vec::new();
    let mut any = false;
    while let Some(v) = queue.pop_front() {
        queued[v] = false;
        let cv = comm[v];
        let kv = strength[v];
        let (t, x) = level.row(v);
        for (&u, &x) in t.iter().zip(x) {
            if u == v {
                continue;
            }
            let c = comm[u];
            if !touched_flag[c] {
                touched_flag[c] = true;
                touched.push(c);
            }
            link[c] += x;
        }
        total[cv] -= kv;
        // Gain (up to a positive factor) of joining group c.
        let gain = |c: usize| link[c] - gamma * total[c] * kv / two_m;
        let mut best = cv;
        let mut best_gain = gain(cv);
        let tol = 1e-12 * kv.max(f64::MIN_POSITIVE);
        for &c in &touched {
            let g = gain(c);
            if g > best_gain + tol {
                best = c;
                best_gain = g;
            }
        }
        total[best] += kv;
        comm[v] = best;
        for &c in &touched {
            link[c] = 0.0;
            touched_flag[c] = false;
        }
        touched.clear();
        if best != cv {
            any = true;
            for &u in t {
                if u != v && !queued[u] && comm[u] != best {
                    queued[u] = true;
                    queue.push_back(u);
                }
            }
        }
    }
    any
}

/// Louvain community detection on an undirected network. See the module
/// documentation for conventions.
pub fn louvain(g: &Graph, resolution: f64, weighted: bool, seed: u64) -> Result<Louvain> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(
            "Louvain community detection is defined here for undirected networks only; \
             build the network as undirected first"
                .into(),
        ));
    }
    check_resolution(resolution)?;
    let w = additive_weights(g, weighted, "modularity")?;
    let n = g.node_count();
    let mut level = Level::from_graph(g, w);
    let two_m: f64 = level.weights.iter().sum();
    let mut membership: Vec<usize> = (0..n).collect();
    let mut levels = 0;
    if two_m > 0.0 {
        loop {
            let mut comm: Vec<usize> = (0..level.n()).collect();
            let mut rng = Stream::new(seed, levels as u64);
            let moved = {
                let _t = crate::linalg::phase::start("louvain: local moves");
                move_nodes(&level, &mut comm, resolution, two_m, &mut rng)
            };
            if !moved {
                break;
            }
            let _t = crate::linalg::phase::start("louvain: aggregation");
            levels += 1;
            let (dense, k) = dense_labels(&comm);
            for m in membership.iter_mut() {
                *m = dense[*m];
            }
            level = level.aggregate(&dense, k);
        }
    }
    let (membership, k) = dense_labels(&membership);
    let modularity = modularity_dense(g, &membership, k, resolution, w);
    Ok(Louvain {
        membership,
        modularity,
        levels,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)], w: Option<&[f64]>) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, w, false).unwrap()
    }

    fn barbell() -> Graph {
        // Two triangles joined by the edge 2-3 (NetworkX's barbell_graph(3, 0)).
        und(
            6,
            &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)],
            None,
        )
    }

    #[test]
    fn known_values() {
        // NetworkX documents Q = 0.35714285714285715 (= 5/14) for this split.
        let q = modularity(&barbell(), &[0, 0, 0, 1, 1, 1], 1.0, false).unwrap();
        assert!((q - 5.0 / 14.0).abs() < 1e-15);
        // One group: Q = 1 - gamma.
        let q = modularity(&barbell(), &[7; 6], 1.0, false).unwrap();
        assert!(q.abs() < 1e-15);
        let q = modularity(&barbell(), &[7; 6], 0.25, false).unwrap();
        assert!((q - 0.75).abs() < 1e-15);
        // Labels are arbitrary.
        let a = modularity(&barbell(), &[0, 0, 1, 1, 2, 2], 1.0, false).unwrap();
        let b = modularity(&barbell(), &[9, 9, 4, 4, 0, 0], 1.0, false).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn loops_and_direction() {
        // Undirected loop: m = 2, degrees (3, 1), one group per node:
        // sum_c [L_c/m - (k_c/2m)^2] = [1/2 - (3/4)^2] + [0 - (1/4)^2] = -1/8.
        let g = und(2, &[(0, 0), (0, 1)], None);
        let q = modularity(&g, &[0, 1], 1.0, false).unwrap();
        assert!((q + 0.125).abs() < 1e-15);
        // Directed: 0 -> 1, 1 -> 0, 1 -> 2, m = 3, k_out = (1, 2, 0),
        // k_in = (1, 1, 1). Groups {0, 1} and {2}:
        // [2/3 - 3 * 2/9] + [0 - 0 * 1/9] = 0.
        let d = Graph::from_edges(3, &[0, 1, 1], &[1, 0, 2], None, true).unwrap();
        let q = modularity(&d, &[0, 0, 1], 1.0, false).unwrap();
        assert!(q.abs() < 1e-15);
    }

    #[test]
    fn errors_and_edge_cases() {
        let g = barbell();
        assert!(modularity(&g, &[0; 5], 1.0, false).is_err());
        assert!(modularity(&g, &[0, 0, 0, 1, 1, -1], 1.0, false).is_err());
        assert!(modularity(&g, &[0; 6], -1.0, false).is_err());
        assert!(modularity(&g, &[0; 6], 1.0, true).is_err());
        let empty = und(3, &[], None);
        assert!(modularity(&empty, &[0, 1, 2], 1.0, false).unwrap().is_nan());
        let neg = und(2, &[(0, 1)], Some(&[-1.0]));
        assert!(modularity(&neg, &[0, 0], 1.0, true).is_err());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(louvain(&d, 1.0, false, 1).is_err());
        let r = louvain(&empty, 1.0, false, 1).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.modularity.is_nan());
    }

    #[test]
    fn louvain_finds_obvious_groups() {
        let r = louvain(&barbell(), 1.0, false, 3).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        assert!((r.modularity - 5.0 / 14.0).abs() < 1e-15);
        // At resolution 0 everything joins one group (Q = 1 - 0 = 1).
        let r = louvain(&barbell(), 0.0, false, 3).unwrap();
        assert_eq!(r.membership, vec![0; 6]);
        // Weights: a heavy bridge pulls 2 and 3 together.
        let w = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 50.0];
        let g = und(
            6,
            &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)],
            Some(&w),
        );
        let r = louvain(&g, 1.0, true, 3).unwrap();
        assert_eq!(r.membership[2], r.membership[3]);
        let q = modularity(
            &g,
            &r.membership.iter().map(|&c| c as i64).collect::<Vec<_>>(),
            1.0,
            true,
        )
        .unwrap();
        assert_eq!(q, r.modularity);
    }

    #[test]
    fn louvain_is_reproducible() {
        // A ring of 8 cliques of size 5.
        let mut e = Vec::new();
        for c in 0..8i64 {
            for i in 0..5 {
                for j in (i + 1)..5 {
                    e.push((5 * c + i, 5 * c + j));
                }
            }
            e.push((5 * c, (5 * (c + 1)) % 40));
        }
        let g = und(40, &e, None);
        let a = louvain(&g, 1.0, false, 11).unwrap();
        let b = louvain(&g, 1.0, false, 11).unwrap();
        assert_eq!(a, b);
        // Each clique is its own group here.
        for c in 0..8 {
            let m = a.membership[5 * c];
            assert!((0..5).all(|i| a.membership[5 * c + i] == m));
        }
    }
}
