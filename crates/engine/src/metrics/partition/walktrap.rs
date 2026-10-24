//! Walktrap: agglomerative communities from random-walk distances.
//!
//! Pons and Latapy (2006, `pons2006computing`) give every node a loop,
//! let P = D^-1 A be the transition matrix of the random walk (D the
//! diagonal of degrees, loops included) and compare nodes by where walks of
//! t steps from them end up:
//!
//! ```text
//! r_ij = || D^-1/2 P^t_i. - D^-1/2 P^t_j. ||,
//! ```
//!
//! with P^t_C. the mean of P^t_i. over the members of a group C. Starting
//! from singletons, they repeatedly merge the two adjacent groups (joined
//! by at least one edge) that least increase the mean squared distance of
//! nodes to their group (Ward's criterion),
//!
//! ```text
//! delta sigma(C1, C2) = (1/n) |C1| |C2| / (|C1| + |C2|) r^2_C1C2,
//! ```
//!
//! and the vector of a merged group is the size-weighted mean of its
//! parts. Here every delta sigma is evaluated from the group vectors
//! directly (rather than through the Lance-Williams-Jambu update, which is
//! the same quantity in exact arithmetic), so every merge is the exact
//! minimum up to rounding. The partition returned is the level of highest
//! modularity along the merges, as in the paper.
//!
//! Conventions:
//! - Undirected networks only, as for the other hierarchical methods here;
//!   igraph's `cluster_walktrap()` ignores edge directions instead.
//! - Parallel edges add their weights, and a self-loop of the network adds
//!   its weight to A_ii. Every node then also gets the loop of the method,
//!   of weight equal to the mean weight of its other incident edges (1
//!   without weights), and an isolated node a loop of weight 1, as igraph
//!   documents for isolated nodes.
//! - Weights (if used) must be positive, as igraph requires.
//! - Groups are named as in igraph's `merges` (leaves 0..n, merge i creates
//!   group n + i); see [`super::dendrogram`]. Among merges with equal
//!   delta sigma the pair with the smaller group names is made first (the
//!   smaller lower name, then the smaller higher name). Networks with
//!   symmetries (structurally equivalent nodes) have many such ties, and
//!   values equal in exact arithmetic can differ in the last bits of
//!   floating point, so values within a relative 1e-9 of the smallest
//!   are taken as tied. igraph documents no tie rule and may merge tied
//!   pairs in a different order.
//! - igraph's `cluster_walktrap()` occasionally merges a pair that is not
//!   the minimum: on some networks its merge sequence skips an adjacent
//!   pair whose exact delta sigma is smaller by far more than rounding
//!   (20% in one validation case) and merges it a few steps later. The
//!   merges here follow the paper's rule exactly, as checked against an
//!   implementation in exact rational arithmetic.
//! - Only adjacent groups merge, so a network with c components gets
//!   n - c merges (an incomplete dendrogram).
//! - Modularity along the merges is [`crate::metrics::community::modularity`]
//!   of the input network (resolution 1; with the same weights), not of the
//!   network with the added loops; the partition is the first level of
//!   highest modularity. Without edges it is NaN, with one group per node.
//! - The walk vectors are dense, so memory grows as n^2 (8 n^2 bytes);
//!   they are computed on `threads` threads, which does not change the
//!   result.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use super::{additive_weights, cut_at_max_modularity, dendrogram_modularity, require_undirected};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel::for_each_item;

/// Result of [`walktrap`].
#[derive(Debug, Clone, PartialEq)]
pub struct Walktrap {
    /// Group of each node in the level of highest modularity.
    pub membership: Vec<usize>,
    /// Its modularity.
    pub modularity: f64,
    /// The merges, in igraph's numbering.
    pub merges: Vec<(usize, usize)>,
    /// Modularity before any merge (entry 0) and after each merge.
    pub path: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    delta: f64,
    a: usize,
    b: usize,
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
    // BinaryHeap is a max-heap: the smallest delta (then names) is greatest.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .delta
            .total_cmp(&self.delta)
            .then_with(|| other.a.cmp(&self.a))
            .then_with(|| other.b.cmp(&self.b))
    }
}

/// Two values of delta sigma equal to a relative 1e-9 (or both below
/// 1e-18) are tied.
fn tied(x: f64, y: f64) -> bool {
    (x - y).abs() <= 1e-9 * x.abs().max(y.abs()) + 1e-18
}

fn squared_distance(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| (a - b) * (a - b)).sum()
}

/// Walktrap communities of an undirected network with walks of `steps`
/// steps. See the module documentation.
pub fn walktrap(g: &Graph, steps: usize, weighted: bool, threads: usize) -> Result<Walktrap> {
    require_undirected(g, "Walktrap")?;
    if steps == 0 {
        return Err(Error::InvalidArgument(
            "the number of random-walk steps must be at least 1".into(),
        ));
    }
    let w = additive_weights(g, weighted)?;
    if let Some(w) = w {
        if let Some(e) = w.iter().position(|&x| x <= 0.0) {
            return Err(Error::NonPositiveWeight {
                edge: e,
                weight: w[e],
            });
        }
    }
    let n = g.node_count();
    // Adjacency rows with weights added up; the method's loop is added
    // below.
    let mut rows: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); n];
    let mut other = vec![(0.0f64, 0usize); n];
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        *rows[u].entry(v).or_insert(0.0) += x;
        if u != v {
            *rows[v].entry(u).or_insert(0.0) += x;
            other[u].0 += x;
            other[u].1 += 1;
            other[v].0 += x;
            other[v].1 += 1;
        }
    }
    for (v, row) in rows.iter_mut().enumerate() {
        let (total, count) = other[v];
        let lw = if count == 0 {
            1.0
        } else {
            total / count as f64
        };
        *row.entry(v).or_insert(0.0) += lw;
    }
    let adj: Vec<Vec<(usize, f64)>> = rows
        .iter()
        .map(|r| r.iter().map(|(&u, &x)| (u, x)).collect())
        .collect();
    let deg: Vec<f64> = adj
        .iter()
        .map(|r| r.iter().map(|&(_, x)| x).sum())
        .collect();
    let inv_sqrt: Vec<f64> = deg.iter().map(|d| 1.0 / d.sqrt()).collect();

    // D^-1/2 P^t_i. for every node i: walks of `steps` steps from i.
    let parts = for_each_item(
        n,
        threads,
        Vec::new,
        |out: &mut Vec<(usize, Vec<f64>)>, i| {
            let mut x = vec![0.0f64; n];
            let mut y = vec![0.0f64; n];
            x[i] = 1.0;
            for _ in 0..steps {
                y.iter_mut().for_each(|v| *v = 0.0);
                for (j, &xj) in x.iter().enumerate() {
                    if xj != 0.0 {
                        let s = xj / deg[j];
                        for &(k, a) in &adj[j] {
                            y[k] += s * a;
                        }
                    }
                }
                std::mem::swap(&mut x, &mut y);
            }
            for (v, s) in x.iter_mut().zip(&inv_sqrt) {
                *v *= s;
            }
            out.push((i, x));
        },
    );
    let mut vecs: Vec<Vec<f64>> = vec![Vec::new(); 2 * n];
    for part in parts {
        for (i, x) in part {
            vecs[i] = x;
        }
    }
    let mut size = vec![0usize; 2 * n];
    size[..n].iter_mut().for_each(|s| *s = 1);
    let nf = n as f64;
    let delta = |vecs: &[Vec<f64>], size: &[usize], a: usize, b: usize| {
        let (sa, sb) = (size[a] as f64, size[b] as f64);
        sa * sb / (sa + sb) * squared_distance(&vecs[a], &vecs[b]) / nf
    };
    // Neighbouring groups of each group.
    let mut nbrs: Vec<Vec<usize>> = vec![Vec::new(); 2 * n];
    let mut heap = BinaryHeap::new();
    for (a, row) in adj.iter().enumerate() {
        for &(b, _) in row {
            if b != a {
                nbrs[a].push(b);
                if b > a {
                    heap.push(Candidate {
                        delta: delta(&vecs, &size, a, b),
                        a,
                        b,
                    });
                }
            }
        }
    }
    let mut alive = vec![false; 2 * n];
    alive[..n].iter_mut().for_each(|s| *s = true);
    let mut merges: Vec<(usize, usize)> = Vec::new();
    let mut near: Vec<Candidate> = Vec::new();
    while let Some(c) = heap.pop() {
        if !alive[c.a] || !alive[c.b] {
            continue;
        }
        // Among the live candidates whose delta sigma equals the smallest
        // to a relative 1e-9, take the one with the smallest names.
        near.clear();
        while let Some(top) = heap.peek() {
            if !tied(top.delta, c.delta) {
                break;
            }
            let top = heap.pop().expect("peeked");
            if alive[top.a] && alive[top.b] {
                near.push(top);
            }
        }
        let mut c = c;
        for x in near.iter_mut() {
            if (x.a, x.b) < (c.a, c.b) {
                std::mem::swap(x, &mut c);
            }
        }
        for x in near.drain(..) {
            heap.push(x);
        }
        let (a, b) = (c.a, c.b);
        let new = n + merges.len();
        merges.push((a, b));
        alive[a] = false;
        alive[b] = false;
        let (sa, sb) = (size[a] as f64, size[b] as f64);
        let va = std::mem::take(&mut vecs[a]);
        let vb = std::mem::take(&mut vecs[b]);
        vecs[new] = va
            .iter()
            .zip(&vb)
            .map(|(x, y)| (sa * x + sb * y) / (sa + sb))
            .collect();
        size[new] = size[a] + size[b];
        let mut around: Vec<usize> = std::mem::take(&mut nbrs[a]);
        around.extend(std::mem::take(&mut nbrs[b]));
        around.retain(|&x| alive[x]);
        around.sort_unstable();
        around.dedup();
        for &x in &around {
            for y in nbrs[x].iter_mut() {
                if *y == a || *y == b {
                    *y = new;
                }
            }
            nbrs[x].sort_unstable();
            nbrs[x].dedup();
            heap.push(Candidate {
                delta: delta(&vecs, &size, x, new),
                a: x,
                b: new,
            });
        }
        nbrs[new] = around;
        alive[new] = true;
    }
    let path = dendrogram_modularity(g, &merges, 1.0, weighted)?;
    let (membership, modularity) = cut_at_max_modularity(g, &merges, 1.0, weighted)?;
    Ok(Walktrap {
        membership,
        modularity,
        merges,
        path,
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
    fn two_triangles() {
        let g = und(6, &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)]);
        let r = walktrap(&g, 4, false, 1).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        assert!((r.modularity - 5.0 / 14.0).abs() < 1e-14);
        assert_eq!(r.merges.len(), 5);
        assert_eq!(r.path.len(), 6);
        assert_eq!(r, walktrap(&g, 4, false, 3).unwrap());
    }

    #[test]
    fn first_merge_by_hand() {
        // Path 0-1-2 with loops: A = [[1,1,0],[1,1,1],[0,1,1]], d = (2,3,2).
        // One step: P_0. = (1/2, 1/2, 0), P_1. = (1/3, 1/3, 1/3),
        // P_2. = (0, 1/2, 1/2). r^2_01 = (1/6)^2/2 + (1/6)^2/3 + (1/3)^2/2,
        // the same as r^2_12 by symmetry and smaller than r^2_02, so the
        // first merge joins the middle node with an end.
        let g = und(3, &[(0, 1), (1, 2)]);
        let r = walktrap(&g, 1, false, 1).unwrap();
        assert!(r.merges[0] == (0, 1) || r.merges[0] == (1, 2));
        assert_eq!(r.merges[1].1, 3);
    }

    #[test]
    fn edge_cases() {
        let r = walktrap(&und(3, &[]), 4, false, 1).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.modularity.is_nan());
        assert!(r.merges.is_empty());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(walktrap(&d, 4, false, 1).is_err());
        assert!(walktrap(&und(2, &[(0, 1)]), 0, false, 1).is_err());
        // Two components: an incomplete dendrogram.
        let r = walktrap(&und(4, &[(0, 1), (2, 3)]), 4, false, 1).unwrap();
        assert_eq!(r.merges.len(), 2);
        assert_eq!(r.membership, vec![0, 0, 1, 1]);
    }
}
