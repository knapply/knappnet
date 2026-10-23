//! Newman's leading eigenvector method.
//!
//! Newman (2006, `newman2006finding`) divides a network in two by the
//! signs of the leading eigenvector of the modularity matrix
//! B = A - k k' / (2m), and divides each group again using the
//! generalized modularity matrix of the group, B(g)_ij = B_ij -
//! delta_ij sum_{l in g} B_il. A group is left undivided when the leading
//! eigenvalue of B(g) is not positive, or when the division it suggests
//! does not increase modularity (s' B(g) s / (4m) <= 0, with s the vector
//! of signs). Division continues until no group can be divided.
//!
//! Conventions (matching igraph's `cluster_leading_eigen()` without a step
//! limit):
//! - Undirected networks only. Weights (if used) are added as given and
//!   must be non-negative; parallel edges add up, and a self-loop of
//!   weight w adds 2w to A_ii and to its node's degree, as in
//!   [`crate::metrics::community::modularity`].
//! - "Positive" means larger than 1e-10 times sum_ij |B(g)_ij| / n_g, a
//!   bound on the size of B(g)'s eigenvalues, so the rule does not depend
//!   on the scale of the weights.
//! - The eigenvector is scaled so that its component of largest absolute
//!   value (the first such) is positive; nodes whose component is
//!   positive form one side and the others, including components within
//!   1e-12 of zero relative to the largest, the other side. The sign
//!   of an eigenvector is arbitrary, so tools may place such nodes
//!   differently.
//! - The final partition does not depend on the order in which groups are
//!   divided. Its modularity is computed on the input network.
//! - Eigenvectors come from a dense symmetric eigendecomposition
//!   (Householder tridiagonalization and implicit QL,
//!   [`crate::linalg::symmetric_eigen`]), which costs O(n^3) time and
//!   O(n^2) memory for the first division. igraph uses the sparse ARPACK
//!   solver instead; results agree up to the conventions above.
//! - A network without edges has modularity NaN and one group per node.

use super::{additive_weights, as_labels, relabel, require_undirected};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::linalg::symmetric_eigen;
use crate::matrix::Dense;
use crate::metrics::community::modularity;

/// Largest network for the dense eigendecomposition.
pub const MAX_NODES: usize = 3_000;

/// Result of [`leading_eigen`].
#[derive(Debug, Clone, PartialEq)]
pub struct LeadingEigen {
    /// Group of each node.
    pub membership: Vec<usize>,
    /// Modularity of the partition.
    pub modularity: f64,
    /// Number of divisions made.
    pub splits: usize,
}

/// Leading eigenvector communities of an undirected network. See the
/// module documentation.
pub fn leading_eigen(g: &Graph, weighted: bool) -> Result<LeadingEigen> {
    require_undirected(g, "The leading eigenvector method")?;
    let n = g.node_count();
    if n > MAX_NODES {
        return Err(Error::InvalidArgument(format!(
            "the leading eigenvector method uses a dense eigendecomposition and is limited \
             to {MAX_NODES} nodes"
        )));
    }
    let w = additive_weights(g, weighted)?;
    // Dense symmetric adjacency with the modularity conventions.
    let mut a = vec![0.0f64; n * n];
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        if u == v {
            a[u * n + u] += 2.0 * x;
        } else {
            a[u * n + v] += x;
            a[v * n + u] += x;
        }
    }
    let k: Vec<f64> = (0..n).map(|i| a[i * n..(i + 1) * n].iter().sum()).collect();
    let two_m: f64 = k.iter().sum();
    let mut membership = vec![0usize; n];
    let mut splits = 0;
    if two_m > 0.0 {
        let mut next_label = 1;
        let mut queue: std::collections::VecDeque<Vec<usize>> =
            std::collections::VecDeque::from([(0..n).collect()]);
        while let Some(group) = queue.pop_front() {
            let ng = group.len();
            if ng < 2 {
                continue;
            }
            let mut b = Dense::zeros(ng, ng);
            for (r, &i) in group.iter().enumerate() {
                for (c, &j) in group.iter().enumerate() {
                    b.set(r, c, a[i * n + j] - k[i] * k[j] / two_m);
                }
            }
            for r in 0..ng {
                let row_sum: f64 = b.row(r).iter().sum();
                b.set(r, r, b.get(r, r) - row_sum);
            }
            let scale = b.as_slice().iter().map(|x| x.abs()).sum::<f64>() / ng as f64;
            let (values, vectors) = symmetric_eigen(&b)?;
            let lambda = values[ng - 1];
            if lambda.is_nan() || lambda <= 1e-10 * scale {
                continue;
            }
            let u: Vec<f64> = (0..ng).map(|r| vectors.get(r, ng - 1)).collect();
            let big = (0..ng)
                .max_by(|&x, &y| u[x].abs().total_cmp(&u[y].abs()).then(y.cmp(&x)))
                .expect("non-empty");
            let sign = if u[big] < 0.0 { -1.0 } else { 1.0 };
            let tiny = 1e-12 * u[big].abs();
            let s: Vec<f64> = u
                .iter()
                .map(|&x| if sign * x > tiny { 1.0 } else { -1.0 })
                .collect();
            let gain: f64 = (0..ng)
                .map(|r| s[r] * (0..ng).map(|c| b.get(r, c) * s[c]).sum::<f64>())
                .sum();
            if gain.is_nan() || gain <= 1e-10 * scale * ng as f64 {
                continue;
            }
            let (plus, minus): (Vec<usize>, Vec<usize>) = (0..ng).partition(|&r| s[r] > 0.0);
            if plus.is_empty() || minus.is_empty() {
                continue;
            }
            let plus: Vec<usize> = plus.into_iter().map(|r| group[r]).collect();
            let minus: Vec<usize> = minus.into_iter().map(|r| group[r]).collect();
            for &v in &minus {
                membership[v] = next_label;
            }
            next_label += 1;
            splits += 1;
            queue.push_back(plus);
            queue.push_back(minus);
        }
    } else {
        membership = (0..n).collect();
    }
    let membership = relabel(&membership);
    let q = modularity(g, &as_labels(&membership), 1.0, weighted)?;
    Ok(LeadingEigen {
        membership,
        modularity: q,
        splits,
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
    fn barbell_and_ring_of_cliques() {
        let g = und(6, &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)]);
        let r = leading_eigen(&g, false).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        assert_eq!(r.splits, 1);
        assert!((r.modularity - 5.0 / 14.0).abs() < 1e-14);
        let mut e = Vec::new();
        for c in 0..4i64 {
            for i in 0..5 {
                for j in (i + 1)..5 {
                    e.push((5 * c + i, 5 * c + j));
                }
            }
            e.push((5 * c, (5 * c + 5) % 20));
        }
        let r = leading_eigen(&und(20, &e), false).unwrap();
        for c in 0..4 {
            assert!((0..5).all(|i| r.membership[5 * c + i] == r.membership[5 * c]));
        }
        assert_eq!(r.splits, 3);
    }

    #[test]
    fn indivisible_and_edge_cases() {
        // A complete graph has no division with positive modularity.
        let mut e = Vec::new();
        for i in 0..5 {
            for j in (i + 1)..5 {
                e.push((i, j));
            }
        }
        let r = leading_eigen(&und(5, &e), false).unwrap();
        assert_eq!(r.membership, vec![0; 5]);
        assert_eq!(r.splits, 0);
        let r = leading_eigen(&und(3, &[]), false).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.modularity.is_nan());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(leading_eigen(&d, false).is_err());
        assert!(leading_eigen(&und(2, &[(0, 1)]), true).is_err());
    }
}
