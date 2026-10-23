//! The continuous core-periphery model.
//!
//! Borgatti and Everett (1999, `borgatti1999models`, section 4) give every
//! node a coreness c_i and take the pattern matrix delta_ij = c_i c_j, so
//! that ties are expected between core nodes, less so between core and
//! peripheral nodes, and least between peripheral nodes. Coreness is
//! chosen to maximize the Pearson correlation between the data and the
//! pattern. The authors maximize it numerically (with a Fletcher-Powell
//! procedure); this module uses the quasi-Newton method of Broyden,
//! Fletcher, Goldfarb, and Shanno with a backtracking line search that
//! enforces the Armijo condition (Nocedal and Wright 2006,
//! `nocedal2006numerical`, chapters 3 and 6), using the exact gradient of
//! the correlation.
//!
//! Conventions:
//! - Pairs and values as in [`super::discrete`]: pairs of distinct nodes,
//!   each unordered pair once for undirected networks and both ordered
//!   pairs for directed networks; binary values, or total weights.
//! - Coreness is constrained to be non-negative, by writing c_i = t_i^2 and
//!   optimizing over t. Borgatti and Everett do not state a sign
//!   constraint, but without one the product pattern also fits two
//!   separate dense groups (with coreness of opposite signs), which is not
//!   a core-periphery pattern.
//! - The correlation does not change when c is multiplied by a positive
//!   number, so the result is scaled to unit Euclidean norm.
//! - Optimization starts from c proportional to each node's total value
//!   plus 1 (degree plus 1 for binary data) and stops when the largest
//!   gradient component is below 1e-9 or no step improves the correlation;
//!   at most `max_iter` iterations are made. The result is a local maximum,
//!   found deterministically. Cost O(n^2) per iteration and O(n^2) memory.
//! - A network whose values are all equal is rejected.

use super::PairRows;
use crate::error::{Error, Result};
use crate::graph::Graph;

/// Result of [`continuous_core_periphery`].
#[derive(Debug, Clone, PartialEq)]
pub struct ContinuousCorePeriphery {
    /// Coreness of each node, non-negative, with unit Euclidean norm.
    pub coreness: Vec<f64>,
    /// Correlation between the data and the pattern c_i c_j.
    pub correlation: f64,
    /// Iterations made.
    pub iterations: usize,
}

struct Objective {
    n: usize,
    directed: bool,
    /// Dense values a_ij (symmetric for undirected networks).
    a: Vec<f64>,
    pairs: f64,
    mean_a: f64,
    sxx: f64,
}

impl Objective {
    fn new(g: &Graph, weighted: bool) -> Result<Objective> {
        let pr = PairRows::new(g, weighted)?;
        let n = g.node_count();
        let mut a = vec![0.0; n * n];
        let directed = g.is_directed();
        let w = if weighted { g.edge_weights() } else { None };
        for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
            if s == t {
                continue;
            }
            let (s, t) = (s as usize, t as usize);
            let x = w.map_or(1.0, |w| w[e]);
            if weighted {
                a[s * n + t] += x;
                if !directed {
                    a[t * n + s] += x;
                }
            } else {
                a[s * n + t] = 1.0;
                if !directed {
                    a[t * n + s] = 1.0;
                }
            }
        }
        let mean_a = pr.sum / pr.pairs;
        let sxx = pr.sum_sq - pr.pairs * mean_a * mean_a;
        Ok(Objective {
            n,
            directed,
            a,
            pairs: pr.pairs,
            mean_a,
            sxx,
        })
    }

    /// Calls f(i, j, a_ij) for every pair in the fit.
    fn for_pairs(&self, mut f: impl FnMut(usize, usize, f64)) {
        let n = self.n;
        for i in 0..n {
            for j in 0..n {
                if i == j || (!self.directed && j < i) {
                    continue;
                }
                f(i, j, self.a[i * n + j]);
            }
        }
    }

    /// Correlation for coreness c, and its gradient with respect to c.
    fn eval(&self, c: &[f64]) -> (f64, Vec<f64>) {
        let (mut sy, mut syy, mut sxy) = (0.0, 0.0, 0.0);
        self.for_pairs(|i, j, a| {
            let y = c[i] * c[j];
            sy += y;
            syy += y * y;
            sxy += (a - self.mean_a) * y;
        });
        let mean_y = sy / self.pairs;
        let vyy = syy - self.pairs * mean_y * mean_y;
        let mut grad = vec![0.0; self.n];
        if !(vyy > 0.0) {
            return (f64::NAN, grad);
        }
        let den = (self.sxx * vyy).sqrt();
        let r = sxy / den;
        // dr/dy_p = (a_p - mean_a) / den - r (y_p - mean_y) / vyy.
        self.for_pairs(|i, j, a| {
            let y = c[i] * c[j];
            let d = (a - self.mean_a) / den - r * (y - mean_y) / vyy;
            grad[i] += d * c[j];
            grad[j] += d * c[i];
        });
        (r, grad)
    }
}

/// Fits the continuous core-periphery model. See the module
/// documentation.
pub fn continuous_core_periphery(
    g: &Graph,
    weighted: bool,
    max_iter: usize,
) -> Result<ContinuousCorePeriphery> {
    let n = g.node_count();
    let obj = Objective::new(g, weighted)?;
    if n < 2 || !(obj.sxx > 0.0) {
        return Err(Error::InvalidArgument(
            "a core-periphery fit needs at least two nodes and pairs with different values \
             (some tied and some not)"
                .into(),
        ));
    }
    // Start: total value plus 1, as t = sqrt(c).
    let mut start = vec![1.0; n];
    obj.for_pairs(|i, j, a| {
        start[i] += a;
        start[j] += a;
    });
    let mut t: Vec<f64> = start.iter().map(|&x| x.sqrt()).collect();
    let norm = t.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in t.iter_mut() {
        *x /= norm;
    }
    // Objective in t (to maximize) and its gradient.
    let f = |t: &[f64]| -> (f64, Vec<f64>) {
        let c: Vec<f64> = t.iter().map(|x| x * x).collect();
        let (r, gc) = obj.eval(&c);
        let gt = gc.iter().zip(t).map(|(g, x)| 2.0 * x * g).collect();
        (r, gt)
    };
    let (mut fx, mut gx) = f(&t);
    if fx.is_nan() {
        return Err(Error::InvalidArgument(
            "the core-periphery correlation is undefined at the starting point".into(),
        ));
    }
    // Inverse Hessian approximation of -f (we minimize -r).
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        h[i * n + i] = 1.0;
    }
    let mut iterations = 0;
    while iterations < max_iter {
        let gmax = gx.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        if gmax < 1e-9 {
            break;
        }
        iterations += 1;
        // Ascent direction p = H g.
        let mut p: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| h[i * n + j] * gx[j]).sum())
            .collect();
        let mut slope: f64 = p.iter().zip(&gx).map(|(a, b)| a * b).sum();
        if !(slope > 0.0) {
            // Not an ascent direction: reset to steepest ascent.
            for i in 0..n {
                for j in 0..n {
                    h[i * n + j] = if i == j { 1.0 } else { 0.0 };
                }
            }
            p.clone_from(&gx);
            slope = gx.iter().map(|x| x * x).sum();
        }
        let mut step = 1.0;
        let mut accepted = None;
        for _ in 0..60 {
            let cand: Vec<f64> = t.iter().zip(&p).map(|(x, d)| x + step * d).collect();
            let (fc, gc) = f(&cand);
            if !fc.is_nan() && fc >= fx + 1e-4 * step * slope {
                accepted = Some((cand, fc, gc));
                break;
            }
            step *= 0.5;
        }
        let Some((cand, fc, gc)) = accepted else {
            break;
        };
        let s: Vec<f64> = cand.iter().zip(&t).map(|(a, b)| a - b).collect();
        // y is the change in the gradient of -f.
        let y: Vec<f64> = gx.iter().zip(&gc).map(|(a, b)| a - b).collect();
        let sy: f64 = s.iter().zip(&y).map(|(a, b)| a * b).sum();
        if sy > 1e-12 * s.iter().map(|x| x * x).sum::<f64>().sqrt() {
            // BFGS update of the inverse Hessian.
            let hy: Vec<f64> = (0..n)
                .map(|i| (0..n).map(|j| h[i * n + j] * y[j]).sum())
                .collect();
            let yhy: f64 = y.iter().zip(&hy).map(|(a, b)| a * b).sum();
            let rho = 1.0 / sy;
            for i in 0..n {
                for j in 0..n {
                    h[i * n + j] +=
                        rho * ((1.0 + rho * yhy) * s[i] * s[j] - hy[i] * s[j] - s[i] * hy[j]);
                }
            }
        }
        let improved = fc - fx;
        t = cand;
        fx = fc;
        gx = gc;
        if improved <= 1e-15 * fx.abs().max(1.0) {
            break;
        }
    }
    let mut c: Vec<f64> = t.iter().map(|x| x * x).collect();
    let norm = c.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in c.iter_mut() {
        *x /= norm;
    }
    let (r, _) = obj.eval(&c);
    Ok(ContinuousCorePeriphery {
        coreness: c,
        correlation: r,
        iterations,
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
    fn gradient_matches_finite_differences() {
        let g = und(5, &[(0, 1), (0, 2), (1, 2), (2, 3), (3, 4)]);
        let obj = Objective::new(&g, false).unwrap();
        let c = [0.9, 0.7, 0.8, 0.3, 0.1];
        let (r, grad) = obj.eval(&c);
        for i in 0..5 {
            let mut up = c;
            up[i] += 1e-6;
            let mut dn = c;
            dn[i] -= 1e-6;
            let fd = (obj.eval(&up).0 - obj.eval(&dn).0) / 2e-6;
            assert!((fd - grad[i]).abs() < 1e-6, "{i}: {fd} vs {}", grad[i]);
        }
        assert!(r > 0.0);
    }

    #[test]
    fn star_core_is_the_center() {
        let g = und(6, &[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]);
        let r = continuous_core_periphery(&g, false, 500).unwrap();
        let top = (0..6)
            .max_by(|&a, &b| r.coreness[a].total_cmp(&r.coreness[b]))
            .unwrap();
        assert_eq!(top, 0);
        assert!(r.coreness.iter().all(|&x| x >= 0.0));
        assert!((r.coreness.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-12);
        // Leaves are symmetric.
        assert!((r.coreness[1] - r.coreness[5]).abs() < 1e-6);
        assert!(continuous_core_periphery(&und(3, &[]), false, 10).is_err());
    }
}
