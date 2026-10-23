//! The continuous core-periphery model.
//!
//! Borgatti and Everett (1999, `borgatti1999models`, section 4) give every
//! node a coreness c_i and take the pattern matrix delta_ij = c_i c_j, so
//! that ties are expected between core nodes, less so between core and
//! peripheral nodes, and least between peripheral nodes. Two criteria fit
//! the coreness:
//!
//! - MINRES (the default, and UCINET's): minimize the sum of squared
//!   residuals sum_{i != j} (a_ij - c_i c_j)^2 over c >= 0, the minimum
//!   residual method of factor analysis (Comrey 1962, `comrey1962minimum`)
//!   with one factor and the diagonal left out, which Borgatti and Everett
//!   suggest and Boyd, Fitzgerald, Mahutga, and Smith (2010,
//!   `boyd2010computing`) describe for this model. Each coreness in turn is
//!   set to the value that minimizes the residuals given the others,
//!   c_i = sum_{j != i} (a_ij + a_ji) c_j / (2 sum_{j != i} c_j^2) (the two
//!   terms coincide for undirected data), clamped at 0; a sweep over all
//!   nodes never raises the residual sum, and sweeps continue until no
//!   coreness changes by more than `tol` (relative to the largest) or
//!   `max_iter` sweeps have been made. The start is c proportional to each
//!   node's total value. With the diagonal included the minimizer would be
//!   the principal eigenvector of the symmetrized data scaled by the square
//!   root of its eigenvalue; leaving the diagonal out changes the solution
//!   slightly, which is the point of MINRES. The result is scaled to unit
//!   Euclidean norm, as UCINET reports it. On data such as a star, whose
//!   leaves share no ties, the residual sum has no minimizer: it keeps
//!   falling as the hub's coreness grows and the leaves' shrink, so the
//!   sweeps stop at `max_iter`; the unit-norm coreness still converges (to
//!   the hub alone) and is reported.
//! - Correlation: maximize the Pearson correlation between the data and the
//!   pattern, the criterion Borgatti and Everett state. The authors
//!   maximize it numerically (with a Fletcher-Powell procedure); this module
//!   uses the quasi-Newton method of Broyden, Fletcher, Goldfarb, and Shanno
//!   with a backtracking line search that enforces the Armijo condition
//!   (Nocedal and Wright 2006, `nocedal2006numerical`, chapters 3 and 6),
//!   using the exact gradient of the correlation. The correlation is
//!   invariant to scaling c, so the result is scaled to unit norm.
//!   Caution: on many networks the supremum of the correlation is only
//!   approached as all coreness values become equal, where the pattern
//!   c_i c_j = (c + d_i)(c + d_j) reduces to the additive pattern d_i + d_j
//!   up to scale; the returned coreness is then nearly uniform and carries
//!   little information, while the correlation is that of the best additive
//!   fit. MINRES does not have this degeneracy.
//!
//! Conventions:
//! - Pairs and values as in [`super::discrete`]: pairs of distinct nodes,
//!   each unordered pair once for undirected networks and both ordered
//!   pairs for directed networks; binary values, or total weights.
//! - Coreness is constrained to be non-negative (by clamping for MINRES,
//!   and by writing c_i = t_i^2 and optimizing over t for the correlation).
//!   Borgatti and Everett do not state a sign constraint, but without one
//!   the product pattern also fits two separate dense groups (with coreness
//!   of opposite signs), which is not a core-periphery pattern.
//! - The correlation fit has several local maxima in general (and the
//!   substitution c = t^2 adds stationary points where some t_i = 0), so
//!   the search runs from `starts` points and keeps the best: the first
//!   start is c proportional to each node's total value plus 1 (degree
//!   plus 1 for binary data), and the others draw t_i uniformly from
//!   (0.05, 1) using [`crate::rng::Stream`] with the given seed, so the
//!   result is reproducible for a seed whatever the thread count. Each
//!   search stops when the largest gradient component is below 1e-9 or no
//!   step improves the correlation, after at most `max_iter` iterations.
//!   Among starts with the same best correlation the earliest wins. Cost
//!   O(n^2) per iteration and O(n^2) memory per thread. MINRES is
//!   deterministic and uses O(n^2) memory and O(n^2) time per sweep.
//! - Both fits report the correlation between the data and the final
//!   pattern, and the residual sum of squares of the unit-norm pattern
//!   scaled to the data (`residual` is the MINRES objective at the
//!   minimizer, before scaling to unit norm, for the MINRES fit, and NaN
//!   for the correlation fit).
//! - A network whose values are all equal is rejected.

use super::PairRows;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;
use crate::rng::Stream;

/// The criterion fitted by [`continuous_core_periphery`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ContinuousFit {
    /// Minimum residual sum of squares (UCINET).
    #[default]
    Minres,
    /// Maximum correlation with the pattern (Borgatti and Everett 1999).
    Correlation,
}

impl ContinuousFit {
    /// Decodes the C interface codes: 0 = MINRES, 1 = correlation.
    pub fn from_code(code: i32) -> Result<ContinuousFit> {
        match code {
            0 => Ok(ContinuousFit::Minres),
            1 => Ok(ContinuousFit::Correlation),
            _ => Err(Error::InvalidArgument(format!(
                "unknown continuous core-periphery fit code {code}"
            ))),
        }
    }
}

/// Options for [`continuous_core_periphery`].
#[derive(Debug, Clone, Copy)]
pub struct ContinuousOptions {
    /// The criterion.
    pub fit: ContinuousFit,
    /// Use total weights between nodes as values instead of binary ties.
    pub weighted: bool,
    /// Iteration limit per start (sweeps, for MINRES).
    pub max_iter: usize,
    /// Convergence tolerance of MINRES: largest change in coreness,
    /// relative to the largest coreness, below which sweeps stop.
    pub tol: f64,
    /// Number of starts (the degree start plus `starts - 1` random ones).
    pub starts: usize,
    /// Seed for the random starts.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for ContinuousOptions {
    fn default() -> Self {
        ContinuousOptions {
            fit: ContinuousFit::Minres,
            weighted: false,
            max_iter: 1000,
            tol: 1e-10,
            starts: 20,
            seed: 1,
            threads: 0,
        }
    }
}

/// Result of [`continuous_core_periphery`].
#[derive(Debug, Clone, PartialEq)]
pub struct ContinuousCorePeriphery {
    /// Coreness of each node, non-negative, with unit Euclidean norm.
    pub coreness: Vec<f64>,
    /// Correlation between the data and the pattern c_i c_j.
    pub correlation: f64,
    /// Residual sum of squares sum_{i != j} (a_ij - c_i c_j)^2 at the
    /// MINRES minimizer (before scaling to unit norm); NaN for the
    /// correlation fit.
    pub residual: f64,
    /// Iterations (sweeps) made.
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
        let mut sxx = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i == j || (!directed && j < i) {
                    continue;
                }
                let d = a[i * n + j] - mean_a;
                sxx += d * d;
            }
        }
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
    /// Sums of squares are taken about the means in a second pass, since
    /// the pattern values c_i c_j can all be nearly equal and the one-pass
    /// formula then loses most of its digits.
    fn eval(&self, c: &[f64]) -> (f64, Vec<f64>) {
        let mut sy = 0.0;
        self.for_pairs(|i, j, _| sy += c[i] * c[j]);
        let mean_y = sy / self.pairs;
        let (mut vyy, mut sxy) = (0.0, 0.0);
        self.for_pairs(|i, j, a| {
            let dy = c[i] * c[j] - mean_y;
            vyy += dy * dy;
            sxy += (a - self.mean_a) * dy;
        });
        let mut grad = vec![0.0; self.n];
        if vyy.is_nan() || vyy <= 0.0 {
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
    opts: ContinuousOptions,
) -> Result<ContinuousCorePeriphery> {
    let n = g.node_count();
    if opts.starts == 0 {
        return Err(Error::InvalidArgument(
            "the number of starts must be at least 1".into(),
        ));
    }
    let obj = Objective::new(g, opts.weighted)?;
    if n < 2 || obj.sxx.is_nan() || obj.sxx <= 0.0 {
        return Err(Error::InvalidArgument(
            "a core-periphery fit needs at least two nodes and pairs with different values \
             (some tied and some not)"
                .into(),
        ));
    }
    if opts.fit == ContinuousFit::Minres {
        return minres(&obj, opts);
    }
    let parts = parallel::for_each_item(opts.starts, opts.threads, Vec::new, |acc, s| {
        let t = if s == 0 {
            // Total value plus 1, as t = sqrt(c).
            let mut start = vec![1.0; n];
            obj.for_pairs(|i, j, a| {
                start[i] += a;
                start[j] += a;
            });
            start.iter().map(|&x| x.sqrt()).collect()
        } else {
            let mut rng = Stream::new(opts.seed, s as u64);
            (0..n).map(|_| 0.05 + 0.95 * rng.uniform()).collect()
        };
        acc.push((s, ascend(&obj, t, opts.max_iter)));
    });
    let mut runs: Vec<(usize, Option<Run>)> = parts.into_iter().flatten().collect();
    runs.sort_by_key(|r| r.0);
    let mut best: Option<Run> = None;
    for (_, run) in runs {
        if let Some(run) = run {
            if best.as_ref().map_or(true, |b| run.1 > b.1) {
                best = Some(run);
            }
        }
    }
    let Some((t, _, iterations)) = best else {
        return Err(Error::InvalidArgument(
            "the core-periphery correlation is undefined at every starting point".into(),
        ));
    };
    let mut c: Vec<f64> = t.iter().map(|x| x * x).collect();
    let norm = c.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in c.iter_mut() {
        *x /= norm;
    }
    let (r, _) = obj.eval(&c);
    Ok(ContinuousCorePeriphery {
        coreness: c,
        correlation: r,
        residual: f64::NAN,
        iterations,
    })
}

/// MINRES: coordinate descent on sum_{i != j} (a_ij - c_i c_j)^2 with
/// c >= 0. See the module documentation.
fn minres(obj: &Objective, opts: ContinuousOptions) -> Result<ContinuousCorePeriphery> {
    if !(opts.tol.is_finite() && opts.tol >= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "tol must be a finite non-negative number, not {}",
            opts.tol
        )));
    }
    let n = obj.n;
    // Symmetrized values s_ij = a_ij + a_ji (2 a_ij for undirected data),
    // so that the update reads c_i = sum_j s_ij c_j / (2 sum_j c_j^2).
    let sym = |i: usize, j: usize| -> f64 {
        if obj.directed {
            obj.a[i * n + j] + obj.a[j * n + i]
        } else {
            2.0 * obj.a[i * n + j]
        }
    };
    // Start: total value of each node (its row and column sums).
    let mut c: Vec<f64> = (0..n)
        .map(|i| (0..n).filter(|&j| j != i).map(|j| sym(i, j)).sum::<f64>())
        .collect();
    let scale = c.iter().cloned().fold(0.0f64, f64::max);
    if scale > 0.0 {
        for x in c.iter_mut() {
            *x /= scale;
        }
    }
    let mut sum_sq: f64 = c.iter().map(|x| x * x).sum();
    let mut iterations = 0;
    while iterations < opts.max_iter {
        iterations += 1;
        let mut change = 0.0f64;
        for i in 0..n {
            let rest = sum_sq - c[i] * c[i];
            let new = if rest > 0.0 {
                let num: f64 = (0..n).filter(|&j| j != i).map(|j| sym(i, j) * c[j]).sum();
                (num / (2.0 * rest)).max(0.0)
            } else {
                0.0
            };
            change = change.max((new - c[i]).abs());
            sum_sq += new * new - c[i] * c[i];
            c[i] = new;
        }
        let top = c.iter().cloned().fold(0.0f64, f64::max);
        if change <= opts.tol * top.max(f64::MIN_POSITIVE) {
            break;
        }
    }
    // Residual sum of squares over the pairs in the fit.
    let mut residual = 0.0;
    obj.for_pairs(|i, j, a| {
        let d = a - c[i] * c[j];
        residual += d * d;
    });
    let (r, _) = obj.eval(&c);
    let norm = c.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 0.0 {
        for x in c.iter_mut() {
            *x /= norm;
        }
    }
    Ok(ContinuousCorePeriphery {
        coreness: c,
        correlation: r,
        residual,
        iterations,
    })
}

/// Final point of one start, its correlation, and the iterations made.
type Run = (Vec<f64>, f64, usize);

/// BFGS ascent of the correlation over t (c = t^2) from `t`, returning the
/// final point, its correlation, and the iterations made, or None when the
/// correlation is undefined at the start (a constant pattern).
fn ascend(obj: &Objective, mut t: Vec<f64>, max_iter: usize) -> Option<Run> {
    let n = obj.n;
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
        return None;
    }
    // Inverse Hessian approximation of -f (we minimize -r).
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        h[i * n + i] = 1.0;
    }
    let mut fresh = true;
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
        if slope.is_nan() || slope <= 0.0 {
            // Not an ascent direction: reset to steepest ascent.
            reset(&mut h, n);
            fresh = true;
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
            // No step along p improves the correlation. A stale inverse
            // Hessian can make p useless while the gradient is still
            // large, so try once more from steepest ascent.
            if fresh {
                break;
            }
            reset(&mut h, n);
            fresh = true;
            continue;
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
            if fresh {
                break;
            }
            reset(&mut h, n);
            fresh = true;
            continue;
        }
        fresh = false;
    }
    Some((t, fx, iterations))
}

fn reset(h: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            h[i * n + j] = if i == j { 1.0 } else { 0.0 };
        }
    }
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
        for fit in [ContinuousFit::Minres, ContinuousFit::Correlation] {
            let opts = ContinuousOptions {
                fit,
                max_iter: 500,
                ..Default::default()
            };
            let r = continuous_core_periphery(&g, opts).unwrap();
            let top = (0..6)
                .max_by(|&a, &b| r.coreness[a].total_cmp(&r.coreness[b]))
                .unwrap();
            assert_eq!(top, 0, "{fit:?}");
            assert!(r.coreness.iter().all(|&x| x >= 0.0));
            assert!((r.coreness.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-12);
            // Leaves are symmetric.
            assert!(
                (r.coreness[1] - r.coreness[5]).abs() < 1e-6,
                "{fit:?} {r:?}"
            );
            assert!(continuous_core_periphery(&und(3, &[]), opts).is_err());
        }
    }

    #[test]
    fn starts_do_not_depend_on_threads() {
        let g = und(
            8,
            &[
                (0, 1),
                (0, 2),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 7),
                (3, 7),
            ],
        );
        let one = continuous_core_periphery(&g, ContinuousOptions::default()).unwrap();
        let opts = ContinuousOptions {
            threads: 3,
            ..Default::default()
        };
        let three = continuous_core_periphery(&g, opts).unwrap();
        assert_eq!(one, three);
        // A regular network has a constant degree start; random starts
        // still give an answer.
        let cycle = und(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
        let opts = ContinuousOptions {
            fit: ContinuousFit::Correlation,
            ..Default::default()
        };
        assert!(continuous_core_periphery(&cycle, opts).is_ok());
    }

    #[test]
    fn minres_is_a_stationary_point_of_the_residuals() {
        // Two triangles joined by a bridge, plus a pendant.
        let g = und(
            7,
            &[
                (0, 1),
                (0, 2),
                (1, 2),
                (2, 3),
                (3, 4),
                (3, 5),
                (4, 5),
                (5, 6),
            ],
        );
        let r = continuous_core_periphery(&g, ContinuousOptions::default()).unwrap();
        let obj = Objective::new(&g, false).unwrap();
        // Rescale the unit-norm coreness to the minimizer and check the
        // first-order condition of each coordinate.
        let scale = {
            let (mut num, mut den) = (0.0, 0.0);
            obj.for_pairs(|i, j, a| {
                let y = r.coreness[i] * r.coreness[j];
                num += a * y;
                den += y * y;
            });
            (num / den).sqrt()
        };
        let c: Vec<f64> = r.coreness.iter().map(|x| x * scale).collect();
        let mut residual = 0.0;
        obj.for_pairs(|i, j, a| residual += (a - c[i] * c[j]).powi(2));
        assert!(
            (residual - r.residual).abs() < 1e-6,
            "{residual} vs {}",
            r.residual
        );
        for i in 0..7 {
            let grad: f64 = (0..7)
                .filter(|&j| j != i)
                .map(|j| -2.0 * (obj.a[i * 7 + j] - c[i] * c[j]) * c[j])
                .sum();
            assert!(
                grad.abs() < 1e-6 || (c[i] == 0.0 && grad > 0.0),
                "{i}: {grad}"
            );
        }
        // The hub of a star has the highest coreness and the fit is exact
        // up to the leaf-leaf pairs.
        let star = und(5, &[(0, 1), (0, 2), (0, 3), (0, 4)]);
        let s = continuous_core_periphery(&star, ContinuousOptions::default()).unwrap();
        assert!(s.coreness[0] > s.coreness[1]);
        assert!(s.correlation > 0.9);
    }
}
