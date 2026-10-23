//! Monte Carlo maximum likelihood estimation.
//!
//! The log-likelihood ratio between theta and a reference theta_0 is
//!
//! ```text
//!     l(theta) - l(theta_0) = (theta - theta_0)' g_obs
//!                             - log E_{theta_0} exp((theta - theta_0)' g(Y)),
//! ```
//!
//! and the expectation can be estimated from a sample drawn at theta_0
//! (Geyer and Thompson 1992, `geyer1992constrained`). Maximizing the
//! estimated ratio gives the Monte Carlo MLE. The estimate is reliable only
//! when theta is close to theta_0, so the fit proceeds in iterations, each
//! sampling at the current estimate.
//!
//! The iterations use the stepping algorithm of Hummel, Hunter, and
//! Handcock (2012, `hummel2012improving`):
//!
//! 1. Sample at theta_t and compute the sample mean g_bar and covariance S.
//! 2. Find the largest gamma in (0, 1] such that the pseudo-observation
//!    gamma g_obs + (1 - gamma) g_bar lies inside the convex hull of the
//!    sampled statistics (shrunk toward g_bar by a small margin, as
//!    `ergm`'s `MCMLE.steplength.margin` does). Membership in the hull is a
//!    linear program, solved here by a phase-one simplex.
//! 3. Update theta_t by one Newton step toward the pseudo-observation under
//!    the lognormal approximation of the ratio, theta + S^{-1} (target -
//!    g_bar), which is exact when the statistics are normal.
//! 4. Once gamma = 1, maximize the full Monte Carlo log-likelihood ratio by
//!    Newton's method with importance weights, sample again at the new
//!    estimate, and stop when the estimating equation g_bar = g_obs holds
//!    within Monte Carlo error (a Hotelling-type statistic of g_bar - g_obs
//!    against its Monte Carlo covariance below a chi-squared critical
//!    value supplied by the caller).
//!
//! Standard errors follow Hunter and Handcock (2006, `hunter2006inference`,
//! Section 3.3): the covariance of the estimate is I^{-1} + I^{-1} V I^{-1},
//! where I is the Fisher information, estimated by the sample covariance of
//! the statistics at the estimate, and V the Monte Carlo covariance of
//! their mean, estimated by batch means (Geyer 1992, `geyer1992practical`).
//! The second term is the contribution of Monte Carlo error, reported
//! separately so it can be shown as a percentage, as `ergm` does.
//!
//! The log-likelihood at the estimate is found by path sampling (Gelman
//! and Meng 1998, `gelman1998simulating`; Hunter and Handcock 2006, Section
//! 4): l(theta) - l(0) = theta' g_obs - integral_0^1 theta' E_{t theta}[g(Y)] dt,
//! the integral estimated by the midpoint rule with a sample at each step,
//! and l(0) = -D log 2 for D free dyads since every network is equally
//! likely when theta = 0.
//!
//! Degeneracy (Handcock 2003, `handcock2003assessing`; Schweinberger 2011,
//! `schweinberger2011instability`): if the sampler returns statistics with
//! no variation (it sits at the empty or the complete network, or another
//! fixed point) the fit stops with an error that says so, since no step
//! can be computed.

use crate::error::{Error, Result};
use crate::matrix::Dense;

use super::linalg::{invert_symmetric, matmul, quadratic_form, solve_spd};
use super::model::Model;
use super::sampler::{simulate, Sample, SampleOptions};

/// Settings for the Monte Carlo MLE.
#[derive(Debug, Clone)]
pub struct McmleOptions {
    /// Sampling settings for each iteration (`stream_offset` is managed
    /// here so every iteration uses fresh streams).
    pub sample: SampleOptions,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Hull shrinkage used when finding the step length, in (0, 1).
    pub margin: f64,
    /// Chi-squared critical value for the convergence test of the
    /// estimating equation, with p degrees of freedom.
    pub critical: f64,
}

impl Default for McmleOptions {
    fn default() -> Self {
        McmleOptions {
            sample: SampleOptions::default(),
            max_iter: 60,
            margin: 0.05,
            critical: 0.0,
        }
    }
}

/// A Monte Carlo maximum likelihood fit.
#[derive(Debug, Clone)]
pub struct McmleFit {
    /// The estimate.
    pub theta: Vec<f64>,
    /// Inverse Fisher information at the estimate.
    pub cov: Dense,
    /// Monte Carlo contribution to the covariance, I^{-1} V I^{-1}.
    pub mcmc_cov: Dense,
    /// Iterations used.
    pub iterations: usize,
    /// Whether the estimating equation test passed.
    pub converged: bool,
    /// The last step length.
    pub step_length: f64,
    /// The Hotelling-type statistic of the final sample.
    pub statistic: f64,
    /// The final sample, drawn at the estimate.
    pub sample: Sample,
}

/// Sample mean and covariance (divisor rows - 1).
fn moments(s: &Sample) -> (Vec<f64>, Dense) {
    let p = s.width;
    let rows = s.rows();
    let mean = s.mean();
    let mut cov = Dense::zeros(p, p);
    for r in 0..rows {
        let g = s.row(r);
        for a in 0..p {
            let da = g[a] - mean[a];
            for b in 0..=a {
                cov.set(a, b, cov.get(a, b) + da * (g[b] - mean[b]));
            }
        }
    }
    let denom = (rows.max(2) - 1) as f64;
    for a in 0..p {
        for b in 0..=a {
            let v = cov.get(a, b) / denom;
            cov.set(a, b, v);
            cov.set(b, a, v);
        }
    }
    (mean, cov)
}

/// Monte Carlo covariance of the overall sample mean, by batch means within
/// each chain, chains taken as independent.
pub fn mcmc_covariance(s: &Sample) -> Dense {
    let p = s.width;
    let mut total = Dense::zeros(p, p);
    let nsim = s.nsim;
    if nsim < 4 {
        // Too few draws for batches: treat every draw as independent.
        let (_, cov) = moments(s);
        let rows = s.rows().max(1) as f64;
        for a in 0..p {
            for c in 0..p {
                total.set(a, c, cov.get(a, c) / rows);
            }
        }
        return total;
    }
    let b = ((nsim as f64).sqrt().floor() as usize).max(2);
    let k = nsim / b;
    for c in 0..s.chains {
        let base = c * nsim;
        let mut means = vec![vec![0.0; p]; b];
        for (t, mean_t) in means.iter_mut().enumerate() {
            for r in 0..k {
                let g = s.row(base + t * k + r);
                for a in 0..p {
                    mean_t[a] += g[a] / k as f64;
                }
            }
        }
        let mut grand = vec![0.0; p];
        for m in &means {
            for a in 0..p {
                grand[a] += m[a] / b as f64;
            }
        }
        // Var(chain mean) = Var(batch mean) / b.
        for m in &means {
            for a in 0..p {
                for c2 in 0..p {
                    let v = (m[a] - grand[a]) * (m[c2] - grand[c2]) / ((b - 1) as f64 * b as f64);
                    total.set(a, c2, total.get(a, c2) + v);
                }
            }
        }
    }
    // The overall mean averages the chain means.
    let chains2 = (s.chains * s.chains) as f64;
    for a in 0..p {
        for c2 in 0..p {
            total.set(a, c2, total.get(a, c2) / chains2);
        }
    }
    total
}

/// Whether `x` lies in the convex hull of the rows of `points` (an m-by-p
/// matrix), by a phase-one simplex on the feasibility problem
/// sum lambda_i s_i = x, sum lambda_i = 1, lambda >= 0.
pub fn in_convex_hull(points: &Dense, x: &[f64]) -> bool {
    let (m, p) = (points.rows(), points.cols());
    if m == 0 {
        return false;
    }
    // Scale each coordinate to unit range for numerical stability; a
    // coordinate without variation must match exactly.
    let mut scale = vec![1.0; p];
    let mut shift = vec![0.0; p];
    for c in 0..p {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for r in 0..m {
            lo = lo.min(points.get(r, c));
            hi = hi.max(points.get(r, c));
        }
        shift[c] = lo;
        if hi > lo {
            scale[c] = hi - lo;
        } else if (x[c] - lo).abs() > 1e-9 * (1.0 + lo.abs()) {
            return false;
        }
    }
    let rows = p + 1;
    let cols = m + rows + 1; // lambda, artificials, right-hand side
    let mut t = vec![0.0; rows * cols];
    let idx = |r: usize, c: usize| r * cols + c;
    for r in 0..p {
        let mut b = (x[r] - shift[r]) / scale[r];
        let sign = if b < 0.0 { -1.0 } else { 1.0 };
        b *= sign;
        for i in 0..m {
            t[idx(r, i)] = sign * (points.get(i, r) - shift[r]) / scale[r];
        }
        t[idx(r, m + r)] = 1.0;
        t[idx(r, cols - 1)] = b;
    }
    for i in 0..m {
        t[idx(p, i)] = 1.0;
    }
    t[idx(p, m + p)] = 1.0;
    t[idx(p, cols - 1)] = 1.0;
    // Reduced costs of phase one: minimize the sum of artificials.
    let mut cost = vec![0.0; cols];
    for i in 0..m {
        cost[i] = -(0..rows).map(|r| t[idx(r, i)]).sum::<f64>();
    }
    cost[cols - 1] = -(0..rows).map(|r| t[idx(r, cols - 1)]).sum::<f64>();
    let mut basis: Vec<usize> = (m..m + rows).collect();
    let eps = 1e-10;
    let max_pivots = 50 * (m + rows);
    for _ in 0..max_pivots {
        // Entering column: most negative reduced cost (Dantzig's rule).
        let mut enter = None;
        let mut best = -eps;
        for (c, &cc) in cost.iter().enumerate().take(cols - 1) {
            if cc < best {
                best = cc;
                enter = Some(c);
            }
        }
        let Some(e) = enter else { break };
        // Leaving row by the ratio test.
        let mut leave = None;
        let mut ratio = f64::INFINITY;
        for r in 0..rows {
            let a = t[idx(r, e)];
            if a > eps {
                let q = t[idx(r, cols - 1)] / a;
                if q < ratio - 1e-14 {
                    ratio = q;
                    leave = Some(r);
                }
            }
        }
        let Some(l) = leave else { break };
        // Pivot.
        let pv = t[idx(l, e)];
        for c in 0..cols {
            t[idx(l, c)] /= pv;
        }
        for r in 0..rows {
            if r == l {
                continue;
            }
            let f = t[idx(r, e)];
            if f != 0.0 {
                for c in 0..cols {
                    t[idx(r, c)] -= f * t[idx(l, c)];
                }
            }
        }
        let f = cost[e];
        for c in 0..cols {
            cost[c] -= f * t[idx(l, c)];
        }
        basis[l] = e;
    }
    // Objective value is -cost[rhs]; feasible when it is (near) zero.
    (-cost[cols - 1]).abs() <= 1e-7
}

/// The largest gamma in [0, 1] such that gamma g_obs + (1 - gamma) g_bar is
/// in the convex hull of the sample shrunk toward g_bar by `margin`.
pub fn step_length(s: &Sample, g_obs: &[f64], g_bar: &[f64], margin: f64) -> f64 {
    let p = s.width;
    let rows = s.rows();
    // Thin very large samples: the hull of a subset is contained in the
    // hull of the whole, so the step length is conservative.
    let keep = 4000usize;
    let stride = (rows / keep).max(1);
    let m = (rows + stride - 1) / stride;
    let mut pts = Dense::zeros(m, p);
    for (r, row) in (0..rows).step_by(stride).enumerate() {
        let g = s.row(row);
        for c in 0..p {
            pts.set(r, c, g_bar[c] + (1.0 - margin) * (g[c] - g_bar[c]));
        }
    }
    let target = |gamma: f64| -> Vec<f64> {
        (0..p)
            .map(|c| gamma * g_obs[c] + (1.0 - gamma) * g_bar[c])
            .collect()
    };
    if in_convex_hull(&pts, &target(1.0)) {
        return 1.0;
    }
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..25 {
        let mid = 0.5 * (lo + hi);
        if in_convex_hull(&pts, &target(mid)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

fn degenerate(cov: &Dense) -> Option<usize> {
    (0..cov.rows()).find(|&a| cov.get(a, a) <= 0.0 || !cov.get(a, a).is_finite())
}

fn degeneracy_error(s: &Sample, stat: usize) -> Error {
    let first = s.row(0)[stat];
    let ties = s.row(0)[0];
    let d = if s.directed {
        s.n * s.n.saturating_sub(1)
    } else {
        s.n * s.n.saturating_sub(1) / 2
    };
    let where_ = if s.row(0).iter().all(|v| *v == 0.0) {
        " (every sampled network is empty)".to_string()
    } else if ties as usize == d {
        " (every sampled network is complete)".to_string()
    } else {
        String::new()
    };
    Error::InvalidArgument(format!(
        "the model is degenerate at the current coefficients: statistic {stat} is {first} \
         in every sampled network{where_}; the sampler cannot move, so no estimate exists. \
         Try a different model specification or starting values"
    ))
}

/// Newton's method on the importance-sampled log-likelihood ratio from the
/// sample's theta_0 (at which the sample was drawn) toward the maximizer.
fn importance_newton(s: &Sample, theta0: &[f64], g_obs: &[f64]) -> Result<Vec<f64>> {
    let p = s.width;
    let rows = s.rows();
    let mut theta = theta0.to_vec();
    let ratio = |theta: &[f64]| -> (f64, Vec<f64>, Dense, f64) {
        // Log ratio, weighted mean, weighted covariance, and effective
        // sample size of the weights.
        let d: Vec<f64> = theta.iter().zip(theta0).map(|(a, b)| a - b).collect();
        let mut lin = Vec::with_capacity(rows);
        let mut max = f64::NEG_INFINITY;
        for r in 0..rows {
            let v: f64 = s.row(r).iter().zip(&d).map(|(a, b)| a * b).sum();
            max = max.max(v);
            lin.push(v);
        }
        let mut w: Vec<f64> = lin.iter().map(|v| (v - max).exp()).collect();
        let sum: f64 = w.iter().sum();
        w.iter_mut().for_each(|x| *x /= sum);
        let log_mean = max + (sum / rows as f64).ln();
        let obs: f64 = g_obs.iter().zip(&d).map(|(a, b)| a * b).sum();
        let mut mean = vec![0.0; p];
        for (r, wr) in w.iter().enumerate() {
            for (m, g) in mean.iter_mut().zip(s.row(r)) {
                *m += wr * g;
            }
        }
        let mut cov = Dense::zeros(p, p);
        for (r, wr) in w.iter().enumerate() {
            let g = s.row(r);
            for a in 0..p {
                for b in 0..=a {
                    cov.set(
                        a,
                        b,
                        cov.get(a, b) + wr * (g[a] - mean[a]) * (g[b] - mean[b]),
                    );
                }
            }
        }
        for a in 0..p {
            for b in 0..a {
                cov.set(b, a, cov.get(a, b));
            }
        }
        let ess = 1.0 / w.iter().map(|x| x * x).sum::<f64>();
        (obs - log_mean, mean, cov, ess)
    };
    let (mut value, mut mean, mut cov, _) = ratio(&theta);
    for _ in 0..100 {
        let score: Vec<f64> = g_obs.iter().zip(&mean).map(|(o, m)| o - m).collect();
        let step = match solve_spd(&cov, &score) {
            Some(st) => st,
            None => return Err(degeneracy_error(s, degenerate(&cov).unwrap_or(0))),
        };
        if step.iter().all(|x| x.abs() < 1e-8) {
            break;
        }
        let mut scale = 1.0;
        let mut moved = false;
        for _ in 0..30 {
            let trial: Vec<f64> = theta
                .iter()
                .zip(&step)
                .map(|(t, st)| t + scale * st)
                .collect();
            let (v, m, c, ess) = ratio(&trial);
            // Keep the importance weights informative: at least a handful
            // of effective draws.
            if v >= value - 1e-12 && ess >= 10.0_f64.min(rows as f64) {
                theta = trial;
                value = v;
                mean = m;
                cov = c;
                moved = true;
                break;
            }
            scale *= 0.5;
        }
        if !moved {
            break;
        }
    }
    Ok(theta)
}

/// Fits the Monte Carlo MLE starting from `theta0` (normally the
/// pseudolikelihood estimate).
pub fn mcmle(model: &Model, theta0: &[f64], opts: &McmleOptions) -> Result<McmleFit> {
    let p = model.width();
    if theta0.len() != p {
        return Err(Error::LengthMismatch {
            what: "starting coefficients",
            expected: p,
            found: theta0.len(),
        });
    }
    if !(0.0..1.0).contains(&opts.margin) {
        return Err(Error::InvalidArgument(
            "the hull margin must be in [0, 1)".into(),
        ));
    }
    if opts.sample.nsim * opts.sample.chains.max(1) < 2 {
        return Err(Error::InvalidArgument(
            "Monte Carlo estimation needs at least two draws per iteration".into(),
        ));
    }
    let g_obs = model.stats();
    let mut theta = theta0.to_vec();
    let mut sopts = opts.sample.clone();
    sopts.keep_networks = false;
    let mut iterations = 0;
    let mut gamma_was_one = false;
    let mut last_gamma = 0.0;
    loop {
        if iterations >= opts.max_iter {
            // Return the current state, flagged as not converged.
            sopts.stream_offset = (iterations as u64 + 1) * sopts.chains as u64;
            let s = simulate(model, &theta, &sopts)?;
            let (mean, cov) = moments(&s);
            if let Some(k) = degenerate(&cov) {
                return Err(degeneracy_error(&s, k));
            }
            let (c, mc, t2) = covariances(&s, &mean, &cov, &g_obs)?;
            return Ok(McmleFit {
                theta,
                cov: c,
                mcmc_cov: mc,
                iterations,
                converged: false,
                step_length: last_gamma,
                statistic: t2,
                sample: s,
            });
        }
        iterations += 1;
        sopts.stream_offset = iterations as u64 * sopts.chains as u64;
        let s = simulate(model, &theta, &sopts)?;
        let (mean, cov) = moments(&s);
        if let Some(k) = degenerate(&cov) {
            return Err(degeneracy_error(&s, k));
        }
        let gamma = step_length(&s, &g_obs, &mean, opts.margin);
        last_gamma = gamma;
        if gamma >= 1.0 {
            if gamma_was_one {
                // The sample at the candidate estimate: test the estimating
                // equation.
                let (c, mc, t2) = covariances(&s, &mean, &cov, &g_obs)?;
                if t2 <= opts.critical {
                    return Ok(McmleFit {
                        theta,
                        cov: c,
                        mcmc_cov: mc,
                        iterations,
                        converged: true,
                        step_length: gamma,
                        statistic: t2,
                        sample: s,
                    });
                }
            }
            theta = importance_newton(&s, &theta, &g_obs)?;
            gamma_was_one = true;
        } else {
            gamma_was_one = false;
            // Lognormal approximation: one Newton step toward the
            // pseudo-observation.
            let target: Vec<f64> = (0..p)
                .map(|c| gamma * g_obs[c] + (1.0 - gamma) * mean[c])
                .collect();
            let diff: Vec<f64> = target.iter().zip(&mean).map(|(t, m)| t - m).collect();
            let step = solve_spd(&cov, &diff)
                .ok_or_else(|| degeneracy_error(&s, degenerate(&cov).unwrap_or(0)))?;
            for (t, st) in theta.iter_mut().zip(&step) {
                *t += st;
            }
        }
        if theta.iter().any(|t| !t.is_finite() || t.abs() > 1e3) {
            return Err(Error::InvalidArgument(
                "the coefficients diverged during estimation; the model is likely \
                 degenerate or the observed statistics lie on the boundary of their \
                 convex hull, so the maximum likelihood estimate does not exist"
                    .into(),
            ));
        }
    }
}

/// Inverse information, the Monte Carlo term, and the Hotelling-type
/// statistic of the estimating equation for a sample at the estimate.
fn covariances(
    s: &Sample,
    mean: &[f64],
    cov: &Dense,
    g_obs: &[f64],
) -> Result<(Dense, Dense, f64)> {
    let inv = invert_symmetric(cov)?;
    let v = mcmc_covariance(s);
    let mc = matmul(&matmul(&inv, &v), &inv);
    let diff: Vec<f64> = mean.iter().zip(g_obs).map(|(m, o)| m - o).collect();
    let t2 = match invert_symmetric(&v) {
        Ok(vinv) => quadratic_form(&vinv, &diff),
        Err(_) => f64::INFINITY,
    };
    Ok((inv, mc, t2))
}

/// Path-sampling estimate of l(theta) - l(0) with `steps` midpoint
/// evaluations, and its Monte Carlo standard error.
pub fn log_likelihood_ratio(
    model: &Model,
    theta: &[f64],
    steps: usize,
    opts: &SampleOptions,
) -> Result<(f64, f64)> {
    let p = model.width();
    if theta.len() != p {
        return Err(Error::LengthMismatch {
            what: "coefficients",
            expected: p,
            found: theta.len(),
        });
    }
    if steps == 0 {
        return Err(Error::InvalidArgument("steps must be at least 1".into()));
    }
    let mut sopts = opts.clone();
    sopts.keep_networks = false;
    // l(theta) - l(0) = theta' g_obs - [log c(theta) - log c(0)], and the
    // bracket is the integral of theta' E_{t theta}[g] over t in [0, 1].
    let g_obs = model.stats();
    let mut total: f64 = theta.iter().zip(&g_obs).map(|(a, b)| a * b).sum();
    let mut var = 0.0;
    for k in 0..steps {
        let t = (k as f64 + 0.5) / steps as f64;
        let at: Vec<f64> = theta.iter().map(|v| v * t).collect();
        sopts.stream_offset = opts.stream_offset + (k as u64 + 1) * sopts.chains as u64;
        let s = simulate(model, &at, &sopts)?;
        let mean = s.mean();
        let dot: f64 = theta.iter().zip(&mean).map(|(a, b)| a * b).sum();
        total -= dot / steps as f64;
        let v = mcmc_covariance(&s);
        var += quadratic_form(&v, theta) / (steps * steps) as f64;
    }
    Ok((total, var.sqrt()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ergm::exact::Enumeration;
    use crate::ergm::state::NetState;
    use crate::ergm::terms::{Covariates, Partners, Term, Terms};
    use crate::ergm::Proposal;

    #[test]
    fn hull_membership() {
        let pts =
            Dense::from_row_major(4, 2, vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap();
        assert!(in_convex_hull(&pts, &[0.5, 0.5]));
        assert!(in_convex_hull(&pts, &[1.0, 0.0]));
        assert!(!in_convex_hull(&pts, &[1.2, 0.5]));
        assert!(!in_convex_hull(&pts, &[0.5, -0.1]));
        // A flat sample: only the plane counts.
        let flat = Dense::from_row_major(3, 2, vec![0.0, 2.0, 1.0, 2.0, 3.0, 2.0]).unwrap();
        assert!(in_convex_hull(&flat, &[2.0, 2.0]));
        assert!(!in_convex_hull(&flat, &[2.0, 2.5]));
    }

    /// On a tiny network the Monte Carlo MLE must agree with the exact MLE
    /// within Monte Carlo error.
    #[test]
    fn recovers_the_exact_mle() {
        let n = 6;
        let y =
            NetState::from_edges(n, &[0, 0, 1, 1, 2, 3, 4], &[1, 2, 2, 3, 4, 4, 5], false).unwrap();
        let terms = Terms::new(
            vec![
                Term::Edges,
                Term::GwPartners {
                    decay: 0.5,
                    partners: Partners::Edgewise,
                },
            ],
            n,
            false,
            Covariates::default(),
        )
        .unwrap();
        let model = Model::new(y, terms, None).unwrap();
        let exact = Enumeration::new(&model)
            .unwrap()
            .mle(&[0.0, 0.0], 100)
            .unwrap();
        let opts = McmleOptions {
            sample: SampleOptions {
                burnin: 2000,
                interval: 20,
                nsim: 2000,
                chains: 2,
                proposal: Proposal::Tnt,
                seed: 11,
                ..Default::default()
            },
            max_iter: 40,
            margin: 0.05,
            critical: 5.99, // chi-squared(2), 0.95
        };
        let fit = mcmle(&model, &[0.0, 0.0], &opts).unwrap();
        for k in 0..2 {
            let se = fit.cov.get(k, k).sqrt();
            assert!(
                (fit.theta[k] - exact.theta[k]).abs() < 0.3 * se.max(0.1) + 0.05,
                "stat {k}: mcmle {} exact {} se {se}",
                fit.theta[k],
                exact.theta[k]
            );
        }
        assert!(fit.converged);
        // Path-sampled log-likelihood vs the exact one.
        let (llr, se) = log_likelihood_ratio(&model, &fit.theta, 16, &opts.sample).unwrap();
        let ll = llr - 15.0 * std::f64::consts::LN_2;
        let exact_ll = Enumeration::new(&model)
            .unwrap()
            .log_likelihood(&fit.theta)
            .unwrap();
        assert!(
            (ll - exact_ll).abs() < 4.0 * se + 0.05,
            "{ll} vs {exact_ll} (se {se})"
        );
    }

    #[test]
    fn degenerate_models_are_reported() {
        let n = 5;
        let y = NetState::from_edges(n, &[0], &[1], false).unwrap();
        let terms = Terms::new(
            vec![Term::Edges, Term::Triangle],
            n,
            false,
            Covariates::default(),
        )
        .unwrap();
        let model = Model::new(y, terms, None).unwrap();
        let opts = McmleOptions {
            sample: SampleOptions {
                burnin: 500,
                interval: 10,
                nsim: 200,
                seed: 1,
                ..Default::default()
            },
            ..Default::default()
        };
        // Huge positive triangle coefficient: the chain sits at the
        // complete network.
        let err = mcmle(&model, &[5.0, 20.0], &opts).unwrap_err();
        assert!(err.to_string().contains("degenerate"), "{err}");
    }
}
