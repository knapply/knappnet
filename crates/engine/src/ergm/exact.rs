//! Exact likelihood by enumeration, for networks with few free dyads.
//!
//! With D free dyads there are 2^D networks. Visiting them in Gray code
//! order (Gray 1953, `gray1953pulse`; Knuth 2011, `knuth2011art`, Section
//! 7.2.1.1) changes one dyad per step, so each network's statistics follow
//! from the previous one's by a single change statistic. The enumeration
//! gives the normalizing constant c(theta) and so the exact log-likelihood,
//! the exact expectation and covariance of the statistics, and the exact
//! maximum likelihood estimate by Newton's method (Geyer and Thompson 1992,
//! `geyer1992constrained`, Section 2, where the exact case motivates the
//! Monte Carlo one). It is the reference against which Monte Carlo
//! estimates are checked, and is limited to D <= 22 (about four million
//! networks; roughly 7 undirected or 5 directed nodes without constraints).
//!
//! The maximum likelihood estimate does not exist when the observed
//! statistics lie on the boundary of the convex hull of attainable
//! statistics (Handcock 2003, `handcock2003assessing`; Rinaldo, Fienberg,
//! and Zhou 2009, `rinaldo2009geometry`); Newton's method then diverges and
//! an error says so.

use crate::error::{Error, Result};
use crate::matrix::Dense;

use super::linalg::{invert_symmetric, solve_spd};
use super::model::Model;

/// Largest number of free dyads enumerated.
pub const MAX_FREE_DYADS: usize = 22;

/// The statistics of every network on the free dyads.
#[derive(Debug, Clone)]
pub struct Enumeration {
    width: usize,
    count: usize,
    /// Row-major: network by network, `width` statistics each.
    stats: Vec<f64>,
    /// The observed statistics.
    observed: Vec<f64>,
}

fn no_mle() -> Error {
    Error::InvalidArgument(
        "the maximum likelihood estimate does not exist: the observed statistics lie on the \
         boundary of their convex hull (for example a count that is as small or as large as it \
         can be)"
            .into(),
    )
}

impl Enumeration {
    /// Enumerates every network that differs from the observed one only on
    /// free dyads.
    pub fn new(model: &Model) -> Result<Enumeration> {
        let dyads = model.free_dyads();
        let d = dyads.len();
        if d > MAX_FREE_DYADS {
            return Err(Error::InvalidArgument(format!(
                "exact enumeration needs at most {MAX_FREE_DYADS} free dyads, but there \
                 are {d}"
            )));
        }
        let p = model.width();
        let count = 1usize << d;
        let terms = model.terms();
        let observed = model.stats();
        // Start from the network with every free dyad empty.
        let mut y = model.state().clone();
        let mut current = observed.clone();
        let mut delta = vec![0.0; p];
        for &(i, j) in &dyads {
            if y.has(i as usize, j as usize) {
                terms.change(&y, i as usize, j as usize, &mut delta);
                for (c, dl) in current.iter_mut().zip(&delta) {
                    *c -= dl;
                }
                y.toggle(i as usize, j as usize);
            }
        }
        let mut stats = Vec::with_capacity(count * p);
        stats.extend_from_slice(&current);
        for k in 1..count {
            // Gray code: flip the dyad indexed by the number of trailing
            // zeros of k.
            let bit = k.trailing_zeros() as usize;
            let (i, j) = (dyads[bit].0 as usize, dyads[bit].1 as usize);
            terms.change(&y, i, j, &mut delta);
            let sign = if y.toggle(i, j) { 1.0 } else { -1.0 };
            for (c, dl) in current.iter_mut().zip(&delta) {
                *c += sign * dl;
            }
            stats.extend_from_slice(&current);
        }
        Ok(Enumeration {
            width: p,
            count,
            stats,
            observed,
        })
    }

    /// Number of networks enumerated.
    pub fn count(&self) -> usize {
        self.count
    }

    /// Statistics of network `k`.
    #[inline]
    pub fn row(&self, k: usize) -> &[f64] {
        &self.stats[k * self.width..(k + 1) * self.width]
    }

    /// The observed statistics.
    pub fn observed(&self) -> &[f64] {
        &self.observed
    }

    /// Log of the normalizing constant, log sum_y exp(theta' g(y)).
    pub fn log_constant(&self, theta: &[f64]) -> f64 {
        let mut max = f64::NEG_INFINITY;
        let mut linear = Vec::with_capacity(self.count);
        for k in 0..self.count {
            let v: f64 = self.row(k).iter().zip(theta).map(|(a, b)| a * b).sum();
            max = max.max(v);
            linear.push(v);
        }
        let sum: f64 = linear.iter().map(|v| (v - max).exp()).sum();
        max + sum.ln()
    }

    /// The exact log-likelihood of the observed network at `theta`.
    pub fn log_likelihood(&self, theta: &[f64]) -> Result<f64> {
        self.check(theta)?;
        let lin: f64 = self.observed.iter().zip(theta).map(|(a, b)| a * b).sum();
        Ok(lin - self.log_constant(theta))
    }

    /// The exact expectation and covariance of the statistics at `theta`.
    pub fn moments(&self, theta: &[f64]) -> Result<(Vec<f64>, Dense)> {
        self.check(theta)?;
        let p = self.width;
        let log_c = self.log_constant(theta);
        let mut mean = vec![0.0; p];
        let mut second = Dense::zeros(p, p);
        for k in 0..self.count {
            let g = self.row(k);
            let v: f64 = g.iter().zip(theta).map(|(a, b)| a * b).sum();
            let w = (v - log_c).exp();
            for a in 0..p {
                mean[a] += w * g[a];
                for b in 0..=a {
                    second.set(a, b, second.get(a, b) + w * g[a] * g[b]);
                }
            }
        }
        let mut cov = Dense::zeros(p, p);
        for a in 0..p {
            for b in 0..=a {
                let c = second.get(a, b) - mean[a] * mean[b];
                cov.set(a, b, c);
                cov.set(b, a, c);
            }
        }
        Ok((mean, cov))
    }

    /// The exact maximum likelihood estimate by Newton's method with step
    /// halving, started from `start`, with the inverse Fisher information
    /// as covariance and the log-likelihood at the estimate.
    pub fn mle(&self, start: &[f64], max_iter: usize) -> Result<ExactMle> {
        self.check(start)?;
        let mut theta = start.to_vec();
        let mut ll = self.log_likelihood(&theta)?;
        let mut converged = false;
        let mut iterations = 0;
        while iterations < max_iter {
            iterations += 1;
            let (mean, cov) = self.moments(&theta)?;
            let score: Vec<f64> = self
                .observed
                .iter()
                .zip(&mean)
                .map(|(o, m)| o - m)
                .collect();
            let step = match solve_spd(&cov, &score) {
                Some(step) => step,
                // As estimates run off toward a boundary the information
                // underflows, so tell that case apart from collinearity.
                None if self.observed_on_boundary() => return Err(no_mle()),
                None => {
                    return Err(Error::InvalidArgument(
                        "the Fisher information is singular: a statistic is constant over \
                         all networks, or two statistics are collinear"
                            .into(),
                    ))
                }
            };
            let mut scale = 1.0;
            let mut accepted = false;
            for _ in 0..30 {
                let trial: Vec<f64> = theta
                    .iter()
                    .zip(&step)
                    .map(|(t, s)| t + scale * s)
                    .collect();
                let ll_trial = self.log_likelihood(&trial)?;
                if ll_trial >= ll - 1e-12 {
                    let gain = ll_trial - ll;
                    theta = trial;
                    ll = ll_trial;
                    accepted = true;
                    if gain.abs() < 1e-12 && step.iter().all(|s| (scale * s).abs() < 1e-9) {
                        converged = true;
                    }
                    break;
                }
                scale *= 0.5;
            }
            if !accepted {
                converged = step.iter().all(|s| s.abs() < 1e-9);
                break;
            }
            if converged {
                break;
            }
            if theta.iter().any(|t| t.abs() > 50.0) {
                return Err(no_mle());
            }
        }
        if !converged {
            return Err(Error::InvalidArgument(
                "Newton's method did not converge to the exact maximum likelihood estimate".into(),
            ));
        }
        let (_, cov) = self.moments(&theta)?;
        let vcov = invert_symmetric(&cov)?;
        Ok(ExactMle {
            theta,
            cov: vcov,
            log_likelihood: ll,
            iterations,
        })
    }

    /// Whether the observed statistics lie on the boundary of the convex
    /// hull of the attainable ones, in which case the MLE does not exist
    /// (Rinaldo, Fienberg & Zhou 2009,
    /// `rinaldo2009geometry`). The test moves the observed point a small
    /// step away from the centroid of the distinct statistic vectors: an
    /// interior point stays inside the hull, a boundary point leaves it.
    fn observed_on_boundary(&self) -> bool {
        let mut seen = std::collections::HashSet::new();
        let mut rows: Vec<&[f64]> = Vec::new();
        for k in 0..self.count {
            let r = self.row(k);
            if seen.insert(r.iter().map(|v| v.to_bits()).collect::<Vec<u64>>()) {
                rows.push(r);
            }
        }
        let p = self.width;
        let mut centroid = vec![0.0; p];
        for r in &rows {
            for (c, v) in centroid.iter_mut().zip(r.iter()) {
                *c += v / rows.len() as f64;
            }
        }
        let x: Vec<f64> = self
            .observed
            .iter()
            .zip(&centroid)
            .map(|(o, c)| o + 1e-4 * (o - c))
            .collect();
        let mut points = Dense::zeros(rows.len(), p);
        for (i, r) in rows.iter().enumerate() {
            for (j, v) in r.iter().enumerate() {
                points.set(i, j, *v);
            }
        }
        !super::mcmle::in_convex_hull(&points, &x)
    }

    fn check(&self, theta: &[f64]) -> Result<()> {
        if theta.len() != self.width {
            return Err(Error::LengthMismatch {
                what: "coefficients",
                expected: self.width,
                found: theta.len(),
            });
        }
        if let Some(k) = theta.iter().position(|t| !t.is_finite()) {
            return Err(Error::InvalidArgument(format!(
                "coefficient {k} is not finite"
            )));
        }
        Ok(())
    }
}

/// An exact maximum likelihood estimate.
#[derive(Debug, Clone)]
pub struct ExactMle {
    /// The estimate.
    pub theta: Vec<f64>,
    /// Inverse Fisher information at the estimate.
    pub cov: Dense,
    /// Log-likelihood at the estimate.
    pub log_likelihood: f64,
    /// Newton iterations used.
    pub iterations: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ergm::state::NetState;
    use crate::ergm::terms::{Covariates, Term, Terms};

    fn model(n: usize, directed: bool, src: &[i64], tgt: &[i64], terms: Vec<Term>) -> Model {
        let y = NetState::from_edges(n, src, tgt, directed).unwrap();
        let t = Terms::new(terms, n, directed, Covariates::default()).unwrap();
        Model::new(y, t, None).unwrap()
    }

    /// An edges-only model is Bernoulli: the MLE is logit(density) and the
    /// log-likelihood is the binomial one.
    #[test]
    fn edges_model_is_bernoulli() {
        let m = model(5, false, &[0, 1, 2], &[1, 2, 3], vec![Term::Edges]);
        let e = Enumeration::new(&m).unwrap();
        assert_eq!(e.count(), 1 << 10);
        let fit = e.mle(&[0.0], 50).unwrap();
        let p: f64 = 3.0 / 10.0;
        assert!((fit.theta[0] - (p / (1.0 - p)).ln()).abs() < 1e-8);
        let ll = 3.0 * p.ln() + 7.0 * (1.0 - p).ln();
        assert!((fit.log_likelihood - ll).abs() < 1e-8);
        assert!((fit.cov.get(0, 0) - 1.0 / (10.0 * p * (1.0 - p))).abs() < 1e-6);
        // theta = 0: every network equally likely.
        assert!((e.log_likelihood(&[0.0]).unwrap() + 10.0 * std::f64::consts::LN_2).abs() < 1e-12);
    }

    #[test]
    fn enumeration_statistics_match_direct_counts() {
        let m = model(
            5,
            true,
            &[0, 1, 2, 0],
            &[1, 2, 0, 3],
            vec![
                Term::Edges,
                Term::Mutual,
                Term::TransitiveTriple,
                Term::CyclicTriple,
            ],
        );
        let e = Enumeration::new(&m).unwrap();
        // The observed network appears in the enumeration with its stats.
        let obs = e.observed().to_vec();
        assert!((0..e.count()).any(|k| e.row(k) == obs.as_slice()));
        // The last Gray code element is the network with only the highest
        // dyad set; its statistics are those of a single arc.
        assert_eq!(e.row(e.count() - 1)[0], 1.0);
        // Expected edges at theta = 0 is D / 2.
        let (mean, _) = e.moments(&[0.0, 0.0, 0.0, 0.0]).unwrap();
        assert!((mean[0] - 10.0).abs() < 1e-9);
        // Expected mutual pairs at theta = 0: 10 pairs, each mutual w.p. 1/4.
        assert!((mean[1] - 2.5).abs() < 1e-9);
    }

    #[test]
    fn boundary_statistics_have_no_mle() {
        // A complete graph: edges at its maximum.
        let m = model(
            4,
            false,
            &[0, 0, 0, 1, 1, 2],
            &[1, 2, 3, 2, 3, 3],
            vec![Term::Edges],
        );
        let e = Enumeration::new(&m).unwrap();
        let err = e.mle(&[0.0], 100).unwrap_err().to_string();
        assert!(err.contains("does not exist"), "{err}");
    }
}
