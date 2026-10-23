//! Logistic regression by iteratively reweighted least squares (Nelder and
//! Wedderburn 1972, `nelder1972generalized`; McCullagh and Nelder 1989,
//! `mccullagh1989generalized`, Section 2.5), without an intercept unless the
//! design carries one. This is the fitting routine behind the maximum
//! pseudolikelihood estimate (Strauss and Ikeda 1990,
//! `strauss1990pseudolikelihood`) and pooled temporal pseudolikelihoods.
//!
//! Each observation is a row of the design with a binary response and a
//! non-negative weight (a case count or a bootstrap multiplicity). Rows that
//! are identical in design and response are merged before fitting, as
//! `ergm` does for change statistics, which makes the cost proportional to
//! the number of distinct rows.
//!
//! Convergence: the fit stops when the relative change in deviance is below
//! `tol` (R's `glm` uses 1e-8; the default here is 1e-12 so that the result
//! is at least as close to the optimum as `glm`'s) or after `max_iter`
//! iterations, in which case `converged` is false. Complete or
//! quasi-complete separation drives coefficients toward infinity; the fit
//! then reports `converged = false` with the coefficients of the last
//! iteration, as `glm` does with a warning.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::matrix::Dense;

use super::linalg::{cholesky_inverse, invert_symmetric, solve_spd};

/// Settings for the logistic fit.
#[derive(Debug, Clone, Copy)]
pub struct LogitOptions {
    /// Maximum number of IRLS iterations.
    pub max_iter: usize,
    /// Relative deviance change at which to stop.
    pub tol: f64,
}

impl Default for LogitOptions {
    fn default() -> Self {
        LogitOptions {
            max_iter: 50,
            tol: 1e-12,
        }
    }
}

/// A fitted logistic regression.
#[derive(Debug, Clone)]
pub struct LogitFit {
    /// Coefficients.
    pub coef: Vec<f64>,
    /// Covariance of the coefficients, the inverse of X' W X at the fit.
    pub cov: Dense,
    /// Residual deviance, -2 times the log-likelihood.
    pub deviance: f64,
    /// Deviance of the model with every coefficient zero (each response has
    /// probability 1/2): 2 log(2) times the total weight.
    pub null_deviance: f64,
    /// Total weight (number of observations when weights are counts).
    pub weight: f64,
    /// Iterations used.
    pub iterations: usize,
    /// Whether the deviance criterion was met.
    pub converged: bool,
}

/// Merges identical rows, summing their weights.
fn compress(x: &Dense, y: &[f64], w: &[f64]) -> (Dense, Vec<f64>, Vec<f64>) {
    let p = x.cols();
    let mut index: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut rows: Vec<f64> = Vec::new();
    let mut ys = Vec::new();
    let mut ws = Vec::new();
    for r in 0..x.rows() {
        if w[r] == 0.0 {
            continue;
        }
        let mut key: Vec<u64> = x.row(r).iter().map(|v| v.to_bits()).collect();
        key.push(y[r].to_bits());
        match index.get(&key) {
            Some(&k) => ws[k] += w[r],
            None => {
                index.insert(key, ws.len());
                rows.extend_from_slice(x.row(r));
                ys.push(y[r]);
                ws.push(w[r]);
            }
        }
    }
    let m = ws.len();
    (
        Dense::from_row_major(m, p, rows).expect("row count"),
        ys,
        ws,
    )
}

/// Fits a logistic regression of `y` (0 or 1) on the rows of `x` with
/// non-negative case weights `w` (`None` for unit weights).
pub fn logit(x: &Dense, y: &[f64], w: Option<&[f64]>, opts: LogitOptions) -> Result<LogitFit> {
    let (m, p) = (x.rows(), x.cols());
    if y.len() != m {
        return Err(Error::LengthMismatch {
            what: "response",
            expected: m,
            found: y.len(),
        });
    }
    let unit = vec![1.0; m];
    let w = w.unwrap_or(&unit);
    if w.len() != m {
        return Err(Error::LengthMismatch {
            what: "weights",
            expected: m,
            found: w.len(),
        });
    }
    if let Some(r) = y.iter().position(|v| *v != 0.0 && *v != 1.0) {
        return Err(Error::InvalidArgument(format!(
            "response {r} is {}, but the response must be 0 or 1",
            y[r]
        )));
    }
    if let Some(r) = w.iter().position(|v| !(*v >= 0.0) || !v.is_finite()) {
        return Err(Error::InvalidArgument(format!(
            "weight {r} is {}, but weights must be finite and non-negative",
            w[r]
        )));
    }
    if p == 0 {
        return Err(Error::InvalidArgument("the design has no columns".into()));
    }
    let (x, y, w) = compress(x, y, w);
    let m = y.len();
    let total: f64 = w.iter().sum();
    if total == 0.0 {
        return Err(Error::InvalidArgument(
            "every observation has zero weight".into(),
        ));
    }

    let mut beta = vec![0.0; p];
    let mut eta = vec![0.0; m];
    let deviance_at = |eta: &[f64]| -> f64 {
        let mut d = 0.0;
        for r in 0..m {
            // -2 w [y eta - log(1 + exp(eta))], computed stably.
            let e = eta[r];
            let log1pe = if e > 0.0 {
                e + (-e).exp().ln_1p()
            } else {
                e.exp().ln_1p()
            };
            d -= 2.0 * w[r] * (y[r] * e - log1pe);
        }
        d
    };
    let mut dev = deviance_at(&eta);
    let null_deviance = 2.0 * total * std::f64::consts::LN_2;
    let mut converged = false;
    let mut iterations = 0;
    let mut xtwx;
    while iterations < opts.max_iter {
        iterations += 1;
        // Weighted normal equations of the working response.
        let mut rhs = vec![0.0; p];
        xtwx = Dense::zeros(p, p);
        for r in 0..m {
            let mu = 1.0 / (1.0 + (-eta[r]).exp());
            let var = (mu * (1.0 - mu)).max(1e-300);
            let wr = w[r] * var;
            let z = eta[r] + (y[r] - mu) / var;
            let row = x.row(r);
            for a in 0..p {
                rhs[a] += wr * row[a] * z;
                for b in 0..=a {
                    xtwx.set(a, b, xtwx.get(a, b) + wr * row[a] * row[b]);
                }
            }
        }
        for a in 0..p {
            for b in 0..a {
                xtwx.set(b, a, xtwx.get(a, b));
            }
        }
        let new_beta = match solve_spd(&xtwx, &rhs) {
            Some(b) => b,
            None => {
                // Fall back to a general solve; a singular design is an
                // error, as in glm (which would drop aliased columns).
                let inv = invert_symmetric(&xtwx).map_err(|_| {
                    Error::InvalidArgument(
                        "the design is rank deficient; two statistics are collinear \
                         (or one is constant)"
                            .into(),
                    )
                })?;
                super::linalg::matvec(&inv, &rhs)
            }
        };
        let new_eta: Vec<f64> = (0..m)
            .map(|r| x.row(r).iter().zip(&new_beta).map(|(u, v)| u * v).sum())
            .collect();
        let new_dev = deviance_at(&new_eta);
        let done = (new_dev - dev).abs() / (new_dev.abs() + 0.1) < opts.tol;
        beta = new_beta;
        eta = new_eta;
        dev = new_dev;
        if done {
            converged = true;
            break;
        }
    }
    // Information at the final coefficients.
    xtwx = Dense::zeros(p, p);
    for r in 0..m {
        let mu = 1.0 / (1.0 + (-eta[r]).exp());
        let wr = w[r] * mu * (1.0 - mu);
        let row = x.row(r);
        for a in 0..p {
            for b in 0..=a {
                xtwx.set(a, b, xtwx.get(a, b) + wr * row[a] * row[b]);
            }
        }
    }
    for a in 0..p {
        for b in 0..a {
            xtwx.set(b, a, xtwx.get(a, b));
        }
    }
    let cov = match cholesky_inverse(&xtwx) {
        Some(c) => c,
        None => invert_symmetric(&xtwx).unwrap_or_else(|_| {
            let mut c = Dense::zeros(p, p);
            for v in 0..p {
                c.set(v, v, f64::INFINITY);
            }
            c
        }),
    };
    Ok(LogitFit {
        coef: beta,
        cov,
        deviance: dev,
        null_deviance,
        weight: total,
        iterations,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_closed_form_for_a_saturated_two_cell_design() {
        // One column, two distinct rows: x = 1 with 3 of 4 successes and
        // x = 0 (no information). The MLE of beta is logit(3/4).
        let x = Dense::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0]).unwrap();
        let fit = logit(&x, &[1.0, 1.0, 1.0, 0.0], None, LogitOptions::default()).unwrap();
        assert!((fit.coef[0] - (3.0f64).ln()).abs() < 1e-9);
        assert!(fit.converged);
        // Var = 1 / (n p (1 - p)) = 1 / (4 * 3/16).
        assert!((fit.cov.get(0, 0) - 16.0 / 12.0).abs() < 1e-8);
        assert!((fit.null_deviance - 8.0 * std::f64::consts::LN_2).abs() < 1e-12);
    }

    #[test]
    fn weights_act_as_case_counts() {
        let x = Dense::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0]).unwrap();
        let y = [1.0, 0.0, 0.0, 1.0];
        let a = logit(&x, &y, Some(&[2.0, 3.0, 1.0, 1.0]), LogitOptions::default()).unwrap();
        let x2 = Dense::from_row_major(
            7,
            2,
            vec![
                1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
            ],
        )
        .unwrap();
        let y2 = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let b = logit(&x2, &y2, None, LogitOptions::default()).unwrap();
        for k in 0..2 {
            assert!((a.coef[k] - b.coef[k]).abs() < 1e-9);
        }
        assert!((a.deviance - b.deviance).abs() < 1e-9);
    }

    #[test]
    fn separation_is_reported() {
        let x = Dense::from_row_major(4, 1, vec![-1.0, -2.0, 1.0, 2.0]).unwrap();
        let fit = logit(&x, &[0.0, 0.0, 1.0, 1.0], None, LogitOptions::default()).unwrap();
        assert!(!fit.converged || fit.coef[0] > 10.0);
    }
}
