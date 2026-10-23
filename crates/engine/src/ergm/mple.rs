//! Maximum pseudolikelihood estimation (Strauss and Ikeda 1990,
//! `strauss1990pseudolikelihood`; Wasserman and Pattison 1996,
//! `wasserman1996logit`).
//!
//! The pseudolikelihood treats the conditional probabilities of the dyads
//! as if they were independent:
//!
//!     logit P(Y_ij = 1 | rest) = theta' delta_ij(y),
//!
//! so theta is the logistic regression of the observed dyad states on the
//! change statistics. When every term is dyad independent the conditional
//! probabilities are the marginal ones and the pseudolikelihood is the
//! likelihood, so the estimate is the exact maximum likelihood estimate.
//! Otherwise the estimate is consistent in some settings but its logistic
//! standard errors are not valid, because the dyads are dependent (van
//! Duijn, Gile, and Handcock 2009, `vanduijn2009framework`); they are
//! reported for completeness and labeled as such by callers.

use crate::error::Result;
use crate::matrix::Dense;

use super::logit::{logit, LogitFit, LogitOptions};
use super::model::Model;

/// A maximum pseudolikelihood fit.
#[derive(Debug, Clone)]
pub struct MpleFit {
    /// The logistic fit on the change statistics.
    pub fit: LogitFit,
    /// Number of free dyads (observations).
    pub dyads: usize,
    /// Whether every term is dyad independent, so the estimate is the
    /// maximum likelihood estimate and its standard errors are valid.
    pub dyad_independent: bool,
}

/// Fits the maximum pseudolikelihood estimate of `model`.
pub fn mple(model: &Model, opts: LogitOptions) -> Result<MpleFit> {
    let (x, y, dyads) = model.change_matrix();
    let fit = logit(&x, &y, None, opts)?;
    Ok(MpleFit {
        fit,
        dyads: dyads.len(),
        dyad_independent: model.terms().dyad_independent(),
    })
}

/// Conditional tie probabilities P(Y_ij = 1 | rest of y) at `theta` for
/// every free dyad of the observed network, in canonical order.
pub fn conditional_probabilities(model: &Model, theta: &[f64]) -> (Vec<f64>, Vec<(u32, u32)>) {
    let (x, _, dyads) = model.change_matrix();
    let p = probabilities_from_design(&x, theta);
    (p, dyads)
}

fn probabilities_from_design(x: &Dense, theta: &[f64]) -> Vec<f64> {
    (0..x.rows())
        .map(|r| {
            let eta: f64 = x.row(r).iter().zip(theta).map(|(a, b)| a * b).sum();
            1.0 / (1.0 + (-eta).exp())
        })
        .collect()
}
