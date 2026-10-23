//! Discrete-time temporal ERGMs by pooled pseudolikelihood.
//!
//! A temporal ERGM (Hanneke, Fu, and Xing 2010, `hanneke2010discrete`)
//! models each network in a sequence conditional on the previous one:
//!
//! ```text
//!     P(Y_t = y_t | y_{t-1}) = exp(theta' g(y_t, y_{t-1})) / c(theta, y_{t-1}),
//! ```
//!
//! with statistics that may depend on the previous network (memory and
//! delayed reciprocity terms, which are dyadic covariates built from
//! y_{t-1}). Desmarais and Cranmer (2012, `desmarais2012statistical`)
//! estimate theta by maximizing the pseudolikelihood pooled over the
//! transitions, and assess uncertainty by a bootstrap over the time steps:
//! transitions are resampled with replacement, the pooled fit repeated, and
//! confidence intervals read from the percentiles of the bootstrap
//! estimates. This is the method of the `btergm` package (Leifeld,
//! Cranmer, and Desmarais 2018, `leifeld2018temporal`), and its
//! pseudolikelihood point estimates agree with `btergm` exactly when the
//! terms agree.
//!
//! The caller assembles the pooled design: the change statistics of every
//! free dyad of every transition (see [`super::model::Model::change_matrix`])
//! stacked, with each row labeled by its transition. Bootstrap replicate r
//! draws from stream r of the seed ([`crate::rng::Stream`]), so results do
//! not depend on the number of threads.

use crate::error::{Error, Result};
use crate::matrix::Dense;
use crate::parallel;
use crate::rng::Stream;

use super::logit::{logit, LogitFit, LogitOptions};

/// A pooled pseudolikelihood fit with bootstrap replicates.
#[derive(Debug, Clone)]
pub struct TergmFit {
    /// The pooled fit over every transition.
    pub fit: LogitFit,
    /// Bootstrap estimates, replicate by replicate (`reps` rows of p
    /// values); NaN where a replicate's fit failed.
    pub boot: Vec<f64>,
    /// Number of replicates.
    pub reps: usize,
    /// Replicates whose fit converged.
    pub valid: usize,
}

/// Fits the pooled pseudolikelihood of `y` on `x`, with `period[r]` the
/// 0-based transition of row `r` among `periods` transitions, and `reps`
/// bootstrap replicates that resample transitions with replacement.
#[allow(clippy::too_many_arguments)]
pub fn pooled_bootstrap(
    x: &Dense,
    y: &[f64],
    period: &[u32],
    periods: usize,
    reps: usize,
    seed: u64,
    threads: usize,
    opts: LogitOptions,
) -> Result<TergmFit> {
    let m = x.rows();
    if period.len() != m {
        return Err(Error::LengthMismatch {
            what: "period labels",
            expected: m,
            found: period.len(),
        });
    }
    if periods == 0 {
        return Err(Error::InvalidArgument(
            "at least one transition is needed".into(),
        ));
    }
    if let Some(r) = period.iter().position(|&t| t as usize >= periods) {
        return Err(Error::InvalidArgument(format!(
            "row {r} is labeled with transition {}, but there are {periods}",
            period[r]
        )));
    }
    let fit = logit(x, y, None, opts)?;
    let p = x.cols();
    let results =
        parallel::for_each_item(reps, threads, Vec::<(usize, Vec<f64>)>::new, |acc, r| {
            let mut rng = Stream::new(seed, r as u64);
            let mut counts = vec![0.0; periods];
            for _ in 0..periods {
                counts[rng.below(periods as u64) as usize] += 1.0;
            }
            let w: Vec<f64> = period.iter().map(|&t| counts[t as usize]).collect();
            let est = match logit(x, y, Some(&w), opts) {
                Ok(f) if f.converged => f.coef,
                _ => vec![f64::NAN; p],
            };
            acc.push((r, est));
        });
    let mut all: Vec<(usize, Vec<f64>)> = results.into_iter().flatten().collect();
    all.sort_by_key(|(r, _)| *r);
    let mut boot = Vec::with_capacity(reps * p);
    let mut valid = 0;
    for (_, est) in all {
        if est.iter().all(|v| v.is_finite()) {
            valid += 1;
        }
        boot.extend(est);
    }
    Ok(TergmFit {
        fit,
        boot,
        reps,
        valid,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bootstrap_is_reproducible_and_centered() {
        // Two transitions with the same design: resampling only reweights
        // whole transitions, so every replicate equals the pooled fit.
        let x = Dense::from_row_major(
            8,
            2,
            vec![
                1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
            ],
        )
        .unwrap();
        let y = [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0];
        let period = [0, 0, 0, 0, 1, 1, 1, 1];
        let a = pooled_bootstrap(&x, &y, &period, 2, 20, 4, 0, LogitOptions::default()).unwrap();
        let b = pooled_bootstrap(&x, &y, &period, 2, 20, 4, 1, LogitOptions::default()).unwrap();
        assert_eq!(a.boot, b.boot);
        assert_eq!(a.valid, 20);
        for r in 0..20 {
            for k in 0..2 {
                assert!((a.boot[r * 2 + k] - a.fit.coef[k]).abs() < 1e-8);
            }
        }
    }
}
