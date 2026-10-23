//! Method-of-moments estimation by stochastic approximation.
//!
//! The estimate solves the moment equation E_theta S = s, where S is the
//! vector of statistics of simulated panels and s its observed value
//! (Snijders 2001, `snijders2001statistical`; Snijders 2017,
//! `snijders2017stochastic`). Statistics: for each effect k, the sum over
//! periods of s_k evaluated on the network at the end of the period; with
//! unconditional estimation also, for each period, the number of tie
//! variables that changed (the rate parameters' statistic). The equation
//! is solved by the Robbins-Monro algorithm (Robbins and Monro 1951,
//! `robbins1951stochastic`) in the three-phase form of Snijders (2001,
//! appendix) with the refinements RSiena documents in `Siena_algorithms`
//! (Ripley et al. 2026, `ripley2026manual`):
//!
//! 1. **Phase 1**: n1 = 7 + 3p simulations at the initial value theta_1
//!    estimate the derivative matrix D = dE S / dtheta by the score-function
//!    method (Schweinberger and Snijders 2007, `schweinberger2007markov`):
//!    D = sum_m cov(S_m, J_m) over periods m, where J_m is the score of
//!    period m. The step matrix is the partial diagonalization
//!    D~ = w diag(D) + (1 - w) D (RSiena's `diagonalize`, default 0.2), and
//!    one step theta = theta_1 - a_1 D~^-1 (mean S - s) is taken.
//! 2. **Phase 2**: `nsub` subphases (default 4). In each iteration one
//!    panel is simulated and theta <- theta - a D~^-1 (S - s), with
//!    deviations truncated at `truncation` standard deviations (estimated
//!    in phase 1). Subphase k runs at least 2.52^k (7 + p) and at most 200
//!    more iterations, and ends early once, for every coordinate, the sum
//!    of products of successive deviations is non-positive (the rule of
//!    Snijders 2001, appendix, after a proposal of Pflug). It ends by
//!    setting theta to its average over the subphase (averaged stochastic
//!    approximation; Polyak and Juditsky 1992, `polyak1992acceleration`)
//!    and multiplying the gain a by `reduceg` (default 0.5, from
//!    `firstg` = 0.2).
//! 3. **Phase 3**: n3 simulations (default 1000) at the final theta give
//!    the mean and covariance Sigma of S, the convergence t-ratios
//!    t_k = (mean S_k - s_k) / sd(S_k), the overall maximum convergence
//!    ratio sqrt((mean S - s)' Sigma^-1 (mean S - s)), the derivative D
//!    again, and the covariance matrix of the estimate by the delta method,
//!    cov(theta) = D^-1 Sigma D^-T, whose diagonal gives the standard
//!    errors. With conditional estimation the rate parameter of each
//!    period is the mean of the simulated stopping times (Snijders 2001,
//!    section 4.2), with their standard deviation reported as its
//!    standard error.
//!
//! Initial values follow the standard initial values of Snijders (2005,
//! `snijders2005models`, equation 11.41) in the form RSiena reports for
//! the out-degree effect: with c+ ties created, c- ties dissolved, and d
//! the mean density, beta_density = log(c+ d / (c- (1 - d))) / 2, every
//! other effect 0, and each rate 2 c_m / n (a heuristic: the changes per
//! actor, doubled because ministeps that undo earlier changes do not count
//! as changes). Phase 1 and phase 3 simulations run in parallel; every
//! simulation has its own random stream, so the estimate does not depend
//! on the thread count except in the last bits of floating-point sums.

use super::data::Panel;
use super::effects::{targets, Model};
use super::gof::{aux_len, aux_statistics};
use super::network::Network;
use super::simulate::{Simulator, Stopping};
use crate::error::{Error, Result};
use crate::linalg::{self, matmul, matvec, pseudo_inverse, quadratic_form, transpose};
use crate::matrix::Dense;
use crate::parallel;
use crate::rng::Stream;

/// Settings of the estimation algorithm.
#[derive(Debug, Clone)]
pub struct EstimateOptions {
    /// Conditional estimation (simulate until the observed number of
    /// changes; the default) or unconditional (simulate one time unit per
    /// period with the rates as parameters).
    pub conditional: bool,
    /// Phase 1 iterations; `None` for 7 + 3p.
    pub n1: Option<usize>,
    /// Number of phase 2 subphases.
    pub nsub: usize,
    /// Phase 3 iterations.
    pub n3: usize,
    /// Initial gain of the Robbins-Monro steps.
    pub firstg: f64,
    /// Factor by which the gain shrinks after each subphase.
    pub reduceg: f64,
    /// Weight of the diagonal in the step matrix, in [0, 1].
    pub diagonalize: f64,
    /// Deviations larger than this many standard deviations are truncated
    /// in phase 2 steps.
    pub truncation: f64,
    /// Estimation stops when a parameter exceeds this in absolute value.
    pub theta_bound: f64,
    /// Random seed.
    pub seed: u64,
    /// Threads for phases 1 and 3; 0 uses all cores.
    pub threads: usize,
    /// Whether to keep the auxiliary goodness-of-fit statistics of every
    /// phase 3 simulation.
    pub keep_aux: bool,
    /// Initial parameters, rates first when unconditional; `None` for the
    /// standard initial values. A NaN entry keeps the standard value.
    pub initial: Option<Vec<f64>>,
}

impl Default for EstimateOptions {
    fn default() -> Self {
        EstimateOptions {
            conditional: true,
            n1: None,
            nsub: 4,
            n3: 1000,
            firstg: 0.2,
            reduceg: 0.5,
            diagonalize: 0.2,
            truncation: 5.0,
            theta_bound: 50.0,
            seed: 0,
            threads: 0,
            keep_aux: true,
            initial: None,
        }
    }
}

/// The result of an estimation.
#[derive(Debug, Clone)]
pub struct Fit {
    /// Whether estimation was conditional.
    pub conditional: bool,
    /// Number of periods.
    pub periods: usize,
    /// Number of effects.
    pub effects: usize,
    /// Estimates: the M - 1 rates first when unconditional, then the
    /// effect parameters.
    pub theta: Vec<f64>,
    /// Standard errors of `theta`.
    pub se: Vec<f64>,
    /// Covariance matrix of `theta` (row-major p by p).
    pub covariance: Vec<f64>,
    /// Estimated derivative matrix dE S / dtheta (row-major p by p).
    pub derivative: Vec<f64>,
    /// Covariance matrix of the statistics in phase 3 (row-major).
    pub statistic_covariance: Vec<f64>,
    /// Convergence t-ratios.
    pub tconv: Vec<f64>,
    /// Overall maximum convergence ratio.
    pub tconv_max: f64,
    /// Observed statistics.
    pub targets: Vec<f64>,
    /// Mean simulated statistics in phase 3.
    pub means: Vec<f64>,
    /// Rate parameter of each period (estimated from the stopping times
    /// when conditional; the corresponding entries of `theta` otherwise).
    pub rates: Vec<f64>,
    /// Standard errors of `rates`.
    pub rate_se: Vec<f64>,
    /// Phase 2 trajectory: for each iteration, the subphase (1-based) and
    /// the parameter vector after the step (row-major, 1 + p columns).
    pub phase2_theta: Vec<f64>,
    /// Phase 2 deviations S - s per iteration (row-major, p columns).
    pub phase2_deviations: Vec<f64>,
    /// Phase 3 deviations S - s per simulation (row-major, p columns).
    pub phase3_deviations: Vec<f64>,
    /// Phase 3 simulated period lengths (row-major, M - 1 columns).
    pub phase3_times: Vec<f64>,
    /// Auxiliary statistics of each phase 3 simulation, summed over
    /// periods (row-major, `aux_len(n)` columns); empty unless kept.
    pub aux_simulated: Vec<f64>,
    /// Auxiliary statistics of the observed waves 2..M, summed.
    pub aux_observed: Vec<f64>,
    /// Iterations in phases 1, 2, and 3.
    pub iterations: [usize; 3],
}

/// One simulated panel reduced to what the algorithm needs.
struct Run {
    /// Statistics S (p values).
    stats: Vec<f64>,
    /// Per period: (S_m, J_m), each p values.
    per_period: Vec<(Vec<f64>, Vec<f64>)>,
    /// Stopping time (or 1) per period.
    times: Vec<f64>,
    /// Auxiliary statistics summed over periods, if requested.
    aux: Option<Vec<f64>>,
}

struct Problem<'a> {
    panel: &'a Panel,
    model: &'a Model,
    conditional: bool,
    rate_count: usize,
    p: usize,
    observed_distances: Vec<f64>,
    seed: u64,
}

impl<'a> Problem<'a> {
    fn run(
        &self,
        sim: &mut Simulator<'_>,
        theta: &[f64],
        rng: &mut Stream,
        aux: bool,
    ) -> Result<Run> {
        let (p, r, k) = (self.p, self.rate_count, self.model.len());
        let beta = &theta[r..];
        let mut stats = vec![0.0; p];
        let mut per_period = Vec::with_capacity(self.panel.period_count());
        let mut times = Vec::with_capacity(self.panel.period_count());
        let mut aux_sum = if aux {
            Some(vec![0.0; aux_len(self.panel.n())])
        } else {
            None
        };
        for m in 0..self.panel.period_count() {
            let stopping = if self.conditional {
                Stopping::Conditional {
                    distance: self.observed_distances[m] as usize,
                }
            } else {
                Stopping::Unconditional { rate: theta[m] }
            };
            let res = sim.period(m, self.panel.wave(m), beta, stopping, rng)?;
            let mut s_m = vec![0.0; p];
            let mut j_m = vec![0.0; p];
            if !self.conditional {
                s_m[m] = res.distance as f64;
                j_m[m] = res.rate_score;
            }
            s_m[r..].copy_from_slice(&res.statistics);
            j_m[r..].copy_from_slice(&res.scores);
            for kk in 0..p {
                stats[kk] += s_m[kk];
            }
            debug_assert_eq!(res.statistics.len(), k);
            per_period.push((s_m, j_m));
            times.push(res.time);
            if let Some(a) = aux_sum.as_mut() {
                for (x, y) in a.iter_mut().zip(aux_statistics(&res.network)) {
                    *x += y;
                }
            }
        }
        Ok(Run {
            stats,
            per_period,
            times,
            aux: aux_sum,
        })
    }

    /// Runs `count` independent simulations at `theta` in parallel; the
    /// simulation of item i draws from stream `base + i` of the seed.
    fn runs(
        &self,
        theta: &[f64],
        count: usize,
        base: u64,
        threads: usize,
        aux: bool,
    ) -> Result<Vec<Run>> {
        let seed = self.seed;
        let parts = parallel::for_each_item(
            count,
            threads,
            || (Simulator::new(self.panel, self.model), Vec::new()),
            |(sim, out): &mut (Simulator<'_>, Vec<(usize, Result<Run>)>), item| {
                let mut rng = Stream::new(seed, base + item as u64);
                out.push((item, self.run(sim, theta, &mut rng, aux)));
            },
        );
        let mut all: Vec<(usize, Result<Run>)> = parts.into_iter().flat_map(|(_, v)| v).collect();
        all.sort_by_key(|(i, _)| *i);
        all.into_iter().map(|(_, r)| r).collect()
    }
}

/// Summary of a set of runs: means, covariance, and derivative.
struct Moments {
    mean: Vec<f64>,
    cov: Vec<f64>,
    derivative: Vec<f64>,
}

fn moments(runs: &[Run], p: usize, periods: usize) -> Moments {
    let count = runs.len() as f64;
    let mut mean = vec![0.0; p];
    for r in runs {
        for (m, s) in mean.iter_mut().zip(&r.stats) {
            *m += s;
        }
    }
    for v in mean.iter_mut() {
        *v /= count;
    }
    let mut cov = vec![0.0; p * p];
    for r in runs {
        for a in 0..p {
            let da = r.stats[a] - mean[a];
            for b in 0..p {
                cov[a * p + b] += da * (r.stats[b] - mean[b]);
            }
        }
    }
    for v in cov.iter_mut() {
        *v /= count;
    }
    // D = sum_m [ mean(S_m J_m') - mean(S_m) mean(J_m)' ].
    let mut derivative = vec![0.0; p * p];
    let mut sum_s = vec![vec![0.0; p]; periods];
    let mut sum_j = vec![vec![0.0; p]; periods];
    for r in runs {
        for (m, (s_m, j_m)) in r.per_period.iter().enumerate() {
            for a in 0..p {
                sum_s[m][a] += s_m[a];
                for b in 0..p {
                    derivative[a * p + b] += s_m[a] * j_m[b];
                }
            }
            for b in 0..p {
                sum_j[m][b] += j_m[b];
            }
        }
    }
    for a in 0..p {
        for b in 0..p {
            let mut v = derivative[a * p + b] / count;
            for m in 0..periods {
                v -= sum_s[m][a] * sum_j[m][b] / (count * count);
            }
            derivative[a * p + b] = v;
        }
    }
    Moments {
        mean,
        cov,
        derivative,
    }
}

/// A row-major `p` by `p` slice as a matrix.
fn square(a: &[f64], p: usize) -> Dense {
    Dense::from_row_major(p, p, a.to_vec()).expect("a p by p matrix")
}

/// The inverse of a `p` by `p` matrix of the estimation algorithm, or an
/// error naming the usual cause when it is singular.
fn invert(a: &[f64], p: usize) -> Result<Dense> {
    linalg::invert(&square(a, p), 1e-13).ok_or_else(|| {
        Error::InvalidArgument(
            "a matrix in the estimation algorithm is singular (the derivative or \
             covariance matrix of the statistics); some effects may be perfectly \
             collinear or have statistics that never vary"
                .into(),
        )
    })
}

/// `w diag(D) + (1 - w) D`, inverted.
fn step_matrix(derivative: &[f64], p: usize, w: f64) -> Result<Dense> {
    let mut d = vec![0.0; p * p];
    for a in 0..p {
        for b in 0..p {
            let full = derivative[a * p + b];
            d[a * p + b] = if a == b { full } else { (1.0 - w) * full };
        }
    }
    for a in 0..p {
        if d[a * p + a].is_nan() || d[a * p + a] <= 0.0 {
            return Err(Error::InvalidArgument(format!(
                "the statistic of parameter {} does not respond to its parameter (estimated \
                 derivative {:.3}); the data carry no information about this effect at the \
                 current values. Remove the effect or change its initial value",
                a + 1,
                d[a * p + a]
            )));
        }
    }
    invert(&d, p)
}

/// Standard initial values (see the module documentation).
pub fn initial_values(panel: &Panel, model: &Model, conditional: bool) -> Vec<f64> {
    let n = panel.n() as f64;
    let periods = panel.period_count();
    let mut theta = Vec::new();
    if !conditional {
        for d in panel.distances() {
            theta.push((2.0 * d / n).max(0.1));
        }
    }
    let (mut up, mut down) = (0usize, 0usize);
    for (u, d) in panel.creations_and_dissolutions() {
        up += u;
        down += d;
    }
    let pairs = (n * (n - 1.0)).max(1.0);
    let density = (0..=periods)
        .map(|w| panel.wave(w).tie_count() as f64 / pairs)
        .sum::<f64>()
        / (periods + 1) as f64;
    let density = density.clamp(1e-6, 1.0 - 1e-6);
    let density_start = if up > 0 && down > 0 {
        0.5 * ((up as f64 * density) / (down as f64 * (1.0 - density))).ln()
    } else {
        (density / (1.0 - density)).ln()
    };
    for e in model.effects() {
        theta.push(if e.kind == super::effects::EffectKind::Density {
            density_start
        } else {
            0.0
        });
    }
    theta
}

fn check_options(opts: &EstimateOptions) -> Result<()> {
    let positive = |v: f64, what: &str| -> Result<()> {
        if !(v > 0.0 && v.is_finite()) {
            return Err(Error::InvalidArgument(format!(
                "{what} must be positive and finite, not {v}"
            )));
        }
        Ok(())
    };
    positive(opts.firstg, "firstg")?;
    positive(opts.truncation, "truncation")?;
    positive(opts.theta_bound, "theta_bound")?;
    if !(opts.reduceg > 0.0 && opts.reduceg <= 1.0) {
        return Err(Error::InvalidArgument(format!(
            "reduceg must be in (0, 1], not {}",
            opts.reduceg
        )));
    }
    if !(0.0..=1.0).contains(&opts.diagonalize) {
        return Err(Error::InvalidArgument(format!(
            "diagonalize must be in [0, 1], not {}",
            opts.diagonalize
        )));
    }
    if opts.n3 < 2 {
        return Err(Error::InvalidArgument(
            "n3 must be at least 2 (at least several hundred for usable standard errors)".into(),
        ));
    }
    Ok(())
}

/// Estimates the model on the panel.
pub fn estimate(panel: &Panel, model: &Model, opts: &EstimateOptions) -> Result<Fit> {
    check_options(opts)?;
    let periods = panel.period_count();
    let k = model.len();
    let rate_count = if opts.conditional { 0 } else { periods };
    let p = rate_count + k;
    let observed_distances = panel.distances();
    if let Some(m) = observed_distances.iter().position(|&d| d == 0.0) {
        return Err(Error::InvalidArgument(format!(
            "no tie changed between waves {} and {}; a period without changes cannot be \
             modelled (its rate would be zero). Merge or drop that wave",
            m + 1,
            m + 2
        )));
    }
    let mut target = vec![0.0; p];
    if !opts.conditional {
        target[..periods].copy_from_slice(&observed_distances);
    }
    let effect_targets = targets(model, panel);
    for m in 0..periods {
        for kk in 0..k {
            target[rate_count + kk] += effect_targets[m * k + kk];
        }
    }
    let problem = Problem {
        panel,
        model,
        conditional: opts.conditional,
        rate_count,
        p,
        observed_distances,
        seed: opts.seed,
    };

    // Initial values.
    let mut theta = initial_values(panel, model, opts.conditional);
    if let Some(init) = &opts.initial {
        if init.len() != p {
            return Err(Error::LengthMismatch {
                what: "initial parameter values",
                expected: p,
                found: init.len(),
            });
        }
        for (t, &v) in theta.iter_mut().zip(init) {
            if !v.is_nan() {
                *t = v;
            }
        }
    }
    if !opts.conditional {
        for (m, &rho) in theta[..periods].iter().enumerate() {
            if rho.is_nan() || rho <= 0.0 {
                return Err(Error::InvalidArgument(format!(
                    "the initial rate of period {} must be positive, not {rho}",
                    m + 1
                )));
            }
        }
    }
    let check_bound = |theta: &[f64]| -> Result<()> {
        if theta
            .iter()
            .any(|v| !v.is_finite() || v.abs() > opts.theta_bound)
        {
            return Err(Error::InvalidArgument(format!(
                "the estimation diverged: a parameter exceeded {} in absolute value. Check \
                 the model specification (collinear effects, effects whose statistics \
                 cannot change), try other initial values, or unconditional estimation",
                opts.theta_bound
            )));
        }
        Ok(())
    };

    // Phase 1.
    let n1 = opts.n1.unwrap_or(7 + 3 * p).max(2);
    let runs1 = problem.runs(&theta, n1, 0, opts.threads, false)?;
    let m1 = moments(&runs1, p, periods);
    let sd1: Vec<f64> = (0..p).map(|a| m1.cov[a * p + a].sqrt()).collect();
    // The step matrix stays as estimated in phase 1 throughout phase 2.
    let step = step_matrix(&m1.derivative, p, opts.diagonalize)?;
    {
        let dev: Vec<f64> = (0..p).map(|a| m1.mean[a] - target[a]).collect();
        let delta = matvec(&step, &dev);
        for a in 0..p {
            theta[a] -= opts.firstg * delta[a];
        }
        check_bound(&theta)?;
    }

    // Phase 2.
    let mut phase2_theta = Vec::new();
    let mut phase2_dev = Vec::new();
    let mut gain = opts.firstg;
    let mut iteration = 0u64;
    let mut sim = Simulator::new(panel, model);
    for sub in 1..=opts.nsub {
        let n2min = (2.52f64.powi(sub as i32) * (7 + p) as f64).round() as usize;
        let n2max = n2min + 200;
        let mut count = 0usize;
        let mut sum_theta = vec![0.0; p];
        let mut ac = vec![0.0; p];
        let mut prev: Option<Vec<f64>> = None;
        loop {
            let mut rng = Stream::new(opts.seed, (1u64 << 32) + iteration);
            iteration += 1;
            let run = problem.run(&mut sim, &theta, &mut rng, false)?;
            let mut dev: Vec<f64> = (0..p).map(|a| run.stats[a] - target[a]).collect();
            phase2_dev.extend_from_slice(&dev);
            for a in 0..p {
                if sd1[a] > 0.0 {
                    let z = dev[a].abs() / sd1[a];
                    if z > opts.truncation {
                        dev[a] *= opts.truncation / z;
                    }
                }
            }
            let delta = matvec(&step, &dev);
            for a in 0..p {
                theta[a] -= gain * delta[a];
            }
            check_bound(&theta)?;
            if !opts.conditional {
                for rho in theta[..periods].iter_mut() {
                    // A non-positive rate has no meaning; keep it small
                    // and positive so the next simulation can run.
                    if *rho < 1e-3 {
                        *rho = 1e-3;
                    }
                }
            }
            count += 1;
            for a in 0..p {
                sum_theta[a] += theta[a];
            }
            if let Some(pv) = &prev {
                for a in 0..p {
                    ac[a] += dev[a] * pv[a];
                }
            }
            prev = Some(dev);
            phase2_theta.push(sub as f64);
            phase2_theta.extend_from_slice(&theta);
            let settled = ac.iter().all(|&v| v <= 0.0);
            if count >= n2max || (count >= n2min && settled) {
                break;
            }
        }
        for a in 0..p {
            theta[a] = sum_theta[a] / count as f64;
        }
        gain *= opts.reduceg;
    }
    let n2_total = iteration as usize;

    // Phase 3.
    let runs3 = problem.runs(&theta, opts.n3, 2u64 << 32, opts.threads, opts.keep_aux)?;
    let m3 = moments(&runs3, p, periods);
    let dev: Vec<f64> = (0..p).map(|a| m3.mean[a] - target[a]).collect();
    let tconv: Vec<f64> = (0..p)
        .map(|a| {
            let sd = m3.cov[a * p + a].sqrt();
            if sd > 0.0 {
                dev[a] / sd
            } else if dev[a] == 0.0 {
                0.0
            } else {
                f64::INFINITY
            }
        })
        .collect();
    let tconv_max = match invert(&m3.cov, p) {
        Ok(inv) => quadratic_form(&inv, &dev).max(0.0).sqrt(),
        Err(_) => {
            let (pinv, _) = pseudo_inverse(&square(&m3.cov, p), 1e-10)?;
            quadratic_form(&pinv, &dev).max(0.0).sqrt()
        }
    };
    let (covariance, se) = match invert(&m3.derivative, p) {
        Ok(dinv) => {
            let cov = matmul(&matmul(&dinv, &square(&m3.cov, p)), &transpose(&dinv));
            let se = (0..p).map(|a| cov.get(a, a).max(0.0).sqrt()).collect();
            (cov.as_slice().to_vec(), se)
        }
        Err(_) => (vec![f64::NAN; p * p], vec![f64::NAN; p]),
    };
    let mut phase3_dev = Vec::with_capacity(opts.n3 * p);
    let mut phase3_times = Vec::with_capacity(opts.n3 * periods);
    let mut aux_simulated = Vec::new();
    for r in &runs3 {
        phase3_dev.extend(r.stats.iter().zip(&target).map(|(s, t)| s - t));
        phase3_times.extend_from_slice(&r.times);
        if let Some(a) = &r.aux {
            aux_simulated.extend_from_slice(a);
        }
    }
    let (rates, rate_se) = if opts.conditional {
        let count = runs3.len() as f64;
        let mut rates = vec![0.0; periods];
        for r in &runs3 {
            for (rate, t) in rates.iter_mut().zip(&r.times) {
                *rate += t;
            }
        }
        for v in rates.iter_mut() {
            *v /= count;
        }
        let mut var = vec![0.0; periods];
        for r in &runs3 {
            for m in 0..periods {
                var[m] += (r.times[m] - rates[m]).powi(2);
            }
        }
        let se = var.iter().map(|v| (v / (count - 1.0)).sqrt()).collect();
        (rates, se)
    } else {
        (theta[..periods].to_vec(), se[..periods].to_vec())
    };
    let mut aux_observed = vec![0.0; aux_len(panel.n())];
    for w in 1..=periods {
        for (x, y) in aux_observed.iter_mut().zip(aux_statistics(panel.wave(w))) {
            *x += y;
        }
    }
    Ok(Fit {
        conditional: opts.conditional,
        periods,
        effects: k,
        theta,
        se,
        covariance,
        derivative: m3.derivative,
        statistic_covariance: m3.cov,
        tconv,
        tconv_max,
        targets: target,
        means: m3.mean,
        rates,
        rate_se,
        phase2_theta,
        phase2_deviations: phase2_dev,
        phase3_deviations: phase3_dev,
        phase3_times,
        aux_simulated,
        aux_observed,
        iterations: [n1, n2_total, opts.n3],
    })
}

/// The network at the end of each period of one simulated panel at
/// `theta`, for callers that want simulated data rather than estimates.
pub fn simulated_networks(
    panel: &Panel,
    model: &Model,
    theta: &[f64],
    conditional: bool,
    seed: u64,
    index: u64,
) -> Result<Vec<Network>> {
    let periods = panel.period_count();
    let rates = if conditional {
        &[][..]
    } else {
        &theta[..periods]
    };
    let beta = if conditional {
        theta
    } else {
        &theta[periods..]
    };
    Ok(
        super::simulate::simulate_panel(panel, model, beta, rates, conditional, seed, index)?
            .networks,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::saom::effects::{Effect, EffectKind};

    /// A first wave with the given density, and later waves that each
    /// toggle about 12% of the tie variables of the wave before.
    fn random_panel(seed: u64, n: usize, waves: usize, density: f64) -> Panel {
        let mut rng = Stream::new(seed, 0);
        let mut adj = Vec::new();
        let mut a = vec![0u8; n * n];
        for i in 0..n {
            for j in 0..n {
                if i != j && rng.uniform() < density {
                    a[i * n + j] = 1;
                }
            }
        }
        adj.push(a.clone());
        for _ in 1..waves {
            for i in 0..n {
                for j in 0..n {
                    if i != j && rng.uniform() < 0.12 {
                        a[i * n + j] ^= 1;
                    }
                }
            }
            adj.push(a.clone());
        }
        Panel::from_adjacency(n, &adj).unwrap()
    }

    #[test]
    fn initial_values_follow_the_formula() {
        // Wave 1: 0 -> 1. Wave 2: 0 -> 1, 1 -> 2 (one creation, no
        // dissolution): the fallback logit of the mean density applies.
        let p = Panel::from_edge_lists(3, &[vec![(0, 1)], vec![(0, 1), (1, 2)]]).unwrap();
        let m = Model::new(&p, &[Effect::structural(EffectKind::Density)]).unwrap();
        let d: f64 = (1.0 / 6.0 + 2.0 / 6.0) / 2.0;
        let theta = initial_values(&p, &m, true);
        assert!((theta[0] - (d / (1.0 - d)).ln()).abs() < 1e-12);
        let theta = initial_values(&p, &m, false);
        assert_eq!(theta.len(), 2);
        assert!((theta[0] - 2.0 / 3.0).abs() < 1e-12);
        // With creations and dissolutions: log(c+ d / (c- (1 - d))) / 2.
        let p = Panel::from_edge_lists(3, &[vec![(0, 1), (2, 1)], vec![(0, 1), (1, 2)]]).unwrap();
        let theta = initial_values(&p, &m, true);
        let d: f64 = 2.0 / 6.0;
        assert!((theta[0] - 0.5 * ((1.0 * d) / (1.0 * (1.0 - d))).ln()).abs() < 1e-12);
    }

    #[test]
    fn estimation_converges_on_a_small_panel_and_is_reproducible() {
        let p = random_panel(21, 15, 3, 0.2);
        let m = Model::new(
            &p,
            &[
                Effect::structural(EffectKind::Density),
                Effect::structural(EffectKind::Recip),
            ],
        )
        .unwrap();
        let opts = EstimateOptions {
            n3: 200,
            nsub: 2,
            seed: 4,
            threads: 1,
            ..Default::default()
        };
        let fit = estimate(&p, &m, &opts).unwrap();
        assert_eq!(fit.theta.len(), 2);
        assert!(fit.se.iter().all(|s| s.is_finite() && *s > 0.0));
        assert!(fit.tconv_max.is_finite());
        assert_eq!(fit.rates.len(), 2);
        assert!(fit.rates.iter().all(|r| *r > 0.0));
        assert_eq!(fit.phase3_deviations.len(), 200 * 2);
        assert_eq!(fit.aux_simulated.len(), 200 * aux_len(15));
        let again = estimate(
            &p,
            &m,
            &EstimateOptions {
                threads: 2,
                ..opts.clone()
            },
        )
        .unwrap();
        for (a, b) in fit.theta.iter().zip(&again.theta) {
            assert!((a - b).abs() < 1e-9, "{a} vs {b}");
        }
        for (a, b) in fit.phase3_deviations.iter().zip(&again.phase3_deviations) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn unconditional_estimation_has_rate_parameters() {
        let p = random_panel(8, 12, 2, 0.25);
        let m = Model::new(&p, &[Effect::structural(EffectKind::Density)]).unwrap();
        let opts = EstimateOptions {
            conditional: false,
            n3: 150,
            nsub: 2,
            seed: 9,
            threads: 1,
            keep_aux: false,
            ..Default::default()
        };
        let fit = estimate(&p, &m, &opts).unwrap();
        assert_eq!(fit.theta.len(), 2);
        assert_eq!(fit.rates, fit.theta[..1].to_vec());
        assert!(fit.aux_simulated.is_empty());
        assert_eq!(fit.targets[0], p.distances()[0]);
    }

    #[test]
    fn rejects_periods_without_change_and_bad_options() {
        let p = Panel::from_edge_lists(3, &[vec![(0, 1)], vec![(0, 1)]]).unwrap();
        let m = Model::new(&p, &[Effect::structural(EffectKind::Density)]).unwrap();
        assert!(estimate(&p, &m, &EstimateOptions::default()).is_err());
        let p = random_panel(3, 8, 2, 0.3);
        let m = Model::new(&p, &[Effect::structural(EffectKind::Density)]).unwrap();
        let bad = EstimateOptions {
            diagonalize: 2.0,
            ..Default::default()
        };
        assert!(estimate(&p, &m, &bad).is_err());
        let bad = EstimateOptions {
            initial: Some(vec![1.0, 2.0]),
            ..Default::default()
        };
        assert!(estimate(&p, &m, &bad).is_err());
    }
}
