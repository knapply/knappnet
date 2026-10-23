//! C ABI for stochastic actor-oriented models (`engine::saom`).
//!
//! A panel is built with `net_saom_new` from adjacency matrices (one per
//! wave, each n-by-n column-major as R stores matrices, waves one after
//! another) and completed with covariates. Effects are passed as three
//! parallel arrays of length k: `codes` (positions in
//! `EffectKind::ALL`), `covariates` (0-based index of the actor or dyadic
//! covariate an effect uses, or -1), and `parameters` (the internal
//! effect parameter; GWESP uses parameter / 100 as its weight).
//!
//! Parameter vectors are ordered with the M - 1 rate parameters first when
//! estimation is unconditional, then the k effect parameters; conditional
//! estimation has only the effect parameters.
//!
//! Estimation returns an opaque `SaomFit`, read with `net_saom_fit_length`
//! and `net_saom_fit_copy` by item number (see `SAOM_FIT_*`), and released
//! with `net_saom_fit_free`. Matrices are row-major unless stated.

use std::ffi::c_char;

use engine::saom::effects::{change_statistics_on_wave, targets, Effect, EffectKind, Model};
use engine::saom::estimate::{estimate, initial_values, EstimateOptions, Fit};
use engine::saom::gof::aux_len;
use engine::saom::simulate::simulate_panel;
use engine::saom::Panel;
use engine::Error;

use crate::{count, deref, guard, slice, slice_mut, thread_count};

/// Opaque panel handle.
pub struct SaomHandle(Panel);

/// Opaque estimation result.
pub struct SaomFit(Fit);

/// Number of doubles in the `options` vector of `net_saom_estimate`:
/// conditional, n1 (0 = default), nsub, n3, firstg, reduceg, diagonalize,
/// truncation, theta_bound, keep_aux.
pub const SAOM_OPTIONS_LEN: usize = 10;

/// Items of a `SaomFit`, for `net_saom_fit_length` and `net_saom_fit_copy`.
pub const SAOM_FIT_THETA: i32 = 0;
/// Standard errors (p).
pub const SAOM_FIT_SE: i32 = 1;
/// Covariance matrix of the estimates (p by p).
pub const SAOM_FIT_COVARIANCE: i32 = 2;
/// Derivative matrix (p by p).
pub const SAOM_FIT_DERIVATIVE: i32 = 3;
/// Covariance matrix of the statistics (p by p).
pub const SAOM_FIT_STATISTIC_COVARIANCE: i32 = 4;
/// Convergence t-ratios (p).
pub const SAOM_FIT_TCONV: i32 = 5;
/// Overall maximum convergence ratio (1).
pub const SAOM_FIT_TCONV_MAX: i32 = 6;
/// Observed statistics (p).
pub const SAOM_FIT_TARGETS: i32 = 7;
/// Mean simulated statistics in phase 3 (p).
pub const SAOM_FIT_MEANS: i32 = 8;
/// Rate parameters (M - 1).
pub const SAOM_FIT_RATES: i32 = 9;
/// Standard errors of the rates (M - 1).
pub const SAOM_FIT_RATE_SE: i32 = 10;
/// Phase 2 trajectory (iterations by 1 + p: subphase, then theta).
pub const SAOM_FIT_PHASE2_THETA: i32 = 11;
/// Phase 2 deviations (iterations by p).
pub const SAOM_FIT_PHASE2_DEVIATIONS: i32 = 12;
/// Phase 3 deviations (n3 by p).
pub const SAOM_FIT_PHASE3_DEVIATIONS: i32 = 13;
/// Phase 3 simulated period lengths (n3 by M - 1).
pub const SAOM_FIT_PHASE3_TIMES: i32 = 14;
/// Auxiliary statistics of phase 3 simulations (n3 by 3n + 16).
pub const SAOM_FIT_AUX_SIMULATED: i32 = 15;
/// Auxiliary statistics of the observed waves 2..M (3n + 16).
pub const SAOM_FIT_AUX_OBSERVED: i32 = 16;
/// Summary: conditional, periods, effects, n1, n2, n3, p, 3n + 16 (8).
pub const SAOM_FIT_META: i32 = 17;

fn panel<'a>(h: *const SaomHandle) -> Result<&'a Panel, Error> {
    // SAFETY: non-null handles are only created by `net_saom_new`.
    unsafe { deref(h, "panel") }.map(|p| &p.0)
}

fn panel_mut<'a>(h: *mut SaomHandle) -> Result<&'a mut Panel, Error> {
    if h.is_null() {
        return Err(Error::InvalidArgument("panel handle is null".into()));
    }
    // SAFETY: non-null handles are only created by `net_saom_new`, and the
    // caller does not use the panel concurrently while adding covariates.
    Ok(unsafe { &mut (*h).0 })
}

fn fit<'a>(h: *const SaomFit) -> Result<&'a Fit, Error> {
    // SAFETY: non-null fit handles are only created by `net_saom_estimate`.
    unsafe { deref(h, "fit") }.map(|f| &f.0)
}

/// # Safety
/// The three arrays must hold `k` elements.
unsafe fn effects(
    k: i64,
    codes: *const i32,
    covariates: *const i32,
    parameters: *const f64,
) -> Result<Vec<Effect>, Error> {
    let k = count(k, "effect count")?;
    let codes = slice(codes, k);
    let covs = slice(covariates, k);
    let pars = slice(parameters, k);
    let mut out = Vec::with_capacity(k);
    for i in 0..k {
        let kind = EffectKind::from_code(codes[i])?;
        let covariate = if covs[i] < 0 {
            None
        } else {
            Some(covs[i] as usize)
        };
        out.push(Effect {
            kind,
            covariate,
            parameter: if kind.has_parameter() { pars[i] } else { 0.0 },
        });
    }
    Ok(out)
}

/// Builds a panel from `waves` adjacency matrices, each n-by-n
/// column-major (entry (i, j) at `j * n + i`), one after another. Any
/// non-zero entry is a tie; the diagonal is ignored. Release the handle
/// with `net_saom_free`.
///
/// # Safety
/// `adjacency` must hold `waves * n * n` doubles; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_saom_new(
    n: i64,
    waves: i64,
    adjacency: *const f64,
    out: *mut *mut SaomHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "actor count")?;
        let waves = count(waves, "wave count")?;
        let values = slice(adjacency, waves * n * n);
        let mut mats = Vec::with_capacity(waves);
        for w in 0..waves {
            let block = &values[w * n * n..(w + 1) * n * n];
            let mut a = vec![0u8; n * n];
            for j in 0..n {
                for i in 0..n {
                    let v = block[j * n + i];
                    if v.is_nan() {
                        return Err(Error::InvalidArgument(format!(
                            "wave {} has a missing tie variable ({}, {}); missing tie data \
                             are not supported",
                            w + 1,
                            i + 1,
                            j + 1
                        )));
                    }
                    a[i * n + j] = (v != 0.0) as u8;
                }
            }
            mats.push(a);
        }
        let p = Panel::from_adjacency(n, &mats)?;
        *out = Box::into_raw(Box::new(SaomHandle(p)));
        Ok(())
    })
}

/// Releases a panel handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_saom_free(h: *mut SaomHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Adds an actor covariate: `len` is n for a constant covariate or
/// n * waves (wave-major, as an n-by-waves column-major matrix) for a
/// changing one. The covariate's 0-based index is the number of actor
/// covariates added before it.
///
/// # Safety
/// `h` must be a live handle not in use elsewhere; `values` must hold
/// `len` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_saom_add_actor_covariate(
    h: *mut SaomHandle,
    values: *const f64,
    len: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel_mut(h)?;
        let len = count(len, "covariate length")?;
        if len != p.n() && len != p.n() * p.wave_count() {
            return Err(Error::InvalidArgument(format!(
                "an actor covariate needs one value per actor ({}) or one per actor and wave \
                 ({}), not {len}",
                p.n(),
                p.n() * p.wave_count()
            )));
        }
        p.add_actor_covariate(slice(values, len))?;
        Ok(())
    })
}

/// Adds a constant dyadic covariate (n-by-n, column-major). Its 0-based
/// index is the number of dyadic covariates added before it.
///
/// # Safety
/// `h` must be a live handle not in use elsewhere; `values` must hold
/// `n * n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_saom_add_dyadic_covariate(
    h: *mut SaomHandle,
    values: *const f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel_mut(h)?;
        let n = p.n();
        let v = slice(values, n * n);
        let mut row_major = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                row_major[i * n + j] = v[j * n + i];
            }
        }
        p.add_dyadic_covariate(&row_major)?;
        Ok(())
    })
}

/// Number of periods (waves - 1), or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_saom_period_count(h: *const SaomHandle) -> i64 {
    panel(h).map(|p| p.period_count() as i64).unwrap_or(-1)
}

/// Number of tie variables that changed in each period (M - 1 doubles).
///
/// # Safety
/// `h` must be a live handle; `out` must hold M - 1 doubles.
#[no_mangle]
pub unsafe extern "C" fn net_saom_distances(
    h: *const SaomHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel(h)?;
        let d = p.distances();
        slice_mut(out, d.len()).copy_from_slice(&d);
        Ok(())
    })
}

/// Target statistics: for each period m, the k effect statistics of wave
/// m + 1 with the covariates of wave m, period-major (`out[m * k + kk]`).
///
/// # Safety
/// `h` must be a live handle; the effect arrays must hold `k` elements;
/// `out` must hold `k * (M - 1)` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_saom_targets(
    h: *const SaomHandle,
    k: i64,
    codes: *const i32,
    covariates: *const i32,
    parameters: *const f64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel(h)?;
        let model = Model::new(p, &effects(k, codes, covariates, parameters)?)?;
        let t = targets(&model, p);
        slice_mut(out, t.len()).copy_from_slice(&t);
        Ok(())
    })
}

/// Change statistics Delta_kij of every effect for actor `actor` (0-based)
/// on wave `wave` (0-based), with the covariates of period `period`
/// (0-based): `n * k` doubles, alter-major (`out[j * k + kk]`), zero for
/// j = actor.
///
/// # Safety
/// `h` must be a live handle; the effect arrays must hold `k` elements;
/// `out` must hold `n * k` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_saom_contributions(
    h: *const SaomHandle,
    k: i64,
    codes: *const i32,
    covariates: *const i32,
    parameters: *const f64,
    wave: i64,
    period: i64,
    actor: i64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel(h)?;
        let model = Model::new(p, &effects(k, codes, covariates, parameters)?)?;
        let wave = count(wave, "wave")?;
        let period = count(period, "period")?;
        let actor = count(actor, "actor")?;
        if wave >= p.wave_count() || period >= p.period_count() || actor >= p.n() {
            return Err(Error::InvalidArgument(
                "wave, period, or actor index out of range".into(),
            ));
        }
        let d = change_statistics_on_wave(&model, p, period, wave, actor);
        slice_mut(out, d.len()).copy_from_slice(&d);
        Ok(())
    })
}

/// Standard initial parameter values (p doubles; see the engine's
/// `saom::estimate` documentation).
///
/// # Safety
/// `h` must be a live handle; the effect arrays must hold `k` elements;
/// `out` must hold `k` doubles, or `k + M - 1` when `conditional == 0`.
#[no_mangle]
pub unsafe extern "C" fn net_saom_initial(
    h: *const SaomHandle,
    k: i64,
    codes: *const i32,
    covariates: *const i32,
    parameters: *const f64,
    conditional: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel(h)?;
        let model = Model::new(p, &effects(k, codes, covariates, parameters)?)?;
        let theta = initial_values(p, &model, conditional != 0);
        slice_mut(out, theta.len()).copy_from_slice(&theta);
        Ok(())
    })
}

/// Simulates every period once from the observed network at its start.
/// `theta` holds the k effect parameters; `rates` the M - 1 rate
/// parameters (ignored, and may be null, when `conditional != 0`). Stream
/// `index` of `seed` supplies the random numbers. `networks` receives the
/// M - 1 simulated end networks as n-by-n column-major 0/1 matrices;
/// `statistics` the k effect statistics per period (period-major); `times`
/// the simulated period lengths (the stopping times under rate 1 when
/// conditional, otherwise 1). Each output pointer may be null.
///
/// # Safety
/// `h` must be a live handle; the arrays must have the stated lengths.
#[no_mangle]
pub unsafe extern "C" fn net_saom_simulate(
    h: *const SaomHandle,
    k: i64,
    codes: *const i32,
    covariates: *const i32,
    parameters: *const f64,
    theta: *const f64,
    rates: *const f64,
    conditional: i32,
    seed: u64,
    index: u64,
    networks: *mut f64,
    statistics: *mut f64,
    times: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel(h)?;
        let model = Model::new(p, &effects(k, codes, covariates, parameters)?)?;
        let beta = slice(theta, model.len());
        let periods = p.period_count();
        let rates = if conditional != 0 {
            &[][..]
        } else {
            slice(rates, periods)
        };
        let sim = simulate_panel(p, &model, beta, rates, conditional != 0, seed, index)?;
        let n = p.n();
        if !networks.is_null() {
            let out = slice_mut(networks, periods * n * n);
            for (m, net) in sim.networks.iter().enumerate() {
                let adj = net.adjacency();
                let block = &mut out[m * n * n..(m + 1) * n * n];
                for i in 0..n {
                    for j in 0..n {
                        block[j * n + i] = adj[i * n + j] as f64;
                    }
                }
            }
        }
        if !statistics.is_null() {
            slice_mut(statistics, sim.statistics.len()).copy_from_slice(&sim.statistics);
        }
        if !times.is_null() {
            slice_mut(times, sim.times.len()).copy_from_slice(&sim.times);
        }
        Ok(())
    })
}

/// Estimates the model by the method of moments. `options` holds
/// `SAOM_OPTIONS_LEN` doubles: conditional (0/1), n1 (0 for the default
/// 7 + 3p), nsub, n3, firstg, reduceg, diagonalize, truncation,
/// theta_bound, keep_aux (0/1). `initial` may be null for the standard
/// initial values; otherwise it holds p values, NaN keeping the standard
/// value of that parameter. `threads == 0` uses all cores. On success
/// `*out` receives a fit to release with `net_saom_fit_free`.
///
/// # Safety
/// `h` must be a live handle; the effect arrays must hold `k` elements;
/// `options` must hold `SAOM_OPTIONS_LEN` doubles; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_saom_estimate(
    h: *const SaomHandle,
    k: i64,
    codes: *const i32,
    covariates: *const i32,
    parameters: *const f64,
    options: *const f64,
    initial: *const f64,
    seed: u64,
    threads: i32,
    out: *mut *mut SaomFit,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let p = panel(h)?;
        let model = Model::new(p, &effects(k, codes, covariates, parameters)?)?;
        let o = slice(options, SAOM_OPTIONS_LEN);
        let whole = |v: f64, what: &str| -> Result<usize, Error> {
            if !(v.is_finite() && v >= 0.0 && v.fract() == 0.0) {
                return Err(Error::InvalidArgument(format!(
                    "{what} must be a non-negative whole number, not {v}"
                )));
            }
            Ok(v as usize)
        };
        let conditional = o[0] != 0.0;
        let n1 = whole(o[1], "n1")?;
        let periods = p.period_count();
        let pdim = model.len() + if conditional { 0 } else { periods };
        let initial = if initial.is_null() {
            None
        } else {
            Some(slice(initial, pdim).to_vec())
        };
        let opts = EstimateOptions {
            conditional,
            n1: if n1 == 0 { None } else { Some(n1) },
            nsub: whole(o[2], "nsub")?,
            n3: whole(o[3], "n3")?,
            firstg: o[4],
            reduceg: o[5],
            diagonalize: o[6],
            truncation: o[7],
            theta_bound: o[8],
            seed,
            threads: thread_count(threads),
            keep_aux: o[9] != 0.0,
            initial,
        };
        let f = estimate(p, &model, &opts)?;
        *out = Box::into_raw(Box::new(SaomFit(f)));
        Ok(())
    })
}

fn item(f: &Fit, which: i32) -> Option<Vec<f64>> {
    let (p, periods) = (f.theta.len(), f.periods);
    let n = (f.aux_observed.len() - 16) / 3;
    Some(match which {
        SAOM_FIT_THETA => f.theta.clone(),
        SAOM_FIT_SE => f.se.clone(),
        SAOM_FIT_COVARIANCE => f.covariance.clone(),
        SAOM_FIT_DERIVATIVE => f.derivative.clone(),
        SAOM_FIT_STATISTIC_COVARIANCE => f.statistic_covariance.clone(),
        SAOM_FIT_TCONV => f.tconv.clone(),
        SAOM_FIT_TCONV_MAX => vec![f.tconv_max],
        SAOM_FIT_TARGETS => f.targets.clone(),
        SAOM_FIT_MEANS => f.means.clone(),
        SAOM_FIT_RATES => f.rates.clone(),
        SAOM_FIT_RATE_SE => f.rate_se.clone(),
        SAOM_FIT_PHASE2_THETA => f.phase2_theta.clone(),
        SAOM_FIT_PHASE2_DEVIATIONS => f.phase2_deviations.clone(),
        SAOM_FIT_PHASE3_DEVIATIONS => f.phase3_deviations.clone(),
        SAOM_FIT_PHASE3_TIMES => f.phase3_times.clone(),
        SAOM_FIT_AUX_SIMULATED => f.aux_simulated.clone(),
        SAOM_FIT_AUX_OBSERVED => f.aux_observed.clone(),
        SAOM_FIT_META => vec![
            f.conditional as u8 as f64,
            periods as f64,
            f.effects as f64,
            f.iterations[0] as f64,
            f.iterations[1] as f64,
            f.iterations[2] as f64,
            p as f64,
            aux_len(n) as f64,
        ],
        _ => return None,
    })
}

/// Length of fit item `which`, or -1 if the item does not exist.
///
/// # Safety
/// `h` must be null or a live fit handle.
#[no_mangle]
pub unsafe extern "C" fn net_saom_fit_length(h: *const SaomFit, which: i32) -> i64 {
    match fit(h) {
        Ok(f) => item(f, which).map(|v| v.len() as i64).unwrap_or(-1),
        Err(_) => -1,
    }
}

/// Copies fit item `which` into `out`.
///
/// # Safety
/// `h` must be a live fit handle; `out` must hold as many doubles as
/// `net_saom_fit_length` reports.
#[no_mangle]
pub unsafe extern "C" fn net_saom_fit_copy(
    h: *const SaomFit,
    which: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = fit(h)?;
        let v = item(f, which)
            .ok_or_else(|| Error::InvalidArgument(format!("fit item {which} does not exist")))?;
        slice_mut(out, v.len()).copy_from_slice(&v);
        Ok(())
    })
}

/// Releases a fit. Null is ignored.
///
/// # Safety
/// `h` must be null or a live fit handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_saom_fit_free(h: *mut SaomFit) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}
