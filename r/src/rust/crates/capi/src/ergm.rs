//! C ABI for exponential random graph models (`engine::ergm`).
//!
//! A model is built once with `net_ergm_new` from a network, its
//! covariates, and a list of terms, and then queried through an opaque
//! `ErgmHandle`. Simulated networks are kept in an `ErgmSample`.
//!
//! Term codes (`term_code`), with the meaning of `term_param` and
//! `term_attr`:
//!
//! | code | term | param | attr |
//! |---|---|---|---|
//! | 0 | edges | | |
//! | 1 | mutual | | |
//! | 2, 3, 4 | nodecov, nodeicov, nodeocov | | numeric attribute |
//! | 5, 6, 7 | nodefactor, nodeifactor, nodeofactor | | categorical attribute; levels |
//! | 8 | nodematch | | categorical attribute |
//! | 9 | nodematch, diff = TRUE | | categorical attribute; levels |
//! | 10 | absdiff | exponent | numeric attribute |
//! | 11 | edgecov | | dyadic covariate |
//! | 12, 13, 14 | kstar, istar, ostar | k | |
//! | 15, 16, 17 | triangle, ttriple, ctriple | | |
//! | 18, 19, 20 | gwesp, gwdsp, gwnsp | decay | |
//! | 21, 22, 23 | gwdegree, gwidegree, gwodegree | decay | |
//! | 24, 25, 26 | degree, idegree, odegree | k | |
//! | 27 | isolates | | |
//!
//! Levels of factor terms are the slice `term_levels[off .. off + len]`
//! given by `term_levels_off` and `term_levels_len`. Numeric attributes are
//! `k_numeric` blocks of `n` doubles; categorical attributes `k_categorical`
//! blocks of `n` level codes; dyadic covariates `k_dyadic` column-major
//! n-by-n matrices. `fixed` may be null, or n*n row-major flags marking
//! dyads that keep their observed state. Matrices written to the caller
//! are column-major, as R stores them.

use std::ffi::c_char;

use engine::ergm::exact::Enumeration;
use engine::ergm::gof::{gof_stat, GofStat};
use engine::ergm::logit::{logit, LogitOptions};
use engine::ergm::mcmle::{log_likelihood_ratio, mcmle, McmleOptions};
use engine::ergm::mple::mple;
use engine::ergm::tergm::pooled_bootstrap;
use engine::ergm::{
    simulate, Constraints, Covariates, DegreeKind, Endpoint, Model, NetState, Partners, Proposal,
    Sample, SampleOptions, Term, Terms,
};
use engine::matrix::Dense;
use engine::Error;

use crate::{count, deref, guard, slice, slice_mut, thread_count};

/// Opaque model handle.
pub struct ErgmHandle(Model);

/// Opaque sample of networks.
pub struct ErgmSample(Sample);

/// Number of values `net_ergm_mple` writes to `fit`.
pub const ERGM_MPLE_FIT_LEN: usize = 6;
/// Number of values `net_logit` and `net_tergm_mple` write to `fit`.
pub const LOGIT_FIT_LEN: usize = 5;
/// Number of values `net_ergm_mcmle` writes to `info`.
pub const ERGM_MCMLE_INFO_LEN: usize = 5;

fn model<'a>(h: *const ErgmHandle) -> Result<&'a Model, Error> {
    // SAFETY: non-null handles are only created by `net_ergm_new`.
    unsafe { deref(h, "model") }.map(|m| &m.0)
}

fn sample<'a>(s: *const ErgmSample) -> Result<&'a Sample, Error> {
    // SAFETY: non-null handles are only created by `net_ergm_simulate`.
    unsafe { deref(s, "sample") }.map(|s| &s.0)
}

fn attr(a: i32) -> Result<usize, Error> {
    if a < 0 {
        return Err(Error::InvalidArgument(
            "a term that needs an attribute was given none".into(),
        ));
    }
    Ok(a as usize)
}

fn whole(v: f64, what: &str) -> Result<usize, Error> {
    if !v.is_finite() || v < 0.0 || v.fract() != 0.0 {
        return Err(Error::InvalidArgument(format!(
            "{what} must be a non-negative whole number"
        )));
    }
    Ok(v as usize)
}

fn decode_term(code: i32, param: f64, attr_index: i32, levels: &[i32]) -> Result<Term, Error> {
    let levels: Vec<u32> = levels
        .iter()
        .map(|&l| {
            if l < 0 {
                Err(Error::InvalidArgument(
                    "level codes must be non-negative".into(),
                ))
            } else {
                Ok(l as u32)
            }
        })
        .collect::<Result<_, _>>()?;
    let endpoint = |c: i32| match c {
        0 => Endpoint::Both,
        1 => Endpoint::Head,
        _ => Endpoint::Tail,
    };
    let kind = |c: i32| match c {
        0 => DegreeKind::Total,
        1 => DegreeKind::In,
        _ => DegreeKind::Out,
    };
    Ok(match code {
        0 => Term::Edges,
        1 => Term::Mutual,
        2..=4 => Term::NodeCov {
            attr: attr(attr_index)?,
            endpoint: endpoint(code - 2),
        },
        5..=7 => Term::NodeFactor {
            attr: attr(attr_index)?,
            endpoint: endpoint(code - 5),
            levels,
        },
        8 => Term::NodeMatch {
            attr: attr(attr_index)?,
            levels: Vec::new(),
        },
        9 => Term::NodeMatch {
            attr: attr(attr_index)?,
            levels,
        },
        10 => Term::AbsDiff {
            attr: attr(attr_index)?,
            pow: param,
        },
        11 => Term::EdgeCov {
            cov: attr(attr_index)?,
        },
        12..=14 => Term::Star {
            k: whole(param, "the star size")?,
            kind: kind(code - 12),
        },
        15 => Term::Triangle,
        16 => Term::TransitiveTriple,
        17 => Term::CyclicTriple,
        18..=20 => Term::GwPartners {
            decay: param,
            partners: match code {
                18 => Partners::Edgewise,
                19 => Partners::Dyadwise,
                _ => Partners::Nonedgewise,
            },
        },
        21..=23 => Term::GwDegree {
            decay: param,
            kind: kind(code - 21),
        },
        24..=26 => Term::Degree {
            k: whole(param, "the degree")?,
            kind: kind(code - 24),
        },
        27 => Term::Isolates,
        _ => return Err(Error::InvalidArgument(format!("unknown term code {code}"))),
    })
}

#[allow(clippy::too_many_arguments)]
fn sample_options(
    burnin: i64,
    interval: i64,
    nsim: i64,
    chains: i32,
    proposal: i32,
    seed: u64,
    threads: i32,
) -> Result<SampleOptions, Error> {
    Ok(SampleOptions {
        burnin: count(burnin, "burnin")?,
        interval: count(interval, "interval")?,
        nsim: count(nsim, "nsim")?,
        chains: count(chains as i64, "chains")?,
        proposal: Proposal::from_code(proposal)?,
        seed,
        stream_offset: 0,
        threads: thread_count(threads),
        keep_networks: false,
        start_empty: false,
    })
}

fn write_matrix(out: *mut f64, m: &Dense) {
    // SAFETY: the caller guarantees `out` holds rows * cols doubles.
    unsafe { slice_mut(out, m.rows() * m.cols()) }.copy_from_slice(&m.to_col_major());
}

/// Writes `rows` rows of `width` values stored row-major in `data` as a
/// column-major matrix.
fn write_rows(out: *mut f64, data: &[f64], rows: usize, width: usize) {
    // SAFETY: the caller guarantees `out` holds rows * width doubles.
    let o = unsafe { slice_mut(out, rows * width) };
    for r in 0..rows {
        for c in 0..width {
            o[c * rows + r] = data[r * width + c];
        }
    }
}

/// Builds a model. See the module documentation for the encoding of terms
/// and covariates. On success `*out` receives a handle that must be
/// released with `net_ergm_free`.
///
/// # Safety
/// `src` and `tgt` must hold `m` elements; the term arrays `nterms`;
/// `term_levels` every slice they reference; `numeric` `k_numeric * n`
/// doubles; `categorical` `k_categorical * n` codes; `dyadic`
/// `k_dyadic * n * n` doubles; `fixed` null or `n * n` flags; `out` must be
/// writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_ergm_new(
    n: i64,
    src: *const i64,
    tgt: *const i64,
    m: i64,
    directed: i32,
    term_code: *const i32,
    term_param: *const f64,
    term_attr: *const i32,
    term_levels_off: *const i32,
    term_levels_len: *const i32,
    term_levels: *const i32,
    nterms: i32,
    numeric: *const f64,
    k_numeric: i32,
    categorical: *const i32,
    k_categorical: i32,
    dyadic: *const f64,
    k_dyadic: i32,
    fixed: *const i32,
    out: *mut *mut ErgmHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "node count")?;
        let m = count(m, "edge count")?;
        let directed = directed != 0;
        let state = NetState::from_edges(n, slice(src, m), slice(tgt, m), directed)?;
        let nterms = count(nterms as i64, "term count")?;
        let codes = slice(term_code, nterms);
        let params = slice(term_param, nterms);
        let attrs = slice(term_attr, nterms);
        let offs = slice(term_levels_off, nterms);
        let lens = slice(term_levels_len, nterms);
        let total_levels: usize = lens.iter().map(|&l| l.max(0) as usize).sum();
        let max_end = offs
            .iter()
            .zip(lens)
            .map(|(&o, &l)| (o.max(0) + l.max(0)) as usize)
            .max()
            .unwrap_or(0);
        let levels = slice(term_levels, max_end.max(total_levels));
        let mut terms = Vec::with_capacity(nterms);
        for t in 0..nterms {
            let (o, l) = (offs[t].max(0) as usize, lens[t].max(0) as usize);
            terms.push(decode_term(
                codes[t],
                params[t],
                attrs[t],
                &levels[o..o + l],
            )?);
        }
        let kn = count(k_numeric as i64, "numeric attribute count")?;
        let kc = count(k_categorical as i64, "categorical attribute count")?;
        let kd = count(k_dyadic as i64, "dyadic covariate count")?;
        let num = slice(numeric, kn * n);
        let cat = slice(categorical, kc * n);
        let dy = slice(dyadic, kd * n * n);
        let cov = Covariates {
            numeric: (0..kn).map(|k| num[k * n..(k + 1) * n].to_vec()).collect(),
            categorical: (0..kc)
                .map(|k| {
                    cat[k * n..(k + 1) * n]
                        .iter()
                        .map(|&c| {
                            if c < 0 {
                                Err(Error::InvalidArgument(
                                    "categorical attributes must not have missing values".into(),
                                ))
                            } else {
                                Ok(c as u32)
                            }
                        })
                        .collect::<Result<Vec<u32>, _>>()
                })
                .collect::<Result<_, _>>()?,
            dyadic: (0..kd)
                .map(|k| Dense::from_col_major(n, n, &dy[k * n * n..(k + 1) * n * n]))
                .collect::<Result<_, _>>()?,
        };
        let terms = Terms::new(terms, n, directed, cov)?;
        let constraints = if fixed.is_null() {
            None
        } else {
            Some(Constraints::new(n, directed, slice(fixed, n * n))?)
        };
        let model = Model::new(state, terms, constraints)?;
        *out = Box::into_raw(Box::new(ErgmHandle(model)));
        Ok(())
    })
}

/// Releases a model handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_free(h: *mut ErgmHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Number of statistics, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_width(h: *const ErgmHandle) -> i64 {
    model(h).map(|m| m.width() as i64).unwrap_or(-1)
}

/// Number of free dyads, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_dyads(h: *const ErgmHandle) -> i64 {
    model(h).map(|m| m.free_dyad_count() as i64).unwrap_or(-1)
}

/// 1 if every term is dyad independent, 0 if not, -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_dyad_independent(h: *const ErgmHandle) -> i32 {
    model(h)
        .map(|m| m.terms().dyad_independent() as i32)
        .unwrap_or(-1)
}

/// The observed statistics; `out` holds `net_ergm_width` doubles.
///
/// # Safety
/// `h` must be a live handle; `out` must hold the statistics.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_stats(
    h: *const ErgmHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let s = m.stats();
        slice_mut(out, s.len()).copy_from_slice(&s);
        Ok(())
    })
}

/// Change statistics of every free dyad of the observed network: `x` is a
/// dyads-by-width column-major matrix, `y` the observed dyad states, and
/// `tail` and `head` the dyads' 0-based endpoints, in canonical order.
///
/// # Safety
/// `h` must be a live handle; `x` must hold `dyads * width` doubles and
/// `y`, `tail`, `head` `dyads` values each.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_changestats(
    h: *const ErgmHandle,
    x: *mut f64,
    y: *mut f64,
    tail: *mut i64,
    head: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let (xm, ys, dyads) = m.change_matrix();
        write_matrix(x, &xm);
        slice_mut(y, ys.len()).copy_from_slice(&ys);
        let t = slice_mut(tail, dyads.len());
        let hd = slice_mut(head, dyads.len());
        for (k, &(a, b)) in dyads.iter().enumerate() {
            t[k] = a as i64;
            hd[k] = b as i64;
        }
        Ok(())
    })
}

fn logit_opts(max_iter: i32, tol: f64) -> LogitOptions {
    let d = LogitOptions::default();
    LogitOptions {
        max_iter: if max_iter > 0 {
            max_iter as usize
        } else {
            d.max_iter
        },
        tol: if tol > 0.0 { tol } else { d.tol },
    }
}

/// Logistic regression of `y` (0/1) on the `m`-by-`p` column-major design
/// `x` with optional case weights `w` (null for unit weights), without an
/// intercept unless the design has one. `max_iter <= 0` and `tol <= 0` use
/// the defaults (50, 1e-12). `coef` receives p values, `cov` the p-by-p
/// covariance, and `fit` 5 values: deviance, null deviance (every
/// coefficient zero), total weight, iterations, and converged (0/1).
///
/// # Safety
/// `x` must hold `m * p` doubles, `y` (and `w` unless null) `m`, `coef`
/// `p`, `cov` `p * p`, and `fit` 5.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_logit(
    x: *const f64,
    y: *const f64,
    w: *const f64,
    m: i64,
    p: i64,
    max_iter: i32,
    tol: f64,
    coef: *mut f64,
    cov: *mut f64,
    fit: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = count(m, "row count")?;
        let p = count(p, "column count")?;
        let xm = Dense::from_col_major(m, p, slice(x, m * p))?;
        let ys = slice(y, m);
        let ws = if w.is_null() { None } else { Some(slice(w, m)) };
        let f = logit(&xm, ys, ws, logit_opts(max_iter, tol))?;
        slice_mut(coef, p).copy_from_slice(&f.coef);
        write_matrix(cov, &f.cov);
        slice_mut(fit, LOGIT_FIT_LEN).copy_from_slice(&[
            f.deviance,
            f.null_deviance,
            f.weight,
            f.iterations as f64,
            f.converged as i32 as f64,
        ]);
        Ok(())
    })
}

/// Maximum pseudolikelihood estimate. `fit` receives 6 values: deviance,
/// null deviance, number of dyads, iterations, converged (0/1), and whether
/// the model is dyad independent (0/1), in which case the estimate is the
/// maximum likelihood estimate.
///
/// # Safety
/// `h` must be a live handle; `coef` must hold `width` doubles, `cov`
/// `width * width`, and `fit` 6.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_mple(
    h: *const ErgmHandle,
    max_iter: i32,
    tol: f64,
    coef: *mut f64,
    cov: *mut f64,
    fit: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let r = mple(m, logit_opts(max_iter, tol))?;
        slice_mut(coef, m.width()).copy_from_slice(&r.fit.coef);
        write_matrix(cov, &r.fit.cov);
        slice_mut(fit, ERGM_MPLE_FIT_LEN).copy_from_slice(&[
            r.fit.deviance,
            r.fit.null_deviance,
            r.dyads as f64,
            r.fit.iterations as f64,
            r.fit.converged as i32 as f64,
            r.dyad_independent as i32 as f64,
        ]);
        Ok(())
    })
}

/// Simulates `chains * nsim` networks at `theta`: `burnin` proposals, then
/// `nsim` samples `interval` proposals apart, per chain; chain c draws from
/// stream c of `seed`. `proposal`: 0 = tie/no-tie, 1 = random dyad.
/// `start_empty != 0` starts every chain from the empty network. `stats`
/// receives a `(chains * nsim)`-by-width column-major matrix of sampled
/// statistics and `rate` two values: proposals accepted and made. If
/// `out` is not null, the sampled networks are kept in a sample handle
/// (`keep_networks != 0`) to be released with `net_ergm_sample_free`.
///
/// # Safety
/// `h` must be a live handle; `theta` must hold `width` doubles, `stats`
/// `chains * nsim * width`, `rate` 2; `out` must be null or writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_ergm_simulate(
    h: *const ErgmHandle,
    theta: *const f64,
    burnin: i64,
    interval: i64,
    nsim: i64,
    chains: i32,
    proposal: i32,
    seed: u64,
    threads: i32,
    keep_networks: i32,
    start_empty: i32,
    stats: *mut f64,
    rate: *mut f64,
    out: *mut *mut ErgmSample,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let mut opts = sample_options(burnin, interval, nsim, chains, proposal, seed, threads)?;
        opts.keep_networks = keep_networks != 0 && !out.is_null();
        opts.start_empty = start_empty != 0;
        let s = simulate(m, slice(theta, m.width()), &opts)?;
        write_rows(stats, &s.stats, s.rows(), s.width);
        slice_mut(rate, 2).copy_from_slice(&[s.accepted as f64, s.proposed as f64]);
        if !out.is_null() {
            *out = Box::into_raw(Box::new(ErgmSample(s)));
        }
        Ok(())
    })
}

/// Releases a sample handle. Null is ignored.
///
/// # Safety
/// `s` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_sample_free(s: *mut ErgmSample) {
    if !s.is_null() {
        drop(Box::from_raw(s));
    }
}

/// Number of networks in a sample (0 if networks were not kept), or -1 for
/// a null handle.
///
/// # Safety
/// `s` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_sample_count(s: *const ErgmSample) -> i64 {
    sample(s)
        .map(|x| x.networks.as_ref().map(|n| n.len()).unwrap_or(0) as i64)
        .unwrap_or(-1)
}

/// Number of ties in network `k` of a sample, or -1 if unavailable.
///
/// # Safety
/// `s` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_sample_ties(s: *const ErgmSample, k: i64) -> i64 {
    sample(s)
        .ok()
        .and_then(|x| x.networks.as_ref())
        .and_then(|n| n.get(k.max(0) as usize))
        .map(|t| t.len() as i64)
        .unwrap_or(-1)
}

/// The ties of network `k` as 0-based `tail` and `head` arrays of
/// `net_ergm_sample_ties(s, k)` elements.
///
/// # Safety
/// `s` must be a live handle; `tail` and `head` must hold the ties.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_sample_edges(
    s: *const ErgmSample,
    k: i64,
    tail: *mut i64,
    head: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let x = sample(s)?;
        let nets = x
            .networks
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("the sample did not keep networks".into()))?;
        let ties = nets
            .get(count(k, "network index")?)
            .ok_or_else(|| Error::InvalidArgument(format!("network index {k} is out of range")))?;
        let t = slice_mut(tail, ties.len());
        let hd = slice_mut(head, ties.len());
        for (i, &(a, b)) in ties.iter().enumerate() {
            t[i] = a as i64;
            hd[i] = b as i64;
        }
        Ok(())
    })
}

/// A goodness-of-fit statistic (`which`: 0 = degree, 1 = in-degree,
/// 2 = out-degree, 3 = edgewise shared partners, 4 = dyadwise shared
/// partners, 5 = geodesic distance, 6 = triad census) of every network in
/// a sample: a count-by-len column-major matrix, where len is n for
/// degrees and distance, n - 1 for shared partners, and 16 for the census.
///
/// # Safety
/// `s` must be a live handle whose networks were kept; `out` must hold
/// `count * len` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_gof_sample(
    s: *const ErgmSample,
    which: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let x = sample(s)?;
        let nets = x
            .networks
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("the sample did not keep networks".into()))?;
        let which = GofStat::from_code(which)?;
        let len = which.len(x.n);
        let o = slice_mut(out, nets.len() * len);
        for (r, net) in nets.iter().enumerate() {
            let v = gof_stat(x.n, x.directed, net, which)?;
            for (c, val) in v.iter().enumerate() {
                o[c * nets.len() + r] = *val;
            }
        }
        Ok(())
    })
}

/// A goodness-of-fit statistic of the observed network (codes as in
/// `net_ergm_gof_sample`); `out` holds len doubles.
///
/// # Safety
/// `h` must be a live handle; `out` must hold the statistic.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_gof_observed(
    h: *const ErgmHandle,
    which: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let which = GofStat::from_code(which)?;
        let y = m.state();
        let v = gof_stat(y.node_count(), y.is_directed(), y.ties(), which)?;
        slice_mut(out, v.len()).copy_from_slice(&v);
        Ok(())
    })
}

/// Monte Carlo maximum likelihood from `theta0`, with the sampling
/// settings of `net_ergm_simulate` for each iteration, at most `max_iter`
/// iterations, hull margin `margin`, and chi-squared critical value
/// `critical` for the convergence test. `theta` receives the estimate,
/// `cov` the inverse Fisher information, `mcmc_cov` the Monte Carlo term
/// of the covariance, `info` 5 values (iterations, converged, last step
/// length, convergence statistic, rows of `stats`), and `stats` the final
/// sample as a `(chains * nsim)`-by-width column-major matrix.
///
/// # Safety
/// `h` must be a live handle; `theta0` and `theta` must hold `width`
/// doubles, `cov` and `mcmc_cov` `width * width`, `info` 5, and `stats`
/// `chains * nsim * width`.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_ergm_mcmle(
    h: *const ErgmHandle,
    theta0: *const f64,
    burnin: i64,
    interval: i64,
    nsim: i64,
    chains: i32,
    proposal: i32,
    seed: u64,
    threads: i32,
    max_iter: i32,
    margin: f64,
    critical: f64,
    theta: *mut f64,
    cov: *mut f64,
    mcmc_cov: *mut f64,
    info: *mut f64,
    stats: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let opts = McmleOptions {
            sample: sample_options(burnin, interval, nsim, chains, proposal, seed, threads)?,
            max_iter: count(max_iter as i64, "max_iter")?,
            margin,
            critical,
        };
        let f = mcmle(m, slice(theta0, m.width()), &opts)?;
        slice_mut(theta, m.width()).copy_from_slice(&f.theta);
        write_matrix(cov, &f.cov);
        write_matrix(mcmc_cov, &f.mcmc_cov);
        slice_mut(info, ERGM_MCMLE_INFO_LEN).copy_from_slice(&[
            f.iterations as f64,
            f.converged as i32 as f64,
            f.step_length,
            f.statistic,
            f.sample.rows() as f64,
        ]);
        write_rows(stats, &f.sample.stats, f.sample.rows(), f.sample.width);
        Ok(())
    })
}

/// Exact log-likelihood at `theta` by enumeration (at most 22 free dyads),
/// with the exact expectation (`mean`, width values) and covariance (`cov`,
/// width-by-width) of the statistics.
///
/// # Safety
/// `h` must be a live handle; `theta` and `mean` must hold `width`
/// doubles, `loglik` 1, and `cov` `width * width`.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_exact(
    h: *const ErgmHandle,
    theta: *const f64,
    loglik: *mut f64,
    mean: *mut f64,
    cov: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let e = Enumeration::new(m)?;
        let th = slice(theta, m.width());
        *loglik = e.log_likelihood(th)?;
        let (mu, c) = e.moments(th)?;
        slice_mut(mean, m.width()).copy_from_slice(&mu);
        write_matrix(cov, &c);
        Ok(())
    })
}

/// Exact maximum likelihood estimate by enumeration and Newton's method
/// from `start`. `cov` receives the inverse Fisher information and `info`
/// two values: the log-likelihood at the estimate and the iterations used.
///
/// # Safety
/// `h` must be a live handle; `start` and `theta` must hold `width`
/// doubles, `cov` `width * width`, and `info` 2.
#[no_mangle]
pub unsafe extern "C" fn net_ergm_exact_mle(
    h: *const ErgmHandle,
    start: *const f64,
    max_iter: i32,
    theta: *mut f64,
    cov: *mut f64,
    info: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let e = Enumeration::new(m)?;
        let f = e.mle(
            slice(start, m.width()),
            count(max_iter as i64, "max_iter")?.max(1),
        )?;
        slice_mut(theta, m.width()).copy_from_slice(&f.theta);
        write_matrix(cov, &f.cov);
        slice_mut(info, 2).copy_from_slice(&[f.log_likelihood, f.iterations as f64]);
        Ok(())
    })
}

/// Path-sampling estimate of the log-likelihood ratio l(theta) - l(0) with
/// `steps` midpoint evaluations, each sampled with the settings of
/// `net_ergm_simulate`. `out` receives the estimate and its Monte Carlo
/// standard error. l(0) is -D log 2 for D free dyads.
///
/// # Safety
/// `h` must be a live handle; `theta` must hold `width` doubles and `out` 2.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_ergm_loglik(
    h: *const ErgmHandle,
    theta: *const f64,
    steps: i32,
    burnin: i64,
    interval: i64,
    nsim: i64,
    chains: i32,
    proposal: i32,
    seed: u64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = model(h)?;
        let opts = sample_options(burnin, interval, nsim, chains, proposal, seed, threads)?;
        let (llr, se) = log_likelihood_ratio(
            m,
            slice(theta, m.width()),
            count(steps as i64, "steps")?,
            &opts,
        )?;
        slice_mut(out, 2).copy_from_slice(&[llr, se]);
        Ok(())
    })
}

/// Pooled pseudolikelihood of a temporal ERGM: logistic regression of `y`
/// on the `m`-by-`p` column-major design `x` (change statistics of every
/// transition stacked), with `period[r]` the 0-based transition of row r
/// among `periods`, and `reps` bootstrap replicates that resample
/// transitions with replacement (replicate r draws from stream r of
/// `seed`). `coef`, `cov`, and `fit` are as in `net_logit`; `boot`
/// receives a reps-by-p column-major matrix of replicate estimates (NaN
/// where a replicate did not converge).
///
/// # Safety
/// `x` must hold `m * p` doubles, `y` and `period` `m`, `coef` `p`, `cov`
/// `p * p`, `fit` 5, and `boot` `reps * p`.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_tergm_mple(
    x: *const f64,
    y: *const f64,
    period: *const i32,
    m: i64,
    p: i64,
    periods: i32,
    reps: i64,
    seed: u64,
    threads: i32,
    max_iter: i32,
    tol: f64,
    coef: *mut f64,
    cov: *mut f64,
    fit: *mut f64,
    boot: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let m = count(m, "row count")?;
        let p = count(p, "column count")?;
        let reps = count(reps, "reps")?;
        let xm = Dense::from_col_major(m, p, slice(x, m * p))?;
        let per: Vec<u32> = slice(period, m)
            .iter()
            .map(|&t| {
                if t < 0 {
                    Err(Error::InvalidArgument(
                        "period labels must be non-negative".into(),
                    ))
                } else {
                    Ok(t as u32)
                }
            })
            .collect::<Result<_, _>>()?;
        let f = pooled_bootstrap(
            &xm,
            slice(y, m),
            &per,
            count(periods as i64, "periods")?,
            reps,
            seed,
            thread_count(threads),
            logit_opts(max_iter, tol),
        )?;
        slice_mut(coef, p).copy_from_slice(&f.fit.coef);
        write_matrix(cov, &f.fit.cov);
        slice_mut(fit, LOGIT_FIT_LEN).copy_from_slice(&[
            f.fit.deviance,
            f.fit.null_deviance,
            f.fit.weight,
            f.fit.iterations as f64,
            f.fit.converged as i32 as f64,
        ]);
        write_rows(boot, &f.boot, reps, p);
        Ok(())
    })
}
