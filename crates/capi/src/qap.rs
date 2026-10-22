//! C ABI for the quadratic assignment procedure (`engine::qap`).
//!
//! Matrices are passed column-major, as R stores them. Counts are returned
//! as doubles, which hold integers exactly up to 2^53.

use std::ffi::c_char;

use engine::matrix::Dense;
use engine::qap::{
    mrqap, qap_cor, MrqapMethod, MrqapOptions, MrqapStatistic, Permutations, QapCorOptions,
    QapMode, Tail,
};
use engine::Error;

use super::{guard, slice, slice_mut};

/// Number of values `net_qap_cor` writes to `summary`.
pub const QAP_COR_SUMMARY_LEN: usize = 7;
/// Number of values per column that `net_mrqap` writes to `coef`.
pub const MRQAP_COEF_COLS: usize = 7;
/// Number of values `net_mrqap` writes to `fit`.
pub const MRQAP_FIT_LEN: usize = 10;

fn permutations(exhaustive: i32, reps: i64, seed: u64) -> Result<Permutations, Error> {
    if exhaustive != 0 {
        return Ok(Permutations::Exhaustive);
    }
    if reps < 0 {
        return Err(Error::InvalidArgument("reps must be non-negative".into()));
    }
    Ok(Permutations::Random {
        reps: reps as usize,
        seed,
    })
}

fn dims(n: i64) -> Result<usize, Error> {
    if n < 0 {
        return Err(Error::InvalidArgument("negative matrix size".into()));
    }
    Ok(n as usize)
}

fn tail_values(t: &Tail) -> [f64; 4] {
    [t.le as f64, t.ge as f64, t.abs_ge as f64, t.valid as f64]
}

/// QAP correlation test between two n-by-n matrices (column-major).
///
/// `mode`: 0 = digraph (all ordered pairs), 1 = graph (lower triangle).
/// With `exhaustive != 0` every one of the n! permutations is evaluated
/// (n <= 10) and `reps` and `seed` are ignored; otherwise `reps` random
/// permutations are drawn, draw i from stream i of `seed`.
///
/// `summary` receives 7 values: the observed correlation, the number of
/// dyads used, the number of draws, and the counts of draws with
/// r <= observed, r >= observed, |r| >= |observed|, and with r defined.
/// `dist` may be null; otherwise it receives every draw's correlation and
/// must hold `reps` doubles (n! when exhaustive).
///
/// # Safety
/// `y` and `x` must hold `n * n` doubles, `summary` 7 doubles, and `dist`
/// (if not null) the number of draws.
#[no_mangle]
pub unsafe extern "C" fn net_qap_cor(
    y: *const f64,
    x: *const f64,
    n: i64,
    mode: i32,
    diag: i32,
    exhaustive: i32,
    reps: i64,
    seed: u64,
    threads: i32,
    summary: *mut f64,
    dist: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = dims(n)?;
        let ym = Dense::from_col_major(n, n, slice(y, n * n))?;
        let xm = Dense::from_col_major(n, n, slice(x, n * n))?;
        let opts = QapCorOptions {
            mode: QapMode::from_code(mode)?,
            diag: diag != 0,
            permutations: permutations(exhaustive, reps, seed)?,
            threads: threads.max(0) as usize,
            keep_draws: !dist.is_null(),
        };
        let res = qap_cor(&ym, &xm, opts)?;
        let t = tail_values(&res.tail);
        slice_mut(summary, QAP_COR_SUMMARY_LEN).copy_from_slice(&[
            res.statistic,
            res.dyads as f64,
            res.draws as f64,
            t[0],
            t[1],
            t[2],
            t[3],
        ]);
        if !dist.is_null() {
            slice_mut(dist, res.distribution.len()).copy_from_slice(&res.distribution);
        }
        Ok(())
    })
}

/// Multiple regression QAP of `y` on `k` predictor matrices.
///
/// `x` holds the k predictors one after another, each n-by-n column-major.
/// The design has p = k + (intercept != 0) columns, the intercept first.
/// `method`: 0 = double semi-partialling, 1 = Y-permutation. `statistic`:
/// 0 = t-statistic, 1 = coefficient. `mode`, `exhaustive`, `reps`, and
/// `seed` are as in `net_qap_cor`.
///
/// `coef` receives a p-by-7 column-major table with columns: coefficient,
/// standard error, t-statistic, and the counts of draws with
/// stat <= observed, stat >= observed, |stat| >= |observed|, and with a
/// defined statistic. `fit` receives 10 values: dyads used, residual
/// degrees of freedom, residual sum of squares, model sum of squares,
/// R-squared, adjusted R-squared, residual standard error, F statistic,
/// the method used (0 or 1; a one-column design always uses 1), and the
/// number of draws. `dist` may be null; otherwise it receives a
/// draws-by-p column-major table of every draw's statistics (NaN where
/// undefined).
///
/// # Safety
/// `y` must hold `n * n` doubles, `x` `k * n * n`, `coef` `7 * p`, `fit` 10,
/// and `dist` (if not null) `draws * p`.
#[no_mangle]
pub unsafe extern "C" fn net_mrqap(
    y: *const f64,
    x: *const f64,
    n: i64,
    k: i64,
    intercept: i32,
    mode: i32,
    diag: i32,
    method: i32,
    statistic: i32,
    tol: f64,
    exhaustive: i32,
    reps: i64,
    seed: u64,
    threads: i32,
    coef: *mut f64,
    fit: *mut f64,
    dist: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = dims(n)?;
        let k = dims(k)?;
        let nn = n * n;
        let ym = Dense::from_col_major(n, n, slice(y, nn))?;
        let xs = slice(x, k * nn);
        let xm = (0..k)
            .map(|i| Dense::from_col_major(n, n, &xs[i * nn..(i + 1) * nn]))
            .collect::<Result<Vec<_>, _>>()?;
        let opts = MrqapOptions {
            mode: QapMode::from_code(mode)?,
            diag: diag != 0,
            intercept: intercept != 0,
            method: MrqapMethod::from_code(method)?,
            statistic: MrqapStatistic::from_code(statistic)?,
            permutations: permutations(exhaustive, reps, seed)?,
            tol,
            threads: threads.max(0) as usize,
            keep_draws: !dist.is_null(),
        };
        let res = mrqap(&ym, &xm, opts)?;
        let p = res.coefficients.len();
        let table = slice_mut(coef, MRQAP_COEF_COLS * p);
        for j in 0..p {
            let t = tail_values(&res.tails[j]);
            let row = [
                res.coefficients[j],
                res.std_errors[j],
                res.t[j],
                t[0],
                t[1],
                t[2],
                t[3],
            ];
            for (c, v) in row.iter().enumerate() {
                table[c * p + j] = *v;
            }
        }
        slice_mut(fit, MRQAP_FIT_LEN).copy_from_slice(&[
            res.dyads as f64,
            res.df_residual as f64,
            res.rss,
            res.mss,
            res.r_squared,
            res.adj_r_squared,
            res.sigma,
            res.f_statistic,
            match res.method {
                MrqapMethod::DoubleSemiPartialling => 0.0,
                MrqapMethod::YPermutation => 1.0,
            },
            res.draws as f64,
        ]);
        if !dist.is_null() {
            let d = res.draws as usize;
            let out = slice_mut(dist, d * p);
            for i in 0..d {
                for j in 0..p {
                    out[j * d + i] = res.distribution[i * p + j];
                }
            }
        }
        Ok(())
    })
}
