//! C ABI for PageRank, Katz centrality, Bonacich power centrality, hubs
//! and authorities, assortativity, mixing matrices, the E-I index, and
//! within-group tie shares.
//!
//! Integer codes used by these functions:
//! - `mode` and `direction`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `loops`: 0 = an undirected self-loop counts twice (igraph), 1 = once
//!   (NetworkX), 2 = not at all.
//! - `dangling` (PageRank): 0 = jump by the personalization vector, 1 =
//!   jump uniformly.
//! - `scale` (hubs and authorities): 0 = largest score 1, 1 = sum 1, 2 =
//!   unit Euclidean norm.
//! - `scope` (E-I): 0 = network, 1 = groups, 2 = nodes.
//! - `counting` (E-I): 0 = Krackhardt and Stern, 1 = the R package
//!   homophily.
//! - Category codes are 0-based; -1 marks a missing value where allowed.
//!
//! `tol <= 0` and `max_iter <= 0` select the defaults (1e-12 and 10000).

use std::ffi::c_char;

use engine::metrics::degree::DegreeMode;
use engine::metrics::eigenvector::LoopWeight;
use engine::metrics::mixing::{
    assortativity_degree, assortativity_nominal, assortativity_values, ei_counts, group_shares,
    mixing_matrix, AssortOptions, EiCounting, EiScope,
};
use engine::metrics::spectral::{
    bonacich_power, hits, katz, pagerank, BonacichOptions, Dangling, HitsOptions, HitsScale,
    KatzOptions, PagerankOptions, DEFAULT_MAX_ITER, DEFAULT_TOL,
};
use engine::Error;

use crate::{guard, handle, slice, slice_mut, GraphHandle};

fn tolerance(tol: f64) -> f64 {
    if tol > 0.0 {
        tol
    } else {
        DEFAULT_TOL
    }
}

fn iterations(max_iter: i64) -> usize {
    if max_iter > 0 {
        max_iter as usize
    } else {
        DEFAULT_MAX_ITER
    }
}

fn count(k: i64, what: &str) -> Result<usize, Error> {
    usize::try_from(k).map_err(|_| Error::InvalidArgument(format!("negative {what}")))
}

/// PageRank. `personalization` may be null (uniform); otherwise it holds n
/// nonnegative values, rescaled to sum to 1. `directed == 0` ignores arc
/// directions.
///
/// # Safety
/// `h` must be a live handle; `personalization` must be null or hold n
/// doubles; `out` must hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_pagerank(
    h: *const GraphHandle,
    damping: f64,
    directed: i32,
    weighted: i32,
    personalization: *const f64,
    dangling: i32,
    loops: i32,
    tol: f64,
    max_iter: i64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let opts = PagerankOptions {
            damping,
            directed: directed != 0,
            weighted: weighted != 0,
            personalization: if personalization.is_null() {
                None
            } else {
                Some(slice(personalization, n).to_vec())
            },
            dangling: Dangling::from_code(dangling)?,
            loops: LoopWeight::from_code(loops)?,
            tol: tolerance(tol),
            max_iter: iterations(max_iter),
            threads: threads.max(0) as usize,
        };
        let values = pagerank(g, &opts)?;
        slice_mut(out, n).copy_from_slice(&values);
        Ok(())
    })
}

/// Katz centrality x = alpha A' x + beta. `beta` holds one value (for
/// every node) or n values (`beta_len`). `normalized != 0` scales to unit
/// Euclidean norm.
///
/// # Safety
/// `h` must be a live handle; `beta` must hold `beta_len` doubles; `out`
/// must hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_katz(
    h: *const GraphHandle,
    alpha: f64,
    beta: *const f64,
    beta_len: i64,
    mode: i32,
    weighted: i32,
    loops: i32,
    normalized: i32,
    tol: f64,
    max_iter: i64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = KatzOptions {
            alpha,
            beta: slice(beta, count(beta_len, "beta length")?).to_vec(),
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            loops: LoopWeight::from_code(loops)?,
            normalized: normalized != 0,
            tol: tolerance(tol),
            max_iter: iterations(max_iter),
            threads: threads.max(0) as usize,
        };
        let values = katz(g, &opts)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Bonacich's power centrality with attenuation `exponent`. `rescale != 0`
/// divides by the sum of the scores instead of scaling the sum of squares
/// to n.
///
/// # Safety
/// `h` must be a live handle; `out` must hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_bonacich_power(
    h: *const GraphHandle,
    exponent: f64,
    mode: i32,
    weighted: i32,
    loops: i32,
    rescale: i32,
    tol: f64,
    max_iter: i64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = BonacichOptions {
            exponent,
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            loops: LoopWeight::from_code(loops)?,
            rescale: rescale != 0,
            tol: tolerance(tol),
            max_iter: iterations(max_iter),
            threads: threads.max(0) as usize,
        };
        let values = bonacich_power(g, opts)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Hub and authority scores. On success `value` receives the largest
/// eigenvalue of A^T A and `unique` 1 if the leading singular value is
/// simple, else 0.
///
/// # Safety
/// `h` must be a live handle; `hubs` and `authorities` must each hold n
/// doubles; `value` and `unique` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_hits(
    h: *const GraphHandle,
    weighted: i32,
    loops: i32,
    scale: i32,
    tol: f64,
    max_iter: i64,
    threads: i32,
    hubs: *mut f64,
    authorities: *mut f64,
    value: *mut f64,
    unique: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let opts = HitsOptions {
            weighted: weighted != 0,
            loops: LoopWeight::from_code(loops)?,
            scale: HitsScale::from_code(scale)?,
            tol: tolerance(tol),
            max_iter: iterations(max_iter),
            threads: threads.max(0) as usize,
        };
        let r = hits(g, opts)?;
        slice_mut(hubs, n).copy_from_slice(&r.hubs);
        slice_mut(authorities, n).copy_from_slice(&r.authorities);
        slice_mut(value, 1)[0] = r.value;
        slice_mut(unique, 1)[0] = r.unique as i32;
        Ok(())
    })
}

fn assort(directed: i32, loops: i32, normalized: i32) -> Result<AssortOptions, Error> {
    Ok(AssortOptions {
        directed: directed != 0,
        loops: LoopWeight::from_code(loops)?,
        normalized: normalized != 0,
    })
}

/// Nominal assortativity of node categories `types` (n codes in 0..k).
/// Writes one double.
///
/// # Safety
/// `h` must be a live handle; `types` must hold n values; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_assortativity_nominal(
    h: *const GraphHandle,
    types: *const i64,
    k: i64,
    directed: i32,
    loops: i32,
    normalized: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = slice(types, g.node_count());
        let r = assortativity_nominal(
            g,
            t,
            count(k, "category count")?,
            assort(directed, loops, normalized)?,
        )?;
        slice_mut(out, 1)[0] = r;
        Ok(())
    })
}

/// Assortativity of numeric node values: `values` at arc sources and
/// `values_in` (null: the same values) at arc targets. Writes one double.
///
/// # Safety
/// `h` must be a live handle; `values` (and `values_in` unless null) must
/// hold n doubles; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_assortativity_values(
    h: *const GraphHandle,
    values: *const f64,
    values_in: *const f64,
    directed: i32,
    loops: i32,
    normalized: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let vin = if values_in.is_null() {
            None
        } else {
            Some(slice(values_in, n))
        };
        let r = assortativity_values(
            g,
            slice(values, n),
            vin,
            assort(directed, loops, normalized)?,
        )?;
        slice_mut(out, 1)[0] = r;
        Ok(())
    })
}

/// Degree assortativity (strengths if `weighted != 0`). Writes one double.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_assortativity_degree(
    h: *const GraphHandle,
    weighted: i32,
    directed: i32,
    loops: i32,
    normalized: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let r = assortativity_degree(g, weighted != 0, assort(directed, loops, normalized)?)?;
        slice_mut(out, 1)[0] = r;
        Ok(())
    })
}

/// Mixing matrix of categories `dim1` (codes 0..k1, -1 missing) by `dim2`
/// (0..k2). `same != 0` says the two are the same attribute. Writes
/// k1 * k2 counts row-major into `out`; `*transposed` is set to 1 when the
/// rows are `dim2` categories and the columns `dim1` (directed networks
/// with direction 2 = in), so the matrix is k2 by k1.
///
/// # Safety
/// `h` must be a live handle; `dim1` and `dim2` must hold n values; `out`
/// must hold k1 * k2 doubles; `transposed` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_mixing_matrix(
    h: *const GraphHandle,
    dim1: *const i64,
    k1: i64,
    dim2: *const i64,
    k2: i64,
    same: i32,
    direction: i32,
    out: *mut f64,
    transposed: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let (k1, k2) = (count(k1, "category count")?, count(k2, "category count")?);
        let m = mixing_matrix(
            g,
            slice(dim1, n),
            k1,
            slice(dim2, n),
            k2,
            same != 0,
            DegreeMode::from_code(direction)?,
        )?;
        slice_mut(out, k1 * k2).copy_from_slice(&m.counts);
        slice_mut(transposed, 1)[0] = m.transposed as i32;
        Ok(())
    })
}

/// E-I index with its external and internal tie counts, for the network
/// (one value), each group (k values), or each node (n values) as `scope`
/// says. `external` and `internal` may be null.
///
/// # Safety
/// `h` must be a live handle; `groups` must hold n codes in 0..k; `out`
/// (and `external`, `internal` unless null) must hold as many doubles as
/// the scope gives values.
#[no_mangle]
pub unsafe extern "C" fn net_ei_index(
    h: *const GraphHandle,
    groups: *const i64,
    k: i64,
    scope: i32,
    direction: i32,
    counting: i32,
    out: *mut f64,
    external: *mut f64,
    internal: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = ei_counts(
            g,
            slice(groups, g.node_count()),
            count(k, "group count")?,
            EiScope::from_code(scope)?,
            DegreeMode::from_code(direction)?,
            EiCounting::from_code(counting)?,
        )?;
        let dst = slice_mut(out, c.len());
        for (d, &(e, i)) in dst.iter_mut().zip(&c) {
            *d = if e + i == 0.0 {
                f64::NAN
            } else {
                (e - i) / (e + i)
            };
        }
        if !external.is_null() {
            for (d, &(e, _)) in slice_mut(external, c.len()).iter_mut().zip(&c) {
                *d = e;
            }
        }
        if !internal.is_null() {
            for (d, &(_, i)) in slice_mut(internal, c.len()).iter_mut().zip(&c) {
                *d = i;
            }
        }
        Ok(())
    })
}

/// Within-group tie shares: for each of the k groups, five values in
/// `out` (row-major k by 5): size, internal ties, external ties, observed
/// internal share, and expected share (size - 1) / (n - 1).
///
/// # Safety
/// `h` must be a live handle; `groups` must hold n codes in 0..k; `out`
/// must hold 5 k doubles.
#[no_mangle]
pub unsafe extern "C" fn net_group_shares(
    h: *const GraphHandle,
    groups: *const i64,
    k: i64,
    direction: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let k = count(k, "group count")?;
        let s = group_shares(
            g,
            slice(groups, g.node_count()),
            k,
            DegreeMode::from_code(direction)?,
        )?;
        let dst = slice_mut(out, 5 * k);
        for (c, share) in s.iter().enumerate() {
            dst[5 * c] = share.size as f64;
            dst[5 * c + 1] = share.internal;
            dst[5 * c + 2] = share.external;
            dst[5 * c + 3] = share.observed;
            dst[5 * c + 4] = share.expected;
        }
        Ok(())
    })
}
