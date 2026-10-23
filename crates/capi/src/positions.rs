//! C ABI for brokerage and positional measures: Burt's structural holes
//! (constraint, effective size, efficiency, and hierarchy),
//! Gould-Fernandez brokerage, structural equivalence, and CONCOR.

use std::ffi::c_char;

use engine::metrics::brokerage::{
    brokerage, constraint, effective_size, efficiency, hierarchy, StructuralHoleOptions,
};
use engine::metrics::equivalence::{
    concor, structural_equivalence, ConcorOptions, EquivalenceMethod, EquivalenceOptions,
};
use engine::Error;

use crate::{guard, handle, slice, slice_mut, GraphHandle};

/// Burt's structural-hole measures. `measure`: 0 = constraint, 1 =
/// effective size, 2 = efficiency, 3 = hierarchy (NaN for nodes with
/// fewer than two contacts). With `weighted != 0`, edge weights are
/// tie strengths. Nodes without contacts get NaN. `threads == 0` uses all
/// cores.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_structural_holes(
    h: *const GraphHandle,
    measure: i32,
    weighted: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = StructuralHoleOptions {
            weighted: weighted != 0,
            threads: threads.max(0) as usize,
        };
        let values = match measure {
            0 => constraint(g, opts)?,
            1 => effective_size(g, opts)?,
            2 => efficiency(g, opts)?,
            3 => hierarchy(g, opts)?,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown structural hole measure code {measure}"
                )))
            }
        };
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Gould-Fernandez brokerage counts. `groups` holds `n_groups` (which must
/// equal the node count) group codes. `out` receives an n-by-6 matrix in
/// column-major order with columns coordinator, itinerant, representative,
/// gatekeeper, liaison, total.
///
/// # Safety
/// `h` must be a live handle; `groups` must hold `n_groups` values; `out`
/// must hold `6 * n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_brokerage(
    h: *const GraphHandle,
    groups: *const i64,
    n_groups: i64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if n_groups < 0 {
            return Err(Error::InvalidArgument("negative group count".into()));
        }
        let m = brokerage(g, slice(groups, n_groups as usize), threads.max(0) as usize)?;
        slice_mut(out, 6 * g.node_count()).copy_from_slice(&m.to_col_major());
        Ok(())
    })
}

/// Structural equivalence of every pair of nodes. `method`: 0 = Euclidean
/// distance, 1 = correlation, 2 = Hamming distance. With `weighted != 0`
/// ties are total edge weights. `out` receives an n-by-n matrix in
/// column-major order (it is symmetric).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n * n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_structural_equivalence(
    h: *const GraphHandle,
    method: i32,
    weighted: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = EquivalenceOptions {
            method: EquivalenceMethod::from_code(method)?,
            weighted: weighted != 0,
            threads: threads.max(0) as usize,
        };
        let m = structural_equivalence(g, opts)?;
        let n = g.node_count();
        slice_mut(out, n * n).copy_from_slice(&m.to_col_major());
        Ok(())
    })
}

/// CONCOR block membership (numbered from 1). `diagnostics` (may be null)
/// receives two counts: splits that reached `max_iter` without converging,
/// and splits in which a correlation with the reference node came within
/// 1e-6 of zero (so rounding may have decided a node's block).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` 32-bit integers;
/// `diagnostics` must be null or hold 2 writable 32-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_concor(
    h: *const GraphHandle,
    depth: i32,
    weighted: i32,
    tolerance: f64,
    max_iter: i32,
    threads: i32,
    out: *mut i32,
    diagnostics: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if depth < 1 || max_iter < 0 {
            return Err(Error::InvalidArgument(
                "CONCOR depth must be at least 1 and max_iter non-negative".into(),
            ));
        }
        let opts = ConcorOptions {
            depth: depth as usize,
            weighted: weighted != 0,
            tolerance,
            max_iter: max_iter as usize,
            threads: threads.max(0) as usize,
        };
        let r = concor(g, opts)?;
        let dst = slice_mut(out, g.node_count());
        for (d, &m) in dst.iter_mut().zip(&r.membership) {
            *d = m as i32;
        }
        if !diagnostics.is_null() {
            let d = slice_mut(diagnostics, 2);
            d[0] = r.unconverged as i32;
            d[1] = r.near_zero as i32;
        }
        Ok(())
    })
}
