//! C ABI for network layouts.
//!
//! Every function writes node coordinates into `x` and `y`, each holding n
//! doubles, aligned with node indices. Integer codes:
//! - `weights_as`: 0 = distance, 1 = strength (only used when `weighted`
//!   is non-zero).
//! - `algorithm` (stress): 0 = automatic by component size, 1 = all pairs,
//!   2 = sparse stress with pivots.
//! - `approximation` (force-directed): 0 = automatic by component size,
//!   1 = exact, 2 = Barnes-Hut.
//! - `kind` (simple layouts): 0 = circle, 1 = star, 2 = grid.
//!
//! Zero (or a negative tolerance or angle) selects the engine's default
//! for any iteration count, pivot count, tolerance, or angle.

use std::ffi::c_char;

use engine::layout::{
    circle_layout, fr_layout, grid_layout, layout_stress, star_layout, stress_layout, tree_layout,
    twomode_layout, FrApproximation, FrOptions, Layout, StressAlgorithm, StressOptions,
};
use engine::{Error, WeightKind};

use crate::{guard, handle, slice, slice_mut, GraphHandle};

fn count(value: i64, what: &str) -> Result<usize, Error> {
    usize::try_from(value)
        .map_err(|_| Error::InvalidArgument(format!("{what} must not be negative")))
}

/// # Safety
/// `x` and `y` must each hold `n` doubles.
unsafe fn write(l: &Layout, x: *mut f64, y: *mut f64) {
    let n = l.x.len();
    slice_mut(x, n).copy_from_slice(&l.x);
    slice_mut(y, n).copy_from_slice(&l.y);
}

/// Stress layout: stress majorization on all pairs, or sparse stress with
/// pivots for large components. `seed` drives the small perturbation of the
/// starting layout and the sample of sparse-stress pivots.
///
/// # Safety
/// `h` must be a live handle; `x` and `y` must each hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_layout_stress(
    h: *const GraphHandle,
    algorithm: i32,
    weighted: i32,
    weights_as: i32,
    pivots: i64,
    max_iter: i64,
    tol: f64,
    seed: u64,
    threads: i32,
    x: *mut f64,
    y: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = StressOptions {
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            algorithm: StressAlgorithm::from_code(algorithm)?,
            pivots: count(pivots, "the number of pivots")?,
            max_iter: count(max_iter, "the iteration cap")?,
            tol,
            seed,
            threads: threads.max(0) as usize,
        };
        write(&stress_layout(g, &opts)?, x, y);
        Ok(())
    })
}

/// Fruchterman-Reingold force-directed layout. `seed` drives the starting
/// positions.
///
/// # Safety
/// `h` must be a live handle; `x` and `y` must each hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_layout_fr(
    h: *const GraphHandle,
    weighted: i32,
    weights_as: i32,
    iterations: i64,
    approximation: i32,
    theta: f64,
    seed: u64,
    threads: i32,
    x: *mut f64,
    y: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = FrOptions {
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            iterations: count(iterations, "the number of iterations")?,
            approximation: FrApproximation::from_code(approximation)?,
            theta,
            seed,
            threads: threads.max(0) as usize,
        };
        write(&fr_layout(g, &opts)?, x, y);
        Ok(())
    })
}

/// Circle (`kind` 0), star (1), or grid (2). `param` is the star's center
/// (a node index, or -1 for the default) or the grid's width (0 for the
/// default); it is ignored for the circle.
///
/// # Safety
/// `h` must be a live handle; `x` and `y` must each hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_layout_simple(
    h: *const GraphHandle,
    kind: i32,
    param: i64,
    x: *mut f64,
    y: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let l = match kind {
            0 => circle_layout(g),
            1 => {
                let center = if param < 0 {
                    None
                } else {
                    Some(param as usize)
                };
                star_layout(g, center)?
            }
            2 => grid_layout(g, count(param, "the grid width")?),
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown simple layout code {kind}"
                )))
            }
        };
        write(&l, x, y);
        Ok(())
    })
}

/// Layered tree layout of each component. `roots` (may be null when
/// `n_roots` is 0) gives at most one root per component; other components
/// use their default root.
///
/// # Safety
/// `h` must be a live handle; `roots` must hold `n_roots` integers; `x` and
/// `y` must each hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_layout_tree(
    h: *const GraphHandle,
    roots: *const i64,
    n_roots: i64,
    x: *mut f64,
    y: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let mut r = Vec::new();
        for &v in slice(roots, count(n_roots, "the number of roots")?) {
            if v < 0 || v as u64 >= n as u64 {
                return Err(Error::InvalidArgument(format!(
                    "root {v} is not a node index (0..{n})"
                )));
            }
            r.push(v as usize);
        }
        write(&tree_layout(g, &r)?, x, y);
        Ok(())
    })
}

/// Two-row layout for a two-mode network: nodes with `types[v] == 0` on
/// top, the others below, ordered by `sweeps` barycenter sweeps (0 for the
/// default of 24).
///
/// # Safety
/// `h` must be a live handle; `types` must hold n integers; `x` and `y`
/// must each hold n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_layout_twomode(
    h: *const GraphHandle,
    types: *const i32,
    sweeps: i64,
    x: *mut f64,
    y: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let t: Vec<bool> = slice(types, n).iter().map(|&v| v != 0).collect();
        let sweeps = match count(sweeps, "the number of sweeps")? {
            0 => 24,
            s => s,
        };
        write(&twomode_layout(g, &t, sweeps)?, x, y);
        Ok(())
    })
}

/// The stress of the layout (`x`, `y`): the sum over pairs of nodes in the
/// same component of (||x_i - x_j|| - d_ij)^2 / d_ij^2, with shortest-path
/// lengths d_ij in the undirected skeleton.
///
/// # Safety
/// `h` must be a live handle; `x` and `y` must each hold n doubles; `out`
/// must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_layout_stress_value(
    h: *const GraphHandle,
    x: *const f64,
    y: *const f64,
    weighted: i32,
    weights_as: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let l = Layout {
            x: slice(x, n).to_vec(),
            y: slice(y, n).to_vec(),
        };
        let s = layout_stress(
            g,
            &l,
            weighted != 0,
            WeightKind::from_code(weights_as)?,
            threads.max(0) as usize,
        )?;
        slice_mut(out, 1)[0] = s;
        Ok(())
    })
}
