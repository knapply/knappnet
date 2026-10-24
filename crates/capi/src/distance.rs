//! C ABI for center, periphery, barycenter, the distance indices (Wiener,
//! Schultz, Gutman, hyper-Wiener, harmonic diameter, global efficiency),
//! local efficiency, and closeness vitality.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength.
//! - `reachable`: 0 = require a (strongly) connected network (NetworkX),
//!   1 = eccentricity over reached nodes (igraph).
//! - `index`: 0 = Wiener, 1 = Schultz, 2 = Gutman, 3 = hyper-Wiener
//!   (published), 4 = hyper-Wiener (NetworkX), 5 = harmonic diameter,
//!   6 = global efficiency.
//! - `tie`: 0 = relative tolerance 1e-10, 1 = exact equality.

use std::ffi::c_char;

use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::PathOptions;
use engine::metrics::distance_measures::{
    barycenter, center, closeness_vitality, global_efficiency, gutman_index, harmonic_diameter,
    hyper_wiener_index, local_efficiency, periphery, schultz_index, wiener_index, HyperWiener,
    Reachable,
};
use engine::paths::Tie;
use engine::{Error, WeightKind};

use crate::{guard, handle, slice_mut, thread_count, GraphHandle};

fn path_options(
    mode: i32,
    weighted: i32,
    weights_as: i32,
    threads: i32,
) -> Result<PathOptions, Error> {
    Ok(PathOptions {
        mode: DegreeMode::from_code(mode)?,
        weighted: weighted != 0,
        weights_as: WeightKind::from_code(weights_as)?,
        threads: thread_count(threads),
        ..Default::default()
    })
}

fn write_nodes(nodes: &[u32], out: *mut i64, count: *mut i64) {
    // SAFETY: the caller guarantees `out` holds n values and `count` one.
    unsafe {
        for (o, &v) in slice_mut(out, nodes.len()).iter_mut().zip(nodes) {
            *o = v as i64;
        }
        slice_mut(count, 1)[0] = nodes.len() as i64;
    }
}

/// The center (`which` 0) or periphery (`which` 1): node indices written to
/// `out`, their number to `*count` (0 when `reachable` is 0 and the network
/// is not strongly connected).
///
/// # Safety
/// `h` must be a live handle; `out` n values; `count` one value.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_center_periphery(
    h: *const GraphHandle,
    which: i32,
    reachable: i32,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    threads: i32,
    out: *mut i64,
    count: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = path_options(mode, weighted, weights_as, threads)?;
        let reach = Reachable::from_code(reachable)?;
        let nodes = match which {
            0 => center(g, opts, reach)?,
            1 => periphery(g, opts, reach)?,
            _ => return Err(Error::InvalidArgument(format!("unknown code {which}"))),
        };
        write_nodes(&nodes, out, count);
        Ok(())
    })
}

/// The barycenter: node indices written to `out` and their number to
/// `*count` (0 unless the network is strongly connected); when not null,
/// `totals` receives every node's total distance (Inf when it does not
/// reach every node).
///
/// # Safety
/// `h` must be a live handle; `out` n values; `count` one value; `totals`
/// null or n doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_barycenter(
    h: *const GraphHandle,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    tie: i32,
    threads: i32,
    out: *mut i64,
    count: *mut i64,
    totals: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let (nodes, sums) = barycenter(
            g,
            path_options(mode, weighted, weights_as, threads)?,
            Tie::from_code(tie)?,
        )?;
        write_nodes(&nodes, out, count);
        if !totals.is_null() {
            slice_mut(totals, g.node_count()).copy_from_slice(&sums);
        }
        Ok(())
    })
}

/// A distance index (see `index` above), one double.
///
/// # Safety
/// `h` must be a live handle; `out` one double.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_distance_index(
    h: *const GraphHandle,
    index: i32,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let o = path_options(mode, weighted, weights_as, threads)?;
        let value = match index {
            0 => wiener_index(g, o)?,
            1 => schultz_index(g, o)?,
            2 => gutman_index(g, o)?,
            3 => hyper_wiener_index(g, o, HyperWiener::Standard)?,
            4 => hyper_wiener_index(g, o, HyperWiener::NetworkX)?,
            5 => harmonic_diameter(g, o)?,
            6 => global_efficiency(g, o)?,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown distance index code {index}"
                )))
            }
        };
        slice_mut(out, 1)[0] = value;
        Ok(())
    })
}

/// Local efficiency (one double); `loop_is_neighbor != 0` counts a node
/// with a self-loop among its own neighbors, as NetworkX does.
///
/// # Safety
/// `h` must be a live handle; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_local_efficiency(
    h: *const GraphHandle,
    loop_is_neighbor: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = local_efficiency(g, loop_is_neighbor != 0, thread_count(threads))?;
        Ok(())
    })
}

/// Closeness vitality; `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `out` n doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_closeness_vitality(
    h: *const GraphHandle,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = closeness_vitality(g, path_options(mode, weighted, weights_as, threads)?)?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}
