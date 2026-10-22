//! C ABI for density, components, fragmentation, centralization,
//! closeness, harmonic centrality, and eigenvector centrality.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength.
//! - `kind` (components): 0 = weak, 1 = strong.
//! - `unreachable` (closeness): 0 = ignore (igraph), 1 = Wasserman-Faust
//!   (NetworkX), 2 = zero (sna), 3 = largest length plus 1 (UCINET).
//! - `measure` (centralization): 0 = degree, 1 = betweenness, 2 = closeness.
//! - `convention`: 0 = igraph, 1 = NetworkX, 2 = sna, 3 = UCINET.
//! - `scale` (eigenvector): 0 = largest score 1, 1 = unit Euclidean norm.
//! - `loops` (eigenvector): 0 = twice, 1 = once, 2 = ignore.

use std::ffi::c_char;

use engine::conventions::Convention;
use engine::metrics::closeness::{closeness, harmonic_centrality, ClosenessOptions, Unreachable};
use engine::metrics::degree::DegreeMode;
use engine::metrics::eigenvector::{
    eigenvector_centrality, EigenScale, EigenvectorOptions, LoopWeight,
};
use engine::metrics::topography::{
    centralization, components, density, fragmentation, simplify, CentralizationMeasure,
    CentralizationOptions, ComponentKind, ParallelWeight,
};
use engine::{Error, WeightKind};

use crate::{guard, handle, slice_mut, GraphHandle};

/// Builds the adjacency-matrix reading of a graph, as sna and UCINET see
/// it: parallel edges collapse into one, keeping the largest (`keep == 0`)
/// or smallest (`keep == 1`) weight, and self-loops are dropped. With
/// `undirected != 0` a directed graph is also symmetrized (a tie wherever
/// either arc exists). On success `*out` receives a new handle that must be
/// released with `net_graph_free`.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_graph_simplify(
    h: *const GraphHandle,
    undirected: i32,
    keep: i32,
    out: *mut *mut GraphHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let s = simplify(g, undirected != 0, ParallelWeight::from_code(keep)?)?;
        *out = Box::into_raw(Box::new(GraphHandle(s)));
        Ok(())
    })
}

/// Edge density. `loops != 0` counts self-loops among the possible edges.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_density(
    h: *const GraphHandle,
    loops: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = density(g, loops != 0);
        Ok(())
    })
}

/// Component membership (0-based labels in order of each component's
/// lowest-numbered node) into `membership`, and the number of components
/// into `count`. Component sizes follow from the membership.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold `n` 32-bit integers
/// and `count` one 64-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_components(
    h: *const GraphHandle,
    kind: i32,
    membership: *mut i32,
    count: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = components(g, ComponentKind::from_code(kind)?);
        let dst = slice_mut(membership, g.node_count());
        for (d, &m) in dst.iter_mut().zip(&c.membership) {
            *d = i32::try_from(m).map_err(|_| Error::TooLarge {
                what: "component count",
                value: m as usize,
            })?;
        }
        slice_mut(count, 1)[0] = c.count() as i64;
        Ok(())
    })
}

/// Fragmentation: the share of ordered pairs of distinct nodes in different
/// (weak or strong) components.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_fragmentation(
    h: *const GraphHandle,
    kind: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = fragmentation(g, ComponentKind::from_code(kind)?);
        Ok(())
    })
}

/// Freeman centralization of degree, betweenness, or closeness.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_centralization(
    h: *const GraphHandle,
    measure: i32,
    mode: i32,
    loops: i32,
    normalized: i32,
    convention: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = CentralizationOptions {
            measure: CentralizationMeasure::from_code(measure)?,
            mode: DegreeMode::from_code(mode)?,
            loops: loops != 0,
            normalized: normalized != 0,
            convention: Convention::from_code(convention)?,
            threads: threads.max(0) as usize,
        };
        slice_mut(out, 1)[0] = centralization(g, opts)?;
        Ok(())
    })
}

/// Closeness (`harmonic == 0`) or harmonic centrality (`harmonic != 0`).
/// `unreachable` only matters for closeness.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_closeness(
    h: *const GraphHandle,
    harmonic: i32,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    normalized: i32,
    unreachable: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = ClosenessOptions {
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            normalized: normalized != 0,
            unreachable: Unreachable::from_code(unreachable)?,
            threads: threads.max(0) as usize,
            ..Default::default()
        };
        let values = if harmonic != 0 {
            harmonic_centrality(g, opts)?
        } else {
            closeness(g, opts)?
        };
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Eigenvector centrality by power iteration. `tol <= 0` and
/// `max_iter <= 0` select the defaults (1e-12 and 10000). On success `eigenvalue` receives the principal eigenvalue and
/// `connected` 1 if the network is (strongly) connected, else 0.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles; `eigenvalue`
/// and `connected` must each be writable.
#[no_mangle]
pub unsafe extern "C" fn net_eigenvector(
    h: *const GraphHandle,
    mode: i32,
    weighted: i32,
    scale: i32,
    loops: i32,
    tol: f64,
    max_iter: i64,
    threads: i32,
    out: *mut f64,
    eigenvalue: *mut f64,
    connected: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let defaults = EigenvectorOptions::default();
        let opts = EigenvectorOptions {
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            scale: EigenScale::from_code(scale)?,
            loops: LoopWeight::from_code(loops)?,
            tol: if tol > 0.0 { tol } else { defaults.tol },
            max_iter: if max_iter > 0 {
                max_iter as usize
            } else {
                defaults.max_iter
            },
            threads: threads.max(0) as usize,
        };
        let e = eigenvector_centrality(g, opts)?;
        slice_mut(out, g.node_count()).copy_from_slice(&e.values);
        slice_mut(eigenvalue, 1)[0] = e.eigenvalue;
        slice_mut(connected, 1)[0] = e.connected as i32;
        Ok(())
    })
}
