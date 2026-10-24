//! C ABI for k-shells, k-crusts, k-coronas, k-trusses, and onion layers.
//!
//! Integer codes used by these functions:
//! - `which` (core subsets): 0 = k-core, 1 = k-shell, 2 = k-crust, 3 =
//!   k-corona.

use std::ffi::c_char;

use engine::metrics::cores::{core_subset, k_truss, onion_layers, truss_numbers, CoreSubset};

use crate::{guard, handle, slice_mut, GraphHandle};

/// The nodes of the k-core, k-shell, k-crust, or k-corona (increasing)
/// into `nodes` and their number into `*count`. `has_k == 0` takes
/// NetworkX's default k (the largest core number; one less for the crust;
/// the corona needs k).
///
/// # Safety
/// `h` must be a live handle; `nodes` must hold n values and `count` one.
#[no_mangle]
pub unsafe extern "C" fn net_core_subset(
    h: *const GraphHandle,
    which: i32,
    k: i64,
    has_k: i32,
    nodes: *mut i64,
    count: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let set = core_subset(g, CoreSubset::from_code(which)?, (has_k != 0).then_some(k))?;
        for (d, &v) in slice_mut(nodes, g.node_count()).iter_mut().zip(&set) {
            *d = v as i64;
        }
        slice_mut(count, 1)[0] = set.len() as i64;
        Ok(())
    })
}

/// Truss number of every edge (0 for self-loops) into `out` (m values).
///
/// # Safety
/// `h` must be a live handle; `out` must hold m values.
#[no_mangle]
pub unsafe extern "C" fn net_truss_numbers(
    h: *const GraphHandle,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = truss_numbers(g);
        for (d, &x) in slice_mut(out, g.edge_count()).iter_mut().zip(&t) {
            *d = x as i64;
        }
        Ok(())
    })
}

/// The k-truss: indices of its edges (increasing) into `edges` with their
/// number in `*n_edges`, and the nodes they touch into `nodes` with their
/// number in `*n_nodes`.
///
/// # Safety
/// `h` must be a live handle; `edges` must hold m values, `nodes` n
/// values, and the counts one each.
#[no_mangle]
pub unsafe extern "C" fn net_k_truss(
    h: *const GraphHandle,
    k: i64,
    edges: *mut i64,
    n_edges: *mut i64,
    nodes: *mut i64,
    n_nodes: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let (e, v) = k_truss(g, k);
        for (d, &x) in slice_mut(edges, g.edge_count()).iter_mut().zip(&e) {
            *d = x as i64;
        }
        for (d, &x) in slice_mut(nodes, g.node_count()).iter_mut().zip(&v) {
            *d = x as i64;
        }
        slice_mut(n_edges, 1)[0] = e.len() as i64;
        slice_mut(n_nodes, 1)[0] = v.len() as i64;
        Ok(())
    })
}

/// Onion layer of every node (from 1) into `out` (n values).
///
/// # Safety
/// `h` must be a live handle; `out` must hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_onion_layers(
    h: *const GraphHandle,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        for (d, &x) in slice_mut(out, g.node_count())
            .iter_mut()
            .zip(&onion_layers(g))
        {
            *d = x as i64;
        }
        Ok(())
    })
}
