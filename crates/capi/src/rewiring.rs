//! C ABI for degree-preserving edge swaps, random and lattice reference
//! graphs, the small-world coefficients sigma and omega, and the
//! rich-club coefficient.
//!
//! Rewired graphs come back as edge lists in caller buffers that hold m
//! values each (a rewired simple graph has at most m edges): (u, v) with
//! u <= v for undirected graphs, sorted; `*count` receives their number.
//! Every function draws from the counter-based stream of `seed`.

use std::ffi::c_char;

use engine::matrix::Dense;
use engine::metrics::rewiring::{
    connected_double_edge_swap, directed_edge_swap, double_edge_swap, lattice_reference, omega,
    random_reference, rich_club_coefficient, sigma,
};
use engine::{Error, Graph};

use crate::{count, guard, handle, slice, slice_mut, thread_count, GraphHandle};

/// # Safety
/// `from` and `to` must hold m values and `n_edges` one.
unsafe fn write_edges(
    g: &Graph,
    e: &[(u32, u32)],
    from: *mut i64,
    to: *mut i64,
    n_edges: *mut i64,
) -> Result<(), Error> {
    let m = g.edge_count();
    if e.len() > m {
        return Err(Error::InvalidArgument(
            "more edges than buffer space".into(),
        ));
    }
    let (f, t) = (slice_mut(from, m), slice_mut(to, m));
    for (i, &(a, b)) in e.iter().enumerate() {
        f[i] = a as i64;
        t[i] = b as i64;
    }
    slice_mut(n_edges, 1)[0] = e.len() as i64;
    Ok(())
}

/// `nswap` double edge swaps of the undirected simple graph (NetworkX's
/// `double_edge_swap`); fails after `max_tries` attempts.
///
/// # Safety
/// `h` must be a live handle; `from` and `to` must hold m values and
/// `n_edges` one.
#[no_mangle]
pub unsafe extern "C" fn net_double_edge_swap(
    h: *const GraphHandle,
    nswap: i64,
    max_tries: i64,
    seed: u64,
    from: *mut i64,
    to: *mut i64,
    n_edges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let e = double_edge_swap(
            g,
            count(nswap, "nswap")? as u64,
            count(max_tries, "max_tries")? as u64,
            seed,
        )?;
        write_edges(g, &e, from, to, n_edges)
    })
}

/// Connectivity-preserving double edge swaps (NetworkX's
/// `connected_double_edge_swap`): `nswap` attempts; `*swaps` receives the
/// number kept.
///
/// # Safety
/// As for `net_double_edge_swap`; `swaps` must hold one value.
#[no_mangle]
pub unsafe extern "C" fn net_connected_double_edge_swap(
    h: *const GraphHandle,
    nswap: i64,
    seed: u64,
    from: *mut i64,
    to: *mut i64,
    n_edges: *mut i64,
    swaps: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let (e, k) = connected_double_edge_swap(g, count(nswap, "nswap")? as u64, seed)?;
        slice_mut(swaps, 1)[0] = k as i64;
        write_edges(g, &e, from, to, n_edges)
    })
}

/// `nswap` directed three-edge swaps (NetworkX's `directed_edge_swap`);
/// fails after `max_tries` attempts.
///
/// # Safety
/// As for `net_double_edge_swap`.
#[no_mangle]
pub unsafe extern "C" fn net_directed_edge_swap(
    h: *const GraphHandle,
    nswap: i64,
    max_tries: i64,
    seed: u64,
    from: *mut i64,
    to: *mut i64,
    n_edges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let e = directed_edge_swap(
            g,
            count(nswap, "nswap")? as u64,
            count(max_tries, "max_tries")? as u64,
            seed,
        )?;
        write_edges(g, &e, from, to, n_edges)
    })
}

/// A random reference graph (NetworkX's `random_reference`).
///
/// # Safety
/// As for `net_double_edge_swap`.
#[no_mangle]
pub unsafe extern "C" fn net_random_reference(
    h: *const GraphHandle,
    niter: i64,
    connectivity: i32,
    seed: u64,
    from: *mut i64,
    to: *mut i64,
    n_edges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let e = random_reference(g, count(niter, "niter")? as u64, connectivity != 0, seed)?;
        write_edges(g, &e, from, to, n_edges)
    })
}

/// A lattice reference graph (NetworkX's `lattice_reference`) towards the
/// n-by-n distance matrix `dist` (column-major; NULL for the ring
/// distance between node positions).
///
/// # Safety
/// As for `net_double_edge_swap`; `dist` must be null or hold n * n values.
#[no_mangle]
pub unsafe extern "C" fn net_lattice_reference(
    h: *const GraphHandle,
    niter: i64,
    dist: *const f64,
    connectivity: i32,
    seed: u64,
    from: *mut i64,
    to: *mut i64,
    n_edges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let d = if dist.is_null() {
            None
        } else {
            Some(Dense::from_col_major(n, n, slice(dist, n * n))?)
        };
        let e = lattice_reference(
            g,
            count(niter, "niter")? as u64,
            d.as_ref(),
            connectivity != 0,
            seed,
        )?;
        write_edges(g, &e, from, to, n_edges)
    })
}

/// The small-world coefficient sigma (NetworkX's `sigma`) into `*out`;
/// reference r uses stream (seed, r).
///
/// # Safety
/// `h` must be a live handle; `out` must hold one value.
#[no_mangle]
pub unsafe extern "C" fn net_small_world_sigma(
    h: *const GraphHandle,
    niter: i64,
    nrand: i64,
    seed: u64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = sigma(
            g,
            count(niter, "niter")? as u64,
            count(nrand, "nrand")?,
            seed,
            thread_count(threads),
        )?;
        Ok(())
    })
}

/// The small-world coefficient omega (NetworkX's `omega`) into `*out`;
/// round r uses streams (seed, 2r) and (seed, 2r + 1).
///
/// # Safety
/// `h` must be a live handle; `out` must hold one value.
#[no_mangle]
pub unsafe extern "C" fn net_small_world_omega(
    h: *const GraphHandle,
    niter: i64,
    nrand: i64,
    seed: u64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = omega(
            g,
            count(niter, "niter")? as u64,
            count(nrand, "nrand")?,
            seed,
            thread_count(threads),
        )?;
        Ok(())
    })
}

/// The rich-club coefficient for k = 0, 1, ... (NetworkX's
/// `rich_club_coefficient`) into `out[0..*count)`; normalized != 0 divides
/// by the mean over n_samples graphs rewired with q * m double edge swaps
/// (sample r from stream (seed, r)).
///
/// # Safety
/// `h` must be a live handle; `out` must hold n values and `count_out`
/// one.
#[no_mangle]
pub unsafe extern "C" fn net_rich_club(
    h: *const GraphHandle,
    normalized: i32,
    q: i64,
    n_samples: i64,
    seed: u64,
    threads: i32,
    out: *mut f64,
    count_out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let rc = rich_club_coefficient(
            g,
            normalized != 0,
            count(q, "q")? as u64,
            count(n_samples, "n_samples")?,
            seed,
            thread_count(threads),
        )?;
        slice_mut(out, g.node_count())[..rc.len()].copy_from_slice(&rc);
        slice_mut(count_out, 1)[0] = rc.len() as i64;
        Ok(())
    })
}
