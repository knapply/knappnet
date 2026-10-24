//! C ABI over the engine.
//!
//! This is the only surface the R package (and later other languages)
//! links against. It follows the pattern of linking a Rust static library
//! into an R package through a thin layer of C (Ooms, `ooms_hellorust`),
//! which keeps the R side free of any Rust-specific tooling.
//!
//! Conventions for every function:
//! - Node indices are 0-based 64-bit integers.
//! - Output buffers are allocated by the caller.
//! - The return value is 0 on success. On failure it is non-zero and a
//!   NUL-terminated message is written to `err` (truncated to `err_len`).
//! - Panics are caught and reported as errors; they never unwind into C.

use std::ffi::c_char;
use std::panic::{catch_unwind, AssertUnwindSafe};

use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::degree::{degree, strength, DegreeMode};
use engine::metrics::triad::triad_census;
use engine::{Error, Graph, WeightKind};

mod assortativity;
mod communities;
mod cores;
mod cut_measures;
mod cyber;
mod descriptive;
mod files;
mod frame;
mod groups;
mod isomorphism;
mod layout;
mod matching;
mod positions;
mod qap;
mod rewiring;
mod roles;
mod similarity;
mod spectral;
mod subgroups;
mod temporal;
mod topography;
mod twomode;
pub use files::*;

/// Opaque graph handle.
pub struct GraphHandle(Graph);

const OK: i32 = 0;
const FAILED: i32 = 1;

fn write_err(err: *mut c_char, err_len: usize, msg: &str) {
    if err.is_null() || err_len == 0 {
        return;
    }
    let bytes = msg.as_bytes();
    let n = bytes.len().min(err_len - 1);
    // SAFETY: the caller guarantees `err` points to `err_len` writable bytes.
    unsafe {
        std::ptr::copy_nonoverlapping(bytes.as_ptr(), err as *mut u8, n);
        *err.add(n) = 0;
    }
}

fn guard(err: *mut c_char, err_len: usize, f: impl FnOnce() -> Result<(), Error>) -> i32 {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(())) => OK,
        Ok(Err(e)) => {
            write_err(err, err_len, &e.to_string());
            FAILED
        }
        Err(_) => {
            write_err(err, err_len, "internal error (panic) in the network engine");
            FAILED
        }
    }
}

/// # Safety
/// `ptr` must be null or point to `len` readable elements.
unsafe fn slice<'a, T>(ptr: *const T, len: usize) -> &'a [T] {
    if len == 0 || ptr.is_null() {
        &[]
    } else {
        std::slice::from_raw_parts(ptr, len)
    }
}

/// # Safety
/// `ptr` must point to `len` writable elements.
unsafe fn slice_mut<'a, T>(ptr: *mut T, len: usize) -> &'a mut [T] {
    if len == 0 {
        &mut []
    } else {
        std::slice::from_raw_parts_mut(ptr, len)
    }
}

/// A non-negative count argument as `usize`, naming it on failure.
fn count(value: i64, what: &str) -> Result<usize, Error> {
    usize::try_from(value)
        .map_err(|_| Error::InvalidArgument(format!("{what} must be non-negative")))
}

/// A thread count argument: 0 (or a negative value) means every core.
fn thread_count(t: i32) -> usize {
    t.max(0) as usize
}

/// Dereferences a handle created by this crate, naming it in the error
/// when it is null.
///
/// # Safety
/// `h` must be null or a live handle created by this crate.
unsafe fn deref<'a, T>(h: *const T, what: &str) -> Result<&'a T, Error> {
    h.as_ref()
        .ok_or_else(|| Error::InvalidArgument(format!("{what} handle is null")))
}

fn handle<'a>(h: *const GraphHandle) -> Result<&'a Graph, Error> {
    // SAFETY: non-null handles are only created by `net_graph_new`.
    unsafe { deref(h, "graph") }.map(|g| &g.0)
}

/// Builds a graph. `weights` may be null. On success `*out` receives a
/// handle that must be released with `net_graph_free`.
///
/// # Safety
/// Pointers must reference `m` elements; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_graph_new(
    n: i64,
    src: *const i64,
    tgt: *const i64,
    weights: *const f64,
    m: i64,
    directed: i32,
    out: *mut *mut GraphHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if n < 0 || m < 0 {
            return Err(Error::InvalidArgument("negative node or edge count".into()));
        }
        let (n, m) = (n as usize, m as usize);
        let w = if weights.is_null() {
            None
        } else {
            Some(slice(weights, m))
        };
        let g = Graph::from_edges(n, slice(src, m), slice(tgt, m), w, directed != 0)?;
        *out = Box::into_raw(Box::new(GraphHandle(g)));
        Ok(())
    })
}

/// Releases a handle created by `net_graph_new`. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_graph_free(h: *mut GraphHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Node count of a graph, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_graph_node_count(h: *const GraphHandle) -> i64 {
    handle(h).map(|g| g.node_count() as i64).unwrap_or(-1)
}

/// Degree (`weighted == 0`) or strength (`weighted != 0`).
/// `mode`: 1 = out, 2 = in, 3 = all. `out` must hold n values.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_degree(
    h: *const GraphHandle,
    mode: i32,
    loops: i32,
    weighted: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let mode = DegreeMode::from_code(mode)?;
        let values = if weighted != 0 {
            strength(g, mode, loops != 0)?
        } else {
            degree(g, mode, loops != 0)
        };
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Shortest-path betweenness. `weights_as` is 0 (distance) or 1 (strength)
/// and only matters when `weighted` is non-zero. `threads == 0` uses all
/// cores.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_betweenness(
    h: *const GraphHandle,
    weighted: i32,
    weights_as: i32,
    normalized: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = BetweennessOptions {
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            normalized: normalized != 0,
            threads: thread_count(threads),
            ..Default::default()
        };
        let values = betweenness(g, opts)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Triad census in MAN order. `out` must hold 16 doubles.
///
/// # Safety
/// `h` must be a live handle; `out` must hold 16 doubles.
#[no_mangle]
pub unsafe extern "C" fn net_triad_census(
    h: *const GraphHandle,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let values = triad_census(g, thread_count(threads));
        slice_mut(out, 16).copy_from_slice(&values);
        Ok(())
    })
}
