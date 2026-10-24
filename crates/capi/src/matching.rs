//! C ABI for matchings, greedy colouring, random maximal independent sets,
//! square clustering, and generalized degree.
//!
//! Integer codes used by these functions:
//! - `strategy` (colouring): 0 = largest first, 1 = random sequential,
//!   2 = smallest last, 3 = independent set, 4 = connected sequential
//!   breadth-first, 5 = depth-first, 6 = saturation largest first (DSATUR).

use std::ffi::c_char;

use engine::metrics::coloring::{greedy_color, maximal_independent_set, ColoringStrategy};
use engine::metrics::matching::{
    check_matching, max_weight_matching, maximal_matching, min_weight_matching,
};
use engine::metrics::squares::{generalized_degree, square_clustering};
use engine::Graph;

use crate::similarity::node_indices;
use crate::{count, guard, handle, slice_mut, thread_count, GraphHandle};

/// # Safety
/// `edges` must hold m values and `found` one.
unsafe fn write_edges(g: &Graph, e: &[usize], edges: *mut i64, found: *mut i64) {
    for (d, &x) in slice_mut(edges, g.edge_count()).iter_mut().zip(e) {
        *d = x as i64;
    }
    slice_mut(found, 1)[0] = e.len() as i64;
}

/// A maximum-weight matching (Edmonds' blossom algorithm; with
/// `max_cardinality != 0`, of largest weight among the largest matchings):
/// matched edge indices, increasing, into `edges[0..*found)`.
///
/// # Safety
/// `h` must be a live handle; `edges` must hold m values and `found` one.
#[no_mangle]
pub unsafe extern "C" fn net_max_weight_matching(
    h: *const GraphHandle,
    weighted: i32,
    max_cardinality: i32,
    edges: *mut i64,
    found: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let e = max_weight_matching(g, weighted != 0, max_cardinality != 0)?;
        write_edges(g, &e, edges, found);
        Ok(())
    })
}

/// A matching of least weight among the largest matchings (NetworkX's
/// `min_weight_matching`), as for `net_max_weight_matching`.
///
/// # Safety
/// As for `net_max_weight_matching`.
#[no_mangle]
pub unsafe extern "C" fn net_min_weight_matching(
    h: *const GraphHandle,
    weighted: i32,
    edges: *mut i64,
    found: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let e = min_weight_matching(g, weighted != 0)?;
        write_edges(g, &e, edges, found);
        Ok(())
    })
}

/// A maximal matching by a greedy scan of the edges in order (NetworkX's
/// `maximal_matching` over `G.edges()`), as for `net_max_weight_matching`.
///
/// # Safety
/// As for `net_max_weight_matching`.
#[no_mangle]
pub unsafe extern "C" fn net_maximal_matching(
    h: *const GraphHandle,
    edges: *mut i64,
    found: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        write_edges(g, &maximal_matching(g), edges, found);
        Ok(())
    })
}

/// Tests the `k` pairs `(us[i], vs[i])`: `out` receives 1 or 0 for "is a
/// matching", "is a maximal matching", "is a perfect matching" (NetworkX's
/// `is_matching`, `is_maximal_matching`, `is_perfect_matching`).
///
/// # Safety
/// `h` must be a live handle; `us` and `vs` must hold `k` values and `out`
/// three.
#[no_mangle]
pub unsafe extern "C" fn net_check_matching(
    h: *const GraphHandle,
    us: *const i64,
    vs: *const i64,
    k: i64,
    out: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let a = node_indices(g, us, k, "us")?;
        let b = node_indices(g, vs, k, "vs")?;
        let c = check_matching(g, &a, &b)?;
        slice_mut(out, 3).copy_from_slice(&[c.matching as i32, c.maximal as i32, c.perfect as i32]);
        Ok(())
    })
}

/// A greedy colouring (colours from 0) into `out` (n values). `seed`
/// drives the random sequential strategy.
///
/// # Safety
/// `h` must be a live handle; `out` must hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_greedy_color(
    h: *const GraphHandle,
    strategy: i32,
    seed: u64,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = greedy_color(g, ColoringStrategy::from_code(strategy)?, seed)?;
        for (d, &x) in slice_mut(out, g.node_count()).iter_mut().zip(&c) {
            *d = x as i64;
        }
        Ok(())
    })
}

/// A random maximal independent set containing `start` (`n_start` nodes;
/// one random node when 0), increasing, into `nodes[0..*found)`. Draws
/// from stream (seed, 0).
///
/// # Safety
/// `h` must be a live handle; `start` must hold `n_start` values, `nodes`
/// n values, and `found` one.
#[no_mangle]
pub unsafe extern "C" fn net_maximal_independent_set(
    h: *const GraphHandle,
    start: *const i64,
    n_start: i64,
    seed: u64,
    nodes: *mut i64,
    found: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let s = node_indices(g, start, n_start, "start")?;
        let set = maximal_independent_set(g, &s, seed)?;
        for (d, &x) in slice_mut(nodes, g.node_count()).iter_mut().zip(&set) {
            *d = x as i64;
        }
        slice_mut(found, 1)[0] = set.len() as i64;
        Ok(())
    })
}

/// Square clustering of every node (NetworkX's `square_clustering`) into
/// `out` (n values).
///
/// # Safety
/// `h` must be a live handle; `out` must hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_square_clustering(
    h: *const GraphHandle,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let s = square_clustering(g, thread_count(threads))?;
        slice_mut(out, g.node_count()).copy_from_slice(&s);
        Ok(())
    })
}

/// Generalized degree (NetworkX's `generalized_degree`): node v's entries
/// are positions offsets[v] .. offsets[v + 1] of `ks` (a triangle count k,
/// increasing) and `counts` (its number of neighbours sharing k triangles
/// with v). `*total` receives the number of entries (at most 2m); when it
/// exceeds `cap` only `*total` is written.
///
/// # Safety
/// `h` must be a live handle; `offsets` must hold n + 1 values, `ks` and
/// `counts` `cap` values, and `total` one.
#[no_mangle]
pub unsafe extern "C" fn net_generalized_degree(
    h: *const GraphHandle,
    threads: i32,
    offsets: *mut i64,
    ks: *mut i64,
    counts: *mut i64,
    cap: i64,
    total: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let cap = count(cap, "cap")?;
        let gd = generalized_degree(g, thread_count(threads));
        let size: usize = gd.iter().map(Vec::len).sum();
        slice_mut(total, 1)[0] = size as i64;
        if size > cap {
            return Ok(());
        }
        let off = slice_mut(offsets, g.node_count() + 1);
        let (k_out, c_out) = (slice_mut(ks, cap), slice_mut(counts, cap));
        let mut pos = 0;
        for (v, row) in gd.iter().enumerate() {
            off[v] = pos as i64;
            for &(k, c) in row {
                k_out[pos] = k as i64;
                c_out[pos] = c as i64;
                pos += 1;
            }
        }
        off[g.node_count()] = pos as i64;
        Ok(())
    })
}
