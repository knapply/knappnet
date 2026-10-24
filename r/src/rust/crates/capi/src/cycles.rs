//! C ABI for cycle bases, minimum cycle bases, simple and chordless
//! cycles, girth, and finding a cycle.
//!
//! Integer codes used by these functions:
//! - `order`: 0 = sorted neighbours, 1 = input order (NetworkX's results
//!   for a graph built from the same edges in the same order).
//! - `length_bound`: the longest cycle (in nodes) to report; -1 for none.
//! - `limit`: the most cycles to report; 0 for all.
//! - `orientation` (find_cycle, directed networks): 0 = original,
//!   1 = reverse, 2 = ignore.

use std::ffi::c_char;

use engine::metrics::cycles::{
    chordless_cycles, cycle_basis, find_cycle, girth, minimum_cycle_basis, simple_cycles,
    Orientation,
};
use engine::order::NeighborOrder;
use engine::Error;

use crate::flow::{give_lists, node, ListsHandle};
use crate::{count, guard, handle, slice_mut, thread_count, GraphHandle};

fn bound(length_bound: i64) -> Result<Option<usize>, Error> {
    if length_bound == -1 {
        Ok(None)
    } else {
        count(length_bound, "the length bound").map(Some)
    }
}

fn limit_of(limit: i64) -> Result<Option<usize>, Error> {
    let l = count(limit, "the limit")?;
    Ok((l > 0).then_some(l))
}

/// Hands node-sequence cycles to the caller with the given values.
///
/// # Safety
/// `out` must be writable.
unsafe fn give_cycles(
    out: *mut *mut ListsHandle,
    cycles: Vec<Vec<u32>>,
    values: Vec<f64>,
) -> Result<(), Error> {
    let lists = cycles
        .into_iter()
        .map(|c| c.into_iter().map(|v| v as i64).collect())
        .collect();
    give_lists(out, lists, values)
}

/// A cycle basis of an undirected network (NetworkX's `cycle_basis` for
/// the same neighbour order; `root` >= 0 starts the first search there,
/// -1 as NetworkX), one node sequence per cycle (value 0).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_cycle_basis(
    h: *const GraphHandle,
    order: i32,
    root: i64,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let root = if root < 0 {
            None
        } else {
            Some(node(root, g.node_count(), "root")?)
        };
        let c = cycle_basis(g, NeighborOrder::from_code(order)?, root)?;
        let values = vec![0.0; c.len()];
        give_cycles(out, c, values)
    })
}

/// A minimum cycle basis of an undirected network (de Pina's algorithm),
/// one node sequence per cycle carrying its weight.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_minimum_cycle_basis(
    h: *const GraphHandle,
    weighted: i32,
    threads: i32,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let b = minimum_cycle_basis(g, weighted != 0, thread_count(threads))?;
        give_cycles(out, b.cycles, b.weights)
    })
}

/// Every simple cycle (Johnson's algorithm), one node sequence each
/// (value 0).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_simple_cycles(
    h: *const GraphHandle,
    length_bound: i64,
    limit: i64,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = simple_cycles(g, bound(length_bound)?, limit_of(limit)?)?;
        let values = vec![0.0; c.len()];
        give_cycles(out, c, values)
    })
}

/// Every chordless cycle, one node sequence each (value 0).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_chordless_cycles(
    h: *const GraphHandle,
    length_bound: i64,
    limit: i64,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = chordless_cycles(g, bound(length_bound)?, limit_of(limit)?)?;
        let values = vec![0.0; c.len()];
        give_cycles(out, c, values)
    })
}

/// Length of a shortest cycle into `out` (INFINITY without one).
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_girth(
    h: *const GraphHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = girth(g);
        Ok(())
    })
}

/// A cycle found by depth-first search over edges (NetworkX's
/// `find_cycle`) from `source` (-1 for every node in turn): its edges'
/// positions into `edges` and 1 for an edge traversed against its
/// direction into `reversed` (both m values, allocated by the caller), and
/// the number of edges into `size` (0 without a cycle).
///
/// # Safety
/// `h` must be a live handle; `edges` must hold m 64-bit and `reversed` m
/// 32-bit integers, `size` one 64-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_find_cycle(
    h: *const GraphHandle,
    source: i64,
    orientation: i32,
    order: i32,
    edges: *mut i64,
    reversed: *mut i32,
    size: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let source = if source < 0 {
            None
        } else {
            Some(node(source, g.node_count(), "source")?)
        };
        let c = find_cycle(
            g,
            source,
            Orientation::from_code(orientation)?,
            NeighborOrder::from_code(order)?,
        )?
        .unwrap_or_default();
        let m = g.edge_count();
        let (e, r) = (slice_mut(edges, m), slice_mut(reversed, m));
        for (i, st) in c.iter().enumerate() {
            e[i] = st.edge as i64;
            r[i] = st.reversed as i32;
        }
        slice_mut(size, 1)[0] = c.len() as i64;
        Ok(())
    })
}
