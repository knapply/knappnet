//! C ABI for cohesive subgroups: cores, cliques, and communities.
//!
//! Results that are collections of node groups (such as maximal cliques)
//! come back in a `GroupsHandle`: the caller asks for its size, allocates
//! two arrays, copies the groups out as flat members plus offsets (the CSR
//! layout used throughout the engine), and frees the handle. The groups are
//! computed once, however many times they are read.

use std::ffi::c_char;

use engine::metrics::cliques::{clique_number, max_cliques};
use engine::metrics::community::{louvain, modularity};
use engine::metrics::cores::coreness;
use engine::metrics::degree::DegreeMode;
use engine::Error;

use crate::{guard, handle, slice, slice_mut, GraphHandle};

/// Opaque collection of node groups.
pub struct GroupsHandle(Vec<Vec<u32>>);

/// Coreness (k-core number) of every node. `mode`: 1 = out, 2 = in,
/// 3 = all; ignored for undirected graphs.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_coreness(
    h: *const GraphHandle,
    mode: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let mode = DegreeMode::from_code(mode)?;
        let values = coreness(g, mode);
        for (o, v) in slice_mut(out, g.node_count()).iter_mut().zip(values) {
            *o = v as f64;
        }
        Ok(())
    })
}

/// Maximal cliques with at least `min_size` members. On success `*out`
/// receives a handle that must be released with `net_groups_free`.
/// `threads == 0` uses all available cores.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_max_cliques(
    h: *const GraphHandle,
    min_size: i64,
    threads: i32,
    out: *mut *mut GroupsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if out.is_null() {
            return Err(Error::InvalidArgument("output pointer is null".into()));
        }
        let groups = max_cliques(g, min_size.max(0) as usize, threads.max(0) as usize);
        *out = Box::into_raw(Box::new(GroupsHandle(groups)));
        Ok(())
    })
}

/// Number of groups, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live groups handle.
#[no_mangle]
pub unsafe extern "C" fn net_groups_count(h: *const GroupsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.len() as i64
}

/// Total number of members over all groups, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live groups handle.
#[no_mangle]
pub unsafe extern "C" fn net_groups_total(h: *const GroupsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.iter().map(|g| g.len() as i64).sum()
}

/// Copies the groups out. Members of group `i` are
/// `members[offsets[i] .. offsets[i + 1]]`, as 0-based node indices in
/// increasing order. `offsets` must hold `count + 1` values and `members`
/// `total` values.
///
/// # Safety
/// `h` must be a live groups handle; the arrays must have the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_groups_copy(
    h: *const GroupsHandle,
    offsets: *mut i64,
    members: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if h.is_null() {
            return Err(Error::InvalidArgument("groups handle is null".into()));
        }
        let groups = &(*h).0;
        let total: usize = groups.iter().map(|g| g.len()).sum();
        let offsets = slice_mut(offsets, groups.len() + 1);
        let members = slice_mut(members, total);
        let mut at = 0;
        offsets[0] = 0;
        for (i, group) in groups.iter().enumerate() {
            for &v in group {
                members[at] = v as i64;
                at += 1;
            }
            offsets[i + 1] = at as i64;
        }
        Ok(())
    })
}

/// Releases a handle created by a function returning groups. Null is
/// ignored.
///
/// # Safety
/// `h` must be null or a live groups handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_groups_free(h: *mut GroupsHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Size of the largest clique, written to `*out`.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_clique_number(
    h: *const GraphHandle,
    threads: i32,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = clique_number(g, threads.max(0) as usize) as i64;
        Ok(())
    })
}

/// Modularity of a partition. `membership` holds `n` non-negative group
/// labels. NaN (not an error) when the network has no edge weight.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold `n` values; `out`
/// must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_modularity(
    h: *const GraphHandle,
    membership: *const i64,
    n: i64,
    resolution: f64,
    weighted: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if n < 0 {
            return Err(Error::InvalidArgument("negative membership length".into()));
        }
        let q = modularity(g, slice(membership, n as usize), resolution, weighted != 0)?;
        slice_mut(out, 1)[0] = q;
        Ok(())
    })
}

/// Louvain community detection on an undirected graph. Writes 0-based
/// group labels (numbered by first appearance) to `membership` and the
/// partition's modularity to `*q`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold `n` values; `q` must
/// be writable.
#[no_mangle]
pub unsafe extern "C" fn net_louvain(
    h: *const GraphHandle,
    resolution: f64,
    weighted: i32,
    seed: u64,
    membership: *mut i64,
    q: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let r = louvain(g, resolution, weighted != 0, seed)?;
        for (o, c) in slice_mut(membership, g.node_count())
            .iter_mut()
            .zip(r.membership)
        {
            *o = c as i64;
        }
        slice_mut(q, 1)[0] = r.modularity;
        Ok(())
    })
}
