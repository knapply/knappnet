//! C ABI for ego networks and two-mode (bipartite) networks.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `types`: one 32-bit integer per node, 0 for the first mode and 1 for
//!   the second.
//! - `measure` (two-mode centrality): 0 = degree, 1 = closeness,
//!   2 = betweenness.
//! - `unreachable` (two-mode closeness): 0 = ignore, 1 = Wasserman-Faust.
//! - `weight` (projection): 0 = shared partners, 1 = Newman's
//!   collaboration weights, 2 = sum of weight products.
//!
//! A projection comes back in a `ProjectionHandle`: ask for its edge
//! count, allocate three arrays, copy, and free the handle.

use std::ffi::c_char;

use engine::metrics::bipartite::{
    project, twomode_betweenness, twomode_closeness, twomode_degree, twomode_density,
    twomode_eigenvector, Projection, ProjectionOptions, ProjectionWeight, TwoModeClosenessOptions,
};
use engine::metrics::closeness::Unreachable;
use engine::metrics::degree::DegreeMode;
use engine::metrics::ego::{
    ego_summary, neighborhood_sizes, neighborhoods, EgoOptions, NeighborhoodOptions,
};
use engine::Error;

use crate::{guard, handle, slice, slice_mut, GraphHandle};

/// Opaque one-mode projection.
pub struct ProjectionHandle(Projection);

fn count(value: i64, what: &str) -> Result<usize, Error> {
    usize::try_from(value)
        .map_err(|_| Error::InvalidArgument(format!("{what} must not be negative")))
}

fn egos_from(ptr: *const i64, k: i64) -> Result<Option<Vec<usize>>, Error> {
    if ptr.is_null() {
        return Ok(None);
    }
    let k = count(k, "the number of egos")?;
    // SAFETY: the caller guarantees `ptr` points to `k` readable values.
    let raw = unsafe { slice(ptr, k) };
    raw.iter()
        .map(|&e| {
            usize::try_from(e)
                .map_err(|_| Error::InvalidArgument(format!("node index {e} is negative")))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

fn types_from(ptr: *const i32, n: usize) -> Result<Vec<bool>, Error> {
    if ptr.is_null() {
        return Err(Error::InvalidArgument("node types are missing".into()));
    }
    // SAFETY: the caller guarantees `ptr` points to `n` readable values.
    let raw = unsafe { slice(ptr, n) };
    raw.iter()
        .enumerate()
        .map(|(v, &t)| match t {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(Error::InvalidArgument(format!(
                "node {v} has type code {t}; use 0 for the first mode and 1 for the second"
            ))),
        })
        .collect()
}

fn neighborhood_options(
    order: i64,
    mode: i32,
    mindist: i64,
    threads: i32,
) -> Result<NeighborhoodOptions, Error> {
    Ok(NeighborhoodOptions {
        order: count(order, "the order")?,
        mode: DegreeMode::from_code(mode)?,
        mindist: count(mindist, "the minimum distance")?,
        threads: threads.max(0) as usize,
    })
}

/// Number of nodes within `order` steps of each ego and at least `mindist`
/// steps away (igraph's `ego_size()`). With `egos == NULL` every node is an
/// ego and `out` holds n values; otherwise `out` holds `k` values.
///
/// # Safety
/// `h` must be a live handle; `egos` must be null or hold `k` values;
/// `out` must hold n (or `k`) doubles.
#[no_mangle]
pub unsafe extern "C" fn net_ego_sizes(
    h: *const GraphHandle,
    egos: *const i64,
    k: i64,
    order: i64,
    mode: i32,
    mindist: i64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let list = egos_from(egos, k)?;
        let opts = neighborhood_options(order, mode, mindist, threads)?;
        let sizes = neighborhood_sizes(g, list.as_deref(), opts)?;
        for (o, s) in slice_mut(out, sizes.len()).iter_mut().zip(sizes) {
            *o = s as f64;
        }
        Ok(())
    })
}

/// Members of the neighborhood of each of `k` egos, in ascending node
/// order: those of ego i are `members[offsets[i] .. offsets[i + 1])`.
/// `offsets` must hold `k + 1` values and `members` `capacity` values;
/// get the total from `net_ego_sizes` with the same arguments. Fails
/// without writing members if `capacity` is too small.
///
/// # Safety
/// `h` must be a live handle; `egos` must hold `k` values, `offsets`
/// `k + 1`, and `members` `capacity`.
#[no_mangle]
pub unsafe extern "C" fn net_ego_members(
    h: *const GraphHandle,
    egos: *const i64,
    k: i64,
    order: i64,
    mode: i32,
    mindist: i64,
    threads: i32,
    offsets: *mut i64,
    members: *mut i64,
    capacity: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let list = egos_from(egos, k)?.unwrap_or_default();
        let opts = neighborhood_options(order, mode, mindist, threads)?;
        let groups = neighborhoods(g, &list, opts)?;
        let total: usize = groups.iter().map(|m| m.len()).sum();
        let capacity = count(capacity, "the capacity")?;
        if total > capacity {
            return Err(Error::InvalidArgument(format!(
                "the neighborhoods hold {total} nodes but room was given for {capacity}"
            )));
        }
        let off = slice_mut(offsets, list.len() + 1);
        let dst = slice_mut(members, total);
        let mut at = 0usize;
        off[0] = 0;
        for (i, group) in groups.iter().enumerate() {
            for &v in group {
                dst[at] = v as i64;
                at += 1;
            }
            off[i + 1] = at as i64;
        }
        Ok(())
    })
}

/// Egocentric measures of every node: number of alters, ties among them,
/// pairs of alters, density, weak components among alters without ego,
/// and (unless `ego_betweenness` is null) ego betweenness. Weights are
/// ignored; ties are distinct ties between distinct nodes.
///
/// # Safety
/// `h` must be a live handle; each output must hold n doubles
/// (`ego_betweenness` may be null).
#[no_mangle]
pub unsafe extern "C" fn net_ego_summary(
    h: *const GraphHandle,
    order: i64,
    mode: i32,
    threads: i32,
    size: *mut f64,
    ties: *mut f64,
    pairs: *mut f64,
    density: *mut f64,
    components: *mut f64,
    ego_betweenness: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let opts = EgoOptions {
            order: count(order, "the order")?,
            mode: DegreeMode::from_code(mode)?,
            betweenness: !ego_betweenness.is_null(),
            threads: threads.max(0) as usize,
        };
        let s = ego_summary(g, opts)?;
        slice_mut(size, n).copy_from_slice(&s.size);
        slice_mut(ties, n).copy_from_slice(&s.ties);
        slice_mut(pairs, n).copy_from_slice(&s.pairs);
        slice_mut(density, n).copy_from_slice(&s.density);
        slice_mut(components, n).copy_from_slice(&s.components);
        if let Some(b) = s.betweenness {
            slice_mut(ego_betweenness, n).copy_from_slice(&b);
        }
        Ok(())
    })
}

/// Two-mode density: distinct ties divided by n1 n2.
///
/// # Safety
/// `h` must be a live handle; `types` must hold n values; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_twomode_density(
    h: *const GraphHandle,
    types: *const i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = types_from(types, g.node_count())?;
        slice_mut(out, 1)[0] = twomode_density(g, &t)?;
        Ok(())
    })
}

/// Degree, closeness, or betweenness normalized for two-mode data
/// (Borgatti and Everett 1997) when `normalized != 0`. `unreachable` is
/// used by closeness only.
///
/// # Safety
/// `h` must be a live handle; `types` must hold n values; `out` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_twomode_centrality(
    h: *const GraphHandle,
    types: *const i32,
    measure: i32,
    normalized: i32,
    unreachable: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = types_from(types, g.node_count())?;
        let threads = threads.max(0) as usize;
        let values = match measure {
            0 => twomode_degree(g, &t, normalized != 0)?,
            1 => {
                let rule = match unreachable {
                    0 => Unreachable::Ignore,
                    1 => Unreachable::WassermanFaust,
                    _ => {
                        return Err(Error::InvalidArgument(format!(
                            "unknown code {unreachable} for unreachable nodes"
                        )))
                    }
                };
                let opts = TwoModeClosenessOptions {
                    normalized: normalized != 0,
                    unreachable: rule,
                    threads,
                };
                twomode_closeness(g, &t, opts)?
            }
            2 => twomode_betweenness(g, &t, normalized != 0, threads)?,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown two-mode centrality code {measure}"
                )))
            }
        };
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Two-mode eigenvector centrality: the leading singular vectors of the
/// binary incidence matrix, each mode scaled to Euclidean length 1.
/// `tol <= 0` and `max_iter <= 0` select the defaults (1e-12 and 10000).
/// `eigenvalue` receives the largest singular value and `connected` 1 if
/// the network is connected.
///
/// # Safety
/// `h` must be a live handle; `types` must hold n values; `out` n doubles;
/// `eigenvalue` and `connected` must each be writable.
#[no_mangle]
pub unsafe extern "C" fn net_twomode_eigenvector(
    h: *const GraphHandle,
    types: *const i32,
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
        let t = types_from(types, g.node_count())?;
        let tol = if tol > 0.0 { tol } else { 1e-12 };
        let max_iter = if max_iter > 0 {
            max_iter as usize
        } else {
            10_000
        };
        let e = twomode_eigenvector(g, &t, tol, max_iter, threads.max(0) as usize)?;
        slice_mut(out, g.node_count()).copy_from_slice(&e.values);
        slice_mut(eigenvalue, 1)[0] = e.eigenvalue;
        slice_mut(connected, 1)[0] = e.connected as i32;
        Ok(())
    })
}

/// One-mode projection onto the first (`second == 0`) or second mode. With
/// `loops != 0` the diagonal is kept as self-loops. On success `*out`
/// receives a handle that must be released with `net_projection_free`.
///
/// # Safety
/// `h` must be a live handle; `types` must hold n values; `out` must be
/// writable.
#[no_mangle]
pub unsafe extern "C" fn net_project(
    h: *const GraphHandle,
    types: *const i32,
    second: i32,
    weight: i32,
    loops: i32,
    threads: i32,
    out: *mut *mut ProjectionHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if out.is_null() {
            return Err(Error::InvalidArgument("output pointer is null".into()));
        }
        let t = types_from(types, g.node_count())?;
        let opts = ProjectionOptions {
            second: second != 0,
            weight: ProjectionWeight::from_code(weight)?,
            loops: loops != 0,
            threads: threads.max(0) as usize,
        };
        let p = project(g, &t, opts)?;
        *out = Box::into_raw(Box::new(ProjectionHandle(p)));
        Ok(())
    })
}

/// Number of edges in a projection, or -1 for a null handle.
///
/// # Safety
/// `p` must be null or a live projection handle.
#[no_mangle]
pub unsafe extern "C" fn net_projection_count(p: *const ProjectionHandle) -> i64 {
    if p.is_null() {
        return -1;
    }
    (*p).0.from.len() as i64
}

/// Copies a projection's edges: 0-based node indices of the original
/// network (`from <= to`) and weights, each `count` values.
///
/// # Safety
/// `p` must be a live projection handle; each output must hold
/// `net_projection_count(p)` values.
#[no_mangle]
pub unsafe extern "C" fn net_projection_copy(
    p: *const ProjectionHandle,
    from: *mut i64,
    to: *mut i64,
    weights: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if p.is_null() {
            return Err(Error::InvalidArgument("projection handle is null".into()));
        }
        let proj = &(*p).0;
        let m = proj.from.len();
        for (d, &s) in slice_mut(from, m).iter_mut().zip(&proj.from) {
            *d = s as i64;
        }
        for (d, &s) in slice_mut(to, m).iter_mut().zip(&proj.to) {
            *d = s as i64;
        }
        slice_mut(weights, m).copy_from_slice(&proj.weights);
        Ok(())
    })
}

/// Releases a projection handle. Null is ignored.
///
/// # Safety
/// `p` must be null or a live projection handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_projection_free(p: *mut ProjectionHandle) {
    if !p.is_null() {
        drop(Box::from_raw(p));
    }
}
