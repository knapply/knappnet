//! C ABI for reach, edge betweenness, information centrality, bridging,
//! regular and automorphic equivalence, blockmodels, classical scaling,
//! and geospatial networks.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength.
//! - `crs`: 0 = longitude and latitude in degrees on a sphere of radius
//!   `radius`, 1 = planar.
//!
//! Matrices cross the interface in column-major order, as R stores them.
//! Proximity networks, whose size is not known in advance, come back in an
//! `EdgeListHandle`: the caller asks for its size, allocates three arrays,
//! copies the edges out, and frees the handle.

use std::ffi::c_char;

use engine::matrix::Dense;
use engine::metrics::blockmodel::{block_model, BlockModelOptions};
use engine::metrics::degree::DegreeMode;
use engine::metrics::reach::{
    bridging_coefficient, edge_betweenness, information_centrality, reach,
    EdgeBetweennessOptions, ReachOptions,
};
use engine::metrics::regular::{automorphic_equivalence, regular_equivalence, RegularOptions};
use engine::metrics::scaling::{classical_mds, stress};
use engine::spatial::{
    edge_lengths, group_geometry, nearest_neighbors, within_distance, Crs, Points,
    ProximityEdges,
};
use engine::{Error, WeightKind};

use crate::{guard, handle, slice, slice_mut, GraphHandle};

/// Opaque list of edges with distances.
pub struct EdgeListHandle(ProximityEdges);

fn count(n: i64, what: &str) -> Result<usize, Error> {
    usize::try_from(n).map_err(|_| Error::InvalidArgument(format!("negative {what}")))
}

/// Number (or with `normalized != 0`, proportion of the other nodes) of
/// nodes within distance `k` of every node: steps, or a length when
/// `weighted != 0`. `k` may be infinite.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_reach(
    h: *const GraphHandle,
    k: f64,
    mode: i32,
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
        let opts = ReachOptions {
            k,
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            normalized: normalized != 0,
            threads: threads.max(0) as usize,
        };
        let values = reach(g, opts)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Shortest-path betweenness of every edge, in input order.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `m` doubles (one per edge).
#[no_mangle]
pub unsafe extern "C" fn net_edge_betweenness(
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
        let opts = EdgeBetweennessOptions {
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            normalized: normalized != 0,
            threads: threads.max(0) as usize,
            ..Default::default()
        };
        let values = edge_betweenness(g, opts)?;
        slice_mut(out, g.edge_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Stephenson-Zelen information centrality. With `weighted != 0` edge
/// weights are tie strengths.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_information_centrality(
    h: *const GraphHandle,
    weighted: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let values = information_centrality(g, weighted != 0)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Hwang et al.'s bridging coefficient (NaN for isolated nodes).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_bridging_coefficient(
    h: *const GraphHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let values = bridging_coefficient(g);
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Maximal regular equivalence by CATREGE. `seed` may be null; otherwise
/// it holds `n_seed` (the node count) labels of a partition to refine.
/// `out` receives 0-based class numbers and `rounds` the number of
/// refinement rounds that split a class.
///
/// # Safety
/// `h` must be a live handle; `seed` must be null or hold `n_seed` values;
/// `out` must hold `n` 32-bit integers; `rounds` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_regular_equivalence(
    h: *const GraphHandle,
    seed: *const i64,
    n_seed: i64,
    weighted: i32,
    out: *mut i32,
    rounds: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let seed = if seed.is_null() {
            None
        } else {
            Some(slice(seed, count(n_seed, "seed length")?))
        };
        let r = regular_equivalence(
            g,
            seed,
            RegularOptions {
                weighted: weighted != 0,
            },
        )?;
        for (o, &c) in slice_mut(out, g.node_count()).iter_mut().zip(&r.membership) {
            *o = c as i32;
        }
        slice_mut(rounds, 1)[0] = r.rounds.min(i32::MAX as usize) as i32;
        Ok(())
    })
}

/// Automorphic equivalence classes (orbits), 0-based. Fails if the search
/// would refine more than `max_steps` colourings (`max_steps <= 0` uses
/// the engine's default).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` 32-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_automorphic_equivalence(
    h: *const GraphHandle,
    max_steps: i64,
    out: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let steps = if max_steps > 0 {
            max_steps as usize
        } else {
            engine::metrics::regular::DEFAULT_MAX_STEPS
        };
        let orbits = automorphic_equivalence(g, steps)?;
        for (o, &c) in slice_mut(out, g.node_count()).iter_mut().zip(&orbits) {
            *o = c as i32;
        }
        Ok(())
    })
}

/// Block densities of a partition. `membership` holds `n_membership` (the
/// node count) labels forming `k` distinct blocks; blocks are numbered in
/// order of first appearance, so passing labels 0..k-1 in that order keeps
/// them. `density`, `sums`, and `pairs` receive k-by-k matrices in
/// column-major order (rows are sending blocks); `summary` receives the
/// overall density and the R-squared fit.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold `n_membership`
/// values; `density`, `sums`, and `pairs` must each hold `k * k` doubles
/// and `summary` two.
#[no_mangle]
pub unsafe extern "C" fn net_block_model(
    h: *const GraphHandle,
    membership: *const i64,
    n_membership: i64,
    k: i64,
    weighted: i32,
    loops: i32,
    density: *mut f64,
    sums: *mut f64,
    pairs: *mut f64,
    summary: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let labels = slice(membership, count(n_membership, "membership length")?);
        let bm = block_model(
            g,
            labels,
            BlockModelOptions {
                weighted: weighted != 0,
                loops: loops != 0,
            },
        )?;
        let k = count(k, "block count")?;
        if bm.blocks() != k {
            return Err(Error::InvalidArgument(format!(
                "the partition has {} blocks, not {k}",
                bm.blocks()
            )));
        }
        slice_mut(density, k * k).copy_from_slice(&bm.density.to_col_major());
        slice_mut(sums, k * k).copy_from_slice(&bm.sums.to_col_major());
        slice_mut(pairs, k * k).copy_from_slice(&bm.pairs.to_col_major());
        let s = slice_mut(summary, 2);
        s[0] = bm.overall;
        s[1] = bm.r_squared;
        Ok(())
    })
}

/// Classical multidimensional scaling of an `n`-by-`n` dissimilarity
/// matrix `d` (column-major) into at most `k` dimensions. `points`
/// receives an n-by-k matrix (column-major) whose first `*dims` columns
/// hold coordinates (the rest NaN); `eigenvalues` all n eigenvalues in
/// decreasing order; `gof` the two goodness-of-fit values.
///
/// # Safety
/// `d` must hold `n * n` doubles, `points` `n * k`, `eigenvalues` `n`,
/// `gof` 2; `dims` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_classical_mds(
    d: *const f64,
    n: i64,
    k: i64,
    points: *mut f64,
    dims: *mut i32,
    eigenvalues: *mut f64,
    gof: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "point count")?;
        let k = count(k, "dimension count")?;
        let m = Dense::from_col_major(n, n, slice(d, n * n))?;
        let r = classical_mds(&m, k)?;
        let dst = slice_mut(points, n * k);
        dst.fill(f64::NAN);
        let cols = r.points.cols();
        for c in 0..cols {
            for i in 0..n {
                dst[c * n + i] = r.points.get(i, c);
            }
        }
        slice_mut(dims, 1)[0] = cols as i32;
        slice_mut(eigenvalues, n).copy_from_slice(&r.eigenvalues);
        slice_mut(gof, 2).copy_from_slice(&r.gof);
        Ok(())
    })
}

/// Kruskal's stress formula 1 of an `n`-by-`p` configuration `points`
/// (column-major) against an `n`-by-`n` dissimilarity matrix `d`.
///
/// # Safety
/// `d` must hold `n * n` doubles, `points` `n * p`; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_mds_stress(
    d: *const f64,
    n: i64,
    points: *const f64,
    p: i64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "point count")?;
        let p = count(p, "dimension count")?;
        let m = Dense::from_col_major(n, n, slice(d, n * n))?;
        let x = Dense::from_col_major(n, p, slice(points, n * p))?;
        slice_mut(out, 1)[0] = stress(&m, &x)?;
        Ok(())
    })
}

/// Length of every edge between its endpoints' coordinates `x`, `y` (one
/// per node, `n_points` of each).
///
/// # Safety
/// `h` must be a live handle; `x` and `y` must hold `n_points` doubles;
/// `out` must hold `m` doubles (one per edge).
#[no_mangle]
pub unsafe extern "C" fn net_edge_lengths(
    h: *const GraphHandle,
    x: *const f64,
    y: *const f64,
    n_points: i64,
    crs: i32,
    radius: f64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = count(n_points, "point count")?;
        let pts = Points::new(slice(x, n), slice(y, n), Crs::from_code(crs, radius)?)?;
        let values = edge_lengths(g, &pts)?;
        slice_mut(out, g.edge_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// A proximity network over `n` located points: with `k > 0`, the `k`
/// nearest neighbours of every point (directed edges); otherwise every
/// pair within distance `r` (undirected edges, smaller index first). On
/// success `*out` receives a handle to release with `net_edgelist_free`.
///
/// # Safety
/// `x` and `y` must hold `n` doubles; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_proximity(
    x: *const f64,
    y: *const f64,
    n: i64,
    crs: i32,
    radius: f64,
    r: f64,
    k: i64,
    threads: i32,
    out: *mut *mut EdgeListHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if out.is_null() {
            return Err(Error::InvalidArgument("output pointer is null".into()));
        }
        let n = count(n, "point count")?;
        let pts = Points::new(slice(x, n), slice(y, n), Crs::from_code(crs, radius)?)?;
        let threads = threads.max(0) as usize;
        let edges = if k > 0 {
            nearest_neighbors(&pts, k as usize, threads)?
        } else {
            within_distance(&pts, r, threads)?
        };
        *out = Box::into_raw(Box::new(EdgeListHandle(edges)));
        Ok(())
    })
}

/// Number of edges in an edge list, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live edge-list handle.
#[no_mangle]
pub unsafe extern "C" fn net_edgelist_size(h: *const EdgeListHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.from.len() as i64
}

/// Copies an edge list out: 0-based endpoints and distances, each array
/// holding `net_edgelist_size(h)` values.
///
/// # Safety
/// `h` must be a live edge-list handle; the arrays must have that size.
#[no_mangle]
pub unsafe extern "C" fn net_edgelist_copy(
    h: *const EdgeListHandle,
    from: *mut i64,
    to: *mut i64,
    distance: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if h.is_null() {
            return Err(Error::InvalidArgument("edge list handle is null".into()));
        }
        let e = &(*h).0;
        let m = e.from.len();
        for (o, &v) in slice_mut(from, m).iter_mut().zip(&e.from) {
            *o = v as i64;
        }
        for (o, &v) in slice_mut(to, m).iter_mut().zip(&e.to) {
            *o = v as i64;
        }
        slice_mut(distance, m).copy_from_slice(&e.distance);
        Ok(())
    })
}

/// Releases an edge-list handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live edge-list handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_edgelist_free(h: *mut EdgeListHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Centroid and dispersion of `k` groups of located points. `groups`
/// holds one label per point forming `k` distinct groups, numbered in
/// order of first appearance. Each output array holds `k` values: group
/// sizes, centroid x and y, mean distance from the centroid, and standard
/// distance.
///
/// # Safety
/// `x`, `y`, and `groups` must hold `n` values; each output must hold `k`
/// doubles.
#[no_mangle]
pub unsafe extern "C" fn net_group_geometry(
    x: *const f64,
    y: *const f64,
    groups: *const i64,
    n: i64,
    crs: i32,
    radius: f64,
    k: i64,
    sizes: *mut f64,
    cx: *mut f64,
    cy: *mut f64,
    mean_distance: *mut f64,
    standard_distance: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "point count")?;
        let k = count(k, "group count")?;
        let pts = Points::new(slice(x, n), slice(y, n), Crs::from_code(crs, radius)?)?;
        let gg = group_geometry(&pts, slice(groups, n))?;
        if gg.groups.len() != k {
            return Err(Error::InvalidArgument(format!(
                "the labels form {} groups, not {k}",
                gg.groups.len()
            )));
        }
        for (o, &s) in slice_mut(sizes, k).iter_mut().zip(&gg.sizes) {
            *o = s as f64;
        }
        slice_mut(cx, k).copy_from_slice(&gg.x);
        slice_mut(cy, k).copy_from_slice(&gg.y);
        slice_mut(mean_distance, k).copy_from_slice(&gg.mean_distance);
        slice_mut(standard_distance, k).copy_from_slice(&gg.standard_distance);
        Ok(())
    })
}
