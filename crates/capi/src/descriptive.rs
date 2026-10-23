//! C ABI for reciprocity, transitivity and clustering, distances and the
//! measures built on them, reachability, cut points and bridges, and the
//! counts a network summary reports.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength.
//! - `reciprocity` mode: 0 = default, 1 = ratio, 2 = dyadic,
//!   3 = correlation.
//! - `weighting`: 0 = binary, 1 = Barrat, 2 = Onnela.
//! - `undefined`: 0 = NaN (igraph), 1 = zero (NetworkX).
//! - `sources`: 0-based node indices; a null pointer or a negative count
//!   means every node, in order.

use std::ffi::c_char;

use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::{
    average_clustering, cut_structure, distance_rows, eccentricity, edge_counts, local_clustering,
    path_summary, reachability, reciprocity, transitivity, ClusteringOptions, PathOptions,
    ReciprocityMode, Undefined, Weighting,
};
use engine::{Error, Graph, WeightKind};

use crate::{guard, handle, slice, slice_mut, thread_count, GraphHandle};

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

fn clustering_options(
    weighting: i32,
    undefined: i32,
    threads: i32,
) -> Result<ClusteringOptions, Error> {
    Ok(ClusteringOptions {
        weighting: Weighting::from_code(weighting)?,
        undefined: Undefined::from_code(undefined)?,
        threads: thread_count(threads),
    })
}

/// # Safety
/// `sources` must be null or hold `n_sources` values.
unsafe fn source_list(
    g: &Graph,
    sources: *const i64,
    n_sources: i64,
) -> Result<Option<Vec<u32>>, Error> {
    if sources.is_null() || n_sources < 0 {
        return Ok(None);
    }
    let n = g.node_count();
    let mut out = Vec::with_capacity(n_sources as usize);
    for &v in slice(sources, n_sources as usize) {
        if v < 0 || v as u64 >= n as u64 {
            return Err(Error::InvalidArgument(format!(
                "source node index {v} is out of range for a network with {n} nodes"
            )));
        }
        out.push(v as u32);
    }
    Ok(Some(out))
}

/// Edge count of a graph (parallel edges and loops included), or -1 for a
/// null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_graph_edge_count(h: *const GraphHandle) -> i64 {
    handle(h).map(|g| g.edge_count() as i64).unwrap_or(-1)
}

/// Reciprocity of the network (one double). See the codes above.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_reciprocity(
    h: *const GraphHandle,
    mode: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = reciprocity(g, ReciprocityMode::from_code(mode)?);
        Ok(())
    })
}

/// Global transitivity (one double). `undefined` sets the value when the
/// network has no connected triples.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one double.
#[no_mangle]
pub unsafe extern "C" fn net_transitivity(
    h: *const GraphHandle,
    undefined: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = clustering_options(0, undefined, threads)?;
        slice_mut(out, 1)[0] = transitivity(g, opts);
        Ok(())
    })
}

/// Local clustering coefficients (`average == 0`, n doubles) or their mean
/// (`average != 0`, one double).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles, or one when
/// `average != 0`.
#[no_mangle]
pub unsafe extern "C" fn net_clustering(
    h: *const GraphHandle,
    average: i32,
    weighting: i32,
    undefined: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = clustering_options(weighting, undefined, threads)?;
        if average != 0 {
            slice_mut(out, 1)[0] = average_clustering(g, opts)?;
        } else {
            let values = local_clustering(g, opts)?;
            slice_mut(out, g.node_count()).copy_from_slice(&values);
        }
        Ok(())
    })
}

/// Shortest-path distances from each source to every node, written to `out`
/// as a column-major matrix with one row per source and one column per
/// node. Unreachable nodes are at distance infinity.
///
/// # Safety
/// `h` must be a live handle; `sources` must be null or hold `n_sources`
/// values; `out` must hold (number of sources) * n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_distances(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
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
        let n = g.node_count();
        let opts = path_options(mode, weighted, weights_as, threads)?;
        let src = source_list(g, sources, n_sources)?;
        let rows = src.as_ref().map_or(n, |s| s.len());
        let dst = slice_mut(out, rows * n);
        distance_rows(g, src.as_deref(), opts, |k, row| {
            for (j, &d) in row.iter().enumerate() {
                dst[j * rows + k] = d;
            }
        })?;
        Ok(())
    })
}

/// Totals over the shortest-path searches from the sources. `out` receives
/// seven doubles: the number of sources, the ordered pairs (source, other
/// node) joined by a path, the sum of their distances, the longest finite
/// distance, the number of sources that do not reach every other node, the
/// mean distance, and the diameter. With `unconnected == 0` the mean and
/// diameter are infinite when some source does not reach every other node.
///
/// # Safety
/// `h` must be a live handle; `sources` must be null or hold `n_sources`
/// values; `out` must hold seven doubles.
#[no_mangle]
pub unsafe extern "C" fn net_path_summary(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    unconnected: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = path_options(mode, weighted, weights_as, threads)?;
        let src = source_list(g, sources, n_sources)?;
        let s = path_summary(g, src.as_deref(), opts)?;
        let unconnected = unconnected != 0;
        slice_mut(out, 7).copy_from_slice(&[
            s.sources as f64,
            s.pairs as f64,
            s.total,
            s.longest,
            s.incomplete as f64,
            s.mean(unconnected),
            s.diameter(unconnected),
        ]);
        Ok(())
    })
}

/// Eccentricity of every node: its largest distance to the nodes it
/// reaches (0 if none).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_eccentricity(
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
        let values = eccentricity(g, path_options(mode, weighted, weights_as, threads)?)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Number of other nodes each node reaches (`mode` 1), is reached from
/// (2), or is joined to ignoring direction (3).
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_reachability(
    h: *const GraphHandle,
    mode: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let values = reachability(g, DegreeMode::from_code(mode)?, thread_count(threads))?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Cut points and bridges of the undirected reading of the network:
/// `cut_point[v]` is 1 if node v is a cut point, and `bridge[e]` is 1 if
/// edge e (in input order) is a bridge, else 0.
///
/// # Safety
/// `h` must be a live handle; `cut_point` must hold `n` and `bridge` `m`
/// 32-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_cut_structure(
    h: *const GraphHandle,
    cut_point: *mut i32,
    bridge: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = cut_structure(g);
        for (d, &x) in slice_mut(cut_point, g.node_count())
            .iter_mut()
            .zip(&c.cut_point)
        {
            *d = x as i32;
        }
        for (d, &x) in slice_mut(bridge, g.edge_count()).iter_mut().zip(&c.bridge) {
            *d = x as i32;
        }
        Ok(())
    })
}

/// Counts for a summary: self-loops, edges repeating an earlier edge
/// between the same nodes, and nodes without an edge to another node.
///
/// # Safety
/// `h` must be a live handle; `out` must hold three doubles.
#[no_mangle]
pub unsafe extern "C" fn net_edge_counts(
    h: *const GraphHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = edge_counts(g);
        slice_mut(out, 3).copy_from_slice(&[c.loops as f64, c.multiple as f64, c.isolates as f64]);
        Ok(())
    })
}
