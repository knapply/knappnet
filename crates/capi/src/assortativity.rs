//! C ABI for degree and attribute mixing as NetworkX computes it: mixing
//! counts, assortativity over node subsets, average neighbour degree, and
//! average degree connectivity.
//!
//! Integer codes used by these functions:
//! - degree kinds (`x`, `y`, `source`, `target`): 1 = out, 2 = in, 3 =
//!   in+out (ignored for undirected graphs).
//! - `coefficient`: 0 = degree, 1 = attribute (nominal), 2 = numeric.
//! - `nodes` (with `n_nodes`): a subset of nodes, or NULL for all.

use std::ffi::c_char;

use engine::metrics::assortativity::{
    attribute_assortativity, attribute_pairs, average_degree_connectivity, average_neighbor_degree,
    degree_assortativity, degree_pairs, mixing_counts, numeric_assortativity, DegreeKind,
    MixingCounts,
};
use engine::{Error, Graph};

use crate::similarity::node_indices;
use crate::{count, guard, handle, slice, slice_mut, GraphHandle};

/// # Safety
/// `nodes` must be null or hold `n_nodes` values.
unsafe fn subset(g: &Graph, nodes: *const i64, n_nodes: i64) -> Result<Option<Vec<u32>>, Error> {
    if nodes.is_null() {
        Ok(None)
    } else {
        node_indices(g, nodes, n_nodes, "nodes").map(Some)
    }
}

/// # Safety
/// `values` must hold `cap`, `counts` `cap * cap`, and `k` one value.
unsafe fn write_mixing(
    m: &MixingCounts,
    cap: i64,
    values: *mut f64,
    counts: *mut f64,
    k: *mut i64,
) -> Result<(), Error> {
    let cap = count(cap, "cap")?;
    let kk = m.values.len();
    slice_mut(k, 1)[0] = kk as i64;
    if kk <= cap {
        slice_mut(values, cap)[..kk].copy_from_slice(&m.values);
        // Row-major k by k into the first k * k slots.
        slice_mut(counts, cap * cap)[..kk * kk].copy_from_slice(&m.counts);
    }
    Ok(())
}

/// Degree mixing counts (NetworkX's `degree_mixing_dict` and
/// `degree_mixing_matrix` before normalization): the distinct degrees met
/// along arcs, increasing, into `values[0..*k)` and the count of
/// (x-degree, y-degree) pairs into `counts[i * k + j]` (row-major k by
/// k). When `*k > cap` only `*k` is written: call again with `cap >= *k`.
///
/// # Safety
/// `h` must be a live handle; `nodes` null or `n_nodes` values; `values`
/// `cap`, `counts` `cap * cap`, and `k` one value.
#[no_mangle]
pub unsafe extern "C" fn net_degree_mixing(
    h: *const GraphHandle,
    x: i32,
    y: i32,
    weighted: i32,
    nodes: *const i64,
    n_nodes: i64,
    cap: i64,
    values: *mut f64,
    counts: *mut f64,
    k: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let sub = subset(g, nodes, n_nodes)?;
        let p = degree_pairs(
            g,
            DegreeKind::from_code(x)?,
            DegreeKind::from_code(y)?,
            weighted != 0,
            sub.as_deref(),
        )?;
        write_mixing(&mixing_counts(&p)?, cap, values, counts, k)
    })
}

/// Attribute mixing counts (`attribute_mixing_dict` and
/// `attribute_mixing_matrix` before normalization) of node values
/// `labels` (n values), as for `net_degree_mixing`.
///
/// # Safety
/// As for `net_degree_mixing`; `labels` must hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_attribute_mixing(
    h: *const GraphHandle,
    labels: *const f64,
    nodes: *const i64,
    n_nodes: i64,
    cap: i64,
    values: *mut f64,
    counts: *mut f64,
    k: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let sub = subset(g, nodes, n_nodes)?;
        let p = attribute_pairs(g, slice(labels, g.node_count()), sub.as_deref())?;
        write_mixing(&mixing_counts(&p)?, cap, values, counts, k)
    })
}

/// Assortativity as NetworkX's coefficient functions over a node subset:
/// `coefficient` 0 = degree (`x`, `y`, `weighted`), 1 = attribute (nominal
/// `values`), 2 = numeric `values`. NaN when undefined.
///
/// # Safety
/// `h` must be a live handle; `values` null (degree) or n values; `nodes`
/// null or `n_nodes` values; `out` one value.
#[no_mangle]
pub unsafe extern "C" fn net_assortativity_subset(
    h: *const GraphHandle,
    coefficient: i32,
    x: i32,
    y: i32,
    weighted: i32,
    values: *const f64,
    nodes: *const i64,
    n_nodes: i64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let sub = subset(g, nodes, n_nodes)?;
        let need = || {
            if values.is_null() {
                Err(Error::InvalidArgument(
                    "this coefficient needs node values".into(),
                ))
            } else {
                Ok(slice(values, g.node_count()))
            }
        };
        let r = match coefficient {
            0 => degree_assortativity(
                g,
                DegreeKind::from_code(x)?,
                DegreeKind::from_code(y)?,
                weighted != 0,
                sub.as_deref(),
            )?,
            1 => attribute_assortativity(g, need()?, sub.as_deref())?,
            2 => numeric_assortativity(g, need()?, sub.as_deref())?,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown coefficient code {coefficient}"
                )))
            }
        };
        slice_mut(out, 1)[0] = r;
        Ok(())
    })
}

/// Average neighbour degree (NetworkX's `average_neighbor_degree`) of the
/// nodes of the subset (all nodes when NULL), in that order, into `out`.
///
/// # Safety
/// `h` must be a live handle; `nodes` null or `n_nodes` values; `out` must
/// hold `n_nodes` values (n when `nodes` is null).
#[no_mangle]
pub unsafe extern "C" fn net_average_neighbor_degree(
    h: *const GraphHandle,
    source: i32,
    target: i32,
    weighted: i32,
    nodes: *const i64,
    n_nodes: i64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let sub = subset(g, nodes, n_nodes)?;
        let r = average_neighbor_degree(
            g,
            DegreeKind::from_code(source)?,
            DegreeKind::from_code(target)?,
            weighted != 0,
            sub.as_deref(),
        )?;
        slice_mut(out, r.len()).copy_from_slice(&r);
        Ok(())
    })
}

/// Average degree connectivity (NetworkX's `average_degree_connectivity`):
/// the distinct source degrees of the subset's nodes, increasing, into
/// `degrees[0..*count)` and their k_nn into `values`.
///
/// # Safety
/// `h` must be a live handle; `nodes` null or `n_nodes` values; `degrees`
/// and `values` must hold n values and `count` one.
#[no_mangle]
pub unsafe extern "C" fn net_average_degree_connectivity(
    h: *const GraphHandle,
    source: i32,
    target: i32,
    weighted: i32,
    nodes: *const i64,
    n_nodes: i64,
    degrees: *mut f64,
    values: *mut f64,
    count_out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let sub = subset(g, nodes, n_nodes)?;
        let (k, v) = average_degree_connectivity(
            g,
            DegreeKind::from_code(source)?,
            DegreeKind::from_code(target)?,
            weighted != 0,
            sub.as_deref(),
        )?;
        let n = g.node_count();
        slice_mut(degrees, n)[..k.len()].copy_from_slice(&k);
        slice_mut(values, n)[..v.len()].copy_from_slice(&v);
        slice_mut(count_out, 1)[0] = k.len() as i64;
        Ok(())
    })
}
