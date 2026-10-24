//! C ABI for node and edge connectivity, minimum cuts, disjoint paths,
//! k-components, k-edge-components, local bridges, and the component
//! structure of directed networks.
//!
//! Integer codes used by these functions:
//! - `kind`: 0 = node (vertex) connectivity, 1 = edge connectivity.
//! - `convention` (node cuts): 0 = engine (the cut nearest the source;
//!   adjacent nodes are an error), 1 = NetworkX (`minimum_st_node_cut`;
//!   empty for adjacent nodes).
//! - `order` (component numbering): 0 = sorted neighbours, 1 = input
//!   order (NetworkX's numbering for a graph built from the same edges).
//! - `s`, `t`: node indices; for the cut functions s = t = -1 asks for a
//!   global cut.

use std::ffi::c_char;

use engine::metrics::connectivity::{
    all_pairs_node_connectivity, attracting_components, average_node_connectivity, condensation,
    edge_connectivity, edge_disjoint_paths, is_semiconnected, k_components, k_edge_components,
    local_bridges, local_edge_connectivity, local_node_connectivity, minimum_edge_cut,
    minimum_node_cut, minimum_st_edge_cut, minimum_st_node_cut, node_connectivity,
    node_disjoint_paths, CutConvention,
};
use engine::order::NeighborOrder;
use engine::Error;

use crate::flow::{give_lists, node, ListsHandle};
use crate::{count, guard, handle, slice_mut, thread_count, GraphHandle};

fn kind_code(kind: i32) -> Result<bool, Error> {
    match kind {
        0 => Ok(true),
        1 => Ok(false),
        _ => Err(Error::InvalidArgument(format!(
            "unknown connectivity kind code {kind}"
        ))),
    }
}

fn labels(values: &[u32], out: &mut [i32]) -> Result<(), Error> {
    for (d, &v) in out.iter_mut().zip(values) {
        *d = i32::try_from(v).map_err(|_| Error::TooLarge {
            what: "component count",
            value: v as usize,
        })?;
    }
    Ok(())
}

/// Local node (`kind == 0`) or edge (`kind == 1`) connectivity between
/// `s` and `t` into `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one 64-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_local_connectivity(
    h: *const GraphHandle,
    s: i64,
    t: i64,
    kind: i32,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let (s, t) = (node(s, n, "source")?, node(t, n, "target")?);
        let v = if kind_code(kind)? {
            local_node_connectivity(g, s, t)?
        } else {
            local_edge_connectivity(g, s, t)?
        };
        slice_mut(out, 1)[0] = v as i64;
        Ok(())
    })
}

/// Global node or edge connectivity into `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one 64-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_connectivity(
    h: *const GraphHandle,
    kind: i32,
    threads: i32,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = thread_count(threads);
        let v = if kind_code(kind)? {
            node_connectivity(g, t)?
        } else {
            edge_connectivity(g, t)?
        };
        slice_mut(out, 1)[0] = v as i64;
        Ok(())
    })
}

/// Local node connectivity of every ordered pair, as an n-by-n
/// column-major matrix with zeros on the diagonal, into `out`; the mean
/// over pairs into `average` (may be null).
///
/// # Safety
/// `h` must be a live handle; `out` must hold n * n doubles and `average`
/// be null or hold one.
#[no_mangle]
pub unsafe extern "C" fn net_all_pairs_node_connectivity(
    h: *const GraphHandle,
    threads: i32,
    out: *mut f64,
    average: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let m = all_pairs_node_connectivity(g, thread_count(threads));
        slice_mut(out, n * n).copy_from_slice(&m.to_col_major());
        if !average.is_null() {
            slice_mut(average, 1)[0] = average_node_connectivity(&m, g.is_directed());
        }
        Ok(())
    })
}

/// A minimum edge cut between `s` and `t`, or a global one when both are
/// -1: its pairs (u, v) into `from` and `to` (at most m each) and their
/// number into `size`.
///
/// # Safety
/// `h` must be a live handle; `from` and `to` must hold m 64-bit integers
/// and `size` one.
#[no_mangle]
pub unsafe extern "C" fn net_min_edge_cut(
    h: *const GraphHandle,
    s: i64,
    t: i64,
    threads: i32,
    from: *mut i64,
    to: *mut i64,
    size: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let cut = if s == -1 && t == -1 {
            minimum_edge_cut(g, thread_count(threads))?
        } else {
            minimum_st_edge_cut(g, node(s, n, "source")?, node(t, n, "target")?)?
        };
        let m = g.edge_count();
        let (f, tt) = (slice_mut(from, m), slice_mut(to, m));
        for (i, &(u, v)) in cut.iter().enumerate() {
            f[i] = u as i64;
            tt[i] = v as i64;
        }
        slice_mut(size, 1)[0] = cut.len() as i64;
        Ok(())
    })
}

/// A minimum node cut between `s` and `t`, or a global one when both are
/// -1: its nodes into `nodes` (at most n) and their number into `size`.
///
/// # Safety
/// `h` must be a live handle; `nodes` must hold n 64-bit integers and
/// `size` one.
#[no_mangle]
pub unsafe extern "C" fn net_min_node_cut(
    h: *const GraphHandle,
    s: i64,
    t: i64,
    convention: i32,
    threads: i32,
    nodes: *mut i64,
    size: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let c = CutConvention::from_code(convention)?;
        let cut = if s == -1 && t == -1 {
            minimum_node_cut(g, c, thread_count(threads))?
        } else {
            minimum_st_node_cut(g, node(s, n, "source")?, node(t, n, "target")?, c)?
        };
        let out = slice_mut(nodes, n);
        for (i, &v) in cut.iter().enumerate() {
            out[i] = v as i64;
        }
        slice_mut(size, 1)[0] = cut.len() as i64;
        Ok(())
    })
}

/// Node-disjoint (`kind == 0`) or edge-disjoint (`kind == 1`) paths from
/// `s` to `t`, as many as possible, each a node sequence (value 0).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_disjoint_paths(
    h: *const GraphHandle,
    s: i64,
    t: i64,
    kind: i32,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let (s, t) = (node(s, n, "source")?, node(t, n, "target")?);
        let paths = if kind_code(kind)? {
            node_disjoint_paths(g, s, t)?
        } else {
            edge_disjoint_paths(g, s, t)?
        };
        let values = vec![0.0; paths.len()];
        let lists = paths
            .into_iter()
            .map(|p| p.into_iter().map(|v| v as i64).collect())
            .collect();
        give_lists(out, lists, values)
    })
}

/// k-components of an undirected network: one list per component (sorted
/// nodes) carrying its k, for k = 1, 2, ... in turn.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_k_components(
    h: *const GraphHandle,
    threads: i32,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let levels = k_components(g, thread_count(threads))?;
        let mut lists = Vec::new();
        let mut values = Vec::new();
        for (i, level) in levels.into_iter().enumerate() {
            for comp in level {
                lists.push(comp.into_iter().map(|v| v as i64).collect());
                values.push((i + 1) as f64);
            }
        }
        give_lists(out, lists, values)
    })
}

/// k-edge-connected components (k >= 1): a class label per node (0-based,
/// in order of each class's lowest node) into `membership`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n 32-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_k_edge_components(
    h: *const GraphHandle,
    k: i64,
    threads: i32,
    membership: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let k = count(k, "k")? as u64;
        let m = k_edge_components(g, k, thread_count(threads))?;
        labels(&m, slice_mut(membership, g.node_count()))
    })
}

/// Local bridges of an undirected network: 1/0 per edge into `flags` and
/// the span of each local bridge (INFINITY without another path, NaN for
/// other edges) into `span` (may be null).
///
/// # Safety
/// `h` must be a live handle; `flags` must hold m 32-bit integers and
/// `span` be null or hold m doubles.
#[no_mangle]
pub unsafe extern "C" fn net_local_bridges(
    h: *const GraphHandle,
    flags: *mut i32,
    span: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let b = local_bridges(g)?;
        let f: Vec<i32> = b.is_local_bridge.iter().map(|&x| x as i32).collect();
        slice_mut(flags, g.edge_count()).copy_from_slice(&f);
        crate::flow::write_opt(span, &b.span);
        Ok(())
    })
}

/// Strongly connected components of a directed network numbered as the
/// nodes of its condensation (the order Tarjan's algorithm completes
/// them, a reverse topological order) into `membership`, and their
/// number into `count`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n 32-bit integers and
/// `count` one 64-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_condensation(
    h: *const GraphHandle,
    order: i32,
    membership: *mut i32,
    count: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let c = condensation(g, NeighborOrder::from_code(order)?)?;
        labels(&c, slice_mut(membership, g.node_count()))?;
        slice_mut(count, 1)[0] = c.iter().map(|&x| x as i64 + 1).max().unwrap_or(0);
        Ok(())
    })
}

/// Attracting components of a directed network: each node's component
/// label (numbered as in `net_condensation`) if its component is
/// attracting, else -1, into `membership`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n 32-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_attracting_components(
    h: *const GraphHandle,
    order: i32,
    membership: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let a = attracting_components(g, NeighborOrder::from_code(order)?)?;
        let out = slice_mut(membership, g.node_count());
        for (d, x) in out.iter_mut().zip(a) {
            *d = match x {
                Some(c) => c as i32,
                None => -1,
            };
        }
        Ok(())
    })
}

/// Whether a directed network is semiconnected (1) or not (0), into `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one 32-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_is_semiconnected(
    h: *const GraphHandle,
    out: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = is_semiconnected(g)? as i32;
        Ok(())
    })
}
