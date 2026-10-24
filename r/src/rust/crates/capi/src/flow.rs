//! C ABI for maximum flows, minimum cuts, minimum-cost flows, and
//! Gomory-Hu trees, and the list handle shared by the connectivity, tree,
//! and cycle functions.
//!
//! Integer codes and conventions used by these functions:
//! - Terminals: `sources` and `sinks` hold 0-based node indices (at least
//!   one each, disjoint, no repeats).
//! - Capacities hold one value per edge in input order, non-negative and
//!   possibly `INFINITY`. For maximum flows and cuts a NULL `capacity`
//!   gives every edge capacity 1; for minimum-cost flows a NULL capacity
//!   is infinite and a NULL cost is 0 (NetworkX's defaults).
//! - `side` (minimum cuts): 0 = smallest source side (nodes reachable from
//!   the sources in the residual network), 1 = largest source side
//!   (NetworkX's `minimum_cut`).
//! - Output arrays marked "may be NULL" are skipped when NULL.
//!
//! A `ListsHandle` holds a sequence of integer lists, each with one double
//! attached (its meaning depends on the function): ask for the sizes,
//! allocate, copy, free.

use std::ffi::c_char;

use engine::flow::{
    gomory_hu_tree, max_flow_min_cost, maximum_flow, min_cost_flow, minimum_cut,
    minimum_node_cut_capacitated, CutSide,
};
use engine::Error;

use crate::{guard, handle, slice, slice_mut, GraphHandle};

/// Opaque sequence of integer lists, each carrying one double.
pub struct ListsHandle {
    pub(crate) lists: Vec<Vec<i64>>,
    pub(crate) values: Vec<f64>,
}

/// Hands a list collection to the caller.
///
/// # Safety
/// `out` must be null or writable.
pub(crate) unsafe fn give_lists(
    out: *mut *mut ListsHandle,
    lists: Vec<Vec<i64>>,
    values: Vec<f64>,
) -> Result<(), Error> {
    if out.is_null() {
        return Err(Error::InvalidArgument("output pointer is null".into()));
    }
    debug_assert_eq!(lists.len(), values.len());
    *out = Box::into_raw(Box::new(ListsHandle { lists, values }));
    Ok(())
}

/// Reads `len` node indices, checking each against `n`.
///
/// # Safety
/// `ptr` must be null or point to `len` readable values.
pub(crate) unsafe fn node_list(
    ptr: *const i64,
    len: i64,
    n: usize,
    what: &str,
) -> Result<Vec<usize>, Error> {
    let len = usize::try_from(len)
        .map_err(|_| Error::InvalidArgument(format!("the number of {what}s is negative")))?;
    if len > 0 && ptr.is_null() {
        return Err(Error::InvalidArgument(format!("{what} pointer is null")));
    }
    slice(ptr, len)
        .iter()
        .map(|&v| {
            if v < 0 || v as u64 >= n as u64 {
                Err(Error::InvalidArgument(format!(
                    "{what} {v} is not a node (n = {n})"
                )))
            } else {
                Ok(v as usize)
            }
        })
        .collect()
}

/// A single node index checked against `n`.
pub(crate) fn node(v: i64, n: usize, what: &str) -> Result<usize, Error> {
    if v < 0 || v as u64 >= n as u64 {
        Err(Error::InvalidArgument(format!(
            "{what} {v} is not a node (n = {n})"
        )))
    } else {
        Ok(v as usize)
    }
}

/// An optional per-item array of `len` doubles.
///
/// # Safety
/// `ptr` must be null or point to `len` readable doubles.
pub(crate) unsafe fn optional<'a>(ptr: *const f64, len: usize) -> Option<&'a [f64]> {
    if ptr.is_null() {
        None
    } else {
        Some(slice(ptr, len))
    }
}

/// Writes `values` to `out` unless it is null.
///
/// # Safety
/// `out` must be null or point to `values.len()` writable elements.
pub(crate) unsafe fn write_opt<T: Copy>(out: *mut T, values: &[T]) {
    if !out.is_null() {
        slice_mut(out, values.len()).copy_from_slice(values);
    }
}

/// Number of lists, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live lists handle.
#[no_mangle]
pub unsafe extern "C" fn net_lists_count(h: *const ListsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).lists.len() as i64
}

/// Total number of entries over all lists, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live lists handle.
#[no_mangle]
pub unsafe extern "C" fn net_lists_total(h: *const ListsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).lists.iter().map(|l| l.len() as i64).sum()
}

/// Copies the lists out: list `i` is `items[offsets[i] .. offsets[i + 1]]`
/// and carries `values[i]`. `offsets` must hold `count + 1` values, `items`
/// `total` values, and `values` (which may be null) `count` values.
///
/// # Safety
/// `h` must be a live lists handle; the arrays must have the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_lists_copy(
    h: *const ListsHandle,
    offsets: *mut i64,
    items: *mut i64,
    values: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if h.is_null() {
            return Err(Error::InvalidArgument("lists handle is null".into()));
        }
        let lists = &(*h).lists;
        let total: usize = lists.iter().map(|l| l.len()).sum();
        let offsets = slice_mut(offsets, lists.len() + 1);
        let items = slice_mut(items, total);
        offsets[0] = 0;
        let mut at = 0;
        for (i, l) in lists.iter().enumerate() {
            items[at..at + l.len()].copy_from_slice(l);
            at += l.len();
            offsets[i + 1] = at as i64;
        }
        write_opt(values, &(*h).values);
        Ok(())
    })
}

/// Releases a lists handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live lists handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_lists_free(h: *mut ListsHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Maximum flow from the sources to the sinks: the value into `value` and
/// the flow on every edge into `edge_flow` (m doubles, may be null; net
/// flow from first to second endpoint for undirected edges).
///
/// # Safety
/// `h` must be a live handle; `sources` and `sinks` must hold `n_sources`
/// and `n_sinks` values, `capacity` null or m values, `value` one double,
/// `edge_flow` null or m doubles.
#[no_mangle]
pub unsafe extern "C" fn net_max_flow(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
    sinks: *const i64,
    n_sinks: i64,
    capacity: *const f64,
    value: *mut f64,
    edge_flow: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let s = node_list(sources, n_sources, n, "source")?;
        let t = node_list(sinks, n_sinks, n, "sink")?;
        let f = maximum_flow(g, &s, &t, optional(capacity, g.edge_count()))?;
        slice_mut(value, 1)[0] = f.value;
        write_opt(edge_flow, &f.edge_flow);
        Ok(())
    })
}

/// Minimum cut between the sources and the sinks: its capacity into
/// `value`, 1 or 0 per node into `source_side` (n values, may be null), and
/// 1 or 0 per edge into `cut_edges` (m values, may be null; 1 for edges of
/// positive capacity crossing from the source side).
///
/// # Safety
/// As for `net_max_flow`; `source_side` null or n and `cut_edges` null or m
/// 32-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_min_cut(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
    sinks: *const i64,
    n_sinks: i64,
    capacity: *const f64,
    side: i32,
    value: *mut f64,
    source_side: *mut i32,
    cut_edges: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let s = node_list(sources, n_sources, n, "source")?;
        let t = node_list(sinks, n_sinks, n, "sink")?;
        let c = minimum_cut(
            g,
            &s,
            &t,
            optional(capacity, g.edge_count()),
            CutSide::from_code(side)?,
        )?;
        slice_mut(value, 1)[0] = c.value;
        let flags: Vec<i32> = c.source_side.iter().map(|&b| b as i32).collect();
        write_opt(source_side, &flags);
        let mut in_cut = vec![0i32; g.edge_count()];
        for &e in &c.edges {
            in_cut[e] = 1;
        }
        write_opt(cut_edges, &in_cut);
        Ok(())
    })
}

/// Minimum cut separating the sources from the sinks when nodes have
/// capacities (`node_capacity`, n values, NULL = 1) as well as edges
/// (`edge_capacity`, m values, NULL = infinite, so only nodes are cut);
/// terminals are never cut. Its capacity into `value` and 1 or 0 flags
/// into `cut_nodes` (n, may be null) and `cut_edges` (m, may be null).
///
/// # Safety
/// As for `net_min_cut`; `node_capacity` null or n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_min_node_cut_capacitated(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
    sinks: *const i64,
    n_sinks: i64,
    node_capacity: *const f64,
    edge_capacity: *const f64,
    side: i32,
    value: *mut f64,
    cut_nodes: *mut i32,
    cut_edges: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let s = node_list(sources, n_sources, n, "source")?;
        let t = node_list(sinks, n_sinks, n, "sink")?;
        let c = minimum_node_cut_capacitated(
            g,
            &s,
            &t,
            optional(node_capacity, n),
            optional(edge_capacity, g.edge_count()),
            CutSide::from_code(side)?,
        )?;
        slice_mut(value, 1)[0] = c.value;
        let mut nodes = vec![0i32; n];
        for &v in &c.nodes {
            nodes[v] = 1;
        }
        write_opt(cut_nodes, &nodes);
        let mut edges = vec![0i32; g.edge_count()];
        for &e in &c.edges {
            edges[e] = 1;
        }
        write_opt(cut_edges, &edges);
        Ok(())
    })
}

/// Minimum-cost flow of a directed network meeting `demand` (n doubles;
/// positive = receives, negative = supplies; summing to 0), with
/// `capacity` (m doubles, NULL = infinite) and `cost` (m doubles, NULL =
/// 0). The cost into `total_cost` and the flow on every edge into
/// `edge_flow` (m doubles, may be null).
///
/// # Safety
/// `h` must be a live handle; `demand` must hold n doubles, `capacity` and
/// `cost` null or m doubles, `total_cost` one double, `edge_flow` null or
/// m doubles.
#[no_mangle]
pub unsafe extern "C" fn net_min_cost_flow(
    h: *const GraphHandle,
    demand: *const f64,
    capacity: *const f64,
    cost: *const f64,
    total_cost: *mut f64,
    edge_flow: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let (n, m) = (g.node_count(), g.edge_count());
        if demand.is_null() {
            return Err(Error::InvalidArgument("demand pointer is null".into()));
        }
        let f = min_cost_flow(
            g,
            slice(demand, n),
            optional(capacity, m),
            optional(cost, m),
        )?;
        slice_mut(total_cost, 1)[0] = f.cost;
        write_opt(edge_flow, &f.edge_flow);
        Ok(())
    })
}

/// Maximum flow from `source` to `sink` of least cost (directed networks;
/// capacity NULL = infinite, cost NULL = 0): the flow value into `value`,
/// its cost into `total_cost`, and edge flows into `edge_flow` (m doubles,
/// may be null).
///
/// # Safety
/// As for `net_min_cost_flow`; `value` one double.
#[no_mangle]
pub unsafe extern "C" fn net_max_flow_min_cost(
    h: *const GraphHandle,
    source: i64,
    sink: i64,
    capacity: *const f64,
    cost: *const f64,
    value: *mut f64,
    total_cost: *mut f64,
    edge_flow: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let (n, m) = (g.node_count(), g.edge_count());
        let s = node(source, n, "source")?;
        let t = node(sink, n, "sink")?;
        let (f, v) = max_flow_min_cost(g, s, t, optional(capacity, m), optional(cost, m))?;
        slice_mut(value, 1)[0] = v;
        slice_mut(total_cost, 1)[0] = f.cost;
        write_opt(edge_flow, &f.edge_flow);
        Ok(())
    })
}

/// Gomory-Hu tree of an undirected network (capacity NULL = 1 per edge),
/// rooted at node 0: node v's tree neighbour toward the root into
/// `parent[v]` (-1 for the root) and the weight of that tree edge into
/// `weight[v]` (0 for the root). With nodes in order this is NetworkX's
/// `gomory_hu_tree`.
///
/// # Safety
/// `h` must be a live handle; `capacity` null or m doubles; `parent` n
/// 64-bit integers and `weight` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_gomory_hu_tree(
    h: *const GraphHandle,
    capacity: *const f64,
    parent: *mut i64,
    weight: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = gomory_hu_tree(g, optional(capacity, g.edge_count()))?;
        let n = g.node_count();
        let p: Vec<i64> = (0..n)
            .map(|v| if v == 0 { -1 } else { t.parent[v] as i64 })
            .collect();
        slice_mut(parent, n).copy_from_slice(&p);
        slice_mut(weight, n).copy_from_slice(&t.weight);
        Ok(())
    })
}
