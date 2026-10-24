//! C ABI for shortest paths: path graphs, single-source, multi-source, and
//! all-pairs searches with parents, predecessor lists, bidirectional and A*
//! search, negative lengths and negative cycles, Johnson and
//! Floyd-Warshall, every shortest path of a pair, and k shortest simple
//! paths.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength (lengths 1/w).
//! - `tie`: 0 = relative tolerance 1e-10, 1 = exact equality (NetworkX).
//! - `method` (negative lengths): 0 = Bellman-Ford, 1 = Goldberg-Radzik.
//! - Node indices are 0-based; -1 marks "no node" in parent arrays and "no
//!   target" where a target is optional.
//! - Matrices are column-major with one row per source.

use std::ffi::c_char;

use engine::metrics::degree::DegreeMode;
use engine::paths::astar::astar;
use engine::paths::bidirectional::bidirectional;
use engine::paths::dense::floyd_warshall;
use engine::paths::enumerate::{all_shortest_paths, count_shortest_paths};
use engine::paths::negative::{
    bellman_ford, find_negative_cycle, goldberg_radzik, has_negative_cycle, johnson,
    tree_predecessors, Outcome,
};
use engine::paths::search::{all_pairs, SearchOptions, Searcher};
use engine::paths::simple::{k_shortest_simple_paths, shortest_path_avoiding};
use engine::paths::{ArcLengths, NodeLists, PathGraph, PathGraphOptions, Tie, NO_NODE};
use engine::{Error, WeightKind};

use crate::{deref, guard, handle, slice, slice_mut, thread_count, GraphHandle};

/// Opaque handle to a path graph.
pub struct PathGraphHandle(PathGraph);

/// Opaque handle to a collection of node lists (paths, predecessor lists,
/// cycles), each with an optional value.
pub struct PathListsHandle(NodeLists);

fn path_graph<'a>(h: *const PathGraphHandle) -> Result<&'a PathGraph, Error> {
    // SAFETY: non-null handles are only created by `net_path_graph_new`.
    unsafe { deref(h, "path graph") }.map(|p| &p.0)
}

fn lists<'a>(h: *const PathListsHandle) -> Result<&'a NodeLists, Error> {
    // SAFETY: non-null handles are only created by this module.
    unsafe { deref(h, "path list") }.map(|p| &p.0)
}

/// Stores a collection in `*out`.
///
/// # Safety
/// `out` must be writable.
unsafe fn emit(out: *mut *mut PathListsHandle, l: NodeLists) -> Result<(), Error> {
    if out.is_null() {
        return Err(Error::InvalidArgument(
            "the output handle pointer is null".into(),
        ));
    }
    *out = Box::into_raw(Box::new(PathListsHandle(l)));
    Ok(())
}

fn node(pg: &PathGraph, v: i64, what: &str) -> Result<u32, Error> {
    if v < 0 || v as u64 >= pg.node_count() as u64 {
        return Err(Error::InvalidArgument(format!(
            "{what} node index {v} is out of range for a network with {} nodes",
            pg.node_count()
        )));
    }
    Ok(v as u32)
}

/// # Safety
/// `nodes` must be null or hold `k` values.
unsafe fn node_list(
    pg: &PathGraph,
    nodes: *const i64,
    k: i64,
    what: &str,
) -> Result<Option<Vec<u32>>, Error> {
    if nodes.is_null() || k < 0 {
        return Ok(None);
    }
    slice(nodes, k as usize)
        .iter()
        .map(|&v| node(pg, v, what))
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

fn index(v: u32) -> i64 {
    if v == NO_NODE {
        -1
    } else {
        v as i64
    }
}

/// Writes a row-major matrix (`rows` by `n`) column-major.
fn transpose_into<T: Copy, U>(src: &[T], rows: usize, n: usize, out: &mut [U], f: impl Fn(T) -> U) {
    for k in 0..rows {
        for v in 0..n {
            out[v * rows + k] = f(src[k * n + v]);
        }
    }
}

/// Builds a path graph: the arcs of `h` oriented for `mode`, with unit
/// lengths (`weighted == 0`), the edge weights read per `weights_as`, or,
/// when `lengths` is not null, one length per edge (`n_lengths` must be
/// the edge count; +Inf hides an edge). `collapse_parallel != 0` keeps only
/// the shortest of parallel arcs. Release with `net_path_graph_free`.
///
/// # Safety
/// `h` must be a live handle; `lengths` null or `n_lengths` values; `out`
/// writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_path_graph_new(
    h: *const GraphHandle,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    lengths: *const f64,
    n_lengths: i64,
    collapse_parallel: i32,
    out: *mut *mut PathGraphHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let per_edge;
        let lengths = if !lengths.is_null() {
            if n_lengths < 0 || n_lengths as usize != g.edge_count() {
                return Err(Error::LengthMismatch {
                    what: "edge lengths",
                    expected: g.edge_count(),
                    found: n_lengths.max(0) as usize,
                });
            }
            per_edge = slice(lengths, n_lengths as usize);
            ArcLengths::PerEdge(per_edge)
        } else if weighted != 0 {
            ArcLengths::Weights(WeightKind::from_code(weights_as)?)
        } else {
            ArcLengths::Unit
        };
        let pg = PathGraph::new(
            g,
            PathGraphOptions {
                mode: DegreeMode::from_code(mode)?,
                lengths,
                collapse_parallel: collapse_parallel != 0,
            },
        )?;
        if out.is_null() {
            return Err(Error::InvalidArgument(
                "the output handle pointer is null".into(),
            ));
        }
        *out = Box::into_raw(Box::new(PathGraphHandle(pg)));
        Ok(())
    })
}

/// Releases a path graph. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_path_graph_free(h: *mut PathGraphHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Number of arcs of a path graph, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_path_graph_arc_count(h: *const PathGraphHandle) -> i64 {
    path_graph(h).map(|p| p.arc_count() as i64).unwrap_or(-1)
}

/// Number of lists, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_path_lists_count(h: *const PathListsHandle) -> i64 {
    lists(h).map(|l| l.len() as i64).unwrap_or(-1)
}

/// Total number of nodes over all lists, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_path_lists_total(h: *const PathListsHandle) -> i64 {
    lists(h).map(|l| l.nodes.len() as i64).unwrap_or(-1)
}

/// Copies the lists: `offsets` (count + 1 values), `nodes` (total values),
/// and, when not null, `values` (count values; NaN for lists without one).
///
/// # Safety
/// `h` must be a live handle; buffers sized as described.
#[no_mangle]
pub unsafe extern "C" fn net_path_lists_copy(
    h: *const PathListsHandle,
    offsets: *mut i64,
    nodes: *mut i64,
    values: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let l = lists(h)?;
        for (o, &x) in slice_mut(offsets, l.offsets.len())
            .iter_mut()
            .zip(&l.offsets)
        {
            *o = x as i64;
        }
        for (o, &x) in slice_mut(nodes, l.nodes.len()).iter_mut().zip(&l.nodes) {
            *o = x as i64;
        }
        if !values.is_null() {
            for (i, o) in slice_mut(values, l.len()).iter_mut().enumerate() {
                *o = l.values.get(i).copied().unwrap_or(f64::NAN);
            }
        }
        Ok(())
    })
}

/// Releases a list collection. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_path_lists_free(h: *mut PathListsHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Shortest paths from the nearest of `n_sources` sources (lengths must be
/// non-negative): `dist` (n doubles; Inf when not reached) and `parent`
/// (n values; -1 for sources and unreached nodes), and, when not null,
/// `parent_edge` (the edge from the parent) and `origin` (the source each
/// node was reached from). `cutoff` bounds the distances (Inf for none),
/// `target >= 0` stops the search once the target is settled, and
/// `reverse != 0` follows arcs backward (paths to the sources).
///
/// # Safety
/// `pg` must be a live handle; `sources` must hold `n_sources` values;
/// output buffers n values each (or null where allowed).
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_shortest_paths(
    pg: *const PathGraphHandle,
    sources: *const i64,
    n_sources: i64,
    cutoff: f64,
    target: i64,
    reverse: i32,
    dist: *mut f64,
    parent: *mut i64,
    parent_edge: *mut i64,
    origin: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let n = pg.node_count();
        let sources = node_list(pg, sources, n_sources, "source")?.unwrap_or_default();
        let target = if target < 0 {
            None
        } else {
            Some(node(pg, target, "target")?)
        };
        let mut ws = Searcher::new(n);
        ws.run(
            pg,
            &sources,
            SearchOptions {
                cutoff,
                target,
                reverse: reverse != 0,
            },
        )?;
        let t = ws.tree(pg);
        slice_mut(dist, n).copy_from_slice(&t.dist);
        for (o, &p) in slice_mut(parent, n).iter_mut().zip(&t.parent) {
            *o = index(p);
        }
        if !parent_edge.is_null() {
            for (o, &p) in slice_mut(parent_edge, n).iter_mut().zip(&t.parent_edge) {
                *o = index(p);
            }
        }
        if !origin.is_null() {
            for (o, &p) in slice_mut(origin, n).iter_mut().zip(&t.origin) {
                *o = index(p);
            }
        }
        Ok(())
    })
}

/// Predecessor lists of one search from the sources (as for
/// `net_shortest_paths`): list v holds every node u with an arc u -> v on
/// a shortest path, in the order the predecessors were settled. `dist`
/// (n doubles) receives the distances. Release `*out` with
/// `net_path_lists_free`.
///
/// # Safety
/// `pg` must be a live handle; `sources` must hold `n_sources` values;
/// `dist` n doubles; `out` writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_shortest_path_predecessors(
    pg: *const PathGraphHandle,
    sources: *const i64,
    n_sources: i64,
    cutoff: f64,
    reverse: i32,
    tie: i32,
    dist: *mut f64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let n = pg.node_count();
        let tie = Tie::from_code(tie)?;
        let sources = node_list(pg, sources, n_sources, "source")?.unwrap_or_default();
        let mut ws = Searcher::new(n);
        let opts = SearchOptions {
            cutoff,
            target: None,
            reverse: reverse != 0,
        };
        ws.run(pg, &sources, opts)?;
        let d = slice_mut(dist, n);
        for (v, x) in d.iter_mut().enumerate() {
            *x = ws.dist(v);
        }
        emit(out, ws.predecessors(pg, tie, reverse != 0))
    })
}

/// Distances and parents from every source (every node when `sources` is
/// null or `n_sources < 0`), as column-major matrices with one row per
/// source and one column per node; lengths must be non-negative. `parent`
/// may be null. `cutoff` bounds the distances.
///
/// # Safety
/// `pg` must be a live handle; `sources` null or `n_sources` values;
/// `dist` and `parent` (number of sources) * n values.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_all_pairs_shortest_paths(
    pg: *const PathGraphHandle,
    sources: *const i64,
    n_sources: i64,
    cutoff: f64,
    threads: i32,
    dist: *mut f64,
    parent: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let n = pg.node_count();
        let sources = node_list(pg, sources, n_sources, "source")?;
        let (d, p) = all_pairs(pg, sources.as_deref(), cutoff, thread_count(threads))?;
        let rows = d.len().checked_div(n).unwrap_or(0);
        transpose_into(&d, rows, n, slice_mut(dist, rows * n), |x| x);
        if !parent.is_null() {
            transpose_into(&p, rows, n, slice_mut(parent, rows * n), index);
        }
        Ok(())
    })
}

/// One shortest path from `s` to `t` by bidirectional search (lengths must
/// be non-negative): `*out` receives zero lists (no path) or one list (the
/// path, with its length as value). Release with `net_path_lists_free`.
///
/// # Safety
/// `pg` must be a live handle; `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_bidirectional_shortest_path(
    pg: *const PathGraphHandle,
    s: i64,
    t: i64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let found = bidirectional(pg, node(pg, s, "source")?, node(pg, t, "target")?)?;
        let mut l = NodeLists::new();
        if let Some((len, path)) = found {
            l.push_with(&path, len);
        }
        emit(out, l)
    })
}

/// A* search from `s` to `t` with one heuristic value per node (n doubles,
/// the estimated remaining length to `t`); `cutoff` drops nodes whose
/// length plus heuristic exceeds it (Inf for none). `*out` receives zero
/// lists or one path with its length.
///
/// # Safety
/// `pg` must be a live handle; `heuristic` n doubles; `out` writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_astar_path(
    pg: *const PathGraphHandle,
    s: i64,
    t: i64,
    heuristic: *const f64,
    cutoff: f64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let h = slice(heuristic, pg.node_count());
        let found = astar(
            pg,
            node(pg, s, "source")?,
            node(pg, t, "target")?,
            h,
            cutoff,
        )?;
        let mut l = NodeLists::new();
        if let Some((len, path)) = found {
            l.push_with(&path, len);
        }
        emit(out, l)
    })
}

/// Shortest paths with lengths of any sign from the nearest of the sources
/// (`method` 1 needs exactly one source). Without a negative cycle
/// reachable from the sources: `dist` and `parent` (n values each) are
/// filled, `*preds` (when `preds` is not null) receives the predecessor
/// lists, and `*cycle` receives zero lists. With one: `*cycle` receives one
/// list (the cycle as a closed walk) and nothing else is written.
///
/// # Safety
/// `pg` must be a live handle; `sources` must hold `n_sources` values;
/// `dist` and `parent` n values; `cycle` writable; `preds` null or
/// writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_negative_shortest_paths(
    pg: *const PathGraphHandle,
    sources: *const i64,
    n_sources: i64,
    method: i32,
    tie: i32,
    dist: *mut f64,
    parent: *mut i64,
    preds: *mut *mut PathListsHandle,
    cycle: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let n = pg.node_count();
        let tie = Tie::from_code(tie)?;
        let sources = node_list(pg, sources, n_sources, "source")?.unwrap_or_default();
        let outcome = match method {
            0 => bellman_ford(pg, &sources)?,
            1 => {
                if sources.len() != 1 {
                    return Err(Error::InvalidArgument(
                        "Goldberg and Radzik's algorithm takes exactly one source".into(),
                    ));
                }
                goldberg_radzik(pg, sources[0])?
            }
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown method code {method} for negative lengths"
                )))
            }
        };
        match outcome {
            Outcome::NegativeCycle(c) => {
                let mut l = NodeLists::new();
                l.push(&c);
                emit(cycle, l)
            }
            Outcome::Tree(t) => {
                slice_mut(dist, n).copy_from_slice(&t.dist);
                for (o, &p) in slice_mut(parent, n).iter_mut().zip(&t.parent) {
                    *o = index(p);
                }
                if !preds.is_null() {
                    emit(preds, tree_predecessors(pg, &t, tie))?;
                }
                emit(cycle, NodeLists::new())
            }
        }
    })
}

/// A negative cycle reachable from `source`, or anywhere when
/// `source < 0`: `*out` receives zero lists (none) or one closed walk.
///
/// # Safety
/// `pg` must be a live handle; `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_negative_cycle(
    pg: *const PathGraphHandle,
    source: i64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let found = if source < 0 {
            has_negative_cycle(pg)?
        } else {
            find_negative_cycle(pg, node(pg, source, "source")?)?
        };
        let mut l = NodeLists::new();
        if let Some(c) = found {
            l.push(&c);
        }
        emit(out, l)
    })
}

/// All-pairs distances and parents by Johnson's algorithm (lengths of any
/// sign), as for `net_all_pairs_shortest_paths`. `*cycle` receives zero
/// lists, or one negative cycle (and the matrices are not written).
///
/// # Safety
/// `pg` must be a live handle; `sources` null or `n_sources` values;
/// `dist` and `parent` (number of sources) * n values; `cycle` writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_johnson(
    pg: *const PathGraphHandle,
    sources: *const i64,
    n_sources: i64,
    threads: i32,
    dist: *mut f64,
    parent: *mut i64,
    cycle: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let n = pg.node_count();
        let sources = node_list(pg, sources, n_sources, "source")?;
        match johnson(pg, sources.as_deref(), thread_count(threads))? {
            Err(c) => {
                let mut l = NodeLists::new();
                l.push(&c);
                emit(cycle, l)
            }
            Ok((d, p)) => {
                let rows = d.len().checked_div(n).unwrap_or(0);
                transpose_into(&d, rows, n, slice_mut(dist, rows * n), |x| x);
                if !parent.is_null() {
                    transpose_into(&p, rows, n, slice_mut(parent, rows * n), index);
                }
                emit(cycle, NodeLists::new())
            }
        }
    })
}

/// All-pairs distances and predecessors by the Floyd-Warshall algorithm:
/// `dist` and `pred` are n-by-n column-major matrices (row u, column v:
/// the distance from u to v and the node before v on that path, -1 when
/// none); `*negative_cycle` is set to 1 when the network has a negative
/// cycle (the matrices are then not shortest-path lengths), else 0.
///
/// # Safety
/// `pg` must be a live handle; `dist` and `pred` n * n values;
/// `negative_cycle` one value.
#[no_mangle]
pub unsafe extern "C" fn net_floyd_warshall(
    pg: *const PathGraphHandle,
    threads: i32,
    dist: *mut f64,
    pred: *mut i64,
    negative_cycle: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let n = pg.node_count();
        let r = floyd_warshall(pg, thread_count(threads));
        transpose_into(&r.dist, n, n, slice_mut(dist, n * n), |x| x);
        if !pred.is_null() {
            transpose_into(&r.pred, n, n, slice_mut(pred, n * n), index);
        }
        slice_mut(negative_cycle, 1)[0] = r.negative_cycle as i32;
        Ok(())
    })
}

/// Number of shortest paths (node sequences) from `s` to every node; `out`
/// holds n doubles.
///
/// # Safety
/// `pg` must be a live handle; `out` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_count_shortest_paths(
    pg: *const PathGraphHandle,
    s: i64,
    tie: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let sigma = count_shortest_paths(pg, node(pg, s, "source")?, Tie::from_code(tie)?)?;
        slice_mut(out, pg.node_count()).copy_from_slice(&sigma);
        Ok(())
    })
}

/// Every shortest path from `s` to `t`, in lexicographic order, each with
/// its length; zero lists when `t` is unreachable. Fails when there are
/// more than `max_paths`.
///
/// # Safety
/// `pg` must be a live handle; `out` writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_all_shortest_paths(
    pg: *const PathGraphHandle,
    s: i64,
    t: i64,
    tie: i32,
    max_paths: i64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let l = all_shortest_paths(
            pg,
            node(pg, s, "source")?,
            node(pg, t, "target")?,
            Tie::from_code(tie)?,
            max_paths.max(0) as usize,
        )?;
        emit(out, l)
    })
}

/// The `k` shortest simple paths from `s` to `t` in order of length (Yen's
/// algorithm), each with its length.
///
/// # Safety
/// `pg` must be a live handle; `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_k_shortest_paths(
    pg: *const PathGraphHandle,
    s: i64,
    t: i64,
    k: i64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let l = k_shortest_simple_paths(
            pg,
            node(pg, s, "source")?,
            node(pg, t, "target")?,
            k.max(0) as usize,
        )?;
        emit(out, l)
    })
}

/// A shortest path from `s` to `t` that avoids `n_nodes` nodes and
/// `n_edges` edges (input indices): zero lists or one path with its
/// length.
///
/// # Safety
/// `pg` must be a live handle; `nodes` and `edges` hold their counts;
/// `out` writable.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_shortest_path_avoiding(
    pg: *const PathGraphHandle,
    s: i64,
    t: i64,
    nodes: *const i64,
    n_nodes: i64,
    edges: *const i64,
    n_edges: i64,
    out: *mut *mut PathListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let pg = path_graph(pg)?;
        let avoid = node_list(pg, nodes, n_nodes, "avoided")?.unwrap_or_default();
        let mut avoid_edges = Vec::new();
        if !edges.is_null() && n_edges > 0 {
            for &e in slice(edges, n_edges as usize) {
                if e < 0 || e as u64 >= pg.edge_count() as u64 {
                    return Err(Error::InvalidArgument(format!(
                        "edge index {e} is out of range for a network with {} edges",
                        pg.edge_count()
                    )));
                }
                avoid_edges.push(e as u32);
            }
        }
        let found = shortest_path_avoiding(
            pg,
            node(pg, s, "source")?,
            node(pg, t, "target")?,
            &avoid,
            &avoid_edges,
        )?;
        let mut l = NodeLists::new();
        if let Some((len, path)) = found {
            l.push_with(&path, len);
        }
        emit(out, l)
    })
}
