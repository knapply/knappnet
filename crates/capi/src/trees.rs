//! C ABI for spanning trees, branchings and arborescences, and directed
//! acyclic graphs (topological orders, transitive closure and reduction,
//! longest paths, dominators, lowest common ancestors).
//!
//! Integer codes used by these functions:
//! - `algorithm` (spanning trees): 0 = Kruskal, 1 = Prim, 2 = Boruvka (all
//!   return the same forest).
//! - `kind` (branchings): 0 = maximum branching, 1 = minimum branching,
//!   2 = minimal branching (NetworkX), 3 = maximum spanning arborescence,
//!   4 = minimum spanning arborescence.
//! - `order`: 0 = sorted neighbours, 1 = input order (NetworkX's results
//!   for a graph built from the same edges in the same order).
//! - `reflexive` (closure): 0 = (v, v) for nodes on cycles (NetworkX's
//!   default), 1 = never, 2 = always.
//! - `root` (lowest common ancestors): -2 = DAG, -1 = tree with its own
//!   root, otherwise the root of a tree.

use std::ffi::c_char;

use engine::metrics::dag::{
    all_topological_sorts, dag_longest_path, dominance_frontiers, immediate_dominators,
    is_directed_acyclic_graph, lexicographical_topological_sort, lowest_common_ancestors,
    topological_generations, transitive_closure, transitive_reduction,
    tree_lowest_common_ancestors, Reflexive,
};
use engine::metrics::trees::{optimum_branching, spanning_tree, BranchingKind, SpanningAlgorithm};
use engine::order::NeighborOrder;

use crate::flow::{give_lists, node, ListsHandle};
use crate::{count, guard, handle, slice, slice_mut, GraphHandle};

/// Writes an edge set: positions to `edges`, their number to `size`, and
/// their weight to `weight`.
///
/// # Safety
/// `edges` must hold at least `set.len()` values; `size` and `weight` one.
unsafe fn write_edges(
    set: &engine::metrics::trees::EdgeSet,
    edges: *mut i64,
    size: *mut i64,
    weight: *mut f64,
    cap: usize,
) {
    let out = slice_mut(edges, cap);
    for (i, &e) in set.edges.iter().enumerate() {
        out[i] = e as i64;
    }
    slice_mut(size, 1)[0] = set.edges.len() as i64;
    slice_mut(weight, 1)[0] = set.weight;
}

/// Minimum (`maximum == 0`) or maximum spanning forest of an undirected
/// network: edge positions in increasing (weight, position) order to
/// `edges` (at most n - 1, allocate n), their number to `size`, and their
/// total weight to `weight`.
///
/// # Safety
/// `h` must be a live handle; `edges` must hold n 64-bit integers, `size`
/// one, and `weight` one double.
#[no_mangle]
pub unsafe extern "C" fn net_spanning_tree(
    h: *const GraphHandle,
    weighted: i32,
    maximum: i32,
    algorithm: i32,
    edges: *mut i64,
    size: *mut i64,
    weight: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let t = spanning_tree(
            g,
            weighted != 0,
            maximum != 0,
            SpanningAlgorithm::from_code(algorithm)?,
        )?;
        write_edges(&t, edges, size, weight, g.node_count());
        Ok(())
    })
}

/// Optimum branching or spanning arborescence of a directed network
/// (`root` >= 0 fixes the arborescence's root; -1 for none): edge
/// positions (increasing) to `edges` (allocate n), their number to
/// `size`, and their total weight to `weight`.
///
/// # Safety
/// As for `net_spanning_tree`.
#[no_mangle]
pub unsafe extern "C" fn net_branching(
    h: *const GraphHandle,
    weighted: i32,
    kind: i32,
    root: i64,
    edges: *mut i64,
    size: *mut i64,
    weight: *mut f64,
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
        let b = optimum_branching(g, weighted != 0, BranchingKind::from_code(kind)?, root)?;
        write_edges(&b, edges, size, weight, g.node_count());
        Ok(())
    })
}

/// A topological order of a DAG to `nodes` (n values) and each node's
/// generation (0 for roots) to `generation` (n values indexed by node, may
/// be null).
///
/// # Safety
/// `h` must be a live handle; `nodes` must hold n 64-bit integers and
/// `generation` be null or hold n.
#[no_mangle]
pub unsafe extern "C" fn net_topological_sort(
    h: *const GraphHandle,
    order: i32,
    nodes: *mut i64,
    generation: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let gens = topological_generations(g, NeighborOrder::from_code(order)?)?;
        let mut flat = Vec::with_capacity(n);
        let mut level = vec![0i64; n];
        for (i, gen) in gens.iter().enumerate() {
            for &v in gen {
                flat.push(v as i64);
                level[v as usize] = i as i64;
            }
        }
        slice_mut(nodes, n).copy_from_slice(&flat);
        crate::flow::write_opt(generation, &level);
        Ok(())
    })
}

/// The topological order that always takes the smallest available node,
/// to `nodes` (n values).
///
/// # Safety
/// `h` must be a live handle; `nodes` must hold n 64-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_lexicographical_topological_sort(
    h: *const GraphHandle,
    nodes: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let o: Vec<i64> = lexicographical_topological_sort(g)?
            .into_iter()
            .map(|v| v as i64)
            .collect();
        slice_mut(nodes, g.node_count()).copy_from_slice(&o);
        Ok(())
    })
}

/// Every topological order of a DAG (at most `limit`; 0 = all), one list
/// each (value 0), in NetworkX's sequence for the input order.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_all_topological_sorts(
    h: *const GraphHandle,
    order: i32,
    limit: i64,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let limit = count(limit, "limit")?;
        let all = all_topological_sorts(
            g,
            NeighborOrder::from_code(order)?,
            (limit > 0).then_some(limit),
        )?;
        let values = vec![0.0; all.len()];
        let lists = all
            .into_iter()
            .map(|o| o.into_iter().map(|v| v as i64).collect())
            .collect();
        give_lists(out, lists, values)
    })
}

/// Whether the network is a directed acyclic graph: 1 or 0 to `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one 32-bit integer.
#[no_mangle]
pub unsafe extern "C" fn net_is_dag(
    h: *const GraphHandle,
    out: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = is_directed_acyclic_graph(g) as i32;
        Ok(())
    })
}

/// Pairs (u, v) as one list per node u of its v's (value 0).
fn pairs_by_node(n: usize, pairs: Vec<(u32, u32)>) -> (Vec<Vec<i64>>, Vec<f64>) {
    let mut lists = vec![Vec::new(); n];
    for (u, v) in pairs {
        lists[u as usize].push(v as i64);
    }
    (lists, vec![0.0; n])
}

/// Transitive closure of a directed network: one list per node u holding
/// every v with (u, v) in the closure, sorted (value 0).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_transitive_closure(
    h: *const GraphHandle,
    reflexive: i32,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let pairs = transitive_closure(g, Reflexive::from_code(reflexive)?)?;
        let (lists, values) = pairs_by_node(g.node_count(), pairs);
        give_lists(out, lists, values)
    })
}

/// Transitive reduction of a DAG: one list per node u of the v with
/// u -> v kept, sorted (value 0).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_transitive_reduction(
    h: *const GraphHandle,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let pairs = transitive_reduction(g)?;
        let (lists, values) = pairs_by_node(g.node_count(), pairs);
        give_lists(out, lists, values)
    })
}

/// A longest path of a DAG: its nodes to `path` (allocate n), their number
/// to `size`, and its length to `length`.
///
/// # Safety
/// `h` must be a live handle; `path` must hold n 64-bit integers, `size`
/// one, and `length` one double.
#[no_mangle]
pub unsafe extern "C" fn net_dag_longest_path(
    h: *const GraphHandle,
    weighted: i32,
    order: i32,
    path: *mut i64,
    size: *mut i64,
    length: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let (p, l) = dag_longest_path(g, weighted != 0, NeighborOrder::from_code(order)?)?;
        let out = slice_mut(path, g.node_count());
        for (i, &v) in p.iter().enumerate() {
            out[i] = v as i64;
        }
        slice_mut(size, 1)[0] = p.len() as i64;
        slice_mut(length, 1)[0] = l;
        Ok(())
    })
}

/// Immediate dominator of every node from `start` to `out` (n values; the
/// start's is itself, -1 for nodes unreachable from the start).
///
/// # Safety
/// `h` must be a live handle; `out` must hold n 64-bit integers.
#[no_mangle]
pub unsafe extern "C" fn net_immediate_dominators(
    h: *const GraphHandle,
    start: i64,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let d = immediate_dominators(g, node(start, n, "start")?)?;
        let o: Vec<i64> = d.iter().map(|x| x.map_or(-1, |v| v as i64)).collect();
        slice_mut(out, n).copy_from_slice(&o);
        Ok(())
    })
}

/// Dominance frontiers from `start`: one sorted list per node, with value
/// 1 for nodes reachable from the start and 0 (empty list) otherwise.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_dominance_frontiers(
    h: *const GraphHandle,
    start: i64,
    out: *mut *mut ListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let f = dominance_frontiers(g, node(start, g.node_count(), "start")?)?;
        let values = f
            .iter()
            .map(|x| if x.is_some() { 1.0 } else { 0.0 })
            .collect();
        let lists = f
            .into_iter()
            .map(|x| {
                x.unwrap_or_default()
                    .into_iter()
                    .map(|v| v as i64)
                    .collect()
            })
            .collect();
        give_lists(out, lists, values)
    })
}

/// Lowest common ancestor of each of `k` pairs (`a[i]`, `b[i]`) to `out`
/// (-1 when a DAG pair has no common ancestor). `root`: -2 for a DAG, -1
/// for a tree with its own root, else the tree's root.
///
/// # Safety
/// `h` must be a live handle; `a`, `b`, and `out` must hold k values.
#[no_mangle]
pub unsafe extern "C" fn net_lowest_common_ancestors(
    h: *const GraphHandle,
    a: *const i64,
    b: *const i64,
    k: i64,
    root: i64,
    order: i32,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let k = count(k, "the number of pairs")?;
        let (sa, sb) = (slice(a, k), slice(b, k));
        let mut pairs = Vec::with_capacity(k);
        for i in 0..k {
            pairs.push((node(sa[i], n, "node")?, node(sb[i], n, "node")?));
        }
        let result: Vec<i64> = match root {
            -2 => lowest_common_ancestors(g, &pairs, NeighborOrder::from_code(order)?)?
                .into_iter()
                .map(|x| x.map_or(-1, |v| v as i64))
                .collect(),
            -1 => tree_lowest_common_ancestors(g, None, &pairs)?
                .into_iter()
                .map(|v| v as i64)
                .collect(),
            r => tree_lowest_common_ancestors(g, Some(node(r, n, "root")?), &pairs)?
                .into_iter()
                .map(|v| v as i64)
                .collect(),
        };
        slice_mut(out, k).copy_from_slice(&result);
        Ok(())
    })
}
