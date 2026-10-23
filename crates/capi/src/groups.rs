//! C ABI for cohesive subsets (bicomponents, n-cliques, n-clans, k-plexes,
//! lambda sets, factions), community partitions (Girvan-Newman, leading
//! eigenvector, greedy modularity), topography (core-periphery,
//! compactness, small worlds), and random graphs.
//!
//! Collections of node sets come back in a `NodeSetsHandle` and batches of
//! random graphs in an `EdgeListsHandle`: ask for the sizes, allocate,
//! copy, and free the handle, as with `GroupsHandle`.
//!
//! Integer codes used by these functions:
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength.
//! - `fit` (factions): 0 = errors, 1 = correlation.
//! - `links` (core-periphery): 0 = missing, 1 = ones, 2 = zeros.

use std::ffi::c_char;

use engine::generate::{self, Degrees, EdgeList};
use engine::metrics::coreperiphery::{
    compactness, continuous_core_periphery, core_periphery, small_world, ContinuousFit,
    ContinuousOptions, CoreLinks, CorePeripheryOptions, SmallWorldOptions,
};
use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::PathOptions;
use engine::metrics::factions::{
    bicomponents, factions, k_plexes, lambda_connectivity, n_clans, n_cliques, FactionFit,
    FactionOptions, KPlexOptions,
};
use engine::metrics::partition::{fast_greedy, girvan_newman, leading_eigen};
use engine::{Error, WeightKind};

use crate::{guard, handle, slice, slice_mut, GraphHandle};

/// Opaque collection of node sets.
pub struct NodeSetsHandle(Vec<Vec<u32>>);

/// Opaque batch of edge lists.
pub struct EdgeListsHandle(Vec<EdgeList>);

fn count(value: i64, what: &str) -> Result<usize, Error> {
    usize::try_from(value)
        .map_err(|_| Error::InvalidArgument(format!("{what} must not be negative")))
}

fn threads(t: i32) -> usize {
    t.max(0) as usize
}

/// # Safety
/// `out` must be writable.
unsafe fn give_sets(out: *mut *mut NodeSetsHandle, sets: Vec<Vec<u32>>) -> Result<(), Error> {
    if out.is_null() {
        return Err(Error::InvalidArgument("output pointer is null".into()));
    }
    *out = Box::into_raw(Box::new(NodeSetsHandle(sets)));
    Ok(())
}

/// # Safety
/// `out` must be writable.
unsafe fn give_edges(
    out: *mut *mut EdgeListsHandle,
    lists: Result<Vec<EdgeList>, Error>,
) -> Result<(), Error> {
    if out.is_null() {
        return Err(Error::InvalidArgument("output pointer is null".into()));
    }
    *out = Box::into_raw(Box::new(EdgeListsHandle(lists?)));
    Ok(())
}

/// Number of node sets, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live node sets handle.
#[no_mangle]
pub unsafe extern "C" fn net_node_sets_count(h: *const NodeSetsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.len() as i64
}

/// Total number of members over all node sets, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live node sets handle.
#[no_mangle]
pub unsafe extern "C" fn net_node_sets_total(h: *const NodeSetsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.iter().map(|s| s.len() as i64).sum()
}

/// Copies the node sets out: members of set `i` are
/// `members[offsets[i] .. offsets[i + 1]]`, 0-based and increasing.
/// `offsets` must hold `count + 1` values and `members` `total` values.
///
/// # Safety
/// `h` must be a live node sets handle; the arrays must have the sizes
/// above.
#[no_mangle]
pub unsafe extern "C" fn net_node_sets_copy(
    h: *const NodeSetsHandle,
    offsets: *mut i64,
    members: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if h.is_null() {
            return Err(Error::InvalidArgument("node sets handle is null".into()));
        }
        let sets = &(*h).0;
        let total: usize = sets.iter().map(|s| s.len()).sum();
        let offsets = slice_mut(offsets, sets.len() + 1);
        let members = slice_mut(members, total);
        offsets[0] = 0;
        let mut at = 0;
        for (i, set) in sets.iter().enumerate() {
            for &v in set {
                members[at] = v as i64;
                at += 1;
            }
            offsets[i + 1] = at as i64;
        }
        Ok(())
    })
}

/// Releases a node sets handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live node sets handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_node_sets_free(h: *mut NodeSetsHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Bicomponents (blocks) of the undirected simple reading of the network.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_bicomponents(
    h: *const GraphHandle,
    out: *mut *mut NodeSetsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        give_sets(out, bicomponents(g))
    })
}

/// n-cliques (`clans == 0`) or n-clans (`clans != 0`) with at least
/// `min_size` members.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_n_cliques(
    h: *const GraphHandle,
    order: i64,
    min_size: i64,
    clans: i32,
    threads_: i32,
    out: *mut *mut NodeSetsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let order = count(order, "n")?;
        let min_size = count(min_size, "the minimum size")?;
        let sets = if clans != 0 {
            n_clans(g, order, min_size, threads(threads_))?
        } else {
            n_cliques(g, order, min_size, threads(threads_))?
        };
        give_sets(out, sets)
    })
}

/// Maximal k-plexes with at least `min_size` members; stops with an error
/// after `max_steps` search steps.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_k_plexes(
    h: *const GraphHandle,
    k: i64,
    min_size: i64,
    max_steps: f64,
    threads_: i32,
    out: *mut *mut NodeSetsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if max_steps.is_nan() || max_steps < 1.0 {
            return Err(Error::InvalidArgument(
                "the step limit must be at least 1".into(),
            ));
        }
        let opts = KPlexOptions {
            k: count(k, "k")?,
            min_size: count(min_size, "the minimum size")?,
            max_steps: max_steps.min(u64::MAX as f64) as u64,
            threads: threads(threads_),
        };
        give_sets(out, k_plexes(g, opts)?)
    })
}

/// Edge connectivity between every pair of nodes: an n-by-n symmetric
/// matrix (zero diagonal) written to `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold n * n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_lambda_connectivity(
    h: *const GraphHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let m = lambda_connectivity(g);
        slice_mut(out, n * n).copy_from_slice(&m.to_col_major());
        Ok(())
    })
}

/// Factions: a partition into `k` groups. Writes 0-based labels to
/// `membership`, and the number of errors and the correlation of the
/// partition. `tenure == 0` chooses the default tabu tenure.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n values; `errors`
/// and `correlation` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_factions(
    h: *const GraphHandle,
    k: i64,
    fit: i32,
    starts: i64,
    max_stall: i64,
    tenure: i64,
    seed: u64,
    threads_: i32,
    membership: *mut i64,
    errors: *mut f64,
    correlation: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = FactionOptions {
            k: count(k, "the number of factions")?,
            fit: FactionFit::from_code(fit)?,
            starts: count(starts, "the number of starts")?,
            max_stall: count(max_stall, "the stall limit")?,
            tenure: count(tenure, "the tabu tenure")?,
            seed,
            threads: threads(threads_),
        };
        let f = factions(g, opts)?;
        for (o, &c) in slice_mut(membership, g.node_count())
            .iter_mut()
            .zip(&f.membership)
        {
            *o = c as i64;
        }
        slice_mut(errors, 1)[0] = f.errors as f64;
        slice_mut(correlation, 1)[0] = f.correlation;
        Ok(())
    })
}

/// Girvan-Newman clustering of an undirected network. Writes 0-based labels
/// of the partition of highest modularity to `membership` and its
/// modularity to `*q`; the removed edges (smaller node first) to
/// `removed_from` / `removed_to` and the modularity before any removal and
/// after each removal to `path`, with the number of removals in
/// `*removed`. The removal buffers must hold `m` values and `path` `m + 1`,
/// where m is the graph's edge count (an upper bound).
///
/// # Safety
/// `h` must be a live handle; buffers must have the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_girvan_newman(
    h: *const GraphHandle,
    threads_: i32,
    membership: *mut i64,
    q: *mut f64,
    removed_from: *mut i64,
    removed_to: *mut i64,
    path: *mut f64,
    removed: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let r = girvan_newman(g, threads(threads_))?;
        let m = g.edge_count();
        for (o, &c) in slice_mut(membership, g.node_count())
            .iter_mut()
            .zip(&r.membership)
        {
            *o = c as i64;
        }
        slice_mut(q, 1)[0] = r.modularity;
        let (from, to) = (slice_mut(removed_from, m), slice_mut(removed_to, m));
        for (i, &(a, b)) in r.removed.iter().enumerate() {
            from[i] = a as i64;
            to[i] = b as i64;
        }
        slice_mut(path, m + 1)[..r.path.len()].copy_from_slice(&r.path);
        slice_mut(removed, 1)[0] = r.removed.len() as i64;
        Ok(())
    })
}

/// Leading eigenvector communities of an undirected network: 0-based
/// labels to `membership`, modularity to `*q`, number of divisions to
/// `*splits`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n values; `q` and
/// `splits` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_leading_eigen(
    h: *const GraphHandle,
    weighted: i32,
    membership: *mut i64,
    q: *mut f64,
    splits: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let r = leading_eigen(g, weighted != 0)?;
        for (o, &c) in slice_mut(membership, g.node_count())
            .iter_mut()
            .zip(&r.membership)
        {
            *o = c as i64;
        }
        slice_mut(q, 1)[0] = r.modularity;
        slice_mut(splits, 1)[0] = r.splits as i64;
        Ok(())
    })
}

/// Greedy modularity communities (Clauset-Newman-Moore) of an undirected
/// network: 0-based labels to `membership`, modularity to `*q`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n values; `q` must be
/// writable.
#[no_mangle]
pub unsafe extern "C" fn net_fast_greedy(
    h: *const GraphHandle,
    resolution: f64,
    weighted: i32,
    membership: *mut i64,
    q: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let r = fast_greedy(g, resolution, weighted != 0)?;
        for (o, &c) in slice_mut(membership, g.node_count())
            .iter_mut()
            .zip(&r.membership)
        {
            *o = c as i64;
        }
        slice_mut(q, 1)[0] = r.modularity;
        Ok(())
    })
}

/// Discrete core-periphery fit: `core[i]` is 1 for core nodes and 0
/// otherwise; the correlation is written to `*correlation`.
///
/// # Safety
/// `h` must be a live handle; `core` must hold n values; `correlation`
/// must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_core_periphery(
    h: *const GraphHandle,
    links: i32,
    weighted: i32,
    starts: i64,
    max_stall: i64,
    tenure: i64,
    seed: u64,
    threads_: i32,
    core: *mut i32,
    correlation: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = CorePeripheryOptions {
            links: CoreLinks::from_code(links)?,
            weighted: weighted != 0,
            starts: count(starts, "the number of starts")?,
            max_stall: count(max_stall, "the stall limit")?,
            tenure: count(tenure, "the tabu tenure")?,
            seed,
            threads: threads(threads_),
        };
        let r = core_periphery(g, opts)?;
        for (o, &c) in slice_mut(core, g.node_count()).iter_mut().zip(&r.core) {
            *o = c as i32;
        }
        slice_mut(correlation, 1)[0] = r.correlation;
        Ok(())
    })
}

/// Continuous core-periphery fit. `fit`: 0 = MINRES (least squares,
/// deterministic, `tol` the convergence tolerance), 1 = maximum correlation
/// from `starts` starting points (the degree start plus seeded random
/// ones). Writes coreness (non-negative, unit norm) to `coreness`, the
/// correlation of the final pattern to `*correlation`, the residual sum of
/// squares to `*residual` (NaN for the correlation fit), and the iterations
/// (sweeps) made to `*iterations`.
///
/// # Safety
/// `h` must be a live handle; `coreness` must hold n values; the scalar
/// outputs must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_core_periphery_continuous(
    h: *const GraphHandle,
    fit: i32,
    weighted: i32,
    max_iter: i64,
    tol: f64,
    starts: i64,
    seed: u64,
    threads_: i32,
    coreness: *mut f64,
    correlation: *mut f64,
    residual: *mut f64,
    iterations: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = ContinuousOptions {
            fit: ContinuousFit::from_code(fit)?,
            weighted: weighted != 0,
            max_iter: count(max_iter, "max_iter")?,
            tol,
            starts: count(starts, "the number of starts")?,
            seed,
            threads: threads(threads_),
        };
        let r = continuous_core_periphery(g, opts)?;
        slice_mut(coreness, g.node_count()).copy_from_slice(&r.coreness);
        slice_mut(correlation, 1)[0] = r.correlation;
        slice_mut(residual, 1)[0] = r.residual;
        slice_mut(iterations, 1)[0] = r.iterations as i64;
        Ok(())
    })
}

/// Compactness (mean reciprocal distance over ordered pairs; 1 minus
/// distance-weighted fragmentation).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_compactness(
    h: *const GraphHandle,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    threads_: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = PathOptions {
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            threads: threads(threads_),
            ..Default::default()
        };
        slice_mut(out, 1)[0] = compactness(g, opts)?;
        Ok(())
    })
}

/// Small-world quotient. Writes five values to `out`: clustering, mean
/// distance, mean random clustering, mean random distance, and the
/// quotient. `local != 0` uses mean local clustering instead of
/// transitivity.
///
/// # Safety
/// `h` must be a live handle; `out` must hold five doubles.
#[no_mangle]
pub unsafe extern "C" fn net_small_world(
    h: *const GraphHandle,
    local: i32,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = SmallWorldOptions {
            local: local != 0,
            reps: count(reps, "the number of random networks")?,
            seed,
            threads: threads(threads_),
        };
        let r = small_world(g, opts)?;
        slice_mut(out, 5).copy_from_slice(&[
            r.clustering,
            r.mean_distance,
            r.random_clustering,
            r.random_mean_distance,
            r.quotient,
        ]);
        Ok(())
    })
}

/// Number of edge lists, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live edge lists handle.
#[no_mangle]
pub unsafe extern "C" fn net_edge_lists_count(h: *const EdgeListsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.len() as i64
}

/// Total number of edges over all lists, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live edge lists handle.
#[no_mangle]
pub unsafe extern "C" fn net_edge_lists_total(h: *const EdgeListsHandle) -> i64 {
    if h.is_null() {
        return -1;
    }
    (*h).0.iter().map(|e| e.len() as i64).sum()
}

/// Copies the edge lists out: the edges of list `i` are positions
/// `offsets[i] .. offsets[i + 1]` of `from` and `to` (0-based nodes).
/// `offsets` must hold `count + 1` values, `from` and `to` `total` values.
///
/// # Safety
/// `h` must be a live edge lists handle; the arrays must have the sizes
/// above.
#[no_mangle]
pub unsafe extern "C" fn net_edge_lists_copy(
    h: *const EdgeListsHandle,
    offsets: *mut i64,
    from: *mut i64,
    to: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if h.is_null() {
            return Err(Error::InvalidArgument("edge lists handle is null".into()));
        }
        let lists = &(*h).0;
        let total: usize = lists.iter().map(|e| e.len()).sum();
        let offsets = slice_mut(offsets, lists.len() + 1);
        let (from, to) = (slice_mut(from, total), slice_mut(to, total));
        offsets[0] = 0;
        let mut at = 0;
        for (i, e) in lists.iter().enumerate() {
            for (&a, &b) in e.from.iter().zip(&e.to) {
                from[at] = a as i64;
                to[at] = b as i64;
                at += 1;
            }
            offsets[i + 1] = at as i64;
        }
        Ok(())
    })
}

/// Releases an edge lists handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live edge lists handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_edge_lists_free(h: *mut EdgeListsHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// `reps` G(n, p) graphs; replicate r uses random stream (seed, r).
///
/// # Safety
/// `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_gnp(
    n: i64,
    p: f64,
    directed: i32,
    loops: i32,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let (n, reps) = (count(n, "n")?, count(reps, "reps")?);
        let lists = generate::gnp(
            n,
            p,
            directed != 0,
            loops != 0,
            reps,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}

/// `reps` G(n, m) graphs.
///
/// # Safety
/// `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_gnm(
    n: i64,
    m: i64,
    directed: i32,
    loops: i32,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let (n, m, reps) = (count(n, "n")?, count(m, "m")?, count(reps, "reps")?);
        let lists = generate::gnm(
            n,
            m as u64,
            directed != 0,
            loops != 0,
            reps,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}

/// `reps` configuration-model graphs. `degrees` holds n degrees for an
/// undirected graph (`in_degrees` null), or n out-degrees with n
/// in-degrees in `in_degrees` for a directed graph. `simple != 0` repeats
/// each matching up to `max_tries` times until it has no loop or parallel
/// edge.
///
/// # Safety
/// `degrees` (and `in_degrees`, if not null) must hold n values; `out` must
/// be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_configuration(
    n: i64,
    degrees: *const i64,
    in_degrees: *const i64,
    simple: i32,
    max_tries: i64,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "n")?;
        let to_u64 = |p: *const i64| -> Result<Vec<u64>, Error> {
            slice(p, n)
                .iter()
                .map(|&d| {
                    u64::try_from(d)
                        .map_err(|_| Error::InvalidArgument("degrees must not be negative".into()))
                })
                .collect()
        };
        let degs = if in_degrees.is_null() {
            Degrees::Undirected(to_u64(degrees)?)
        } else {
            Degrees::Directed(to_u64(degrees)?, to_u64(in_degrees)?)
        };
        let lists = generate::configuration(
            &degs,
            simple != 0,
            count(max_tries, "max_tries")?,
            count(reps, "reps")?,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}

/// `reps` U|MAN graphs with `mutual` mutual and `asym` asymmetric dyads.
///
/// # Safety
/// `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_uman(
    n: i64,
    mutual: i64,
    asym: i64,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let lists = generate::uman(
            count(n, "n")?,
            count(mutual, "the number of mutual dyads")? as u64,
            count(asym, "the number of asymmetric dyads")? as u64,
            count(reps, "reps")?,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}

/// `reps` Watts-Strogatz small-world graphs.
///
/// # Safety
/// `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_small_world(
    n: i64,
    k: i64,
    p: f64,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let lists = generate::small_world(
            count(n, "n")?,
            count(k, "k")?,
            p,
            count(reps, "reps")?,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}

/// `reps` Barabasi-Albert preferential attachment graphs.
///
/// # Safety
/// `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_preferential(
    n: i64,
    m: i64,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let lists = generate::preferential(
            count(n, "n")?,
            count(m, "m")?,
            count(reps, "reps")?,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}

/// `reps` stochastic block model graphs: `k` blocks of `sizes[r]` nodes
/// (numbered block by block) and tie probabilities `probs[r * k + s]`.
///
/// # Safety
/// `sizes` must hold k values and `probs` k * k; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_random_sbm(
    k: i64,
    sizes: *const i64,
    probs: *const f64,
    directed: i32,
    loops: i32,
    reps: i64,
    seed: u64,
    threads_: i32,
    out: *mut *mut EdgeListsHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let k = count(k, "the number of blocks")?;
        let sizes: Vec<usize> = slice(sizes, k)
            .iter()
            .map(|&s| count(s, "block sizes"))
            .collect::<Result<_, _>>()?;
        let probs = slice(probs, k * k);
        let lists = generate::sbm(
            &sizes,
            probs,
            directed != 0,
            loops != 0,
            count(reps, "reps")?,
            seed,
            threads(threads_),
        );
        give_edges(out, lists)
    })
}
