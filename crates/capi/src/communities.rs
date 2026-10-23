//! C ABI for community detection: Walktrap, label propagation, Infomap
//! and the map equation, Leiden, spinglass, exact maximum modularity,
//! fluid communities, Voronoi communities, partition comparison, and
//! dendrograms (conversion and cutting).
//!
//! Memberships are written as 0-based labels numbered by first
//! appearance. Dendrograms use igraph's `merges` numbering: leaves 0..n,
//! merge i creates cluster n + i; a dendrogram is passed as two arrays
//! `merge_a`, `merge_b` of `n_merges` values.
//!
//! Integer codes used by these functions:
//! - `variant` (label propagation): 0 = dominance, 1 = retention,
//!   2 = fast, 3 = semi-synchronous.
//! - `mode` (label propagation, Voronoi): 1 = out, 2 = in, 3 = all.
//! - `objective` (Leiden): 0 = Constant Potts Model, 1 = modularity.
//! - `null_model` (spinglass): 0 = simple, 1 = config.
//! - `method` (comparison): 0 = variation of information, 1 = normalized
//!   mutual information, 2 = split-join distance, 3 = Rand index,
//!   4 = adjusted Rand index, 5 = adjusted mutual information.

use std::ffi::c_char;

use engine::metrics::partition::{
    compare_partitions, cut_at_groups, cut_at_max_modularity, cut_dendrogram,
    dendrogram_modularity, fluid_communities, infomap, label_propagation, leiden, map_equation,
    merges_from_edge_removals, merges_from_node_pairs, optimal_modularity, spinglass,
    split_join_distance, voronoi, walktrap, CompareMethod, LabelPropagationOptions,
    LeidenObjective, LeidenOptions, LpaMode, LpaVariant, SpinglassNull, SpinglassOptions,
    VoronoiMode,
};
use engine::Error;

use super::{count, guard, handle, slice, slice_mut, thread_count, GraphHandle};

fn write_membership(out: &mut [i64], membership: &[usize]) {
    for (o, &c) in out.iter_mut().zip(membership) {
        *o = c as i64;
    }
}

/// An optional array: `None` when the pointer is null.
///
/// # Safety
/// `ptr` must be null or point to `len` readable elements.
unsafe fn optional<'a, T>(ptr: *const T, len: usize) -> Option<&'a [T]> {
    if ptr.is_null() {
        None
    } else {
        Some(slice(ptr, len))
    }
}

/// A dendrogram from two arrays of cluster indices.
///
/// # Safety
/// `a` and `b` must point to `len` readable values (or be null with
/// `len` 0).
unsafe fn merges(a: *const i64, b: *const i64, len: i64) -> Result<Vec<(usize, usize)>, Error> {
    let len = count(len, "the number of merges")?;
    let (a, b) = (slice(a, len), slice(b, len));
    a.iter()
        .zip(b)
        .map(|(&x, &y)| {
            Ok((
                count(x, "a merge's cluster index")?,
                count(y, "a merge's cluster index")?,
            ))
        })
        .collect()
}

fn write_merges(list: &[(usize, usize)], a: &mut [i64], b: &mut [i64]) {
    for (i, &(x, y)) in list.iter().enumerate() {
        a[i] = x as i64;
        b[i] = y as i64;
    }
}

/// Walktrap communities of an undirected network with walks of `steps`
/// steps. Writes the labels of the level of highest modularity to
/// `membership` (n values) and its modularity to `*q`; the merges to
/// `merge_a` / `merge_b` (n - 1 values each, of which `*n_merges` are
/// used) and the modularity before any merge and after each to `path`
/// (n values, of which `*n_merges + 1` are used).
///
/// # Safety
/// `h` must be a live handle; buffers must have the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_walktrap(
    h: *const GraphHandle,
    steps: i64,
    weighted: i32,
    threads_: i32,
    membership: *mut i64,
    q: *mut f64,
    merge_a: *mut i64,
    merge_b: *mut i64,
    path: *mut f64,
    n_merges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let r = walktrap(
            g,
            count(steps, "the number of steps")?,
            weighted != 0,
            thread_count(threads_),
        )?;
        write_membership(slice_mut(membership, n), &r.membership);
        slice_mut(q, 1)[0] = r.modularity;
        let cap = n.saturating_sub(1);
        write_merges(&r.merges, slice_mut(merge_a, cap), slice_mut(merge_b, cap));
        slice_mut(path, n)[..r.path.len()].copy_from_slice(&r.path);
        slice_mut(n_merges, 1)[0] = r.merges.len() as i64;
        Ok(())
    })
}

/// Label propagation communities. `initial` (n labels, negative =
/// unlabeled) and `fixed` (n flags, nonzero = fixed) may be null. Writes
/// labels to `membership` (n values).
///
/// # Safety
/// `h` must be a live handle; `initial` and `fixed` must be null or hold
/// n values; `membership` must hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_label_propagation(
    h: *const GraphHandle,
    variant: i32,
    mode: i32,
    weighted: i32,
    initial: *const i64,
    fixed: *const i32,
    seed: u64,
    membership: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let fixed: Option<Vec<bool>> =
            optional(fixed, n).map(|f| f.iter().map(|&x| x != 0).collect());
        let opts = LabelPropagationOptions {
            variant: LpaVariant::from_code(variant)?,
            mode: LpaMode::from_code(mode)?,
            weighted: weighted != 0,
            initial: optional(initial, n),
            fixed: fixed.as_deref(),
            seed,
        };
        let m = label_propagation(g, &opts)?;
        write_membership(slice_mut(membership, n), &m);
        Ok(())
    })
}

/// Infomap communities: the best of `trials` searches. `node_weights`
/// (n values, teleportation weights for directed networks) may be null.
/// Writes labels to `membership` (n values) and the codelength in bits to
/// `*codelength`.
///
/// # Safety
/// `h` must be a live handle; `node_weights` must be null or hold n
/// values; `membership` must hold n values; `codelength` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_infomap(
    h: *const GraphHandle,
    weighted: i32,
    node_weights: *const f64,
    trials: i64,
    seed: u64,
    membership: *mut i64,
    codelength: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let r = infomap(
            g,
            weighted != 0,
            optional(node_weights, n),
            count(trials, "the number of trials")?,
            seed,
        )?;
        write_membership(slice_mut(membership, n), &r.membership);
        slice_mut(codelength, 1)[0] = r.codelength;
        Ok(())
    })
}

/// The map equation (bits per step) of the partition `membership` (n
/// non-negative labels). `node_weights` may be null.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold `n` values;
/// `node_weights` must be null or hold `n` values; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_map_equation(
    h: *const GraphHandle,
    membership: *const i64,
    n: i64,
    weighted: i32,
    node_weights: *const f64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = count(n, "the membership length")?;
        slice_mut(out, 1)[0] = map_equation(
            g,
            slice(membership, n),
            weighted != 0,
            optional(node_weights, n),
        )?;
        Ok(())
    })
}

/// Leiden communities of an undirected network. `node_weights` and
/// `initial` (n values each) may be null; `n_iterations` < 0 runs until
/// stable. Writes labels to `membership` (n values), the quality to
/// `*quality`, and the number of iterations run to `*iterations`.
///
/// # Safety
/// `h` must be a live handle; optional arrays must be null or hold n
/// values; outputs must be writable with the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_leiden(
    h: *const GraphHandle,
    objective: i32,
    resolution: f64,
    beta: f64,
    n_iterations: i64,
    weighted: i32,
    node_weights: *const f64,
    initial: *const i64,
    seed: u64,
    membership: *mut i64,
    quality: *mut f64,
    iterations: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let opts = LeidenOptions {
            objective: LeidenObjective::from_code(objective)?,
            resolution,
            beta,
            n_iterations,
            weighted: weighted != 0,
            node_weights: optional(node_weights, n),
            initial: optional(initial, n),
            seed,
        };
        let r = leiden(g, &opts)?;
        write_membership(slice_mut(membership, n), &r.membership);
        slice_mut(quality, 1)[0] = r.quality;
        slice_mut(iterations, 1)[0] = r.iterations as i64;
        Ok(())
    })
}

/// Spinglass communities of a connected network by simulated annealing.
/// Writes labels to `membership` (n values), the modularity to
/// `*modularity`, and the final energy to `*energy`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n values;
/// `modularity` and `energy` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_spinglass(
    h: *const GraphHandle,
    spins: i64,
    parallel: i32,
    start_temp: f64,
    stop_temp: f64,
    cool_fact: f64,
    null_model: i32,
    gamma: f64,
    signed_: i32,
    gamma_minus: f64,
    weighted: i32,
    seed: u64,
    membership: *mut i64,
    modularity: *mut f64,
    energy: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = SpinglassOptions {
            spins: count(spins, "the number of spins")?,
            parallel: parallel != 0,
            start_temp,
            stop_temp,
            cool_fact,
            null: SpinglassNull::from_code(null_model)?,
            gamma,
            signed: signed_ != 0,
            gamma_minus,
            weighted: weighted != 0,
            seed,
        };
        let r = spinglass(g, &opts)?;
        write_membership(slice_mut(membership, g.node_count()), &r.membership);
        slice_mut(modularity, 1)[0] = r.modularity;
        slice_mut(energy, 1)[0] = r.energy;
        Ok(())
    })
}

/// A partition of maximum modularity (exact; components of at most 60
/// nodes). Writes labels to `membership` (n values) and the modularity to
/// `*q`.
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n values; `q` must be
/// writable.
#[no_mangle]
pub unsafe extern "C" fn net_optimal_modularity(
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
        let r = optimal_modularity(g, resolution, weighted != 0)?;
        write_membership(slice_mut(membership, g.node_count()), &r.membership);
        slice_mut(q, 1)[0] = r.modularity;
        Ok(())
    })
}

/// Fluid communities (k groups) of a connected undirected network. Writes
/// labels to `membership` (n values).
///
/// # Safety
/// `h` must be a live handle; `membership` must hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_fluid_communities(
    h: *const GraphHandle,
    k: i64,
    seed: u64,
    membership: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let m = fluid_communities(g, count(k, "the number of communities")?, seed)?;
        write_membership(slice_mut(membership, g.node_count()), &m);
        Ok(())
    })
}

/// Voronoi communities. `lengths` (one per edge) may be null for unit
/// lengths; a negative or NaN `radius` chooses the radius of highest
/// modularity. Writes labels to `membership` (n values), the generators
/// to `generators` (n values, of which `*n_generators` are used), the
/// modularity to `*modularity`, and the radius used to `*radius_used`.
///
/// # Safety
/// `h` must be a live handle; `lengths` must be null or hold m values;
/// outputs must be writable with the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_voronoi(
    h: *const GraphHandle,
    lengths: *const f64,
    weighted: i32,
    mode: i32,
    radius: f64,
    threads_: i32,
    membership: *mut i64,
    generators: *mut i64,
    n_generators: *mut i64,
    modularity: *mut f64,
    radius_used: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let r = voronoi(
            g,
            optional(lengths, g.edge_count()),
            weighted != 0,
            VoronoiMode::from_code(mode)?,
            (radius >= 0.0).then_some(radius),
            thread_count(threads_),
        )?;
        write_membership(slice_mut(membership, n), &r.membership);
        write_membership(slice_mut(generators, n), &r.generators);
        slice_mut(n_generators, 1)[0] = r.generators.len() as i64;
        slice_mut(modularity, 1)[0] = r.modularity;
        slice_mut(radius_used, 1)[0] = r.radius;
        Ok(())
    })
}

/// Compares two partitions of `n` nodes (arbitrary integer labels) with
/// `method`; writes the value to `*out`.
///
/// # Safety
/// `x` and `y` must hold `n` values; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_compare_partitions(
    x: *const i64,
    y: *const i64,
    n: i64,
    method: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "the partition length")?;
        slice_mut(out, 1)[0] =
            compare_partitions(slice(x, n), slice(y, n), CompareMethod::from_code(method)?)?;
        Ok(())
    })
}

/// The two projection distances of van Dongen: of `x` from `y` to
/// `*x_to_y` and of `y` from `x` to `*y_to_x`.
///
/// # Safety
/// `x` and `y` must hold `n` values; outputs must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_split_join_distance(
    x: *const i64,
    y: *const i64,
    n: i64,
    x_to_y: *mut i64,
    y_to_x: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "the partition length")?;
        let (a, b) = split_join_distance(slice(x, n), slice(y, n))?;
        slice_mut(x_to_y, 1)[0] = a as i64;
        slice_mut(y_to_x, 1)[0] = b as i64;
        Ok(())
    })
}

/// Membership of `n` leaves after the first `steps` merges of a
/// dendrogram, or (when `groups` > 0) cut at exactly `groups` groups.
///
/// # Safety
/// `merge_a` / `merge_b` must hold `n_merges` values; `membership` must
/// hold `n` values.
#[no_mangle]
pub unsafe extern "C" fn net_cut_dendrogram(
    n: i64,
    merge_a: *const i64,
    merge_b: *const i64,
    n_merges: i64,
    steps: i64,
    groups: i64,
    membership: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "the number of leaves")?;
        let list = merges(merge_a, merge_b, n_merges)?;
        let m = if groups > 0 {
            cut_at_groups(n, &list, groups as usize)?
        } else {
            cut_dendrogram(n, &list, count(steps, "the number of steps")?)?
        };
        write_membership(slice_mut(membership, n), &m);
        Ok(())
    })
}

/// Modularity of an undirected network after 0, 1, ..., `n_merges`
/// merges of a dendrogram on its nodes, to `path` (`n_merges + 1`
/// values); the labels of the first level of highest modularity to
/// `membership` (n values) and its modularity to `*q`.
///
/// # Safety
/// `h` must be a live handle; buffers must have the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_dendrogram_modularity(
    h: *const GraphHandle,
    merge_a: *const i64,
    merge_b: *const i64,
    n_merges: i64,
    resolution: f64,
    weighted: i32,
    path: *mut f64,
    membership: *mut i64,
    q: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let list = merges(merge_a, merge_b, n_merges)?;
        let p = dendrogram_modularity(g, &list, resolution, weighted != 0)?;
        slice_mut(path, list.len() + 1).copy_from_slice(&p);
        let (m, best) = cut_at_max_modularity(g, &list, resolution, weighted != 0)?;
        write_membership(slice_mut(membership, g.node_count()), &m);
        slice_mut(q, 1)[0] = best;
        Ok(())
    })
}

/// Converts a merge list whose groups are named by any of their nodes
/// (`node_a`, `node_b`; `len` pairs, as `net_fast_greedy`-style methods
/// report them) or, when `removals` is nonzero, a list of removed edges in
/// order of removal (as `net_girvan_newman` reports them) into a
/// dendrogram on `n` leaves: `merge_a` / `merge_b` receive up to n - 1
/// merges and `*n_merges` their number.
///
/// # Safety
/// `node_a` / `node_b` must hold `len` values; `merge_a` / `merge_b` must
/// hold n - 1 values; `n_merges` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_merges_from_pairs(
    n: i64,
    node_a: *const i64,
    node_b: *const i64,
    len: i64,
    removals: i32,
    merge_a: *mut i64,
    merge_b: *mut i64,
    n_merges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "the number of leaves")?;
        let len = count(len, "the number of pairs")?;
        let node = |x: i64| -> Result<u32, Error> {
            if x < 0 || x as u64 >= n as u64 {
                Err(Error::InvalidArgument(format!(
                    "node {x} is out of range for {n} nodes"
                )))
            } else {
                Ok(x as u32)
            }
        };
        let pairs: Vec<(u32, u32)> = slice(node_a, len)
            .iter()
            .zip(slice(node_b, len))
            .map(|(&a, &b)| Ok((node(a)?, node(b)?)))
            .collect::<Result<_, Error>>()?;
        let list = if removals != 0 {
            merges_from_edge_removals(n, &pairs)?
        } else {
            merges_from_node_pairs(n, &pairs)?
        };
        let cap = n.saturating_sub(1);
        write_merges(&list, slice_mut(merge_a, cap), slice_mut(merge_b, cap));
        slice_mut(n_merges, 1)[0] = list.len() as i64;
        Ok(())
    })
}
