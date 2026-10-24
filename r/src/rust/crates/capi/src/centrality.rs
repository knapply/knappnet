//! C ABI for betweenness for sources and targets, percolation and load
//! centrality, current-flow measures, resistance distance, Kemeny's
//! constant, second-order centrality, communicability, dispersion,
//! VoteRank, Laplacian centrality, trophic levels, and reaching
//! centrality.
//!
//! Integer codes used by these functions:
//! - `weights_as`: 0 = distance, 1 = strength. Path-based measures read
//!   weights as lengths; current-flow measures read a distance as a
//!   resistance (conductance 1/w) and a strength as a conductance.
//! - `rescale`: 0 = unordered pairs (igraph), 1 = NetworkX's
//!   `betweenness_centrality` (sampled sources given as `sources`),
//!   2 = NetworkX's `*_subset` functions.
//! - `parallel`: 0 = parallel edges are distinct paths (igraph),
//!   1 = adjacent nodes joined once (NetworkX).
//! - `edge_split`: 0 = proportional to path counts, 1 = NetworkX's subset
//!   edge betweenness.
//! - `tie`: 0 = relative tolerance 1e-10, 1 = exact equality (NetworkX).
//! - Node lists: a null pointer or a negative count means every node.

use std::ffi::c_char;

use engine::metrics::communicability::{communicability, communicability_betweenness};
use engine::metrics::current_flow::{
    current_flow_betweenness, current_flow_betweenness_approx, current_flow_betweenness_subset,
    current_flow_closeness, effective_graph_resistance, kemeny_constant, resistance_distance,
    second_order_centrality, FlowOptions, Sample,
};
use engine::metrics::hierarchy::{
    global_reaching_centrality, local_reaching_centrality, trophic_differences,
    trophic_incoherence, trophic_levels, ReachingOptions,
};
use engine::metrics::path_centrality::{
    betweenness_subset, edge_load, load_centrality, percolation_centrality, EdgeSplit, LoadOptions,
    Parallel, PercolationOptions, Rescale, SubsetOptions,
};
use engine::metrics::spreading::{dispersion, laplacian_centrality, voterank, DispersionOptions};
use engine::paths::Tie;
use engine::{Error, Graph, WeightKind};

use crate::{guard, handle, slice, slice_mut, thread_count, GraphHandle};

/// # Safety
/// `nodes` must be null or hold `k` values.
unsafe fn node_list(g: &Graph, nodes: *const i64, k: i64) -> Result<Option<Vec<u32>>, Error> {
    if nodes.is_null() || k < 0 {
        return Ok(None);
    }
    let n = g.node_count();
    slice(nodes, k as usize)
        .iter()
        .map(|&v| {
            if v < 0 || v as u64 >= n as u64 {
                Err(Error::InvalidArgument(format!(
                    "node index {v} is out of range for a network with {n} nodes"
                )))
            } else {
                Ok(v as u32)
            }
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

/// # Safety
/// `pairs` must be null or hold `2 * k` values.
unsafe fn pair_list(g: &Graph, pairs: *const i64, k: i64) -> Result<Vec<(u32, u32)>, Error> {
    if pairs.is_null() || k <= 0 {
        return Ok(Vec::new());
    }
    let flat = node_list(g, pairs, 2 * k)?.unwrap_or_default();
    Ok(flat.chunks(2).map(|p| (p[0], p[1])).collect())
}

fn flow_options(
    weighted: i32,
    weights_as: i32,
    normalized: i32,
    threads: i32,
) -> Result<FlowOptions, Error> {
    Ok(FlowOptions {
        weighted: weighted != 0,
        weights_as: WeightKind::from_code(weights_as)?,
        normalized: normalized != 0,
        threads: thread_count(threads),
    })
}

/// Betweenness for sources and targets: `nodes` receives n values and, when
/// not null, `edges` one value per edge. See the codes above;
/// `endpoints != 0` counts path endpoints in node values.
///
/// # Safety
/// `h` must be a live handle; `sources` and `targets` null or holding
/// their counts; `nodes` n doubles; `edges` null or m doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_betweenness_subset(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
    targets: *const i64,
    n_targets: i64,
    weighted: i32,
    weights_as: i32,
    endpoints: i32,
    normalized: i32,
    rescale: i32,
    parallel: i32,
    edge_split: i32,
    tie: i32,
    threads: i32,
    nodes: *mut f64,
    edges: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let s = node_list(g, sources, n_sources)?;
        let t = node_list(g, targets, n_targets)?;
        let r = betweenness_subset(
            g,
            SubsetOptions {
                sources: s.as_deref(),
                targets: t.as_deref(),
                weighted: weighted != 0,
                weights_as: WeightKind::from_code(weights_as)?,
                endpoints: endpoints != 0,
                normalized: normalized != 0,
                rescale: Rescale::from_code(rescale)?,
                parallel: Parallel::from_code(parallel)?,
                edge_split: EdgeSplit::from_code(edge_split)?,
                tie: Tie::from_code(tie)?,
                threads: thread_count(threads),
            },
        )?;
        slice_mut(nodes, g.node_count()).copy_from_slice(&r.nodes);
        if !edges.is_null() {
            slice_mut(edges, g.edge_count()).copy_from_slice(&r.edges);
        }
        Ok(())
    })
}

/// Percolation centrality with one percolation state per node (`states`,
/// n doubles); `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `states` and `out` n doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_percolation_centrality(
    h: *const GraphHandle,
    states: *const f64,
    weighted: i32,
    weights_as: i32,
    parallel: i32,
    tie: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        let v = percolation_centrality(
            g,
            slice(states, n),
            PercolationOptions {
                weighted: weighted != 0,
                weights_as: WeightKind::from_code(weights_as)?,
                parallel: Parallel::from_code(parallel)?,
                tie: Tie::from_code(tie)?,
                threads: thread_count(threads),
            },
        )?;
        slice_mut(out, n).copy_from_slice(&v);
        Ok(())
    })
}

/// Load centrality; `cutoff` bounds path lengths (Inf for none) and
/// `networkx_ties != 0` reproduces NetworkX's handling of the source among
/// several predecessors. `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `out` n doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_load_centrality(
    h: *const GraphHandle,
    weighted: i32,
    weights_as: i32,
    normalized: i32,
    cutoff: f64,
    networkx_ties: i32,
    tie: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = load_centrality(
            g,
            LoadOptions {
                weighted: weighted != 0,
                weights_as: WeightKind::from_code(weights_as)?,
                normalized: normalized != 0,
                cutoff,
                networkx_ties: networkx_ties != 0,
                tie: Tie::from_code(tie)?,
                threads: thread_count(threads),
            },
        )?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// NetworkX's edge load: `out` receives 2m values, for each edge the
/// ordered pairs (source, target) and (target, source).
///
/// # Safety
/// `h` must be a live handle; `out` 2m doubles.
#[no_mangle]
pub unsafe extern "C" fn net_edge_load(
    h: *const GraphHandle,
    cutoff: f64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = edge_load(g, cutoff, thread_count(threads))?;
        slice_mut(out, 2 * g.edge_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// Current-flow betweenness: `nodes` n values and, when not null, `edges`
/// m values.
///
/// # Safety
/// `h` must be a live handle; `nodes` n doubles; `edges` null or m doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_current_flow_betweenness(
    h: *const GraphHandle,
    weighted: i32,
    weights_as: i32,
    normalized: i32,
    threads: i32,
    nodes: *mut f64,
    edges: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let r =
            current_flow_betweenness(g, flow_options(weighted, weights_as, normalized, threads)?)?;
        slice_mut(nodes, g.node_count()).copy_from_slice(&r.nodes);
        if !edges.is_null() {
            slice_mut(edges, g.edge_count()).copy_from_slice(&r.edges);
        }
        Ok(())
    })
}

/// Current-flow betweenness for sources and targets (`networkx != 0`:
/// NetworkX's formula); outputs as for `net_current_flow_betweenness`.
///
/// # Safety
/// `h` must be a live handle; `sources` and `targets` hold their counts;
/// `nodes` n doubles; `edges` null or m doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_current_flow_betweenness_subset(
    h: *const GraphHandle,
    sources: *const i64,
    n_sources: i64,
    targets: *const i64,
    n_targets: i64,
    networkx: i32,
    weighted: i32,
    weights_as: i32,
    normalized: i32,
    threads: i32,
    nodes: *mut f64,
    edges: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let all: Vec<u32> = (0..g.node_count() as u32).collect();
        let s = node_list(g, sources, n_sources)?.unwrap_or_else(|| all.clone());
        let t = node_list(g, targets, n_targets)?.unwrap_or(all);
        let r = current_flow_betweenness_subset(
            g,
            &s,
            &t,
            networkx != 0,
            flow_options(weighted, weights_as, normalized, threads)?,
        )?;
        slice_mut(nodes, g.node_count()).copy_from_slice(&r.nodes);
        if !edges.is_null() {
            slice_mut(edges, g.edge_count()).copy_from_slice(&r.edges);
        }
        Ok(())
    })
}

/// Approximate current-flow betweenness from `n_pairs` pairs (`pairs`
/// holds s0, t0, s1, t1, ...), or, when `n_pairs < 0`, from random pairs
/// drawn with `seed`, as many as `epsilon` and `sample_weight` require
/// (an error beyond `kmax`). `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `pairs` null or 2 * `n_pairs` values; `out`
/// n doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_current_flow_betweenness_approx(
    h: *const GraphHandle,
    pairs: *const i64,
    n_pairs: i64,
    epsilon: f64,
    kmax: i64,
    sample_weight: f64,
    seed: u64,
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
        let opts = flow_options(weighted, weights_as, normalized, threads)?;
        let given = pair_list(g, pairs, n_pairs)?;
        let sample = if n_pairs >= 0 {
            Sample::Pairs(&given)
        } else {
            Sample::Random {
                epsilon,
                kmax: kmax.max(0) as usize,
                sample_weight,
                seed,
            }
        };
        let v = current_flow_betweenness_approx(g, sample, opts)?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// Current-flow closeness; `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `out` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_current_flow_closeness(
    h: *const GraphHandle,
    weighted: i32,
    weights_as: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = current_flow_closeness(g, flow_options(weighted, weights_as, 0, threads)?)?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// Resistance distances as an n-by-n matrix (symmetric; Inf between
/// components).
///
/// # Safety
/// `h` must be a live handle; `out` n * n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_resistance_distance(
    h: *const GraphHandle,
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
        let r = resistance_distance(g, flow_options(weighted, weights_as, 0, threads)?)?;
        slice_mut(out, n * n).copy_from_slice(&r.to_col_major());
        Ok(())
    })
}

/// Effective graph resistance (one double).
///
/// # Safety
/// `h` must be a live handle; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_effective_graph_resistance(
    h: *const GraphHandle,
    weighted: i32,
    weights_as: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] =
            effective_graph_resistance(g, flow_options(weighted, weights_as, 0, threads)?)?;
        Ok(())
    })
}

/// Kemeny's constant (one double); weights are used as given.
///
/// # Safety
/// `h` must be a live handle; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_kemeny_constant(
    h: *const GraphHandle,
    weighted: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = kemeny_constant(g, weighted != 0)?;
        Ok(())
    })
}

/// Second-order centrality; `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `out` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_second_order_centrality(
    h: *const GraphHandle,
    weighted: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = second_order_centrality(g, weighted != 0)?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// The communicability matrix (n by n, symmetric).
///
/// # Safety
/// `h` must be a live handle; `out` n * n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_communicability(
    h: *const GraphHandle,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let n = g.node_count();
        slice_mut(out, n * n).copy_from_slice(&communicability(g)?.to_col_major());
        Ok(())
    })
}

/// Communicability betweenness; `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `out` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_communicability_betweenness(
    h: *const GraphHandle,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = communicability_betweenness(g, thread_count(threads))?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// Dispersion of `n_pairs` pairs (`pairs` holds u0, v0, u1, v1, ...);
/// normalized != 0 uses (dispersion + b)^alpha / (embeddedness + c).
/// `out` receives `n_pairs` values.
///
/// # Safety
/// `h` must be a live handle; `pairs` 2 * `n_pairs` values; `out`
/// `n_pairs` doubles.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn net_dispersion(
    h: *const GraphHandle,
    pairs: *const i64,
    n_pairs: i64,
    normalized: i32,
    alpha: f64,
    b: f64,
    c: f64,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let p = pair_list(g, pairs, n_pairs)?;
        let v = dispersion(
            g,
            &p,
            DispersionOptions {
                normalized: normalized != 0,
                alpha,
                b,
                c,
            },
        )?;
        slice_mut(out, p.len()).copy_from_slice(&v);
        Ok(())
    })
}

/// VoteRank: up to `count` nodes (every node when `count < 0`), in the
/// order chosen, written to `out`; `*n_out` receives how many.
///
/// # Safety
/// `h` must be a live handle; `out` n values; `n_out` one value.
#[no_mangle]
pub unsafe extern "C" fn net_voterank(
    h: *const GraphHandle,
    count: i64,
    out: *mut i64,
    n_out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let chosen = voterank(
            g,
            if count < 0 {
                None
            } else {
                Some(count as usize)
            },
        );
        for (o, &v) in slice_mut(out, chosen.len()).iter_mut().zip(&chosen) {
            *o = v as i64;
        }
        slice_mut(n_out, 1)[0] = chosen.len() as i64;
        Ok(())
    })
}

/// Laplacian centrality (undirected reading); `out` receives n values.
///
/// # Safety
/// `h` must be a live handle; `out` n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_laplacian_centrality(
    h: *const GraphHandle,
    weighted: i32,
    normalized: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let v = laplacian_centrality(g, weighted != 0, normalized != 0)?;
        slice_mut(out, g.node_count()).copy_from_slice(&v);
        Ok(())
    })
}

/// Trophic measures of a directed network: `what` 0 = levels (n values),
/// 1 = differences (m values, one per edge), 2 = incoherence parameter (one
/// value; `cannibalism != 0` keeps self-loops).
///
/// # Safety
/// `h` must be a live handle; `out` sized as described.
#[no_mangle]
pub unsafe extern "C" fn net_trophic(
    h: *const GraphHandle,
    what: i32,
    weighted: i32,
    cannibalism: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let w = weighted != 0;
        match what {
            0 => slice_mut(out, g.node_count()).copy_from_slice(&trophic_levels(g, w)?),
            1 => slice_mut(out, g.edge_count()).copy_from_slice(&trophic_differences(g, w)?),
            2 => slice_mut(out, 1)[0] = trophic_incoherence(g, w, cannibalism != 0)?,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown trophic measure code {what}"
                )))
            }
        }
        Ok(())
    })
}

/// Local reaching centrality (`local`, n values; may be null) and global
/// reaching centrality (`global`, one value; may be null).
///
/// # Safety
/// `h` must be a live handle; `local` null or n doubles; `global` null or
/// one double.
#[no_mangle]
pub unsafe extern "C" fn net_reaching_centrality(
    h: *const GraphHandle,
    weighted: i32,
    normalized: i32,
    threads: i32,
    local: *mut f64,
    global: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = ReachingOptions {
            weighted: weighted != 0,
            normalized: normalized != 0,
            threads: thread_count(threads),
        };
        if !local.is_null() {
            slice_mut(local, g.node_count()).copy_from_slice(&local_reaching_centrality(g, opts)?);
        }
        if !global.is_null() {
            slice_mut(global, 1)[0] = global_reaching_centrality(g, opts)?;
        }
        Ok(())
    })
}
