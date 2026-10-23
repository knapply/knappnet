//! C ABI for key players, group centrality, and flow betweenness.
//!
//! Integer codes used by these functions:
//! - `measure` (key players): 0 = fragmentation, 1 = distance-weighted
//!   fragmentation, 2 = m-reach, 3 = distance-weighted reach.
//! - `measure` (group centrality): 0 = degree, 1 = closeness,
//!   2 = betweenness.
//! - `mode`: 1 = out, 2 = in, 3 = all (as for `net_degree`).
//! - `weights_as`: 0 = distance, 1 = strength.
//! - `method`: 0 = greedy, 1 = exact.
//! - `scale` (flow betweenness): 0 = raw, 1 = normalized, 2 = fractional.
//!
//! Node sets are 0-based indices without repeats.

use std::ffi::c_char;

use engine::metrics::degree::DegreeMode;
use engine::metrics::resilience::{
    flow_betweenness, group_centrality, key_players, kpp_score, FlowScale, GroupMeasure,
    GroupOptions, KppMeasure, KppMethod, KppOptions,
};
use engine::{Error, Graph, WeightKind};

use crate::{guard, handle, slice, slice_mut, GraphHandle};

fn node_set(g: &Graph, nodes: *const i64, k: i64) -> Result<Vec<u32>, Error> {
    if k < 0 {
        return Err(Error::InvalidArgument("negative set size".into()));
    }
    // SAFETY: the caller guarantees `nodes` holds `k` values.
    let raw = unsafe { slice(nodes, k as usize) };
    raw.iter()
        .map(|&v| {
            if v < 0 || v as u64 >= g.node_count() as u64 {
                Err(Error::InvalidArgument(format!(
                    "node index {v} is out of range for a network with {} nodes",
                    g.node_count()
                )))
            } else {
                Ok(v as u32)
            }
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn kpp_options(
    measure: i32,
    m: f64,
    mode: i32,
    weighted: i32,
    weights_as: i32,
) -> Result<KppOptions, Error> {
    Ok(KppOptions {
        measure: KppMeasure::from_code(measure)?,
        m,
        mode: DegreeMode::from_code(mode)?,
        weighted: weighted != 0,
        weights_as: WeightKind::from_code(weights_as)?,
        ..Default::default()
    })
}

/// Value of a key player criterion for the set `set` of `k` nodes. `m` is
/// the distance limit for m-reach and is ignored otherwise.
///
/// # Safety
/// `h` must be a live handle; `set` must hold `k` values; `out` one double.
#[no_mangle]
pub unsafe extern "C" fn net_kpp_score(
    h: *const GraphHandle,
    set: *const i64,
    k: i64,
    measure: i32,
    m: f64,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let set = node_set(g, set, k)?;
        let opts = kpp_options(measure, m, mode, weighted, weights_as)?;
        slice_mut(out, 1)[0] = kpp_score(g, &set, opts)?;
        Ok(())
    })
}

/// Finds `k` key players by greedy search (`restarts` random starts, start
/// r drawn from stream r of `seed`) or exact enumeration. `set` receives
/// the members in increasing order, `score` the criterion's value, and
/// `evaluations` the number of sets evaluated.
///
/// # Safety
/// `h` must be a live handle; `set` must hold `k` values; `score` and
/// `evaluations` one value each.
#[no_mangle]
pub unsafe extern "C" fn net_key_players(
    h: *const GraphHandle,
    k: i64,
    measure: i32,
    m: f64,
    mode: i32,
    weighted: i32,
    weights_as: i32,
    method: i32,
    restarts: i64,
    seed: u64,
    threads: i32,
    set: *mut i64,
    score: *mut f64,
    evaluations: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        if k < 1 {
            return Err(Error::InvalidArgument(
                "the number of key players must be at least 1".into(),
            ));
        }
        let opts = KppOptions {
            method: KppMethod::from_code(method)?,
            restarts: restarts.max(1) as usize,
            seed,
            threads: threads.max(0) as usize,
            ..kpp_options(measure, m, mode, weighted, weights_as)?
        };
        let found = key_players(g, k as usize, opts)?;
        for (d, &v) in slice_mut(set, k as usize).iter_mut().zip(&found.set) {
            *d = v as i64;
        }
        slice_mut(score, 1)[0] = found.score;
        slice_mut(evaluations, 1)[0] = found.evaluations as f64;
        Ok(())
    })
}

/// Group degree, closeness, or betweenness of the `k` nodes in `group`.
///
/// # Safety
/// `h` must be a live handle; `group` must hold `k` values; `out` one
/// double.
#[no_mangle]
pub unsafe extern "C" fn net_group_centrality(
    h: *const GraphHandle,
    group: *const i64,
    k: i64,
    measure: i32,
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
        let group = node_set(g, group, k)?;
        let opts = GroupOptions {
            mode: DegreeMode::from_code(mode)?,
            weighted: weighted != 0,
            weights_as: WeightKind::from_code(weights_as)?,
            normalized: normalized != 0,
            threads: threads.max(0) as usize,
        };
        slice_mut(out, 1)[0] =
            group_centrality(g, &group, GroupMeasure::from_code(measure)?, opts)?;
        Ok(())
    })
}

/// Flow betweenness of every node. With `capacity != 0` edge weights are
/// capacities; otherwise every edge has capacity 1.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_flow_betweenness(
    h: *const GraphHandle,
    capacity: i32,
    scale: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let values = flow_betweenness(
            g,
            capacity != 0,
            FlowScale::from_code(scale)?,
            threads.max(0) as usize,
        )?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    })
}
