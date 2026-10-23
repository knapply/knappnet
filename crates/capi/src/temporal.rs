//! C ABI for temporal networks: building a network of timed contacts,
//! time-respecting paths, reachability, temporal closeness, slicing into
//! windows, aggregation, and inter-event statistics.
//!
//! Integer codes used by these functions:
//! - `criterion`: 0 = earliest arrival, 1 = latest departure, 2 = fastest,
//!   3 = shortest (fewest contacts).
//! - `direction`: 1 = walks from the roots, 2 = walks into the roots.
//! - `traversal`: 0 = the constant `traversal_time` for every contact,
//!   1 = each contact's duration.
//! - `scope`: 0 = nodes, 1 = node pairs; `mode` (nodes of directed
//!   networks): 1 = out, 2 = in, 3 = all.
//!
//! Window bounds `start` and `end` may be NaN to use the earliest start
//! time and the latest arrival of any contact; `max_wait` may be infinite.
//!
//! Results whose length is not known in advance are returned as a
//! `TemporalResult` holding numbered columns of doubles, read with
//! `net_temporal_result_length` and `net_temporal_result_copy` and released
//! with `net_temporal_result_free`.

use std::ffi::c_char;

use engine::metrics::degree::DegreeMode;
use engine::temporal::{
    aggregate, burstiness, contacts_at, contacts_in_windows, inter_event_times, temporal_closeness,
    temporal_paths, temporal_reach_counts, Criterion, Direction, EventScope, PathOptions,
    TemporalNetwork, Traversal,
};
use engine::Error;

use crate::{deref, guard, slice, slice_mut, thread_count};

/// Opaque temporal network handle.
pub struct TemporalHandle(TemporalNetwork);

/// Opaque result: numbered columns of doubles.
pub struct TemporalResult(Vec<Vec<f64>>);

fn temporal<'a>(h: *const TemporalHandle) -> Result<&'a TemporalNetwork, Error> {
    // SAFETY: non-null handles are only created by `net_temporal_new`.
    unsafe { deref(h, "temporal network") }.map(|t| &t.0)
}

fn result<'a>(h: *const TemporalResult) -> Result<&'a [Vec<f64>], Error> {
    // SAFETY: non-null result handles are only created by this module.
    unsafe { deref(h, "result") }.map(|r| &r.0[..])
}

#[allow(clippy::too_many_arguments)]
fn options(
    start: f64,
    end: f64,
    min_wait: f64,
    max_wait: f64,
    strict: i32,
    traversal: i32,
    traversal_time: f64,
    threads: i32,
) -> Result<PathOptions, Error> {
    let traversal = match traversal {
        0 => Traversal::Constant(traversal_time),
        1 => Traversal::Duration,
        _ => {
            return Err(Error::InvalidArgument(format!(
                "unknown traversal code {traversal}"
            )))
        }
    };
    Ok(PathOptions {
        start: (!start.is_nan()).then_some(start),
        end: (!end.is_nan()).then_some(end),
        min_wait,
        max_wait,
        strict: strict != 0,
        traversal,
        threads: thread_count(threads),
    })
}

fn scope(scope: i32, mode: i32) -> Result<EventScope, Error> {
    match scope {
        0 => Ok(EventScope::Node(DegreeMode::from_code(mode)?)),
        1 => Ok(EventScope::Pair),
        _ => Err(Error::InvalidArgument(format!(
            "unknown event scope code {scope}"
        ))),
    }
}

fn indices(v: &[usize]) -> Vec<f64> {
    v.iter().map(|&i| i as f64).collect()
}

/// Builds a temporal network from `m` contacts with 0-based endpoints and
/// start times; `end` may be null for instantaneous contacts. On success
/// `*out` receives a handle that must be released with
/// `net_temporal_free`.
///
/// # Safety
/// `src`, `tgt`, `start` (and `end` unless null) must reference `m`
/// elements; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_new(
    n: i64,
    src: *const i64,
    tgt: *const i64,
    start: *const f64,
    end: *const f64,
    m: i64,
    directed: i32,
    out: *mut *mut TemporalHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        if n < 0 || m < 0 {
            return Err(Error::InvalidArgument(
                "negative node or contact count".into(),
            ));
        }
        let (n, m) = (n as usize, m as usize);
        let e = if end.is_null() {
            None
        } else {
            Some(slice(end, m))
        };
        let net = TemporalNetwork::new(
            n,
            slice(src, m),
            slice(tgt, m),
            slice(start, m),
            e,
            directed != 0,
        )?;
        *out = Box::into_raw(Box::new(TemporalHandle(net)));
        Ok(())
    })
}

/// Releases a handle created by `net_temporal_new`. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_free(h: *mut TemporalHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Node count of a temporal network, or -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_node_count(h: *const TemporalHandle) -> i64 {
    temporal(h).map(|t| t.node_count() as i64).unwrap_or(-1)
}

/// Optimal time-respecting walks between `roots` and every node:
/// `out[r * n + j]` describes walks from `roots[r]` to j (`direction` 1)
/// or from j to `roots[r]` (`direction` 2). Unreachable pairs get +inf
/// (-inf for latest departure).
///
/// # Safety
/// `h` must be a live handle; `roots` must hold `nroots` values and `out`
/// `nroots * n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_paths(
    h: *const TemporalHandle,
    criterion: i32,
    direction: i32,
    roots: *const i64,
    nroots: i64,
    start: f64,
    end: f64,
    min_wait: f64,
    max_wait: f64,
    strict: i32,
    traversal: i32,
    traversal_time: f64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let opts = options(
            start,
            end,
            min_wait,
            max_wait,
            strict,
            traversal,
            traversal_time,
            threads,
        )?;
        let k = nroots.max(0) as usize;
        let mut r = Vec::with_capacity(k);
        for &x in slice(roots, k) {
            if x < 0 || x as u64 >= net.node_count() as u64 {
                return Err(Error::InvalidArgument(format!(
                    "node index {x} is out of range; valid indices are 0..{}",
                    net.node_count()
                )));
            }
            r.push(x as usize);
        }
        let values = temporal_paths(
            net,
            Criterion::from_code(criterion)?,
            Direction::from_code(direction)?,
            &r,
            &opts,
        )?;
        slice_mut(out, k * net.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// For every node, the number of other nodes it reaches (`direction` 1)
/// or that reach it (`direction` 2) by time-respecting walks.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_reach(
    h: *const TemporalHandle,
    direction: i32,
    start: f64,
    end: f64,
    min_wait: f64,
    max_wait: f64,
    strict: i32,
    traversal: i32,
    traversal_time: f64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let opts = options(
            start,
            end,
            min_wait,
            max_wait,
            strict,
            traversal,
            traversal_time,
            threads,
        )?;
        let values = temporal_reach_counts(net, Direction::from_code(direction)?, &opts)?;
        slice_mut(out, net.node_count()).copy_from_slice(&values);
        Ok(())
    })
}

/// Pan and Saramäki's temporal closeness into `closeness` (n values) and,
/// unless `distances` is null, the average temporal distances (n * n
/// values, row = origin).
///
/// # Safety
/// `h` must be a live handle; `closeness` must hold `n` doubles and
/// `distances`, unless null, `n * n`.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_closeness(
    h: *const TemporalHandle,
    start: f64,
    end: f64,
    min_wait: f64,
    max_wait: f64,
    strict: i32,
    traversal: i32,
    traversal_time: f64,
    threads: i32,
    closeness: *mut f64,
    distances: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let opts = options(
            start,
            end,
            min_wait,
            max_wait,
            strict,
            traversal,
            traversal_time,
            threads,
        )?;
        let n = net.node_count();
        let c = temporal_closeness(net, &opts, !distances.is_null())?;
        slice_mut(closeness, n).copy_from_slice(&c.closeness);
        if !distances.is_null() {
            slice_mut(distances, n * n).copy_from_slice(&c.distances);
        }
        Ok(())
    })
}

/// Contacts active in each window [from[k], to[k]). Result columns:
/// 0 = offsets (k + 1 values), 1 = 0-based contact indices, window k
/// holding entries offsets[k] to offsets[k + 1] - 1.
///
/// # Safety
/// `h` must be a live handle; `from` and `to` must hold `k` values; `out`
/// must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_windows(
    h: *const TemporalHandle,
    from: *const f64,
    to: *const f64,
    k: i64,
    out: *mut *mut TemporalResult,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let k = k.max(0) as usize;
        let w = contacts_in_windows(net, slice(from, k), slice(to, k))?;
        *out = Box::into_raw(Box::new(TemporalResult(vec![
            indices(&w.offsets),
            indices(&w.indices),
        ])));
        Ok(())
    })
}

/// Contacts active at time `t`. Result column 0: 0-based contact indices.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_active_at(
    h: *const TemporalHandle,
    t: f64,
    out: *mut *mut TemporalResult,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let c = contacts_at(net, t)?;
        *out = Box::into_raw(Box::new(TemporalResult(vec![indices(&c)])));
        Ok(())
    })
}

/// Contacts grouped by node pair; with `windowed != 0` only contacts
/// active in [from, to) count and durations are clipped to it. Result
/// columns: 0 = from, 1 = to (0-based), 2 = count, 3 = duration,
/// 4 = first start, 5 = last end, 6 = pair of each contact (-1 when
/// outside the window).
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_aggregate(
    h: *const TemporalHandle,
    windowed: i32,
    from: f64,
    to: f64,
    out: *mut *mut TemporalResult,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let a = aggregate(net, (windowed != 0).then_some((from, to)))?;
        let pair_of = a
            .pair_of
            .iter()
            .map(|&p| if p == usize::MAX { -1.0 } else { p as f64 })
            .collect();
        *out = Box::into_raw(Box::new(TemporalResult(vec![
            a.from.iter().map(|&x| x as f64).collect(),
            a.to.iter().map(|&x| x as f64).collect(),
            a.count,
            a.duration,
            a.first,
            a.last,
            pair_of,
        ])));
        Ok(())
    })
}

/// Inter-event times of every node (`scope` 0) or node pair (`scope` 1).
/// Result columns: 0 = group (node, or pair as numbered by
/// `net_temporal_aggregate` without a window) of each inter-event time,
/// 1 = the inter-event times, 2 = number of events per group.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_inter_event(
    h: *const TemporalHandle,
    scope_code: i32,
    mode: i32,
    out: *mut *mut TemporalResult,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let ie = inter_event_times(net, scope(scope_code, mode)?)?;
        *out = Box::into_raw(Box::new(TemporalResult(vec![
            ie.group.iter().map(|&g| g as f64).collect(),
            ie.gap,
            ie.events,
        ])));
        Ok(())
    })
}

/// Goh and Barabási's burstiness and memory of every node or pair.
/// Result columns: 0 = number of events, 1 = burstiness B, 2 = memory M.
///
/// # Safety
/// `h` must be a live handle; `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_burstiness(
    h: *const TemporalHandle,
    scope_code: i32,
    mode: i32,
    out: *mut *mut TemporalResult,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = temporal(h)?;
        let b = burstiness(net, scope(scope_code, mode)?)?;
        *out = Box::into_raw(Box::new(TemporalResult(vec![
            b.events,
            b.burstiness,
            b.memory,
        ])));
        Ok(())
    })
}

/// Length of result column `col`, or -1 for a null handle or a column
/// that does not exist.
///
/// # Safety
/// `h` must be null or a live result handle.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_result_length(h: *const TemporalResult, col: i32) -> i64 {
    match result(h) {
        Ok(cols) if col >= 0 && (col as usize) < cols.len() => cols[col as usize].len() as i64,
        _ => -1,
    }
}

/// Copies result column `col` into `out`.
///
/// # Safety
/// `h` must be a live result handle; `out` must hold as many doubles as
/// `net_temporal_result_length` reports.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_result_copy(
    h: *const TemporalResult,
    col: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let cols = result(h)?;
        if col < 0 || col as usize >= cols.len() {
            return Err(Error::InvalidArgument(format!(
                "result column {col} does not exist"
            )));
        }
        let c = &cols[col as usize];
        slice_mut(out, c.len()).copy_from_slice(c);
        Ok(())
    })
}

/// Releases a result. Null is ignored.
///
/// # Safety
/// `h` must be null or a live result handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_temporal_result_free(h: *mut TemporalResult) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}
