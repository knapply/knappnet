//! Minimal deterministic parallelism built on `std::thread::scope`.
//!
//! [`for_each_item`] assigns work items to threads in a fixed strided
//! pattern and returns the partial results in thread order. For a given
//! thread count the result is therefore bit-for-bit reproducible.
//! Floating-point sums of per-thread partials can differ in the last few
//! bits across different thread counts because addition is not
//! associative; integer results never do.
//!
//! [`sum_in_order`] removes that dependence for sums over many items: items
//! are split into a fixed number of contiguous groups that depends only on
//! the item count, each group is summed in item order into its own
//! accumulator, and the group sums are added to the total in group order,
//! so the result is the same for every thread count. Threads claim groups
//! dynamically, which also balances items of very different cost (the
//! hubs of a scale-free network).

/// Resolves a requested thread count. Zero means "use all available cores".
pub fn resolve_threads(requested: usize, items: usize) -> usize {
    // WebAssembly targets without thread support cannot spawn threads, so
    // all work runs on the calling thread there.
    if cfg!(target_family = "wasm") {
        return 1;
    }
    let available = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let t = if requested == 0 { available } else { requested };
    t.clamp(1, items.max(1))
}

/// Runs `work(state, item)` for every item in `0..items`, where each thread
/// owns one `state` created by `init`. Returns the per-thread states in
/// thread order so the caller can merge them deterministically.
pub fn for_each_item<S, I, W>(items: usize, threads: usize, init: I, work: W) -> Vec<S>
where
    S: Send,
    I: Fn() -> S + Sync,
    W: Fn(&mut S, usize) + Sync,
{
    let threads = resolve_threads(threads, items);
    if threads == 1 {
        let mut state = init();
        for item in 0..items {
            work(&mut state, item);
        }
        return vec![state];
    }
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let init = &init;
                let work = &work;
                scope.spawn(move || {
                    let mut state = init();
                    let mut item = t;
                    while item < items {
                        work(&mut state, item);
                        item += threads;
                    }
                    state
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("worker thread panicked"))
            .collect()
    })
}

/// Number of groups [`sum_in_order`] splits its items into (fewer when
/// there are fewer items).
pub const SUM_GROUPS: usize = 64;

/// Runs `work(state, acc, item)` for every item in `0..items` and calls
/// `merge(acc)` once per group of items, in group order, so the caller can
/// add each group's accumulator to its total and reset it to empty.
///
/// The items are split into at most [`SUM_GROUPS`] contiguous groups whose
/// boundaries depend only on `items`. Each group is processed by one
/// thread in increasing item order into an accumulator that is empty when
/// the group starts (a new one from `fresh`, or one that `merge` has
/// reset). Each thread owns one `state` created by `init`. The sequence of
/// operations on every accumulator, and the order of the merges, are
/// therefore the same for every thread count, and so is the result.
///
/// At most about two accumulators per thread are alive at once: a thread
/// that gets more than `2 * threads` groups ahead of the next group to be
/// merged waits for it.
pub fn sum_in_order<S, A, I, F, W, M>(
    items: usize,
    threads: usize,
    init: I,
    fresh: F,
    work: W,
    merge: M,
) where
    S: Send,
    A: Send,
    I: Fn() -> S + Sync,
    F: Fn() -> A + Sync,
    W: Fn(&mut S, &mut A, usize) + Sync,
    M: FnMut(&mut A) + Send,
{
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Condvar, Mutex};

    if items == 0 {
        return;
    }
    let size = (items + SUM_GROUPS - 1) / SUM_GROUPS;
    let groups = (items + size - 1) / size;
    let range = |g: usize| g * size..((g + 1) * size).min(items);
    let threads = resolve_threads(threads, groups);
    let mut merge = merge;
    if threads == 1 {
        let mut state = init();
        let mut acc = fresh();
        for g in 0..groups {
            for item in range(g) {
                work(&mut state, &mut acc, item);
            }
            merge(&mut acc);
        }
        return;
    }

    struct Ordered<A, M> {
        /// The next group to merge.
        next: usize,
        /// Finished groups waiting for their turn.
        done: Vec<Option<A>>,
        /// Accumulators already merged and reset, for reuse.
        pool: Vec<A>,
        merge: M,
    }
    let claim = AtomicUsize::new(0);
    let ordered = Mutex::new(Ordered {
        next: 0,
        done: (0..groups).map(|_| None).collect(),
        pool: Vec::new(),
        merge,
    });
    let turn = Condvar::new();
    let window = 2 * threads;
    let lock = || ordered.lock().unwrap_or_else(|e| e.into_inner());
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let (init, fresh, work) = (&init, &fresh, &work);
                let (claim, turn, lock) = (&claim, &turn, &lock);
                scope.spawn(move || {
                    let mut state = init();
                    loop {
                        let g = claim.fetch_add(1, Ordering::Relaxed);
                        if g >= groups {
                            break;
                        }
                        let mut acc = {
                            let mut o = lock();
                            while g >= o.next + window {
                                o = turn.wait(o).unwrap_or_else(|e| e.into_inner());
                            }
                            o.pool.pop()
                        }
                        .unwrap_or_else(fresh);
                        for item in range(g) {
                            work(&mut state, &mut acc, item);
                        }
                        let mut o = lock();
                        o.done[g] = Some(acc);
                        let mut merged = false;
                        while o.next < groups {
                            let k = o.next;
                            let Some(mut a) = o.done[k].take() else {
                                break;
                            };
                            (o.merge)(&mut a);
                            o.pool.push(a);
                            o.next += 1;
                            merged = true;
                        }
                        drop(o);
                        if merged {
                            turn.notify_all();
                        }
                    }
                })
            })
            .collect();
        for h in handles {
            h.join().expect("worker thread panicked");
        }
    });
}

/// Coarse phase timings for profiling expensive measures.
///
/// Off by default, when each probe costs one relaxed atomic load. After
/// [`probe::enable`]`(true)`, instrumented measures record how long each of
/// their phases took on the calling thread (for example setup, traversal,
/// accumulation, and copying the result), and [`probe::take`] returns and
/// clears the record. The benchmark example uses this to show where the
/// time of a measure goes; results never depend on it.
pub mod probe {
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::Mutex;
    use std::time::Instant;

    static ENABLED: AtomicBool = AtomicBool::new(false);
    static RECORD: Mutex<Vec<(&'static str, f64)>> = Mutex::new(Vec::new());

    /// Turns phase recording on or off.
    pub fn enable(on: bool) {
        ENABLED.store(on, Ordering::Relaxed);
    }

    /// Returns the phases recorded since the last call, as (name, seconds)
    /// in the order they ended, and clears the record.
    pub fn take() -> Vec<(&'static str, f64)> {
        std::mem::take(&mut *RECORD.lock().unwrap_or_else(|e| e.into_inner()))
    }

    /// A running phase clock; `None` inside when recording is off.
    #[derive(Debug)]
    pub struct Clock(Option<Instant>);

    impl Clock {
        /// Starts the clock if recording is on.
        #[inline]
        pub fn start() -> Clock {
            Clock(if ENABLED.load(Ordering::Relaxed) {
                Some(Instant::now())
            } else {
                None
            })
        }

        /// Records the time since the last mark (or the start) under
        /// `name`, and restarts the clock.
        #[inline]
        pub fn mark(&mut self, name: &'static str) {
            if let Some(t) = self.0 {
                let now = Instant::now();
                RECORD
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .push((name, (now - t).as_secs_f64()));
                self.0 = Some(now);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sums of values whose floating-point sum depends on the order.
    fn ordered_sum(items: usize, threads: usize) -> (u64, Vec<usize>) {
        let mut total = 0.0f64;
        let mut seen = Vec::new();
        sum_in_order(
            items,
            threads,
            || (),
            || (0.0f64, Vec::new()),
            |_, acc: &mut (f64, Vec<usize>), i| {
                acc.0 += 1.0 / (1.0 + i as f64) * if i % 3 == 0 { 1e8 } else { 1.0 };
                acc.1.push(i);
            },
            |acc| {
                total += acc.0;
                seen.append(&mut acc.1);
                acc.0 = 0.0;
            },
        );
        (total.to_bits(), seen)
    }

    #[test]
    fn sums_in_order_for_every_thread_count() {
        for items in [0usize, 1, 5, 63, 64, 65, 1000, 4097] {
            let (one, seen) = ordered_sum(items, 1);
            assert_eq!(seen, (0..items).collect::<Vec<_>>());
            for t in [2, 3, 8, 100] {
                assert_eq!(
                    ordered_sum(items, t),
                    (one, seen.clone()),
                    "{items} items, {t} threads"
                );
            }
        }
    }
}
