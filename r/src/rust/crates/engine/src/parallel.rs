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
//!
//! Threads are started for each parallel call and joined at its end, which
//! costs some microseconds per thread. A call whose total work is small
//! therefore runs faster on fewer threads, or on one. [`threads_for`] turns
//! an estimate of the work into a thread count: at most one thread per
//! [`GRAIN`] elementary steps (an arc scanned, a multiply-add, a
//! comparison), so every thread has enough to do to repay its start. The
//! kernels in this crate pass their thread count through it. Every result
//! that does not depend on the thread count (all of those documented as
//! such) is unchanged by it.

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

/// Smallest number of elementary steps worth a thread of its own: about
/// 50 to 150 microseconds of work, against some microseconds to start and
/// join a thread.
pub const GRAIN: u64 = 1 << 17;

/// Threads to use for `items` independent items that take about `work`
/// elementary steps in all: the requested count (0 for all cores), at most
/// one per item, and at most one per [`GRAIN`] steps of work, so small
/// problems run on the calling thread alone.
pub fn threads_for(requested: usize, items: usize, work: u64) -> usize {
    let by_work = (work / GRAIN).clamp(1, usize::MAX as u64) as usize;
    resolve_threads(requested, items).min(by_work)
}

/// Work of `searches` graph searches over a network with `n` nodes and
/// `arcs` arcs, each of which may visit every node and arc once.
pub fn search_work(searches: usize, n: usize, arcs: usize) -> u64 {
    (searches as u64).saturating_mul(n as u64 + arcs as u64)
}

/// Steps of one search stopped after `levels` levels in a network with `n`
/// nodes and `arcs` arcs: with average degree d it reaches about d^levels
/// nodes and scans their arcs, and never more than the whole network.
pub fn ball_steps(n: usize, arcs: usize, levels: usize) -> u64 {
    let d = arcs as f64 / n.max(1) as f64;
    let ball = d.max(1.0).powf(levels.min(64) as f64) * (d + 1.0);
    ball.min((n + arcs) as f64) as u64
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

/// Runs `work(state, item)` for every item in `0..items`, where each thread
/// owns one `state` created by `init`, and returns the per-thread states.
///
/// Unlike [`for_each_item`], threads claim items dynamically, which
/// balances items of very different cost, but which thread processes which
/// item varies from run to run. A claim takes a strided chunk (items c,
/// c + k, c + 2k, ... for k chunks), so costly items that are numbered
/// together (the hubs of a preferential-attachment network come first) are
/// spread over many claims.
/// Use it only when the merged result cannot depend on that: integer sums,
/// maxima, results keyed by item, collections sorted afterwards.
pub fn for_each_claimed<S, I, W>(items: usize, threads: usize, init: I, work: W) -> Vec<S>
where
    S: Send,
    I: Fn() -> S + Sync,
    W: Fn(&mut S, usize) + Sync,
{
    use std::sync::atomic::{AtomicUsize, Ordering};
    let threads = resolve_threads(threads, items);
    if threads == 1 {
        let mut state = init();
        for item in 0..items {
            work(&mut state, item);
        }
        return vec![state];
    }
    let chunks = items.min(threads * 64);
    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let (init, work, next) = (&init, &work, &next);
                scope.spawn(move || {
                    let mut state = init();
                    loop {
                        let c = next.fetch_add(1, Ordering::Relaxed);
                        if c >= chunks {
                            break;
                        }
                        for item in (c..items).step_by(chunks) {
                            work(&mut state, item);
                        }
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

/// Returns `f(state, item)` for every item in `0..items`, in item order,
/// where each thread owns one `state` created by `init`. Items are claimed
/// dynamically as in [`for_each_claimed`]; the output is the same for every
/// thread count as long as each result depends only on its item.
pub fn map_items<S, T, I, F>(items: usize, threads: usize, init: I, f: F) -> Vec<T>
where
    T: Send,
    S: Send,
    I: Fn() -> S + Sync,
    F: Fn(&mut S, usize) -> T + Sync,
{
    let parts = for_each_claimed(
        items,
        threads,
        || (init(), Vec::new()),
        |(state, out): &mut (S, Vec<(usize, T)>), item| out.push((item, f(state, item))),
    );
    let mut slots: Vec<Option<T>> = (0..items).map(|_| None).collect();
    for (_, out) in parts {
        for (item, t) in out {
            slots[item] = Some(t);
        }
    }
    slots
        .into_iter()
        .map(|t| t.expect("every item is computed"))
        .collect()
}

/// Calls `f(state, c, chunk)` for every chunk `c` of `chunk` consecutive
/// elements of `out` (the last may be shorter), where each thread owns one
/// `state` created by `init`. Threads claim chunks dynamically and write
/// their results in place, so a large output is filled (and its memory
/// first touched) in parallel without intermediate copies; each element is
/// written by exactly one call, so the result does not depend on the
/// thread count as long as each chunk's values depend only on the chunk.
pub fn fill_chunks<T, S, I, F>(out: &mut [T], chunk: usize, threads: usize, init: I, f: F)
where
    T: Send,
    S: Send,
    I: Fn() -> S + Sync,
    F: Fn(&mut S, usize, &mut [T]) + Sync,
{
    let chunk = chunk.max(1);
    let count = (out.len() + chunk - 1) / chunk;
    let threads = resolve_threads(threads, count);
    if threads <= 1 {
        let mut state = init();
        for (c, part) in out.chunks_mut(chunk).enumerate() {
            f(&mut state, c, part);
        }
        return;
    }
    let queue = std::sync::Mutex::new(out.chunks_mut(chunk).enumerate());
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let (queue, init, f) = (&queue, &init, &f);
            scope.spawn(move || {
                let mut state = init();
                loop {
                    let next = queue.lock().unwrap_or_else(|e| e.into_inner()).next();
                    match next {
                        Some((c, part)) => f(&mut state, c, part),
                        None => break,
                    }
                }
            });
        }
    });
}

/// One member of a [`crew`]: its index, the crew's size, and the barrier
/// all members wait at between phases.
pub struct Crew<'a> {
    index: usize,
    size: usize,
    barrier: Option<&'a SpinBarrier>,
}

impl Crew<'_> {
    /// This member's index, from 0 to `size() - 1`.
    #[inline]
    pub fn index(&self) -> usize {
        self.index
    }

    /// Number of members.
    #[inline]
    pub fn size(&self) -> usize {
        self.size
    }

    /// This member's share of `items` items: one of `size()` contiguous
    /// ranges of nearly equal length, in member order.
    #[inline]
    pub fn share(&self, items: usize) -> std::ops::Range<usize> {
        (items * self.index / self.size)..(items * (self.index + 1) / self.size)
    }

    /// Waits until every member has reached this point. Everything a member
    /// wrote before waiting is visible to every member after it.
    pub fn wait(&self) {
        if let Some(b) = self.barrier {
            b.wait();
        }
    }
}

/// Barrier for [`crew`]: members spin briefly, then yield, so phases of a
/// few microseconds do not pay for putting threads to sleep and waking them.
pub struct SpinBarrier {
    size: usize,
    arrived: std::sync::atomic::AtomicUsize,
    generation: std::sync::atomic::AtomicUsize,
    /// Set when a member panicked, so the others stop waiting for it.
    poisoned: std::sync::atomic::AtomicBool,
}

impl SpinBarrier {
    fn new(size: usize) -> SpinBarrier {
        SpinBarrier {
            size,
            arrived: std::sync::atomic::AtomicUsize::new(0),
            generation: std::sync::atomic::AtomicUsize::new(0),
            poisoned: std::sync::atomic::AtomicBool::new(false),
        }
    }

    fn wait(&self) {
        use std::sync::atomic::Ordering;
        let gen = self.generation.load(Ordering::Acquire);
        if self.arrived.fetch_add(1, Ordering::AcqRel) + 1 == self.size {
            self.arrived.store(0, Ordering::Relaxed);
            self.generation
                .store(gen.wrapping_add(1), Ordering::Release);
            return;
        }
        let mut spins = 0u32;
        while self.generation.load(Ordering::Acquire) == gen {
            if self.poisoned.load(Ordering::Relaxed) {
                panic!("another member of the crew panicked");
            }
            if spins < 1 << 12 {
                spins += 1;
                std::hint::spin_loop();
            } else {
                std::thread::yield_now();
            }
        }
    }
}

/// Runs `body` on `threads` threads at once (0 uses all available cores),
/// each with its own [`Crew`] member, and returns when all have finished.
///
/// The members run the same code on different parts of the data and meet
/// at [`Crew::wait`] between phases, so an iterative method starts its
/// threads once instead of once per phase. The caller decides the
/// partition of work (for example with [`Crew::share`]); results that do not
/// depend on that partition do not depend on the thread count. With one
/// thread `body` runs on the calling thread and waiting costs nothing.
pub fn crew<F>(threads: usize, body: F)
where
    F: Fn(&Crew) + Sync,
{
    let size = resolve_threads(threads, usize::MAX);
    if size == 1 {
        body(&Crew {
            index: 0,
            size: 1,
            barrier: None,
        });
        return;
    }
    let barrier = SpinBarrier::new(size);
    /// Marks the barrier poisoned if its member unwinds.
    struct Guard<'a>(&'a SpinBarrier);
    impl Drop for Guard<'_> {
        fn drop(&mut self) {
            if std::thread::panicking() {
                self.0
                    .poisoned
                    .store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }
    }
    std::thread::scope(|scope| {
        let handles: Vec<_> = (1..size)
            .map(|index| {
                let (body, barrier) = (&body, &barrier);
                scope.spawn(move || {
                    let _guard = Guard(barrier);
                    body(&Crew {
                        index,
                        size,
                        barrier: Some(barrier),
                    });
                })
            })
            .collect();
        {
            let _guard = Guard(&barrier);
            body(&Crew {
                index: 0,
                size,
                barrier: Some(&barrier),
            });
        }
        for h in handles {
            h.join().expect("worker thread panicked");
        }
    });
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
    fn claimed_items_are_each_done_once() {
        for items in [0usize, 1, 7, 1000, 5000] {
            for t in [1, 3, 16] {
                let parts = for_each_claimed(items, t, Vec::new, |v: &mut Vec<usize>, i| v.push(i));
                let mut all: Vec<usize> = parts.into_iter().flatten().collect();
                all.sort_unstable();
                assert_eq!(all, (0..items).collect::<Vec<_>>());
                let squares = map_items(items, t, || (), |_, i| i * i);
                assert_eq!(squares, (0..items).map(|i| i * i).collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn crews_meet_at_barriers() {
        use std::sync::atomic::{AtomicU64, Ordering};
        for t in [1usize, 2, 5] {
            let cells: Vec<AtomicU64> = (0..100).map(|_| AtomicU64::new(0)).collect();
            crew(t, |c| {
                // Each phase reads what every member wrote in the last one.
                for phase in 1..=50u64 {
                    for i in c.share(cells.len()) {
                        assert_eq!(cells[(i + 99) % 100].load(Ordering::Relaxed), phase - 1);
                    }
                    c.wait();
                    for i in c.share(cells.len()) {
                        cells[i].store(phase, Ordering::Relaxed);
                    }
                    c.wait();
                }
            });
            assert!(cells.iter().all(|x| x.load(Ordering::Relaxed) == 50));
        }
    }

    #[test]
    fn chunks_are_filled_once() {
        for t in [1usize, 4] {
            let mut v = vec![0usize; 1003];
            fill_chunks(
                &mut v,
                10,
                t,
                || (),
                |_, c, part| {
                    for (i, x) in part.iter_mut().enumerate() {
                        *x += c * 10 + i;
                    }
                },
            );
            assert_eq!(v, (0..1003).collect::<Vec<_>>());
        }
    }

    #[test]
    fn threads_follow_the_work() {
        assert_eq!(threads_for(8, 100, 10), 1);
        assert_eq!(threads_for(8, 100, 3 * GRAIN), 3);
        assert_eq!(threads_for(8, 2, 100 * GRAIN), 2);
        assert_eq!(threads_for(8, 100, u64::MAX), 8);
        assert_eq!(threads_for(1, 100, u64::MAX), 1);
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
