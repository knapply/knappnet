//! Minimal deterministic parallelism built on `std::thread::scope`.
//!
//! Work items are assigned to threads in a fixed strided pattern and partial
//! results are merged in thread order. For a given thread count the result
//! is therefore bit-for-bit reproducible. Floating-point results can differ
//! in the last few bits across different thread counts because addition is
//! not associative; integer results never do.

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
