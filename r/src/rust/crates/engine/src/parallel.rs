//! Minimal deterministic parallelism built on `std::thread::scope`.
//!
//! Work items are assigned to threads in a fixed strided pattern and partial
//! results are merged in thread order. For a given thread count the result
//! is therefore bit-for-bit reproducible. Floating-point results can differ
//! in the last few bits across different thread counts because addition is
//! not associative; integer results never do.

/// Resolves a requested thread count. Zero means "use all available cores".
pub fn resolve_threads(requested: usize, items: usize) -> usize {
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
