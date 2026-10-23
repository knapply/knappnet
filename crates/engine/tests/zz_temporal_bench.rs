extern crate knappnet as engine;
use engine::rng::Stream;
use engine::temporal::*;
use std::time::Instant;

#[test]
#[ignore]
fn bench() {
    let (n, m) = (50_000u64, 5_000_000usize);
    let mut r = Stream::new(42, 0);
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    let mut st = Vec::with_capacity(m);
    for _ in 0..m {
        src.push(r.below(n) as i64);
        tgt.push(r.below(n) as i64);
        st.push(r.below(30 * 86400) as f64);
    }
    let t0 = Instant::now();
    let net = TemporalNetwork::new(n as usize, &src, &tgt, &st, None, true).unwrap();
    println!("construct {:?}", t0.elapsed());
    let t0 = Instant::now();
    let (a, b) = regular_windows(0.0, 30.0 * 86400.0, 86400.0, 3600.0 * 6.0).unwrap();
    let w = contacts_in_windows(&net, &a, &b).unwrap();
    println!(
        "snapshots {} windows, {} memberships {:?}",
        a.len(),
        w.indices.len(),
        t0.elapsed()
    );
    let t0 = Instant::now();
    let agg = aggregate(&net, None).unwrap();
    println!("aggregate {} pairs {:?}", agg.from.len(), t0.elapsed());
    let roots: Vec<usize> = (0..100).map(|i| i * 500).collect();
    for threads in [1, 0] {
        let t0 = Instant::now();
        let o = PathOptions {
            threads,
            ..Default::default()
        };
        let x =
            temporal_paths(&net, Criterion::EarliestArrival, Direction::Out, &roots, &o).unwrap();
        let reached = x.iter().filter(|v| v.is_finite()).count();
        println!(
            "EA 100 sources threads={threads}: {:?} (reached {reached})",
            t0.elapsed()
        );
    }
    let t0 = Instant::now();
    let o = PathOptions {
        threads: 0,
        max_wait: 3600.0,
        ..Default::default()
    };
    let _ = temporal_paths(&net, Criterion::Fastest, Direction::Out, &roots, &o).unwrap();
    println!("fastest max_wait 100 sources: {:?}", t0.elapsed());
    let t0 = Instant::now();
    let _ = burstiness(
        &net,
        EventScope::Node(engine::metrics::degree::DegreeMode::All),
    )
    .unwrap();
    println!("burstiness: {:?}", t0.elapsed());
}
