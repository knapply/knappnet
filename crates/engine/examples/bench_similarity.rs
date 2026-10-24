//! Benchmarks the link prediction, hashing, similarity, isomorphism,
//! core, mixing, cut, matching, colouring, and rewiring kernels on random
//! networks, on one thread and on every core.
//!
//!     cargo run --release --offline --example bench_similarity
//!     BENCH=truss,wl REPS=5 cargo run --release --offline --example bench_similarity
//!
//! Environment variables: `BENCH` (comma-separated substrings of workload
//! names), `REPS` (runs per measurement, default 3; the median is
//! reported), `SCALE` (multiplies the sizes, default 1). Prints one line
//! per workload: size, time on one thread, time on every core, and the
//! peak of live heap bytes during the one-thread run. Allocations are
//! counted only during the one-thread runs, so the shared counter does
//! not slow the parallel ones.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

use knappnet::generate::{gnm, preferential};
use knappnet::metrics::assortativity::{average_neighbor_degree, degree_assortativity, DegreeKind};
use knappnet::metrics::coloring::{greedy_color, ColoringStrategy};
use knappnet::metrics::cores::{onion_layers, truss_numbers};
use knappnet::metrics::cut_measures::cut_measures;
use knappnet::metrics::isomorphism::{count_mappings, Matcher, Problem};
use knappnet::metrics::linkpred::{
    link_scores, link_scores_non_edges, non_edge_count, LinkMethod, LinkOptions,
};
use knappnet::metrics::matching::max_weight_matching;
use knappnet::metrics::rewiring::{double_edge_swap, random_reference};
use knappnet::metrics::similarity::{simrank, SimRankOptions};
use knappnet::metrics::squares::square_clustering;
use knappnet::metrics::wl::{wl_graph_hash, WlOptions};
use knappnet::rng::Stream;
use knappnet::Graph;

struct Counting;
static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static COUNT: AtomicBool = AtomicBool::new(false);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = System.alloc(l);
        if !p.is_null() && COUNT.load(Ordering::Relaxed) {
            let now = LIVE.fetch_add(l.size(), Ordering::Relaxed) + l.size();
            PEAK.fetch_max(now, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l);
        if COUNT.load(Ordering::Relaxed) {
            LIVE.fetch_sub(l.size(), Ordering::Relaxed);
        }
    }
}

#[global_allocator]
static ALLOC: Counting = Counting;

fn graph_of(e: &knappnet::generate::EdgeList, weights: bool) -> Graph {
    let s: Vec<i64> = e.from.iter().map(|&x| x as i64).collect();
    let t: Vec<i64> = e.to.iter().map(|&x| x as i64).collect();
    let w: Option<Vec<f64>> = weights.then(|| {
        let mut rng = Stream::new(9, 9);
        (0..s.len()).map(|_| 1.0 + rng.below(100) as f64).collect()
    });
    Graph::from_edges(e.n, &s, &t, w.as_deref(), e.directed).unwrap()
}

fn sparse(n: usize, deg: usize, directed: bool, weights: bool) -> Graph {
    let m = (n * deg / 2) as u64;
    graph_of(&gnm(n, m, directed, false, 1, 1, 0).unwrap()[0], weights)
}

fn scale_free(n: usize, m: usize) -> Graph {
    graph_of(&preferential(n, m, 1, 2, 0).unwrap()[0], false)
}

fn time(reps: usize, mut f: impl FnMut()) -> (f64, usize) {
    let mut times = Vec::new();
    let mut peak = 0;
    for _ in 0..reps {
        PEAK.store(LIVE.load(Ordering::Relaxed), Ordering::Relaxed);
        let base = LIVE.load(Ordering::Relaxed);
        let t0 = Instant::now();
        f();
        times.push(t0.elapsed().as_secs_f64());
        peak = peak.max(PEAK.load(Ordering::Relaxed).saturating_sub(base));
    }
    times.sort_by(|a, b| a.total_cmp(b));
    (times[times.len() / 2], peak)
}

fn main() {
    let filter: Vec<String> = std::env::var("BENCH")
        .map(|s| s.split(',').map(str::to_string).collect())
        .unwrap_or_default();
    let reps: usize = std::env::var("REPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let k: usize = std::env::var("SCALE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);
    let want = |name: &str| filter.is_empty() || filter.iter().any(|f| name.contains(f.as_str()));
    println!("| workload | size | 1 thread (s) | all cores (s) | peak heap (MB) |");
    println!("|---|---|---|---|---|");
    let run = |name: &str, size: String, f: &dyn Fn(usize)| {
        if !want(name) {
            return;
        }
        COUNT.store(true, Ordering::Relaxed);
        let (t1, peak) = time(reps, || f(1));
        COUNT.store(false, Ordering::Relaxed);
        let (tn, _) = time(reps, || f(0));
        println!(
            "| {name} | {size} | {t1:.4} | {tn:.4} | {:.1} |",
            peak as f64 / 1e6
        );
    };

    let g = sparse(100_000 * k, 10, false, false);
    let size = format!("n {} m {}", g.node_count(), g.edge_count());
    let mut rng = Stream::new(3, 3);
    let pairs: Vec<(u32, u32)> = (0..1_000_000)
        .map(|_| {
            (
                rng.below(g.node_count() as u64) as u32,
                rng.below(g.node_count() as u64) as u32,
            )
        })
        .collect();
    let (us, vs): (Vec<u32>, Vec<u32>) = pairs.iter().copied().unzip();
    run("linkpred_ra_pairs_1M", size.clone(), &|t| {
        let o = LinkOptions {
            method: LinkMethod::ResourceAllocation,
            threads: t,
            ..Default::default()
        };
        link_scores(&g, &us, &vs, None, o).unwrap();
    });
    let small = sparse(3_000 * k, 40, false, false);
    run(
        "linkpred_jaccard_all_non_edges",
        format!("n {} m {}", small.node_count(), small.edge_count()),
        &|t| {
            let o = LinkOptions {
                method: LinkMethod::Jaccard,
                threads: t,
                ..Default::default()
            };
            let c = non_edge_count(&small);
            let (mut a, mut b, mut s) = (vec![0; c], vec![0; c], vec![0.0; c]);
            link_scores_non_edges(&small, None, o, &mut a, &mut b, &mut s).unwrap();
        },
    );
    run("wl_graph_hash", size.clone(), &|t| {
        wl_graph_hash(
            &g,
            WlOptions {
                threads: t,
                ..Default::default()
            },
        )
        .unwrap();
    });
    let sim = sparse(2_000 * k, 10, true, false);
    run(
        "simrank",
        format!("n {} m {}", sim.node_count(), sim.edge_count()),
        &|t| {
            let _ = simrank(
                &sim,
                SimRankOptions {
                    threads: t,
                    ..Default::default()
                },
            );
        },
    );
    let pattern = Graph::from_edges(4, &[0, 1, 2, 3], &[1, 2, 3, 0], None, false).unwrap();
    run("vf2pp_mono_4cycles_cap1e6", size.clone(), &|_| {
        let mut m = Matcher::new(&g, &pattern, None, None, Problem::Monomorphism, false).unwrap();
        count_mappings(&mut m, 1_000_000);
    });
    let sf = scale_free(200_000 * k, 5);
    let sfs = format!("BA n {} m {}", sf.node_count(), sf.edge_count());
    run("truss_numbers", sfs.clone(), &|_| {
        truss_numbers(&sf);
    });
    run("onion_layers", sfs.clone(), &|_| {
        onion_layers(&sf);
    });
    run("square_clustering", sfs.clone(), &|t| {
        square_clustering(&sf, t).unwrap();
    });
    run("degree_assortativity", sfs.clone(), &|_| {
        degree_assortativity(&sf, DegreeKind::Out, DegreeKind::In, false, None).unwrap();
    });
    run("average_neighbor_degree", sfs.clone(), &|_| {
        average_neighbor_degree(&sf, DegreeKind::Total, DegreeKind::Total, false, None).unwrap();
    });
    let half: Vec<u32> = (0..sf.node_count() as u32 / 2).collect();
    run("cut_measures", sfs.clone(), &|_| {
        cut_measures(&sf, &half, None, false).unwrap();
    });
    run("greedy_color_dsatur", sfs.clone(), &|_| {
        greedy_color(&sf, ColoringStrategy::SaturationLargestFirst, 0).unwrap();
    });
    run("greedy_color_smallest_last", sfs.clone(), &|_| {
        greedy_color(&sf, ColoringStrategy::SmallestLast, 0).unwrap();
    });
    let wm = sparse(20_000 * k, 6, false, true);
    run(
        "max_weight_matching",
        format!("n {} m {}", wm.node_count(), wm.edge_count()),
        &|_| {
            max_weight_matching(&wm, true, false).unwrap();
        },
    );
    run("double_edge_swap_1M", size.clone(), &|_| {
        double_edge_swap(&g, 1_000_000, 100_000_000, 1).unwrap();
    });
    let rr = sparse(2_000 * k, 8, false, false);
    run(
        "random_reference_niter1",
        format!("n {} m {}", rr.node_count(), rr.edge_count()),
        &|_| {
            random_reference(&rr, 1, true, 1).unwrap();
        },
    );
}
