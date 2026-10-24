//! Benchmarks the flow, connectivity, tree, DAG, and cycle kernels: wall
//! time (median of `REPS` runs), peak live heap, and allocations, on
//! random networks of growing size, on one thread and on every core for
//! the parallel kernels.
//!
//!     cargo run --release --offline --example bench_flows
//!     BENCH=flow,cycle REPS=5 cargo run --release --offline --example bench_flows
//!
//! Environment variables: `BENCH` (comma-separated substrings of workload
//! names to run), `REPS` (runs per measurement, default 3), `MAX` (largest
//! size multiplier, default 1; 0 runs only the smallest sizes).
//! Results are printed as a Markdown table.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

use knappnet::flow::{gomory_hu_tree, max_flow_min_cost, maximum_flow, minimum_cut, CutSide};
use knappnet::metrics::connectivity::{
    all_pairs_node_connectivity, edge_connectivity, k_components, k_edge_components,
    node_connectivity,
};
use knappnet::metrics::cycles::{chordless_cycles, girth, minimum_cycle_basis, simple_cycles};
use knappnet::metrics::dag::{
    immediate_dominators, topological_sort, transitive_closure, transitive_reduction, Reflexive,
};
use knappnet::metrics::trees::{
    optimum_branching, spanning_tree, BranchingKind, SpanningAlgorithm,
};
use knappnet::order::NeighborOrder;
use knappnet::rng::Stream;
use knappnet::Graph;

struct Counting;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static ALLOCS: AtomicU64 = AtomicU64::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        if !p.is_null() {
            let now = LIVE.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            PEAK.fetch_max(now, Ordering::Relaxed);
            ALLOCS.fetch_add(1, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        System.dealloc(ptr, layout);
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
    }
}

#[global_allocator]
static GLOBAL: Counting = Counting;

/// A random network with `n` nodes and about `n * degree / 2` edges
/// (arcs if directed), capacities or weights 1..=10.
fn random(n: usize, degree: f64, directed: bool, seed: u64) -> Graph {
    let mut rng = Stream::new(seed, 0);
    let m = (n as f64 * degree / if directed { 1.0 } else { 2.0 }) as usize;
    let mut s = Vec::with_capacity(m);
    let mut t = Vec::with_capacity(m);
    let mut w = Vec::with_capacity(m);
    while s.len() < m {
        let a = rng.below(n as u64) as i64;
        let b = rng.below(n as u64) as i64;
        if a != b {
            s.push(a);
            t.push(b);
            w.push((1 + rng.below(10)) as f64);
        }
    }
    Graph::from_edges(n, &s, &t, Some(&w), directed).expect("valid")
}

/// A random DAG: arcs from lower to higher index.
fn random_dag(n: usize, degree: f64, seed: u64) -> Graph {
    let g = random(n, degree, true, seed);
    let (s, t): (Vec<i64>, Vec<i64>) = g
        .edge_sources()
        .iter()
        .zip(g.edge_targets())
        .map(|(&a, &b)| (a.min(b) as i64, a.max(b) as i64))
        .unzip();
    Graph::from_edges(n, &s, &t, g.edge_weights(), true).expect("valid")
}

struct Row {
    name: String,
    size: String,
    threads: String,
    seconds: f64,
    heap_mb: f64,
    allocs: u64,
}

fn measure(reps: usize, mut f: impl FnMut()) -> (f64, f64, u64) {
    let mut times = Vec::new();
    let mut heap = 0usize;
    let mut allocs = 0u64;
    for _ in 0..reps {
        let base = LIVE.load(Ordering::Relaxed);
        PEAK.store(base, Ordering::Relaxed);
        let a0 = ALLOCS.load(Ordering::Relaxed);
        let t0 = Instant::now();
        f();
        times.push(t0.elapsed().as_secs_f64());
        heap = heap.max(PEAK.load(Ordering::Relaxed) - base);
        allocs = ALLOCS.load(Ordering::Relaxed) - a0;
    }
    times.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    (times[times.len() / 2], heap as f64 / 1e6, allocs)
}

fn main() {
    let filter: Vec<String> = std::env::var("BENCH")
        .ok()
        .map(|s| s.split(',').map(|x| x.trim().to_string()).collect())
        .unwrap_or_default();
    let reps: usize = std::env::var("REPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let scale: usize = std::env::var("MAX")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);
    let wanted = |name: &str| filter.is_empty() || filter.iter().any(|f| name.contains(f));
    let cores = std::thread::available_parallelism()
        .map(|c| c.get())
        .unwrap_or(1);
    let mut rows: Vec<Row> = Vec::new();
    let mut run = |name: &str, size: String, threads: usize, f: &mut dyn FnMut()| {
        if !wanted(name) {
            return;
        }
        let (s, h, a) = measure(reps, f);
        eprintln!("{name} {size} t={threads}: {s:.4} s");
        rows.push(Row {
            name: name.into(),
            size,
            threads: if threads == 0 {
                format!("{cores}")
            } else {
                threads.to_string()
            },
            seconds: s,
            heap_mb: h,
            allocs: a,
        });
    };
    // Sizes, with a size listed twice (when MAX = 0 caps it) run once.
    let sizes = |list: &[usize]| -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        for &x in list {
            if !out.contains(&x) {
                out.push(x);
            }
        }
        out
    };
    let big = |small: usize, large: usize| if scale == 0 { small } else { large };

    // Flows.
    for n in sizes(&[1_000usize, 10_000, big(10_000, 100_000)]) {
        let g = random(n, 8.0, true, 1);
        let cap = g.edge_weights().expect("weights").to_vec();
        run(
            "flow_max_flow",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                maximum_flow(&g, &[0], &[n - 1], Some(&cap)).expect("flow");
            },
        );
        run(
            "flow_min_cut",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                minimum_cut(&g, &[0], &[n - 1], Some(&cap), CutSide::Largest).expect("cut");
            },
        );
    }
    for n in sizes(&[1_000usize, big(1_000, 10_000)]) {
        let g = random(n, 8.0, true, 2);
        let cap = g.edge_weights().expect("weights").to_vec();
        let cost: Vec<f64> = (0..g.edge_count()).map(|e| (e % 7) as f64).collect();
        run(
            "flow_min_cost",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                max_flow_min_cost(&g, 0, n - 1, Some(&cap), Some(&cost)).expect("cost");
            },
        );
    }
    for n in sizes(&[300usize, 1_000, big(1_000, 3_000)]) {
        let g = random(n, 8.0, false, 3);
        let cap = g.edge_weights().expect("weights").to_vec();
        run(
            "flow_gomory_hu",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                gomory_hu_tree(&g, Some(&cap)).expect("tree");
            },
        );
    }

    // Connectivity.
    for n in sizes(&[300usize, 1_000, big(1_000, 3_000)]) {
        let g = random(n, 20.0, false, 4);
        for t in [1, 0] {
            run(
                "conn_node",
                format!("n={n} m={}", g.edge_count()),
                t,
                &mut || {
                    node_connectivity(&g, t).expect("kappa");
                },
            );
            run(
                "conn_edge",
                format!("n={n} m={}", g.edge_count()),
                t,
                &mut || {
                    edge_connectivity(&g, t).expect("lambda");
                },
            );
        }
    }
    for n in sizes(&[100usize, 300]) {
        let g = random(n, 8.0, false, 5);
        for t in [1, 0] {
            run(
                "conn_all_pairs_node",
                format!("n={n} m={}", g.edge_count()),
                t,
                &mut || {
                    all_pairs_node_connectivity(&g, t);
                },
            );
        }
        run(
            "conn_k_components",
            format!("n={n} m={}", g.edge_count()),
            0,
            &mut || {
                k_components(&g, 0).expect("k-components");
            },
        );
    }
    for n in sizes(&[1_000usize, big(1_000, 5_000)]) {
        let g = random(n, 8.0, false, 6);
        run(
            "conn_k_edge_components",
            format!("n={n} k=3"),
            0,
            &mut || {
                k_edge_components(&g, 3, 0).expect("classes");
            },
        );
    }

    // Trees and DAGs.
    for n in sizes(&[10_000usize, big(10_000, 1_000_000)]) {
        let g = random(n, 8.0, false, 7);
        for (label, alg) in [
            ("tree_kruskal", SpanningAlgorithm::Kruskal),
            ("tree_prim", SpanningAlgorithm::Prim),
            ("tree_boruvka", SpanningAlgorithm::Boruvka),
        ] {
            run(label, format!("n={n} m={}", g.edge_count()), 1, &mut || {
                spanning_tree(&g, true, false, alg).expect("tree");
            });
        }
        let d = random(n, 8.0, true, 8);
        run(
            "tree_arborescence",
            format!("n={n} m={}", d.edge_count()),
            1,
            &mut || {
                let _ = optimum_branching(&d, true, BranchingKind::MaximumBranching, None);
            },
        );
        run(
            "dag_dominators",
            format!("n={n} m={}", d.edge_count()),
            1,
            &mut || {
                immediate_dominators(&d, 0).expect("idom");
            },
        );
        let dag = random_dag(n, 8.0, 9);
        run(
            "dag_topological_sort",
            format!("n={n} m={}", dag.edge_count()),
            1,
            &mut || {
                topological_sort(&dag, NeighborOrder::Sorted).expect("order");
            },
        );
    }
    for n in sizes(&[1_000usize, big(1_000, 10_000)]) {
        let dag = random_dag(n, 8.0, 10);
        run(
            "dag_closure",
            format!("n={n} m={}", dag.edge_count()),
            1,
            &mut || {
                transitive_closure(&dag, Reflexive::Never).expect("closure");
            },
        );
        run(
            "dag_reduction",
            format!("n={n} m={}", dag.edge_count()),
            1,
            &mut || {
                transitive_reduction(&dag).expect("reduction");
            },
        );
    }

    // Cycles.
    for n in sizes(&[1_000usize, big(1_000, 10_000)]) {
        let g = random(n, 6.0, false, 11);
        run(
            "cycle_girth",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                girth(&g);
            },
        );
        run(
            "cycle_simple_bound5",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                simple_cycles(&g, Some(5), None).expect("cycles");
            },
        );
        run(
            "cycle_chordless_bound6",
            format!("n={n} m={}", g.edge_count()),
            1,
            &mut || {
                chordless_cycles(&g, Some(6), None).expect("cycles");
            },
        );
    }
    for n in sizes(&[100usize, 300, big(300, 1_000)]) {
        let g = random(n, 4.0, false, 12);
        for t in [1, 0] {
            run(
                "cycle_min_basis",
                format!("n={n} m={}", g.edge_count()),
                t,
                &mut || {
                    minimum_cycle_basis(&g, true, t).expect("basis");
                },
            );
        }
    }

    println!("| workload | size | threads | seconds | peak heap MB | allocations |");
    println!("|---|---|---|---|---|---|");
    for r in &rows {
        println!(
            "| {} | {} | {} | {:.4} | {:.1} | {} |",
            r.name, r.size, r.threads, r.seconds, r.heap_mb, r.allocs
        );
    }
}
