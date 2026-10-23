//! Benchmarks and profiles the community detection, subgroup, position,
//! layout, permutation test, random graph, and dense linear algebra
//! workloads, as a tool for finding what to optimize.
//!
//!     cargo run --release --offline --example bench_models
//!     MODE=all REPS=5 BENCH=louvain,qap cargo run --release --offline --example bench_models
//!
//! Environment variables:
//! - `MODE`: a comma-separated list of `fixed` (the default: every workload
//!   on sparse, medium, and dense networks, on one thread and on every
//!   core), `ladder` (empirical complexity: each workload over a geometric
//!   ladder of sizes at a fixed average degree, and at a fixed size over a
//!   ladder of densities, with log-log slopes of time and memory compared
//!   with the expected exponents), `scaling` (speedup and efficiency at 1,
//!   2, 4, ... threads up to the number of cores), or `all`.
//! - `BENCH`: comma-separated substrings; only workloads whose name
//!   contains one of them run.
//! - `REPS`: runs per measurement (default 3); the median is reported.
//! - `OUT`: output directory (default `target/bench_models`).
//! - `PERF=1`: also run the selected workloads under `perf stat` for
//!   hardware counters (skipped with a note when `perf` is absent).
//!
//! Every measurement records wall time, the peak of live heap bytes and
//! the number of allocations during the run (through a counting wrapper
//! around the system allocator), and the process's peak resident set size
//! (`VmHWM`, Linux only; it never decreases, so it bounds the run from
//! above). Long-running algorithms report per-phase times through
//! `knappnet::linalg::phase`. Results go to `results.csv`, `phases.csv`,
//! `ladder.csv`, `slopes.csv`, `scaling.csv`, and a Markdown `summary.md`
//! listing the largest opportunities.
//!
//! Density classes: sparse (average degree 8), medium (density 2%), and
//! dense (density 30%).

use std::alloc::{GlobalAlloc, Layout, System};
use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

use knappnet::generate::gnp;
use knappnet::layout::{stress_layout, StressAlgorithm, StressOptions};
use knappnet::linalg::{cholesky, matmul, phase, symmetric_eigen};
use knappnet::matrix::Dense;
use knappnet::metrics::cliques::max_cliques;
use knappnet::metrics::community::louvain;
use knappnet::metrics::coreperiphery::discrete::{core_periphery, CorePeripheryOptions};
use knappnet::metrics::equivalence::{concor, ConcorOptions};
use knappnet::metrics::factions::kgroups::{factions, FactionOptions};
use knappnet::metrics::partition::fast_greedy::fast_greedy;
use knappnet::metrics::partition::girvan_newman::girvan_newman;
use knappnet::metrics::partition::leading_eigen::leading_eigen;
use knappnet::metrics::regular::regular_equivalence;
use knappnet::metrics::scaling::classical_mds;
use knappnet::qap::{
    mrqap, qap_cor, MrqapMethod, MrqapOptions, Permutations, QapCorOptions, QapMode,
};
use knappnet::rng::Stream;
use knappnet::Graph;

// ---- Counting allocator --------------------------------------------------------

struct Counting;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static ALLOCS: AtomicU64 = AtomicU64::new(0);

fn grow(by: usize) {
    let now = LIVE.fetch_add(by, Ordering::Relaxed) + by;
    PEAK.fetch_max(now, Ordering::Relaxed);
}

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        if !p.is_null() {
            grow(layout.size());
            ALLOCS.fetch_add(1, Ordering::Relaxed);
        }
        p
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc_zeroed(layout);
        if !p.is_null() {
            grow(layout.size());
            ALLOCS.fetch_add(1, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        System.dealloc(ptr, layout);
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let p = System.realloc(ptr, layout, new_size);
        if !p.is_null() {
            ALLOCS.fetch_add(1, Ordering::Relaxed);
            if new_size >= layout.size() {
                grow(new_size - layout.size());
            } else {
                LIVE.fetch_sub(layout.size() - new_size, Ordering::Relaxed);
            }
        }
        p
    }
}

#[global_allocator]
static GLOBAL: Counting = Counting;

/// Peak resident set size in kB (Linux), or 0.
fn vm_hwm_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1).and_then(|x| x.parse().ok()))
        })
        .unwrap_or(0)
}

#[derive(Clone, Copy, Default)]
struct Sample {
    secs: f64,
    /// Peak live heap bytes above the level at the start of the run.
    peak: usize,
    allocs: u64,
}

fn measure(f: &mut dyn FnMut()) -> Sample {
    let base = LIVE.load(Ordering::Relaxed);
    PEAK.store(base, Ordering::Relaxed);
    let a0 = ALLOCS.load(Ordering::Relaxed);
    let t = Instant::now();
    f();
    let secs = t.elapsed().as_secs_f64();
    Sample {
        secs,
        peak: PEAK.load(Ordering::Relaxed).saturating_sub(base),
        allocs: ALLOCS.load(Ordering::Relaxed) - a0,
    }
}

/// Median time over `reps` runs (memory from the median run).
fn median(reps: usize, f: &mut dyn FnMut()) -> Sample {
    let mut s: Vec<Sample> = (0..reps.max(1)).map(|_| measure(f)).collect();
    s.sort_by(|a, b| a.secs.total_cmp(&b.secs));
    s[s.len() / 2]
}

// ---- Inputs -------------------------------------------------------------------------

/// Density of a class at size n.
fn class_density(class: &str, n: usize) -> f64 {
    match class {
        "sparse" => (8.0 / (n as f64 - 1.0)).min(1.0),
        "medium" => 0.02,
        "dense" => 0.30,
        _ => panic!("unknown class {class}"),
    }
}

fn random_graph(n: usize, p: f64, directed: bool, weighted: bool, seed: u64) -> Graph {
    let e = gnp(n, p, directed, false, 1, seed, 1).unwrap().remove(0);
    let src: Vec<i64> = e.from.iter().map(|&v| v as i64).collect();
    let tgt: Vec<i64> = e.to.iter().map(|&v| v as i64).collect();
    let w: Option<Vec<f64>> = weighted.then(|| {
        let mut r = Stream::new(seed, 99);
        (0..src.len()).map(|_| 1.0 + r.below(5) as f64).collect()
    });
    Graph::from_edges(n, &src, &tgt, w.as_deref(), directed).unwrap()
}

fn random_matrix(n: usize, seed: u64, symmetric: bool) -> Dense {
    let mut r = Stream::new(seed, 7);
    let mut m = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            if symmetric && j > i {
                continue;
            }
            let x = r.uniform();
            m.set(i, j, x);
            if symmetric {
                m.set(j, i, x);
            }
        }
    }
    m
}

fn spd_matrix(n: usize, seed: u64) -> Dense {
    let a = random_matrix(n, seed, false);
    let mut m = matmul(&knappnet::linalg::transpose(&a), &a);
    for i in 0..n {
        m.set(i, i, m.get(i, i) + n as f64);
    }
    m
}

fn distance_matrix(n: usize) -> Dense {
    let mut d = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let x = ((i as f64 - j as f64).abs() + 1.0).ln();
            d.set(i, j, if i == j { 0.0 } else { x });
        }
    }
    d
}

// ---- Workloads ----------------------------------------------------------------------

/// A prepared workload: its "size" m (edges, or matrix entries) and a
/// runner taking the thread count.
struct Prepared {
    m: usize,
    run: Box<dyn FnMut(usize)>,
}

/// A workload family.
struct Spec {
    name: &'static str,
    /// Whether the workload takes a thread count.
    parallel: bool,
    /// Whether the workload depends on density (a network input).
    graph: bool,
    /// Fixed-mode cases: (class, n).
    cases: &'static [(&'static str, usize)],
    /// Ladder at average degree 8 (or plain sizes for matrix workloads).
    ladder: &'static [usize],
    /// Density ladder: fixed n and densities.
    density_n: usize,
    densities: &'static [f64],
    /// Expected exponents of time in n (fixed average degree) and in m
    /// (fixed n, growing density); NaN when no single exponent applies.
    exp_n: f64,
    exp_m: f64,
    /// Why those exponents are expected.
    expect: &'static str,
    setup: fn(usize, f64) -> Prepared,
}

fn graph_job(
    n: usize,
    p: f64,
    directed: bool,
    weighted: bool,
    seed: u64,
    f: fn(&Graph, usize),
) -> Prepared {
    let g = random_graph(n, p, directed, weighted, seed);
    Prepared {
        m: g.edge_count(),
        run: Box::new(move |t| f(&g, t)),
    }
}

const DENS: &[f64] = &[0.01, 0.02, 0.05, 0.1, 0.2, 0.4];
const NONE: &[f64] = &[];

fn specs() -> Vec<Spec> {
    vec![
        Spec {
            name: "louvain",
            parallel: false,
            graph: true,
            cases: &[("sparse", 20_000), ("medium", 3_000), ("dense", 1_000)],
            ladder: &[2_500, 5_000, 10_000, 20_000, 40_000],
            density_n: 1_000,
            densities: DENS,
            exp_n: 1.0,
            exp_m: 1.0,
            expect: "near-linear in m per level",
            setup: |n, p| {
                graph_job(n, p, false, false, 1, |g, _| {
                    louvain(g, 1.0, false, 7).unwrap();
                })
            },
        },
        Spec {
            name: "louvain weighted",
            parallel: false,
            graph: true,
            cases: &[("sparse", 20_000), ("medium", 3_000), ("dense", 1_000)],
            ladder: &[],
            density_n: 0,
            densities: NONE,
            exp_n: 1.0,
            exp_m: 1.0,
            expect: "near-linear in m per level",
            setup: |n, p| {
                graph_job(n, p, false, true, 1, |g, _| {
                    louvain(g, 1.0, true, 7).unwrap();
                })
            },
        },
        Spec {
            name: "fast greedy",
            parallel: false,
            graph: true,
            cases: &[("sparse", 10_000), ("medium", 2_000), ("dense", 600)],
            ladder: &[1_250, 2_500, 5_000, 10_000],
            density_n: 600,
            densities: DENS,
            exp_n: 1.3,
            exp_m: 1.3,
            expect: "O(m d log n), Clauset et al. 2004",
            setup: |n, p| {
                graph_job(n, p, false, false, 2, |g, _| {
                    fast_greedy(g, 1.0, false).unwrap();
                })
            },
        },
        Spec {
            name: "leading eigenvector",
            parallel: false,
            graph: true,
            cases: &[("sparse", 1_000), ("medium", 1_000), ("dense", 1_000)],
            ladder: &[125, 250, 500, 1_000],
            density_n: 500,
            densities: DENS,
            exp_n: 3.0,
            exp_m: 0.0,
            expect: "dense eigendecomposition, O(n^3)",
            setup: |n, p| {
                graph_job(n, p, false, false, 3, |g, _| {
                    leading_eigen(g, false).unwrap();
                })
            },
        },
        Spec {
            name: "girvan-newman",
            parallel: true,
            graph: true,
            cases: &[("sparse", 300), ("dense", 60)],
            ladder: &[50, 100, 200],
            density_n: 60,
            densities: &[0.05, 0.1, 0.2, 0.4],
            exp_n: 3.0,
            exp_m: 2.0,
            expect: "m removals of O(nm) betweenness",
            setup: |n, p| {
                graph_job(n, p, false, false, 4, |g, t| {
                    girvan_newman(g, t).unwrap();
                })
            },
        },
        Spec {
            name: "maximal cliques",
            parallel: true,
            graph: true,
            cases: &[("sparse", 20_000), ("medium", 3_000), ("dense", 150)],
            ladder: &[2_500, 5_000, 10_000, 20_000, 40_000],
            density_n: 150,
            densities: &[0.05, 0.1, 0.2, 0.3],
            exp_n: 1.0,
            exp_m: f64::NAN,
            expect: "linear for bounded degeneracy; output-sensitive in density",
            setup: |n, p| {
                graph_job(n, p, false, false, 5, |g, t| {
                    max_cliques(g, 3, t);
                })
            },
        },
        Spec {
            name: "factions k=3 undirected",
            parallel: true,
            graph: true,
            cases: &[("sparse", 400), ("dense", 200)],
            ladder: &[100, 200, 400],
            density_n: 200,
            densities: &[0.02, 0.05, 0.1, 0.2, 0.4],
            exp_n: f64::NAN,
            exp_m: f64::NAN,
            expect: "tabu search; iterations depend on the landscape",
            setup: |n, p| {
                graph_job(n, p, false, false, 6, |g, t| {
                    let opts = FactionOptions {
                        k: 3,
                        threads: t,
                        ..Default::default()
                    };
                    factions(g, opts).unwrap();
                })
            },
        },
        Spec {
            name: "factions k=3 directed",
            parallel: true,
            graph: true,
            cases: &[("sparse", 400), ("dense", 200)],
            ladder: &[],
            density_n: 0,
            densities: NONE,
            exp_n: f64::NAN,
            exp_m: f64::NAN,
            expect: "tabu search",
            setup: |n, p| {
                graph_job(n, p, true, false, 6, |g, t| {
                    let opts = FactionOptions {
                        k: 3,
                        threads: t,
                        ..Default::default()
                    };
                    factions(g, opts).unwrap();
                })
            },
        },
        Spec {
            name: "core-periphery",
            parallel: true,
            graph: true,
            cases: &[("sparse", 400), ("dense", 200)],
            ladder: &[100, 200, 400],
            density_n: 200,
            densities: &[0.02, 0.05, 0.1, 0.2, 0.4],
            exp_n: f64::NAN,
            exp_m: f64::NAN,
            expect: "tabu search",
            setup: |n, p| {
                graph_job(n, p, false, false, 7, |g, t| {
                    let opts = CorePeripheryOptions {
                        threads: t,
                        ..Default::default()
                    };
                    core_periphery(g, opts).unwrap();
                })
            },
        },
        Spec {
            name: "concor depth 3 directed",
            parallel: true,
            graph: true,
            cases: &[("sparse", 2_000), ("dense", 1_000)],
            ladder: &[250, 500, 1_000],
            density_n: 500,
            densities: DENS,
            exp_n: 3.0,
            exp_m: 0.0,
            expect: "column correlations, O(n^3) per iteration",
            setup: |n, p| {
                graph_job(n, p, true, false, 8, |g, t| {
                    let opts = ConcorOptions {
                        depth: 3,
                        threads: t,
                        ..Default::default()
                    };
                    concor(g, opts).unwrap();
                })
            },
        },
        Spec {
            name: "regular equivalence directed",
            parallel: false,
            graph: true,
            cases: &[("sparse", 2_000), ("dense", 1_000)],
            ladder: &[250, 500, 1_000],
            density_n: 300,
            densities: DENS,
            exp_n: 2.0,
            exp_m: 2.0,
            expect: "REGE, O(n^2 d^2) per iteration",
            setup: |n, p| {
                graph_job(n, p, true, false, 8, |g, _| {
                    regular_equivalence(g, None, Default::default()).unwrap();
                })
            },
        },
        Spec {
            name: "stress layout all pairs",
            parallel: true,
            graph: true,
            cases: &[("sparse", 1_000), ("dense", 1_000)],
            ladder: &[250, 500, 1_000, 2_000],
            density_n: 1_000,
            densities: DENS,
            exp_n: 2.0,
            exp_m: 1.0,
            expect: "O(n^2) per iteration plus O(nm) distances",
            setup: |n, p| {
                graph_job(n, p, false, false, 9, |g, t| {
                    let opts = StressOptions {
                        algorithm: StressAlgorithm::Full,
                        threads: t,
                        ..Default::default()
                    };
                    stress_layout(g, &opts).unwrap();
                })
            },
        },
        Spec {
            name: "sparse stress",
            parallel: true,
            graph: true,
            cases: &[("sparse", 10_000)],
            ladder: &[1_250, 2_500, 5_000, 10_000, 20_000],
            density_n: 2_000,
            densities: &[0.005, 0.01, 0.02, 0.05],
            exp_n: 1.0,
            exp_m: 1.0,
            expect: "O(k m) pivot distances and O(m + kn) per iteration",
            setup: |n, p| {
                graph_job(n, p, false, false, 10, |g, t| {
                    let opts = StressOptions {
                        algorithm: StressAlgorithm::Sparse,
                        threads: t,
                        ..Default::default()
                    };
                    stress_layout(g, &opts).unwrap();
                })
            },
        },
        Spec {
            name: "qap correlation digraph 1000 permutations",
            parallel: true,
            graph: false,
            cases: &[("matrix", 100)],
            ladder: &[50, 100, 200, 400],
            density_n: 0,
            densities: NONE,
            exp_n: 2.0,
            exp_m: f64::NAN,
            expect: "O(n^2) per permutation",
            setup: |n, _| {
                let y = random_matrix(n, 11, false);
                let x = random_matrix(n, 12, false);
                Prepared {
                    m: n * n,
                    run: Box::new(move |t| {
                        let opts = QapCorOptions {
                            mode: QapMode::Digraph,
                            permutations: Permutations::Random {
                                reps: 1000,
                                seed: 1,
                            },
                            threads: t,
                            ..Default::default()
                        };
                        qap_cor(&y, &x, opts).unwrap();
                    }),
                }
            },
        },
        Spec {
            name: "mrqap dsp digraph 3 predictors 1000 permutations",
            parallel: true,
            graph: false,
            cases: &[("matrix", 100)],
            ladder: &[50, 100, 200],
            density_n: 0,
            densities: NONE,
            exp_n: 2.0,
            exp_m: f64::NAN,
            expect: "O(p^2 n^2) per permutation",
            setup: |n, _| mrqap_job(n, MrqapMethod::DoubleSemiPartialling),
        },
        Spec {
            name: "mrqap y-perm digraph 3 predictors 1000 permutations",
            parallel: true,
            graph: false,
            cases: &[("matrix", 100)],
            ladder: &[50, 100, 200],
            density_n: 0,
            densities: NONE,
            exp_n: 2.0,
            exp_m: f64::NAN,
            expect: "O(p^2 n^2) per permutation",
            setup: |n, _| mrqap_job(n, MrqapMethod::YPermutation),
        },
        Spec {
            name: "gnp undirected 20 draws",
            parallel: true,
            graph: true,
            cases: &[("sparse", 2_000), ("dense", 2_000)],
            ladder: &[2_000, 8_000, 32_000, 128_000],
            density_n: 2_000,
            densities: DENS,
            exp_n: 1.0,
            exp_m: 1.0,
            expect: "geometric skipping, O(n + m)",
            setup: |n, p| Prepared {
                m: (p * (n * (n - 1) / 2) as f64) as usize,
                run: Box::new(move |t| {
                    gnp(n, p, false, false, 20, 3, t).unwrap();
                }),
            },
        },
        Spec {
            name: "gnp directed 1000 draws",
            parallel: true,
            graph: true,
            cases: &[("sparse", 100), ("dense", 100)],
            ladder: &[],
            density_n: 0,
            densities: NONE,
            exp_n: 1.0,
            exp_m: 1.0,
            expect: "geometric skipping, O(n + m)",
            setup: |n, p| Prepared {
                m: (p * (n * (n - 1)) as f64) as usize,
                run: Box::new(move |t| {
                    gnp(n, p, true, false, 1000, 3, t).unwrap();
                }),
            },
        },
        Spec {
            name: "symmetric eigendecomposition",
            parallel: false,
            graph: false,
            cases: &[("matrix", 300), ("matrix", 1_000)],
            ladder: &[125, 250, 500, 1_000],
            density_n: 0,
            densities: NONE,
            exp_n: 3.0,
            exp_m: f64::NAN,
            expect: "O(n^3)",
            setup: |n, _| {
                let a = random_matrix(n, 21, true);
                Prepared {
                    m: n * n,
                    run: Box::new(move |_| {
                        symmetric_eigen(&a).unwrap();
                    }),
                }
            },
        },
        Spec {
            name: "cholesky",
            parallel: false,
            graph: false,
            cases: &[("matrix", 300), ("matrix", 1_000)],
            ladder: &[125, 250, 500, 1_000],
            density_n: 0,
            densities: NONE,
            exp_n: 3.0,
            exp_m: f64::NAN,
            expect: "O(n^3)",
            setup: |n, _| {
                let s = spd_matrix(n, 22);
                Prepared {
                    m: n * n,
                    run: Box::new(move |_| {
                        cholesky(&s).unwrap();
                    }),
                }
            },
        },
        Spec {
            name: "matmul",
            parallel: false,
            graph: false,
            cases: &[("matrix", 300), ("matrix", 1_000)],
            ladder: &[125, 250, 500, 1_000],
            density_n: 0,
            densities: NONE,
            exp_n: 3.0,
            exp_m: f64::NAN,
            expect: "O(n^3)",
            setup: |n, _| {
                let a = random_matrix(n, 23, false);
                let b = random_matrix(n, 24, false);
                Prepared {
                    m: n * n,
                    run: Box::new(move |_| {
                        matmul(&a, &b);
                    }),
                }
            },
        },
        Spec {
            name: "classical mds k=2",
            parallel: false,
            graph: false,
            cases: &[("matrix", 300), ("matrix", 1_000)],
            ladder: &[125, 250, 500, 1_000],
            density_n: 0,
            densities: NONE,
            exp_n: 3.0,
            exp_m: f64::NAN,
            expect: "dense eigendecomposition, O(n^3)",
            setup: |n, _| {
                let d = distance_matrix(n);
                Prepared {
                    m: n * n,
                    run: Box::new(move |_| {
                        classical_mds(&d, 2).unwrap();
                    }),
                }
            },
        },
    ]
}

fn mrqap_job(n: usize, method: MrqapMethod) -> Prepared {
    let y = random_matrix(n, 11, false);
    let xs = vec![
        random_matrix(n, 12, false),
        random_matrix(n, 13, false),
        random_matrix(n, 14, false),
    ];
    Prepared {
        m: n * n,
        run: Box::new(move |t| {
            let opts = MrqapOptions {
                mode: QapMode::Digraph,
                method,
                permutations: Permutations::Random {
                    reps: 1000,
                    seed: 1,
                },
                threads: t,
                ..Default::default()
            };
            mrqap(&y, &xs, opts).unwrap();
        }),
    }
}

// ---- Driver -------------------------------------------------------------------------

struct Ctx {
    reps: usize,
    cores: usize,
    out: std::path::PathBuf,
    results: String,
    phases: String,
    ladder: String,
    slopes: String,
    scaling: String,
    /// (workload, class, n, m, 1-thread seconds, peak bytes)
    fixed_rows: Vec<(String, String, usize, usize, f64, usize)>,
    phase_rows: Vec<(String, String, String, f64)>,
    flags: Vec<String>,
    scale_rows: Vec<(String, usize, f64)>,
}

fn p_of(spec: &Spec, class: &str, n: usize) -> f64 {
    if spec.graph {
        class_density(class, n)
    } else {
        0.0
    }
}

fn run_fixed(ctx: &mut Ctx, spec: &Spec) {
    for &(class, n) in spec.cases {
        let p = p_of(spec, class, n);
        let mut job = (spec.setup)(n, p);
        let m = job.m;
        let threads: Vec<usize> = if spec.parallel && ctx.cores > 1 {
            vec![1, ctx.cores]
        } else {
            vec![1]
        };
        let mut line = format!("{:<52} {:<6} n={:<6} m={:<8}", spec.name, class, n, m);
        for &t in &threads {
            if t == 1 {
                phase::take();
                phase::enable(true);
            }
            let s = median(ctx.reps, &mut || (job.run)(t));
            if t == 1 {
                phase::enable(false);
                for (ph, secs, _) in phase::take() {
                    let per = secs / ctx.reps as f64;
                    let _ = writeln!(
                        ctx.phases,
                        "{},{},{},{},{:.6},{:.4}",
                        spec.name,
                        class,
                        n,
                        ph,
                        per,
                        per / s.secs
                    );
                    ctx.phase_rows.push((
                        spec.name.into(),
                        format!("{class} n={n}"),
                        ph.into(),
                        per / s.secs,
                    ));
                }
                ctx.fixed_rows
                    .push((spec.name.into(), class.into(), n, m, s.secs, s.peak));
            }
            let hwm = vm_hwm_kb();
            let _ = writeln!(
                ctx.results,
                "{},{},{},{},{},{:.6},{},{},{:.1},{:.1},{}",
                spec.name,
                class,
                n,
                m,
                t,
                s.secs,
                s.peak,
                s.allocs,
                s.peak as f64 / n as f64,
                s.peak as f64 / m.max(1) as f64,
                hwm
            );
            let _ = write!(
                line,
                "  t={t}: {:>8.3}s {:>9.1} MB {:>9} allocs",
                s.secs,
                s.peak as f64 / 1e6,
                s.allocs
            );
        }
        println!("{line}");
    }
}

/// Least-squares slope of log y on log x.
fn slope(points: &[(f64, f64)]) -> f64 {
    let pts: Vec<(f64, f64)> = points
        .iter()
        .filter(|p| p.0 > 0.0 && p.1 > 0.0)
        .map(|p| (p.0.ln(), p.1.ln()))
        .collect();
    if pts.len() < 2 {
        return f64::NAN;
    }
    let k = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / k;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / k;
    let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
    sxy / sxx
}

fn run_ladder(ctx: &mut Ctx, spec: &Spec) {
    let mut ladders: Vec<(&str, Vec<(usize, f64)>)> = Vec::new();
    if !spec.ladder.is_empty() {
        ladders.push((
            "n",
            spec.ladder
                .iter()
                .map(|&n| (n, p_of(spec, "sparse", n)))
                .collect(),
        ));
    }
    if !spec.densities.is_empty() {
        ladders.push((
            "density",
            spec.densities
                .iter()
                .map(|&p| (spec.density_n, p))
                .collect(),
        ));
    }
    for (kind, points) in ladders {
        let mut time_n = Vec::new();
        let mut time_m = Vec::new();
        let mut mem_n = Vec::new();
        let mut mem_m = Vec::new();
        for (n, p) in points {
            let mut job = (spec.setup)(n, p);
            let s = median(ctx.reps, &mut || (job.run)(1));
            let _ = writeln!(
                ctx.ladder,
                "{},{},{},{:.5},{},{:.6},{},{}",
                spec.name, kind, n, p, job.m, s.secs, s.peak, s.allocs
            );
            println!(
                "  ladder {:<44} {kind:<7} n={n:<7} p={p:<8.4} m={:<9} {:>8.3}s {:>9.1} MB",
                spec.name,
                job.m,
                s.secs,
                s.peak as f64 / 1e6
            );
            // Runs under a millisecond are timer noise.
            if s.secs >= 1e-3 {
                time_n.push((n as f64, s.secs));
                time_m.push((job.m as f64, s.secs));
            }
            mem_n.push((n as f64, s.peak as f64));
            mem_m.push((job.m as f64, s.peak as f64));
        }
        let (st, sm, expected) = if kind == "n" {
            (slope(&time_n), slope(&mem_n), spec.exp_n)
        } else {
            (slope(&time_m), slope(&mem_m), spec.exp_m)
        };
        let flagged = expected.is_finite() && st.is_finite() && st > expected + 0.25;
        let _ = writeln!(
            ctx.slopes,
            "{},{},{:.3},{:.3},{},{}",
            spec.name,
            kind,
            st,
            sm,
            if expected.is_finite() {
                format!("{expected:.2}")
            } else {
                String::new()
            },
            flagged
        );
        println!(
            "  slope  {:<44} time ~ {}^{st:.2} (expected {}), memory ~ {}^{sm:.2}{}",
            spec.name,
            if kind == "n" { "n" } else { "m" },
            if expected.is_finite() {
                format!("{expected:.2}")
            } else {
                "-".into()
            },
            if kind == "n" { "n" } else { "m" },
            if flagged {
                "  <-- above expectation"
            } else {
                ""
            }
        );
        if flagged {
            ctx.flags.push(format!(
                "{} ({} ladder): time grows as {}^{st:.2}, expected {expected:.2} ({})",
                spec.name,
                kind,
                if kind == "n" { "n" } else { "m" },
                spec.expect
            ));
        }
    }
}

fn run_scaling(ctx: &mut Ctx, spec: &Spec) {
    if !spec.parallel {
        return;
    }
    let mut ts = vec![1usize];
    while ts.last().unwrap() * 2 <= ctx.cores {
        ts.push(ts.last().unwrap() * 2);
    }
    if *ts.last().unwrap() != ctx.cores {
        ts.push(ctx.cores);
    }
    let &(class, n) = spec.cases.last().unwrap();
    let mut job = (spec.setup)(n, p_of(spec, class, n));
    let mut base = 0.0;
    for t in ts {
        let s = median(ctx.reps, &mut || (job.run)(t));
        if t == 1 {
            base = s.secs;
        }
        let speedup = base / s.secs;
        let eff = speedup / t as f64;
        let _ = writeln!(
            ctx.scaling,
            "{},{},{},{},{:.6},{:.3},{:.3}",
            spec.name, class, n, t, s.secs, speedup, eff
        );
        println!(
            "  scaling {:<44} {class} n={n} t={t}: {:.3}s speedup {speedup:.2} efficiency {eff:.2}",
            spec.name, s.secs
        );
        if t > 1 {
            ctx.scale_rows
                .push((format!("{} ({class} n={n})", spec.name), t, eff));
        }
    }
}

fn run_perf(filter: &str) -> String {
    let ok = std::process::Command::new("perf")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !ok {
        return "`perf` is not available on this machine; hardware counters skipped.\n".into();
    }
    let exe = std::env::current_exe().unwrap();
    let out = std::process::Command::new("perf")
        .args([
            "stat",
            "-x,",
            "-e",
            "cycles,instructions,cache-references,cache-misses,branch-misses",
        ])
        .arg(&exe)
        .env("PERF", "0")
        .env("MODE", "fixed")
        .env("REPS", "1")
        .env("BENCH", filter)
        .output();
    match out {
        Ok(o) => String::from_utf8_lossy(&o.stderr).into_owned(),
        Err(e) => format!("perf stat failed: {e}\n"),
    }
}

fn main() {
    let filter: Vec<String> = std::env::var("BENCH")
        .ok()
        .map(|s| s.split(',').map(|x| x.trim().to_string()).collect())
        .unwrap_or_default();
    let modes: Vec<String> = std::env::var("MODE")
        .unwrap_or_else(|_| "fixed".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let has = |m: &str| modes.iter().any(|x| x == m || x == "all");
    let reps: usize = std::env::var("REPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    let out = std::path::PathBuf::from(
        std::env::var("OUT").unwrap_or_else(|_| "target/bench_models".into()),
    );
    std::fs::create_dir_all(&out).expect("output directory");
    let mut ctx = Ctx {
        reps,
        cores,
        out,
        results: "workload,class,n,m,threads,seconds,peak_heap_bytes,allocations,bytes_per_node,bytes_per_edge,vm_hwm_kb\n".into(),
        phases: "workload,class,n,phase,seconds_per_run,share_of_wall\n".into(),
        ladder: "workload,ladder,n,density,m,seconds,peak_heap_bytes,allocations\n".into(),
        slopes: "workload,ladder,time_slope,memory_slope,expected_time_slope,above_expectation\n".into(),
        scaling: "workload,class,n,threads,seconds,speedup,efficiency\n".into(),
        fixed_rows: Vec::new(),
        phase_rows: Vec::new(),
        flags: Vec::new(),
        scale_rows: Vec::new(),
    };
    println!("{cores} cores, {reps} runs per measurement (median)");
    let specs: Vec<Spec> = specs()
        .into_iter()
        .filter(|s| filter.is_empty() || filter.iter().any(|f| s.name.contains(f.as_str())))
        .collect();
    for spec in &specs {
        if has("fixed") {
            run_fixed(&mut ctx, spec);
        }
        if has("ladder") {
            run_ladder(&mut ctx, spec);
        }
        if has("scaling") {
            run_scaling(&mut ctx, spec);
        }
    }
    let perf = if std::env::var("PERF").is_ok_and(|v| v == "1") {
        run_perf(&filter.join(","))
    } else {
        String::new()
    };

    // ---- Summary ------------------------------------------------------------------
    let mut md = String::from("# Model workload profile\n\n");
    let _ = writeln!(md, "{cores} cores, median of {reps} runs.\n");
    if !ctx.fixed_rows.is_empty() {
        md.push_str("## Slowest workloads (one thread)\n\n| workload | class | n | m | seconds | peak heap MB | bytes/edge |\n|---|---|---|---|---|---|---|\n");
        let mut rows = ctx.fixed_rows.clone();
        rows.sort_by(|a, b| b.4.total_cmp(&a.4));
        for r in rows.iter().take(12) {
            let _ = writeln!(
                md,
                "| {} | {} | {} | {} | {:.3} | {:.1} | {:.0} |",
                r.0,
                r.1,
                r.2,
                r.3,
                r.4,
                r.5 as f64 / 1e6,
                r.5 as f64 / r.3.max(1) as f64
            );
        }
        md.push('\n');
    }
    if !ctx.phase_rows.is_empty() {
        md.push_str("## Dominant phases (share of wall time, one thread)\n\n");
        let mut ph = ctx.phase_rows.clone();
        ph.sort_by(|a, b| b.3.total_cmp(&a.3));
        for r in ph.iter().filter(|r| r.3 >= 0.3).take(15) {
            let _ = writeln!(md, "- {} ({}): {} {:.0}%", r.0, r.1, r.2, 100.0 * r.3);
        }
        md.push('\n');
    }
    if has("ladder") {
        md.push_str("## Empirical complexity above expectation\n\n");
        if ctx.flags.is_empty() {
            md.push_str("None.\n");
        }
        for f in &ctx.flags {
            let _ = writeln!(md, "- {f}");
        }
        md.push('\n');
    }
    if !ctx.scale_rows.is_empty() {
        md.push_str("## Parallel efficiency below 60%\n\n");
        let mut any = false;
        for r in ctx.scale_rows.iter().filter(|r| r.2 < 0.6) {
            any = true;
            let _ = writeln!(md, "- {} at {} threads: {:.0}%", r.0, r.1, 100.0 * r.2);
        }
        if !any {
            md.push_str("None.\n");
        }
        md.push('\n');
    }
    if !perf.is_empty() {
        md.push_str("## Hardware counters\n\n```\n");
        md.push_str(&perf);
        md.push_str("```\n");
    }
    let write = |name: &str, body: &str| {
        std::fs::write(ctx.out.join(name), body).expect("write output");
    };
    write("results.csv", &ctx.results);
    write("phases.csv", &ctx.phases);
    write("ladder.csv", &ctx.ladder);
    write("slopes.csv", &ctx.slopes);
    write("scaling.csv", &ctx.scaling);
    write("summary.md", &md);
    println!("\n{md}");
    println!("Wrote {}", ctx.out.display());
}
