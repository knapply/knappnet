//! Times the shortest-path API, the centrality variants, and the distance
//! measures on random networks, on one thread and on every core, for the
//! scale notes (docs/drafts/scale_notes_paths.md).
//!
//!     cargo run --release --offline --example bench_paths_extra
//!     BENCH=johnson,floyd cargo run --release --offline --example bench_paths_extra
//!
//! `BENCH` is a comma-separated list of substrings; only workloads whose
//! name contains one of them run. `SCALE` multiplies the sizes (default 1).

use std::time::Instant;

use knappnet::metrics::communicability::communicability_betweenness;
use knappnet::metrics::current_flow::{
    current_flow_betweenness, current_flow_betweenness_approx, FlowOptions, Sample,
};
use knappnet::metrics::descriptive::PathOptions;
use knappnet::metrics::distance_measures::{closeness_vitality, local_efficiency, wiener_index};
use knappnet::metrics::path_centrality::{
    betweenness_subset, load_centrality, LoadOptions, SubsetOptions,
};
use knappnet::paths::dense::floyd_warshall;
use knappnet::paths::negative::{bellman_ford, goldberg_radzik, johnson};
use knappnet::paths::search::{all_pairs, single_source, SearchOptions};
use knappnet::paths::simple::k_shortest_simple_paths;
use knappnet::paths::{ArcLengths, PathGraph, PathGraphOptions};
use knappnet::{Graph, WeightKind};

struct Xor(u64);
impl Xor {
    fn below(&mut self, k: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % k
    }
}

/// A random network with `n` nodes, about `deg` edges per node, integer
/// weights 1..8, and a spanning path so that it is connected.
fn network(n: usize, deg: usize, directed: bool, seed: u64) -> Graph {
    let mut r = Xor(seed);
    let mut s: Vec<i64> = (0..n as i64 - 1).collect();
    let mut t: Vec<i64> = (1..n as i64).collect();
    if directed {
        s.push(n as i64 - 1);
        t.push(0);
    }
    for _ in 0..n * deg {
        s.push(r.below(n as u64) as i64);
        t.push(r.below(n as u64) as i64);
    }
    let w: Vec<f64> = (0..s.len()).map(|_| (1 + r.below(8)) as f64).collect();
    Graph::from_edges(n, &s, &t, Some(&w), directed).unwrap()
}

fn path_graph(g: &Graph, lengths: ArcLengths) -> PathGraph {
    PathGraph::new(
        g,
        PathGraphOptions {
            lengths,
            ..Default::default()
        },
    )
    .unwrap()
}

fn time(label: &str, wanted: &[String], mut f: impl FnMut(usize)) {
    if !wanted.is_empty() && !wanted.iter().any(|w| label.contains(w.as_str())) {
        return;
    }
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut line = format!("{label:<48}");
    for threads in [1, cores] {
        let t = Instant::now();
        f(threads);
        line += &format!("{:>10.3} s ({threads:>2} thr)", t.elapsed().as_secs_f64());
    }
    println!("{line}");
}

fn main() {
    let wanted: Vec<String> = std::env::var("BENCH")
        .map(|b| b.split(',').map(|s| s.trim().to_string()).collect())
        .unwrap_or_default();
    let scale: f64 = std::env::var("SCALE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1.0);
    let sz = |n: usize| ((n as f64) * scale) as usize;
    let w = ArcLengths::Weights(WeightKind::Distance);

    let big = network(sz(100_000), 4, true, 1);
    let pg = path_graph(&big, w);
    let unit = path_graph(&big, ArcLengths::Unit);
    println!(
        "sparse directed network: n = {}, m = {}",
        big.node_count(),
        big.edge_count()
    );
    time("single-source Dijkstra (n = 100k)", &wanted, |_| {
        single_source(&pg, 0, SearchOptions::default()).unwrap();
    });
    let sources: Vec<u32> = (0..500).collect();
    time("all pairs, 500 sources, BFS (n = 100k)", &wanted, |t| {
        all_pairs(&unit, Some(&sources), f64::INFINITY, t).unwrap();
    });
    time(
        "all pairs, 500 sources, Dijkstra (n = 100k)",
        &wanted,
        |t| {
            all_pairs(&pg, Some(&sources), f64::INFINITY, t).unwrap();
        },
    );
    // Negative lengths from a potential p(v) = v mod 7.
    let (src, tgt) = (big.edge_sources(), big.edge_targets());
    let wts = big.edge_weights().unwrap();
    let neg: Vec<f64> = (0..big.edge_count())
        .map(|e| wts[e] + (src[e] % 7) as f64 - (tgt[e] % 7) as f64)
        .collect();
    let npg = path_graph(&big, ArcLengths::PerEdge(&neg));
    time("Bellman-Ford, negative lengths (n = 100k)", &wanted, |_| {
        bellman_ford(&npg, &[0]).unwrap();
    });
    time(
        "Goldberg-Radzik, negative lengths (n = 100k)",
        &wanted,
        |_| {
            goldberg_radzik(&npg, 0).unwrap();
        },
    );
    time(
        "k shortest simple paths, k = 50 (n = 100k)",
        &wanted,
        |_| {
            k_shortest_simple_paths(&pg, 0, (big.node_count() / 2) as u32, 50).unwrap();
        },
    );
    let samples: Vec<u32> = (0..200).map(|i| i * 37 % big.node_count() as u32).collect();
    time(
        "betweenness from 200 sampled sources (n = 100k)",
        &wanted,
        |t| {
            betweenness_subset(
                &big,
                SubsetOptions {
                    sources: Some(&samples),
                    weighted: true,
                    threads: t,
                    ..Default::default()
                },
            )
            .unwrap();
        },
    );

    let mid = network(sz(10_000), 4, true, 2);
    let mpg = path_graph(&mid, w);
    let (src, tgt) = (mid.edge_sources(), mid.edge_targets());
    let wts = mid.edge_weights().unwrap();
    let neg: Vec<f64> = (0..mid.edge_count())
        .map(|e| wts[e] + (src[e] % 7) as f64 - (tgt[e] % 7) as f64)
        .collect();
    let mneg = path_graph(&mid, ArcLengths::PerEdge(&neg));
    println!(
        "medium directed network: n = {}, m = {}",
        mid.node_count(),
        mid.edge_count()
    );
    time("all pairs, Dijkstra (n = 10k)", &wanted, |t| {
        all_pairs(&mpg, None, f64::INFINITY, t).unwrap();
    });
    time("Johnson, negative lengths (n = 10k)", &wanted, |t| {
        johnson(&mneg, None, t).unwrap().unwrap();
    });
    time("load centrality (n = 10k)", &wanted, |t| {
        load_centrality(
            &mid,
            LoadOptions {
                threads: t,
                ..Default::default()
            },
        )
        .unwrap();
    });
    time("Wiener index (n = 10k)", &wanted, |t| {
        wiener_index(
            &mid,
            PathOptions {
                threads: t,
                ..Default::default()
            },
        )
        .unwrap();
    });

    let small = network(sz(2_000), 4, false, 3);
    let spg = path_graph(&small, w);
    println!(
        "small undirected network: n = {}, m = {}",
        small.node_count(),
        small.edge_count()
    );
    time("Floyd-Warshall (n = 2k)", &wanted, |t| {
        floyd_warshall(&spg, t);
    });
    time("current-flow betweenness, exact (n = 2k)", &wanted, |t| {
        current_flow_betweenness(
            &small,
            FlowOptions {
                threads: t,
                ..Default::default()
            },
        )
        .unwrap();
    });
    let pairs: Vec<(u32, u32)> = (0..400u32)
        .map(|i| (i * 7 % 2000, (i * 13 + 1) % 2000))
        .filter(|p| p.0 != p.1)
        .collect();
    time(
        "current-flow betweenness, 400 pairs (n = 2k)",
        &wanted,
        |t| {
            current_flow_betweenness_approx(
                &small,
                Sample::Pairs(&pairs),
                FlowOptions {
                    threads: t,
                    ..Default::default()
                },
            )
            .unwrap();
        },
    );
    time("closeness vitality (n = 2k)", &wanted, |t| {
        closeness_vitality(
            &small,
            PathOptions {
                threads: t,
                ..Default::default()
            },
        )
        .unwrap();
    });
    let huge = network(sz(200_000), 4, false, 5);
    time("local efficiency (n = 200k)", &wanted, |t| {
        local_efficiency(&huge, false, t).unwrap();
    });
    let tiny = network(sz(200), 3, false, 4);
    time("communicability betweenness (n = 200)", &wanted, |t| {
        communicability_betweenness(&tiny, t).unwrap();
    });
}
