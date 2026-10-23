//! Randomized property tests for the traversal kernels behind the
//! path-based measures: the bit-parallel multi-source breadth-first search
//! (`then2014more`) must give exactly what one plain breadth-first search
//! per source gives, the indexed-heap Dijkstra (`johnson1975priority`)
//! must agree with breadth-first search on unit lengths, and results must
//! not depend on the thread count.

extern crate knappnet as engine;

use engine::bitset::{for_each_bit, MultiBfs, Step};
use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::closeness::{closeness, harmonic_centrality, ClosenessOptions};
use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::{
    distances, eccentricity, path_summary, reachability, PathOptions,
};
use engine::metrics::ego::{neighborhood_sizes, NeighborhoodOptions};
use engine::metrics::reach::{reach, ReachOptions};
use engine::{generate, Graph, PathAlgorithm};

/// Random networks of awkward sizes (not multiples of 64, several
/// components, isolates) in both directions and at several densities.
fn cases() -> Vec<Graph> {
    let mut out = Vec::new();
    let mut seed = 1;
    for &(n, m) in &[
        (1usize, 0u64),
        (7, 6),
        (70, 60),
        (150, 300),
        (130, 2000),
        (200, 150),
    ] {
        for directed in [false, true] {
            seed += 1;
            let l = generate::gnm(n, m, directed, false, 1, seed, 1)
                .unwrap()
                .remove(0);
            let from: Vec<i64> = l.from.iter().map(|&v| v as i64).collect();
            let to: Vec<i64> = l.to.iter().map(|&v| v as i64).collect();
            out.push(Graph::from_edges(n, &from, &to, None, directed).unwrap());
        }
    }
    out
}

/// Hop counts from `s` along out-arcs; `u32::MAX` when unreachable.
fn plain_bfs(g: &Graph, s: usize) -> Vec<u32> {
    let adj = g.out();
    let mut d = vec![u32::MAX; g.node_count()];
    d[s] = 0;
    let mut q = vec![s as u32];
    let mut head = 0;
    while head < q.len() {
        let v = q[head] as usize;
        head += 1;
        for &w in adj.neighbors(v) {
            if d[w as usize] == u32::MAX {
                d[w as usize] = d[v] + 1;
                q.push(w);
            }
        }
    }
    d
}

/// The same network with every edge weighted 1.
fn unit_weighted(g: &Graph) -> Graph {
    let from: Vec<i64> = g.edge_sources().iter().map(|&v| v as i64).collect();
    let to: Vec<i64> = g.edge_targets().iter().map(|&v| v as i64).collect();
    let w = vec![1.0; from.len()];
    Graph::from_edges(g.node_count(), &from, &to, Some(&w), g.is_directed()).unwrap()
}

#[test]
fn multi_source_search_equals_plain_search() {
    for g in cases() {
        let n = g.node_count();
        let mut ms = MultiBfs::new(n);
        // Batches of every size up to 64, including repeated sources.
        let mut sources: Vec<u32> = (0..n as u32).rev().collect();
        sources.extend(0..(n as u32).min(5));
        for chunk in sources.chunks(37).chain(sources.chunks(64)) {
            let mut got = vec![vec![u32::MAX; n]; chunk.len()];
            ms.run(g.out(), chunk, u32::MAX, |step| {
                if let Step::Reached {
                    level,
                    node,
                    sources,
                } = step
                {
                    for_each_bit(sources, |i| {
                        assert_eq!(got[i][node as usize], u32::MAX, "reported twice");
                        got[i][node as usize] = level;
                    });
                }
            });
            for (i, &s) in chunk.iter().enumerate() {
                let mut want = plain_bfs(&g, s as usize);
                want[s as usize] = u32::MAX;
                assert_eq!(got[i], want, "n = {n}, source {s}");
            }
        }
    }
}

#[test]
fn distances_and_reach_match_plain_search() {
    for g in cases() {
        let n = g.node_count();
        let d = distances(&g, None, PathOptions::default()).unwrap();
        let reach_counts = reachability(&g, DegreeMode::Out, 1).unwrap();
        let sizes = neighborhood_sizes(
            &g,
            None,
            NeighborhoodOptions {
                order: 2,
                mode: DegreeMode::Out,
                ..Default::default()
            },
        )
        .unwrap();
        let reach2 = reach(
            &g,
            ReachOptions {
                k: 2.0,
                ..Default::default()
            },
        )
        .unwrap();
        for s in 0..n {
            let want = plain_bfs(&g, s);
            for t in 0..n {
                let w = if want[t] == u32::MAX {
                    f64::INFINITY
                } else {
                    want[t] as f64
                };
                assert_eq!(d.get(s, t), w);
            }
            let reached = want.iter().filter(|&&x| x != u32::MAX).count() - 1;
            if g.is_directed() {
                assert_eq!(reach_counts[s], reached as f64);
            }
            let within2 = want.iter().filter(|&&x| x <= 2).count();
            assert_eq!(sizes[s], within2, "neighborhood of {s}");
            assert_eq!(reach2[s], (within2 - 1) as f64, "reach of {s}");
        }
    }
}

#[test]
fn search_kernels_agree_on_unit_lengths() {
    for g in cases() {
        let u = unit_weighted(&g);
        let bfs_c = ClosenessOptions::default();
        let dij_c = ClosenessOptions {
            weighted: true,
            algorithm: PathAlgorithm::Dijkstra,
            ..Default::default()
        };
        let a = closeness(&g, bfs_c).unwrap();
        let b = closeness(&u, dij_c).unwrap();
        for (x, y) in a.iter().zip(&b) {
            assert!(x == y || (x.is_nan() && y.is_nan()), "{x} vs {y}");
        }
        // Harmonic sums are formed per level by the multi-source search and
        // per node by Dijkstra; they agree to rounding.
        let a = harmonic_centrality(&g, bfs_c).unwrap();
        let b = harmonic_centrality(&u, dij_c).unwrap();
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() <= 1e-12 * x.abs().max(1.0), "{x} vs {y}");
        }
        let bfs_p = PathOptions::default();
        let dij_p = PathOptions {
            weighted: true,
            algorithm: PathAlgorithm::Dijkstra,
            ..Default::default()
        };
        assert_eq!(
            eccentricity(&g, bfs_p).unwrap(),
            eccentricity(&u, dij_p).unwrap()
        );
        let (p, q) = (
            path_summary(&g, None, bfs_p).unwrap(),
            path_summary(&u, None, dij_p).unwrap(),
        );
        assert_eq!((p.pairs, p.total, p.longest), (q.pairs, q.total, q.longest));
        let bfs_b = BetweennessOptions::default();
        let dij_b = BetweennessOptions {
            weighted: true,
            algorithm: PathAlgorithm::Dijkstra,
            ..Default::default()
        };
        let (x, y) = (
            betweenness(&g, bfs_b).unwrap(),
            betweenness(&u, dij_b).unwrap(),
        );
        for (a, b) in x.iter().zip(&y) {
            assert!((a - b).abs() <= 1e-9 * a.abs().max(1.0), "{a} vs {b}");
        }
    }
}

#[test]
fn batched_searches_do_not_depend_on_thread_count() {
    for g in cases() {
        let run = |t: usize| {
            let c = ClosenessOptions {
                threads: t,
                ..Default::default()
            };
            let p = PathOptions {
                threads: t,
                ..Default::default()
            };
            let mut v = closeness(&g, c).unwrap();
            v.extend(harmonic_centrality(&g, c).unwrap());
            v.extend(eccentricity(&g, p).unwrap());
            v.extend(distances(&g, None, p).unwrap().as_slice());
            v.extend(reachability(&g, DegreeMode::Out, t).unwrap());
            v.into_iter().map(f64::to_bits).collect::<Vec<u64>>()
        };
        let one = run(1);
        assert_eq!(one, run(2));
        assert_eq!(one, run(3));
    }
}
