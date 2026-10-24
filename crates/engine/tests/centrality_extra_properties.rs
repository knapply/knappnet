//! Randomized property tests for the centrality variants and distance
//! measures: results do not depend on the thread count, relabeling nodes
//! permutes them, special cases reduce to known measures (sources and
//! targets equal to every node give ordinary betweenness; on trees load
//! equals twice betweenness; resistance distance on a tree is the path
//! length; the approximation over every ordered pair is exact current-flow
//! betweenness).

#![allow(clippy::needless_range_loop)]

extern crate knappnet as engine;

use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::communicability::{communicability, communicability_betweenness};
use engine::metrics::current_flow::{
    current_flow_betweenness, current_flow_betweenness_approx, current_flow_betweenness_subset,
    current_flow_closeness, kemeny_constant, resistance_distance, second_order_centrality,
    FlowOptions, Sample,
};
use engine::metrics::descriptive::PathOptions;
use engine::metrics::hierarchy::{local_reaching_centrality, trophic_levels, ReachingOptions};
use engine::metrics::path_centrality::{
    betweenness_subset, edge_load, load_centrality, percolation_centrality, LoadOptions,
    PercolationOptions, SubsetOptions,
};
use engine::metrics::spreading::{laplacian_centrality, voterank};
use engine::Graph;

/// SplitMix64 (Steele, Lea, and Flood 2014, `steele2014fast`).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    w: Vec<f64>,
    directed: bool,
}

fn random_case(rng: &mut Rng, directed: bool) -> Case {
    let n = 2 + rng.below(18) as usize;
    let m = rng.below((n * n / 3 + 2) as u64) as usize;
    Case {
        n,
        src: (0..m).map(|_| rng.below(n as u64) as i64).collect(),
        tgt: (0..m).map(|_| rng.below(n as u64) as i64).collect(),
        w: (0..m).map(|_| (1 + rng.below(8)) as f64 / 2.0).collect(),
        directed,
    }
}

/// A random connected undirected network: a random tree plus extra edges.
fn connected_case(rng: &mut Rng) -> Case {
    let n = 3 + rng.below(14) as usize;
    let mut src = Vec::new();
    let mut tgt = Vec::new();
    for v in 1..n {
        src.push(rng.below(v as u64) as i64);
        tgt.push(v as i64);
    }
    for _ in 0..rng.below(n as u64) {
        src.push(rng.below(n as u64) as i64);
        tgt.push(rng.below(n as u64) as i64);
    }
    let m = src.len();
    Case {
        n,
        w: (0..m).map(|_| (1 + rng.below(8)) as f64 / 2.0).collect(),
        src,
        tgt,
        directed: false,
    }
}

fn graph(c: &Case) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, Some(&c.w), c.directed).unwrap()
}

fn relabeled(c: &Case, p: &[usize]) -> Case {
    Case {
        n: c.n,
        src: c.src.iter().map(|&v| p[v as usize] as i64).collect(),
        tgt: c.tgt.iter().map(|&v| p[v as usize] as i64).collect(),
        w: c.w.clone(),
        directed: c.directed,
    }
}

fn permutation(rng: &mut Rng, n: usize) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        p.swap(i, rng.below(i as u64 + 1) as usize);
    }
    p
}

fn close(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || a == b || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

#[test]
fn path_variants_are_thread_independent_and_permute() {
    let mut rng = Rng(5);
    for i in 0..60 {
        let c = random_case(&mut rng, i % 2 == 0);
        let g = graph(&c);
        let p = permutation(&mut rng, c.n);
        let h = graph(&relabeled(&c, &p));
        let opts = |threads| SubsetOptions {
            weighted: true,
            threads,
            ..Default::default()
        };
        let one = betweenness_subset(&g, opts(1)).unwrap();
        assert_eq!(one, betweenness_subset(&g, opts(4)).unwrap());
        let other = betweenness_subset(&h, opts(1)).unwrap();
        for v in 0..c.n {
            assert!(close(one.nodes[v], other.nodes[p[v]]));
        }
        assert_eq!(one.edges.len(), other.edges.len());
        for e in 0..c.src.len() {
            assert!(close(one.edges[e], other.edges[e]));
        }
        let bc = betweenness(
            &g,
            BetweennessOptions {
                weighted: true,
                ..Default::default()
            },
        )
        .unwrap();
        for v in 0..c.n {
            assert!(close(one.nodes[v], bc[v]));
        }
        let states: Vec<f64> = (0..c.n).map(|v| (v % 4) as f64 / 3.0).collect();
        let pc = |t| {
            percolation_centrality(
                &g,
                &states,
                PercolationOptions {
                    threads: t,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        assert_eq!(bits(pc(1)), bits(pc(3)));
        let lc = |t| {
            load_centrality(
                &g,
                LoadOptions {
                    threads: t,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        assert_eq!(lc(1), lc(5));
        assert_eq!(
            edge_load(&g, f64::INFINITY, 1).unwrap(),
            edge_load(&g, f64::INFINITY, 4).unwrap()
        );
    }
}

#[test]
fn load_on_trees_is_twice_betweenness() {
    // Shortest paths are unique in a tree, so splitting equally and in
    // proportion agree; load sums over ordered pairs.
    let mut rng = Rng(9);
    for _ in 0..40 {
        let n = 2 + rng.below(20) as usize;
        let src: Vec<i64> = (1..n).map(|v| rng.below(v as u64) as i64).collect();
        let tgt: Vec<i64> = (1..n as i64).collect();
        let g = Graph::from_edges(n, &src, &tgt, None, false).unwrap();
        let load = load_centrality(&g, LoadOptions::default()).unwrap();
        let bc = betweenness(&g, BetweennessOptions::default()).unwrap();
        for v in 0..n {
            assert!(close(load[v], 2.0 * bc[v]));
        }
        // Resistance distance on a tree is the path length.
        let r = resistance_distance(&g, FlowOptions::default()).unwrap();
        let d = engine::metrics::descriptive::distances(&g, None, PathOptions::default()).unwrap();
        for a in 0..n {
            for b in 0..n {
                assert!(close(r.get(a, b), d.get(a, b)));
            }
        }
    }
}

#[test]
fn current_flow_identities() {
    let mut rng = Rng(13);
    for _ in 0..40 {
        let c = connected_case(&mut rng);
        let g = graph(&c);
        let n = c.n;
        let opts = FlowOptions {
            weighted: true,
            threads: 1,
            ..Default::default()
        };
        let full = current_flow_betweenness(&g, opts).unwrap();
        let many = current_flow_betweenness(&g, FlowOptions { threads: 4, ..opts }).unwrap();
        assert_eq!(full, many);
        let all: Vec<u32> = (0..n as u32).collect();
        let sub = current_flow_betweenness_subset(&g, &all, &all, false, opts).unwrap();
        for v in 0..n {
            assert!(close(sub.nodes[v], full.nodes[v]));
        }
        for e in 0..c.src.len() {
            assert!(close(sub.edges[e], full.edges[e]));
        }
        let pairs: Vec<(u32, u32)> = (0..n as u32)
            .flat_map(|s| (0..n as u32).filter(move |&t| t != s).map(move |t| (s, t)))
            .collect();
        let approx = current_flow_betweenness_approx(&g, Sample::Pairs(&pairs), opts).unwrap();
        for v in 0..n {
            assert!(close(approx[v], full.nodes[v]));
        }
        // Closeness is 1 / sum of resistances, and resistances are a metric.
        let r = resistance_distance(&g, opts).unwrap();
        let cl = current_flow_closeness(&g, opts).unwrap();
        for a in 0..n {
            let s: f64 = (0..n).map(|b| r.get(a, b)).sum();
            assert!(close(cl[a], 1.0 / s));
            for b in 0..n {
                assert!(close(r.get(a, b), r.get(b, a)));
                for x in 0..n {
                    assert!(r.get(a, b) <= r.get(a, x) + r.get(x, b) + 1e-9);
                }
            }
        }
        // Relabeling permutes node values and leaves Kemeny's constant.
        let p = permutation(&mut rng, n);
        let h = graph(&relabeled(&c, &p));
        let other = current_flow_betweenness(&h, opts).unwrap();
        for v in 0..n {
            assert!(close(full.nodes[v], other.nodes[p[v]]));
        }
        assert!(close(
            kemeny_constant(&g, true).unwrap(),
            kemeny_constant(&h, true).unwrap()
        ));
        let so = second_order_centrality(&g, true).unwrap();
        let so2 = second_order_centrality(&h, true).unwrap();
        for v in 0..n {
            assert!(close(so[v], so2[p[v]]));
        }
    }
}

#[test]
fn communicability_is_symmetric_and_permutes() {
    let mut rng = Rng(17);
    for _ in 0..30 {
        let c = random_case(&mut rng, false);
        let g = graph(&c);
        let n = c.n;
        let m = communicability(&g).unwrap();
        for a in 0..n {
            for b in 0..n {
                assert!(close(m.get(a, b), m.get(b, a)));
                assert!(m.get(a, b) >= -1e-12);
            }
        }
        let p = permutation(&mut rng, n);
        let h = graph(&relabeled(&c, &p));
        let one = communicability_betweenness(&g, 1).unwrap();
        assert_eq!(one, communicability_betweenness(&g, 3).unwrap());
        let other = communicability_betweenness(&h, 1).unwrap();
        for v in 0..n {
            assert!(close(one[v], other[p[v]]));
        }
    }
}

#[test]
fn local_measures_permute() {
    let mut rng = Rng(21);
    for i in 0..60 {
        let directed = i % 2 == 0;
        let c = random_case(&mut rng, directed);
        let g = graph(&c);
        let n = c.n;
        let p = permutation(&mut rng, n);
        let h = graph(&relabeled(&c, &p));
        let a = laplacian_centrality(&g, true, false).unwrap();
        let b = laplacian_centrality(&h, true, false).unwrap();
        for v in 0..n {
            assert!(close(a[v], b[p[v]]));
        }
        let ro = ReachingOptions::default();
        let a = local_reaching_centrality(&g, ro).unwrap();
        let b = local_reaching_centrality(&h, ro).unwrap();
        for v in 0..n {
            assert!(close(a[v], b[p[v]]));
        }
        if directed {
            let a = trophic_levels(&g, true).unwrap();
            let b = trophic_levels(&h, true).unwrap();
            for v in 0..n {
                assert!(close(a[v], b[p[v]]));
            }
        }
        // VoteRank never picks a node twice.
        let mut chosen = voterank(&g, None);
        let k = chosen.len();
        chosen.sort_unstable();
        chosen.dedup();
        assert_eq!(chosen.len(), k);
    }
}
