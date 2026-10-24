//! Randomized property tests for the shortest-path module
//! (`engine::paths`): every algorithm agrees with every other on the
//! lengths they share, relabeling nodes permutes the results, unit weights
//! reproduce unweighted results, strengths w give the results of distances
//! 1/w, the distances agree with the engine's distance matrix
//! (`metrics::descriptive::distances`), and nothing depends on the thread
//! count.

#![allow(clippy::needless_range_loop)]

extern crate knappnet as engine;

use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::{distances, PathOptions};
use engine::paths::bidirectional::bidirectional;
use engine::paths::dense::floyd_warshall;
use engine::paths::enumerate::{all_shortest_paths, count_shortest_paths};
use engine::paths::negative::{bellman_ford, goldberg_radzik, johnson, Outcome};
use engine::paths::search::{all_pairs, single_source, SearchOptions};
use engine::paths::simple::k_shortest_simple_paths;
use engine::paths::{ArcLengths, PathGraph, PathGraphOptions, Tie, NO_NODE};
use engine::{Graph, WeightKind};

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

fn random_case(rng: &mut Rng) -> Case {
    let n = 1 + rng.below(30) as usize;
    let m = rng.below((n * n / 3 + 1) as u64) as usize;
    let src = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    let tgt = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    // Quarter-integers, so sums are exact and ties are real ties.
    let w = (0..m).map(|_| (1 + rng.below(16)) as f64 / 4.0).collect();
    Case {
        n,
        src,
        tgt,
        w,
        directed: rng.below(2) == 0,
    }
}

fn graph(c: &Case) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, Some(&c.w), c.directed).unwrap()
}

fn path_graph(g: &Graph, mode: DegreeMode, lengths: ArcLengths) -> PathGraph {
    PathGraph::new(
        g,
        PathGraphOptions {
            mode,
            lengths,
            collapse_parallel: false,
        },
    )
    .unwrap()
}

const WEIGHTS: ArcLengths<'static> = ArcLengths::Weights(WeightKind::Distance);

#[test]
fn every_algorithm_gives_the_same_distances() {
    let mut rng = Rng(11);
    for _ in 0..120 {
        let c = random_case(&mut rng);
        let g = graph(&c);
        let n = c.n;
        for mode in [DegreeMode::Out, DegreeMode::In, DegreeMode::All] {
            for (lengths, weighted) in [(ArcLengths::Unit, false), (WEIGHTS, true)] {
                let pg = path_graph(&g, mode, lengths);
                let (d, _) = all_pairs(&pg, None, f64::INFINITY, 3).unwrap();
                // The engine's distance matrix.
                let reference = distances(
                    &g,
                    None,
                    PathOptions {
                        mode,
                        weighted,
                        threads: 1,
                        ..Default::default()
                    },
                )
                .unwrap();
                let (dj, _) = johnson(&pg, None, 2).unwrap().unwrap();
                let fw = floyd_warshall(&pg, 2);
                assert!(!fw.negative_cycle);
                for s in 0..n {
                    let bf = bellman_ford(&pg, &[s as u32]).unwrap().tree().unwrap();
                    let gr = goldberg_radzik(&pg, s as u32).unwrap().tree().unwrap();
                    for t in 0..n {
                        let want = reference.get(s, t);
                        assert_eq!(d[s * n + t], want);
                        assert_eq!(dj[s * n + t], want);
                        assert_eq!(fw.dist[s * n + t], want);
                        assert_eq!(bf.dist[t], want);
                        assert_eq!(gr.dist[t], want);
                        let bd = bidirectional(&pg, s as u32, t as u32).unwrap();
                        assert_eq!(bd.map_or(f64::INFINITY, |x| x.0), want);
                    }
                }
            }
        }
    }
}

#[test]
fn relabeling_permutes_distances_and_counts() {
    let mut rng = Rng(23);
    for _ in 0..80 {
        let c = random_case(&mut rng);
        let n = c.n;
        // A random permutation p: node v becomes p[v].
        let mut p: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            p.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let relabeled = Case {
            n,
            src: c.src.iter().map(|&v| p[v as usize] as i64).collect(),
            tgt: c.tgt.iter().map(|&v| p[v as usize] as i64).collect(),
            w: c.w.clone(),
            directed: c.directed,
        };
        let (g, h) = (graph(&c), graph(&relabeled));
        for lengths in [ArcLengths::Unit, WEIGHTS] {
            let (pa, pb) = (
                path_graph(&g, DegreeMode::Out, lengths),
                path_graph(&h, DegreeMode::Out, lengths),
            );
            let (da, _) = all_pairs(&pa, None, f64::INFINITY, 1).unwrap();
            let (db, _) = all_pairs(&pb, None, f64::INFINITY, 1).unwrap();
            for s in 0..n {
                let ca = count_shortest_paths(&pa, s as u32, Tie::Exact).unwrap();
                let cb = count_shortest_paths(&pb, p[s] as u32, Tie::Exact).unwrap();
                for t in 0..n {
                    assert_eq!(da[s * n + t], db[p[s] * n + p[t]]);
                    assert_eq!(ca[t], cb[p[t]]);
                }
            }
        }
    }
}

#[test]
fn unit_weights_and_strengths() {
    let mut rng = Rng(37);
    for _ in 0..80 {
        let mut c = random_case(&mut rng);
        let n = c.n;
        let g = graph(&c);
        // Strengths w give the distances of lengths 1/w.
        let inv: Vec<f64> = c.w.iter().map(|w| 1.0 / w).collect();
        let a = path_graph(
            &g,
            DegreeMode::Out,
            ArcLengths::Weights(WeightKind::Strength),
        );
        let b = path_graph(&g, DegreeMode::Out, ArcLengths::PerEdge(&inv));
        assert_eq!(
            all_pairs(&a, None, f64::INFINITY, 1).unwrap(),
            all_pairs(&b, None, f64::INFINITY, 1).unwrap()
        );
        // Unit weights reproduce unweighted searches.
        c.w = vec![1.0; c.w.len()];
        let g = graph(&c);
        let a = path_graph(&g, DegreeMode::Out, WEIGHTS);
        let b = path_graph(&g, DegreeMode::Out, ArcLengths::Unit);
        let (da, pa) = all_pairs(&a, None, f64::INFINITY, 1).unwrap();
        let (db, pb) = all_pairs(&b, None, f64::INFINITY, 1).unwrap();
        assert_eq!(da, db);
        // Parents may differ between the two searches (both are shortest).
        for s in 0..n {
            for t in 0..n {
                assert_eq!(pa[s * n + t] == NO_NODE, pb[s * n + t] == NO_NODE);
            }
        }
    }
}

#[test]
fn thread_counts_do_not_matter() {
    let mut rng = Rng(41);
    for _ in 0..30 {
        let c = random_case(&mut rng);
        let g = graph(&c);
        let pg = path_graph(&g, DegreeMode::Out, WEIGHTS);
        let one = all_pairs(&pg, None, f64::INFINITY, 1).unwrap();
        let many = all_pairs(&pg, None, f64::INFINITY, 7).unwrap();
        assert_eq!(one, many);
        assert_eq!(
            johnson(&pg, None, 1).unwrap(),
            johnson(&pg, None, 5).unwrap()
        );
        assert_eq!(floyd_warshall(&pg, 1), floyd_warshall(&pg, 6));
    }
}

#[test]
fn shortest_path_lists_count_and_rank() {
    let mut rng = Rng(53);
    for _ in 0..80 {
        let c = random_case(&mut rng);
        let n = c.n;
        let g = graph(&c);
        for lengths in [ArcLengths::Unit, WEIGHTS] {
            let pg = path_graph(&g, DegreeMode::Out, lengths);
            let sigma = count_shortest_paths(&pg, 0, Tie::Exact).unwrap();
            let tree = single_source(&pg, 0, SearchOptions::default()).unwrap();
            for t in 0..n {
                let all = all_shortest_paths(&pg, 0, t as u32, Tie::Exact, 1 << 20).unwrap();
                assert_eq!(all.len() as f64, sigma[t]);
                for i in 0..all.len() {
                    assert_eq!(pg.path_length(all.get(i)), Some(tree.dist[t]));
                }
                // The first of the k shortest simple paths is a shortest path.
                let k = k_shortest_simple_paths(&pg, 0, t as u32, 3).unwrap();
                if tree.dist[t].is_finite() {
                    assert_eq!(k.values[0], tree.dist[t]);
                    assert!(k.values.windows(2).all(|w| w[0] <= w[1]));
                } else {
                    assert!(k.is_empty());
                }
            }
        }
    }
}

#[test]
fn negative_lengths_from_potentials_keep_shortest_paths() {
    // Reweighting by a potential p (w + p(u) - p(v)) changes every path
    // from s to t by p(s) - p(t) and creates no negative cycle.
    let mut rng = Rng(67);
    for _ in 0..80 {
        let mut c = random_case(&mut rng);
        c.directed = true;
        let n = c.n;
        let g = graph(&c);
        let pot: Vec<f64> = (0..n).map(|_| rng.below(9) as f64).collect();
        let lens: Vec<f64> = (0..c.src.len())
            .map(|e| c.w[e] + pot[c.src[e] as usize] - pot[c.tgt[e] as usize])
            .collect();
        let base = path_graph(&g, DegreeMode::Out, WEIGHTS);
        let neg = path_graph(&g, DegreeMode::Out, ArcLengths::PerEdge(&lens));
        let (d, _) = all_pairs(&base, None, f64::INFINITY, 1).unwrap();
        let (dj, _) = johnson(&neg, None, 1).unwrap().unwrap();
        let fw = floyd_warshall(&neg, 1);
        for s in 0..n {
            let bf = match bellman_ford(&neg, &[s as u32]).unwrap() {
                Outcome::Tree(t) => t,
                Outcome::NegativeCycle(cy) => panic!("unexpected cycle {cy:?}"),
            };
            for t in 0..n {
                let want = if d[s * n + t].is_finite() {
                    d[s * n + t] + pot[s] - pot[t]
                } else {
                    f64::INFINITY
                };
                assert_eq!(bf.dist[t], want);
                assert_eq!(dj[s * n + t], want);
                assert_eq!(fw.dist[s * n + t], want);
            }
        }
    }
}
