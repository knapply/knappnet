//! Randomized property tests for key players, group centrality, and flow
//! betweenness (`metrics::resilience`).

extern crate knappnet as engine;

use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::resilience::{
    flow_betweenness, group_centrality, key_players, kpp_score, FlowScale, GroupMeasure,
    GroupOptions, KppMeasure, KppMethod, KppOptions,
};
use engine::metrics::topography::{fragmentation, ComponentKind};
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
    weights: Vec<f64>,
    directed: bool,
}

fn random_case(rng: &mut Rng, max_n: u64) -> Case {
    let n = 2 + rng.below(max_n - 1) as usize;
    let m = rng.below((n * n / 2 + 1) as u64) as usize;
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    let mut weights = Vec::with_capacity(m);
    for _ in 0..m {
        src.push(rng.below(n as u64) as i64);
        tgt.push(rng.below(n as u64) as i64);
        weights.push(1.0 + rng.below(4) as f64);
    }
    Case {
        n,
        src,
        tgt,
        weights,
        directed: rng.below(2) == 0,
    }
}

impl Case {
    fn graph(&self, weighted: bool) -> Graph {
        let w = weighted.then_some(self.weights.as_slice());
        Graph::from_edges(self.n, &self.src, &self.tgt, w, self.directed).unwrap()
    }

    /// The same network with node v renamed perm[v].
    fn relabeled(&self, perm: &[usize]) -> Graph {
        let s: Vec<i64> = self.src.iter().map(|&v| perm[v as usize] as i64).collect();
        let t: Vec<i64> = self.tgt.iter().map(|&v| perm[v as usize] as i64).collect();
        Graph::from_edges(self.n, &s, &t, None, self.directed).unwrap()
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
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

const CASES: u64 = 150;

#[test]
fn kpp_scores_permute_with_node_labels() {
    let mut rng = Rng(1);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 14);
        let g = c.graph(false);
        let perm = permutation(&mut rng, c.n);
        let h = c.relabeled(&perm);
        for measure in [
            KppMeasure::Fragmentation,
            KppMeasure::DistanceFragmentation,
            KppMeasure::MReach,
            KppMeasure::DistanceReach,
        ] {
            let opts = KppOptions {
                measure,
                m: 2.0,
                ..Default::default()
            };
            let set = [0u32, (c.n - 1) as u32];
            let set = if c.n > 2 { &set[..] } else { &set[..1] };
            let mapped: Vec<u32> = set.iter().map(|&v| perm[v as usize] as u32).collect();
            let a = kpp_score(&g, set, opts).unwrap();
            let b = kpp_score(&h, &mapped, opts).unwrap();
            assert!(close(a, b), "{measure:?}: {a} vs {b}");
        }
    }
}

#[test]
fn fragmentation_of_the_empty_set_matches_components() {
    let mut rng = Rng(2);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 20);
        let g = c.graph(false);
        let f = kpp_score(&g, &[], KppOptions::default()).unwrap();
        if !c.directed {
            assert!(close(f, fragmentation(&g, ComponentKind::Weak)));
        } else {
            // Ordered pairs that cannot reach each other lie between the
            // shares of pairs in different strong and weak components.
            assert!(f <= fragmentation(&g, ComponentKind::Strong) + 1e-12);
            assert!(f >= fragmentation(&g, ComponentKind::Weak) - 1e-12);
        }
    }
}

#[test]
fn greedy_never_beats_exact_and_ignores_threads() {
    let mut rng = Rng(3);
    for case in 0..60 {
        let c = random_case(&mut rng, 12);
        let g = c.graph(false);
        let k = 1 + (case % 2).min(c.n - 1);
        for measure in [KppMeasure::DistanceFragmentation, KppMeasure::DistanceReach] {
            let opts = KppOptions {
                measure,
                restarts: 3,
                seed: case as u64,
                ..Default::default()
            };
            let one = key_players(&g, k, KppOptions { threads: 1, ..opts }).unwrap();
            let many = key_players(&g, k, KppOptions { threads: 4, ..opts }).unwrap();
            // Compared field by field: NaN scores (too few nodes left) are
            // not equal to themselves.
            assert_eq!(one.set, many.set);
            assert_eq!(one.evaluations, many.evaluations);
            assert!(close(one.score, many.score));
            let exact = key_players(
                &g,
                k,
                KppOptions {
                    method: KppMethod::Exact,
                    threads: 3,
                    ..opts
                },
            )
            .unwrap();
            if !exact.score.is_nan() {
                assert!(one.score <= exact.score + 1e-12);
            }
            let rescored = kpp_score(&g, &one.set, opts).unwrap();
            assert!(close(rescored, one.score));
        }
    }
}

#[test]
fn group_betweenness_of_one_node_is_its_betweenness() {
    let mut rng = Rng(4);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 16);
        for weighted in [false, true] {
            let g = c.graph(weighted);
            let b = betweenness(
                &g,
                BetweennessOptions {
                    weighted,
                    ..Default::default()
                },
            )
            .unwrap();
            let v = (c.n - 1) as u32;
            let gb = group_centrality(
                &g,
                &[v],
                GroupMeasure::Betweenness,
                GroupOptions {
                    weighted,
                    normalized: false,
                    threads: 2,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(close(gb, b[v as usize]), "{gb} vs {}", b[v as usize]);
        }
    }
}

#[test]
fn flow_betweenness_is_consistent() {
    let mut rng = Rng(5);
    for _ in 0..60 {
        let c = random_case(&mut rng, 10);
        let g = c.graph(false);
        let raw = flow_betweenness(&g, false, FlowScale::Raw, 1).unwrap();
        let raw4 = flow_betweenness(&g, false, FlowScale::Raw, 4).unwrap();
        assert_eq!(raw, raw4);
        let frac = flow_betweenness(&g, false, FlowScale::Fractional, 2).unwrap();
        let norm = flow_betweenness(&g, false, FlowScale::Normalized, 2).unwrap();
        for v in 0..c.n {
            // Unit capacities give integer flows, and each pair's share is
            // at most 1 and at most its loss.
            assert_eq!(raw[v], raw[v].round());
            assert!(frac[v] <= raw[v] + 1e-12);
            assert!(norm[v].is_nan() || (0.0..=1.0 + 1e-12).contains(&norm[v]));
        }
        let perm = permutation(&mut rng, c.n);
        let h = c.relabeled(&perm);
        let moved = flow_betweenness(&h, false, FlowScale::Raw, 2).unwrap();
        for v in 0..c.n {
            assert!(close(raw[v], moved[perm[v]]));
        }
        // Doubling every capacity doubles every flow.
        let w2 = vec![2.0; c.src.len()];
        let g2 = Graph::from_edges(c.n, &c.src, &c.tgt, Some(&w2), c.directed).unwrap();
        let doubled = flow_betweenness(&g2, true, FlowScale::Raw, 2).unwrap();
        for v in 0..c.n {
            assert!(close(doubled[v], 2.0 * raw[v]));
        }
    }
}
