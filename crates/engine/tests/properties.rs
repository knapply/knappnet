//! Randomized property tests. These check invariants that must hold for
//! every graph, complementing the oracle comparisons in `validation/`.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::degree::{degree, DegreeMode};
use engine::metrics::triad::triad_census;
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
    directed: bool,
}

fn random_case(rng: &mut Rng) -> Case {
    let n = 1 + rng.below(40) as usize;
    let m = rng.below((n * n / 2 + 1) as u64) as usize;
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    for _ in 0..m {
        src.push(rng.below(n as u64) as i64);
        tgt.push(rng.below(n as u64) as i64);
    }
    Case {
        n,
        src,
        tgt,
        directed: rng.below(2) == 0,
    }
}

fn close(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| (x - y).abs() <= 1e-9 * (1.0 + x.abs().max(y.abs())))
}

const CASES: u64 = 300;

#[test]
fn relabeling_nodes_permutes_results() {
    let mut rng = Rng(1);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        // Random permutation via Fisher-Yates.
        let mut perm: Vec<usize> = (0..c.n).collect();
        for i in (1..c.n).rev() {
            perm.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let ps: Vec<i64> = c.src.iter().map(|&v| perm[v as usize] as i64).collect();
        let pt: Vec<i64> = c.tgt.iter().map(|&v| perm[v as usize] as i64).collect();
        let g = Graph::from_edges(c.n, &c.src, &c.tgt, None, c.directed).unwrap();
        let h = Graph::from_edges(c.n, &ps, &pt, None, c.directed).unwrap();

        let bg = betweenness(&g, BetweennessOptions::default()).unwrap();
        let bh = betweenness(&h, BetweennessOptions::default()).unwrap();
        let mapped: Vec<f64> = (0..c.n).map(|v| bh[perm[v]]).collect();
        assert!(close(&bg, &mapped));

        let dg = degree(&g, DegreeMode::All, true);
        let dh = degree(&h, DegreeMode::All, true);
        assert!((0..c.n).all(|v| dg[v] == dh[perm[v]]));

        assert_eq!(triad_census(&g, 1), triad_census(&h, 1));
    }
}

#[test]
fn unit_weights_match_unweighted() {
    let mut rng = Rng(2);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let w = vec![1.0; c.src.len()];
        let g = Graph::from_edges(c.n, &c.src, &c.tgt, Some(&w), c.directed).unwrap();
        let a = betweenness(&g, BetweennessOptions::default()).unwrap();
        let b = betweenness(
            &g,
            BetweennessOptions {
                weighted: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(close(&a, &b), "{a:?} vs {b:?}");
    }
}

#[test]
fn results_do_not_depend_on_thread_count() {
    let mut rng = Rng(3);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = Graph::from_edges(c.n, &c.src, &c.tgt, None, c.directed).unwrap();
        let one = betweenness(
            &g,
            BetweennessOptions {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let many = betweenness(
            &g,
            BetweennessOptions {
                threads: 4,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(close(&one, &many));
        assert_eq!(triad_census(&g, 1), triad_census(&g, 4));
    }
}

#[test]
fn degree_and_census_totals() {
    let mut rng = Rng(4);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = Graph::from_edges(c.n, &c.src, &c.tgt, None, c.directed).unwrap();
        let total: f64 = degree(&g, DegreeMode::All, true).iter().sum();
        assert_eq!(total, 2.0 * c.src.len() as f64);
        let n = c.n as f64;
        let triples = if c.n < 3 {
            0.0
        } else {
            n * (n - 1.0) * (n - 2.0) / 6.0
        };
        assert_eq!(triad_census(&g, 1).iter().sum::<f64>(), triples);
    }
}

#[test]
fn undirected_betweenness_matches_symmetric_directed() {
    // An undirected edge is a pair of opposite arcs; undirected betweenness
    // counts each unordered pair once, so it is half the directed value.
    let mut rng = Rng(5);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = Graph::from_edges(c.n, &c.src, &c.tgt, None, false).unwrap();
        let (mut s2, mut t2) = (c.src.clone(), c.tgt.clone());
        for (s, t) in c.src.iter().zip(&c.tgt) {
            if s != t {
                s2.push(*t);
                t2.push(*s);
            }
        }
        let d = Graph::from_edges(c.n, &s2, &t2, None, true).unwrap();
        let bu = betweenness(&g, BetweennessOptions::default()).unwrap();
        let bd: Vec<f64> = betweenness(&d, BetweennessOptions::default())
            .unwrap()
            .into_iter()
            .map(|x| x / 2.0)
            .collect();
        assert!(close(&bu, &bd), "{bu:?} vs {bd:?}");
    }
}
