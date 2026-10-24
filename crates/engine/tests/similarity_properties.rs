//! Randomized property tests for Weisfeiler-Lehman hashes and colours,
//! SimRank, Panther, and graph edit distance: relabeling, the thread
//! count, symmetries, and brute-force edit distances on small networks.

extern crate knappnet as engine;

use engine::metrics::edit_distance::{graph_edit_distance, EditOptions};
use engine::metrics::similarity::{
    panther_similarity, simrank, PantherOptions, SimRankOptions, SimRankStop,
};
use engine::metrics::wl::{wl_colors, wl_graph_hash, wl_subgraph_hashes, Refiner, WlOptions};
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

fn random_case(rng: &mut Rng, max_n: u64) -> Case {
    let n = 1 + rng.below(max_n) as usize;
    let m = rng.below((n * n / 2 + 1) as u64) as usize;
    Case {
        n,
        src: (0..m).map(|_| rng.below(n as u64) as i64).collect(),
        tgt: (0..m).map(|_| rng.below(n as u64) as i64).collect(),
        directed: rng.below(2) == 0,
    }
}

fn permuted(c: &Case, rng: &mut Rng) -> (Vec<usize>, Graph) {
    let mut perm: Vec<usize> = (0..c.n).collect();
    for i in (1..c.n).rev() {
        perm.swap(i, rng.below(i as u64 + 1) as usize);
    }
    // Edges in a shuffled order too.
    let mut order: Vec<usize> = (0..c.src.len()).collect();
    for i in (1..order.len()).rev() {
        order.swap(i, rng.below(i as u64 + 1) as usize);
    }
    let s: Vec<i64> = order
        .iter()
        .map(|&e| perm[c.src[e] as usize] as i64)
        .collect();
    let t: Vec<i64> = order
        .iter()
        .map(|&e| perm[c.tgt[e] as usize] as i64)
        .collect();
    (
        perm,
        Graph::from_edges(c.n, &s, &t, None, c.directed).unwrap(),
    )
}

fn graph(c: &Case) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, None, c.directed).unwrap()
}

#[test]
fn wl_hashes_are_invariant_and_thread_independent() {
    let mut rng = Rng(3);
    for _ in 0..40 {
        let c = random_case(&mut rng, 25);
        let g = graph(&c);
        let (perm, h) = permuted(&c, &mut rng);
        let o = WlOptions {
            threads: 1,
            ..Default::default()
        };
        let o4 = WlOptions { threads: 4, ..o };
        assert_eq!(wl_graph_hash(&g, o).unwrap(), wl_graph_hash(&h, o).unwrap());
        assert_eq!(
            wl_graph_hash(&g, o).unwrap(),
            wl_graph_hash(&g, o4).unwrap()
        );
        let a = wl_subgraph_hashes(&g, o).unwrap();
        let b = wl_subgraph_hashes(&h, o4).unwrap();
        for v in 0..c.n {
            assert_eq!(a[v], b[perm[v]]);
        }
        // Colour classes are a relabeling-invariant partition, and each
        // round refines the one before.
        let mut r = Refiner::new();
        let ca = wl_colors(&g, None, 3, &mut r).unwrap();
        let mut r = Refiner::new();
        let cb = wl_colors(&h, None, 3, &mut r).unwrap();
        for round in 0..3 {
            for u in 0..c.n {
                for v in 0..c.n {
                    let same_a = ca[round][u] == ca[round][v];
                    assert_eq!(same_a, cb[round][perm[u]] == cb[round][perm[v]]);
                    if round > 0 && same_a {
                        assert_eq!(ca[round - 1][u], ca[round - 1][v]);
                    }
                }
            }
        }
    }
}

#[test]
fn simrank_is_symmetric_bounded_and_invariant() {
    let mut rng = Rng(9);
    for _ in 0..30 {
        let c = random_case(&mut rng, 20);
        let g = graph(&c);
        let (perm, h) = permuted(&c, &mut rng);
        for stop in [SimRankStop::Allclose, SimRankStop::Relative] {
            let o = SimRankOptions {
                stop,
                threads: 1,
                ..Default::default()
            };
            let a = simrank(&g, o).unwrap();
            let b = simrank(&h, o).unwrap();
            let t = simrank(&g, SimRankOptions { threads: 3, ..o }).unwrap();
            for u in 0..c.n {
                assert_eq!(a.get(u, u), 1.0);
                for v in 0..c.n {
                    let x = a.get(u, v);
                    assert!((0.0..=1.0).contains(&x));
                    assert!((x - a.get(v, u)).abs() < 1e-12);
                    assert!((x - b.get(perm[u], perm[v])).abs() < 1e-12);
                    assert_eq!(x.to_bits(), t.get(u, v).to_bits());
                }
            }
        }
    }
}

#[test]
fn panther_is_reproducible_and_well_formed() {
    let mut rng = Rng(21);
    let mut checked = 0;
    while checked < 20 {
        let mut c = random_case(&mut rng, 20);
        c.directed = false;
        let g = graph(&c);
        let Some(source) = (0..c.n).find(|&v| g.out().degree(v) > 0) else {
            continue;
        };
        let kept = (0..c.n).filter(|&v| g.out().degree(v) > 0).count();
        let o = PantherOptions {
            k: kept.min(3),
            eps: Some(0.1),
            seed: 5,
            threads: 1,
            ..Default::default()
        };
        let a = panther_similarity(&g, source, o).unwrap();
        let b = panther_similarity(&g, source, PantherOptions { threads: 4, ..o }).unwrap();
        assert_eq!(a, b);
        assert!(!a.nodes.contains(&(source as u32)));
        assert!(a.scores.windows(2).all(|w| w[0] >= w[1]));
        assert!(a.scores.iter().all(|&x| (0.0..=1.0).contains(&x)));
        let other = panther_similarity(&g, source, PantherOptions { seed: 6, ..o }).unwrap();
        assert_eq!(other.nodes.len(), a.nodes.len());
        checked += 1;
    }
}

/// Edit distance by trying every partial injective mapping.
fn brute_ged(a: &Case, b: &Case) -> f64 {
    let key = |x: usize, y: usize, d: bool| if d || x <= y { (x, y) } else { (y, x) };
    let e1: std::collections::BTreeSet<(usize, usize)> = a
        .src
        .iter()
        .zip(&a.tgt)
        .map(|(&s, &t)| key(s as usize, t as usize, a.directed))
        .collect();
    let e2: std::collections::BTreeSet<(usize, usize)> = b
        .src
        .iter()
        .zip(&b.tgt)
        .map(|(&s, &t)| key(s as usize, t as usize, b.directed))
        .collect();
    fn rec(
        n1: usize,
        f: &mut Vec<Option<usize>>,
        used: &mut Vec<bool>,
        best: &mut f64,
        eval: &dyn Fn(&[Option<usize>]) -> f64,
    ) {
        if f.len() == n1 {
            *best = best.min(eval(f));
            return;
        }
        f.push(None);
        rec(n1, f, used, best, eval);
        f.pop();
        for x in 0..used.len() {
            if !used[x] {
                used[x] = true;
                f.push(Some(x));
                rec(n1, f, used, best, eval);
                f.pop();
                used[x] = false;
            }
        }
    }
    let eval = |f: &[Option<usize>]| -> f64 {
        let mut c = f.iter().filter(|x| x.is_none()).count();
        c += b.n - f.iter().filter(|x| x.is_some()).count();
        let mut covered = std::collections::BTreeSet::new();
        for &(s, t) in &e1 {
            match (f[s], f[t]) {
                (Some(x), Some(y)) if e2.contains(&key(x, y, a.directed)) => {
                    covered.insert(key(x, y, a.directed));
                }
                _ => c += 1,
            }
        }
        (c + e2.len() - covered.len()) as f64
    };
    let mut f = Vec::with_capacity(a.n);
    let mut used = vec![false; b.n];
    let mut best = f64::INFINITY;
    rec(a.n, &mut f, &mut used, &mut best, &eval);
    best
}

#[test]
fn edit_distance_matches_brute_force() {
    let mut rng = Rng(17);
    for _ in 0..60 {
        let a = random_case(&mut rng, 5);
        let mut b = random_case(&mut rng, 5);
        b.directed = a.directed;
        let (ga, gb) = (graph(&a), graph(&b));
        let d = graph_edit_distance(&ga, &gb, EditOptions::default())
            .unwrap()
            .unwrap();
        assert_eq!(
            d,
            brute_ged(&a, &b),
            "{:?} {:?}",
            (&a.src, &a.tgt),
            (&b.src, &b.tgt)
        );
        let back = graph_edit_distance(&gb, &ga, EditOptions::default())
            .unwrap()
            .unwrap();
        assert_eq!(d, back);
        let (_, pa) = permuted(&a, &mut rng);
        let zero = graph_edit_distance(&ga, &pa, EditOptions::default())
            .unwrap()
            .unwrap();
        assert_eq!(zero, 0.0);
    }
}
