//! Randomized property tests for greedy colourings, maximal independent
//! sets, square clustering, generalized degree, and cut measures: validity,
//! seeds, thread counts, and relabeling.

extern crate knappnet as engine;

use engine::metrics::coloring::{greedy_color, maximal_independent_set, ColoringStrategy};
use engine::metrics::cut_measures::cut_measures;
use engine::metrics::squares::{generalized_degree, square_clustering};
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

fn random(rng: &mut Rng, directed: bool) -> (usize, Vec<i64>, Vec<i64>, Graph) {
    let n = 1 + rng.below(30) as usize;
    let m = rng.below((n * n / 3 + 1) as u64) as usize;
    let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
    (n, s, t, g)
}

#[test]
fn colourings_are_proper_and_sets_independent() {
    let mut rng = Rng(19);
    for _ in 0..80 {
        let (n, s, t, g) = random(&mut rng, false);
        for code in 0..7 {
            let c = greedy_color(&g, ColoringStrategy::from_code(code).unwrap(), 5).unwrap();
            assert_eq!(c.len(), n);
            for (&a, &b) in s.iter().zip(&t) {
                if a != b {
                    assert_ne!(c[a as usize], c[b as usize], "strategy {code}");
                }
            }
            // Greedy colourings never use more than max degree + 1 colours.
            let maxdeg = (0..n).map(|v| g.out().degree(v)).max().unwrap_or(0);
            assert!(c.iter().all(|&x| x as usize <= maxdeg));
        }
        let a = maximal_independent_set(&g, &[], 3).unwrap();
        assert_eq!(a, maximal_independent_set(&g, &[], 3).unwrap());
        let inset: Vec<bool> = (0..n).map(|v| a.contains(&(v as u32))).collect();
        for (&x, &y) in s.iter().zip(&t) {
            if x != y {
                assert!(!(inset[x as usize] && inset[y as usize]));
            }
        }
        for v in 0..n {
            let covered = inset[v]
                || s.iter().zip(&t).any(|(&x, &y)| {
                    (x as usize == v && inset[y as usize]) || (y as usize == v && inset[x as usize])
                });
            assert!(covered, "node {v} could be added");
        }
    }
}

#[test]
fn squares_and_cuts_are_thread_and_label_independent() {
    let mut rng = Rng(23);
    for round in 0..60 {
        let (n, s, t, g) = random(&mut rng, round % 2 == 0);
        let one = square_clustering(&g, 1).unwrap();
        assert_eq!(one, square_clustering(&g, 4).unwrap());
        assert!(one.iter().all(|&x| (0.0..=1.0).contains(&x)));
        if !g.is_directed() {
            assert_eq!(generalized_degree(&g, 1), generalized_degree(&g, 3));
        }
        // Relabeling permutes square clustering and leaves cuts unchanged.
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            perm.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let ps: Vec<i64> = s.iter().map(|&x| perm[x as usize] as i64).collect();
        let pt: Vec<i64> = t.iter().map(|&x| perm[x as usize] as i64).collect();
        let h = Graph::from_edges(n, &ps, &pt, None, g.is_directed()).unwrap();
        let sq = square_clustering(&h, 1).unwrap();
        for v in 0..n {
            assert!((one[v] - sq[perm[v]]).abs() < 1e-12);
        }
        let set: Vec<u32> = (0..n as u32).filter(|v| v % 2 == 0).collect();
        let pset: Vec<u32> = set.iter().map(|&v| perm[v as usize] as u32).collect();
        let a = cut_measures(&g, &set, None, false).unwrap();
        let b = cut_measures(&h, &pset, None, false).unwrap();
        assert_eq!(a.cut_size, b.cut_size);
        assert_eq!(a.volume_s, b.volume_s);
        assert_eq!(a.node_expansion, b.node_expansion);
    }
}
