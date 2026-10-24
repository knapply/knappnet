//! Randomized property tests for the distance measures: the Wiener index
//! agrees with the engine's distance summary, efficiency stays in [0, 1],
//! and results do not depend on the thread count.

#![allow(clippy::needless_range_loop)]

extern crate knappnet as engine;

use engine::metrics::descriptive::{path_summary, PathOptions};
use engine::metrics::distance_measures::{
    closeness_vitality, global_efficiency, local_efficiency, wiener_index,
};
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

fn graph(c: &Case) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, Some(&c.w), c.directed).unwrap()
}

fn close(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || a == b || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

#[test]
fn distance_indices_agree_with_the_distance_summary() {
    let mut rng = Rng(29);
    for i in 0..60 {
        let c = random_case(&mut rng, i % 2 == 0);
        let g = graph(&c);
        let n = c.n;
        for weighted in [false, true] {
            let opts = PathOptions {
                weighted,
                threads: 1,
                ..Default::default()
            };
            let s = path_summary(&g, None, opts).unwrap();
            let w = wiener_index(&g, opts).unwrap();
            if s.pairs == (n * (n - 1)) as u64 {
                let want = if c.directed { s.total } else { s.total / 2.0 };
                assert!(close(w, want));
            } else {
                assert!(w.is_infinite());
            }
            let e = global_efficiency(&g, opts).unwrap();
            assert!((0.0..=1.0 + 1e-12).contains(&e) || weighted);
            let v1 = closeness_vitality(&g, opts).unwrap();
            let v4 = closeness_vitality(&g, PathOptions { threads: 4, ..opts }).unwrap();
            for v in 0..n {
                assert!(close(v1[v], v4[v]));
            }
        }
        if !c.directed {
            let a = local_efficiency(&g, false, 1).unwrap();
            assert_eq!(a, local_efficiency(&g, false, 4).unwrap());
            assert!((0.0..=1.0).contains(&a));
        }
    }
}
