//! Randomized property tests for NetworkX-style mixing: agreement with the
//! igraph-style coefficients where the conventions coincide, relabeling,
//! subsets, and the mixing counts.

extern crate knappnet as engine;

use engine::metrics::assortativity::{
    attribute_assortativity, average_degree_connectivity, average_neighbor_degree,
    degree_assortativity, degree_pairs, mixing_counts, numeric_assortativity, DegreeKind,
};
use engine::metrics::eigenvector::LoopWeight;
use engine::metrics::mixing::{
    assortativity_degree, assortativity_nominal, assortativity_values, AssortOptions,
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

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-10 * a.abs().max(1.0)
}

#[test]
fn coefficients_agree_with_the_igraph_style_ones() {
    let mut rng = Rng(31);
    for _ in 0..80 {
        let n = 2 + rng.below(25) as usize;
        let m = rng.below((n * n / 2 + 1) as u64) as usize;
        let mut s = Vec::new();
        let mut t = Vec::new();
        for _ in 0..m {
            let (a, b) = (rng.below(n as u64) as i64, rng.below(n as u64) as i64);
            // NetworkX counts an undirected loop once, igraph twice: keep
            // loops out where the two are compared.
            if a != b {
                s.push(a);
                t.push(b);
            }
        }
        let directed = rng.below(2) == 0;
        let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
        let o = AssortOptions {
            directed: true,
            loops: LoopWeight::Once,
            normalized: true,
        };
        let a = degree_assortativity(&g, DegreeKind::Out, DegreeKind::In, false, None).unwrap();
        let b = assortativity_degree(&g, false, o).unwrap();
        assert!(same(a, b), "{a} {b}");
        let types: Vec<i64> = (0..n).map(|_| rng.below(3) as i64).collect();
        let labels: Vec<f64> = types.iter().map(|&x| x as f64).collect();
        let a = attribute_assortativity(&g, &labels, None).unwrap();
        let b = assortativity_nominal(&g, &types, 3, o).unwrap();
        if s.is_empty() {
            assert_eq!(a, 0.0);
        } else {
            assert!(same(a, b), "{a} {b}");
        }
        let a = numeric_assortativity(&g, &labels, None).unwrap();
        let b = assortativity_values(&g, &labels, None, o).unwrap();
        assert!(same(a, b));
        // Mixing counts add up to the number of pairs.
        let p = degree_pairs(&g, DegreeKind::Out, DegreeKind::In, false, None).unwrap();
        let c = mixing_counts(&p).unwrap();
        assert_eq!(c.counts.iter().sum::<f64>() as usize, p.len());
        // A subset of every node is the whole.
        let all: Vec<u32> = (0..n as u32).collect();
        let t3 = DegreeKind::Total;
        assert_eq!(
            average_neighbor_degree(&g, t3, t3, false, Some(&all)).unwrap(),
            average_neighbor_degree(&g, t3, t3, false, None).unwrap()
        );
        let (k, v) = average_degree_connectivity(&g, t3, t3, false, None).unwrap();
        assert!(k.windows(2).all(|w| w[0] < w[1]));
        assert!(v.iter().all(|x| x.is_finite() && *x >= 0.0));
    }
}
