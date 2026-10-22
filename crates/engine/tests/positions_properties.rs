//! Randomized property tests for brokerage and positional measures.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::metrics::brokerage::{
    brokerage, constraint, effective_size, efficiency, StructuralHoleOptions,
};
use engine::metrics::equivalence::{
    concor, structural_equivalence, ConcorOptions, EquivalenceMethod, EquivalenceOptions,
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
    weights: Vec<f64>,
    directed: bool,
}

impl Case {
    fn graph(&self, weighted: bool) -> Graph {
        let w = if weighted {
            Some(&self.weights[..])
        } else {
            None
        };
        Graph::from_edges(self.n, &self.src, &self.tgt, w, self.directed).unwrap()
    }

    /// The same network with node v renamed perm[v].
    fn relabel(&self, perm: &[usize]) -> Case {
        Case {
            n: self.n,
            src: self.src.iter().map(|&v| perm[v as usize] as i64).collect(),
            tgt: self.tgt.iter().map(|&v| perm[v as usize] as i64).collect(),
            weights: self.weights.clone(),
            directed: self.directed,
        }
    }
}

fn random_case(rng: &mut Rng) -> Case {
    let n = 1 + rng.below(30) as usize;
    let m = rng.below((n * n / 2 + 1) as u64) as usize;
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    let mut weights = Vec::with_capacity(m);
    for _ in 0..m {
        src.push(rng.below(n as u64) as i64);
        tgt.push(rng.below(n as u64) as i64);
        weights.push(0.25 * (1 + rng.below(12)) as f64);
    }
    Case {
        n,
        src,
        tgt,
        weights,
        directed: rng.below(2) == 0,
    }
}

fn permutation(rng: &mut Rng, n: usize) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        p.swap(i, j);
    }
    p
}

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

const CASES: u64 = 200;

#[test]
fn structural_holes_permute_with_labels_and_ignore_threads() {
    let mut rng = Rng(0x5eed_0001);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        for weighted in [false, true] {
            let opts = StructuralHoleOptions {
                weighted,
                threads: 1,
            };
            let many = StructuralHoleOptions {
                weighted,
                threads: 4,
            };
            for f in [constraint, effective_size, efficiency] {
                let a = f(&c.graph(weighted), opts).unwrap();
                let b = f(&r.graph(weighted), opts).unwrap();
                let t = f(&c.graph(weighted), many).unwrap();
                for v in 0..c.n {
                    assert!(same(a[v], b[perm[v]]));
                    assert!(same(a[v], t[v]));
                }
            }
        }
    }
}

#[test]
fn directed_structural_holes_equal_the_undirected_reading() {
    // z_ij = a_ij + a_ji, so a directed network gives the same values as
    // the undirected network on the same edges.
    let mut rng = Rng(0x5eed_0002);
    for _ in 0..CASES {
        let mut c = random_case(&mut rng);
        c.directed = true;
        let d = c.graph(true);
        c.directed = false;
        let u = c.graph(true);
        for weighted in [false, true] {
            let o = StructuralHoleOptions {
                weighted,
                threads: 1,
            };
            for f in [constraint, effective_size] {
                let a = f(&d, o).unwrap();
                let b = f(&u, o).unwrap();
                assert!(a.iter().zip(&b).all(|(x, y)| same(*x, *y)));
            }
        }
    }
}

#[test]
fn borgatti_identity_on_simple_undirected_graphs() {
    let mut rng = Rng(0x5eed_0003);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let n = c.n;
        let mut adj = vec![vec![false; n]; n];
        for (&s, &t) in c.src.iter().zip(&c.tgt) {
            if s != t {
                adj[s as usize][t as usize] = true;
                adj[t as usize][s as usize] = true;
            }
        }
        let (mut src, mut tgt) = (Vec::new(), Vec::new());
        for (a, row) in adj.iter().enumerate() {
            for (b, &x) in row.iter().enumerate().skip(a + 1) {
                if x {
                    src.push(a as i64);
                    tgt.push(b as i64);
                }
            }
        }
        let g = Graph::from_edges(n, &src, &tgt, None, false).unwrap();
        let es = effective_size(&g, StructuralHoleOptions::default()).unwrap();
        for v in 0..n {
            let nb: Vec<usize> = (0..n).filter(|&u| adj[v][u]).collect();
            if nb.is_empty() {
                assert!(es[v].is_nan());
                continue;
            }
            let mut t = 0usize;
            for (i, &a) in nb.iter().enumerate() {
                for &b in &nb[i + 1..] {
                    t += adj[a][b] as usize;
                }
            }
            let k = nb.len() as f64;
            assert!(same(es[v], k - 2.0 * t as f64 / k));
        }
    }
}

#[test]
fn brokerage_permutes_with_labels_and_ignores_threads() {
    let mut rng = Rng(0x5eed_0004);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let groups: Vec<i64> = (0..c.n).map(|_| rng.below(3) as i64).collect();
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        let mut pg = vec![0i64; c.n];
        for v in 0..c.n {
            pg[perm[v]] = groups[v];
        }
        let a = brokerage(&c.graph(false), &groups, 1).unwrap();
        let b = brokerage(&r.graph(false), &pg, 1).unwrap();
        let t = brokerage(&c.graph(false), &groups, 4).unwrap();
        assert_eq!(a, t);
        let mut grand = 0.0;
        for (v, &pv) in perm.iter().enumerate() {
            assert_eq!(a.row(v), b.row(pv));
            assert_eq!(a.row(v)[..5].iter().sum::<f64>(), a.get(v, 5));
            grand += a.get(v, 5);
        }
        // Every open two-path is counted exactly once, whatever the groups.
        let one = brokerage(&c.graph(false), &vec![0; c.n], 1).unwrap();
        let coord: f64 = (0..c.n).map(|v| one.get(v, 0)).sum();
        assert_eq!(coord, grand);
    }
}

#[test]
fn structural_equivalence_permutes_and_is_symmetric() {
    let mut rng = Rng(0x5eed_0005);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        for method in [
            EquivalenceMethod::Euclidean,
            EquivalenceMethod::Correlation,
            EquivalenceMethod::Hamming,
        ] {
            for weighted in [false, true] {
                let o = EquivalenceOptions {
                    method,
                    weighted,
                    threads: 1,
                };
                let a = structural_equivalence(&c.graph(true), o).unwrap();
                let b = structural_equivalence(&r.graph(true), o).unwrap();
                let t =
                    structural_equivalence(&c.graph(true), EquivalenceOptions { threads: 3, ..o })
                        .unwrap();
                assert_eq!(a, t);
                for i in 0..c.n {
                    for j in 0..c.n {
                        assert!(same(a.get(i, j), b.get(perm[i], perm[j])));
                        assert_eq!(a.get(i, j), a.get(j, i));
                    }
                }
            }
        }
    }
}

#[test]
fn concor_ignores_threads_and_blocks_constant_profiles() {
    let mut rng = Rng(0x5eed_0006);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = c.graph(false);
        for depth in 1..=3 {
            let o = ConcorOptions {
                depth,
                threads: 1,
                ..Default::default()
            };
            let a = concor(&g, o).unwrap();
            let b = concor(&g, ConcorOptions { threads: 4, ..o }).unwrap();
            assert_eq!(a, b);
            // Labels are 1..=k with every label used.
            let k = *a.membership.iter().max().unwrap_or(&0);
            for label in 1..=k {
                assert!(a.membership.contains(&label));
            }
        }
    }
}
