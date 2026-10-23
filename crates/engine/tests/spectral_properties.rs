//! Randomized property tests for the Lanczos eigenvector solver, PageRank,
//! Katz and Bonacich power centrality, hubs and authorities, assortativity,
//! mixing matrices, and the E-I index.

extern crate knappnet as engine;

use engine::metrics::degree::DegreeMode;
use engine::metrics::eigenvector::{
    eigenvector_centrality, EigenScale, EigenvectorOptions, LoopWeight,
};
use engine::metrics::mixing::{
    assortativity_degree, assortativity_nominal, assortativity_values, ei_counts, mixing_matrix,
    AssortOptions, EiCounting, EiScope,
};
use engine::metrics::spectral::{
    bonacich_power, hits, katz, pagerank, BonacichOptions, HitsOptions, KatzOptions,
    PagerankOptions,
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
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    w: Vec<f64>,
    directed: bool,
}

fn random_case(rng: &mut Rng, max_n: u64) -> Case {
    let n = 2 + rng.below(max_n) as usize;
    let m = rng.below((n * n / 3 + 2) as u64) as usize;
    let mut c = Case {
        n,
        src: Vec::with_capacity(m),
        tgt: Vec::with_capacity(m),
        w: Vec::with_capacity(m),
        directed: rng.below(2) == 0,
    };
    for _ in 0..m {
        c.src.push(rng.below(n as u64) as i64);
        c.tgt.push(rng.below(n as u64) as i64);
        c.w.push(0.25 + 3.0 * rng.unit());
    }
    c
}

fn graph(c: &Case, weighted: bool) -> Graph {
    let w = if weighted { Some(&c.w[..]) } else { None };
    Graph::from_edges(c.n, &c.src, &c.tgt, w, c.directed).unwrap()
}

/// The same network with node v renamed perm[v].
fn relabeled(c: &Case, perm: &[usize], weighted: bool) -> Graph {
    let s: Vec<i64> = c.src.iter().map(|&v| perm[v as usize] as i64).collect();
    let t: Vec<i64> = c.tgt.iter().map(|&v| perm[v as usize] as i64).collect();
    let w = if weighted { Some(&c.w[..]) } else { None };
    Graph::from_edges(c.n, &s, &t, w, c.directed).unwrap()
}

fn permutation(rng: &mut Rng, n: usize) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        p.swap(i, rng.below(i as u64 + 1) as usize);
    }
    p
}

/// Dense matrix whose row i lists arcs into i (`into`) or out of i, with
/// undirected edges both ways and an undirected loop `loop_factor` times.
fn dense(c: &Case, weighted: bool, into: bool, loop_factor: f64) -> Vec<Vec<f64>> {
    let mut a = vec![vec![0.0; c.n]; c.n];
    for (k, (&s, &t)) in c.src.iter().zip(&c.tgt).enumerate() {
        let (s, t) = (s as usize, t as usize);
        let w = if weighted { c.w[k] } else { 1.0 };
        if s == t {
            // A directed loop counts once unless loops are ignored.
            let f = if c.directed {
                loop_factor.min(1.0)
            } else {
                loop_factor
            };
            a[s][s] += f * w;
        } else if c.directed {
            if into {
                a[t][s] += w;
            } else {
                a[s][t] += w;
            }
        } else {
            a[s][t] += w;
            a[t][s] += w;
        }
    }
    a
}

fn times(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|r| r.iter().zip(x).map(|(p, q)| p * q).sum())
        .collect()
}

fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(x, y)| {
            (x - y).abs() <= tol * (1.0 + x.abs().max(y.abs())) || (x.is_nan() && y.is_nan())
        })
}

#[test]
fn lanczos_eigenvectors_satisfy_the_eigen_equation() {
    let mut rng = Rng(11);
    for _ in 0..200 {
        let mut c = random_case(&mut rng, 40);
        c.directed = false;
        let weighted = rng.below(2) == 0;
        let g = graph(&c, weighted);
        let opts = EigenvectorOptions {
            weighted,
            scale: EigenScale::Unit,
            threads: 1,
            ..Default::default()
        };
        let e = eigenvector_centrality(&g, opts).unwrap();
        if e.eigenvalue == 0.0 {
            continue;
        }
        let a = dense(&c, weighted, true, 2.0);
        let y = times(&a, &e.values);
        for (yi, xi) in y.iter().zip(&e.values) {
            assert!(*xi >= 0.0);
            assert!((yi - e.eigenvalue * xi).abs() <= 1e-9 * (1.0 + e.eigenvalue));
        }
        // Relabeling the nodes permutes the scores.
        let perm = permutation(&mut rng, c.n);
        let h = relabeled(&c, &perm, weighted);
        let f = eigenvector_centrality(&h, opts).unwrap();
        let back: Vec<f64> = (0..c.n).map(|v| f.values[perm[v]]).collect();
        assert!(close(&e.values, &back, 1e-9));
    }
}

#[test]
fn lanczos_results_do_not_depend_on_thread_count() {
    // Large enough for parallel matrix products and several blocks, with
    // many components of different sizes.
    let mut rng = Rng(5);
    let n = 30_000usize;
    let (mut s, mut t) = (Vec::new(), Vec::new());
    for _ in 0..80_000 {
        s.push(rng.below(n as u64) as i64);
        t.push(rng.below(n as u64) as i64);
    }
    let g = Graph::from_edges(n, &s, &t, None, false).unwrap();
    let run = |threads| {
        eigenvector_centrality(
            &g,
            EigenvectorOptions {
                threads,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let one = run(1);
    assert_eq!(one, run(2));
    assert_eq!(one, run(4));
    assert!(!one.connected);
}

#[test]
fn pagerank_is_a_stationary_distribution() {
    let mut rng = Rng(21);
    for _ in 0..150 {
        let c = random_case(&mut rng, 30);
        let weighted = rng.below(2) == 0;
        let g = graph(&c, weighted);
        let opts = PagerankOptions {
            weighted,
            threads: 1,
            ..Default::default()
        };
        let x = pagerank(&g, &opts).unwrap();
        assert!((x.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        // x = d (W^T D^-1 x + dangling mass / n) + (1 - d) / n.
        let a = dense(&c, weighted, false, 2.0);
        let n = c.n as f64;
        let out: Vec<f64> = a.iter().map(|r| r.iter().sum()).collect();
        let stuck: f64 = (0..c.n).filter(|&j| out[j] == 0.0).map(|j| x[j]).sum();
        for i in 0..c.n {
            let flow: f64 = (0..c.n)
                .filter(|&j| out[j] > 0.0)
                .map(|j| a[j][i] / out[j] * x[j])
                .sum();
            let rhs = 0.85 * (flow + stuck / n) + 0.15 / n;
            assert!((x[i] - rhs).abs() < 1e-10);
        }
        let perm = permutation(&mut rng, c.n);
        let y = pagerank(&relabeled(&c, &perm, weighted), &opts).unwrap();
        let back: Vec<f64> = (0..c.n).map(|v| y[perm[v]]).collect();
        assert!(close(&x, &back, 1e-9));
    }
}

#[test]
fn katz_and_power_centrality_solve_their_linear_systems() {
    let mut rng = Rng(31);
    for _ in 0..150 {
        let c = random_case(&mut rng, 30);
        let g = graph(&c, false);
        // Spectral radius of the matrix without loops, from the engine.
        let lambda = eigenvector_centrality(
            &g,
            EigenvectorOptions {
                loops: LoopWeight::Ignore,
                ..Default::default()
            },
        )
        .map(|e| e.eigenvalue)
        .unwrap_or(0.0);
        let alpha = if lambda > 0.0 { 0.7 / lambda } else { 0.5 };
        let k = katz(
            &g,
            &KatzOptions {
                alpha,
                loops: LoopWeight::Ignore,
                normalized: false,
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let a_in = dense(&c, false, true, 0.0);
        let ax = times(&a_in, &k);
        for i in 0..c.n {
            assert!((k[i] - (1.0 + alpha * ax[i])).abs() < 1e-9 * (1.0 + k[i].abs()));
        }
        // Above 1 / lambda_max Katz refuses.
        if lambda > 0.0 {
            let bad = KatzOptions {
                alpha: 1.01 / lambda,
                loops: LoopWeight::Ignore,
                ..Default::default()
            };
            assert!(katz(&g, &bad).is_err());
        }
        for beta in [alpha, -alpha, 1.0] {
            let p = match bonacich_power(
                &g,
                BonacichOptions {
                    exponent: beta,
                    threads: 1,
                    ..Default::default()
                },
            ) {
                Ok(p) => p,
                Err(_) => continue, // singular for beta = 1
            };
            if p.iter().any(|v| v.is_nan()) {
                continue; // no ties
            }
            let sq: f64 = p.iter().map(|v| v * v).sum();
            assert!((sq - c.n as f64).abs() < 1e-8 * c.n as f64);
            // (I - beta A) p is proportional to A 1.
            let a_out = dense(&c, false, false, 0.0);
            let ap = times(&a_out, &p);
            let a1 = times(&a_out, &vec![1.0; c.n]);
            let lhs: Vec<f64> = (0..c.n).map(|i| p[i] - beta * ap[i]).collect();
            let j = (0..c.n).max_by(|&x, &y| a1[x].total_cmp(&a1[y])).unwrap();
            let scale = lhs[j] / a1[j];
            for i in 0..c.n {
                assert!((lhs[i] - scale * a1[i]).abs() < 1e-7 * (1.0 + lhs[i].abs()));
            }
        }
    }
}

#[test]
fn authorities_are_the_leading_eigenvector_of_ata() {
    let mut rng = Rng(41);
    for _ in 0..150 {
        let c = random_case(&mut rng, 30);
        let g = graph(&c, false);
        let r = hits(
            &g,
            HitsOptions {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        if r.value == 0.0 {
            assert!(r.authorities.iter().all(|v| v.is_nan()));
            continue;
        }
        let a = dense(&c, false, false, 2.0);
        // A^T A a = value a, up to the scaling of a.
        let aa = times(&a, &r.authorities);
        let ata: Vec<f64> = (0..c.n)
            .map(|j| (0..c.n).map(|i| a[i][j] * aa[i]).sum())
            .collect();
        for (x, a) in ata.iter().zip(&r.authorities) {
            assert!((x - r.value * a).abs() < 1e-8 * (1.0 + r.value));
        }
        // Hubs are proportional to A a.
        let top = aa.iter().copied().fold(0.0f64, f64::max);
        for (h, x) in r.hubs.iter().zip(&aa) {
            assert!((h - x / top).abs() < 1e-9);
        }
    }
}

#[test]
fn mixing_measures_are_consistent() {
    let mut rng = Rng(51);
    for _ in 0..300 {
        let c = random_case(&mut rng, 25);
        let g = graph(&c, false);
        let k = 1 + rng.below(4) as usize;
        let groups: Vec<i64> = (0..c.n).map(|_| rng.below(k as u64) as i64).collect();
        let m = c.src.len() as f64;
        // Mixing matrix totals: each arc once, or each edge twice.
        let mm = mixing_matrix(&g, &groups, k, &groups, k, true, DegreeMode::All).unwrap();
        let total: f64 = mm.counts.iter().sum();
        assert_eq!(total, if c.directed { m } else { 2.0 * m });
        // E-I counts: the network total, and each tie counted for both of
        // its groups (or once, inside one group).
        let global = ei_counts(
            &g,
            &groups,
            k,
            EiScope::Global,
            DegreeMode::All,
            EiCounting::Standard,
        )
        .unwrap()[0];
        assert_eq!(global.0 + global.1, m);
        let by_group = ei_counts(
            &g,
            &groups,
            k,
            EiScope::Group,
            DegreeMode::All,
            EiCounting::Standard,
        )
        .unwrap();
        let internal: f64 = by_group.iter().map(|x| x.1).sum();
        let external: f64 = by_group.iter().map(|x| x.0).sum();
        assert_eq!(internal, global.1);
        assert_eq!(external, 2.0 * global.0);
        // Assortativity does not depend on node labels or category labels,
        // and Pearson correlation ignores positive affine maps.
        let opts = AssortOptions::default();
        let r = assortativity_nominal(&g, &groups, k, opts).unwrap();
        let swapped: Vec<i64> = groups.iter().map(|&x| (k as i64 - 1) - x).collect();
        let r2 = assortativity_nominal(&g, &swapped, k, opts).unwrap();
        assert!(close(&[r], &[r2], 1e-12));
        let vals: Vec<f64> = (0..c.n).map(|_| rng.unit()).collect();
        let affine: Vec<f64> = vals.iter().map(|v| 3.0 * v + 7.0).collect();
        let rv = assortativity_values(&g, &vals, None, opts).unwrap();
        let ra = assortativity_values(&g, &affine, None, opts).unwrap();
        assert!(close(&[rv], &[ra], 1e-9));
        let perm = permutation(&mut rng, c.n);
        let h = relabeled(&c, &perm, false);
        let mut pg = vec![0i64; c.n];
        for v in 0..c.n {
            pg[perm[v]] = groups[v];
        }
        let rp = assortativity_nominal(&h, &pg, k, opts).unwrap();
        assert!(close(&[r], &[rp], 1e-12));
        let d1 = assortativity_degree(&g, false, opts).unwrap();
        let d2 = assortativity_degree(&h, false, opts).unwrap();
        assert!(close(&[d1], &[d2], 1e-12));
    }
}
