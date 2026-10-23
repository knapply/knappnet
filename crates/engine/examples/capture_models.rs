//! Prints the results of the model, layout, permutation test, generator,
//! and linear algebra workloads on fixed random inputs, so that a change
//! to their code can be checked for identical output: run before and after
//! and compare the files.
//!
//!     cargo run --release --offline --example capture_models > before.txt

use knappnet::generate::{gnm, gnp, sbm};
use knappnet::layout::{fr_layout, stress_layout, FrOptions, StressAlgorithm, StressOptions};
use knappnet::linalg::{
    cholesky, cholesky_inverse, invert, matmul, pseudo_inverse, solve, solve_spd, symmetric_eigen,
};
use knappnet::matrix::Dense;
use knappnet::metrics::cliques::max_cliques;
use knappnet::metrics::community::louvain;
use knappnet::metrics::coreperiphery::discrete::{core_periphery, CorePeripheryOptions};
use knappnet::metrics::equivalence::{concor, ConcorOptions};
use knappnet::metrics::factions::kgroups::{factions, FactionOptions};
use knappnet::metrics::partition::fast_greedy::fast_greedy;
use knappnet::metrics::partition::girvan_newman::girvan_newman;
use knappnet::metrics::partition::leading_eigen::leading_eigen;
use knappnet::metrics::regular::regular_equivalence;
use knappnet::metrics::scaling::classical_mds;
use knappnet::qap::{
    mrqap, qap_cor, MrqapMethod, MrqapOptions, Permutations, QapCorOptions, QapMode,
};
use knappnet::rng::Stream;
use knappnet::Graph;

fn graph(n: usize, p: f64, directed: bool, weighted: bool, seed: u64) -> Graph {
    let e = gnp(n, p, directed, false, 1, seed, 1).unwrap().remove(0);
    let src: Vec<i64> = e.from.iter().map(|&v| v as i64).collect();
    let tgt: Vec<i64> = e.to.iter().map(|&v| v as i64).collect();
    let w: Option<Vec<f64>> = weighted.then(|| {
        let mut r = Stream::new(seed, 99);
        (0..src.len()).map(|_| 1.0 + r.below(5) as f64).collect()
    });
    Graph::from_edges(n, &src, &tgt, w.as_deref(), directed).unwrap()
}

fn matrix(n: usize, seed: u64, symmetric: bool, nan: bool) -> Dense {
    let mut r = Stream::new(seed, 7);
    let mut m = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            if symmetric && j > i {
                continue;
            }
            let mut x = r.uniform();
            if nan && r.below(50) == 0 {
                x = f64::NAN;
            }
            m.set(i, j, x);
            if symmetric {
                m.set(j, i, x);
            }
        }
    }
    m
}

fn bits(v: &[f64]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for x in v {
        for b in x.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    format!("{h:016x}")
}

fn main() {
    for (n, p, weighted) in [
        (300, 0.03, false),
        (300, 0.03, true),
        (200, 0.3, false),
        (2000, 0.004, false),
    ] {
        let g = graph(n, p, false, weighted, 1);
        let l = louvain(&g, 1.0, weighted, 7).unwrap();
        println!(
            "louvain {n} {p} {weighted}: {:?} {} {}",
            l.membership, l.modularity, l.levels
        );
        let f = fast_greedy(&g, 1.0, weighted).unwrap();
        println!(
            "fast_greedy {n} {p} {weighted}: {:?} {} {:?}",
            f.membership, f.modularity, f.merges
        );
        let f = fast_greedy(&g, 0.7, weighted).unwrap();
        println!(
            "fast_greedy 0.7 {n} {p} {weighted}: {:?} {} {:?}",
            f.membership, f.modularity, f.merges
        );
        if n <= 300 {
            let e = leading_eigen(&g, weighted).unwrap();
            println!(
                "leading_eigen {n} {p} {weighted}: {:?} {} {}",
                e.membership, e.modularity, e.splits
            );
        }
        if n <= 200 {
            for t in [1, 2] {
                let r = girvan_newman(&g, t).unwrap();
                println!(
                    "girvan_newman {n} {p} t={t}: {:?} {} {:?} {}",
                    r.membership,
                    r.modularity,
                    r.removed,
                    bits(&r.path)
                );
            }
        }
    }
    // A graph with repeated merges and ties (ring of cliques) for fast greedy.
    {
        let mut e = Vec::new();
        for c in 0..30i64 {
            for i in 0..5 {
                for j in (i + 1)..5 {
                    e.push((5 * c + i, 5 * c + j));
                }
            }
            e.push((5 * c, (5 * (c + 1)) % 150));
        }
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        let g = Graph::from_edges(150, &s, &t, None, false).unwrap();
        let f = fast_greedy(&g, 1.0, false).unwrap();
        println!(
            "fast_greedy ring: {:?} {} {:?}",
            f.membership, f.modularity, f.merges
        );
        let r = girvan_newman(&g, 2).unwrap();
        println!(
            "girvan_newman ring: {:?} {} {:?} {}",
            r.membership,
            r.modularity,
            r.removed,
            bits(&r.path)
        );
        let l = leading_eigen(&g, false).unwrap();
        println!("leading_eigen ring: {:?} {}", l.membership, l.modularity);
    }
    for (n, p) in [(500, 0.02), (100, 0.3)] {
        let g = graph(n, p, false, false, 5);
        let mut c = max_cliques(&g, 3, 2);
        c.sort();
        println!(
            "cliques {n} {p}: {} {:?}",
            c.len(),
            c.iter().take(20).collect::<Vec<_>>()
        );
        for directed in [false, true] {
            let g = graph(n, p, directed, false, 6);
            for t in [1, 2] {
                let f = factions(
                    &g,
                    FactionOptions {
                        k: 3,
                        threads: t,
                        starts: 6,
                        ..Default::default()
                    },
                )
                .unwrap();
                println!(
                    "factions {n} {p} {directed} t={t}: {:?} {} {}",
                    f.membership, f.errors, f.correlation
                );
                let cp = core_periphery(
                    &g,
                    CorePeripheryOptions {
                        threads: t,
                        starts: 6,
                        ..Default::default()
                    },
                )
                .unwrap();
                println!(
                    "coreperiphery {n} {p} {directed} t={t}: {:?} {}",
                    cp.core, cp.correlation
                );
                let c = concor(
                    &g,
                    ConcorOptions {
                        depth: 3,
                        threads: t,
                        ..Default::default()
                    },
                )
                .unwrap();
                println!(
                    "concor {n} {p} {directed} t={t}: {:?} {} {}",
                    c.membership, c.unconverged, c.near_zero
                );
            }
            let r = regular_equivalence(&g, None, Default::default()).unwrap();
            println!(
                "regular {n} {p} {directed}: {:?} {}",
                r.membership, r.rounds
            );
        }
    }
    for (n, p, weighted) in [
        (300, 0.02, false),
        (300, 0.02, true),
        (150, 0.3, false),
        (3000, 0.003, false),
    ] {
        let g = graph(n, p, false, weighted, 9);
        for alg in [StressAlgorithm::Full, StressAlgorithm::Sparse] {
            if alg == StressAlgorithm::Full && n > 2000 {
                continue;
            }
            for t in [1, 2] {
                let opts = StressOptions {
                    algorithm: alg,
                    weighted,
                    threads: t,
                    ..Default::default()
                };
                let l = stress_layout(&g, &opts).unwrap();
                println!(
                    "stress {n} {p} {weighted} {alg:?} t={t}: {} {}",
                    bits(&l.x),
                    bits(&l.y)
                );
            }
        }
        for t in [1, 2] {
            let opts = FrOptions {
                weighted,
                threads: t,
                ..Default::default()
            };
            let l = fr_layout(&g, &opts).unwrap();
            println!("fr {n} {p} {weighted} t={t}: {} {}", bits(&l.x), bits(&l.y));
        }
    }
    for nan in [false, true] {
        let n = 60;
        let y = matrix(n, 11, false, nan);
        let x1 = matrix(n, 12, false, nan);
        let x2 = matrix(n, 13, false, false);
        let x3 = matrix(n, 14, true, false);
        for mode in [QapMode::Digraph, QapMode::Graph] {
            for diag in [false, true] {
                for t in [1, 2] {
                    let opts = QapCorOptions {
                        mode,
                        diag,
                        permutations: Permutations::Random { reps: 300, seed: 1 },
                        threads: t,
                        keep_draws: true,
                    };
                    let r = qap_cor(&y, &x1, opts).unwrap();
                    println!(
                        "qap_cor nan={nan} {mode:?} {diag} t={t}: {} {:?} {}",
                        r.statistic,
                        r.tail,
                        bits(&r.distribution)
                    );
                    for method in [
                        MrqapMethod::DoubleSemiPartialling,
                        MrqapMethod::YPermutation,
                    ] {
                        for intercept in [true, false] {
                            let opts = MrqapOptions {
                                mode,
                                diag,
                                method,
                                intercept,
                                permutations: Permutations::Random { reps: 300, seed: 1 },
                                threads: t,
                                keep_draws: true,
                                ..Default::default()
                            };
                            let r = mrqap(&y, &[x1.clone(), x2.clone(), x3.clone()], opts).unwrap();
                            println!("mrqap nan={nan} {mode:?} {diag} {method:?} {intercept} t={t}: {:?} {:?} {:?} {:?} {} {} {}", r.coefficients, r.std_errors, r.tails, r.method, r.rss, r.f_statistic, bits(&r.distribution));
                        }
                    }
                }
            }
        }
    }
    {
        let n = 7;
        let y = matrix(n, 31, false, false);
        let x1 = matrix(n, 32, false, false);
        let opts = QapCorOptions {
            permutations: Permutations::Exhaustive,
            threads: 2,
            keep_draws: true,
            ..Default::default()
        };
        let r = qap_cor(&y, &x1, opts).unwrap();
        println!(
            "qap_cor exhaustive: {} {:?} {}",
            r.statistic,
            r.tail,
            bits(&r.distribution)
        );
        let opts = MrqapOptions {
            permutations: Permutations::Exhaustive,
            threads: 2,
            keep_draws: true,
            ..Default::default()
        };
        let r = mrqap(&y, std::slice::from_ref(&x1), opts).unwrap();
        println!(
            "mrqap exhaustive: {:?} {:?} {}",
            r.coefficients,
            r.tails,
            bits(&r.distribution)
        );
    }
    for (n, p) in [
        (50, 0.1),
        (50, 0.6),
        (500, 0.01),
        (500, 0.4),
        (3, 1.0),
        (1, 0.5),
        (0, 0.5),
    ] {
        for directed in [false, true] {
            for loops in [false, true] {
                let e = gnp(n, p, directed, loops, 5, 3, 2).unwrap();
                let s: Vec<String> = e
                    .iter()
                    .map(|e| {
                        format!(
                            "{:?}/{:?}",
                            bits(&e.from.iter().map(|&v| v as f64).collect::<Vec<_>>()),
                            bits(&e.to.iter().map(|&v| v as f64).collect::<Vec<_>>())
                        )
                    })
                    .collect();
                println!("gnp {n} {p} {directed} {loops}: {} {:?}", e[0].len(), s);
                let total = e[0].len() as u64;
                let e = gnm(n, total, directed, loops, 3, 4, 2).unwrap();
                let s: Vec<String> = e
                    .iter()
                    .map(|e| {
                        format!(
                            "{:?}/{:?}",
                            bits(&e.from.iter().map(|&v| v as f64).collect::<Vec<_>>()),
                            bits(&e.to.iter().map(|&v| v as f64).collect::<Vec<_>>())
                        )
                    })
                    .collect();
                println!("gnm {n} {p} {directed} {loops}: {:?}", s);
            }
        }
    }
    {
        let e = sbm(
            &[20, 30, 10],
            &[0.5, 0.05, 0.1, 0.05, 0.4, 0.02, 0.1, 0.02, 0.6],
            false,
            false,
            3,
            8,
            2,
        )
        .unwrap();
        let s: Vec<String> = e
            .iter()
            .map(|e| {
                format!(
                    "{}:{}/{}",
                    e.len(),
                    bits(&e.from.iter().map(|&v| v as f64).collect::<Vec<_>>()),
                    bits(&e.to.iter().map(|&v| v as f64).collect::<Vec<_>>())
                )
            })
            .collect();
        println!("sbm: {s:?}");
    }
    for n in [1, 2, 3, 10, 40, 150, 600] {
        let a = matrix(n, 21, true, false);
        let (vals, vecs) = symmetric_eigen(&a).unwrap();
        println!("eigen {n}: {} {}", bits(&vals), bits(vecs.as_slice()));
        let b = matrix(n, 23, false, false);
        let mut s = matmul(&knappnet::linalg::transpose(&b), &b);
        for i in 0..n {
            s.set(i, i, s.get(i, i) + 1.0);
        }
        let l = cholesky(&s).unwrap();
        println!("cholesky {n}: {}", bits(l.as_slice()));
        println!(
            "cholesky_inverse {n}: {}",
            bits(cholesky_inverse(&s).unwrap().as_slice())
        );
        let rhs: Vec<f64> = (0..n).map(|i| i as f64 * 0.5 - 1.0).collect();
        println!("solve_spd {n}: {}", bits(&solve_spd(&s, &rhs).unwrap()));
        println!("solve {n}: {}", bits(&solve(&b, &rhs, 1e-12).unwrap()));
        println!(
            "invert {n}: {}",
            bits(invert(&b, 1e-12).unwrap().as_slice())
        );
        println!("matmul {n}: {}", bits(matmul(&a, &b).as_slice()));
        let (pi, rank) = pseudo_inverse(&a, 1e-10).unwrap();
        println!("pinv {n}: {rank} {}", bits(pi.as_slice()));
        if n >= 3 {
            let mut d = Dense::zeros(n, n);
            for i in 0..n {
                for j in 0..n {
                    d.set(
                        i,
                        j,
                        if i == j {
                            0.0
                        } else {
                            ((i as f64 - j as f64).abs() + 1.0).ln()
                        },
                    );
                }
            }
            let m = classical_mds(&d, 2).unwrap();
            println!("mds {n}: {} {:?}", bits(m.points.as_slice()), m.gof);
        }
    }
}
