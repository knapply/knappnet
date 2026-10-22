//! Randomized property tests for the quadratic assignment procedure.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::matrix::Dense;
use engine::qap::{
    mrqap, qap_cor, MrqapMethod, MrqapOptions, Permutations, QapCorOptions, QapMode,
};
use engine::rng::Stream;

fn random_matrix(s: &mut Stream, n: usize, binary: bool) -> Dense {
    let mut m = Dense::zeros(n, n);
    for r in 0..n {
        for c in 0..n {
            let u = s.uniform();
            m.set(
                r,
                c,
                if binary {
                    (u < 0.4) as u8 as f64
                } else {
                    u * 4.0 - 2.0
                },
            );
        }
    }
    m
}

fn relabel(m: &Dense, perm: &[usize]) -> Dense {
    let n = m.rows();
    let mut out = Dense::zeros(n, n);
    for r in 0..n {
        for c in 0..n {
            out.set(r, c, m.get(perm[r], perm[c]));
        }
    }
    out
}

fn affine(m: &Dense, a: f64, b: f64) -> Dense {
    let n = m.rows();
    let mut out = Dense::zeros(n, n);
    for r in 0..n {
        for c in 0..n {
            out.set(r, c, a * m.get(r, c) + b);
        }
    }
    out
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// Exact p-values do not depend on how the nodes are labeled: relabeling
/// every matrix the same way leaves the observed statistics unchanged and
/// permutes the reference distribution onto itself.
#[test]
fn exact_results_are_invariant_to_relabeling() {
    for case in 0..12u64 {
        let mut s = Stream::new(31, case);
        let n = 4 + (case % 3) as usize;
        let mode = if case % 2 == 0 {
            QapMode::Digraph
        } else {
            QapMode::Graph
        };
        let y = random_matrix(&mut s, n, case % 4 == 1);
        let x1 = random_matrix(&mut s, n, false);
        let x2 = random_matrix(&mut s, n, true);
        let mut perm = vec![0; n];
        s.permutation(&mut perm);
        // Graph mode reads the lower triangle, which relabeling does not
        // preserve unless the matrices are symmetric.
        let sym = |m: &Dense| {
            let mut o = m.clone();
            for r in 0..n {
                for c in 0..r {
                    o.set(c, r, m.get(r, c));
                }
            }
            o
        };
        let prep = |m: &Dense| {
            if mode == QapMode::Graph {
                sym(m)
            } else {
                m.clone()
            }
        };
        let (y, x1, x2) = (prep(&y), prep(&x1), prep(&x2));
        let (yp, x1p, x2p) = (relabel(&y, &perm), relabel(&x1, &perm), relabel(&x2, &perm));

        let copts = QapCorOptions {
            mode,
            permutations: Permutations::Exhaustive,
            ..Default::default()
        };
        let a = qap_cor(&y, &x1, copts).unwrap();
        let b = qap_cor(&yp, &x1p, copts).unwrap();
        assert!(close(a.statistic, b.statistic));
        assert_eq!(a.tail, b.tail, "case {case}");

        for method in [
            MrqapMethod::DoubleSemiPartialling,
            MrqapMethod::YPermutation,
        ] {
            let opts = MrqapOptions {
                mode,
                method,
                permutations: Permutations::Exhaustive,
                ..Default::default()
            };
            let a = mrqap(&y, &[x1.clone(), x2.clone()], opts).unwrap();
            let b = mrqap(&yp, &[x1p.clone(), x2p.clone()], opts).unwrap();
            for k in 0..3 {
                assert!(close(a.coefficients[k], b.coefficients[k]));
                assert_eq!(a.tails[k], b.tails[k], "case {case} {method:?} column {k}");
            }
        }
    }
}

/// t-statistics, and therefore p-values, are unchanged by a positive
/// affine transformation of y; a negative scale swaps the tails.
#[test]
fn t_tests_are_invariant_to_rescaling_y() {
    for case in 0..8u64 {
        let mut s = Stream::new(77, case);
        let n = 5;
        let y = random_matrix(&mut s, n, false);
        let xs = [
            random_matrix(&mut s, n, false),
            random_matrix(&mut s, n, true),
        ];
        for method in [
            MrqapMethod::DoubleSemiPartialling,
            MrqapMethod::YPermutation,
        ] {
            let opts = MrqapOptions {
                method,
                permutations: Permutations::Exhaustive,
                ..Default::default()
            };
            let a = mrqap(&y, &xs, opts).unwrap();
            let b = mrqap(&affine(&y, 3.5, 0.0), &xs, opts).unwrap();
            let c = mrqap(&affine(&y, -2.0, 0.0), &xs, opts).unwrap();
            for k in 0..3 {
                assert!(close(a.t[k], b.t[k]));
                assert_eq!(a.tails[k], b.tails[k]);
                assert!(close(a.t[k], -c.t[k]));
                assert_eq!(a.tails[k].le, c.tails[k].ge);
                assert_eq!(a.tails[k].abs_ge, c.tails[k].abs_ge);
            }
            // A shift only moves the intercept.
            let d = mrqap(&affine(&y, 1.0, 10.0), &xs, opts).unwrap();
            for k in 1..3 {
                assert!(close(a.t[k], d.t[k]));
                assert_eq!(a.tails[k], d.tails[k]);
            }
        }
    }
}

/// With one predictor and an intercept, the Y-permutation t-test and the
/// QAP correlation test are the same test: t is a monotone function of r.
#[test]
fn simple_regression_agrees_with_correlation() {
    for case in 0..10u64 {
        let mut s = Stream::new(5, case);
        let n = 5 + (case % 2) as usize;
        let y = random_matrix(&mut s, n, case % 3 == 0);
        let x = random_matrix(&mut s, n, false);
        let c = qap_cor(
            &y,
            &x,
            QapCorOptions {
                permutations: Permutations::Exhaustive,
                ..Default::default()
            },
        )
        .unwrap();
        let m = mrqap(
            &y,
            std::slice::from_ref(&x),
            MrqapOptions {
                method: MrqapMethod::YPermutation,
                permutations: Permutations::Exhaustive,
                ..Default::default()
            },
        )
        .unwrap();
        // t = r sqrt(N - 2) / sqrt(1 - r^2).
        let r = c.statistic;
        let df = (c.dyads - 2) as f64;
        assert!(close(m.t[1], r * df.sqrt() / (1.0 - r * r).sqrt()));
        assert!(close(m.r_squared, r * r));
        assert_eq!(m.tails[1], c.tail, "case {case}");
    }
}

/// Monte Carlo p-values converge to the exact ones.
#[test]
fn monte_carlo_is_close_to_exact() {
    let mut s = Stream::new(2, 0);
    let y = random_matrix(&mut s, 6, false);
    let x = affine(&y, 0.3, 0.0);
    let noise = random_matrix(&mut s, 6, false);
    let mut xn = Dense::zeros(6, 6);
    for r in 0..6 {
        for c in 0..6 {
            xn.set(r, c, x.get(r, c) + noise.get(r, c));
        }
    }
    let exact = qap_cor(
        &y,
        &xn,
        QapCorOptions {
            permutations: Permutations::Exhaustive,
            ..Default::default()
        },
    )
    .unwrap();
    let reps = 50_000;
    let mc = qap_cor(
        &y,
        &xn,
        QapCorOptions {
            permutations: Permutations::Random { reps, seed: 9 },
            ..Default::default()
        },
    )
    .unwrap();
    for (e, m) in [
        (exact.tail.p_le(), mc.tail.p_le()),
        (exact.tail.p_ge(), mc.tail.p_ge()),
        (exact.tail.p_abs(), mc.tail.p_abs()),
    ] {
        let se = (e * (1.0 - e) / reps as f64).sqrt();
        assert!((e - m).abs() <= 5.0 * se + 1e-4, "{e} vs {m}");
    }
}
