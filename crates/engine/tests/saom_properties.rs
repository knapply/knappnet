//! Randomized checks of the stochastic actor-oriented model: invariance of
//! statistics under relabeling, reproducibility across thread counts, and
//! recovery of known parameters from simulated panels.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::rng::Stream;
use engine::saom::effects::{change_statistics_on_wave, targets, Effect, EffectKind, Model};
use engine::saom::simulate::simulate_panel;
use engine::saom::{estimate, EstimateOptions, Panel};

fn random_adjacency(rng: &mut Stream, n: usize, density: f64) -> Vec<u8> {
    let mut a = vec![0u8; n * n];
    for i in 0..n {
        for j in 0..n {
            if i != j && rng.uniform() < density {
                a[i * n + j] = 1;
            }
        }
    }
    a
}

/// A panel whose later waves toggle a share of the previous wave's ties.
fn drifting_panel(seed: u64, n: usize, waves: usize, density: f64, flip: f64) -> Panel {
    let mut rng = Stream::new(seed, 0);
    let mut adj = vec![random_adjacency(&mut rng, n, density)];
    for w in 1..waves {
        let mut a = adj[w - 1].clone();
        for i in 0..n {
            for j in 0..n {
                if i != j && rng.uniform() < flip {
                    a[i * n + j] ^= 1;
                }
            }
        }
        adj.push(a);
    }
    let mut p = Panel::from_adjacency(n, &adj).unwrap();
    let v: Vec<f64> = (0..n).map(|i| ((i * 7) % 5) as f64).collect();
    p.add_actor_covariate(&v).unwrap();
    let w: Vec<f64> = (0..n * n).map(|i| ((i * 3) % 7) as f64).collect();
    p.add_dyadic_covariate(&w).unwrap();
    p
}

fn every_effect() -> Vec<Effect> {
    EffectKind::ALL
        .iter()
        .map(|&k| {
            if k.needs_actor_covariate() || k.needs_dyadic_covariate() {
                Effect::with_covariate(k, 0)
            } else {
                Effect::structural(k)
            }
        })
        .collect()
}

fn permute_panel(p: &Panel, perm: &[usize]) -> Panel {
    let n = p.n();
    let waves: Vec<Vec<u8>> = (0..p.wave_count())
        .map(|w| {
            let a = p.wave(w).adjacency();
            let mut b = vec![0u8; n * n];
            for i in 0..n {
                for j in 0..n {
                    b[perm[i] * n + perm[j]] = a[i * n + j];
                }
            }
            b
        })
        .collect();
    let mut q = Panel::from_adjacency(n, &waves).unwrap();
    let v: Vec<f64> = (0..n).map(|i| ((i * 7) % 5) as f64).collect();
    let mut pv = vec![0.0; n];
    for i in 0..n {
        pv[perm[i]] = v[i];
    }
    q.add_actor_covariate(&pv).unwrap();
    let w: Vec<f64> = (0..n * n).map(|i| ((i * 3) % 7) as f64).collect();
    let mut pw = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            pw[perm[i] * n + perm[j]] = w[i * n + j];
        }
    }
    q.add_dyadic_covariate(&pw).unwrap();
    q
}

#[test]
fn relabeling_actors_leaves_targets_unchanged_and_permutes_change_statistics() {
    for seed in 0..5 {
        let n = 9 + seed as usize;
        let p = drifting_panel(seed, n, 3, 0.25, 0.1);
        let mut perm: Vec<usize> = (0..n).collect();
        Stream::new(seed, 99).permutation(&mut perm);
        let q = permute_panel(&p, &perm);
        let m_p = Model::new(&p, &every_effect()).unwrap();
        let m_q = Model::new(&q, &every_effect()).unwrap();
        let tp = targets(&m_p, &p);
        let tq = targets(&m_q, &q);
        for (a, b) in tp.iter().zip(&tq) {
            assert!((a - b).abs() < 1e-9, "seed {seed}: {a} vs {b}");
        }
        let k = m_p.len();
        for i in 0..n {
            let dp = change_statistics_on_wave(&m_p, &p, 1, 1, i);
            let dq = change_statistics_on_wave(&m_q, &q, 1, 1, perm[i]);
            for j in 0..n {
                for kk in 0..k {
                    let a = dp[j * k + kk];
                    let b = dq[perm[j] * k + kk];
                    assert!((a - b).abs() < 1e-9, "seed {seed} i {i} j {j} effect {kk}");
                }
            }
        }
    }
}

#[test]
fn simulations_are_reproducible_and_streams_differ() {
    let p = drifting_panel(3, 14, 2, 0.2, 0.1);
    let m = Model::new(
        &p,
        &[
            Effect::structural(EffectKind::Density),
            Effect::structural(EffectKind::Recip),
            Effect::structural(EffectKind::TransTrip),
        ],
    )
    .unwrap();
    let beta = [-1.8, 1.0, 0.2];
    let a = simulate_panel(&p, &m, &beta, &[], true, 5, 0).unwrap();
    let b = simulate_panel(&p, &m, &beta, &[], true, 5, 0).unwrap();
    assert_eq!(a.networks, b.networks);
    assert_eq!(a.statistics, b.statistics);
    let c = simulate_panel(&p, &m, &beta, &[], true, 5, 1).unwrap();
    assert_ne!(a.networks, c.networks);
    // The statistics reported are those of the returned networks.
    let s = m.statistics(&p, 0, &a.networks[0]);
    assert_eq!(s, a.statistics[..3].to_vec());
}

/// Simulates waves 2 and 3 from wave 1 under known parameters.
fn simulated_panel(seed: u64, n: usize, rate: f64, beta: &[f64], effects: &[Effect]) -> Panel {
    let mut rng = Stream::new(seed, 0);
    let w1 = random_adjacency(&mut rng, n, 0.08);
    let mut waves = vec![w1];
    for step in 0..2 {
        let last = waves.last().unwrap().clone();
        let scaffold = Panel::from_adjacency(n, &[last.clone(), last]).unwrap();
        let m = Model::new(&scaffold, effects).unwrap();
        let sim = simulate_panel(&scaffold, &m, beta, &[rate], false, seed, 100 + step).unwrap();
        waves.push(sim.networks[0].adjacency().to_vec());
    }
    Panel::from_adjacency(n, &waves).unwrap()
}

#[test]
fn estimates_recover_the_generating_parameters() {
    let effects = [
        Effect::structural(EffectKind::Density),
        Effect::structural(EffectKind::Recip),
        Effect::structural(EffectKind::TransTrip),
    ];
    let truth = [-2.2, 1.5, 0.4];
    let rate = 5.0;
    let mut covered = 0;
    let mut total = 0;
    let mut z_sum = [0.0f64; 3];
    let replicates = 6;
    for r in 0..replicates {
        let p = simulated_panel(40 + r, 30, rate, &truth, &effects);
        let m = Model::new(&p, &effects).unwrap();
        let opts = EstimateOptions {
            n3: 400,
            nsub: 3,
            seed: 1000 + r,
            threads: 2,
            keep_aux: false,
            ..Default::default()
        };
        let fit = estimate(&p, &m, &opts).unwrap();
        assert!(
            fit.tconv_max < 0.5,
            "replicate {r}: tconv.max {}",
            fit.tconv_max
        );
        for k in 0..3 {
            let z = (fit.theta[k] - truth[k]) / fit.se[k];
            z_sum[k] += z;
            total += 1;
            if z.abs() < 2.5 {
                covered += 1;
            }
        }
        // The rate is estimated from stopping times; the period length is
        // 1, so it should be near the generating rate.
        for &rho in &fit.rates {
            assert!((rho - rate).abs() < 0.4 * rate, "replicate {r}: rate {rho}");
        }
    }
    // Most 2.5-sigma intervals cover the truth, and the standardized
    // errors do not all lean the same way.
    assert!(covered >= total - 2, "{covered} of {total} covered");
    for (k, z) in z_sum.iter().enumerate() {
        let mean_z = z / replicates as f64;
        assert!(
            mean_z.abs() < 1.2,
            "effect {k}: mean standardized error {mean_z}"
        );
    }
}

#[test]
fn unconditional_and_conditional_estimates_agree() {
    let effects = [
        Effect::structural(EffectKind::Density),
        Effect::structural(EffectKind::Recip),
    ];
    let p = simulated_panel(7, 25, 4.0, &[-2.0, 1.2], &effects);
    let m = Model::new(&p, &effects).unwrap();
    let cond = estimate(
        &p,
        &m,
        &EstimateOptions {
            n3: 300,
            nsub: 3,
            seed: 2,
            threads: 2,
            keep_aux: false,
            ..Default::default()
        },
    )
    .unwrap();
    let uncond = estimate(
        &p,
        &m,
        &EstimateOptions {
            conditional: false,
            n3: 300,
            nsub: 3,
            seed: 3,
            threads: 2,
            keep_aux: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(uncond.theta.len(), 4);
    for k in 0..2 {
        let a = cond.theta[k];
        let b = uncond.theta[2 + k];
        let tol = 3.0 * (cond.se[k].powi(2) + uncond.se[2 + k].powi(2)).sqrt();
        assert!((a - b).abs() < tol, "effect {k}: {a} vs {b} (tol {tol})");
    }
    for m in 0..2 {
        let a = cond.rates[m];
        let b = uncond.rates[m];
        assert!(
            (a - b).abs() < 0.5 * a.max(b),
            "period {m}: rates {a} vs {b}"
        );
    }
}
