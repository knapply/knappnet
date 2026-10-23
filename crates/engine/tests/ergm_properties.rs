//! Randomized property tests for exponential random graph models.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::ergm::terms::{Covariates, DegreeKind, Endpoint, Partners, Term, Terms};
use engine::ergm::{Constraints, Model, NetState, Proposal, SampleOptions};
use engine::matrix::Dense;
use engine::rng::Stream;

fn random_state(s: &mut Stream, n: usize, directed: bool, density: f64) -> NetState {
    let mut y = NetState::empty(n, directed);
    for i in 0..n {
        for j in 0..n {
            if i == j || (!directed && j < i) {
                continue;
            }
            if s.uniform() < density {
                y.toggle(i, j);
            }
        }
    }
    y
}

fn random_covariates(s: &mut Stream, n: usize) -> Covariates {
    let numeric = vec![
        (0..n).map(|_| (s.below(10) as f64) / 2.0).collect(),
        (0..n).map(|_| s.uniform() * 3.0 - 1.0).collect(),
    ];
    let categorical = vec![(0..n).map(|_| s.below(3) as u32).collect()];
    let mut w = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            w.set(i, j, s.uniform() * 2.0 - 1.0);
        }
    }
    Covariates {
        numeric,
        categorical,
        dyadic: vec![w],
    }
}

fn all_terms(directed: bool) -> Vec<Term> {
    let mut t = vec![
        Term::Edges,
        Term::NodeCov {
            attr: 0,
            endpoint: Endpoint::Both,
        },
        Term::NodeFactor {
            attr: 0,
            endpoint: Endpoint::Both,
            levels: vec![1, 2],
        },
        Term::NodeMatch {
            attr: 0,
            levels: vec![],
        },
        Term::NodeMatch {
            attr: 0,
            levels: vec![0, 1, 2],
        },
        Term::AbsDiff { attr: 1, pow: 1.0 },
        Term::AbsDiff { attr: 1, pow: 2.0 },
        Term::EdgeCov { cov: 0 },
        Term::Triangle,
        Term::GwPartners {
            decay: 0.5,
            partners: Partners::Edgewise,
        },
        Term::GwPartners {
            decay: 0.25,
            partners: Partners::Dyadwise,
        },
        Term::GwPartners {
            decay: 0.7,
            partners: Partners::Nonedgewise,
        },
        Term::Isolates,
    ];
    if directed {
        t.extend([
            Term::Mutual,
            Term::NodeCov {
                attr: 1,
                endpoint: Endpoint::Tail,
            },
            Term::NodeCov {
                attr: 1,
                endpoint: Endpoint::Head,
            },
            Term::NodeFactor {
                attr: 0,
                endpoint: Endpoint::Tail,
                levels: vec![0, 2],
            },
            Term::NodeFactor {
                attr: 0,
                endpoint: Endpoint::Head,
                levels: vec![1],
            },
            Term::TransitiveTriple,
            Term::CyclicTriple,
            Term::Star {
                k: 2,
                kind: DegreeKind::In,
            },
            Term::Star {
                k: 3,
                kind: DegreeKind::Out,
            },
            Term::GwDegree {
                decay: 0.3,
                kind: DegreeKind::In,
            },
            Term::GwDegree {
                decay: 1.1,
                kind: DegreeKind::Out,
            },
            Term::Degree {
                k: 2,
                kind: DegreeKind::In,
            },
            Term::Degree {
                k: 0,
                kind: DegreeKind::Out,
            },
        ]);
    } else {
        t.extend([
            Term::Star {
                k: 2,
                kind: DegreeKind::Total,
            },
            Term::Star {
                k: 3,
                kind: DegreeKind::Total,
            },
            Term::GwDegree {
                decay: 0.3,
                kind: DegreeKind::Total,
            },
            Term::Degree {
                k: 1,
                kind: DegreeKind::Total,
            },
            Term::Degree {
                k: 3,
                kind: DegreeKind::Total,
            },
        ]);
    }
    t
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// The definition of a change statistic: the difference between the
/// statistics of the network with the tie and without it, every other tie
/// as it is. Checked for every term on random networks, from states where
/// the dyad is present and where it is absent.
#[test]
fn change_statistics_match_differences_of_whole_network_statistics() {
    for case in 0..24u64 {
        let mut s = Stream::new(101, case);
        let directed = case % 2 == 1;
        let n = 5 + (case % 5) as usize;
        let density = 0.15 + 0.1 * (case % 4) as f64;
        let y0 = random_state(&mut s, n, directed, density);
        let cov = random_covariates(&mut s, n);
        let terms = Terms::new(all_terms(directed), n, directed, cov).unwrap();
        let p = terms.width();
        let mut delta = vec![0.0; p];
        for i in 0..n {
            for j in 0..n {
                if i == j || (!directed && j < i) {
                    continue;
                }
                let mut y = y0.clone();
                // From the current state.
                terms.change(&y, i, j, &mut delta);
                let before = terms.stats(&y);
                let now_present = y.toggle(i, j);
                let after = terms.stats(&y);
                for k in 0..p {
                    let diff = if now_present {
                        after[k] - before[k]
                    } else {
                        before[k] - after[k]
                    };
                    assert!(
                        close(delta[k], diff),
                        "case {case} dyad ({i}, {j}) term {:?}: change {} vs difference {}",
                        terms.terms()[term_index(&terms, k)],
                        delta[k],
                        diff
                    );
                }
                // From the toggled state the change statistic is the same
                // vector (it is defined on the rest of the network).
                let mut delta2 = vec![0.0; p];
                terms.change(&y, i, j, &mut delta2);
                for k in 0..p {
                    assert!(
                        close(delta[k], delta2[k]),
                        "case {case} dyad ({i}, {j}) stat {k}"
                    );
                }
            }
        }
    }
}

fn term_index(terms: &Terms, stat: usize) -> usize {
    let mut acc = 0;
    for (t, term) in terms.terms().iter().enumerate() {
        acc += term.width();
        if stat < acc {
            return t;
        }
    }
    unreachable!()
}

/// Relabeling the nodes (and permuting covariates the same way) leaves the
/// statistics unchanged.
#[test]
fn statistics_are_invariant_to_relabeling() {
    for case in 0..8u64 {
        let mut s = Stream::new(202, case);
        let directed = case % 2 == 0;
        let n = 7;
        let y = random_state(&mut s, n, directed, 0.3);
        let cov = random_covariates(&mut s, n);
        let mut perm = vec![0; n];
        s.permutation(&mut perm);
        let mut yp = NetState::empty(n, directed);
        for &(a, b) in y.ties() {
            yp.toggle(perm[a as usize], perm[b as usize]);
        }
        let mut covp = cov.clone();
        for x in covp.numeric.iter_mut() {
            let old = x.clone();
            for v in 0..n {
                x[perm[v]] = old[v];
            }
        }
        for a in covp.categorical.iter_mut() {
            let old = a.clone();
            for v in 0..n {
                a[perm[v]] = old[v];
            }
        }
        let mut w = Dense::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let v = cov.dyadic[0].get(i, j);
                // Undirected edgecov reads the (min, max) cell, so the
                // covariate must be symmetric for relabeling to commute.
                w.set(
                    perm[i],
                    perm[j],
                    if directed {
                        v
                    } else {
                        v + cov.dyadic[0].get(j, i)
                    },
                );
                if !directed {
                    w.set(perm[j], perm[i], v + cov.dyadic[0].get(j, i));
                }
            }
        }
        let mut cov_sym = cov.clone();
        if !directed {
            let mut ws = Dense::zeros(n, n);
            for i in 0..n {
                for j in 0..n {
                    ws.set(i, j, cov.dyadic[0].get(i, j) + cov.dyadic[0].get(j, i));
                }
            }
            cov_sym.dyadic = vec![ws];
        }
        covp.dyadic = vec![w];
        let t1 = Terms::new(all_terms(directed), n, directed, cov_sym).unwrap();
        let t2 = Terms::new(all_terms(directed), n, directed, covp).unwrap();
        let (a, b) = (t1.stats(&y), t2.stats(&yp));
        for k in 0..a.len() {
            assert!(
                close(a[k], b[k]),
                "case {case} stat {k}: {} vs {}",
                a[k],
                b[k]
            );
        }
    }
}

/// The sampler never touches a fixed dyad, and its running statistics stay
/// equal to the statistics of the network it holds.
#[test]
fn sampler_respects_constraints_and_tracks_statistics() {
    let mut s = Stream::new(303, 0);
    let n = 8;
    for directed in [false, true] {
        let y = random_state(&mut s, n, directed, 0.3);
        let mut fixed = vec![0; n * n];
        for i in 0..n {
            for j in 0..n {
                if s.uniform() < 0.3 {
                    fixed[i * n + j] = 1;
                }
            }
        }
        let c = Constraints::new(n, directed, &fixed).unwrap();
        let cov = random_covariates(&mut s, n);
        let terms = Terms::new(
            vec![
                Term::Edges,
                Term::Triangle,
                Term::GwPartners {
                    decay: 0.5,
                    partners: Partners::Edgewise,
                },
                Term::NodeCov {
                    attr: 0,
                    endpoint: Endpoint::Both,
                },
            ],
            n,
            directed,
            cov,
        )
        .unwrap();
        let model = Model::new(y.clone(), terms, Some(c.clone())).unwrap();
        for proposal in [Proposal::Tnt, Proposal::Random] {
            let sample = engine::ergm::simulate(
                &model,
                &[-0.5, 0.3, 0.2, 0.1],
                &SampleOptions {
                    burnin: 200,
                    interval: 7,
                    nsim: 30,
                    chains: 2,
                    proposal,
                    seed: 5,
                    keep_networks: true,
                    ..Default::default()
                },
            )
            .unwrap();
            let nets = sample.networks.as_ref().unwrap();
            for (r, net) in nets.iter().enumerate() {
                let mut z = NetState::empty(n, directed);
                for &(a, b) in net {
                    z.toggle(a as usize, b as usize);
                }
                for i in 0..n {
                    for j in 0..n {
                        if i != j && c.is_fixed(i, j) {
                            assert_eq!(z.has(i, j), y.has(i, j), "fixed dyad ({i}, {j}) changed");
                        }
                    }
                }
                let stats = model.terms().stats(&z);
                for k in 0..stats.len() {
                    assert!(close(stats[k], sample.row(r)[k]), "row {r} stat {k}");
                }
            }
        }
    }
}
