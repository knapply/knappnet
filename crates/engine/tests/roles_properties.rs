//! Randomized property tests for reach, edge betweenness, information
//! centrality, bridging, hierarchy, regular and automorphic equivalence,
//! blockmodels, classical scaling, and geospatial networks.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::matrix::Dense;
use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::blockmodel::{block_model, BlockModelOptions};
use engine::metrics::brokerage::{constraint, hierarchy, StructuralHoleOptions};
use engine::metrics::closeness::{closeness, ClosenessOptions, Unreachable};
use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::{distance_rows, PathOptions};
use engine::metrics::equivalence::{structural_equivalence, EquivalenceMethod, EquivalenceOptions};
use engine::metrics::reach::{
    bridging_coefficient, edge_betweenness, information_centrality, reach, EdgeBetweennessOptions,
    ReachOptions,
};
use engine::metrics::regular::{automorphic_equivalence, regular_equivalence, RegularOptions};
use engine::metrics::scaling::{classical_mds, stress};
use engine::spatial::{
    edge_lengths, group_geometry, haversine, nearest_neighbors, within_distance, Crs, Points,
};
use engine::{Graph, PathAlgorithm, WeightKind};

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
    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
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

    /// The same network with its edges listed in a different order.
    fn shuffle_edges(&self, order: &[usize]) -> Case {
        Case {
            n: self.n,
            src: order.iter().map(|&e| self.src[e]).collect(),
            tgt: order.iter().map(|&e| self.tgt[e]).collect(),
            weights: order.iter().map(|&e| self.weights[e]).collect(),
            directed: self.directed,
        }
    }
}

fn random_case(rng: &mut Rng, max_n: usize) -> Case {
    let n = 1 + rng.below(max_n as u64) as usize;
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

/// Whether two labelings describe the same partition.
fn same_partition(a: &[u32], b: &[u32]) -> bool {
    let canon = |v: &[u32]| {
        let mut seen = std::collections::HashMap::new();
        v.iter()
            .map(|&x| {
                let len = seen.len() as u32;
                *seen.entry(x).or_insert(len)
            })
            .collect::<Vec<u32>>()
    };
    canon(a) == canon(b)
}

const CASES: u64 = 150;

#[test]
fn reach_is_monotone_in_k_and_permutes_with_labels() {
    let mut rng = Rng(0x5eed_1001);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 30);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        for weighted in [false, true] {
            let g = c.graph(weighted);
            let mut previous = vec![0.0; c.n];
            for k in [0.0, 1.0, 2.0, 3.5, f64::INFINITY] {
                let opts = ReachOptions {
                    k,
                    weighted,
                    threads: 1,
                    ..Default::default()
                };
                let a = reach(&g, opts).unwrap();
                let b = reach(&r.graph(weighted), ReachOptions { threads: 3, ..opts }).unwrap();
                for v in 0..c.n {
                    assert!(a[v] >= previous[v]);
                    assert!(a[v] <= (c.n - 1) as f64);
                    assert!(same(a[v], b[perm[v]]));
                }
                previous = a;
            }
            // k = 1 unweighted counts distinct neighbours; k = infinity
            // counts every node reached.
            let one = reach(
                &g,
                ReachOptions {
                    k: 1.0,
                    weighted: false,
                    threads: 1,
                    ..Default::default()
                },
            )
            .unwrap();
            let mut nb = vec![std::collections::HashSet::new(); c.n];
            for (&s, &t) in c.src.iter().zip(&c.tgt) {
                if s != t {
                    nb[s as usize].insert(t);
                    if !c.directed {
                        nb[t as usize].insert(s);
                    }
                }
            }
            for v in 0..c.n {
                assert_eq!(one[v], nb[v].len() as f64);
            }
        }
    }
}

#[test]
fn reach_with_in_mode_is_the_transpose_count() {
    let mut rng = Rng(0x5eed_1002);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        let g = c.graph(false);
        let k = 2.0;
        let out = reach(
            &g,
            ReachOptions {
                k,
                mode: DegreeMode::Out,
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let inn = reach(
            &g,
            ReachOptions {
                k,
                mode: DegreeMode::In,
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        // The number of ordered pairs within k is the same from either end.
        assert!(same(out.iter().sum(), inn.iter().sum()));
        if !c.directed {
            for v in 0..c.n {
                assert!(same(out[v], inn[v]));
            }
        }
    }
}

#[test]
fn proximity_prestige_lies_in_the_unit_interval() {
    // Closeness with the Wasserman-Faust rule over paths to each node is
    // proximity prestige: (I/(n-1)) / (mean distance from I).
    let mut rng = Rng(0x5eed_1003);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        let g = c.graph(false);
        let p = closeness(
            &g,
            ClosenessOptions {
                mode: DegreeMode::In,
                normalized: true,
                unreachable: Unreachable::WassermanFaust,
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let n = c.n;
        let mut rows: Vec<Vec<f64>> = Vec::new();
        distance_rows(
            &g,
            None,
            PathOptions {
                mode: DegreeMode::In,
                ..Default::default()
            },
            |_, row| rows.push(row.to_vec()),
        )
        .unwrap();
        for v in 0..n {
            assert!((0.0..=1.0).contains(&p[v]) || n == 1, "{}", p[v]);
            let dom: Vec<f64> = (0..n)
                .filter(|&j| j != v && rows[v][j].is_finite())
                .map(|j| rows[v][j])
                .collect();
            let expect = if dom.is_empty() {
                0.0
            } else {
                let size = dom.len() as f64;
                (size / (n as f64 - 1.0)) / (dom.iter().sum::<f64>() / size)
            };
            assert!(same(p[v], expect));
        }
    }
}

#[test]
fn edge_betweenness_follows_edge_order_and_conserves_flow() {
    let mut rng = Rng(0x5eed_1004);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        let m = c.src.len();
        let order = permutation(&mut rng, m);
        let shuffled = c.shuffle_edges(&order);
        for weighted in [false, true] {
            let opts = EdgeBetweennessOptions {
                weighted,
                threads: 1,
                ..Default::default()
            };
            let g = c.graph(weighted);
            let a = edge_betweenness(&g, opts).unwrap();
            let b = edge_betweenness(
                &shuffled.graph(weighted),
                EdgeBetweennessOptions { threads: 4, ..opts },
            )
            .unwrap();
            for (k, &e) in order.iter().enumerate() {
                assert!(same(a[e], b[k]));
            }
            // Dial and Dijkstra agree on integer-valued lengths.
            if weighted {
                let integer = Case {
                    n: c.n,
                    src: c.src.clone(),
                    tgt: c.tgt.clone(),
                    weights: c.weights.iter().map(|w| w.ceil()).collect(),
                    directed: c.directed,
                };
                let gi = integer.graph(true);
                let dial = edge_betweenness(
                    &gi,
                    EdgeBetweennessOptions {
                        algorithm: PathAlgorithm::Dial,
                        weights_as: WeightKind::Distance,
                        ..opts
                    },
                )
                .unwrap();
                let dijkstra = edge_betweenness(
                    &gi,
                    EdgeBetweennessOptions {
                        algorithm: PathAlgorithm::Dijkstra,
                        weights_as: WeightKind::Distance,
                        ..opts
                    },
                )
                .unwrap();
                for e in 0..m {
                    assert!(same(dial[e], dijkstra[e]));
                }
            }
            // Flow conservation. Directed: the edge betweenness entering a
            // node sums to its betweenness plus the number of nodes that
            // reach it (the pairs ending there). Undirected: every shortest
            // path crosses one more edge than it has intermediate nodes, so
            // the edge values sum to the node values plus the number of
            // connected unordered pairs.
            let nb = betweenness(
                &g,
                BetweennessOptions {
                    weighted,
                    ..Default::default()
                },
            )
            .unwrap();
            let reached = reach(
                &g,
                ReachOptions {
                    k: f64::INFINITY,
                    mode: DegreeMode::In,
                    weighted,
                    threads: 1,
                    ..Default::default()
                },
            )
            .unwrap();
            if c.directed {
                let mut inflow = vec![0.0; c.n];
                for (e, &t) in c.tgt.iter().enumerate() {
                    inflow[t as usize] += a[e];
                }
                for v in 0..c.n {
                    assert!(same(inflow[v], nb[v] + reached[v]), "node {v}");
                }
            } else {
                let total_edges: f64 = a.iter().sum();
                let total_nodes: f64 = nb.iter().sum();
                let pairs: f64 = reached.iter().sum::<f64>() / 2.0;
                assert!(same(total_edges, total_nodes + pairs));
            }
            // Normalization is a constant factor.
            let norm = edge_betweenness(
                &g,
                EdgeBetweennessOptions {
                    normalized: true,
                    ..opts
                },
            )
            .unwrap();
            let n = c.n as f64;
            let pairs = if c.directed {
                n * (n - 1.0)
            } else {
                n * (n - 1.0) / 2.0
            };
            for e in 0..m {
                if c.n >= 2 {
                    assert!(same(norm[e], a[e] / pairs));
                } else {
                    assert!(norm[e].is_nan());
                }
            }
        }
    }
}

#[test]
fn information_centrality_is_symmetric_scaled_and_permutes() {
    let mut rng = Rng(0x5eed_1005);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 20);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        for weighted in [false, true] {
            let a = information_centrality(&c.graph(weighted), weighted).unwrap();
            let b = information_centrality(&r.graph(weighted), weighted).unwrap();
            for v in 0..c.n {
                assert!(same(a[v], b[perm[v]]));
                assert!(a[v].is_nan() || a[v] >= 0.0);
            }
            // Without weights, direction does not matter: the symmetrized
            // ties are the same either way. (With weights a directed pair
            // takes the larger direction's total, an undirected one the sum
            // of both, so the readings differ.)
            if !weighted {
                let flipped = Case {
                    n: c.n,
                    src: c.src.clone(),
                    tgt: c.tgt.clone(),
                    weights: c.weights.clone(),
                    directed: !c.directed,
                };
                let other = information_centrality(&flipped.graph(false), false).unwrap();
                for v in 0..c.n {
                    assert!(same(a[v], other[v]));
                }
            }
            // Doubling weights doubles the scores.
            if weighted {
                let doubled = Case {
                    n: c.n,
                    src: c.src.clone(),
                    tgt: c.tgt.clone(),
                    weights: c.weights.iter().map(|w| 2.0 * w).collect(),
                    directed: c.directed,
                };
                let d = information_centrality(&doubled.graph(true), true).unwrap();
                for v in 0..c.n {
                    assert!(same(d[v], 2.0 * a[v]));
                }
            }
        }
    }
}

#[test]
fn bridging_coefficient_and_hierarchy_permute_and_lie_in_range() {
    let mut rng = Rng(0x5eed_1006);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 30);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        let a = bridging_coefficient(&c.graph(false));
        let b = bridging_coefficient(&r.graph(false));
        for v in 0..c.n {
            assert!(same(a[v], b[perm[v]]));
            assert!(a[v].is_nan() || a[v] > 0.0);
        }
        for weighted in [false, true] {
            let opts = StructuralHoleOptions {
                weighted,
                threads: 1,
            };
            let h = hierarchy(&c.graph(weighted), opts).unwrap();
            let hp = hierarchy(
                &r.graph(weighted),
                StructuralHoleOptions {
                    weighted,
                    threads: 4,
                },
            )
            .unwrap();
            let cons = constraint(&c.graph(weighted), opts).unwrap();
            for v in 0..c.n {
                assert!(same(h[v], hp[perm[v]]));
                if h[v].is_nan() {
                    continue;
                }
                assert!((-1e-12..=1.0 + 1e-12).contains(&h[v]), "{}", h[v]);
                assert!(cons[v].is_finite());
            }
        }
    }
}

#[test]
fn regular_equivalence_is_regular_and_refines_the_seed() {
    let mut rng = Rng(0x5eed_1007);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        let seed: Vec<i64> = (0..c.n).map(|_| rng.below(3) as i64).collect();
        let seed_r: Vec<i64> = {
            let mut s = vec![0; c.n];
            for v in 0..c.n {
                s[perm[v]] = seed[v];
            }
            s
        };
        for weighted in [false, true] {
            let opts = RegularOptions { weighted };
            for (s, sr) in [(None, None), (Some(&seed[..]), Some(&seed_r[..]))] {
                let a = regular_equivalence(&c.graph(weighted), s, opts).unwrap();
                let b = regular_equivalence(&r.graph(weighted), sr, opts).unwrap();
                let mut bp = vec![0; c.n];
                for v in 0..c.n {
                    bp[v] = b.membership[perm[v]];
                }
                assert!(same_partition(&a.membership, &bp));
                assert_eq!(a.rounds, b.rounds);
                // Refines the seed.
                if let Some(s) = s {
                    for i in 0..c.n {
                        for j in 0..c.n {
                            if a.membership[i] == a.membership[j] {
                                assert_eq!(s[i], s[j]);
                            }
                        }
                    }
                }
                // Regularity: equivalent nodes have the same set of
                // (tie type, alter class) pairs.
                let g = c.graph(weighted);
                let w = g.edge_weights();
                let mut ties: std::collections::HashMap<(usize, usize), f64> =
                    std::collections::HashMap::new();
                for (e, (&s, &t)) in c.src.iter().zip(&c.tgt).enumerate() {
                    if s == t {
                        continue;
                    }
                    let x = if weighted { w.unwrap()[e] } else { 1.0 };
                    for (a_, b_) in if c.directed {
                        vec![(s as usize, t as usize)]
                    } else {
                        vec![(s as usize, t as usize), (t as usize, s as usize)]
                    } {
                        let v = ties.entry((a_, b_)).or_insert(0.0);
                        *v = if weighted { *v + x } else { 1.0 };
                    }
                }
                let profile = |i: usize| {
                    let mut p: Vec<(u64, u64, u32)> = (0..c.n)
                        .filter(|&k| k != i)
                        .filter_map(|k| {
                            let out = ties.get(&(i, k)).copied().unwrap_or(0.0);
                            let inc = ties.get(&(k, i)).copied().unwrap_or(0.0);
                            if out == 0.0 && inc == 0.0 {
                                None
                            } else {
                                Some((out.to_bits(), inc.to_bits(), a.membership[k]))
                            }
                        })
                        .collect();
                    p.sort_unstable();
                    p.dedup();
                    p
                };
                let profiles: Vec<_> = (0..c.n).map(profile).collect();
                for i in 0..c.n {
                    for j in 0..c.n {
                        if a.membership[i] == a.membership[j] {
                            assert_eq!(profiles[i], profiles[j]);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn automorphic_equivalence_refines_regular_and_is_refined_by_structural() {
    let mut rng = Rng(0x5eed_1008);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 14);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabel(&perm);
        let g = c.graph(false);
        let orbits = automorphic_equivalence(&g, 1_000_000).unwrap();
        let orbits_r = automorphic_equivalence(&r.graph(false), 1_000_000).unwrap();
        let mut back = vec![0; c.n];
        for v in 0..c.n {
            back[v] = orbits_r[perm[v]];
        }
        assert!(same_partition(&orbits, &back));
        // Automorphically equivalent nodes are regularly equivalent (the
        // maximal regular equivalence is coarser).
        let reg = regular_equivalence(&g, None, RegularOptions::default()).unwrap();
        for i in 0..c.n {
            for j in 0..c.n {
                if orbits[i] == orbits[j] {
                    assert_eq!(reg.membership[i], reg.membership[j]);
                }
            }
        }
        // Structurally equivalent nodes (distance 0 with loops as node
        // properties) are automorphically equivalent: swapping them is an
        // automorphism.
        let se = structural_equivalence(
            &g,
            EquivalenceOptions {
                method: EquivalenceMethod::Hamming,
                weighted: false,
                threads: 1,
            },
        )
        .unwrap();
        let has_loop: Vec<bool> = (0..c.n)
            .map(|v| {
                c.src
                    .iter()
                    .zip(&c.tgt)
                    .any(|(&s, &t)| s == t && s as usize == v)
            })
            .collect();
        for i in 0..c.n {
            for j in 0..c.n {
                let tie_ij = c.src.iter().zip(&c.tgt).any(|(&s, &t)| {
                    (s as usize == i && t as usize == j)
                        || (!c.directed && s as usize == j && t as usize == i)
                });
                let tie_ji = c.src.iter().zip(&c.tgt).any(|(&s, &t)| {
                    (s as usize == j && t as usize == i)
                        || (!c.directed && s as usize == i && t as usize == j)
                });
                if i != j && se.get(i, j) == 0.0 && has_loop[i] == has_loop[j] && tie_ij == tie_ji {
                    assert_eq!(orbits[i], orbits[j], "{i} {j}");
                }
            }
        }
    }
}

#[test]
fn block_model_densities_are_weighted_means_and_fit_is_a_correlation() {
    let mut rng = Rng(0x5eed_1009);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        let labels: Vec<i64> = (0..c.n).map(|_| 10 + rng.below(4) as i64).collect();
        for weighted in [false, true] {
            for loops in [false, true] {
                let opts = BlockModelOptions { weighted, loops };
                let bm = block_model(&c.graph(weighted), &labels, opts).unwrap();
                let k = bm.blocks();
                // Overall density is the pair-weighted mean of block densities.
                let mut total = 0.0;
                let mut pairs = 0.0;
                for r in 0..k {
                    for s in 0..k {
                        let p = bm.pairs.get(r, s);
                        if p > 0.0 {
                            assert!(same(bm.density.get(r, s), bm.sums.get(r, s) / p));
                            total += bm.sums.get(r, s);
                            pairs += p;
                        } else {
                            assert!(bm.density.get(r, s).is_nan());
                        }
                    }
                }
                if pairs > 0.0 {
                    assert!(same(bm.overall, total / pairs));
                    assert!(bm.r_squared.is_nan() || (-1e-9..=1.0 + 1e-9).contains(&bm.r_squared));
                }
                // The image at the overall density has a 1 wherever the
                // density is at least that.
                let img = bm.image(bm.overall);
                for r in 0..k {
                    for s in 0..k {
                        let d = bm.density.get(r, s);
                        let i = img.get(r, s);
                        assert!(
                            (d.is_nan() && i.is_nan())
                                || (d >= bm.overall && i == 1.0)
                                || (d < bm.overall && i == 0.0)
                        );
                    }
                }
                // Undirected: symmetric densities.
                if !c.directed {
                    for r in 0..k {
                        for s in 0..k {
                            assert!(same(bm.density.get(r, s), bm.density.get(s, r)));
                        }
                    }
                }
                // One block: density is the network density, fit undefined.
                let one = block_model(&c.graph(weighted), &vec![7; c.n], opts).unwrap();
                assert!(same(one.density.get(0, 0), one.overall));
                assert!(one.r_squared.is_nan());
            }
        }
    }
}

#[test]
fn classical_mds_recovers_euclidean_configurations() {
    let mut rng = Rng(0x5eed_100a);
    for _ in 0..60 {
        let n = 2 + rng.below(20) as usize;
        let p = 1 + rng.below(3) as usize;
        let pts: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..p).map(|_| 10.0 * rng.uniform() - 5.0).collect())
            .collect();
        let mut d = Dense::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let dist = pts[i]
                    .iter()
                    .zip(&pts[j])
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>()
                    .sqrt();
                d.set(i, j, dist);
            }
        }
        let k = (n - 1).min(p);
        let mds = classical_mds(&d, k).unwrap();
        assert_eq!(mds.eigenvalues.len(), n);
        for w in mds.eigenvalues.windows(2) {
            assert!(w[0] >= w[1] - 1e-9);
        }
        // Distances are reproduced when enough dimensions are kept.
        if p < n && mds.points.cols() == p {
            let s = stress(&d, &mds.points).unwrap();
            assert!(s < 1e-7, "stress {s}");
            assert!(mds.gof[0] > 1.0 - 1e-7 && mds.gof[1] > 1.0 - 1e-7);
        }
        // Scaling the dissimilarities scales the coordinates.
        let mut d2 = Dense::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                d2.set(i, j, 3.0 * d.get(i, j));
            }
        }
        let mds2 = classical_mds(&d2, k).unwrap();
        assert_eq!(mds2.points.cols(), mds.points.cols());
        for c in 0..mds.points.cols() {
            for i in 0..n {
                assert!(same(mds2.points.get(i, c), 3.0 * mds.points.get(i, c)));
            }
        }
        for j in 0..n {
            assert!(same(mds2.eigenvalues[j], 9.0 * mds.eigenvalues[j]));
        }
        // Stress of the exact configuration is 0 and of a collapsed one 1.
        let exact = Dense::from_row_major(n, p, pts.iter().flatten().copied().collect()).unwrap();
        assert!(stress(&d, &exact).unwrap() < 1e-12);
    }
}

#[test]
fn spatial_distances_are_metric_and_proximity_matches_brute_force() {
    let mut rng = Rng(0x5eed_100b);
    for _ in 0..40 {
        let n = 1 + rng.below(120) as usize;
        let planar = rng.below(2) == 0;
        let (x, y): (Vec<f64>, Vec<f64>) = if planar {
            (
                (0..n).map(|_| 100.0 * rng.uniform()).collect(),
                (0..n).map(|_| 100.0 * rng.uniform()).collect(),
            )
        } else {
            (
                (0..n).map(|_| 360.0 * rng.uniform() - 180.0).collect(),
                (0..n)
                    .map(|_| (2.0 * rng.uniform() - 1.0).asin().to_degrees())
                    .collect(),
            )
        };
        let crs = if planar {
            Crs::Planar
        } else {
            Crs::LonLat { radius: 1.0 }
        };
        let pts = Points::new(&x, &y, crs).unwrap();
        // Symmetry, zero on the diagonal, triangle inequality.
        for i in 0..n {
            assert_eq!(pts.distance(i, i), 0.0);
            for j in 0..n {
                assert!(same(pts.distance(i, j), pts.distance(j, i)));
                if !planar {
                    assert!(same(
                        pts.distance(i, j),
                        haversine(x[i], y[i], x[j], y[j], 1.0)
                    ));
                    assert!(pts.distance(i, j) <= std::f64::consts::PI + 1e-12);
                }
                for l in 0..n.min(8) {
                    assert!(pts.distance(i, j) <= pts.distance(i, l) + pts.distance(l, j) + 1e-9);
                }
            }
        }
        let r = if planar {
            15.0 * rng.uniform()
        } else {
            0.5 * rng.uniform()
        };
        for threads in [1, 3] {
            let w = within_distance(&pts, r, threads).unwrap();
            let mut expect = Vec::new();
            for i in 0..n {
                for j in i + 1..n {
                    if pts.distance(i, j) <= r {
                        expect.push((i as u32, j as u32));
                    }
                }
            }
            let got: Vec<(u32, u32)> = w.from.iter().copied().zip(w.to.iter().copied()).collect();
            assert_eq!(got, expect);
            for (e, &d) in w.distance.iter().enumerate() {
                assert!(same(d, pts.distance(w.from[e] as usize, w.to[e] as usize)));
            }
            let k = 1 + rng.below(4) as usize;
            let nn = nearest_neighbors(&pts, k, threads).unwrap();
            let kk = k.min(n.saturating_sub(1));
            assert_eq!(nn.from.len(), n * kk);
            for i in 0..n {
                let mut cand: Vec<(f64, u32)> = (0..n)
                    .filter(|&j| j != i)
                    .map(|j| (pts.distance(i, j), j as u32))
                    .collect();
                cand.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                for (slot, &(d, j)) in cand.iter().take(kk).enumerate() {
                    let e = i * kk + slot;
                    assert_eq!(nn.from[e], i as u32);
                    assert_eq!(nn.to[e], j);
                    assert!(same(nn.distance[e], d));
                }
            }
        }
        // Edge lengths of a random network over the points.
        let m = rng.below(3 * n as u64) as usize;
        let src: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
        let tgt: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
        let g = Graph::from_edges(n, &src, &tgt, None, true).unwrap();
        let lens = edge_lengths(&g, &pts).unwrap();
        for e in 0..m {
            assert!(same(
                lens[e],
                pts.distance(src[e] as usize, tgt[e] as usize)
            ));
        }
        // Group geometry: sizes, mean and root-mean-square distances from
        // the reported centroid, and the planar centroid as the coordinate
        // mean.
        let labels: Vec<i64> = (0..n).map(|_| rng.below(3) as i64).collect();
        let gg = group_geometry(&pts, &labels).unwrap();
        for (gi, &label) in gg.groups.iter().enumerate() {
            let members: Vec<usize> = (0..n).filter(|&v| labels[v] == label).collect();
            assert_eq!(gg.sizes[gi], members.len());
            if gg.x[gi].is_nan() {
                continue;
            }
            let cx = [gg.x[gi]];
            let cy = [gg.y[gi]];
            let mut xs: Vec<f64> = members.iter().map(|&v| x[v]).collect();
            let mut ys: Vec<f64> = members.iter().map(|&v| y[v]).collect();
            xs.push(cx[0]);
            ys.push(cy[0]);
            let joined = Points::new(&xs, &ys, crs).unwrap();
            let last = members.len();
            let dists: Vec<f64> = (0..last).map(|i| joined.distance(i, last)).collect();
            let mean = dists.iter().sum::<f64>() / last as f64;
            let rms = (dists.iter().map(|d| d * d).sum::<f64>() / last as f64).sqrt();
            assert!(same(gg.mean_distance[gi], mean));
            assert!(same(gg.standard_distance[gi], rms));
            assert!(gg.standard_distance[gi] >= gg.mean_distance[gi] - 1e-9);
            if planar {
                let mx = members.iter().map(|&v| x[v]).sum::<f64>() / last as f64;
                let my = members.iter().map(|&v| y[v]).sum::<f64>() / last as f64;
                assert!(same(gg.x[gi], mx) && same(gg.y[gi], my));
            }
        }
    }
}
