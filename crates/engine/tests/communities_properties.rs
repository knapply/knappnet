//! Randomized property tests for community detection, partition
//! comparison, and dendrograms (`metrics::partition`).

extern crate knappnet as engine;

use engine::metrics::community::{louvain, modularity};
use engine::metrics::partition::{
    check_merges, compare_partitions, cut_at_groups, cut_at_max_modularity, dendrogram_modularity,
    fast_greedy, fluid_communities, infomap, label_propagation, leiden, map_equation,
    optimal_modularity, spinglass, split_join_distance, voronoi, walktrap, CompareMethod,
    LabelPropagationOptions, LeidenObjective, LeidenOptions, LpaMode, LpaVariant, SpinglassOptions,
    VoronoiMode,
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
    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    weights: Vec<f64>,
}

/// A random undirected network on 2..=max_n nodes with edge probability
/// p and real weights in [0.5, 2.5).
fn random_case(rng: &mut Rng, max_n: u64, p: f64) -> Case {
    let n = 2 + rng.below(max_n - 1) as usize;
    let (mut src, mut tgt, mut weights) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        for j in (i + 1)..n {
            if rng.uniform() < p {
                src.push(i as i64);
                tgt.push(j as i64);
                weights.push(0.5 + 2.0 * rng.uniform());
            }
        }
    }
    Case {
        n,
        src,
        tgt,
        weights,
    }
}

/// A connected random network: a random spanning tree plus random edges.
fn connected_case(rng: &mut Rng, max_n: u64, p: f64) -> Case {
    let mut c = random_case(rng, max_n, p);
    for v in 1..c.n {
        let u = rng.below(v as u64) as i64;
        c.src.push(u);
        c.tgt.push(v as i64);
        c.weights.push(0.5 + 2.0 * rng.uniform());
    }
    c
}

fn graph(c: &Case) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, Some(&c.weights), false).unwrap()
}

/// The same network with node v renamed perm[v].
fn permuted(c: &Case, perm: &[usize]) -> Graph {
    let s: Vec<i64> = c.src.iter().map(|&v| perm[v as usize] as i64).collect();
    let t: Vec<i64> = c.tgt.iter().map(|&v| perm[v as usize] as i64).collect();
    Graph::from_edges(c.n, &s, &t, Some(&c.weights), false).unwrap()
}

fn shuffle(rng: &mut Rng, n: usize) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        p.swap(i, rng.below(i as u64 + 1) as usize);
    }
    p
}

fn labels(m: &[usize]) -> Vec<i64> {
    m.iter().map(|&c| c as i64).collect()
}

/// Whether two memberships define the same partition.
fn same_partition(a: &[usize], b: &[usize]) -> bool {
    a.len() == b.len()
        && (0..a.len()).all(|i| (0..a.len()).all(|j| (a[i] == a[j]) == (b[i] == b[j])))
}

fn components_of_groups(g: &Graph, m: &[usize]) -> usize {
    let n = g.node_count();
    let adj = g.simple_undirected(false);
    let mut seen = vec![false; n];
    let mut count = 0;
    for s in 0..n {
        if seen[s] {
            continue;
        }
        count += 1;
        seen[s] = true;
        let mut stack = vec![s];
        while let Some(v) = stack.pop() {
            for &u in adj.neighbors(v) {
                let u = u as usize;
                if !seen[u] && m[u] == m[v] {
                    seen[u] = true;
                    stack.push(u);
                }
            }
        }
    }
    count
}

fn groups(m: &[usize]) -> usize {
    m.iter().copied().max().map_or(0, |x| x + 1)
}

#[test]
fn deterministic_methods_are_invariant_under_relabelling() {
    let mut rng = Rng(1);
    for _ in 0..30 {
        let c = connected_case(&mut rng, 14, 0.25);
        let g = graph(&c);
        let perm = shuffle(&mut rng, c.n);
        let h = permuted(&c, &perm);
        let back = |m: &[usize]| -> Vec<usize> { (0..c.n).map(|v| m[perm[v]]).collect() };
        // Real weights make ties (and so tie rules) practically impossible.
        let a = walktrap(&g, 4, true, 1).unwrap();
        let b = walktrap(&h, 4, true, 1).unwrap();
        assert!(same_partition(&a.membership, &back(&b.membership)));
        assert!((a.modularity - b.modularity).abs() < 1e-12);
        let a = voronoi(&g, None, true, VoronoiMode::All, None, 1).unwrap();
        let b = voronoi(&h, None, true, VoronoiMode::All, None, 1).unwrap();
        assert!(same_partition(&a.membership, &back(&b.membership)));
        let a = optimal_modularity(&g, 1.0, true).unwrap();
        let b = optimal_modularity(&h, 1.0, true).unwrap();
        assert!((a.modularity - b.modularity).abs() < 1e-12);
    }
}

#[test]
fn seeded_methods_are_reproducible() {
    let mut rng = Rng(2);
    for trial in 0..20u64 {
        let c = connected_case(&mut rng, 25, 0.15);
        let g = graph(&c);
        let lo = LeidenOptions {
            objective: LeidenObjective::Modularity,
            resolution: 1.0,
            beta: 0.01,
            n_iterations: -1,
            weighted: true,
            node_weights: None,
            initial: None,
            seed: trial,
        };
        assert_eq!(leiden(&g, &lo).unwrap(), leiden(&g, &lo).unwrap());
        assert_eq!(
            infomap(&g, true, None, 3, trial).unwrap(),
            infomap(&g, true, None, 3, trial).unwrap()
        );
        for variant in [
            LpaVariant::Dominance,
            LpaVariant::Retention,
            LpaVariant::Fast,
        ] {
            let o = LabelPropagationOptions {
                variant,
                mode: LpaMode::All,
                weighted: false,
                initial: None,
                fixed: None,
                seed: trial,
            };
            assert_eq!(
                label_propagation(&g, &o).unwrap(),
                label_propagation(&g, &o).unwrap()
            );
        }
        let so = SpinglassOptions {
            seed: trial,
            cool_fact: 0.95,
            ..Default::default()
        };
        assert_eq!(spinglass(&g, &so).unwrap(), spinglass(&g, &so).unwrap());
        assert_eq!(
            fluid_communities(&g, 2, trial).unwrap(),
            fluid_communities(&g, 2, trial).unwrap()
        );
        // Thread counts do not change the parallel methods.
        assert_eq!(
            walktrap(&g, 3, true, 1).unwrap(),
            walktrap(&g, 3, true, 4).unwrap()
        );
        assert_eq!(
            voronoi(&g, None, true, VoronoiMode::All, None, 1).unwrap(),
            voronoi(&g, None, true, VoronoiMode::All, None, 3).unwrap()
        );
    }
}

#[test]
fn optimal_modularity_bounds_every_heuristic() {
    let mut rng = Rng(3);
    for trial in 0..25u64 {
        let c = random_case(&mut rng, 16, 0.3);
        if c.src.is_empty() {
            continue;
        }
        let g = graph(&c);
        for weighted in [false, true] {
            let best = optimal_modularity(&g, 1.0, weighted).unwrap().modularity;
            let mut others = vec![
                louvain(&g, 1.0, weighted, trial).unwrap().modularity,
                fast_greedy(&g, 1.0, weighted).unwrap().modularity,
                walktrap(&g, 4, weighted, 1).unwrap().modularity,
            ];
            let lo = LeidenOptions {
                objective: LeidenObjective::Modularity,
                resolution: 1.0,
                beta: 0.01,
                n_iterations: -1,
                weighted,
                node_weights: None,
                initial: None,
                seed: trial,
            };
            others.push(leiden(&g, &lo).unwrap().quality);
            for q in others {
                assert!(q <= best + 1e-12, "{q} > optimum {best}");
            }
        }
    }
}

#[test]
fn leiden_groups_are_connected() {
    let mut rng = Rng(4);
    for trial in 0..40u64 {
        let c = random_case(&mut rng, 40, 0.08);
        let g = graph(&c);
        for (objective, resolution) in [
            (LeidenObjective::Modularity, 1.0),
            (LeidenObjective::Modularity, 2.0),
            (LeidenObjective::Cpm, 0.05),
        ] {
            let r = leiden(
                &g,
                &LeidenOptions {
                    objective,
                    resolution,
                    beta: 0.01,
                    n_iterations: 2,
                    weighted: trial % 2 == 0,
                    node_weights: None,
                    initial: None,
                    seed: trial,
                },
            )
            .unwrap();
            assert_eq!(
                components_of_groups(&g, &r.membership),
                groups(&r.membership)
            );
        }
    }
}

#[test]
fn label_propagation_ends_at_fixed_points() {
    let mut rng = Rng(5);
    for trial in 0..30u64 {
        let c = random_case(&mut rng, 30, 0.12);
        let g = graph(&c);
        let adj = g.simple_undirected(false);
        for variant in [
            LpaVariant::Dominance,
            LpaVariant::Retention,
            LpaVariant::Fast,
            LpaVariant::SemiSynchronous,
        ] {
            let m = label_propagation(
                &g,
                &LabelPropagationOptions {
                    variant,
                    mode: LpaMode::All,
                    weighted: false,
                    initial: None,
                    fixed: None,
                    seed: trial,
                },
            )
            .unwrap();
            for v in 0..c.n {
                let mut count = std::collections::HashMap::new();
                for &u in adj.neighbors(v) {
                    *count.entry(m[u as usize]).or_insert(0usize) += 1;
                }
                if let Some(&top) = count.values().max() {
                    assert_eq!(count.get(&m[v]).copied().unwrap_or(0), top, "{variant:?}");
                }
            }
        }
    }
}

#[test]
fn group_counts_and_codelengths() {
    let mut rng = Rng(6);
    for trial in 0..25u64 {
        let c = connected_case(&mut rng, 30, 0.1);
        let g = graph(&c);
        let k = 1 + (trial as usize % c.n.min(4));
        assert_eq!(groups(&fluid_communities(&g, k, trial).unwrap()), k);
        let so = SpinglassOptions {
            spins: 3,
            seed: trial,
            cool_fact: 0.95,
            ..Default::default()
        };
        assert!(groups(&spinglass(&g, &so).unwrap().membership) <= 3);
        let r = infomap(&g, true, None, 2, trial).unwrap();
        let l = map_equation(&g, &labels(&r.membership), true, None).unwrap();
        assert!((r.codelength - l).abs() < 1e-9);
        // One group is never better than the partition found.
        let one = map_equation(&g, &vec![0; c.n], true, None).unwrap();
        assert!(r.codelength <= one + 1e-9);
    }
}

#[test]
fn comparison_measures_behave_as_metrics_and_similarities() {
    let mut rng = Rng(7);
    let m = |rng: &mut Rng, n: usize, k: u64| -> Vec<i64> {
        (0..n).map(|_| rng.below(k) as i64).collect()
    };
    for _ in 0..200 {
        let n = 2 + rng.below(30) as usize;
        let (x, y, z) = (m(&mut rng, n, 4), m(&mut rng, n, 3), m(&mut rng, n, 5));
        let vi = |a: &[i64], b: &[i64]| {
            compare_partitions(a, b, CompareMethod::VariationOfInformation).unwrap()
        };
        assert!(vi(&x, &x).abs() < 1e-12);
        assert!((vi(&x, &y) - vi(&y, &x)).abs() < 1e-12);
        assert!(vi(&x, &z) <= vi(&x, &y) + vi(&y, &z) + 1e-12);
        // Renaming the groups changes nothing.
        let renamed: Vec<i64> = x.iter().map(|&c| 100 - 7 * c).collect();
        for method in [
            CompareMethod::VariationOfInformation,
            CompareMethod::NormalizedMutualInformation,
            CompareMethod::SplitJoin,
            CompareMethod::Rand,
            CompareMethod::AdjustedRand,
            CompareMethod::AdjustedMutualInformation,
        ] {
            let a = compare_partitions(&x, &y, method).unwrap();
            let b = compare_partitions(&renamed, &y, method).unwrap();
            let c = compare_partitions(&y, &x, method).unwrap();
            assert!((a - b).abs() < 1e-12 || (a.is_nan() && b.is_nan()));
            assert!(
                (a - c).abs() < 1e-9 || (a.is_nan() && c.is_nan()),
                "{method:?}"
            );
        }
        for method in [
            CompareMethod::NormalizedMutualInformation,
            CompareMethod::Rand,
            CompareMethod::AdjustedRand,
        ] {
            assert!((compare_partitions(&x, &x, method).unwrap() - 1.0).abs() < 1e-12);
        }
        let r = compare_partitions(&x, &y, CompareMethod::Rand).unwrap();
        assert!((0.0..=1.0).contains(&r));
        assert_eq!(split_join_distance(&x, &x).unwrap(), (0, 0));
    }
}

#[test]
fn dendrogram_cuts_are_consistent() {
    let mut rng = Rng(8);
    for _ in 0..30 {
        let c = random_case(&mut rng, 20, 0.2);
        let g = graph(&c);
        let w = walktrap(&g, 4, true, 1).unwrap();
        check_merges(c.n, &w.merges).unwrap();
        let comps = c.n - w.merges.len();
        for k in comps..=c.n {
            assert_eq!(groups(&cut_at_groups(c.n, &w.merges, k).unwrap()), k);
        }
        if comps > 1 {
            assert!(cut_at_groups(c.n, &w.merges, comps - 1).is_err());
        }
        let path = dendrogram_modularity(&g, &w.merges, 1.0, true).unwrap();
        let (m, q) = cut_at_max_modularity(&g, &w.merges, 1.0, true).unwrap();
        if !c.src.is_empty() {
            let direct = modularity(&g, &labels(&m), 1.0, true).unwrap();
            assert!((q - direct).abs() < 1e-12);
            assert!(path.iter().all(|&p| p <= q + 1e-12));
        }
    }
}
