//! Randomized property tests for ego networks and two-mode networks.
//!
//! Each property compares the engine with a direct evaluation written here
//! (neighborhoods by repeated expansion, ties checked pair by pair,
//! projections as incidence-matrix products) or with an invariant: thread
//! count and node relabeling do not change results, and Borgatti and
//! Everett's betweenness maximum is attained and never exceeded.

extern crate knappnet as engine;

use std::collections::{BTreeMap, BTreeSet};

use engine::metrics::betweenness::{betweenness, BetweennessOptions};
use engine::metrics::bipartite::{
    betweenness_max, check_types, project, twomode_betweenness, twomode_closeness, twomode_degree,
    twomode_density, twomode_eigenvector, ProjectionOptions, ProjectionWeight,
    TwoModeClosenessOptions,
};
use engine::metrics::closeness::Unreachable;
use engine::metrics::degree::DegreeMode;
use engine::metrics::ego::{
    ego_summary, neighborhood, neighborhood_sizes, neighborhoods, EgoOptions, NeighborhoodOptions,
};
use engine::rng::Stream;
use engine::Graph;

const CASES: u64 = 200;
const MODES: [DegreeMode; 3] = [DegreeMode::Out, DegreeMode::In, DegreeMode::All];

struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    weights: Vec<f64>,
    directed: bool,
}

impl Case {
    fn graph(&self, weighted: bool) -> Graph {
        Graph::from_edges(
            self.n,
            &self.src,
            &self.tgt,
            weighted.then_some(&self.weights[..]),
            self.directed,
        )
        .unwrap()
    }

    /// The same network with nodes renamed by `perm` (old v becomes
    /// perm[v]).
    fn relabeled(&self, perm: &[usize]) -> Case {
        Case {
            n: self.n,
            src: self.src.iter().map(|&v| perm[v as usize] as i64).collect(),
            tgt: self.tgt.iter().map(|&v| perm[v as usize] as i64).collect(),
            weights: self.weights.clone(),
            directed: self.directed,
        }
    }
}

/// A random multigraph with self-loops.
fn random_case(rng: &mut Stream) -> Case {
    let n = 1 + rng.below(25) as usize;
    let m = rng.below((n * n / 3 + 1) as u64) as usize;
    let mut c = Case {
        n,
        src: Vec::new(),
        tgt: Vec::new(),
        weights: Vec::new(),
        directed: rng.below(2) == 0,
    };
    for _ in 0..m {
        c.src.push(rng.below(n as u64) as i64);
        c.tgt.push(rng.below(n as u64) as i64);
        c.weights.push(1.0 + rng.below(4) as f64);
    }
    c
}

/// A random two-mode multigraph with interleaved modes, random edge
/// orientation, repeated edges, and integer weights.
fn random_twomode(rng: &mut Stream) -> (Case, Vec<bool>) {
    let n = 2 + rng.below(20) as usize;
    let types: Vec<bool> = (0..n).map(|_| rng.below(3) == 0).collect();
    let first: Vec<usize> = (0..n).filter(|&v| !types[v]).collect();
    let second: Vec<usize> = (0..n).filter(|&v| types[v]).collect();
    let mut c = Case {
        n,
        src: Vec::new(),
        tgt: Vec::new(),
        weights: Vec::new(),
        directed: rng.below(2) == 0,
    };
    if !first.is_empty() && !second.is_empty() {
        let m = rng.below((first.len() * second.len()) as u64 + 1) as usize;
        for _ in 0..m {
            let a = first[rng.below(first.len() as u64) as usize] as i64;
            let b = second[rng.below(second.len() as u64) as usize] as i64;
            let (s, t) = if rng.below(2) == 0 { (a, b) } else { (b, a) };
            c.src.push(s);
            c.tgt.push(t);
            c.weights.push(1.0 + rng.below(5) as f64);
        }
    }
    (c, types)
}

fn permutation(rng: &mut Stream, n: usize) -> Vec<usize> {
    let mut p = vec![0; n];
    rng.permutation(&mut p);
    p
}

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn close(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| same(*x, *y))
}

fn permuted(values: &[f64], perm: &[usize]) -> Vec<f64> {
    let mut out = vec![0.0; values.len()];
    for (v, &x) in values.iter().enumerate() {
        out[perm[v]] = x;
    }
    out
}

/// Nodes within `order` steps of ego by repeated expansion.
fn naive_within(c: &Case, ego: usize, order: usize, mode: DegreeMode) -> BTreeMap<usize, usize> {
    let mut steps = BTreeMap::new();
    steps.insert(ego, 0);
    let mut frontier = vec![ego];
    for k in 1..=order {
        let mut next = Vec::new();
        for &v in &frontier {
            for (&s, &t) in c.src.iter().zip(&c.tgt) {
                let (s, t) = (s as usize, t as usize);
                let forward = !c.directed || mode != DegreeMode::In;
                let backward = !c.directed || mode != DegreeMode::Out;
                for (from, to, ok) in [(s, t, forward), (t, s, backward)] {
                    if ok && from == v && !steps.contains_key(&to) {
                        steps.insert(to, k);
                        next.push(to);
                    }
                }
            }
        }
        frontier = next;
    }
    steps
}

/// Distinct arcs between distinct nodes (both directions if undirected).
fn arc_set(c: &Case) -> BTreeSet<(usize, usize)> {
    let mut arcs = BTreeSet::new();
    for (&s, &t) in c.src.iter().zip(&c.tgt) {
        if s != t {
            arcs.insert((s as usize, t as usize));
            if !c.directed {
                arcs.insert((t as usize, s as usize));
            }
        }
    }
    arcs
}

#[test]
fn neighborhoods_match_repeated_expansion() {
    for case in 0..CASES {
        let mut rng = Stream::new(0xE60, case);
        let c = random_case(&mut rng);
        let g = c.graph(false);
        for mode in MODES {
            for order in 0..3 {
                for mindist in 0..2 {
                    let opts = NeighborhoodOptions {
                        order,
                        mode,
                        mindist,
                        threads: 1 + (case % 3) as usize,
                    };
                    let egos: Vec<usize> = (0..c.n).collect();
                    let all = neighborhoods(&g, &egos, opts).unwrap();
                    let sizes = neighborhood_sizes(&g, None, opts).unwrap();
                    for ego in 0..c.n {
                        let want: Vec<u32> = naive_within(&c, ego, order, mode)
                            .into_iter()
                            .filter(|&(_, k)| k >= mindist)
                            .map(|(v, _)| v as u32)
                            .collect();
                        assert_eq!(all[ego], want, "case {case} ego {ego}");
                        assert_eq!(sizes[ego], want.len());
                        assert_eq!(neighborhood(&g, ego, opts).unwrap(), want);
                    }
                }
            }
        }
    }
}

#[test]
fn ego_summary_matches_direct_evaluation() {
    for case in 0..CASES {
        let mut rng = Stream::new(0xE61, case);
        let c = random_case(&mut rng);
        let g = c.graph(true);
        let arcs = arc_set(&c);
        for mode in MODES {
            for order in 1..3 {
                let opts = EgoOptions {
                    order,
                    mode,
                    betweenness: true,
                    threads: 1 + (case % 4) as usize,
                };
                let s = ego_summary(&g, opts).unwrap();
                for ego in 0..c.n {
                    let members = naive_within(&c, ego, order, mode);
                    let k = members.len() - 1;
                    let among: Vec<(usize, usize)> = arcs
                        .iter()
                        .copied()
                        .filter(|(a, b)| {
                            members.contains_key(a)
                                && members.contains_key(b)
                                && *a != ego
                                && *b != ego
                        })
                        .collect();
                    let ties = if c.directed {
                        among.len()
                    } else {
                        among.len() / 2
                    };
                    let pairs = if c.directed {
                        k * k.saturating_sub(1)
                    } else {
                        k * k.saturating_sub(1) / 2
                    };
                    assert_eq!(s.size[ego], k as f64);
                    assert_eq!(s.ties[ego], ties as f64);
                    assert_eq!(s.pairs[ego], pairs as f64);
                    if pairs == 0 {
                        assert!(s.density[ego].is_nan());
                    } else {
                        assert_eq!(s.density[ego], ties as f64 / pairs as f64);
                    }

                    // Components by labels spread to a fixed point.
                    let alters: Vec<usize> =
                        members.keys().copied().filter(|&v| v != ego).collect();
                    let mut label: BTreeMap<usize, usize> =
                        alters.iter().map(|&v| (v, v)).collect();
                    loop {
                        let mut changed = false;
                        for &(a, b) in &among {
                            let m = label[&a].min(label[&b]);
                            for v in [a, b] {
                                if label[&v] != m {
                                    label.insert(v, m);
                                    changed = true;
                                }
                            }
                        }
                        if !changed {
                            break;
                        }
                    }
                    let comps: BTreeSet<usize> = label.values().copied().collect();
                    assert_eq!(s.components[ego], comps.len() as f64, "case {case}");

                    // Ego betweenness is the betweenness of ego in the
                    // induced simple ego network.
                    let nodes: Vec<usize> = members.keys().copied().collect();
                    let local = |v: usize| nodes.iter().position(|&u| u == v).unwrap() as i64;
                    let sub: Vec<(i64, i64)> = arcs
                        .iter()
                        .filter(|(a, b)| {
                            members.contains_key(a)
                                && members.contains_key(b)
                                && (c.directed || a < b)
                        })
                        .map(|&(a, b)| (local(a), local(b)))
                        .collect();
                    let (ss, tt): (Vec<i64>, Vec<i64>) = sub.into_iter().unzip();
                    let h = Graph::from_edges(nodes.len(), &ss, &tt, None, c.directed).unwrap();
                    let b = betweenness(&h, BetweennessOptions::default()).unwrap();
                    let got = s.betweenness.as_ref().unwrap()[ego];
                    assert!(same(got, b[local(ego) as usize]), "case {case} ego {ego}");
                }
            }
        }
    }
}

#[test]
fn ego_summary_is_invariant_to_threads_and_relabeling() {
    for case in 0..CASES {
        let mut rng = Stream::new(0xE62, case);
        let c = random_case(&mut rng);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabeled(&perm);
        for mode in MODES {
            let opts = |threads| EgoOptions {
                order: 1,
                mode,
                betweenness: true,
                threads,
            };
            let one = ego_summary(&c.graph(false), opts(1)).unwrap();
            let many = ego_summary(&c.graph(false), opts(5)).unwrap();
            let moved = ego_summary(&r.graph(false), opts(3)).unwrap();
            assert_eq!(one.size, many.size);
            assert_eq!(one.ties, many.ties);
            assert_eq!(one.components, many.components);
            assert_eq!(
                one.betweenness.as_ref().unwrap(),
                many.betweenness.as_ref().unwrap()
            );
            assert_eq!(permuted(&one.size, &perm), moved.size);
            assert_eq!(permuted(&one.ties, &perm), moved.ties);
            assert_eq!(permuted(&one.components, &perm), moved.components);
            assert!(close(&permuted(&one.density, &perm), &moved.density));
            assert!(close(
                &permuted(one.betweenness.as_ref().unwrap(), &perm),
                moved.betweenness.as_ref().unwrap()
            ));
        }
    }
}

/// Incidence-matrix projection: weights of every pair (i <= j) with a
/// shared partner, from the binary or summed-weight matrix.
fn naive_projection(
    c: &Case,
    types: &[bool],
    second: bool,
    weight: ProjectionWeight,
    loops: bool,
) -> Vec<(u32, u32, f64)> {
    let mut cell: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for (e, (&s, &t)) in c.src.iter().zip(&c.tgt).enumerate() {
        let (s, t) = (s as usize, t as usize);
        let (node, partner) = if types[s] == second { (s, t) } else { (t, s) };
        *cell.entry((node, partner)).or_insert(0.0) += c.weights[e];
    }
    let mut partners_of: BTreeMap<usize, usize> = BTreeMap::new();
    for &(_, p) in cell.keys() {
        *partners_of.entry(p).or_insert(0) += 1;
    }
    let mut out = Vec::new();
    let nodes: Vec<usize> = (0..c.n).filter(|&v| types[v] == second).collect();
    for (a, &i) in nodes.iter().enumerate() {
        for &j in &nodes[a..] {
            if i == j && !loops {
                continue;
            }
            let mut any = false;
            let mut w = 0.0;
            for (&(node, p), &wi) in cell.range((i, 0)..(i + 1, 0)) {
                debug_assert_eq!(node, i);
                if let Some(&wj) = cell.get(&(j, p)) {
                    any = true;
                    w += match weight {
                        ProjectionWeight::Count => 1.0,
                        ProjectionWeight::Newman => 1.0 / (partners_of[&p] - 1) as f64,
                        ProjectionWeight::Product => wi * wj,
                    };
                }
            }
            if any {
                out.push((i as u32, j as u32, w));
            }
        }
    }
    out
}

#[test]
fn projections_match_incidence_products() {
    for case in 0..CASES {
        let mut rng = Stream::new(0xB10, case);
        let (c, types) = random_twomode(&mut rng);
        let g = c.graph(true);
        for second in [false, true] {
            for weight in [
                ProjectionWeight::Count,
                ProjectionWeight::Newman,
                ProjectionWeight::Product,
            ] {
                for loops in [false, true] {
                    if loops && weight == ProjectionWeight::Newman {
                        continue;
                    }
                    let opts = ProjectionOptions {
                        second,
                        weight,
                        loops,
                        threads: 1 + (case % 4) as usize,
                    };
                    let p = project(&g, &types, opts).unwrap();
                    let want = naive_projection(&c, &types, second, weight, loops);
                    assert_eq!(p.from.len(), want.len(), "case {case}");
                    for (k, &(i, j, w)) in want.iter().enumerate() {
                        assert_eq!((p.from[k], p.to[k]), (i, j));
                        assert!(same(p.weights[k], w), "case {case}");
                    }
                }
            }
        }
    }
}

#[test]
fn twomode_results_are_invariant_to_threads_orientation_and_relabeling() {
    for case in 0..CASES {
        let mut rng = Stream::new(0xB11, case);
        let (c, types) = random_twomode(&mut rng);
        let perm = permutation(&mut rng, c.n);
        let r = c.relabeled(&perm);
        let mut rtypes = vec![false; c.n];
        for v in 0..c.n {
            rtypes[perm[v]] = types[v];
        }
        // Reversed edges in the other directedness.
        let flipped = Case {
            n: c.n,
            src: c.tgt.clone(),
            tgt: c.src.clone(),
            weights: c.weights.clone(),
            directed: !c.directed,
        };
        let g = c.graph(false);
        let measures = |g: &Graph, t: &[bool], threads: usize| -> Vec<Vec<f64>> {
            let wf = TwoModeClosenessOptions {
                threads,
                ..Default::default()
            };
            vec![
                vec![twomode_density(g, t).unwrap()],
                twomode_degree(g, t, true).unwrap(),
                twomode_closeness(g, t, wf).unwrap(),
                twomode_betweenness(g, t, true, threads).unwrap(),
                twomode_betweenness(g, t, false, threads).unwrap(),
            ]
        };
        let base = measures(&g, &types, 1);
        let threaded = measures(&g, &types, 4);
        let other = measures(&flipped.graph(false), &types, 2);
        let moved = measures(&r.graph(false), &rtypes, 3);
        for k in 0..base.len() {
            assert!(close(&base[k], &threaded[k]), "case {case} measure {k}");
            assert!(close(&base[k], &other[k]), "case {case} measure {k}");
            let expect = if k == 0 {
                base[0].clone()
            } else {
                permuted(&base[k], &perm)
            };
            assert!(close(&expect, &moved[k]), "case {case} measure {k}");
        }
        for values in &base[1..] {
            assert!(values.iter().all(|&x| x.is_nan() || x >= 0.0));
        }
        // Normalized degree and betweenness never exceed 1.
        for k in [1, 3] {
            assert!(base[k].iter().all(|&x| x.is_nan() || x <= 1.0 + 1e-12));
        }
    }
}

#[test]
fn connected_closeness_rules_agree_and_eigenvectors_are_singular_vectors() {
    let mut tested = 0;
    for case in 0..CASES {
        let mut rng = Stream::new(0xB12, case);
        let (c, types) = random_twomode(&mut rng);
        let g = c.graph(false);
        let e = twomode_eigenvector(&g, &types, 1e-13, 1_000_000, 1).unwrap();
        if !e.connected || c.src.is_empty() {
            continue;
        }
        tested += 1;
        let ignore = TwoModeClosenessOptions {
            unreachable: Unreachable::Ignore,
            ..Default::default()
        };
        assert!(close(
            &twomode_closeness(&g, &types, TwoModeClosenessOptions::default()).unwrap(),
            &twomode_closeness(&g, &types, ignore).unwrap()
        ));
        // lambda x_v = sum of the scores of v's distinct partners, with each
        // mode of unit length.
        let mut partners: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); c.n];
        for (&s, &t) in c.src.iter().zip(&c.tgt) {
            partners[s as usize].insert(t as usize);
            partners[t as usize].insert(s as usize);
        }
        for (v, nbrs) in partners.iter().enumerate() {
            let sum: f64 = nbrs.iter().map(|&u| e.values[u]).sum();
            assert!(
                (e.eigenvalue * e.values[v] - sum).abs() < 1e-8,
                "case {case} node {v}"
            );
            assert!(e.values[v] > 0.0);
        }
        for mode in [false, true] {
            let norm: f64 = (0..c.n)
                .filter(|&v| types[v] == mode)
                .map(|v| e.values[v] * e.values[v])
                .sum();
            assert!((norm - 1.0).abs() < 1e-12);
        }
    }
    assert!(tested > 10, "too few connected cases ({tested})");
}

/// Borgatti and Everett's maximum is the largest betweenness of a node in
/// any bipartite graph with the given mode sizes: checked against every
/// such graph for small sizes.
#[test]
fn betweenness_maximum_is_attained_and_never_exceeded() {
    for own in 1..=5usize {
        for other in 1..=4usize {
            if own * other > 12 {
                continue;
            }
            let n = own + other;
            let types: Vec<bool> = (0..n).map(|v| v >= own).collect();
            let pairs: Vec<(i64, i64)> = (0..own)
                .flat_map(|a| (own..n).map(move |b| (a as i64, b as i64)))
                .collect();
            let mut best = 0.0f64;
            for mask in 0u32..(1 << pairs.len()) {
                let (s, t): (Vec<i64>, Vec<i64>) = pairs
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| mask >> k & 1 == 1)
                    .map(|(_, &p)| p)
                    .unzip();
                let g = Graph::from_edges(n, &s, &t, None, false).unwrap();
                check_types(&g, &types).unwrap();
                let b = betweenness(&g, BetweennessOptions::default()).unwrap();
                best = best.max(b[0]);
            }
            assert_eq!(
                best,
                betweenness_max(own, other),
                "own {own}, other {other}"
            );
        }
    }
}
