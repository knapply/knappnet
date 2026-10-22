//! Randomized property tests for density, components, fragmentation,
//! closeness, harmonic centrality, eigenvector centrality, and
//! centralization.

extern crate knappnet as engine;

use engine::conventions::Convention;
use engine::metrics::closeness::{closeness, harmonic_centrality, ClosenessOptions, Unreachable};
use engine::metrics::degree::DegreeMode;
use engine::metrics::eigenvector::{eigenvector_centrality, EigenScale, EigenvectorOptions};
use engine::metrics::topography::{
    centralization, components, density, fragmentation, simplify, CentralizationMeasure,
    CentralizationOptions, ComponentKind, ParallelWeight,
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
}

struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    directed: bool,
}

fn random_case(rng: &mut Rng) -> Case {
    let n = 1 + rng.below(30) as usize;
    let m = rng.below((n * n / 3 + 1) as u64) as usize;
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    for _ in 0..m {
        src.push(rng.below(n as u64) as i64);
        tgt.push(rng.below(n as u64) as i64);
    }
    Case {
        n,
        src,
        tgt,
        directed: rng.below(2) == 0,
    }
}

fn graph(c: &Case, w: Option<&[f64]>) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, w, c.directed).unwrap()
}

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn close(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| same(*x, *y))
}

const CASES: u64 = 300;
const RULES: [Unreachable; 4] = [
    Unreachable::Ignore,
    Unreachable::WassermanFaust,
    Unreachable::Zero,
    Unreachable::MaxPlusOne,
];
const MODES: [DegreeMode; 3] = [DegreeMode::Out, DegreeMode::In, DegreeMode::All];

fn copts(mode: DegreeMode, unreachable: Unreachable) -> ClosenessOptions {
    ClosenessOptions {
        mode,
        unreachable,
        normalized: true,
        ..Default::default()
    }
}

/// Reachability by breadth-first search, independent of the engine.
fn reach(c: &Case) -> Vec<Vec<bool>> {
    let mut adj = vec![Vec::new(); c.n];
    for (&s, &t) in c.src.iter().zip(&c.tgt) {
        adj[s as usize].push(t as usize);
        if !c.directed {
            adj[t as usize].push(s as usize);
        }
    }
    (0..c.n)
        .map(|s| {
            let mut seen = vec![false; c.n];
            seen[s] = true;
            let mut queue = vec![s];
            while let Some(v) = queue.pop() {
                for &w in &adj[v] {
                    if !seen[w] {
                        seen[w] = true;
                        queue.push(w);
                    }
                }
            }
            seen
        })
        .collect()
}

#[test]
fn relabeling_nodes_permutes_results() {
    let mut rng = Rng(11);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let mut perm: Vec<usize> = (0..c.n).collect();
        for i in (1..c.n).rev() {
            perm.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let p = Case {
            n: c.n,
            src: c.src.iter().map(|&v| perm[v as usize] as i64).collect(),
            tgt: c.tgt.iter().map(|&v| perm[v as usize] as i64).collect(),
            directed: c.directed,
        };
        let (g, h) = (graph(&c, None), graph(&p, None));
        let mapped = |x: &[f64]| -> Vec<f64> { (0..c.n).map(|v| x[perm[v]]).collect() };
        for mode in MODES {
            for rule in RULES {
                let o = copts(mode, rule);
                let a = closeness(&g, o).unwrap();
                assert!(close(&a, &mapped(&closeness(&h, o).unwrap())));
            }
            let o = copts(mode, Unreachable::Ignore);
            let a = harmonic_centrality(&g, o).unwrap();
            assert!(close(&a, &mapped(&harmonic_centrality(&h, o).unwrap())));
        }
        // Undirected graphs always converge, connected or not; the result
        // on a disconnected graph is defined component by component and
        // must not depend on node order either.
        let o = EigenvectorOptions {
            mode: DegreeMode::All,
            ..Default::default()
        };
        let e = eigenvector_centrality(&g, o).unwrap();
        let f = eigenvector_centrality(&h, o).unwrap();
        let fm = mapped(&f.values);
        assert!(e.values.iter().zip(&fm).all(|(x, y)| (x - y).abs() < 1e-8));
        if let (Ok(e), Ok(f)) = (
            eigenvector_centrality(&g, Default::default()),
            eigenvector_centrality(&h, Default::default()),
        ) {
            if e.connected {
                let fm = mapped(&f.values);
                assert!(e.values.iter().zip(&fm).all(|(x, y)| (x - y).abs() < 1e-8));
            }
        }
        for kind in [ComponentKind::Weak, ComponentKind::Strong] {
            let a = components(&g, kind).membership;
            let b = components(&h, kind).membership;
            for u in 0..c.n {
                for v in 0..c.n {
                    assert_eq!(a[u] == a[v], b[perm[u]] == b[perm[v]]);
                }
            }
            assert!(same(fragmentation(&g, kind), fragmentation(&h, kind)));
        }
        assert!(same(density(&g, false), density(&h, false)));
    }
}

#[test]
fn components_match_reachability() {
    let mut rng = Rng(12);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = graph(&c, None);
        let r = reach(&c);
        let strong = components(&g, ComponentKind::Strong);
        let weak = components(&g, ComponentKind::Weak);
        let sym = Case {
            n: c.n,
            src: c.src.clone(),
            tgt: c.tgt.clone(),
            directed: false,
        };
        let rs = reach(&sym);
        let mut ordered_apart = 0usize;
        for u in 0..c.n {
            for v in 0..c.n {
                assert_eq!(
                    strong.membership[u] == strong.membership[v],
                    r[u][v] && r[v][u]
                );
                assert_eq!(weak.membership[u] == weak.membership[v], rs[u][v]);
                if u != v && !rs[u][v] {
                    ordered_apart += 1;
                }
            }
        }
        assert_eq!(weak.sizes.iter().sum::<usize>(), c.n);
        assert!(strong.count() >= weak.count());
        if c.n > 1 {
            let f = fragmentation(&g, ComponentKind::Weak);
            assert!(same(f, ordered_apart as f64 / (c.n * (c.n - 1)) as f64));
        }
    }
}

#[test]
fn results_do_not_depend_on_thread_count() {
    let mut rng = Rng(13);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = graph(&c, None);
        for rule in RULES {
            let one = ClosenessOptions {
                threads: 1,
                ..copts(DegreeMode::Out, rule)
            };
            let many = ClosenessOptions { threads: 4, ..one };
            let (a, b) = (closeness(&g, one).unwrap(), closeness(&g, many).unwrap());
            assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
        }
    }
}

#[test]
fn weights_behave_as_lengths() {
    let mut rng = Rng(14);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let plain = graph(&c, None);
        let k = (1 + rng.below(5)) as f64 * 0.5;
        let uniform = graph(&c, Some(&vec![k; c.src.len()]));
        let w: Vec<f64> = (0..c.src.len())
            .map(|_| (1 + rng.below(6)) as f64)
            .collect();
        let weighted = graph(&c, Some(&w));
        let inverse: Vec<f64> = w.iter().map(|x| 1.0 / x).collect();
        let inverted = graph(&c, Some(&inverse));
        for mode in MODES {
            let base = copts(mode, Unreachable::Ignore);
            let wt = ClosenessOptions {
                weighted: true,
                normalized: false,
                ..base
            };
            // Uniform lengths k scale closeness and harmonic centrality by 1/k.
            let a: Vec<f64> = closeness(
                &plain,
                ClosenessOptions {
                    normalized: false,
                    ..base
                },
            )
            .unwrap()
            .iter()
            .map(|x| x / k)
            .collect();
            assert!(close(&a, &closeness(&uniform, wt).unwrap()));
            let h: Vec<f64> = harmonic_centrality(&plain, base)
                .unwrap()
                .iter()
                .map(|x| x / k)
                .collect();
            let hw = ClosenessOptions {
                normalized: true,
                ..wt
            };
            assert!(close(&h, &harmonic_centrality(&uniform, hw).unwrap()));
            // Every algorithm agrees, and strengths w are distances 1/w.
            let run = |g: &Graph, algorithm, weights_as| {
                let o = ClosenessOptions {
                    algorithm,
                    weights_as,
                    ..wt
                };
                (closeness(g, o).unwrap(), harmonic_centrality(g, o).unwrap())
            };
            let dial = run(&weighted, PathAlgorithm::Dial, WeightKind::Distance);
            let heap = run(&weighted, PathAlgorithm::Dijkstra, WeightKind::Distance);
            assert!(close(&dial.0, &heap.0) && close(&dial.1, &heap.1));
            let s = run(&weighted, PathAlgorithm::Auto, WeightKind::Strength);
            let d = run(&inverted, PathAlgorithm::Auto, WeightKind::Distance);
            assert!(close(&s.0, &d.0) && close(&s.1, &d.1));
        }
    }
}

#[test]
fn closeness_rules_agree_when_strongly_connected() {
    let mut rng = Rng(15);
    let mut checked = 0;
    for _ in 0..3 * CASES {
        let c = random_case(&mut rng);
        let g = graph(&c, None);
        if c.n < 2 || components(&g, ComponentKind::Strong).count() != 1 {
            continue;
        }
        checked += 1;
        let a = closeness(&g, copts(DegreeMode::Out, Unreachable::Ignore)).unwrap();
        for rule in RULES {
            assert!(close(
                &a,
                &closeness(&g, copts(DegreeMode::Out, rule)).unwrap()
            ));
        }
    }
    assert!(checked > 20);
}

#[test]
fn eigenvector_satisfies_its_equation() {
    let mut rng = Rng(16);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let w: Vec<f64> = (0..c.src.len())
            .map(|_| 0.5 + rng.below(8) as f64 / 4.0)
            .collect();
        let g = graph(&c, Some(&w));
        for mode in MODES {
            let o = EigenvectorOptions {
                mode,
                weighted: true,
                scale: EigenScale::Unit,
                ..Default::default()
            };
            let symmetric = !c.directed || mode == DegreeMode::All;
            let strongly = components(&g, ComponentKind::Strong).count() == 1;
            let e = match eigenvector_centrality(&g, o) {
                Ok(e) => e,
                // Power iteration may fail only when the largest eigenvalue
                // can be shared by strongly connected components.
                Err(err) => {
                    assert!(!symmetric && !strongly, "{err}");
                    continue;
                }
            };
            let norm: f64 = e.values.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!((norm - 1.0).abs() < 1e-12 || c.n == 0);
            assert!(e.values.iter().all(|&x| x >= 0.0));
            if e.eigenvalue == 0.0 {
                continue;
            }
            // Dense check of M x = lambda x, where row i of M sums what
            // node i receives in this mode (loops twice when undirected).
            let n = c.n;
            let mut m = vec![vec![0.0; n]; n];
            for ((&s, &t), &x) in c.src.iter().zip(&c.tgt).zip(&w) {
                let (s, t) = (s as usize, t as usize);
                let undirected = !c.directed || mode == DegreeMode::All;
                if s == t {
                    m[s][s] += if undirected { 2.0 * x } else { x };
                } else if undirected {
                    m[s][t] += x;
                    m[t][s] += x;
                } else if mode == DegreeMode::In {
                    m[t][s] += x;
                } else {
                    m[s][t] += x;
                }
            }
            for (i, row) in m.iter().enumerate() {
                let lhs: f64 = row.iter().zip(&e.values).map(|(a, b)| a * b).sum();
                assert!(
                    (lhs - e.eigenvalue * e.values[i]).abs() < 1e-7 * (1.0 + e.eigenvalue),
                    "row {i}: {lhs} vs {}",
                    e.eigenvalue * e.values[i]
                );
            }
        }
    }
}

#[test]
fn centralization_is_a_proportion_on_simple_graphs() {
    let mut rng = Rng(17);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        if c.n < 3 {
            continue;
        }
        let g = simplify(&graph(&c, None), false, ParallelWeight::Max).unwrap();
        let mut cases = vec![(CentralizationMeasure::Betweenness, DegreeMode::All)];
        if c.directed {
            cases.push((CentralizationMeasure::Degree, DegreeMode::In));
            cases.push((CentralizationMeasure::Degree, DegreeMode::Out));
        } else {
            cases.push((CentralizationMeasure::Degree, DegreeMode::All));
            if components(&g, ComponentKind::Weak).count() == 1 {
                cases.push((CentralizationMeasure::Closeness, DegreeMode::All));
            }
        }
        for (measure, mode) in cases {
            for convention in [Convention::Igraph, Convention::Sna, Convention::Ucinet] {
                let opts = CentralizationOptions {
                    measure,
                    mode,
                    loops: false,
                    convention,
                    ..Default::default()
                };
                let x = centralization(&g, opts).unwrap();
                assert!(
                    (-1e-12..=1.0 + 1e-12).contains(&x),
                    "{measure:?} {mode:?} {convention:?}: {x}"
                );
            }
        }
    }
}
