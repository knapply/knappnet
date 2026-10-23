//! Randomized property tests for layouts. Layouts are not unique, so these
//! check invariants rather than values: coordinates are finite, results
//! are identical for any thread count and reproducible from the seed,
//! components do not overlap, and weights are read consistently.

// The crate is published under the project name; `engine` keeps that name
// out of the code so renaming the project only touches Cargo.toml.
extern crate knappnet as engine;

use engine::layout::{
    circle_layout, fr_layout, grid_layout, layout_stress, star_layout, stress_layout, tree_layout,
    twomode_layout, FrApproximation, FrOptions, Layout, StressAlgorithm, StressOptions,
};
use engine::metrics::topography::{components, ComponentKind};
use engine::{Graph, WeightKind};

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
    w: Vec<f64>,
    directed: bool,
}

impl Case {
    fn graph(&self, weighted: bool) -> Graph {
        let w = if weighted { Some(&self.w[..]) } else { None };
        Graph::from_edges(self.n, &self.src, &self.tgt, w, self.directed).unwrap()
    }
}

/// Sparse random graphs with loops, parallel edges, isolates, and several
/// components, and positive integer or real weights.
fn random_case(rng: &mut Rng, max_n: u64) -> Case {
    let n = 1 + rng.below(max_n) as usize;
    let m = rng.below(2 * n as u64 + 1) as usize;
    let integer = rng.below(2) == 0;
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
        c.w.push(if integer {
            1.0 + rng.below(5) as f64
        } else {
            0.25 + rng.below(1000) as f64 / 250.0
        });
    }
    c
}

fn finite(l: &Layout, n: usize) -> bool {
    l.x.len() == n && l.y.len() == n && l.x.iter().chain(&l.y).all(|v| v.is_finite())
}

fn dist(l: &Layout, a: usize, b: usize) -> f64 {
    ((l.x[a] - l.x[b]).powi(2) + (l.y[a] - l.y[b]).powi(2)).sqrt()
}

/// Nodes of different components are at least `gap` apart.
fn separated(g: &Graph, l: &Layout, gap: f64) -> bool {
    let comp = components(g, ComponentKind::Weak).membership;
    let n = g.node_count();
    for a in 0..n {
        for b in a + 1..n {
            if comp[a] != comp[b] && dist(l, a, b) < gap * (1.0 - 1e-9) {
                return false;
            }
        }
    }
    true
}

/// Mean shortest merged edge length, the packing gap of stress layouts.
fn gap(c: &Case, weighted: bool, kind: WeightKind) -> f64 {
    let mut best: std::collections::BTreeMap<(i64, i64), f64> = Default::default();
    for e in 0..c.src.len() {
        let (a, b) = (c.src[e].min(c.tgt[e]), c.src[e].max(c.tgt[e]));
        if a == b {
            continue;
        }
        let len = match (weighted, kind) {
            (false, _) => 1.0,
            (true, WeightKind::Distance) => c.w[e],
            (true, WeightKind::Strength) => 1.0 / c.w[e],
        };
        let slot = best.entry((a, b)).or_insert(f64::INFINITY);
        *slot = slot.min(len);
    }
    if best.is_empty() {
        1.0
    } else {
        best.values().sum::<f64>() / best.len() as f64
    }
}

#[test]
fn stress_layouts_are_finite_reproducible_and_thread_independent() {
    let mut rng = Rng(11);
    for _ in 0..60 {
        let c = random_case(&mut rng, 60);
        let weighted = rng.below(2) == 0;
        let g = c.graph(weighted);
        for algorithm in [StressAlgorithm::Full, StressAlgorithm::Sparse] {
            let opts = |threads: usize| StressOptions {
                weighted,
                algorithm,
                pivots: if algorithm == StressAlgorithm::Sparse {
                    7
                } else {
                    0
                },
                seed: 5,
                threads,
                ..Default::default()
            };
            let one = stress_layout(&g, &opts(1)).unwrap();
            assert!(finite(&one, c.n));
            assert_eq!(one, stress_layout(&g, &opts(1)).unwrap());
            assert_eq!(one, stress_layout(&g, &opts(2)).unwrap());
            assert_eq!(one, stress_layout(&g, &opts(5)).unwrap());
            assert!(separated(&g, &one, gap(&c, weighted, WeightKind::Distance)));
        }
    }
}

#[test]
fn large_components_use_parallel_loops_without_changing_results() {
    // A 300-node ring with chords exercises the parallel inner loops
    // (components of 256 nodes or more), and an isolate the packing.
    let n = 301;
    let mut s: Vec<i64> = (0..300).collect();
    let mut t: Vec<i64> = (0..300).map(|i| (i + 1) % 300).collect();
    for i in (0..300).step_by(7) {
        s.push(i);
        t.push((i * 13 + 5) % 300);
    }
    let g = Graph::from_edges(n, &s, &t, None, false).unwrap();
    for algorithm in [StressAlgorithm::Full, StressAlgorithm::Sparse] {
        let opts = |threads: usize| StressOptions {
            algorithm,
            pivots: 30,
            threads,
            ..Default::default()
        };
        let one = stress_layout(&g, &opts(1)).unwrap();
        assert!(finite(&one, n));
        assert_eq!(one, stress_layout(&g, &opts(2)).unwrap());
        assert_eq!(one, stress_layout(&g, &opts(3)).unwrap());
    }
    for approximation in [FrApproximation::Exact, FrApproximation::BarnesHut] {
        let opts = |threads: usize| FrOptions {
            approximation,
            iterations: 60,
            threads,
            ..Default::default()
        };
        let one = fr_layout(&g, &opts(1)).unwrap();
        assert!(finite(&one, n));
        assert_eq!(one, fr_layout(&g, &opts(2)).unwrap());
        assert_eq!(one, fr_layout(&g, &opts(4)).unwrap());
    }
}

#[test]
fn seeds_change_random_starts_only() {
    let mut rng = Rng(12);
    let c = loop {
        let c = random_case(&mut rng, 40);
        if c.n > 10 && c.src.len() > 10 {
            break c;
        }
    };
    let g = c.graph(false);
    let fr = |seed: u64| {
        fr_layout(
            &g,
            &FrOptions {
                seed,
                iterations: 50,
                ..Default::default()
            },
        )
        .unwrap()
    };
    assert_eq!(fr(3), fr(3));
    assert_ne!(fr(3), fr(4));
}

#[test]
fn force_layouts_are_finite_reproducible_and_separated() {
    let mut rng = Rng(13);
    for _ in 0..40 {
        let c = random_case(&mut rng, 50);
        let weighted = rng.below(2) == 0;
        let g = c.graph(weighted);
        for approximation in [FrApproximation::Exact, FrApproximation::BarnesHut] {
            let opts = |threads: usize| FrOptions {
                weighted,
                approximation,
                iterations: 80,
                seed: 9,
                threads,
                ..Default::default()
            };
            let one = fr_layout(&g, &opts(1)).unwrap();
            assert!(finite(&one, c.n));
            assert_eq!(one, fr_layout(&g, &opts(3)).unwrap());
            assert!(separated(&g, &one, 1.0));
        }
    }
}

#[test]
fn strengths_are_reciprocal_distances_and_uniform_weights_scale() {
    let mut rng = Rng(14);
    for _ in 0..40 {
        let c = random_case(&mut rng, 40);
        let g = c.graph(true);
        let recip = Case {
            n: c.n,
            src: c.src.clone(),
            tgt: c.tgt.clone(),
            w: c.w.iter().map(|w| 1.0 / w).collect(),
            directed: c.directed,
        };
        let gr = recip.graph(true);
        let strength = StressOptions {
            weighted: true,
            weights_as: WeightKind::Strength,
            ..Default::default()
        };
        let distance = StressOptions {
            weighted: true,
            ..Default::default()
        };
        assert_eq!(
            stress_layout(&g, &strength).unwrap(),
            stress_layout(&gr, &distance).unwrap()
        );

        // All weights 2.5: the unweighted layout scaled by 2.5, up to
        // rotation (symmetric drawings have no principal axis), so compare
        // distances within components.
        let uniform = Case {
            w: vec![2.5; c.src.len()],
            ..recip
        };
        // Rounding differs between the two, so the iterations can stop a
        // step apart; distances agree to well within the tolerance.
        let tight = StressOptions {
            weighted: true,
            tol: 1e-10,
            max_iter: 5_000,
            ..Default::default()
        };
        let a = stress_layout(&uniform.graph(true), &tight).unwrap();
        let unweighted = StressOptions {
            weighted: false,
            ..tight
        };
        let b = stress_layout(&uniform.graph(false), &unweighted).unwrap();
        let comp = components(&g, ComponentKind::Weak).membership;
        for v in 0..c.n {
            for u in v + 1..c.n {
                if comp[u] == comp[v] {
                    let (da, db) = (dist(&a, u, v), dist(&b, u, v));
                    assert!((da - 2.5 * db).abs() < 1e-4 * (1.0 + da), "{da} vs {db}");
                }
            }
        }
    }
}

#[test]
fn stress_is_lower_than_simple_layouts() {
    // On connected graphs, stress majorization beats the circle and grid
    // drawings, each scaled to its own best size.
    let mut rng = Rng(15);
    let mut checked = 0;
    while checked < 20 {
        let c = random_case(&mut rng, 40);
        let g = c.graph(false);
        if c.n < 4 || components(&g, ComponentKind::Weak).count() != 1 {
            continue;
        }
        checked += 1;
        let best_scaled = |l: &Layout| -> f64 {
            // stress(sX) is a quadratic in s; take its minimum over s by a
            // fine search, since layout_stress has no scale argument.
            (1..=400)
                .map(|k| {
                    let s = k as f64 * 0.05;
                    let scaled = Layout {
                        x: l.x.iter().map(|v| v * s).collect(),
                        y: l.y.iter().map(|v| v * s).collect(),
                    };
                    layout_stress(&g, &scaled, false, WeightKind::Distance, 1).unwrap()
                })
                .fold(f64::INFINITY, f64::min)
        };
        let s = layout_stress(
            &g,
            &stress_layout(&g, &StressOptions::default()).unwrap(),
            false,
            WeightKind::Distance,
            1,
        )
        .unwrap();
        assert!(s <= best_scaled(&circle_layout(&g)) + 1e-9);
        assert!(s <= best_scaled(&grid_layout(&g, 0)) + 1e-9);
    }
}

#[test]
fn simple_layouts_are_well_formed() {
    let mut rng = Rng(16);
    for _ in 0..60 {
        let c = random_case(&mut rng, 40);
        let g = c.graph(false);
        assert!(finite(&circle_layout(&g), c.n));
        assert!(finite(&grid_layout(&g, 0), c.n));
        assert!(finite(&star_layout(&g, None).unwrap(), c.n));
        let tree = tree_layout(&g, &[]).unwrap();
        assert!(finite(&tree, c.n));
        // Depths are whole numbers: heights differ by whole numbers within
        // a component.
        let comp = components(&g, ComponentKind::Weak).membership;
        for v in 0..c.n {
            for u in 0..c.n {
                if comp[u] == comp[v] {
                    let d = tree.y[v] - tree.y[u];
                    assert!((d - d.round()).abs() < 1e-9);
                }
            }
        }
        let types: Vec<bool> = (0..c.n).map(|_| rng.below(2) == 0).collect();
        let two = twomode_layout(&g, &types, 12).unwrap();
        assert!(finite(&two, c.n));
        if types.iter().any(|&t| t) && types.iter().any(|&t| !t) {
            for v in 0..c.n {
                for u in 0..c.n {
                    if !types[v] && types[u] {
                        assert!(two.y[v] > two.y[u]);
                    }
                }
            }
        }
    }
}
