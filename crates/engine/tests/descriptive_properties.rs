//! Randomized property tests for reciprocity, transitivity and clustering,
//! distances and the measures built on them, reachability, and cut points
//! and bridges.

extern crate knappnet as engine;

use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::{
    average_clustering, bridges, cut_points, diameter, distances, eccentricity, local_clustering,
    mean_distance, path_summary, radius, reachability, reciprocity, transitivity, triangles,
    ClusteringOptions, PathOptions, ReciprocityMode, Undefined, Weighting,
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
    let n = 1 + rng.below(25) as usize;
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
    (a.is_nan() && b.is_nan()) || a == b || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

const CASES: u64 = 300;
const MODES: [DegreeMode; 3] = [DegreeMode::Out, DegreeMode::In, DegreeMode::All];

/// Floyd-Warshall on lengths, independent of the engine.
fn floyd(c: &Case, lengths: &[f64], mode: DegreeMode) -> Vec<Vec<f64>> {
    let n = c.n;
    let mut d = vec![vec![f64::INFINITY; n]; n];
    for (v, row) in d.iter_mut().enumerate() {
        row[v] = 0.0;
    }
    for (e, (&a, &b)) in c.src.iter().zip(&c.tgt).enumerate() {
        let (a, b) = (a as usize, b as usize);
        let both = !c.directed || mode == DegreeMode::All;
        let mut arcs = vec![];
        if both || mode == DegreeMode::Out {
            arcs.push((a, b));
        }
        if both || mode == DegreeMode::In {
            arcs.push((b, a));
        }
        for (x, y) in arcs {
            if x != y && lengths[e] < d[x][y] {
                d[x][y] = lengths[e];
            }
        }
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let alt = d[i][k] + d[k][j];
                if alt < d[i][j] {
                    d[i][j] = alt;
                }
            }
        }
    }
    d
}

#[test]
fn distances_match_floyd_warshall_with_every_algorithm() {
    let mut rng = Rng(101);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let ints: Vec<f64> = c.src.iter().map(|_| (1 + rng.below(9)) as f64).collect();
        let g = graph(&c, Some(&ints));
        for mode in MODES {
            let hops = floyd(&c, &vec![1.0; c.src.len()], mode);
            let weighted = floyd(&c, &ints, mode);
            let o = PathOptions {
                mode,
                threads: 1,
                ..Default::default()
            };
            let d = distances(&g, None, o).unwrap();
            for (i, row) in hops.iter().enumerate() {
                assert_eq!(d.row(i), &row[..]);
            }
            for algorithm in [
                PathAlgorithm::Auto,
                PathAlgorithm::Dial,
                PathAlgorithm::Dijkstra,
            ] {
                let w = PathOptions {
                    weighted: true,
                    algorithm,
                    ..o
                };
                let d = distances(&g, None, w).unwrap();
                for (i, row) in weighted.iter().enumerate() {
                    assert_eq!(d.row(i), &row[..], "{algorithm:?} {mode:?}");
                }
            }
        }
    }
}

#[test]
fn strengths_are_reciprocal_distances() {
    let mut rng = Rng(102);
    for _ in 0..CASES / 3 {
        let c = random_case(&mut rng);
        let w: Vec<f64> = c
            .src
            .iter()
            .map(|_| 0.1 + rng.below(1000) as f64 / 100.0)
            .collect();
        let inv: Vec<f64> = w.iter().map(|x| 1.0 / x).collect();
        let base = PathOptions {
            weighted: true,
            ..Default::default()
        };
        let a = distances(
            &graph(&c, Some(&w)),
            None,
            PathOptions {
                weights_as: WeightKind::Strength,
                ..base
            },
        )
        .unwrap();
        let b = distances(&graph(&c, Some(&inv)), None, base).unwrap();
        assert_eq!(a, b);
    }
}

#[test]
fn summaries_follow_from_the_distance_matrix() {
    let mut rng = Rng(103);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let w: Vec<f64> = c
            .src
            .iter()
            .map(|_| 0.5 + rng.below(40) as f64 / 8.0)
            .collect();
        let g = graph(&c, Some(&w));
        for mode in MODES {
            for weighted in [false, true] {
                let o = PathOptions {
                    mode,
                    weighted,
                    threads: 1,
                    ..Default::default()
                };
                let d = distances(&g, None, o).unwrap();
                let n = c.n;
                let mut finite = vec![];
                let mut complete = true;
                let mut ecc = vec![0.0f64; n];
                for (i, e) in ecc.iter_mut().enumerate() {
                    for j in 0..n {
                        let x = d.get(i, j);
                        if i != j && x.is_finite() {
                            finite.push(x);
                            *e = e.max(x);
                        } else if i != j {
                            complete = false;
                        }
                    }
                }
                let mean = if finite.is_empty() {
                    f64::NAN
                } else {
                    finite.iter().sum::<f64>() / finite.len() as f64
                };
                let longest = finite.iter().copied().fold(0.0, f64::max);
                assert!(same(mean_distance(&g, o, true).unwrap(), mean));
                assert!(same(diameter(&g, o, true).unwrap(), longest));
                let strict = mean_distance(&g, o, false).unwrap();
                assert!(if complete {
                    same(strict, mean)
                } else {
                    strict == f64::INFINITY
                });
                let got = eccentricity(&g, o).unwrap();
                assert!(got.iter().zip(&ecc).all(|(a, b)| same(*a, *b)));
                let r = radius(&g, o).unwrap();
                assert!(same(r, ecc.iter().copied().fold(f64::INFINITY, f64::min)));
                if !weighted {
                    let reach = reachability(&g, mode, 1).unwrap();
                    for (i, &x) in reach.iter().enumerate() {
                        let count = (0..n)
                            .filter(|&j| j != i && d.get(i, j).is_finite())
                            .count();
                        assert_eq!(x, count as f64);
                    }
                }
            }
        }
    }
}

#[test]
fn results_do_not_depend_on_threads() {
    let mut rng = Rng(104);
    for _ in 0..CASES / 3 {
        let c = random_case(&mut rng);
        let w: Vec<f64> = c
            .src
            .iter()
            .map(|_| 0.1 + rng.below(997) as f64 / 37.0)
            .collect();
        let g = graph(&c, Some(&w));
        let o1 = PathOptions {
            weighted: true,
            threads: 1,
            ..Default::default()
        };
        let o4 = PathOptions { threads: 4, ..o1 };
        assert_eq!(
            distances(&g, None, o1).unwrap(),
            distances(&g, None, o4).unwrap()
        );
        assert_eq!(
            path_summary(&g, None, o1).unwrap(),
            path_summary(&g, None, o4).unwrap()
        );
        assert_eq!(triangles(&g, 1), triangles(&g, 4));
        for weighting in [Weighting::None, Weighting::Barrat, Weighting::Onnela] {
            let a = ClusteringOptions {
                weighting,
                threads: 1,
                ..Default::default()
            };
            let b = ClusteringOptions { threads: 4, ..a };
            let (x, y) = (
                local_clustering(&g, a).unwrap(),
                local_clustering(&g, b).unwrap(),
            );
            assert!(x.iter().zip(&y).all(|(p, q)| p.to_bits() == q.to_bits()));
            assert_eq!(
                average_clustering(&g, a).unwrap().to_bits(),
                average_clustering(&g, b).unwrap().to_bits()
            );
        }
    }
}

/// The simple undirected reading as adjacency sets.
fn neighbor_sets(c: &Case) -> Vec<Vec<bool>> {
    let mut a = vec![vec![false; c.n]; c.n];
    for (&s, &t) in c.src.iter().zip(&c.tgt) {
        if s != t {
            a[s as usize][t as usize] = true;
            a[t as usize][s as usize] = true;
        }
    }
    a
}

#[test]
fn triangles_and_clustering_match_enumeration() {
    let mut rng = Rng(105);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = graph(&c, None);
        let a = neighbor_sets(&c);
        let n = c.n;
        let t = triangles(&g, 0);
        let local = local_clustering(&g, ClusteringOptions::default()).unwrap();
        let (mut closed, mut triples) = (0u64, 0u64);
        for v in 0..n {
            let nb: Vec<usize> = (0..n).filter(|&u| a[v][u]).collect();
            let mut tv = 0u64;
            for (i, &x) in nb.iter().enumerate() {
                for &y in &nb[i + 1..] {
                    tv += a[x][y] as u64;
                }
            }
            assert_eq!(t[v], tv);
            let k = nb.len() as u64;
            closed += tv;
            triples += k * k.saturating_sub(1) / 2;
            if k < 2 {
                assert!(local[v].is_nan());
            } else {
                assert!(same(local[v], tv as f64 / (k * (k - 1) / 2) as f64));
            }
        }
        assert_eq!(closed % 3, 0);
        let global = transitivity(&g, ClusteringOptions::default());
        if triples == 0 {
            assert!(global.is_nan());
        } else {
            assert!(same(global, closed as f64 / triples as f64));
        }
    }
}

#[test]
fn equal_weights_reduce_weighted_clustering_to_binary() {
    let mut rng = Rng(106);
    for _ in 0..CASES / 3 {
        let mut c = random_case(&mut rng);
        // Barrat's and Onnela's coefficients reduce to the binary one when
        // every pair of adjacent nodes has the same total weight, so keep
        // one undirected edge per pair.
        c.directed = false;
        let mut seen = std::collections::HashSet::new();
        let (mut s, mut t) = (vec![], vec![]);
        for (&a, &b) in c.src.iter().zip(&c.tgt) {
            if a != b && seen.insert((a.min(b), a.max(b))) {
                s.push(a);
                t.push(b);
            }
        }
        c.src = s;
        c.tgt = t;
        let w = vec![2.5; c.src.len()];
        let g = graph(&c, Some(&w));
        let binary = local_clustering(&g, ClusteringOptions::default()).unwrap();
        for weighting in [Weighting::Barrat, Weighting::Onnela] {
            let o = ClusteringOptions {
                weighting,
                ..Default::default()
            };
            let x = local_clustering(&g, o).unwrap();
            assert!(x.iter().zip(&binary).all(|(a, b)| same(*a, *b)));
        }
    }
}

fn components_without(c: &Case, node: Option<usize>, edge: Option<usize>) -> usize {
    let mut parent: Vec<usize> = (0..c.n).collect();
    fn find(p: &mut [usize], mut v: usize) -> usize {
        while p[v] != v {
            v = p[v];
        }
        v
    }
    for e in 0..c.src.len() {
        let (a, b) = (c.src[e] as usize, c.tgt[e] as usize);
        if Some(e) == edge || node == Some(a) || node == Some(b) {
            continue;
        }
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..c.n)
        .filter(|&v| Some(v) != node && find(&mut parent, v) == v)
        .count()
}

#[test]
fn cut_points_and_bridges_match_deletion() {
    let mut rng = Rng(107);
    for _ in 0..CASES {
        let mut c = random_case(&mut rng);
        // Sparser networks have more cut points and bridges.
        let keep = c.src.len() / 3 + 1;
        c.src.truncate(keep);
        c.tgt.truncate(keep);
        let g = graph(&c, None);
        let base = components_without(&c, None, None);
        let isolated = |v: usize| {
            c.src
                .iter()
                .zip(&c.tgt)
                .all(|(&a, &b)| a == b || (a as usize != v && b as usize != v))
        };
        let expect: Vec<u32> = (0..c.n)
            .filter(|&v| {
                let after = components_without(&c, Some(v), None);
                after > base - usize::from(isolated(v))
            })
            .map(|v| v as u32)
            .collect();
        assert_eq!(cut_points(&g), expect);
        let expect: Vec<usize> = (0..c.src.len())
            .filter(|&e| c.src[e] != c.tgt[e] && components_without(&c, None, Some(e)) > base)
            .collect();
        assert_eq!(bridges(&g), expect);
    }
}

#[test]
fn relabeling_permutes_results() {
    let mut rng = Rng(108);
    for _ in 0..CASES / 3 {
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
        for mode in [
            ReciprocityMode::Default,
            ReciprocityMode::Ratio,
            ReciprocityMode::Dyadic,
            ReciprocityMode::Correlation,
        ] {
            assert!(same(reciprocity(&g, mode), reciprocity(&h, mode)));
        }
        let o = ClusteringOptions {
            undefined: Undefined::Zero,
            ..Default::default()
        };
        assert!(same(transitivity(&g, o), transitivity(&h, o)));
        let (a, b) = (
            local_clustering(&g, o).unwrap(),
            local_clustering(&h, o).unwrap(),
        );
        assert!((0..c.n).all(|v| same(a[v], b[perm[v]])));
        let mut cg: Vec<u32> = cut_points(&g)
            .iter()
            .map(|&v| perm[v as usize] as u32)
            .collect();
        cg.sort_unstable();
        assert_eq!(cg, cut_points(&h));
        assert_eq!(bridges(&g), bridges(&h));
        let (e1, e2) = (
            eccentricity(&g, PathOptions::default()).unwrap(),
            eccentricity(&h, PathOptions::default()).unwrap(),
        );
        assert!((0..c.n).all(|v| e1[v] == e2[perm[v]]));
    }
}

#[test]
fn reciprocity_bounds() {
    let mut rng = Rng(109);
    for _ in 0..CASES {
        let c = random_case(&mut rng);
        let g = graph(&c, None);
        for mode in [
            ReciprocityMode::Default,
            ReciprocityMode::Ratio,
            ReciprocityMode::Dyadic,
        ] {
            let r = reciprocity(&g, mode);
            assert!(r.is_nan() || (0.0..=1.0).contains(&r));
        }
        let rho = reciprocity(&g, ReciprocityMode::Correlation);
        assert!(rho.is_nan() || (-1.0 - 1e-12..=1.0 + 1e-12).contains(&rho));
        // The ratio and default indices are tied by r_ratio = r / (2 - r).
        let r = reciprocity(&g, ReciprocityMode::Default);
        assert!(same(reciprocity(&g, ReciprocityMode::Ratio), r / (2.0 - r)));
    }
}
