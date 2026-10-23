//! Randomized property tests for cohesive subsets (bicomponents,
//! n-cliques, n-clans, k-plexes, lambda sets, factions), community
//! partitions (Girvan-Newman, leading eigenvector, greedy modularity),
//! topography (core-periphery fits, compactness, small worlds), and the
//! random graph generators.
//!
//! The checks here are independent of the engine's algorithms: sets are
//! tested against their definitions with plain breadth-first searches and
//! exhaustive enumeration on small networks, and fits against direct
//! recomputation of the criteria.

extern crate knappnet as engine;

use std::collections::{HashSet, VecDeque};

use engine::generate::{self, Degrees, EdgeList};
use engine::metrics::cliques::max_cliques;
use engine::metrics::community::modularity;
use engine::metrics::coreperiphery::{
    compactness, continuous_core_periphery, core_periphery, small_world, ContinuousOptions,
    CoreLinks, CorePeripheryOptions, SmallWorldOptions,
};
use engine::metrics::degree::DegreeMode;
use engine::metrics::descriptive::PathOptions;
use engine::metrics::factions::{
    bicomponents, factions, k_plexes, lambda_connectivity, n_clans, n_cliques, FactionFit,
    FactionOptions, KPlexOptions,
};
use engine::metrics::partition::{fast_greedy, girvan_newman, leading_eigen};
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
}

struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    directed: bool,
}

fn random_case(rng: &mut Rng, max_n: usize) -> Case {
    let n = 1 + rng.below(max_n as u64) as usize;
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

fn graph(c: &Case) -> Graph {
    Graph::from_edges(c.n, &c.src, &c.tgt, None, c.directed).unwrap()
}

/// Undirected simple adjacency: an arc in either direction joins two
/// distinct nodes.
fn simple_adj(c: &Case) -> Vec<HashSet<usize>> {
    let mut adj = vec![HashSet::new(); c.n];
    for (&s, &t) in c.src.iter().zip(&c.tgt) {
        if s != t {
            adj[s as usize].insert(t as usize);
            adj[t as usize].insert(s as usize);
        }
    }
    adj
}

/// Distances from `s` within `allowed` (every node when None); usize::MAX
/// means unreachable.
fn bfs(adj: &[HashSet<usize>], s: usize, allowed: Option<&HashSet<usize>>) -> Vec<usize> {
    let mut d = vec![usize::MAX; adj.len()];
    d[s] = 0;
    let mut q = VecDeque::from([s]);
    while let Some(v) = q.pop_front() {
        for &u in &adj[v] {
            if allowed.map_or(true, |a| a.contains(&u)) && d[u] == usize::MAX {
                d[u] = d[v] + 1;
                q.push_back(u);
            }
        }
    }
    d
}

/// Equality that treats NaN as equal to NaN (via Debug output).
fn same_debug<T: std::fmt::Debug>(a: &T, b: &T) -> bool {
    format!("{a:?}") == format!("{b:?}")
}

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn sets_to_usize(sets: &[Vec<u32>]) -> Vec<Vec<usize>> {
    sets.iter()
        .map(|s| s.iter().map(|&v| v as usize).collect())
        .collect()
}

fn is_sorted_canonical(sets: &[Vec<u32>]) -> bool {
    sets.iter().all(|s| s.windows(2).all(|w| w[0] < w[1])) && sets.windows(2).all(|w| w[0] < w[1])
}

const CASES: u64 = 200;

// ---- Bicomponents ---------------------------------------------------------

#[test]
fn bicomponents_cover_every_edge_once_and_stay_connected_without_any_node() {
    let mut rng = Rng(11);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        let g = graph(&c);
        let adj = simple_adj(&c);
        let blocks = bicomponents(&g);
        assert!(is_sorted_canonical(&blocks));
        let blocks = sets_to_usize(&blocks);
        let sets: Vec<HashSet<usize>> =
            blocks.iter().map(|b| b.iter().copied().collect()).collect();
        // Every edge lies in exactly one block.
        for (v, nbrs) in adj.iter().enumerate() {
            for &u in nbrs {
                if u > v {
                    let count = sets
                        .iter()
                        .filter(|b| b.contains(&u) && b.contains(&v))
                        .count();
                    assert_eq!(count, 1, "edge {v}-{u} in {count} blocks");
                }
            }
        }
        // Blocks have at least two nodes and share at most one.
        for (i, a) in sets.iter().enumerate() {
            assert!(a.len() >= 2);
            for b in &sets[i + 1..] {
                assert!(a.intersection(b).count() <= 1);
            }
            // Removing any one node leaves the block connected.
            for &x in a {
                let rest: HashSet<usize> = a.iter().copied().filter(|&y| y != x).collect();
                let start = *rest.iter().next().unwrap();
                let d = bfs(&adj, start, Some(&rest));
                assert!(
                    rest.iter().all(|&y| d[y] != usize::MAX),
                    "block {a:?} cut by {x}"
                );
            }
        }
        // Non-isolated nodes are covered; isolated nodes are not.
        for (v, nbrs) in adj.iter().enumerate() {
            let covered = sets.iter().any(|b| b.contains(&v));
            assert_eq!(covered, !nbrs.is_empty());
        }
    }
}

// ---- n-cliques, n-clans, k-plexes -------------------------------------------

#[test]
fn n_cliques_are_maximal_sets_within_distance_n_and_clans_have_small_diameter() {
    let mut rng = Rng(23);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 14);
        let g = graph(&c);
        let adj = simple_adj(&c);
        let dist: Vec<Vec<usize>> = (0..c.n).map(|v| bfs(&adj, v, None)).collect();
        // 1-cliques are the maximal cliques.
        assert_eq!(n_cliques(&g, 1, 0, 2).unwrap(), max_cliques(&g, 0, 1));
        for order in 1..=3 {
            let cl = n_cliques(&g, order, 0, 2).unwrap();
            assert!(is_sorted_canonical(&cl));
            let cl = sets_to_usize(&cl);
            let within = |a: usize, b: usize| dist[a][b] <= order;
            for s in &cl {
                assert!(s.iter().all(|&a| s.iter().all(|&b| within(a, b))), "{s:?}");
                // Maximal: no outside node is within n of every member.
                for v in 0..c.n {
                    if !s.contains(&v) {
                        assert!(!s.iter().all(|&a| within(a, v)), "{s:?} + {v}");
                    }
                }
            }
            // Complete: every maximal such set is listed (exhaustive for
            // small n, by growing sets greedily from each node subset).
            let listed: HashSet<Vec<usize>> = cl.iter().cloned().collect();
            if c.n <= 10 {
                for bits in 1u32..(1 << c.n) {
                    let s: Vec<usize> = (0..c.n).filter(|&v| bits >> v & 1 == 1).collect();
                    let ok = s.iter().all(|&a| s.iter().all(|&b| within(a, b)));
                    let maximal =
                        ok && (0..c.n).all(|v| s.contains(&v) || !s.iter().all(|&a| within(a, v)));
                    if maximal {
                        assert!(listed.contains(&s), "missing {order}-clique {s:?}");
                    }
                }
            }
            // n-clans: the n-cliques with inner diameter at most n.
            let clans = sets_to_usize(&n_clans(&g, order, 0, 2).unwrap());
            let expected: Vec<Vec<usize>> = cl
                .iter()
                .filter(|s| {
                    let allowed: HashSet<usize> = s.iter().copied().collect();
                    s.iter().all(|&a| {
                        let d = bfs(&adj, a, Some(&allowed));
                        s.iter().all(|&b| d[b] <= order)
                    })
                })
                .cloned()
                .collect();
            assert_eq!(clans, expected);
            // min_size only filters.
            let big = sets_to_usize(&n_cliques(&g, order, 3, 2).unwrap());
            let filtered: Vec<Vec<usize>> = cl.iter().filter(|s| s.len() >= 3).cloned().collect();
            assert_eq!(big, filtered);
        }
    }
}

#[test]
fn k_plexes_are_maximal_and_complete() {
    let mut rng = Rng(37);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 11);
        let g = graph(&c);
        let adj = simple_adj(&c);
        let is_plex = |s: &[usize], k: usize| {
            s.iter().all(|&v| {
                s.iter()
                    .filter(|&&u| u != v && !adj[v].contains(&u))
                    .count()
                    < k
            })
        };
        for k in 1..=3 {
            for min_size in [1usize, 2 * k - 1, 2 * k] {
                let opts = KPlexOptions {
                    k,
                    min_size,
                    max_steps: 1_000_000,
                    threads: 2,
                };
                let found = k_plexes(&g, opts).unwrap();
                assert!(is_sorted_canonical(&found));
                let found = sets_to_usize(&found);
                let listed: HashSet<Vec<usize>> = found.iter().cloned().collect();
                for s in &found {
                    assert!(s.len() >= min_size);
                    assert!(is_plex(s, k), "{s:?} is not a {k}-plex");
                }
                // Every maximal k-plex of at least min_size nodes is listed,
                // and nothing else.
                let mut expected = HashSet::new();
                for bits in 1u32..(1 << c.n) {
                    let s: Vec<usize> = (0..c.n).filter(|&v| bits >> v & 1 == 1).collect();
                    if s.len() < min_size || !is_plex(&s, k) {
                        continue;
                    }
                    let maximal = (0..c.n).all(|v| {
                        if s.contains(&v) {
                            return true;
                        }
                        let mut t = s.clone();
                        t.push(v);
                        t.sort_unstable();
                        !is_plex(&t, k)
                    });
                    if maximal {
                        expected.insert(s);
                    }
                }
                assert_eq!(listed, expected, "k = {k}, min_size = {min_size}");
                if k == 1 {
                    let cl: HashSet<Vec<usize>> = sets_to_usize(&max_cliques(&g, min_size, 1))
                        .into_iter()
                        .collect();
                    assert_eq!(listed, cl);
                }
            }
        }
    }
}

// ---- Lambda sets ----------------------------------------------------------

/// Edge-disjoint paths between s and t by augmenting paths on unit
/// capacities (Menger), independent of the engine's Gusfield tree.
fn edge_disjoint_paths(adj: &[HashSet<usize>], s: usize, t: usize) -> u64 {
    let n = adj.len();
    let mut cap = vec![vec![0u8; n]; n];
    for v in 0..n {
        for &u in &adj[v] {
            cap[v][u] = 1;
        }
    }
    let mut flow = 0;
    loop {
        let mut prev = vec![usize::MAX; n];
        prev[s] = s;
        let mut q = VecDeque::from([s]);
        while let Some(v) = q.pop_front() {
            for u in 0..n {
                if cap[v][u] > 0 && prev[u] == usize::MAX {
                    prev[u] = v;
                    q.push_back(u);
                }
            }
        }
        if prev[t] == usize::MAX {
            return flow;
        }
        let mut v = t;
        while v != s {
            let p = prev[v];
            cap[p][v] -= 1;
            cap[v][p] += 1;
            v = p;
        }
        flow += 1;
    }
}

#[test]
fn lambda_connectivity_matches_edge_disjoint_paths() {
    let mut rng = Rng(41);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 12);
        let g = graph(&c);
        let adj = simple_adj(&c);
        let m = lambda_connectivity(&g);
        assert_eq!((m.rows(), m.cols()), (c.n, c.n));
        for a in 0..c.n {
            assert_eq!(m.get(a, a), 0.0);
            for b in 0..c.n {
                assert_eq!(m.get(a, b), m.get(b, a));
                if a != b {
                    assert_eq!(
                        m.get(a, b),
                        edge_disjoint_paths(&adj, a, b) as f64,
                        "{a} {b}"
                    );
                    assert!(m.get(a, b) <= adj[a].len().min(adj[b].len()) as f64);
                }
                for x in 0..c.n {
                    // lambda(a, b) >= min(lambda(a, x), lambda(x, b)).
                    if a != b && x != a && x != b {
                        assert!(m.get(a, b) >= m.get(a, x).min(m.get(x, b)));
                    }
                }
            }
        }
    }
}

// ---- Factions -------------------------------------------------------------

/// Binary adjacency as the fits read it: arc a -> b (both ways when
/// undirected), loops dropped.
fn arcs(c: &Case) -> Vec<Vec<bool>> {
    let mut arc = vec![vec![false; c.n]; c.n];
    for (&s, &t) in c.src.iter().zip(&c.tgt) {
        if s != t {
            arc[s as usize][t as usize] = true;
            if !c.directed {
                arc[t as usize][s as usize] = true;
            }
        }
    }
    arc
}

/// Errors and correlation of a partition, from the definitions.
fn faction_fit(c: &Case, arc: &[Vec<bool>], labels: &[usize]) -> (u64, f64) {
    let (mut pairs, mut ties, mut within, mut tied_within) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for a in 0..c.n {
        for b in 0..c.n {
            if a == b || (!c.directed && b < a) {
                continue;
            }
            pairs += 1.0;
            let tie = arc[a][b];
            let same = labels[a] == labels[b];
            if tie {
                ties += 1.0;
            }
            if same {
                within += 1.0;
                if tie {
                    tied_within += 1.0;
                }
            }
        }
    }
    let errors = (within - tied_within) + (ties - tied_within);
    let num = pairs * tied_within - ties * within;
    let den = ((pairs * ties - ties * ties) * (pairs * within - within * within)).sqrt();
    (errors as u64, num / den)
}

#[test]
fn factions_report_their_own_fit_and_find_the_optimum_on_small_networks() {
    let mut rng = Rng(53);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 9);
        let g = graph(&c);
        let arc = arcs(&c);
        let pairs = if c.directed {
            c.n * (c.n - 1)
        } else {
            c.n * (c.n - 1) / 2
        };
        let ties = (0..c.n)
            .map(|a| {
                (0..c.n)
                    .filter(|&b| arc[a][b] && (c.directed || b > a))
                    .count()
            })
            .sum::<usize>();
        for k in 1..=3.min(c.n) {
            for fit in [FactionFit::Errors, FactionFit::Correlation] {
                let opts = FactionOptions {
                    k,
                    fit,
                    starts: 10,
                    max_stall: 50,
                    tenure: 0,
                    seed: 7,
                    threads: 2,
                };
                let r = factions(&g, opts);
                if fit == FactionFit::Correlation
                    && (k == 1 || k == c.n || ties == 0 || ties == pairs)
                {
                    assert!(r.is_err(), "{r:?}");
                    continue;
                }
                let r = r.unwrap();
                // Exactly k groups, numbered by first appearance.
                let mut seen = Vec::new();
                for &l in &r.membership {
                    if !seen.contains(&l) {
                        assert_eq!(l, seen.len());
                        seen.push(l);
                    }
                }
                assert_eq!(seen.len(), k);
                let (errors, corr) = faction_fit(&c, &arc, &r.membership);
                assert_eq!(r.errors, errors);
                assert!(same(r.correlation, corr), "{} vs {corr}", r.correlation);
                // Thread count does not matter.
                let again = factions(&g, FactionOptions { threads: 1, ..opts }).unwrap();
                assert!(same_debug(&r, &again));
                // Optimum over every partition into k non-empty groups.
                let mut best_err = u64::MAX;
                let mut best_corr = f64::NEG_INFINITY;
                let mut labels = vec![0usize; c.n];
                let total = (k as u64).pow(c.n as u32);
                for code in 0..total {
                    let mut x = code;
                    for l in labels.iter_mut() {
                        *l = (x % k as u64) as usize;
                        x /= k as u64;
                    }
                    if (0..k).any(|grp| !labels.contains(&grp)) {
                        continue;
                    }
                    let (e, cr) = faction_fit(&c, &arc, &labels);
                    best_err = best_err.min(e);
                    if !cr.is_nan() {
                        best_corr = best_corr.max(cr);
                    }
                }
                match fit {
                    FactionFit::Errors => assert_eq!(r.errors, best_err),
                    FactionFit::Correlation => {
                        if best_corr.is_finite() {
                            assert!(same(r.correlation, best_corr) || r.correlation > best_corr);
                        }
                    }
                }
            }
        }
    }
}

// ---- Partitions -----------------------------------------------------------

fn undirected_case(rng: &mut Rng, max_n: usize) -> Case {
    let mut c = random_case(rng, max_n);
    c.directed = false;
    c
}

fn canonical_labels(labels: &[usize]) -> bool {
    let mut next = 0;
    for &l in labels {
        if l == next {
            next += 1;
        } else if l > next {
            return false;
        }
    }
    true
}

fn as_i64(labels: &[usize]) -> Vec<i64> {
    labels.iter().map(|&l| l as i64).collect()
}

#[test]
fn partitions_report_the_modularity_of_their_membership() {
    let mut rng = Rng(67);
    for _ in 0..CASES {
        let c = undirected_case(&mut rng, 20);
        let g = graph(&c);
        let simple = {
            let adj = simple_adj(&c);
            let mut s = Vec::new();
            let mut t = Vec::new();
            for (v, nbrs) in adj.iter().enumerate() {
                for &u in nbrs {
                    if u > v {
                        s.push(v as i64);
                        t.push(u as i64);
                    }
                }
            }
            Graph::from_edges(c.n, &s, &t, None, false).unwrap()
        };
        let gn = girvan_newman(&g, 2).unwrap();
        assert!(canonical_labels(&gn.membership));
        assert!(same(
            gn.modularity,
            modularity(&simple, &as_i64(&gn.membership), 1.0, false).unwrap()
        ));
        assert!(same_debug(&gn, &girvan_newman(&g, 1).unwrap()));
        // The path records the modularity after each removal, and the
        // chosen partition is its maximum.
        assert_eq!(gn.path.len(), gn.removed.len() + 1);
        if simple.edge_count() > 0 {
            let best = gn.path.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            assert!(same(gn.modularity, best));
            assert_eq!(gn.removed.len(), simple.edge_count());
        }
        let le = leading_eigen(&g, false).unwrap();
        assert!(canonical_labels(&le.membership));
        assert!(same(
            le.modularity,
            modularity(&g, &as_i64(&le.membership), 1.0, false).unwrap()
        ));
        for resolution in [0.5, 1.0, 2.0] {
            let fg = fast_greedy(&g, resolution, false).unwrap();
            assert!(canonical_labels(&fg.membership));
            assert!(same(
                fg.modularity,
                modularity(&g, &as_i64(&fg.membership), resolution, false).unwrap()
            ));
            // Groups never span components (merges only join connected
            // groups), so every group is connected in the simple graph.
            let adj = simple_adj(&c);
            let k = fg.membership.iter().max().map_or(0, |&m| m + 1);
            for grp in 0..k {
                let members: HashSet<usize> =
                    (0..c.n).filter(|&v| fg.membership[v] == grp).collect();
                let start = *members.iter().next().unwrap();
                let d = bfs(&adj, start, Some(&members));
                assert!(members.iter().all(|&v| d[v] != usize::MAX));
            }
        }
        // Directed networks are refused.
        let d = Graph::from_edges(c.n, &c.src, &c.tgt, None, true).unwrap();
        assert!(girvan_newman(&d, 1).is_err());
        assert!(leading_eigen(&d, false).is_err());
        assert!(fast_greedy(&d, 1.0, false).is_err());
    }
}

// ---- Core-periphery, compactness, small worlds -------------------------------

/// Pearson correlation.
fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    let syy: f64 = y.iter().map(|b| (b - my) * (b - my)).sum();
    sxy / (sxx * syy).sqrt()
}

#[test]
fn core_periphery_fits_report_their_correlation() {
    let mut rng = Rng(71);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 12);
        let g = graph(&c);
        let adj = simple_adj(&c);
        let pairs: Vec<(usize, usize)> = (0..c.n)
            .flat_map(|a| (0..c.n).map(move |b| (a, b)))
            .filter(|&(a, b)| a != b && (c.directed || a < b))
            .collect();
        let arc = arcs(&c);
        let x: Vec<f64> = pairs.iter().map(|&(a, b)| arc[a][b] as u8 as f64).collect();
        let ties = x.iter().sum::<f64>();
        let degenerate = c.n < 2 || ties == 0.0 || ties == x.len() as f64;
        for links in [CoreLinks::Missing, CoreLinks::Ones, CoreLinks::Zeros] {
            let opts = CorePeripheryOptions {
                links,
                weighted: false,
                starts: 5,
                max_stall: 50,
                tenure: 0,
                seed: 3,
                threads: 2,
            };
            let r = core_periphery(&g, opts);
            if degenerate {
                assert!(r.is_err());
                continue;
            }
            let r = r.unwrap();
            let again = core_periphery(&g, CorePeripheryOptions { threads: 1, ..opts }).unwrap();
            assert!(same_debug(&r, &again));
            // Recompute the correlation of the reported core.
            let (mut xs, mut ys) = (Vec::new(), Vec::new());
            for (i, &(a, b)) in pairs.iter().enumerate() {
                let (ca, cb) = (r.core[a], r.core[b]);
                let d = match (ca && cb, !ca && !cb) {
                    (true, _) => 1.0,
                    (_, true) => 0.0,
                    _ => match links {
                        CoreLinks::Missing => continue,
                        CoreLinks::Ones => 1.0,
                        CoreLinks::Zeros => 0.0,
                    },
                };
                xs.push(x[i]);
                ys.push(d);
            }
            let expected = pearson(&xs, &ys);
            assert!(
                same(r.correlation, expected),
                "{} vs {expected}",
                r.correlation
            );
        }
        let opts = ContinuousOptions {
            weighted: false,
            max_iter: 500,
            starts: 3,
            seed: 3,
            threads: 2,
        };
        let r = continuous_core_periphery(&g, opts);
        if degenerate {
            assert!(r.is_err());
            continue;
        }
        // Two nodes give a single pair, whose pattern has no variance.
        let Ok(r) = r else {
            assert!(c.n == 2, "{r:?}");
            continue;
        };
        let again =
            continuous_core_periphery(&g, ContinuousOptions { threads: 1, ..opts }).unwrap();
        assert!(same_debug(&r, &again));
        assert!(r.coreness.iter().all(|&v| v >= 0.0));
        assert!((r.coreness.iter().map(|v| v * v).sum::<f64>() - 1.0).abs() < 1e-9);
        let y: Vec<f64> = pairs
            .iter()
            .map(|&(a, b)| r.coreness[a] * r.coreness[b])
            .collect();
        let expected = pearson(&x, &y);
        assert!(same(r.correlation, expected) || (r.correlation.is_nan() && expected.is_nan()));
        // At least as good as the best single-node and degree patterns,
        // which are feasible points.
        if !r.correlation.is_nan() {
            let deg: Vec<f64> = (0..c.n).map(|v| adj[v].len() as f64 + 1.0).collect();
            let yd: Vec<f64> = pairs.iter().map(|&(a, b)| deg[a] * deg[b]).collect();
            let rd = pearson(&x, &yd);
            if !rd.is_nan() {
                assert!(
                    r.correlation + 1e-9 >= rd,
                    "{} < degree start {rd}",
                    r.correlation
                );
            }
        }
    }
}

#[test]
fn compactness_is_the_mean_reciprocal_distance() {
    let mut rng = Rng(83);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 20);
        let g = graph(&c);
        let opts = PathOptions {
            mode: DegreeMode::Out,
            weighted: false,
            threads: 2,
            ..Default::default()
        };
        let v = compactness(&g, opts).unwrap();
        if c.n < 2 {
            assert!(v.is_nan());
            continue;
        }
        // Directed distances by breadth-first search.
        let mut out = vec![HashSet::new(); c.n];
        for (&s, &t) in c.src.iter().zip(&c.tgt) {
            if s != t {
                out[s as usize].insert(t as usize);
                if !c.directed {
                    out[t as usize].insert(s as usize);
                }
            }
        }
        let mut total = 0.0;
        for s in 0..c.n {
            for (t, &d) in bfs(&out, s, None).iter().enumerate() {
                if t != s && d != usize::MAX {
                    total += 1.0 / d as f64;
                }
            }
        }
        let expected = total / (c.n * (c.n - 1)) as f64;
        assert!(same(v, expected), "{v} vs {expected}");
        assert!((0.0..=1.0 + 1e-12).contains(&v));
        // Ignoring direction can only shorten paths.
        let all = compactness(
            &g,
            PathOptions {
                mode: DegreeMode::All,
                ..opts
            },
        )
        .unwrap();
        assert!(all + 1e-12 >= v);
    }
}

#[test]
fn small_world_quotient_is_the_ratio_of_ratios_and_is_reproducible() {
    let mut rng = Rng(97);
    for _ in 0..40 {
        let c = random_case(&mut rng, 25);
        let g = graph(&c);
        for local in [false, true] {
            let opts = SmallWorldOptions {
                local,
                reps: 5,
                seed: 11,
                threads: 2,
            };
            let r = small_world(&g, opts).unwrap();
            let again = small_world(&g, SmallWorldOptions { threads: 1, ..opts }).unwrap();
            assert!(same_debug(&r, &again));
            let q =
                (r.clustering / r.random_clustering) / (r.mean_distance / r.random_mean_distance);
            assert!(same(r.quotient, q) || (r.quotient.is_nan() && q.is_nan()));
            // Another seed gives other random networks but the same
            // observed values.
            let other = small_world(&g, SmallWorldOptions { seed: 12, ..opts }).unwrap();
            assert!(same(other.clustering, r.clustering));
            assert!(same(other.mean_distance, r.mean_distance));
        }
    }
}

// ---- Random graphs ----------------------------------------------------------

fn check_edges(e: &EdgeList, n: usize, simple: bool) {
    assert_eq!(e.n, n);
    assert_eq!(e.from.len(), e.to.len());
    let pairs: Vec<(u32, u32)> = e.from.iter().copied().zip(e.to.iter().copied()).collect();
    assert!(pairs
        .iter()
        .all(|&(a, b)| (a as usize) < n && (b as usize) < n));
    if simple {
        assert!(
            pairs.windows(2).all(|w| w[0] < w[1]),
            "unsorted or repeated edge"
        );
        assert!(pairs.iter().all(|&(a, b)| a != b), "self-loop");
    } else {
        assert!(pairs.windows(2).all(|w| w[0] <= w[1]));
    }
    if !e.directed {
        assert!(pairs.iter().all(|&(a, b)| a <= b));
    }
}

fn degrees(e: &EdgeList) -> (Vec<u64>, Vec<u64>) {
    let (mut dout, mut din) = (vec![0u64; e.n], vec![0u64; e.n]);
    for (&a, &b) in e.from.iter().zip(&e.to) {
        dout[a as usize] += 1;
        din[b as usize] += 1;
    }
    (dout, din)
}

#[test]
fn generators_keep_their_exact_constraints_and_reproduce_by_seed() {
    let mut rng = Rng(101);
    for _ in 0..60 {
        let n = 2 + rng.below(30) as usize;
        let directed = rng.below(2) == 0;
        let loops = rng.below(2) == 0;
        let seed = rng.next();
        let possible = if directed {
            n * (n - 1)
        } else {
            n * (n - 1) / 2
        } + if loops { n } else { 0 };
        let m = rng.below(possible as u64 + 1);
        // G(n, m): exactly m distinct pairs.
        let lists = generate::gnm(n, m, directed, loops, 3, seed, 2).unwrap();
        assert_eq!(lists.len(), 3);
        for e in &lists {
            check_edges(e, n, !loops);
            assert_eq!(e.len() as u64, m);
            assert_eq!(e.directed, directed);
            if !loops {
                assert!(e.from.iter().zip(&e.to).all(|(a, b)| a != b));
            }
        }
        assert_eq!(
            lists,
            generate::gnm(n, m, directed, loops, 3, seed, 1).unwrap()
        );
        assert_eq!(
            lists[1],
            generate::gnm(n, m, directed, loops, 2, seed, 1).unwrap()[1]
        );
        // G(n, p): simple unless loops; p = 1 gives every pair.
        let p = rng.below(101) as f64 / 100.0;
        let lists = generate::gnp(n, p, directed, loops, 2, seed, 2).unwrap();
        for e in &lists {
            check_edges(e, n, !loops);
            assert!(e.len() <= possible);
        }
        assert_eq!(
            lists,
            generate::gnp(n, p, directed, loops, 2, seed, 1).unwrap()
        );
        let full = generate::gnp(n, 1.0, directed, loops, 1, seed, 1).unwrap();
        assert_eq!(full[0].len(), possible);
        assert!(generate::gnp(n, 0.0, directed, loops, 1, seed, 1).unwrap()[0].is_empty());

        // Configuration model: degrees exact, multigraph and simple.
        let base = &generate::gnm(n, m.min((n * (n - 1) / 2) as u64), false, false, 1, seed, 1)
            .unwrap()[0];
        let (d, _) = degrees(base);
        let (_, din) = degrees(base);
        let und: Vec<u64> = d.iter().zip(&din).map(|(a, b)| a + b).collect();
        let lists =
            generate::configuration(&Degrees::Undirected(und.clone()), false, 1, 3, seed, 2)
                .unwrap();
        for e in &lists {
            check_edges(e, n, false);
            let (a, b) = degrees(e);
            let total: Vec<u64> = a.iter().zip(&b).map(|(x, y)| x + y).collect();
            assert_eq!(total, und);
        }
        assert_eq!(
            lists,
            generate::configuration(&Degrees::Undirected(und.clone()), false, 1, 3, seed, 1)
                .unwrap()
        );
        if let Ok(simple) =
            generate::configuration(&Degrees::Undirected(und.clone()), true, 200, 2, seed, 2)
        {
            for e in &simple {
                check_edges(e, n, true);
                let (a, b) = degrees(e);
                let total: Vec<u64> = a.iter().zip(&b).map(|(x, y)| x + y).collect();
                assert_eq!(total, und);
            }
        }
        let dir =
            &generate::gnm(n, m.min((n * (n - 1)) as u64), true, false, 1, seed, 1).unwrap()[0];
        let (dout, din) = degrees(dir);
        let lists = generate::configuration(
            &Degrees::Directed(dout.clone(), din.clone()),
            false,
            1,
            2,
            seed,
            2,
        )
        .unwrap();
        for e in &lists {
            assert!(e.directed);
            let (a, b) = degrees(e);
            assert_eq!((a, b), (dout.clone(), din.clone()));
        }

        // U|MAN: the dyad census is exact.
        let dyads = (n * (n - 1) / 2) as u64;
        let mutual = rng.below(dyads + 1);
        let asym = rng.below(dyads - mutual + 1);
        let lists = generate::uman(n, mutual, asym, 2, seed, 2).unwrap();
        for e in &lists {
            check_edges(e, n, true);
            assert!(e.directed);
            let set: HashSet<(u32, u32)> =
                e.from.iter().copied().zip(e.to.iter().copied()).collect();
            let (mut mm, mut aa) = (0, 0);
            for a in 0..n as u32 {
                for b in a + 1..n as u32 {
                    match (set.contains(&(a, b)), set.contains(&(b, a))) {
                        (true, true) => mm += 1,
                        (false, false) => {}
                        _ => aa += 1,
                    }
                }
            }
            assert_eq!((mm, aa), (mutual, asym));
        }
        assert_eq!(lists, generate::uman(n, mutual, asym, 2, seed, 1).unwrap());

        // Watts-Strogatz: n k / 2 edges, simple; p = 0 is the lattice.
        let k = 2 * (1 + rng.below(((n - 1) / 2).max(1) as u64) as usize);
        if k < n {
            let lists = generate::small_world(n, k, p, 2, seed, 2).unwrap();
            for e in &lists {
                check_edges(e, n, true);
                assert_eq!(e.len(), n * k / 2);
                assert!(!e.directed);
            }
            assert_eq!(lists, generate::small_world(n, k, p, 2, seed, 1).unwrap());
            let lattice = &generate::small_world(n, k, 0.0, 1, seed, 1).unwrap()[0];
            let (a, b) = degrees(lattice);
            assert!(a.iter().zip(&b).all(|(x, y)| (x + y) as usize == k));
        }
        // Preferential attachment: m (n - m) edges, simple.
        let ma = 1 + rng.below((n - 1) as u64) as usize;
        let lists = generate::preferential(n, ma, 2, seed, 2).unwrap();
        for e in &lists {
            check_edges(e, n, true);
            assert_eq!(e.len(), ma * (n - ma));
            let (a, b) = degrees(e);
            assert!((0..n).all(|v| a[v] + b[v] >= 1));
        }
        assert_eq!(lists, generate::preferential(n, ma, 2, seed, 1).unwrap());
        // Stochastic block model: p = 1 within blocks and 0 between gives
        // disjoint cliques.
        let sizes = [n / 2, n - n / 2];
        let probs = [1.0, 0.0, 0.0, 1.0];
        let e = &generate::sbm(&sizes, &probs, directed, loops, 1, seed, 1).unwrap()[0];
        check_edges(e, n, !loops);
        let within = |v: u32| (v as usize) < sizes[0];
        assert!(e
            .from
            .iter()
            .zip(&e.to)
            .all(|(&a, &b)| within(a) == within(b)));
        let expect_block = |s: usize| if directed { s * (s - 1) } else { s * (s - 1) / 2 } + if loops { s } else { 0 };
        assert_eq!(e.len(), expect_block(sizes[0]) + expect_block(sizes[1]));
        let probs = [0.3, 0.1, 0.1, 0.3];
        let lists = generate::sbm(&sizes, &probs, directed, loops, 2, seed, 2).unwrap();
        assert_eq!(
            lists,
            generate::sbm(&sizes, &probs, directed, loops, 2, seed, 1).unwrap()
        );
    }
    // Invalid requests are errors.
    assert!(generate::gnm(5, 11, false, false, 1, 1, 1).is_err());
    assert!(generate::gnp(5, 1.5, false, false, 1, 1, 1).is_err());
    assert!(generate::uman(4, 4, 3, 1, 1, 1).is_err());
    assert!(generate::small_world(10, 3, 0.1, 1, 1, 1).is_err());
    assert!(generate::preferential(5, 5, 1, 1, 1).is_err());
    assert!(
        generate::configuration(&Degrees::Undirected(vec![1, 1, 1]), false, 1, 1, 1, 1).is_err()
    );
    assert!(generate::sbm(&[2, 2], &[0.5, 0.2, 0.3, 0.5], false, false, 1, 1, 1).is_err());
}
