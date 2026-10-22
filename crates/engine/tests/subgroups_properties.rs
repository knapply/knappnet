//! Randomized property tests for cores, cliques, and communities.

extern crate knappnet as engine;

use std::collections::BTreeSet;

use engine::metrics::cliques::{clique_number, max_cliques};
use engine::metrics::community::{louvain, modularity};
use engine::metrics::cores::coreness;
use engine::metrics::degree::DegreeMode;
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
    weights: Vec<f64>,
    directed: bool,
}

fn random_case(rng: &mut Rng, max_n: u64) -> Case {
    let n = 1 + rng.below(max_n) as usize;
    let m = rng.below((n * n / 3 + 1) as u64) as usize;
    let mut c = Case {
        n,
        src: Vec::with_capacity(m),
        tgt: Vec::with_capacity(m),
        weights: Vec::with_capacity(m),
        directed: rng.below(2) == 0,
    };
    for _ in 0..m {
        c.src.push(rng.below(n as u64) as i64);
        c.tgt.push(rng.below(n as u64) as i64);
        c.weights.push(1.0 + rng.below(5) as f64 / 2.0);
    }
    c
}

fn graph(c: &Case, weighted: bool, directed: bool) -> Graph {
    let w = if weighted { Some(&c.weights[..]) } else { None };
    Graph::from_edges(c.n, &c.src, &c.tgt, w, directed).unwrap()
}

fn permutation(rng: &mut Rng, n: usize) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        p.swap(i, rng.below(i as u64 + 1) as usize);
    }
    p
}

fn permuted(c: &Case, perm: &[usize]) -> Case {
    Case {
        n: c.n,
        src: c.src.iter().map(|&v| perm[v as usize] as i64).collect(),
        tgt: c.tgt.iter().map(|&v| perm[v as usize] as i64).collect(),
        weights: c.weights.clone(),
        directed: c.directed,
    }
}

/// Simple undirected adjacency without loops.
fn adjacency(c: &Case) -> Vec<BTreeSet<usize>> {
    let mut adj = vec![BTreeSet::new(); c.n];
    for (&a, &b) in c.src.iter().zip(&c.tgt) {
        if a != b {
            adj[a as usize].insert(b as usize);
            adj[b as usize].insert(a as usize);
        }
    }
    adj
}

const CASES: u64 = 300;

#[test]
fn every_core_is_a_core() {
    // The nodes with coreness >= k induce a subgraph of minimum degree >= k,
    // and the nodes with coreness exactly k cannot all be kept in a
    // (k+1)-core (checked by maximality: the k-core is the largest such set).
    let mut rng = Rng(11);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 30);
        let g = graph(&c, false, c.directed);
        for mode in [DegreeMode::All, DegreeMode::Out, DegreeMode::In] {
            let core = coreness(&g, mode);
            let kmax = core.iter().copied().max().unwrap_or(0);
            for k in 0..=kmax + 1 {
                let inside: Vec<bool> = core.iter().map(|&x| x >= k).collect();
                let mut deg = vec![0usize; c.n];
                for (&a, &b) in c.src.iter().zip(&c.tgt) {
                    let (a, b) = (a as usize, b as usize);
                    if !(inside[a] && inside[b]) {
                        continue;
                    }
                    let all = !c.directed || mode == DegreeMode::All;
                    if all || mode == DegreeMode::Out {
                        deg[a] += 1;
                    }
                    if all || mode == DegreeMode::In {
                        deg[b] += 1;
                    }
                }
                for v in 0..c.n {
                    if inside[v] {
                        assert!(
                            deg[v] >= k,
                            "node {v} in the {k}-core has degree {}",
                            deg[v]
                        );
                    }
                }
            }
            // Peeling from scratch: nodes outside the k-core must be
            // removable by repeatedly deleting nodes of degree < k.
            let k = kmax;
            let mut alive = vec![true; c.n];
            loop {
                let mut deg = vec![0usize; c.n];
                for (&a, &b) in c.src.iter().zip(&c.tgt) {
                    let (a, b) = (a as usize, b as usize);
                    if alive[a] && alive[b] {
                        let all = !c.directed || mode == DegreeMode::All;
                        if all || mode == DegreeMode::Out {
                            deg[a] += 1;
                        }
                        if all || mode == DegreeMode::In {
                            deg[b] += 1;
                        }
                    }
                }
                let drop: Vec<usize> = (0..c.n).filter(|&v| alive[v] && deg[v] < k).collect();
                if drop.is_empty() {
                    break;
                }
                for v in drop {
                    alive[v] = false;
                }
            }
            for v in 0..c.n {
                assert_eq!(alive[v], core[v] >= k);
            }
        }
    }
}

#[test]
fn cliques_are_maximal_complete_and_exhaustive() {
    let mut rng = Rng(12);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 24);
        let g = graph(&c, false, c.directed);
        let adj = adjacency(&c);
        let cliques = max_cliques(&g, 0, 1 + rng.below(4) as usize);
        let mut covered = vec![false; c.n];
        for k in &cliques {
            assert!(k.windows(2).all(|w| w[0] < w[1]));
            for (i, &a) in k.iter().enumerate() {
                covered[a as usize] = true;
                for &b in &k[i + 1..] {
                    assert!(adj[a as usize].contains(&(b as usize)));
                }
            }
            // No outside node extends the clique.
            for (v, nbrs) in adj.iter().enumerate() {
                if !k.contains(&(v as u32)) {
                    assert!(!k.iter().all(|&a| nbrs.contains(&(a as usize))));
                }
            }
        }
        // Every node and every edge lies in some maximal clique.
        assert!(covered.iter().all(|&x| x));
        for (a, nbrs) in adj.iter().enumerate() {
            for &b in nbrs {
                assert!(cliques
                    .iter()
                    .any(|k| k.contains(&(a as u32)) && k.contains(&(b as u32))));
            }
        }
        let set: BTreeSet<_> = cliques.iter().collect();
        assert_eq!(set.len(), cliques.len(), "duplicate cliques");
        let largest = cliques.iter().map(|k| k.len()).max().unwrap_or(0);
        assert_eq!(clique_number(&g, 0), largest);
        // min_size filters and nothing else.
        let big: Vec<_> = cliques.iter().filter(|k| k.len() >= 3).cloned().collect();
        assert_eq!(max_cliques(&g, 3, 0), big);
    }
}

#[test]
fn cliques_do_not_depend_on_labels_or_threads() {
    let mut rng = Rng(13);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 30);
        let perm = permutation(&mut rng, c.n);
        let g = graph(&c, false, c.directed);
        let h = graph(&permuted(&c, &perm), false, c.directed);
        let a = max_cliques(&g, 0, 1);
        assert_eq!(a, max_cliques(&g, 0, 3));
        let mut mapped: Vec<Vec<u32>> = a
            .iter()
            .map(|k| {
                let mut m: Vec<u32> = k.iter().map(|&v| perm[v as usize] as u32).collect();
                m.sort_unstable();
                m
            })
            .collect();
        mapped.sort_unstable();
        assert_eq!(mapped, max_cliques(&h, 0, 2));
        let core_g = coreness(&g, DegreeMode::All);
        let core_h = coreness(&h, DegreeMode::All);
        for v in 0..c.n {
            assert_eq!(core_g[v], core_h[perm[v]]);
        }
    }
}

/// Modularity by the double-sum definition, independently of the engine.
fn modularity_by_definition(c: &Case, memb: &[i64], gamma: f64, weighted: bool) -> f64 {
    let n = c.n;
    let mut a = vec![vec![0.0; n]; n];
    for (e, (&s, &t)) in c.src.iter().zip(&c.tgt).enumerate() {
        let w = if weighted { c.weights[e] } else { 1.0 };
        let (s, t) = (s as usize, t as usize);
        if c.directed {
            a[s][t] += w;
        } else if s == t {
            a[s][s] += 2.0 * w;
        } else {
            a[s][t] += w;
            a[t][s] += w;
        }
    }
    let kout: Vec<f64> = a.iter().map(|r| r.iter().sum()).collect();
    let kin: Vec<f64> = (0..n).map(|j| (0..n).map(|i| a[i][j]).sum()).collect();
    let total: f64 = kout.iter().sum();
    let mut q = 0.0;
    for i in 0..n {
        for j in 0..n {
            if memb[i] == memb[j] {
                q += a[i][j] - gamma * kout[i] * kin[j] / total;
            }
        }
    }
    q / total
}

#[test]
fn modularity_matches_the_double_sum() {
    let mut rng = Rng(14);
    for _ in 0..CASES {
        let c = random_case(&mut rng, 25);
        if c.src.is_empty() {
            continue;
        }
        let groups = 1 + rng.below(5);
        let memb: Vec<i64> = (0..c.n).map(|_| rng.below(groups) as i64).collect();
        let gamma = rng.below(5) as f64 / 2.0;
        for weighted in [false, true] {
            let g = graph(&c, weighted, c.directed);
            let q = modularity(&g, &memb, gamma, weighted).unwrap();
            let d = modularity_by_definition(&c, &memb, gamma, weighted);
            assert!((q - d).abs() <= 1e-12 * (1.0 + d.abs()), "{q} vs {d}");
        }
    }
}

#[test]
fn louvain_is_seeded_and_locally_optimal() {
    let mut rng = Rng(15);
    for _ in 0..CASES {
        let mut c = random_case(&mut rng, 30);
        c.directed = false;
        let weighted = rng.below(2) == 0;
        let g = graph(&c, weighted, false);
        let gamma = [0.5, 1.0, 1.5][rng.below(3) as usize];
        let seed = rng.next();
        let r = louvain(&g, gamma, weighted, seed).unwrap();
        let again = louvain(&g, gamma, weighted, seed).unwrap();
        assert_eq!(r.membership, again.membership);
        assert_eq!(r.modularity.to_bits(), again.modularity.to_bits());
        // Labels are 0..k by first appearance.
        let mut next = 0;
        for &x in &r.membership {
            assert!(x <= next);
            if x == next {
                next += 1;
            }
        }
        if c.src.is_empty() {
            assert!(r.modularity.is_nan());
            continue;
        }
        let memb: Vec<i64> = r.membership.iter().map(|&x| x as i64).collect();
        assert_eq!(
            r.modularity,
            modularity(&g, &memb, gamma, weighted).unwrap()
        );
        // The last pass, run on the network of final groups, moved nothing:
        // merging two adjacent groups never raises Q. (Moving a single
        // node can, since Louvain only moves single nodes on the first
        // level.)
        let adj = adjacency(&c);
        for v in 0..c.n {
            for &u in &adj[v] {
                if memb[u] == memb[v] {
                    continue;
                }
                let merged: Vec<i64> = memb
                    .iter()
                    .map(|&x| if x == memb[v] { memb[u] } else { x })
                    .collect();
                let q = modularity(&g, &merged, gamma, weighted).unwrap();
                assert!(q <= r.modularity + 1e-9, "merging groups raises Q");
            }
        }
        // Never worse than putting every node in its own group.
        let singletons: Vec<i64> = (0..c.n as i64).collect();
        assert!(r.modularity >= modularity(&g, &singletons, gamma, weighted).unwrap() - 1e-12);
    }
}
