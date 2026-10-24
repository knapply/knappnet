//! Randomized property tests for VF2++ and tree isomorphism: every
//! mapping is valid and distinct, counts match brute force on small
//! graphs (multigraphs, loops, labels, both directions), and relabeling
//! preserves counts.

extern crate knappnet as engine;

use std::collections::{BTreeMap, BTreeSet};

use engine::metrics::isomorphism::{
    could_be_isomorphic, rooted_tree_isomorphism, tree_isomorphism, Matcher, Problem,
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
}

#[derive(Clone)]
struct Case {
    n: usize,
    edges: Vec<(i64, i64)>,
    directed: bool,
}

impl Case {
    fn graph(&self) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = self.edges.iter().copied().unzip();
        Graph::from_edges(self.n, &s, &t, None, self.directed).unwrap()
    }
    fn counts(&self) -> BTreeMap<(usize, usize), u32> {
        let mut c = BTreeMap::new();
        for &(a, b) in &self.edges {
            let (a, b) = (a as usize, b as usize);
            let k = if self.directed || a <= b {
                (a, b)
            } else {
                (b, a)
            };
            *c.entry(k).or_insert(0) += 1;
        }
        c
    }
}

fn random_case(rng: &mut Rng, max_n: u64, directed: bool) -> Case {
    let n = 1 + rng.below(max_n) as usize;
    let m = rng.below((n * n / 2 + 2) as u64) as usize;
    Case {
        n,
        edges: (0..m)
            .map(|_| (rng.below(n as u64) as i64, rng.below(n as u64) as i64))
            .collect(),
        directed,
    }
}

fn valid(m: &[u32], sg: &Case, fg: &Case, p: Problem, ls: &[i64], lf: &[i64]) -> bool {
    let (cs, cf) = (sg.counts(), fg.counts());
    let mut seen = BTreeSet::new();
    for (u, &v) in m.iter().enumerate() {
        if !seen.insert(v) || ls[u] != lf[v as usize] {
            return false;
        }
    }
    for u in 0..sg.n {
        for w in 0..sg.n {
            if !sg.directed && w < u {
                continue;
            }
            let (x, y) = (m[u] as usize, m[w] as usize);
            let kf = if fg.directed || x <= y {
                (x, y)
            } else {
                (y, x)
            };
            let s = cs.get(&(u, w)).copied().unwrap_or(0);
            let f = cf.get(&kf).copied().unwrap_or(0);
            let ok = match p {
                Problem::Monomorphism => s <= f,
                _ => s == f,
            };
            if !ok {
                return false;
            }
        }
    }
    true
}

/// Every injective map, checked against the definition.
fn brute(sg: &Case, fg: &Case, p: Problem, ls: &[i64], lf: &[i64]) -> usize {
    fn rec(
        m: &mut Vec<u32>,
        sg: &Case,
        fg: &Case,
        p: Problem,
        ls: &[i64],
        lf: &[i64],
        total: &mut usize,
    ) {
        if m.len() == sg.n {
            if valid(m, sg, fg, p, ls, lf) {
                *total += 1;
            }
            return;
        }
        for x in 0..fg.n as u32 {
            if !m.contains(&x) {
                m.push(x);
                rec(m, sg, fg, p, ls, lf, total);
                m.pop();
            }
        }
    }
    let mut total = 0;
    rec(&mut Vec::new(), sg, fg, p, ls, lf, &mut total);
    total
}

#[test]
fn mappings_match_brute_force() {
    let mut rng = Rng(42);
    for round in 0..150 {
        let directed = round % 2 == 0;
        let fg = random_case(&mut rng, 6, directed);
        let mut sg = random_case(&mut rng, 4, directed);
        if round % 5 == 0 {
            // A relabeled copy, so isomorphisms exist.
            let mut perm: Vec<i64> = (0..fg.n as i64).collect();
            for i in (1..fg.n).rev() {
                perm.swap(i, rng.below(i as u64 + 1) as usize);
            }
            sg = Case {
                n: fg.n,
                edges: fg
                    .edges
                    .iter()
                    .map(|&(a, b)| (perm[a as usize], perm[b as usize]))
                    .collect(),
                directed,
            };
        }
        let labelled = round % 3 == 0;
        let lf: Vec<i64> = (0..fg.n)
            .map(|v| if labelled { (v % 2) as i64 } else { 0 })
            .collect();
        let ls: Vec<i64> = (0..sg.n)
            .map(|v| if labelled { (v % 2) as i64 } else { 0 })
            .collect();
        let (gf, gs) = (fg.graph(), sg.graph());
        for p in [
            Problem::Isomorphism,
            Problem::InducedSubgraph,
            Problem::Monomorphism,
        ] {
            let mut m = Matcher::new(&gf, &gs, Some(&lf), Some(&ls), p, false).unwrap();
            let mut found = BTreeSet::new();
            while let Some(x) = m.next() {
                assert!(
                    valid(&x, &sg, &fg, p, &ls, &lf),
                    "{p:?} invalid mapping {x:?}"
                );
                assert!(found.insert(x), "{p:?} repeated mapping");
            }
            if p == Problem::Isomorphism && fg.n != sg.n {
                assert!(found.is_empty());
                continue;
            }
            assert_eq!(found.len(), brute(&sg, &fg, p, &ls, &lf), "{p:?}");
        }
    }
}

#[test]
fn tree_isomorphisms_are_valid() {
    let mut rng = Rng(8);
    for _ in 0..200 {
        let n = 1 + rng.below(12) as usize;
        let parent: Vec<i64> = (1..n).map(|v| rng.below(v as u64) as i64).collect();
        let edges: Vec<(i64, i64)> = (1..n).map(|v| (parent[v - 1], v as i64)).collect();
        let mut perm: Vec<i64> = (0..n as i64).collect();
        for i in (1..n).rev() {
            perm.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let a = Case {
            n,
            edges: edges.clone(),
            directed: false,
        };
        let b = Case {
            n,
            edges: edges
                .iter()
                .map(|&(x, y)| (perm[y as usize], perm[x as usize]))
                .collect(),
            directed: false,
        };
        let (ga, gb) = (a.graph(), b.graph());
        let zero = vec![0; n];
        let m = tree_isomorphism(&ga, &gb).unwrap().expect("isomorphic");
        assert!(valid(&m, &a, &b, Problem::Isomorphism, &zero, &zero));
        let r = rooted_tree_isomorphism(&ga, 0, &gb, perm[0] as usize)
            .unwrap()
            .expect("rooted");
        assert_eq!(r[0], perm[0] as u32);
        assert!(valid(&r, &a, &b, Problem::Isomorphism, &zero, &zero));
        // Agreement with VF2++ on a random other tree.
        let parent2: Vec<i64> = (1..n).map(|v| rng.below(v as u64) as i64).collect();
        let c = Case {
            n,
            edges: (1..n).map(|v| (parent2[v - 1], v as i64)).collect(),
            directed: false,
        };
        let gc = c.graph();
        let vf = Matcher::new(&ga, &gc, None, None, Problem::Isomorphism, true)
            .unwrap()
            .next()
            .is_some();
        assert_eq!(tree_isomorphism(&ga, &gc).unwrap().is_some(), vf);
        assert!(could_be_isomorphic(&ga, &gb, true, true, true, 1).unwrap());
    }
}
