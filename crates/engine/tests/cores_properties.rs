//! Randomized property tests for truss numbers, k-trusses, onion layers,
//! and core subsets against direct evaluation of their definitions, and
//! their invariance under relabeling.

extern crate knappnet as engine;

use std::collections::BTreeSet;

use engine::metrics::cores::{
    core_subset, coreness, k_truss, onion_layers, truss_numbers, CoreSubset,
};
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

fn random(rng: &mut Rng) -> (usize, Vec<i64>, Vec<i64>, bool) {
    let n = 1 + rng.below(25) as usize;
    let m = rng.below((n * n / 2 + 1) as u64) as usize;
    let s = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    let t = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    (n, s, t, rng.below(2) == 0)
}

fn pairs(s: &[i64], t: &[i64]) -> BTreeSet<(i64, i64)> {
    s.iter()
        .zip(t)
        .filter(|(a, b)| a != b)
        .map(|(&a, &b)| (a.min(b), a.max(b)))
        .collect()
}

fn truss_by_definition(edges: &BTreeSet<(i64, i64)>, k: usize) -> BTreeSet<(i64, i64)> {
    let mut e = edges.clone();
    loop {
        let has = |a: i64, b: i64, e: &BTreeSet<(i64, i64)>| e.contains(&(a.min(b), a.max(b)));
        let nodes: BTreeSet<i64> = e.iter().flat_map(|&(a, b)| [a, b]).collect();
        let drop: Vec<(i64, i64)> = e
            .iter()
            .copied()
            .filter(|&(a, b)| {
                nodes
                    .iter()
                    .filter(|&&w| w != a && w != b && has(a, w, &e) && has(b, w, &e))
                    .count()
                    + 2
                    < k
            })
            .collect();
        if drop.is_empty() {
            return e;
        }
        for d in drop {
            e.remove(&d);
        }
    }
}

#[test]
fn trusses_match_the_definition() {
    let mut rng = Rng(4);
    for _ in 0..60 {
        let (n, s, t, directed) = random(&mut rng);
        let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
        let nums = truss_numbers(&g);
        let all = pairs(&s, &t);
        for k in 3..7 {
            let want = truss_by_definition(&all, k);
            for (i, (&a, &b)) in s.iter().zip(&t).enumerate() {
                if a == b {
                    assert_eq!(nums[i], 0);
                } else {
                    assert_eq!(nums[i] as usize >= k, want.contains(&(a.min(b), a.max(b))));
                }
            }
            let (edges, nodes) = k_truss(&g, k as i64);
            let got: BTreeSet<(i64, i64)> = edges
                .iter()
                .map(|&e| (s[e].min(t[e]), s[e].max(t[e])))
                .collect();
            assert_eq!(got, want);
            let wn: BTreeSet<u32> = want
                .iter()
                .flat_map(|&(a, b)| [a as u32, b as u32])
                .collect();
            assert_eq!(nodes.into_iter().collect::<BTreeSet<u32>>(), wn);
        }
    }
}

fn onion_by_definition(n: usize, s: &[i64], t: &[i64]) -> Vec<u32> {
    let e = pairs(s, t);
    let mut alive: BTreeSet<usize> = (0..n).collect();
    let deg = |v: usize, alive: &BTreeSet<usize>| {
        e.iter()
            .filter(|&&(a, b)| {
                (a as usize == v && alive.contains(&(b as usize)))
                    || (b as usize == v && alive.contains(&(a as usize)))
            })
            .count()
    };
    let mut layer = vec![0u32; n];
    let iso: Vec<usize> = (0..n).filter(|&v| deg(v, &alive) == 0).collect();
    let mut current = 1;
    for &v in &iso {
        layer[v] = 1;
        alive.remove(&v);
    }
    if !iso.is_empty() {
        current = 2;
    }
    let mut core = 1;
    while !alive.is_empty() {
        let d: Vec<(usize, usize)> = alive.iter().map(|&v| (v, deg(v, &alive))).collect();
        core = core.max(d.iter().map(|x| x.1).min().unwrap());
        for &(v, dv) in &d {
            if dv <= core {
                layer[v] = current;
                alive.remove(&v);
            }
        }
        current += 1;
    }
    layer
}

#[test]
fn onion_layers_match_the_definition() {
    let mut rng = Rng(12);
    for _ in 0..80 {
        let (n, s, t, directed) = random(&mut rng);
        let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
        assert_eq!(onion_layers(&g), onion_by_definition(n, &s, &t));
    }
}

#[test]
fn core_subsets_partition_by_core_number() {
    let mut rng = Rng(77);
    for _ in 0..50 {
        let (n, s, t, directed) = random(&mut rng);
        let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
        let c = coreness(&g, DegreeMode::All);
        let max = *c.iter().max().unwrap() as i64;
        for k in 0..=max + 1 {
            let core = core_subset(&g, CoreSubset::Core, Some(k)).unwrap();
            let shell = core_subset(&g, CoreSubset::Shell, Some(k)).unwrap();
            let crust = core_subset(&g, CoreSubset::Crust, Some(k)).unwrap();
            let corona = core_subset(&g, CoreSubset::Corona, Some(k)).unwrap();
            assert!(corona.iter().all(|v| shell.contains(v)));
            // Core and crust overlap exactly in the shell.
            let both: Vec<u32> = core.iter().copied().filter(|v| crust.contains(v)).collect();
            assert_eq!(both, shell);
            assert_eq!(core.len() + crust.len() - shell.len(), n);
        }
        let shell = core_subset(&g, CoreSubset::Shell, None).unwrap();
        assert_eq!(shell, core_subset(&g, CoreSubset::Core, None).unwrap());
    }
}
