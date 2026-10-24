//! Randomized property tests for rewiring: every swap keeps degrees and
//! simplicity (and connectivity where required), results are reproducible
//! for a seed, and the means over random graphs do not depend on the
//! thread count.

extern crate knappnet as engine;

use std::collections::{BTreeMap, BTreeSet};

use engine::metrics::rewiring::{
    connected_double_edge_swap, directed_edge_swap, double_edge_swap, lattice_reference, omega,
    random_reference, rich_club_coefficient, sigma, SwapGraph,
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

/// A connected random graph: a random tree plus extra edges.
fn connected_graph(rng: &mut Rng, directed: bool) -> Graph {
    let n = 5 + rng.below(20) as usize;
    let mut s = Vec::new();
    let mut t = Vec::new();
    for v in 1..n {
        s.push(rng.below(v as u64) as i64);
        t.push(v as i64);
    }
    for _ in 0..rng.below(3 * n as u64) {
        let (a, b) = (rng.below(n as u64) as i64, rng.below(n as u64) as i64);
        if a != b {
            s.push(a);
            t.push(b);
        }
    }
    Graph::from_edges(n, &s, &t, None, directed).unwrap()
}

fn degree_map(edges: &[(u32, u32)], directed: bool) -> BTreeMap<(u32, bool), usize> {
    let mut d = BTreeMap::new();
    for &(a, b) in edges {
        *d.entry((a, false)).or_insert(0) += 1;
        *d.entry((b, directed)).or_insert(0) += 1;
    }
    d
}

fn connected(n: usize, edges: &[(u32, u32)]) -> bool {
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in edges {
        adj[a as usize].push(b as usize);
        adj[b as usize].push(a as usize);
    }
    let mut seen = vec![false; n];
    let mut stack = vec![0];
    seen[0] = true;
    while let Some(v) = stack.pop() {
        for &w in &adj[v] {
            if !seen[w] {
                seen[w] = true;
                stack.push(w);
            }
        }
    }
    seen.iter().all(|&x| x)
}

#[test]
fn swaps_keep_degrees_simplicity_and_connectivity() {
    let mut rng = Rng(7);
    for round in 0..40 {
        let g = connected_graph(&mut rng, false);
        let n = g.node_count();
        let base = SwapGraph::new(&g).edge_list();
        let want = degree_map(&base, false);
        let check = |e: &[(u32, u32)], must_connect: bool| {
            assert_eq!(degree_map(e, false), want);
            assert_eq!(e.iter().collect::<BTreeSet<_>>().len(), e.len());
            if must_connect {
                assert!(connected(n, e));
            }
        };
        if let Ok(e) = double_edge_swap(&g, 3, 3000, round) {
            check(&e, false);
            assert_eq!(e, double_edge_swap(&g, 3, 3000, round).unwrap());
        }
        let (e, k) = connected_double_edge_swap(&g, 30, round).unwrap();
        check(&e, true);
        assert!(k <= 30);
        check(&random_reference(&g, 1, true, round).unwrap(), true);
        check(&lattice_reference(&g, 1, None, true, round).unwrap(), true);
        let d = connected_graph(&mut rng, true);
        let base = SwapGraph::new(&d).edge_list();
        if let Ok(e) = directed_edge_swap(&d, 3, 3000, round) {
            assert_eq!(degree_map(&e, true), degree_map(&base, true));
            assert_eq!(e.iter().collect::<BTreeSet<_>>().len(), e.len());
        }
    }
}

#[test]
fn random_means_do_not_depend_on_threads() {
    let mut rng = Rng(99);
    for _ in 0..4 {
        let g = connected_graph(&mut rng, false);
        let a = sigma(&g, 2, 4, 5, 1).unwrap();
        let b = sigma(&g, 2, 4, 5, 3).unwrap();
        assert!(a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()));
        let a = omega(&g, 1, 3, 5, 1).unwrap();
        let b = omega(&g, 1, 3, 5, 4).unwrap();
        assert!(a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()));
        if let (Ok(a), Ok(b)) = (
            rich_club_coefficient(&g, true, 5, 3, 1, 1),
            rich_club_coefficient(&g, true, 5, 3, 1, 4),
        ) {
            assert_eq!(a, b);
        }
    }
}
