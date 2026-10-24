//! Randomized checks of connectivity, cuts, disjoint paths, k-components,
//! and k-edge-components against brute force on small networks, and of
//! thread-count and relabelling invariance.

extern crate knappnet as engine;

use engine::metrics::connectivity::{
    all_pairs_node_connectivity, condensation, edge_connectivity, edge_disjoint_paths,
    is_semiconnected, k_components, k_edge_components, local_edge_connectivity,
    local_node_connectivity, minimum_edge_cut, minimum_node_cut, minimum_st_edge_cut,
    minimum_st_node_cut, node_connectivity, node_disjoint_paths, CutConvention,
};
use engine::order::NeighborOrder;
use engine::rng::Stream;
use engine::Graph;

fn random_graph(rng: &mut Stream, n: usize, p: f64, directed: bool) -> Graph {
    let mut s = Vec::new();
    let mut t = Vec::new();
    for a in 0..n {
        for b in 0..n {
            if (!directed && b < a) || rng.uniform() >= p {
                continue;
            }
            s.push(a as i64);
            t.push(b as i64);
            if rng.uniform() < 0.1 {
                s.push(a as i64);
                t.push(b as i64);
            }
        }
    }
    Graph::from_edges(n, &s, &t, None, directed).unwrap()
}

/// Simple adjacency as bit masks of successors (neighbours if undirected).
fn masks(g: &Graph) -> Vec<u32> {
    let n = g.node_count();
    let mut m = vec![0u32; n];
    for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
        if a != b {
            m[a as usize] |= 1 << b;
            if !g.is_directed() {
                m[b as usize] |= 1 << a;
            }
        }
    }
    m
}

/// Nodes reachable from `s` avoiding `removed`, without the arc s -> t
/// when `skip` is set.
fn reach(adj: &[u32], s: usize, removed: u32, skip: Option<usize>) -> u32 {
    let mut seen = 1u32 << s;
    let mut frontier = seen;
    while frontier != 0 {
        let mut next = 0u32;
        let mut f = frontier;
        while f != 0 {
            let v = f.trailing_zeros() as usize;
            f &= f - 1;
            let mut out = adj[v] & !removed;
            if v == s {
                if let Some(t) = skip {
                    out &= !(1 << t);
                }
            }
            next |= out;
        }
        next &= !seen;
        seen |= next;
        frontier = next;
    }
    seen
}

fn brute_kappa(adj: &[u32], n: usize, s: usize, t: usize) -> u64 {
    let adjacent = adj[s] >> t & 1 == 1;
    let mut best = u32::MAX;
    for mask in 0u32..(1 << n) {
        if mask >> s & 1 == 1 || mask >> t & 1 == 1 || mask.count_ones() >= best {
            continue;
        }
        if reach(adj, s, mask, Some(t)) >> t & 1 == 0 {
            best = mask.count_ones();
        }
    }
    best as u64 + adjacent as u64
}

fn brute_lambda(g: &Graph, s: usize, t: usize) -> u64 {
    let n = g.node_count();
    let adj = masks(g);
    let mut best = u64::MAX;
    for mask in 0u32..(1 << n) {
        if mask >> s & 1 == 0 || mask >> t & 1 == 1 {
            continue;
        }
        let c: u64 = (0..n)
            .filter(|&v| mask >> v & 1 == 1)
            .map(|v| (adj[v] & !mask).count_ones() as u64)
            .sum();
        best = best.min(c);
    }
    best
}

fn strongly_connected(adj: &[u32], n: usize, removed: u32) -> bool {
    let alive = ((1u64 << n) - 1) as u32 & !removed;
    if alive.count_ones() <= 1 {
        return true;
    }
    let r = alive.trailing_zeros() as usize;
    let mut v = alive;
    while v != 0 {
        let x = v.trailing_zeros() as usize;
        v &= v - 1;
        if reach(adj, r, removed, None) >> x & 1 == 0 || reach(adj, x, removed, None) >> r & 1 == 0
        {
            return false;
        }
    }
    true
}

fn brute_global_kappa(adj: &[u32], n: usize) -> u64 {
    let complete = (0..n).all(|v| (0..n).all(|w| v == w || adj[v] >> w & 1 == 1));
    if complete {
        return n as u64 - 1;
    }
    (0u32..(1 << n))
        .filter(|&m| !strongly_connected(adj, n, m))
        .map(|m| m.count_ones() as u64)
        .min()
        .unwrap()
}

#[test]
fn connectivity_matches_brute_force() {
    let mut rng = Stream::new(201, 0);
    for round in 0..150 {
        let n = 2 + round % 8;
        let directed = round % 2 == 0;
        let g = random_graph(&mut rng, n, 0.25 + 0.05 * (round % 8) as f64, directed);
        let adj = masks(&g);
        let (s, t) = (0, n - 1);
        let kappa = brute_kappa(&adj, n, s, t);
        assert_eq!(
            local_node_connectivity(&g, s, t).unwrap(),
            kappa,
            "round {round}"
        );
        let lambda = brute_lambda(&g, s, t);
        assert_eq!(local_edge_connectivity(&g, s, t).unwrap(), lambda);
        let global = brute_global_kappa(&adj, n);
        for threads in [1, 3] {
            assert_eq!(
                node_connectivity(&g, threads).unwrap(),
                global,
                "round {round}"
            );
        }
        let global_lambda = (0..n)
            .flat_map(|a| (0..n).map(move |b| (a, b)))
            .filter(|(a, b)| a != b)
            .map(|(a, b)| brute_lambda(&g, a, b))
            .min()
            .unwrap_or(0);
        assert_eq!(edge_connectivity(&g, 2).unwrap(), global_lambda);
        // Menger: as many disjoint paths as the connectivity, and valid.
        let paths = node_disjoint_paths(&g, s, t).unwrap();
        assert_eq!(paths.len() as u64, kappa);
        let mut inner = 0u32;
        for p in &paths {
            assert_eq!((p[0], *p.last().unwrap()), (s, t));
            for w in p.windows(2) {
                assert!(adj[w[0]] >> w[1] & 1 == 1);
            }
            for &x in &p[1..p.len() - 1] {
                assert!(inner >> x & 1 == 0);
                inner |= 1 << x;
            }
        }
        assert_eq!(edge_disjoint_paths(&g, s, t).unwrap().len() as u64, lambda);
        // Cuts: valid and minimum.
        let cut = minimum_st_edge_cut(&g, s, t).unwrap();
        assert_eq!(cut.len() as u64, lambda);
        for conv in [CutConvention::Engine, CutConvention::Networkx] {
            match minimum_st_node_cut(&g, s, t, conv) {
                Ok(c) if adj[s] >> t & 1 == 0 => {
                    assert_eq!(c.len() as u64, kappa);
                    let removed = c.iter().fold(0u32, |a, &v| a | 1 << v);
                    assert_eq!(reach(&adj, s, removed, None) >> t & 1, 0);
                }
                Ok(c) => {
                    assert_eq!(conv, CutConvention::Networkx);
                    assert!(c.is_empty());
                }
                Err(_) => assert!(conv == CutConvention::Engine && adj[s] >> t & 1 == 1),
            }
            let c = minimum_node_cut(&g, conv, 2).unwrap();
            let removed = c.iter().fold(0u32, |a, &v| a | 1 << v);
            assert_eq!(c.len() as u64, global, "round {round}");
            assert!(!strongly_connected(&adj, n, removed) || c.len() == n - 1);
        }
        let ec = minimum_edge_cut(&g, 1).unwrap();
        assert_eq!(ec.len() as u64, global_lambda);
    }
}

#[test]
fn all_pairs_is_thread_and_label_independent() {
    let mut rng = Stream::new(202, 0);
    for round in 0..20 {
        let n = 3 + round % 7;
        let directed = round % 2 == 1;
        let g = random_graph(&mut rng, n, 0.4, directed);
        let a = all_pairs_node_connectivity(&g, 1);
        let b = all_pairs_node_connectivity(&g, 4);
        assert_eq!(a.as_slice(), b.as_slice());
        let mut perm: Vec<usize> = (0..n).collect();
        rng.permutation(&mut perm);
        let s: Vec<i64> = g
            .edge_sources()
            .iter()
            .map(|&v| perm[v as usize] as i64)
            .collect();
        let t: Vec<i64> = g
            .edge_targets()
            .iter()
            .map(|&v| perm[v as usize] as i64)
            .collect();
        let h = Graph::from_edges(n, &s, &t, None, directed).unwrap();
        let c = all_pairs_node_connectivity(&h, 2);
        for u in 0..n {
            for v in 0..n {
                assert_eq!(a.get(u, v), c.get(perm[u], perm[v]));
            }
        }
    }
}

#[test]
fn k_components_are_maximal_k_connected_sets() {
    let mut rng = Stream::new(203, 0);
    for round in 0..60 {
        let n = 2 + round % 8;
        let g = random_graph(&mut rng, n, 0.3 + 0.07 * (round % 6) as f64, false);
        let adj = masks(&g);
        // kappa of every induced subgraph with at least two nodes.
        let mut kappa = vec![0u64; 1 << n];
        for set in 1u32..(1 << n) {
            if set.count_ones() < 2 {
                continue;
            }
            let nodes: Vec<usize> = (0..n).filter(|&v| set >> v & 1 == 1).collect();
            let sub: Vec<u32> = (0..n).map(|v| adj[v] & set).collect();
            let complete = nodes
                .iter()
                .all(|&v| nodes.iter().all(|&w| v == w || sub[v] >> w & 1 == 1));
            kappa[set as usize] = if complete {
                nodes.len() as u64 - 1
            } else {
                (0u32..(1 << n))
                    .filter(|&m| m & !set == 0)
                    .filter(|&m| !strongly_connected(&sub, n, m | !set))
                    .map(|m| m.count_ones() as u64)
                    .min()
                    .unwrap()
            };
        }
        let got = k_components(&g, 2).unwrap();
        let kmax = kappa.iter().copied().max().unwrap_or(0) as usize;
        assert_eq!(got.len(), kmax, "round {round}");
        for (i, level) in got.iter().enumerate() {
            let k = (i + 1) as u64;
            let sets: Vec<u32> = (1u32..(1 << n))
                .filter(|&s| kappa[s as usize] >= k && s.count_ones() as u64 > k)
                .collect();
            let mut want: Vec<Vec<u32>> = sets
                .iter()
                .filter(|&&s| !sets.iter().any(|&o| o != s && o & s == s))
                .map(|&s| (0..n as u32).filter(|&v| s >> v & 1 == 1).collect())
                .collect();
            want.sort();
            assert_eq!(level, &want, "round {round}, k = {k}");
        }
    }
}

#[test]
fn k_edge_components_match_pairwise_connectivity() {
    let mut rng = Stream::new(204, 0);
    for round in 0..80 {
        let n = 1 + round % 8;
        let directed = round % 2 == 0;
        let g = random_graph(&mut rng, n, 0.45, directed);
        for k in 1..4u64 {
            let labels = k_edge_components(&g, k, 2).unwrap();
            for u in 0..n {
                for v in (u + 1)..n {
                    let together = brute_lambda(&g, u, v) >= k && brute_lambda(&g, v, u) >= k;
                    assert_eq!(labels[u] == labels[v], together, "round {round} k {k}");
                }
            }
        }
    }
}

#[test]
fn condensation_orders_and_semiconnectedness() {
    let mut rng = Stream::new(205, 0);
    for round in 0..60 {
        let n = 1 + round % 9;
        let g = random_graph(&mut rng, n, 0.25, true);
        let adj = masks(&g);
        let reach_all: Vec<u32> = (0..n).map(|v| reach(&adj, v, 0, None)).collect();
        for order in [NeighborOrder::Sorted, NeighborOrder::Input] {
            let c = condensation(&g, order).unwrap();
            for u in 0..n {
                for v in 0..n {
                    let mutual = reach_all[u] >> v & 1 == 1 && reach_all[v] >> u & 1 == 1;
                    assert_eq!(c[u] == c[v], mutual);
                    // Reverse topological numbering: an arc between
                    // components goes to a lower label.
                    if adj[u] >> v & 1 == 1 && c[u] != c[v] {
                        assert!(c[u] > c[v]);
                    }
                }
            }
        }
        let semi = (0..n)
            .all(|u| (0..n).all(|v| reach_all[u] >> v & 1 == 1 || reach_all[v] >> u & 1 == 1));
        assert_eq!(is_semiconnected(&g).unwrap(), semi);
    }
}
