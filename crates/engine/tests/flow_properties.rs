//! Randomized checks of maximum flows, minimum cuts, minimum-cost flows,
//! and Gomory-Hu trees against brute force and optimality certificates.

extern crate knappnet as engine;

use engine::flow::{
    gomory_hu_tree, max_flow_min_cost, maximum_flow, min_cost_flow, minimum_cut,
    minimum_node_cut_capacitated, CutSide,
};
use engine::rng::Stream;
use engine::Graph;

struct Net {
    g: Graph,
    cap: Vec<f64>,
}

fn random_net(rng: &mut Stream, n: usize, p: f64, directed: bool, integer: bool) -> Net {
    let mut s = Vec::new();
    let mut t = Vec::new();
    let mut cap = Vec::new();
    for a in 0..n {
        for b in 0..n {
            if (!directed && b < a) || rng.uniform() >= p {
                continue;
            }
            // Occasional loops and parallel edges.
            let reps = if rng.uniform() < 0.1 { 2 } else { 1 };
            for _ in 0..reps {
                s.push(a as i64);
                t.push(b as i64);
                cap.push(if integer {
                    rng.below(6) as f64
                } else {
                    0.25 * (1 + rng.below(12)) as f64
                });
            }
        }
    }
    let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
    Net { g, cap }
}

/// Least cut capacity over all subsets containing every source and no
/// sink, and the union and intersection of the minimizing subsets.
fn brute_cuts(net: &Net, sources: &[usize], sinks: &[usize]) -> (f64, u32, u32) {
    let n = net.g.node_count();
    let (src, tgt) = (net.g.edge_sources(), net.g.edge_targets());
    let mut best = f64::INFINITY;
    let mut sets = Vec::new();
    for mask in 0u32..(1 << n) {
        if sources.iter().any(|&v| mask >> v & 1 == 0) || sinks.iter().any(|&v| mask >> v & 1 == 1)
        {
            continue;
        }
        let mut c = 0.0;
        for e in 0..src.len() {
            let (a, b) = (mask >> src[e] & 1 == 1, mask >> tgt[e] & 1 == 1);
            let crosses = if net.g.is_directed() { a && !b } else { a != b };
            if crosses && src[e] != tgt[e] {
                c += net.cap[e];
            }
        }
        if c < best - 1e-12 {
            best = c;
            sets.clear();
        }
        if (c - best).abs() <= 1e-12 {
            sets.push(mask);
        }
    }
    let union = sets.iter().fold(0, |a, &b| a | b);
    let inter = sets.iter().fold(u32::MAX, |a, &b| a & b);
    (best, union, inter)
}

fn mask_of(side: &[bool]) -> u32 {
    side.iter()
        .enumerate()
        .filter(|(_, &b)| b)
        .fold(0, |a, (i, _)| a | 1 << i)
}

#[test]
fn max_flow_equals_min_cut_and_flows_are_feasible() {
    let mut rng = Stream::new(101, 0);
    for round in 0..120 {
        let n = 2 + (round % 8);
        let directed = round % 2 == 0;
        let net = random_net(&mut rng, n, 0.35, directed, round % 3 != 0);
        let multi = n >= 4 && round % 4 == 1;
        let (sources, sinks): (Vec<usize>, Vec<usize>) = if multi {
            (vec![0, 1], vec![n - 2, n - 1])
        } else {
            (vec![0], vec![n - 1])
        };
        let (best, union, inter) = brute_cuts(&net, &sources, &sinks);
        let f = maximum_flow(&net.g, &sources, &sinks, Some(&net.cap)).unwrap();
        assert!((f.value - best).abs() < 1e-9, "round {round}");
        // Feasibility and conservation.
        let mut excess = vec![0.0f64; n];
        let (src, tgt) = (net.g.edge_sources(), net.g.edge_targets());
        for e in 0..src.len() {
            let x = f.edge_flow[e];
            if directed {
                assert!(x >= -1e-12 && x <= net.cap[e] + 1e-12);
            } else {
                assert!(x.abs() <= net.cap[e] + 1e-12);
            }
            excess[src[e] as usize] -= x;
            excess[tgt[e] as usize] += x;
        }
        for (v, x) in excess.iter().enumerate() {
            if !sources.contains(&v) && !sinks.contains(&v) {
                assert!(x.abs() < 1e-9, "conservation at {v}");
            }
        }
        let out: f64 = sources.iter().map(|&s| -excess[s]).sum();
        assert!((out - f.value).abs() < 1e-9);
        let small =
            minimum_cut(&net.g, &sources, &sinks, Some(&net.cap), CutSide::Smallest).unwrap();
        let large =
            minimum_cut(&net.g, &sources, &sinks, Some(&net.cap), CutSide::Largest).unwrap();
        assert_eq!(mask_of(&small.source_side), inter, "round {round}");
        assert_eq!(mask_of(&large.source_side), union, "round {round}");
    }
}

#[test]
fn relabeling_permutes_cuts() {
    let mut rng = Stream::new(102, 0);
    for round in 0..40 {
        let n = 3 + round % 7;
        let directed = round % 2 == 1;
        let net = random_net(&mut rng, n, 0.4, directed, true);
        let mut perm: Vec<usize> = (0..n).collect();
        rng.permutation(&mut perm);
        let s: Vec<i64> = net
            .g
            .edge_sources()
            .iter()
            .map(|&v| perm[v as usize] as i64)
            .collect();
        let t: Vec<i64> = net
            .g
            .edge_targets()
            .iter()
            .map(|&v| perm[v as usize] as i64)
            .collect();
        let h = Graph::from_edges(n, &s, &t, None, directed).unwrap();
        let a = minimum_cut(&net.g, &[0], &[n - 1], Some(&net.cap), CutSide::Largest).unwrap();
        let b = minimum_cut(
            &h,
            &[perm[0]],
            &[perm[n - 1]],
            Some(&net.cap),
            CutSide::Largest,
        )
        .unwrap();
        assert_eq!(a.value, b.value);
        for (v, &p) in perm.iter().enumerate() {
            assert_eq!(a.source_side[v], b.source_side[p]);
        }
    }
}

#[test]
fn node_cuts_match_brute_force() {
    let mut rng = Stream::new(103, 0);
    for round in 0..80 {
        let n = 3 + round % 7;
        let directed = round % 2 == 0;
        let net = random_net(&mut rng, n, 0.35, directed, true);
        let w: Vec<f64> = (0..n).map(|_| (1 + rng.below(4)) as f64).collect();
        let (s, t) = (0, n - 1);
        let (src, tgt) = (net.g.edge_sources(), net.g.edge_targets());
        let adjacent = (0..src.len()).any(|e| {
            (src[e] as usize == s && tgt[e] as usize == t)
                || (!directed && src[e] as usize == t && tgt[e] as usize == s)
        });
        let got =
            minimum_node_cut_capacitated(&net.g, &[s], &[t], Some(&w), None, CutSide::Smallest);
        if adjacent {
            assert!(got.is_err());
            continue;
        }
        let got = got.unwrap();
        // Brute force: least weight of a node set (not s or t) leaving no
        // path from s to t.
        let mut best = f64::INFINITY;
        for mask in 0u32..(1 << n) {
            if mask >> s & 1 == 1 || mask >> t & 1 == 1 {
                continue;
            }
            let cost: f64 = (0..n).filter(|&v| mask >> v & 1 == 1).map(|v| w[v]).sum();
            if cost >= best {
                continue;
            }
            let mut seen = vec![false; n];
            seen[s] = true;
            let mut stack = vec![s];
            while let Some(v) = stack.pop() {
                for e in 0..src.len() {
                    let (a, b) = (src[e] as usize, tgt[e] as usize);
                    let next = if a == v {
                        Some(b)
                    } else if !directed && b == v {
                        Some(a)
                    } else {
                        None
                    };
                    if let Some(x) = next {
                        if !seen[x] && mask >> x & 1 == 0 {
                            seen[x] = true;
                            stack.push(x);
                        }
                    }
                }
            }
            if !seen[t] {
                best = cost;
            }
        }
        assert_eq!(got.value, best, "round {round}");
        let listed: f64 = got.nodes.iter().map(|&v| w[v]).sum();
        assert_eq!(listed, best);
    }
}

/// Bellman-Ford over the residual network: whether some cycle of
/// negative cost has residual capacity (an optimal flow has none).
fn has_negative_residual_cycle(g: &Graph, cap: &[f64], cost: &[f64], flow: &[f64]) -> bool {
    let n = g.node_count();
    let mut arcs = Vec::new();
    for e in 0..g.edge_count() {
        let (a, b) = (g.edge_sources()[e] as usize, g.edge_targets()[e] as usize);
        if a == b {
            continue;
        }
        if flow[e] < cap[e] - 1e-9 {
            arcs.push((a, b, cost[e]));
        }
        if flow[e] > 1e-9 {
            arcs.push((b, a, -cost[e]));
        }
    }
    let mut dist = vec![0.0f64; n];
    for _ in 0..n {
        let mut changed = false;
        for &(a, b, c) in &arcs {
            if dist[a] + c < dist[b] - 1e-9 {
                dist[b] = dist[a] + c;
                changed = true;
            }
        }
        if !changed {
            return false;
        }
    }
    true
}

#[test]
fn min_cost_flows_are_optimal() {
    let mut rng = Stream::new(104, 0);
    for round in 0..150 {
        let n = 2 + round % 9;
        let net = random_net(&mut rng, n, 0.35, true, true);
        let m = net.g.edge_count();
        let cost: Vec<f64> = (0..m).map(|_| rng.below(9) as f64 - 2.0).collect();
        let (f, value) = max_flow_min_cost(&net.g, 0, n - 1, Some(&net.cap), Some(&cost)).unwrap();
        let mf = maximum_flow(&net.g, &[0], &[n - 1], Some(&net.cap)).unwrap();
        assert_eq!(value, mf.value);
        // Feasible: bounds and demands.
        let mut excess = vec![0.0f64; n];
        let mut total = 0.0;
        for (e, &x) in f.edge_flow.iter().enumerate() {
            assert!(x >= 0.0 && x <= net.cap[e], "round {round}");
            let (a, b) = (
                net.g.edge_sources()[e] as usize,
                net.g.edge_targets()[e] as usize,
            );
            excess[a] -= x;
            excess[b] += x;
            total += x * cost[e];
            if a == b {
                assert_eq!(x, if cost[e] < 0.0 { net.cap[e] } else { 0.0 });
            }
        }
        assert_eq!(total, f.cost);
        assert!(excess[1..n - 1].iter().all(|&x| x == 0.0));
        assert_eq!(excess[n - 1], value);
        // Optimal: no negative cycle in the residual network.
        assert!(
            !has_negative_residual_cycle(&net.g, &net.cap, &cost, &f.edge_flow),
            "round {round}"
        );
    }
}

#[test]
fn min_cost_flow_with_several_supplies() {
    let mut rng = Stream::new(105, 0);
    for round in 0..60 {
        let n = 4 + round % 6;
        let net = random_net(&mut rng, n, 0.5, true, true);
        let m = net.g.edge_count();
        let cost: Vec<f64> = (0..m).map(|_| rng.below(7) as f64).collect();
        // Supplies at 0 and 1, demands at the last two nodes, sized by what
        // the network can carry.
        let f1 = maximum_flow(&net.g, &[0, 1], &[n - 2, n - 1], Some(&net.cap))
            .unwrap()
            .value;
        let mut demand = vec![0.0; n];
        let half = (f1 / 2.0).floor();
        demand[0] = -half;
        demand[n - 1] = half;
        match min_cost_flow(&net.g, &demand, Some(&net.cap), Some(&cost)) {
            Ok(f) => assert!(!has_negative_residual_cycle(
                &net.g,
                &net.cap,
                &cost,
                &f.edge_flow
            )),
            Err(_) => {
                // Infeasible only if node 0 cannot send `half` to n - 1.
                let v = maximum_flow(&net.g, &[0], &[n - 1], Some(&net.cap))
                    .unwrap()
                    .value;
                assert!(v < half);
            }
        }
    }
}

#[test]
fn gomory_hu_tree_gives_every_minimum_cut() {
    let mut rng = Stream::new(106, 0);
    for round in 0..40 {
        let n = 1 + round % 9;
        let net = random_net(&mut rng, n, 0.4, false, round % 2 == 0);
        let tree = gomory_hu_tree(&net.g, Some(&net.cap)).unwrap();
        for u in 0..n {
            for v in (u + 1)..n {
                let (best, _, _) = brute_cuts(&net, &[u], &[v]);
                assert!((tree.min_cut_value(u, v) - best).abs() < 1e-9);
            }
        }
    }
}
