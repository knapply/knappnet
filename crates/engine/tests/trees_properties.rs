//! Randomized checks of spanning trees, branchings, arborescences, and DAG
//! functions against brute force on small networks.

extern crate knappnet as engine;

use engine::metrics::dag::{
    all_topological_sorts, dag_longest_path, dominance_frontiers, immediate_dominators,
    lexicographical_topological_sort, lowest_common_ancestors, topological_generations,
    topological_sort, transitive_closure, transitive_reduction, tree_lowest_common_ancestors,
    Reflexive,
};
use engine::metrics::trees::{optimum_branching, spanning_tree, BranchingKind, SpanningAlgorithm};
use engine::order::NeighborOrder;
use engine::rng::Stream;
use engine::Graph;

struct Net {
    g: Graph,
    s: Vec<usize>,
    t: Vec<usize>,
    w: Vec<f64>,
}

fn random_net(rng: &mut Stream, n: usize, p: f64, directed: bool, dag: bool) -> Net {
    let mut s = Vec::new();
    let mut t = Vec::new();
    let mut w = Vec::new();
    for a in 0..n {
        for b in 0..n {
            if (dag && b <= a) || (!directed && b < a) || rng.uniform() >= p {
                continue;
            }
            let reps = if rng.uniform() < 0.1 { 2 } else { 1 };
            for _ in 0..reps {
                s.push(a);
                t.push(b);
                // Small integer weights make ties common.
                w.push(rng.below(4) as f64 - 1.0);
            }
        }
    }
    let si: Vec<i64> = s.iter().map(|&x| x as i64).collect();
    let ti: Vec<i64> = t.iter().map(|&x| x as i64).collect();
    let g = Graph::from_edges(n, &si, &ti, Some(&w), directed).unwrap();
    Net { g, s, t, w }
}

/// Every set of edges that is a spanning forest (as many edges as n minus
/// the number of components), with the best total weight.
fn brute_forest(net: &Net, n: usize, maximum: bool) -> f64 {
    let m = net.s.len();
    let mut comps = engine::union_find::UnionFind::new(n);
    for e in 0..m {
        comps.union(net.s[e], net.t[e]);
    }
    let c = (0..n).filter(|&v| comps.find(v) == v).count();
    let need = n - c;
    let mut best: Option<f64> = None;
    for mask in 0u32..(1 << m) {
        if mask.count_ones() as usize != need {
            continue;
        }
        let mut uf = engine::union_find::UnionFind::new(n);
        let mut ok = true;
        let mut total = 0.0;
        for e in 0..m {
            if mask >> e & 1 == 1 {
                if uf.union(net.s[e], net.t[e]).is_none() {
                    ok = false;
                    break;
                }
                total += net.w[e];
            }
        }
        if ok {
            best = Some(match best {
                None => total,
                Some(b) if maximum => b.max(total),
                Some(b) => b.min(total),
            });
        }
    }
    best.unwrap_or(0.0)
}

#[test]
fn spanning_forests_are_optimal_and_agree() {
    let mut rng = Stream::new(301, 0);
    for round in 0..120 {
        let n = 1 + round % 7;
        let net = random_net(&mut rng, n, 0.5, false, false);
        if net.s.len() > 16 {
            continue;
        }
        for maximum in [false, true] {
            let k = spanning_tree(&net.g, true, maximum, SpanningAlgorithm::Kruskal).unwrap();
            for alg in [SpanningAlgorithm::Prim, SpanningAlgorithm::Boruvka] {
                assert_eq!(spanning_tree(&net.g, true, maximum, alg).unwrap(), k);
            }
            assert_eq!(k.weight, brute_forest(&net, n, maximum), "round {round}");
        }
    }
}

/// Best branching (or spanning arborescence) by trying every choice of at
/// most one incoming arc per node.
fn brute_branching(net: &Net, n: usize, kind: BranchingKind) -> Option<f64> {
    let mut choices: Vec<Vec<Option<usize>>> = vec![vec![None]; n];
    for e in 0..net.s.len() {
        if net.s[e] != net.t[e] {
            choices[net.t[e]].push(Some(e));
        }
    }
    let arbo = matches!(
        kind,
        BranchingKind::MaximumArborescence | BranchingKind::MinimumArborescence
    );
    let minimize = matches!(
        kind,
        BranchingKind::MinimumBranching | BranchingKind::MinimumArborescence
    );
    let mut best: Option<f64> = None;
    let mut pick = vec![0usize; n];
    loop {
        let parent: Vec<Option<usize>> = (0..n)
            .map(|v| choices[v][pick[v]].map(|e| net.s[e]))
            .collect();
        let roots = parent.iter().filter(|p| p.is_none()).count();
        let acyclic = (0..n).all(|v| {
            let mut x = Some(v);
            let mut steps = 0;
            while let Some(y) = x {
                steps += 1;
                if steps > n {
                    return false;
                }
                x = parent[y];
            }
            true
        });
        if acyclic && (!arbo || roots == 1) {
            let total: f64 = (0..n)
                .filter_map(|v| choices[v][pick[v]])
                .map(|e| net.w[e])
                .sum();
            best = Some(match best {
                None => total,
                Some(b) if minimize => b.min(total),
                Some(b) => b.max(total),
            });
        }
        // Next choice.
        let mut i = 0;
        loop {
            if i == n {
                return best;
            }
            pick[i] += 1;
            if pick[i] < choices[i].len() {
                break;
            }
            pick[i] = 0;
            i += 1;
        }
    }
}

#[test]
fn branchings_are_optimal() {
    let mut rng = Stream::new(302, 0);
    for round in 0..300 {
        let n = 1 + round % 6;
        let net = random_net(&mut rng, n, 0.4, true, false);
        for kind in [
            BranchingKind::MaximumBranching,
            BranchingKind::MinimumBranching,
            BranchingKind::MaximumArborescence,
            BranchingKind::MinimumArborescence,
        ] {
            let got = optimum_branching(&net.g, true, kind, None);
            let want = brute_branching(&net, n, kind);
            match (got, want) {
                (Ok(b), Some(w)) => {
                    assert_eq!(b.weight, w, "round {round} {kind:?}");
                    // A branching: one arc into a node at most, no cycle.
                    let mut into = vec![usize::MAX; n];
                    for &e in &b.edges {
                        assert_eq!(into[net.t[e]], usize::MAX);
                        into[net.t[e]] = net.s[e];
                    }
                }
                (Err(_), None) => {}
                (got, want) => panic!("round {round} {kind:?}: {got:?} vs {want:?}"),
            }
        }
    }
}

fn reach(n: usize, s: &[usize], t: &[usize]) -> Vec<Vec<bool>> {
    let mut r = vec![vec![false; n]; n];
    for (u, row) in r.iter_mut().enumerate() {
        let mut stack = vec![u];
        while let Some(v) = stack.pop() {
            for e in 0..s.len() {
                if s[e] == v && !row[t[e]] {
                    row[t[e]] = true;
                    stack.push(t[e]);
                }
            }
        }
    }
    r
}

#[test]
fn dag_functions_match_brute_force() {
    let mut rng = Stream::new(303, 0);
    for round in 0..80 {
        let n = 1 + round % 7;
        let net = random_net(&mut rng, n, 0.4, true, true);
        let r = reach(n, &net.s, &net.t);
        for order in [NeighborOrder::Sorted, NeighborOrder::Input] {
            let topo = topological_sort(&net.g, order).unwrap();
            let mut pos = vec![0; n];
            for (i, &v) in topo.iter().enumerate() {
                pos[v as usize] = i;
            }
            for e in 0..net.s.len() {
                assert!(pos[net.s[e]] < pos[net.t[e]]);
            }
            let gens = topological_generations(&net.g, order).unwrap();
            assert_eq!(gens.iter().map(|g| g.len()).sum::<usize>(), n);
            // Every linear extension exactly once.
            let all = all_topological_sorts(&net.g, order, None).unwrap();
            let mut count = 0usize;
            let mut perm: Vec<usize> = (0..n).collect();
            permutations(&mut perm, 0, &mut |p| {
                let mut at = vec![0; n];
                for (i, &v) in p.iter().enumerate() {
                    at[v] = i;
                }
                if (0..net.s.len()).all(|e| at[net.s[e]] < at[net.t[e]]) {
                    count += 1;
                }
            });
            let mut distinct = all.clone();
            distinct.sort();
            distinct.dedup();
            assert_eq!(distinct.len(), all.len());
            assert_eq!(all.len(), count, "round {round}");
        }
        let lex = lexicographical_topological_sort(&net.g).unwrap();
        assert_eq!(lex.len(), n);
        let closure = transitive_closure(&net.g, Reflexive::Cycles).unwrap();
        let want: Vec<(u32, u32)> = (0..n)
            .flat_map(|u| (0..n).map(move |v| (u, v)))
            .filter(|&(u, v)| r[u][v])
            .map(|(u, v)| (u as u32, v as u32))
            .collect();
        assert_eq!(closure, want);
        let red = transitive_reduction(&net.g).unwrap();
        for &(u, v) in &red {
            let (u, v) = (u as usize, v as usize);
            assert!((0..net.s.len()).any(|e| net.s[e] == u && net.t[e] == v));
            assert!(!(0..net.s.len()).any(|e| net.s[e] == u && net.t[e] != v && r[net.t[e]][v]));
        }
        // Reduction preserves reachability.
        let (rs, rt): (Vec<usize>, Vec<usize>) =
            red.iter().map(|&(a, b)| (a as usize, b as usize)).unzip();
        assert_eq!(reach(n, &rs, &rt), r);
        // Longest path against enumeration (weights as lengths, never
        // extending a path below zero, as NetworkX).
        let (path, len) = dag_longest_path(&net.g, true, NeighborOrder::Input).unwrap();
        let mut best = 0.0f64;
        for v in 0..n {
            longest_from(v, 0.0, &net, &mut best);
        }
        assert_eq!(len, best, "round {round}");
        let mut total = 0.0;
        for p in path.windows(2) {
            total += (0..net.s.len())
                .filter(|&e| net.s[e] == p[0] as usize && net.t[e] == p[1] as usize)
                .map(|e| net.w[e])
                .fold(f64::NEG_INFINITY, f64::max);
        }
        assert_eq!(total, len);
        // Lowest common ancestors are common ancestors without a child that
        // is one.
        let pairs: Vec<(usize, usize)> = (0..n).flat_map(|u| (u..n).map(move |v| (u, v))).collect();
        let lca = lowest_common_ancestors(&net.g, &pairs, NeighborOrder::Input).unwrap();
        for (&(u, v), x) in pairs.iter().zip(lca) {
            let anc = |a: usize, b: usize| a == b || r[a][b];
            let common: Vec<usize> = (0..n).filter(|&c| anc(c, u) && anc(c, v)).collect();
            match x {
                None => assert!(common.is_empty()),
                Some(x) => {
                    let x = x as usize;
                    assert!(common.contains(&x));
                    assert!(!(0..net.s.len()).any(|e| net.s[e] == x && common.contains(&net.t[e])));
                }
            }
        }
    }
}

fn longest_from(v: usize, total: f64, net: &Net, best: &mut f64) {
    *best = best.max(total);
    for e in 0..net.s.len() {
        if net.s[e] == v && total + net.w[e] >= 0.0 {
            longest_from(net.t[e], total + net.w[e], net, best);
        }
    }
}

fn permutations(p: &mut Vec<usize>, k: usize, f: &mut impl FnMut(&[usize])) {
    if k == p.len() {
        f(p);
        return;
    }
    for i in k..p.len() {
        p.swap(k, i);
        permutations(p, k + 1, f);
        p.swap(k, i);
    }
}

#[test]
fn dominators_match_the_definition() {
    let mut rng = Stream::new(304, 0);
    for round in 0..80 {
        let n = 1 + round % 8;
        let net = random_net(&mut rng, n, 0.3, true, false);
        let idom = immediate_dominators(&net.g, 0).unwrap();
        let df = dominance_frontiers(&net.g, 0).unwrap();
        // d dominates v: v unreachable from 0 once d is removed.
        let reachable = |removed: Option<usize>| -> Vec<bool> {
            let mut seen = vec![false; n];
            if removed == Some(0) {
                return seen;
            }
            seen[0] = true;
            let mut stack = vec![0];
            while let Some(v) = stack.pop() {
                for e in 0..net.s.len() {
                    let w = net.t[e];
                    if net.s[e] == v && !seen[w] && Some(w) != removed {
                        seen[w] = true;
                        stack.push(w);
                    }
                }
            }
            seen
        };
        let base = reachable(None);
        let dom = |d: usize, v: usize| d == v || d == 0 || !reachable(Some(d))[v];
        for v in 0..n {
            if !base[v] {
                assert_eq!(idom[v], None);
                assert_eq!(df[v], None);
                continue;
            }
            if v == 0 {
                assert_eq!(idom[v], Some(0));
            } else {
                let strict: Vec<usize> =
                    (0..n).filter(|&d| d != v && base[d] && dom(d, v)).collect();
                let want = strict
                    .iter()
                    .copied()
                    .find(|&d| strict.iter().all(|&o| dom(o, d)))
                    .unwrap();
                assert_eq!(idom[v], Some(want as u32), "round {round}");
            }
            let want: Vec<u32> = (0..n)
                .filter(|&y| base[y])
                .filter(|&y| {
                    let preds_dom = (0..net.s.len())
                        .any(|e| net.t[e] == y && base[net.s[e]] && dom(v, net.s[e]));
                    preds_dom && !(dom(v, y) && v != y)
                })
                .map(|y| y as u32)
                .collect();
            assert_eq!(df[v], Some(want), "round {round} node {v}");
        }
    }
}

#[test]
fn tree_ancestors() {
    let mut rng = Stream::new(305, 0);
    for round in 0..40 {
        let n = 1 + round % 12;
        // A random rooted tree: each node v > 0 hangs from a node below v.
        let parent: Vec<usize> = (0..n)
            .map(|v| {
                if v == 0 {
                    0
                } else {
                    rng.below(v as u64) as usize
                }
            })
            .collect();
        let s: Vec<i64> = (1..n).map(|v| parent[v] as i64).collect();
        let t: Vec<i64> = (1..n).map(|v| v as i64).collect();
        let g = Graph::from_edges(n, &s, &t, None, true).unwrap();
        let pairs: Vec<(usize, usize)> = (0..n).flat_map(|u| (0..n).map(move |v| (u, v))).collect();
        let got = tree_lowest_common_ancestors(&g, None, &pairs).unwrap();
        let chain = |mut v: usize| {
            let mut c = vec![v];
            while v != 0 {
                v = parent[v];
                c.push(v);
            }
            c
        };
        for (&(u, v), x) in pairs.iter().zip(got) {
            let cv = chain(v);
            let want = chain(u).into_iter().find(|a| cv.contains(a)).unwrap();
            assert_eq!(x as usize, want);
        }
    }
}
