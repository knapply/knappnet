//! Randomized checks of cycle enumeration, cycle bases, girth, and
//! find_cycle against brute force on small multigraphs.

extern crate knappnet as engine;

use std::collections::BTreeSet;

use engine::metrics::cycles::{
    chordless_cycles, cycle_basis, find_cycle, girth, minimum_cycle_basis, simple_cycles,
    Orientation,
};
use engine::order::NeighborOrder;
use engine::rng::Stream;
use engine::Graph;

struct Net {
    g: Graph,
    s: Vec<usize>,
    t: Vec<usize>,
}

fn random_net(rng: &mut Stream, n: usize, p: f64, directed: bool) -> Net {
    let mut s = Vec::new();
    let mut t = Vec::new();
    for a in 0..n {
        for b in 0..n {
            if (!directed && b < a) || rng.uniform() >= p {
                continue;
            }
            let reps = if rng.uniform() < 0.12 {
                1 + rng.below(3) as usize
            } else {
                1
            };
            for _ in 0..reps {
                s.push(a);
                t.push(b);
            }
        }
    }
    let w: Vec<f64> = (0..s.len()).map(|_| (1 + rng.below(3)) as f64).collect();
    let si: Vec<i64> = s.iter().map(|&x| x as i64).collect();
    let ti: Vec<i64> = t.iter().map(|&x| x as i64).collect();
    let g = Graph::from_edges(n, &si, &ti, Some(&w), directed).unwrap();
    Net { g, s, t }
}

fn canon(mut c: Vec<u32>, directed: bool) -> Vec<u32> {
    let i = (0..c.len()).min_by_key(|&i| c[i]).unwrap();
    c.rotate_left(i);
    if !directed && c.len() > 2 && c[1] > c[c.len() - 1] {
        c[1..].reverse();
    }
    c
}

/// Node-sequence cycles of the multigraph by brute force (each from its
/// lowest node; undirected ones in one orientation).
fn brute_cycles(net: &Net, n: usize, bound: Option<usize>) -> BTreeSet<Vec<u32>> {
    let directed = net.g.is_directed();
    let mut succ = vec![BTreeSet::new(); n];
    let mut out = BTreeSet::new();
    let mut mult = std::collections::HashMap::new();
    for e in 0..net.s.len() {
        let (a, b) = (net.s[e], net.t[e]);
        if a == b {
            out.insert(vec![a as u32]);
            continue;
        }
        succ[a].insert(b);
        if !directed {
            succ[b].insert(a);
        }
        let key = if directed {
            (a, b)
        } else {
            (a.min(b), a.max(b))
        };
        *mult.entry(key).or_insert(0) += 1;
    }
    if bound == Some(1) {
        return out;
    }
    if !directed {
        for (&(a, b), &k) in &mult {
            if k > 1 {
                out.insert(vec![a as u32, b as u32]);
            }
        }
    }
    for s in 0..n {
        let mut stack = vec![vec![s]];
        while let Some(path) = stack.pop() {
            let v = *path.last().unwrap();
            for &w in &succ[v] {
                if w == s {
                    let ok = if directed {
                        path.len() >= 2
                    } else {
                        path.len() >= 3 && path[1] < path[path.len() - 1]
                    };
                    if ok {
                        out.insert(path.iter().map(|&x| x as u32).collect());
                    }
                } else if w > s && !path.contains(&w) && bound.map_or(true, |b| path.len() < b) {
                    let mut p = path.clone();
                    p.push(w);
                    stack.push(p);
                }
            }
        }
    }
    out
}

fn induced_edges(net: &Net, nodes: &[u32]) -> usize {
    (0..net.s.len())
        .filter(|&e| nodes.contains(&(net.s[e] as u32)) && nodes.contains(&(net.t[e] as u32)))
        .count()
}

#[test]
fn simple_and_chordless_cycles_match_brute_force() {
    let mut rng = Stream::new(401, 0);
    for round in 0..200 {
        let n = 1 + round % 7;
        let directed = round % 2 == 0;
        let net = random_net(&mut rng, n, 0.3 + 0.05 * (round % 5) as f64, directed);
        for bound in [None, Some(2), Some(3), Some(4)] {
            let want = brute_cycles(&net, n, bound);
            let got: BTreeSet<Vec<u32>> = simple_cycles(&net.g, bound, None)
                .unwrap()
                .into_iter()
                .map(|c| canon(c, directed))
                .collect();
            assert_eq!(got, want, "round {round} bound {bound:?}");
            let got_len = simple_cycles(&net.g, bound, None).unwrap().len();
            assert_eq!(got_len, want.len(), "duplicates in round {round}");
            let chordless_want: BTreeSet<Vec<u32>> = want
                .iter()
                .filter(|c| induced_edges(&net, c) == c.len())
                .cloned()
                .collect();
            let chordless_got: Vec<Vec<u32>> = chordless_cycles(&net.g, bound, None)
                .unwrap()
                .into_iter()
                .map(|c| canon(c, directed))
                .collect();
            let set: BTreeSet<Vec<u32>> = chordless_got.iter().cloned().collect();
            assert_eq!(
                set.len(),
                chordless_got.len(),
                "duplicates in round {round}"
            );
            assert_eq!(set, chordless_want, "round {round} bound {bound:?}");
        }
        // The girth is the shortest cycle.
        let all = brute_cycles(&net, n, None);
        let want = all.iter().map(|c| c.len()).min();
        let g = girth(&net.g);
        match want {
            Some(l) => assert_eq!(g, l as f64, "round {round}"),
            None => assert!(g.is_infinite()),
        }
        // find_cycle reports a closed walk of edges when a cycle exists.
        let has = !all.is_empty();
        let found = find_cycle(&net.g, None, Orientation::Original, NeighborOrder::Input).unwrap();
        assert_eq!(found.is_some(), has, "round {round}");
        if let Some(steps) = found {
            let ends: Vec<(usize, usize)> = steps
                .iter()
                .map(|st| {
                    let (a, b) = (net.s[st.edge], net.t[st.edge]);
                    if st.reversed {
                        (b, a)
                    } else {
                        (a, b)
                    }
                })
                .collect();
            if directed {
                for w in ends.windows(2) {
                    assert_eq!(w[0].1, w[1].0);
                }
                assert_eq!(ends.last().unwrap().1, ends[0].0);
            } else {
                let mut distinct: Vec<usize> = steps.iter().map(|s| s.edge).collect();
                distinct.sort_unstable();
                distinct.dedup();
                assert_eq!(distinct.len(), steps.len());
            }
        }
    }
}

/// GF(2) rank of edge-incidence vectors.
fn rank(vectors: &[u64]) -> usize {
    let mut basis: Vec<u64> = Vec::new();
    for &v in vectors {
        let mut x = v;
        for &b in &basis {
            x = x.min(x ^ b);
        }
        if x != 0 {
            basis.push(x);
            basis.sort_unstable_by(|a, b| b.cmp(a));
        }
    }
    basis.len()
}

/// Edge vector of a node-sequence cycle (parallel edges told apart by the
/// occurrence count in `taken`).
fn vector(
    net: &Net,
    c: &[u32],
    taken: &mut std::collections::HashMap<(usize, usize), usize>,
) -> u64 {
    let edges_between = |a: usize, b: usize| -> Vec<usize> {
        (0..net.s.len())
            .filter(|&e| (net.s[e].min(net.t[e]), net.s[e].max(net.t[e])) == (a.min(b), a.max(b)))
            .collect()
    };
    if c.len() <= 2 {
        let (a, b) = (c[0] as usize, *c.last().unwrap() as usize);
        let list = edges_between(a, b);
        let k = taken.entry((a.min(b), a.max(b))).or_insert(0);
        *k += 1;
        return if c.len() == 1 {
            1 << list[*k - 1]
        } else {
            (1 << list[0]) | (1 << list[*k])
        };
    }
    let mut v = 0u64;
    for i in 0..c.len() {
        let (a, b) = (c[i] as usize, c[(i + 1) % c.len()] as usize);
        v ^= 1 << edges_between(a, b)[0];
    }
    v
}

#[test]
fn bases_have_full_rank_and_minimum_weight() {
    let mut rng = Stream::new(402, 0);
    for round in 0..120 {
        let n = 1 + round % 7;
        let net = random_net(&mut rng, n, 0.45, false);
        let m = net.s.len();
        if m > 60 {
            continue;
        }
        let mut uf = engine::union_find::UnionFind::new(n);
        for e in 0..m {
            uf.union(net.s[e], net.t[e]);
        }
        let comps = (0..n).filter(|&v| uf.find(v) == v).count();
        let dim = m + comps - n;
        for order in [NeighborOrder::Sorted, NeighborOrder::Input] {
            let b = cycle_basis(&net.g, order, None).unwrap();
            assert_eq!(b.len(), dim);
            let mut taken = std::collections::HashMap::new();
            let vecs: Vec<u64> = b.iter().map(|c| vector(&net, c, &mut taken)).collect();
            assert_eq!(rank(&vecs), dim, "round {round}");
        }
        for weighted in [false, true] {
            let mb1 = minimum_cycle_basis(&net.g, weighted, 1).unwrap();
            let mb3 = minimum_cycle_basis(&net.g, weighted, 3).unwrap();
            assert_eq!(mb1, mb3);
            assert_eq!(mb1.cycles.len(), dim);
            // Greedy over every cycle (a matroid) gives the minimum
            // weights; parallel edges and loops get their own cycles.
            let w = net.g.edge_weights().unwrap();
            let weight_of = |e: usize| if weighted { w[e] } else { 1.0 };
            // Every cycle of the multigraph as an edge set: simple cycles
            // with every choice of parallel edges.
            let mut cands: Vec<(f64, u64)> = Vec::new();
            for c in brute_cycles(&net, n, None) {
                let mut sets: Vec<(f64, u64)> = vec![(0.0, 0)];
                let len = c.len();
                if len == 1 {
                    sets = (0..m)
                        .filter(|&e| net.s[e] == c[0] as usize && net.t[e] == c[0] as usize)
                        .map(|e| (weight_of(e), 1u64 << e))
                        .collect();
                } else if len == 2 {
                    let list: Vec<usize> = (0..m)
                        .filter(|&e| {
                            (net.s[e].min(net.t[e]), net.s[e].max(net.t[e]))
                                == (c[0] as usize, c[1] as usize)
                        })
                        .collect();
                    sets.clear();
                    for i in 0..list.len() {
                        for j in (i + 1)..list.len() {
                            sets.push((
                                weight_of(list[i]) + weight_of(list[j]),
                                (1 << list[i]) | (1 << list[j]),
                            ));
                        }
                    }
                } else {
                    for i in 0..len {
                        let (a, b) = (c[i] as usize, c[(i + 1) % len] as usize);
                        let list: Vec<usize> = (0..m)
                            .filter(|&e| {
                                (net.s[e].min(net.t[e]), net.s[e].max(net.t[e]))
                                    == (a.min(b), a.max(b))
                            })
                            .collect();
                        let mut next = Vec::new();
                        for &(sw, sb) in &sets {
                            for &e in &list {
                                next.push((sw + weight_of(e), sb | (1 << e)));
                            }
                        }
                        sets = next;
                    }
                }
                cands.extend(sets);
            }
            cands.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let mut chosen: Vec<u64> = Vec::new();
            let mut want = Vec::new();
            for (cw, cb) in cands {
                if chosen.len() == dim {
                    break;
                }
                let mut trial = chosen.clone();
                trial.push(cb);
                if rank(&trial) == trial.len() {
                    chosen = trial;
                    want.push(cw);
                }
            }
            let mut got = mb1.weights.clone();
            got.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(got, want, "round {round} weighted {weighted}");
        }
    }
}
