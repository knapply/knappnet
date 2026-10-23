//! Randomized checks that the optimized dense kernels and the lazily
//! evaluated greedy modularity merger reproduce the plain textbook
//! computations exactly (bit for bit, or merge for merge).

use knappnet::linalg::{matmul, symmetric_eigen};
use knappnet::matrix::Dense;
use knappnet::metrics::partition::fast_greedy::fast_greedy;
use knappnet::rng::Stream;
use knappnet::Graph;

fn random_dense(rows: usize, cols: usize, r: &mut Stream, zeros: bool) -> Dense {
    let mut m = Dense::zeros(rows, cols);
    for i in 0..rows {
        for j in 0..cols {
            let x = if zeros && r.below(4) == 0 {
                0.0
            } else {
                r.uniform() * 2.0 - 1.0
            };
            m.set(i, j, x);
        }
    }
    m
}

/// The plain i-k-j product, as the blocked kernel must reproduce.
fn naive_matmul(a: &Dense, b: &Dense) -> Dense {
    let mut out = Dense::zeros(a.rows(), b.cols());
    for r in 0..a.rows() {
        for k in 0..a.cols() {
            let v = a.get(r, k);
            if v == 0.0 {
                continue;
            }
            for c in 0..b.cols() {
                out.set(r, c, out.get(r, c) + v * b.get(k, c));
            }
        }
    }
    out
}

#[test]
fn blocked_matmul_is_bitwise_the_plain_product() {
    let mut r = Stream::new(1, 0);
    for case in 0..40 {
        // Sizes straddle the block edges (64 inner, 512 columns).
        let dims = [1usize, 2, 7, 63, 64, 65, 130, 513, 600];
        let pick = |r: &mut Stream| dims[r.below(dims.len() as u64) as usize];
        let (p, q, s) = (pick(&mut r), pick(&mut r), pick(&mut r));
        if p * q * s > 30_000_000 {
            continue;
        }
        let a = random_dense(p, q, &mut r, case % 2 == 0);
        let b = random_dense(q, s, &mut r, false);
        let fast = matmul(&a, &b);
        let slow = naive_matmul(&a, &b);
        let fb: Vec<u64> = fast.as_slice().iter().map(|x| x.to_bits()).collect();
        let sb: Vec<u64> = slow.as_slice().iter().map(|x| x.to_bits()).collect();
        assert_eq!(fb, sb, "{p}x{q} times {q}x{s}");
    }
}

#[test]
fn symmetric_eigen_decomposes() {
    let mut r = Stream::new(2, 0);
    for n in [1usize, 2, 3, 5, 17, 64, 101] {
        let mut a = random_dense(n, n, &mut r, n > 5);
        for i in 0..n {
            for j in 0..i {
                a.set(i, j, a.get(j, i));
            }
        }
        let (vals, vecs) = symmetric_eigen(&a).unwrap();
        assert!(vals.windows(2).all(|w| w[0] <= w[1]));
        let scale = a.as_slice().iter().fold(1.0f64, |m, x| m.max(x.abs())) * n as f64;
        for (j, &val) in vals.iter().enumerate() {
            for i in 0..n {
                let av: f64 = (0..n).map(|k| a.get(i, k) * vecs.get(k, j)).sum();
                assert!((av - val * vecs.get(i, j)).abs() < 1e-11 * scale);
            }
            for l in 0..n {
                let dot: f64 = (0..n).map(|k| vecs.get(k, j) * vecs.get(k, l)).sum();
                let want = if j == l { 1.0 } else { 0.0 };
                assert!((dot - want).abs() < 1e-11 * n as f64);
            }
        }
    }
}

/// Reference greedy merger (Clauset, Newman, and Moore 2004) by brute
/// force: at every step evaluate the gain 2m W_IJ - gamma K_I K_J of every
/// connected pair of groups and merge the best, ties to the pair with the
/// smallest lowest-numbered members.
fn reference_merges(n: usize, edges: &[(usize, usize, f64)], gamma: f64) -> Vec<(u32, u32)> {
    let mut w = vec![vec![0.0f64; n]; n];
    let mut k = vec![0.0f64; n];
    for &(u, v, x) in edges {
        k[u] += x;
        k[v] += x;
        if u != v {
            w[u][v] += x;
            w[v][u] += x;
        }
    }
    let two_m: f64 = k.iter().sum();
    let mut alive = vec![true; n];
    let mut connected = vec![vec![false; n]; n];
    for &(u, v, _) in edges {
        if u != v {
            connected[u][v] = true;
            connected[v][u] = true;
        }
    }
    let mut merges = Vec::new();
    loop {
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..n {
            for b in a + 1..n {
                if !(alive[a] && alive[b] && connected[a][b]) {
                    continue;
                }
                let g = two_m * w[a][b] - gamma * k[a] * k[b];
                let better = match best {
                    None => true,
                    Some((bg, _, _)) => g.total_cmp(&bg).is_gt(),
                };
                if better {
                    best = Some((g, a, b));
                }
            }
        }
        let Some((_, a, b)) = best else { break };
        merges.push((a as u32, b as u32));
        for x in 0..n {
            if x != a && x != b {
                w[a][x] += w[b][x];
                w[x][a] = w[a][x];
                connected[a][x] |= connected[b][x];
                connected[x][a] = connected[a][x];
            }
        }
        k[a] += k[b];
        alive[b] = false;
    }
    merges
}

#[test]
fn lazy_greedy_merges_match_brute_force() {
    let mut r = Stream::new(3, 0);
    for case in 0..150 {
        let n = 2 + r.below(28) as usize;
        let p = 0.05 + r.uniform() * 0.5;
        let weighted = case % 3 == 1;
        let gamma = if case % 5 == 4 { 0.5 } else { 1.0 };
        let mut edges = Vec::new();
        for u in 0..n {
            for v in u..n {
                if (u != v || r.below(10) == 0) && r.uniform() < p {
                    let x = if weighted {
                        1.0 + r.below(4) as f64
                    } else {
                        1.0
                    };
                    edges.push((u, v, x));
                }
            }
        }
        let s: Vec<i64> = edges.iter().map(|e| e.0 as i64).collect();
        let t: Vec<i64> = edges.iter().map(|e| e.1 as i64).collect();
        let x: Vec<f64> = edges.iter().map(|e| e.2).collect();
        let g = Graph::from_edges(n, &s, &t, Some(&x), false).unwrap();
        let fast = fast_greedy(&g, gamma, weighted).unwrap();
        let unit: Vec<(usize, usize, f64)> = if weighted {
            edges.clone()
        } else {
            edges.iter().map(|&(u, v, _)| (u, v, 1.0)).collect()
        };
        assert_eq!(
            fast.merges,
            reference_merges(n, &unit, gamma),
            "case {case}"
        );
    }
}
