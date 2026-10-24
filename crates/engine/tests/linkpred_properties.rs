//! Randomized property tests for link prediction: brute-force scores on
//! small networks (dense ones use the bit-row kernels, sparse ones the
//! merged lists), agreement of the pair and non-edge entry points,
//! relabeling, and the thread count.

extern crate knappnet as engine;

use engine::metrics::linkpred::{
    link_scores, link_scores_non_edges, non_edge_count, LinkMethod, LinkOptions,
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

const METHODS: [LinkMethod; 9] = [
    LinkMethod::CommonNeighbors,
    LinkMethod::Jaccard,
    LinkMethod::AdamicAdar,
    LinkMethod::ResourceAllocation,
    LinkMethod::PreferentialAttachment,
    LinkMethod::CommonNeighborCentrality,
    LinkMethod::CnSoundarajanHopcroft,
    LinkMethod::RaSoundarajanHopcroft,
    LinkMethod::WithinInterCluster,
];

fn random_graph(rng: &mut Rng) -> (usize, Vec<i64>, Vec<i64>, bool) {
    let n = 1 + rng.below(40) as usize;
    let dense = rng.below(2) == 0;
    let m = if dense {
        rng.below((n * n) as u64 + 1) as usize
    } else {
        rng.below(2 * n as u64 + 1) as usize
    };
    let src = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    let tgt = (0..m).map(|_| rng.below(n as u64) as i64).collect();
    (n, src, tgt, rng.below(2) == 0)
}

fn brute(
    n: usize,
    src: &[i64],
    tgt: &[i64],
    comm: &[i64],
    method: LinkMethod,
    u: usize,
    v: usize,
) -> f64 {
    let mut a = vec![vec![false; n]; n];
    for (&s, &t) in src.iter().zip(tgt) {
        a[s as usize][t as usize] = true;
        a[t as usize][s as usize] = true;
    }
    let deg: Vec<f64> = (0..n)
        .map(|x| (a[x].iter().filter(|&&b| b).count() + a[x][x] as usize) as f64)
        .collect();
    let cn: Vec<usize> = (0..n)
        .filter(|&w| w != u && w != v && a[u][w] && a[v][w])
        .collect();
    // Hop distance by repeated expansion.
    let mut dist = vec![usize::MAX; n];
    dist[u] = 0;
    for step in 0..n {
        for x in 0..n {
            if dist[x] == step {
                for y in 0..n {
                    if a[x][y] && dist[y] == usize::MAX {
                        dist[y] = step + 1;
                    }
                }
            }
        }
    }
    let same: Vec<usize> = cn.iter().copied().filter(|&w| comm[w] == comm[u]).collect();
    match method {
        LinkMethod::CommonNeighbors => cn.len() as f64,
        LinkMethod::Jaccard => {
            let union = (0..n).filter(|&w| a[u][w] || a[v][w]).count();
            if union == 0 {
                0.0
            } else {
                cn.len() as f64 / union as f64
            }
        }
        LinkMethod::AdamicAdar => cn.iter().map(|&w| 1.0 / deg[w].ln()).sum(),
        LinkMethod::ResourceAllocation => cn.iter().map(|&w| 1.0 / deg[w]).sum(),
        LinkMethod::PreferentialAttachment => deg[u] * deg[v],
        LinkMethod::CommonNeighborCentrality => {
            let d = if dist[v] == usize::MAX {
                f64::INFINITY
            } else {
                dist[v] as f64
            };
            0.8 * cn.len() as f64 + 0.2 * n as f64 / d
        }
        LinkMethod::CnSoundarajanHopcroft => {
            (cn.len() + if comm[u] == comm[v] { same.len() } else { 0 }) as f64
        }
        LinkMethod::RaSoundarajanHopcroft if comm[u] == comm[v] => {
            same.iter().map(|&w| 1.0 / deg[w]).sum()
        }
        LinkMethod::WithinInterCluster if comm[u] == comm[v] => {
            same.len() as f64 / ((cn.len() - same.len()) as f64 + 0.001)
        }
        _ => 0.0,
    }
}

fn opts(method: LinkMethod, threads: usize) -> LinkOptions {
    LinkOptions {
        method,
        threads,
        ..Default::default()
    }
}

#[test]
fn scores_match_brute_force_and_entry_points_agree() {
    let mut rng = Rng(11);
    for _ in 0..60 {
        let (n, src, tgt, directed) = random_graph(&mut rng);
        let g = Graph::from_edges(n, &src, &tgt, None, directed).unwrap();
        let comm: Vec<i64> = (0..n).map(|_| rng.below(3) as i64).collect();
        let k = non_edge_count(&g);
        let mut us: Vec<u32> = Vec::new();
        let mut vs: Vec<u32> = Vec::new();
        for u in 0..n as u32 {
            for v in u + 1..n as u32 {
                us.push(u);
                vs.push(v);
            }
        }
        for &method in &METHODS {
            let all = link_scores(&g, &us, &vs, Some(&comm), opts(method, 1)).unwrap();
            for (i, (&u, &v)) in us.iter().zip(&vs).enumerate() {
                let b = brute(n, &src, &tgt, &comm, method, u as usize, v as usize);
                assert!(
                    (all[i] - b).abs() <= 1e-12 * b.abs().max(1.0),
                    "{method:?} ({u}, {v}): {} vs {b}",
                    all[i]
                );
            }
            let (mut a, mut b, mut s) = (vec![0; k], vec![0; k], vec![0.0; k]);
            link_scores_non_edges(&g, Some(&comm), opts(method, 3), &mut a, &mut b, &mut s)
                .unwrap();
            for i in 0..k {
                let j = us
                    .iter()
                    .zip(&vs)
                    .position(|(&x, &y)| x == a[i] && y == b[i])
                    .unwrap();
                assert_eq!(s[i].to_bits(), all[j].to_bits());
            }
        }
    }
}

#[test]
fn relabeling_permutes_and_threads_do_not_matter() {
    let mut rng = Rng(5);
    for _ in 0..30 {
        let (n, src, tgt, directed) = random_graph(&mut rng);
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            perm.swap(i, rng.below(i as u64 + 1) as usize);
        }
        let g = Graph::from_edges(n, &src, &tgt, None, directed).unwrap();
        let ps: Vec<i64> = src.iter().map(|&x| perm[x as usize] as i64).collect();
        let pt: Vec<i64> = tgt.iter().map(|&x| perm[x as usize] as i64).collect();
        let h = Graph::from_edges(n, &ps, &pt, None, directed).unwrap();
        let comm: Vec<i64> = (0..n).map(|v| (v % 3) as i64).collect();
        let mut pcomm = vec![0; n];
        for v in 0..n {
            pcomm[perm[v]] = comm[v];
        }
        let us: Vec<u32> = (0..n as u32)
            .flat_map(|u| (0..n as u32).map(move |_| u))
            .collect();
        let vs: Vec<u32> = (0..n as u32).flat_map(|_| 0..n as u32).collect();
        let keep: Vec<usize> = (0..us.len()).filter(|&i| us[i] != vs[i]).collect();
        let (us, vs): (Vec<u32>, Vec<u32>) = keep.iter().map(|&i| (us[i], vs[i])).unzip();
        let pu: Vec<u32> = us.iter().map(|&x| perm[x as usize] as u32).collect();
        let pv: Vec<u32> = vs.iter().map(|&x| perm[x as usize] as u32).collect();
        for &method in &METHODS {
            let a = link_scores(&g, &us, &vs, Some(&comm), opts(method, 1)).unwrap();
            let b = link_scores(&h, &pu, &pv, Some(&pcomm), opts(method, 1)).unwrap();
            let c = link_scores(&g, &us, &vs, Some(&comm), opts(method, 4)).unwrap();
            for i in 0..a.len() {
                assert!((a[i] - b[i]).abs() <= 1e-12 * a[i].abs().max(1.0));
                assert_eq!(a[i].to_bits(), c[i].to_bits());
            }
        }
    }
}

#[test]
fn invalid_arguments_are_errors() {
    let g = Graph::from_edges(3, &[0, 1], &[1, 2], None, false).unwrap();
    let cnc = opts(LinkMethod::CommonNeighborCentrality, 1);
    assert!(link_scores(&g, &[1], &[1], None, cnc).is_err());
    assert!(link_scores(&g, &[0], &[3], None, cnc).is_err());
    let sh = opts(LinkMethod::CnSoundarajanHopcroft, 1);
    assert!(link_scores(&g, &[0], &[2], None, sh).is_err());
    assert!(link_scores(&g, &[0], &[2], Some(&[0, 1]), sh).is_err());
    let mut wic = opts(LinkMethod::WithinInterCluster, 1);
    wic.delta = 0.0;
    assert!(link_scores(&g, &[0], &[2], Some(&[0, 0, 0]), wic).is_err());
}
