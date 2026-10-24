//! Randomized property tests for matchings: maximum-weight (with and
//! without maximum cardinality) and minimum-weight matchings against
//! exhaustive search on small graphs, and greedy maximal matchings.

extern crate knappnet as engine;

use engine::metrics::matching::{
    check_matching, max_weight_matching, maximal_matching, min_weight_matching,
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

/// Best (size, weight) over all matchings by recursion over the lowest
/// unmatched vertex; `key` orders the candidates.
fn best(
    n: usize,
    w: &[Vec<Option<f64>>],
    used: &mut Vec<bool>,
    key: &dyn Fn(usize, f64) -> (i64, f64),
) -> (usize, f64) {
    let Some(v) = (0..n).find(|&v| !used[v]) else {
        return (0, 0.0);
    };
    used[v] = true;
    let mut top = best(n, w, used, key);
    for u in v + 1..n {
        if let (false, Some(x)) = (used[u], w[v][u]) {
            used[u] = true;
            let (c, s) = best(n, w, used, key);
            used[u] = false;
            let cand = (c + 1, s + x);
            if key(cand.0, cand.1) > key(top.0, top.1) {
                top = cand;
            }
        }
    }
    used[v] = false;
    top
}

#[test]
fn optimal_matchings_match_exhaustive_search() {
    let mut rng = Rng(2024);
    for round in 0..400 {
        let n = 1 + rng.below(10) as usize;
        let m = rng.below((n * n / 2 + 1) as u64) as usize;
        let integer = round % 2 == 0;
        let mut s = Vec::new();
        let mut t = Vec::new();
        let mut wts = Vec::new();
        for _ in 0..m {
            s.push(rng.below(n as u64) as i64);
            t.push(rng.below(n as u64) as i64);
            let x = if integer {
                rng.below(20) as f64 - 4.0
            } else {
                (rng.below(1000) as f64) / 37.0 - 3.0
            };
            wts.push(x);
        }
        let g = Graph::from_edges(n, &s, &t, Some(&wts), false).unwrap();
        // Best and worst weight of each pair (parallel edges).
        let mut hi = vec![vec![None; n]; n];
        let mut lo = vec![vec![None; n]; n];
        for i in 0..m {
            let (a, b) = (s[i] as usize, t[i] as usize);
            if a == b {
                continue;
            }
            for (tab, better) in [(&mut hi, true), (&mut lo, false)] {
                for (x, y) in [(a, b), (b, a)] {
                    let cur: Option<f64> = tab[x][y];
                    if cur.map_or(true, |c| if better { wts[i] > c } else { wts[i] < c }) {
                        tab[x][y] = Some(wts[i]);
                    }
                }
            }
        }
        let total = |e: &[usize]| e.iter().map(|&i| wts[i]).sum::<f64>();
        let mut used = vec![false; n];
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
        let valid = |e: &[usize]| {
            let us: Vec<u32> = e.iter().map(|&i| s[i] as u32).collect();
            let vs: Vec<u32> = e.iter().map(|&i| t[i] as u32).collect();
            check_matching(&g, &us, &vs).unwrap().matching
        };
        let mw = max_weight_matching(&g, true, false).unwrap();
        assert!(valid(&mw));
        let (_, want) = best(n, &hi, &mut used, &|_, w| (0, w));
        assert!(
            close(total(&mw), want),
            "max weight {} vs {want}",
            total(&mw)
        );
        let mc = max_weight_matching(&g, true, true).unwrap();
        assert!(valid(&mc));
        let (c, w2) = best(n, &hi, &mut used, &|c, w| (c as i64, w));
        assert_eq!(mc.len(), c);
        assert!(close(total(&mc), w2), "max card {} vs {w2}", total(&mc));
        let mn = min_weight_matching(&g, true).unwrap();
        assert!(valid(&mn));
        let (c3, w3) = best(n, &lo, &mut used, &|c, w| (c as i64, -w));
        assert_eq!(mn.len(), c3);
        assert!(close(total(&mn), w3), "min weight {} vs {w3}", total(&mn));
        let mm = maximal_matching(&g);
        let us: Vec<u32> = mm.iter().map(|&i| s[i] as u32).collect();
        let vs: Vec<u32> = mm.iter().map(|&i| t[i] as u32).collect();
        let chk = check_matching(&g, &us, &vs).unwrap();
        assert!(chk.matching && chk.maximal);
    }
}
