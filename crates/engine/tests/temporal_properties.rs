//! Randomized property tests for temporal networks: every path criterion
//! against a brute-force search over walks, time-reversal symmetry,
//! window monotonicity, the static limit, slicing, and event statistics.

extern crate knappnet as engine;

use engine::temporal::{
    temporal_closeness, temporal_paths, temporal_reach_counts, Criterion, Direction, PathOptions,
    TemporalNetwork, Traversal,
};

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

#[derive(Debug, Clone)]
struct Case {
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    start: Vec<f64>,
    end: Vec<f64>,
    directed: bool,
}

impl Case {
    fn net(&self) -> TemporalNetwork {
        TemporalNetwork::new(
            self.n,
            &self.src,
            &self.tgt,
            &self.start,
            Some(&self.end),
            self.directed,
        )
        .unwrap()
    }
}

fn random_case(rng: &mut Rng) -> Case {
    let n = 1 + rng.below(7) as usize;
    let m = rng.below(16) as usize;
    let half = rng.below(3) == 0;
    let mut c = Case {
        n,
        src: Vec::new(),
        tgt: Vec::new(),
        start: Vec::new(),
        end: Vec::new(),
        directed: rng.below(2) == 0,
    };
    for _ in 0..m {
        c.src.push(rng.below(n as u64) as i64);
        c.tgt.push(rng.below(n as u64) as i64);
        let mut s = rng.below(7) as f64;
        if half {
            s /= 2.0;
        }
        let d = if rng.below(2) == 0 {
            0.0
        } else {
            rng.below(4) as f64
        };
        c.start.push(s);
        c.end.push(s + d);
    }
    c
}

fn random_options(rng: &mut Rng) -> PathOptions {
    let traversal = match rng.below(4) {
        0 | 1 => Traversal::Constant(0.0),
        2 => Traversal::Constant(1.0),
        _ => Traversal::Duration,
    };
    let min_wait = if rng.below(3) == 0 { 1.0 } else { 0.0 };
    let max_wait = match rng.below(4) {
        0 => min_wait,
        1 => min_wait + 1.0,
        2 => min_wait + 2.5,
        _ => f64::INFINITY,
    };
    let (start, end) = if rng.below(2) == 0 {
        (None, None)
    } else {
        let a = rng.below(4) as f64;
        (Some(a), Some(a + rng.below(6) as f64))
    };
    PathOptions {
        start,
        end,
        min_wait,
        max_wait,
        strict: rng.below(3) == 0,
        traversal,
        threads: 1 + rng.below(3) as usize,
    }
}

/// Temporal arcs built directly from the contacts.
struct Arc {
    from: usize,
    to: usize,
    dep: f64,
    arr: f64,
}

/// Brute force over walks: the graph whose nodes are arcs, with an edge
/// from a to b when b can follow a. Shares no code with the engine.
struct Brute {
    arcs: Vec<Arc>,
    next: Vec<Vec<usize>>,
    lo: f64,
    hi: f64,
}

impl Brute {
    fn new(c: &Case, o: &PathOptions) -> Brute {
        let mut arcs = Vec::new();
        for e in 0..c.src.len() {
            let (u, v) = (c.src[e] as usize, c.tgt[e] as usize);
            if u == v {
                continue;
            }
            let lambda = match o.traversal {
                Traversal::Constant(x) => x,
                Traversal::Duration => c.end[e] - c.start[e],
            };
            let mut push = |a: usize, b: usize| {
                arcs.push(Arc {
                    from: a,
                    to: b,
                    dep: c.start[e],
                    arr: c.start[e] + lambda,
                })
            };
            push(u, v);
            if !c.directed {
                push(v, u);
            }
        }
        let next = arcs
            .iter()
            .map(|a| {
                (0..arcs.len())
                    .filter(|&j| {
                        let b = &arcs[j];
                        let w = b.dep - a.arr;
                        b.from == a.to
                            && (if o.strict {
                                w > o.min_wait
                            } else {
                                w >= o.min_wait
                            })
                            && w <= o.max_wait
                    })
                    .collect()
            })
            .collect();
        let lo = o.start.unwrap_or_else(|| {
            c.start
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min)
                .min(if c.start.is_empty() {
                    0.0
                } else {
                    f64::INFINITY
                })
        });
        let last = c
            .start
            .iter()
            .zip(&c.end)
            .map(|(&s, &e)| match o.traversal {
                Traversal::Constant(x) => s + x,
                Traversal::Duration => e,
            })
            .fold(f64::NEG_INFINITY, f64::max);
        let hi = o.end.unwrap_or(if c.start.is_empty() { lo } else { last });
        Brute { arcs, next, lo, hi }
    }

    /// Breadth-first search over arcs from `starts`: hop count of every
    /// reachable arc (1 for the starts).
    fn bfs(&self, starts: &[usize]) -> Vec<Option<usize>> {
        let mut depth = vec![None; self.arcs.len()];
        let mut queue = std::collections::VecDeque::new();
        for &s in starts {
            if depth[s].is_none() {
                depth[s] = Some(1);
                queue.push_back(s);
            }
        }
        while let Some(a) = queue.pop_front() {
            let d = depth[a].unwrap();
            for &b in &self.next[a] {
                if depth[b].is_none() {
                    depth[b] = Some(d + 1);
                    queue.push_back(b);
                }
            }
        }
        depth
    }

    fn starts(&self, from: usize) -> Vec<usize> {
        (0..self.arcs.len())
            .filter(|&a| self.arcs[a].from == from && self.arcs[a].dep >= self.lo)
            .collect()
    }

    /// Value for walks from i to j under a criterion.
    fn value(&self, c: Criterion, i: usize, j: usize) -> f64 {
        if i == j {
            return match c {
                Criterion::EarliestArrival => self.lo,
                Criterion::LatestDeparture => self.hi,
                _ => 0.0,
            };
        }
        let ends = |b: usize| self.arcs[b].to == j && self.arcs[b].arr <= self.hi;
        match c {
            Criterion::EarliestArrival | Criterion::Shortest => {
                let depth = self.bfs(&self.starts(i));
                let mut best = f64::INFINITY;
                for (b, d) in depth.iter().enumerate() {
                    if let Some(d) = *d {
                        if ends(b) {
                            let v = if c == Criterion::Shortest {
                                d as f64
                            } else {
                                self.arcs[b].arr
                            };
                            best = best.min(v);
                        }
                    }
                }
                best
            }
            Criterion::Fastest | Criterion::LatestDeparture => {
                let mut best = if c == Criterion::Fastest {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                };
                for s in self.starts(i) {
                    let depth = self.bfs(&[s]);
                    for (b, d) in depth.iter().enumerate() {
                        if d.is_some() && ends(b) {
                            if c == Criterion::Fastest {
                                best = best.min(self.arcs[b].arr - self.arcs[s].dep);
                            } else {
                                best = best.max(self.arcs[s].dep);
                            }
                        }
                    }
                }
                best
            }
        }
    }

    /// Pan and Saramäki's average temporal distance, integrated exactly
    /// over the pieces between consecutive walk departures.
    fn tau(&self, i: usize, j: usize) -> f64 {
        let mut pairs = Vec::new();
        for s in self.starts(i) {
            let depth = self.bfs(&[s]);
            for (b, d) in depth.iter().enumerate() {
                if d.is_some() && self.arcs[b].to == j && self.arcs[b].arr <= self.hi {
                    pairs.push((self.arcs[s].dep, self.arcs[b].arr));
                }
            }
        }
        if pairs.is_empty() {
            return f64::INFINITY;
        }
        let first_arrival = pairs.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let span = self.hi - self.lo;
        let mut cuts: Vec<f64> = pairs.iter().map(|p| p.0).collect();
        cuts.push(self.lo);
        cuts.push(self.hi);
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let mut total = 0.0;
        for w in cuts.windows(2) {
            let (a, b) = (w[0], w[1]);
            // On (a, b] the earliest arrival comes from walks leaving at
            // or after b, or from the first walk of the next period.
            let arrive = pairs
                .iter()
                .filter(|p| p.0 >= b)
                .map(|p| p.1)
                .fold(first_arrival + span, f64::min);
            total += (b - a) * (arrive - (a + b) / 2.0);
        }
        total / span
    }
}

fn same(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

#[test]
fn every_criterion_matches_brute_force() {
    let mut rng = Rng(20261023);
    let crits = [
        Criterion::EarliestArrival,
        Criterion::LatestDeparture,
        Criterion::Fastest,
        Criterion::Shortest,
    ];
    for _ in 0..600 {
        let c = random_case(&mut rng);
        let o = random_options(&mut rng);
        let net = c.net();
        let brute = Brute::new(&c, &o);
        let roots: Vec<usize> = (0..c.n).collect();
        for &crit in &crits {
            for dir in [Direction::Out, Direction::In] {
                let got = temporal_paths(&net, crit, dir, &roots, &o).unwrap();
                for r in 0..c.n {
                    for j in 0..c.n {
                        let want = match dir {
                            Direction::Out => brute.value(crit, r, j),
                            Direction::In => brute.value(crit, j, r),
                        };
                        let g = got[r * c.n + j];
                        assert!(
                            same(g, want),
                            "{crit:?} {dir:?} root {r} node {j}: got {g}, want {want}\n{c:?}\n{o:?}"
                        );
                    }
                }
            }
        }
        for dir in [Direction::Out, Direction::In] {
            let counts = temporal_reach_counts(&net, dir, &o).unwrap();
            for (r, &count) in counts.iter().enumerate() {
                let want = (0..c.n)
                    .filter(|&j| {
                        j != r
                            && match dir {
                                Direction::Out => brute.value(Criterion::Shortest, r, j),
                                Direction::In => brute.value(Criterion::Shortest, j, r),
                            }
                            .is_finite()
                    })
                    .count();
                assert_eq!(count, want as f64, "{dir:?} reach of {r}\n{c:?}\n{o:?}");
            }
        }
        if brute.hi > brute.lo && c.n > 1 {
            let tc = temporal_closeness(&net, &o, true).unwrap();
            for i in 0..c.n {
                let mut sum = 0.0;
                for j in 0..c.n {
                    let want = if i == j { 0.0 } else { brute.tau(i, j) };
                    let g = tc.distances[i * c.n + j];
                    assert!(
                        same(g, want),
                        "tau {i} {j}: got {g}, want {want}\n{c:?}\n{o:?}"
                    );
                    if i != j {
                        sum += 1.0 / want;
                    }
                }
                assert!(same(tc.closeness[i], sum / (c.n - 1) as f64));
            }
        }
    }
}

/// Reversing time maps walks from i to j onto walks from j to i, so
/// earliest arrival becomes minus the latest departure, and fastest and
/// shortest walks are unchanged.
#[test]
fn time_reversal_symmetry() {
    let mut rng = Rng(7);
    for _ in 0..300 {
        let c = random_case(&mut rng);
        let mut o = random_options(&mut rng);
        if o.traversal == Traversal::Constant(1.0) {
            o.traversal = Traversal::Constant(0.0);
        }
        let mut r = c.clone();
        std::mem::swap(&mut r.src, &mut r.tgt);
        if o.traversal == Traversal::Duration {
            r.start = c.end.iter().map(|e| -e).collect();
            r.end = c.start.iter().map(|s| -s).collect();
        } else {
            // A zero traversal time ignores durations: reverse start times.
            r.start = c.start.iter().map(|s| -s).collect();
            r.end = r.start.clone();
        }
        let mut ro = o;
        let brute = Brute::new(&c, &o);
        ro.start = Some(-brute.hi);
        ro.end = Some(-brute.lo);
        let (a, b) = (c.net(), r.net());
        let roots: Vec<usize> = (0..c.n).collect();
        let ea =
            temporal_paths(&a, Criterion::EarliestArrival, Direction::Out, &roots, &o).unwrap();
        let ld =
            temporal_paths(&b, Criterion::LatestDeparture, Direction::In, &roots, &ro).unwrap();
        for k in 0..ea.len() {
            assert!(same(ea[k], -ld[k]), "{c:?} {o:?}");
        }
        for crit in [Criterion::Fastest, Criterion::Shortest] {
            let x = temporal_paths(&a, crit, Direction::Out, &roots, &o).unwrap();
            let y = temporal_paths(&b, crit, Direction::In, &roots, &ro).unwrap();
            for k in 0..x.len() {
                assert!(same(x[k], y[k]), "{crit:?} {c:?} {o:?}");
            }
        }
    }
}

/// A larger window never loses a walk.
#[test]
fn window_monotonicity() {
    let mut rng = Rng(99);
    for _ in 0..300 {
        let c = random_case(&mut rng);
        let net = c.net();
        let a = rng.below(4) as f64;
        let b = a + rng.below(5) as f64;
        let inner = PathOptions {
            start: Some(a),
            end: Some(b),
            ..Default::default()
        };
        let outer = PathOptions {
            start: Some(a - 1.0),
            end: Some(b + 1.0),
            ..inner
        };
        for dir in [Direction::Out, Direction::In] {
            let x = temporal_reach_counts(&net, dir, &inner).unwrap();
            let y = temporal_reach_counts(&net, dir, &outer).unwrap();
            assert!(x.iter().zip(&y).all(|(p, q)| p <= q));
        }
        let roots: Vec<usize> = (0..c.n).collect();
        let x = temporal_paths(&net, Criterion::Shortest, Direction::Out, &roots, &inner).unwrap();
        let y = temporal_paths(&net, Criterion::Shortest, Direction::Out, &roots, &outer).unwrap();
        assert!(x.iter().zip(&y).all(|(p, q)| q <= p));
    }
}

/// With every contact at the same time and zero waiting allowed, temporal
/// reachability is static reachability and the fewest contacts is the
/// static shortest-path length.
#[test]
fn equal_times_give_static_paths() {
    let mut rng = Rng(5);
    for _ in 0..200 {
        let mut c = random_case(&mut rng);
        let t = rng.below(5) as f64;
        c.start.iter_mut().for_each(|s| *s = t);
        c.end = c.start.clone();
        let net = c.net();
        let n = c.n;
        let roots: Vec<usize> = (0..n).collect();
        let sp = temporal_paths(
            &net,
            Criterion::Shortest,
            Direction::Out,
            &roots,
            &PathOptions::default(),
        )
        .unwrap();
        // Static BFS.
        let mut adj = vec![Vec::new(); n];
        for e in 0..c.src.len() {
            let (u, v) = (c.src[e] as usize, c.tgt[e] as usize);
            adj[u].push(v);
            if !c.directed {
                adj[v].push(u);
            }
        }
        for r in 0..n {
            let mut d = vec![f64::INFINITY; n];
            d[r] = 0.0;
            let mut q = std::collections::VecDeque::from([r]);
            while let Some(u) = q.pop_front() {
                for &v in &adj[u] {
                    if d[v].is_infinite() {
                        d[v] = d[u] + 1.0;
                        q.push_back(v);
                    }
                }
            }
            assert_eq!(&sp[r * n..(r + 1) * n], &d[..]);
        }
        // Strictly increasing times allow single contacts only.
        let strict = PathOptions {
            strict: true,
            ..Default::default()
        };
        let sp =
            temporal_paths(&net, Criterion::Shortest, Direction::Out, &roots, &strict).unwrap();
        assert!(sp.iter().all(|&x| x <= 1.0 || x.is_infinite()));
    }
}

#[test]
fn results_do_not_depend_on_threads() {
    let mut rng = Rng(11);
    for _ in 0..50 {
        let c = random_case(&mut rng);
        let net = c.net();
        let roots: Vec<usize> = (0..c.n).collect();
        let mut o = random_options(&mut rng);
        o.threads = 1;
        let a = temporal_paths(&net, Criterion::Fastest, Direction::Out, &roots, &o).unwrap();
        o.threads = 3;
        let b = temporal_paths(&net, Criterion::Fastest, Direction::Out, &roots, &o).unwrap();
        assert_eq!(a, b);
    }
}

/// The one-pass window sweep returns exactly the contacts active in each
/// window, and aggregation conserves counts and durations.
#[test]
fn windows_and_aggregation_match_definitions() {
    use engine::temporal::{aggregate, contacts_in_windows, regular_windows};
    let mut rng = Rng(3);
    for _ in 0..300 {
        let c = random_case(&mut rng);
        let net = c.net();
        let width = 0.5 * (1 + rng.below(6)) as f64;
        let step = 0.5 * (1 + rng.below(6)) as f64;
        let (a, b) = regular_windows(-1.0, 7.0, width, step).unwrap();
        let w = contacts_in_windows(&net, &a, &b).unwrap();
        for k in 0..a.len() {
            let want: Vec<usize> = (0..c.src.len())
                .filter(|&e| {
                    let (s, t) = (c.start[e], c.end[e]);
                    if s == t {
                        a[k] <= s && s < b[k]
                    } else {
                        s < b[k] && t > a[k]
                    }
                })
                .collect();
            assert_eq!(&w.indices[w.offsets[k]..w.offsets[k + 1]], &want[..]);
        }
        let agg = aggregate(&net, None).unwrap();
        assert_eq!(agg.count.iter().sum::<f64>(), c.src.len() as f64);
        let total: f64 = c.start.iter().zip(&c.end).map(|(s, e)| e - s).sum();
        assert!(same(agg.duration.iter().sum::<f64>(), total));
        // Clipped durations over a tumbling series add up to the total.
        let (a, b) = regular_windows(0.0, 10.0, 1.0, 1.0).unwrap();
        let mut clipped = 0.0;
        for k in 0..a.len() {
            clipped += aggregate(&net, Some((a[k], b[k])))
                .unwrap()
                .duration
                .iter()
                .sum::<f64>();
        }
        assert!(same(clipped, total));
    }
}

/// Earliest arrival without a waiting limit runs 64 roots per pass; a
/// waiting limit too large to matter forces the one-root sweep, which
/// must agree, across batch boundaries and in both directions.
#[test]
fn batched_sweeps_match_single_root_sweeps() {
    let mut rng = Rng(17);
    for _ in 0..20 {
        let n = 60 + rng.below(100) as usize;
        let m = rng.below(4 * n as u64) as usize;
        let mut c = Case {
            n,
            src: (0..m).map(|_| rng.below(n as u64) as i64).collect(),
            tgt: (0..m).map(|_| rng.below(n as u64) as i64).collect(),
            start: (0..m).map(|_| rng.below(20) as f64).collect(),
            end: Vec::new(),
            directed: rng.below(2) == 0,
        };
        c.end = c.start.iter().map(|s| s + rng.below(3) as f64).collect();
        let net = c.net();
        let mut o = random_options(&mut rng);
        o.max_wait = f64::INFINITY;
        let mut huge = o;
        huge.max_wait = 1e300;
        let roots: Vec<usize> = (0..n).rev().collect();
        for (crit, dir) in [
            (Criterion::EarliestArrival, Direction::Out),
            (Criterion::LatestDeparture, Direction::In),
            (Criterion::EarliestArrival, Direction::In),
        ] {
            let a = temporal_paths(&net, crit, dir, &roots, &o).unwrap();
            let b = temporal_paths(&net, crit, dir, &roots, &huge).unwrap();
            assert_eq!(a, b, "{crit:?} {dir:?} {o:?}");
        }
        for dir in [Direction::Out, Direction::In] {
            let a = temporal_reach_counts(&net, dir, &o).unwrap();
            let b = temporal_reach_counts(&net, dir, &huge).unwrap();
            assert_eq!(a, b);
        }
    }
}
