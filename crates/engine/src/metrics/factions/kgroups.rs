//! Factions: a partition of the nodes into a given number of groups that
//! best matches the ideal of cohesive, mutually separate groups.
//!
//! UCINET's factions routine (Borgatti, Everett, and Freeman 2002,
//! `borgatti2002ucinet`) compares the network with an idealized pattern in
//! which every pair of nodes in the same group is tied and no pair in
//! different groups is, and counts as errors the absent ties within groups
//! plus the present ties between groups; the final number of errors
//! measures how well the partition fits (Hanneman and Riddle 2005,
//! `hanneman2005introduction`, chapter 11). This is the Hamming distance
//! between the adjacency matrix and the ideal block pattern, the criterion
//! of a binary blockmodel with complete diagonal and null off-diagonal
//! blocks (Doreian, Batagelj, and Ferligoj 2005, `doreian2005generalized`).
//! The alternative criterion is the Pearson correlation between the
//! adjacency matrix and the ideal pattern, which Borgatti and Everett
//! (1999, `borgatti1999models`) use to fit core-periphery patterns.
//!
//! Search. Finding the best partition is a combinatorial optimization
//! problem. UCINET's documentation does not specify its search procedure
//! in enough detail to reproduce it, so this module uses its own: tabu
//! search (Glover 1989, `glover1989tabu`) over moves of single nodes
//! between groups, from several random starting partitions. At each step
//! the best move that is not tabu is made (a tabu move is allowed only if
//! it improves on the best partition found so far); moving a node away
//! from a group forbids moving it back for `tenure` steps; a start ends
//! after `max_stall` steps without improving on its best partition. The
//! best partition over all starts is returned, the earliest start winning
//! ties. Like any heuristic this can miss the global optimum; the
//! validation suite compares the fit with exhaustive search on small
//! networks.
//!
//! Conventions:
//! - The adjacency matrix is binary: parallel edges count once, self-loops
//!   are ignored, and weights are not used. Undirected networks compare
//!   each unordered pair once; directed networks compare both ordered
//!   pairs.
//! - Every group is non-empty (k groups exactly, 1 <= k <= n). Groups are
//!   labelled 0, 1, ... in order of their first node.
//! - With m ties among P pairs, W pairs within groups, and T ties within
//!   groups, errors = (W - T) + (m - T) and the correlation is
//!   (P T - m W) / sqrt((P m - m^2)(P W - W^2)), which is undefined (NaN)
//!   when there are no ties, every pair is tied, or the partition puts all
//!   pairs within groups (k = 1) or none (k = n). The correlation
//!   criterion therefore requires 1 < k < n and 0 < m < P.
//! - Starts use streams (seed, start) of [`crate::rng::Stream`] and run in
//!   parallel, so results do not depend on the number of threads.

use super::simple_rows;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::community::dense_labels;
use crate::parallel;
use crate::rng::Stream;

/// The criterion a faction partition optimizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FactionFit {
    /// Fewest errors: absent ties within groups plus ties between groups.
    #[default]
    Errors,
    /// Largest correlation between the adjacency matrix and the ideal
    /// pattern.
    Correlation,
}

impl FactionFit {
    /// Decodes the C interface codes: 0 = errors, 1 = correlation.
    pub fn from_code(code: i32) -> Result<FactionFit> {
        match code {
            0 => Ok(FactionFit::Errors),
            1 => Ok(FactionFit::Correlation),
            _ => Err(Error::InvalidArgument(format!(
                "unknown faction fit code {code}"
            ))),
        }
    }
}

/// Options for [`factions`].
#[derive(Debug, Clone, Copy)]
pub struct FactionOptions {
    /// Number of groups.
    pub k: usize,
    /// The criterion.
    pub fit: FactionFit,
    /// Number of random starting partitions.
    pub starts: usize,
    /// Steps without improvement after which a start ends.
    pub max_stall: usize,
    /// Steps for which a node may not return to a group it left; 0 chooses
    /// max(7, n / 4).
    pub tenure: usize,
    /// Seed for the starting partitions.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for FactionOptions {
    fn default() -> Self {
        FactionOptions {
            k: 2,
            fit: FactionFit::Errors,
            starts: 20,
            max_stall: 200,
            tenure: 0,
            seed: 1,
            threads: 0,
        }
    }
}

/// Result of [`factions`].
#[derive(Debug, Clone, PartialEq)]
pub struct Factions {
    /// Group of each node, 0..k, numbered by first appearance.
    pub membership: Vec<usize>,
    /// Number of errors of the partition.
    pub errors: u64,
    /// Correlation between the adjacency matrix and the ideal pattern
    /// (NaN where undefined).
    pub correlation: f64,
}

/// Counts that determine both criteria.
struct Counts {
    pairs: f64,
    ties: f64,
    /// 1 for undirected (unordered pairs), 2 for directed (ordered pairs).
    w: i64,
}

impl Counts {
    fn errors(&self, within: i64, tied_within: i64) -> i64 {
        within - tied_within + self.ties as i64 - tied_within
    }

    fn correlation(&self, within: i64, tied_within: i64) -> f64 {
        let (p, m, w, t) = (self.pairs, self.ties, within as f64, tied_within as f64);
        let den = ((p * m - m * m) * (p * w - w * w)).sqrt();
        if den > 0.0 {
            (p * t - m * w) / den
        } else {
            f64::NAN
        }
    }
}

/// Binary adjacency in the form the search needs: for undirected networks
/// the simple rows; for directed networks, for each node, its out- and
/// in-neighbors (each ordered pair counts on its own).
struct Ties {
    rows: Vec<Vec<u32>>,
}

impl Ties {
    fn new(g: &Graph) -> Ties {
        if !g.is_directed() {
            return Ties {
                rows: simple_rows(g),
            };
        }
        let n = g.node_count();
        let rows = (0..n)
            .map(|v| {
                let mut out: Vec<u32> = g.out().neighbors(v).to_vec();
                out.dedup();
                out.retain(|&u| u as usize != v);
                let mut inc: Vec<u32> = g.inc().neighbors(v).to_vec();
                inc.dedup();
                inc.retain(|&u| u as usize != v);
                // A mutual pair appears twice: once per arc.
                out.extend(inc);
                out
            })
            .collect();
        Ties { rows }
    }
}

struct State {
    group: Vec<usize>,
    size: Vec<i64>,
    /// link[v * k + c]: ties (arcs, for directed networks) between v and
    /// members of group c.
    link: Vec<i64>,
    within: i64,
    tied_within: i64,
}

impl State {
    fn new(ties: &Ties, group: Vec<usize>, k: usize, w: i64) -> State {
        let n = group.len();
        let mut size = vec![0i64; k];
        for &c in &group {
            size[c] += 1;
        }
        let mut link = vec![0i64; n * k];
        let mut tied = 0i64;
        for v in 0..n {
            for &u in &ties.rows[v] {
                link[v * k + group[u as usize]] += 1;
                if group[u as usize] == group[v] {
                    tied += 1;
                }
            }
        }
        let within = size.iter().map(|&s| w * s * (s - 1) / 2).sum();
        // Every within tie was counted from both ends (both arcs of a
        // directed pair are listed at each end, so the same holds).
        State {
            group,
            size,
            link,
            within,
            tied_within: tied / 2,
        }
    }

    /// Changes in (within pairs, ties within) if v moves to group h.
    fn delta(&self, v: usize, h: usize, k: usize, w: i64) -> (i64, i64) {
        let g = self.group[v];
        let dw = w * (self.size[h] - (self.size[g] - 1));
        let dt = self.link[v * k + h] - self.link[v * k + g];
        (dw, dt)
    }

    fn apply(&mut self, ties: &Ties, v: usize, h: usize, k: usize, w: i64) {
        let (dw, dt) = self.delta(v, h, k, w);
        let g = self.group[v];
        self.within += dw;
        self.tied_within += dt;
        self.size[g] -= 1;
        self.size[h] += 1;
        self.group[v] = h;
        for &u in &ties.rows[v] {
            let u = u as usize;
            self.link[u * k + g] -= 1;
            self.link[u * k + h] += 1;
        }
    }
}

/// Score to maximize, exact for errors.
#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct Score(f64);

fn score(fit: FactionFit, c: &Counts, within: i64, tied: i64) -> Score {
    match fit {
        FactionFit::Errors => Score(-(c.errors(within, tied) as f64)),
        FactionFit::Correlation => {
            let r = c.correlation(within, tied);
            Score(if r.is_nan() { f64::NEG_INFINITY } else { r })
        }
    }
}

fn random_partition(n: usize, k: usize, rng: &mut Stream) -> Vec<usize> {
    let mut order: Vec<usize> = vec![0; n];
    rng.permutation(&mut order);
    let mut group = vec![0usize; n];
    for (i, &v) in order.iter().enumerate() {
        // The first k nodes in the random order seed one group each.
        group[v] = if i < k {
            i
        } else {
            rng.below(k as u64) as usize
        };
    }
    group
}

fn tabu_search(
    ties: &Ties,
    c: &Counts,
    opts: &FactionOptions,
    tenure: usize,
    start: usize,
) -> (Score, Vec<usize>) {
    let n = ties.rows.len();
    let k = opts.k;
    let mut rng = Stream::new(opts.seed, start as u64);
    let mut st = State::new(ties, random_partition(n, k, &mut rng), k, c.w);
    let mut best = score(opts.fit, c, st.within, st.tied_within);
    let mut best_group = st.group.clone();
    // tabu[v * k + g]: the step until which v may not return to g.
    let mut tabu = vec![0usize; n * k];
    let mut stall = 0;
    let mut step = 0usize;
    while stall < opts.max_stall {
        step += 1;
        let mut chosen: Option<(Score, usize, usize)> = None;
        for v in 0..n {
            let g = st.group[v];
            if st.size[g] == 1 {
                continue;
            }
            for h in 0..k {
                if h == g {
                    continue;
                }
                let (dw, dt) = st.delta(v, h, k, c.w);
                let s = score(opts.fit, c, st.within + dw, st.tied_within + dt);
                let allowed = tabu[v * k + h] < step || s > best;
                if allowed && chosen.map_or(true, |(b, _, _)| s > b) {
                    chosen = Some((s, v, h));
                }
            }
        }
        let Some((s, v, h)) = chosen else { break };
        let g = st.group[v];
        st.apply(ties, v, h, k, c.w);
        tabu[v * k + g] = step + tenure;
        if s > best {
            best = s;
            best_group.clone_from(&st.group);
            stall = 0;
        } else {
            stall += 1;
        }
    }
    (best, best_group)
}

/// Factions of `g`. See the module documentation.
pub fn factions(g: &Graph, opts: FactionOptions) -> Result<Factions> {
    let n = g.node_count();
    let k = opts.k;
    if k == 0 || k > n {
        return Err(Error::InvalidArgument(format!(
            "the number of factions must be between 1 and the number of nodes ({n}), not {k}"
        )));
    }
    if opts.starts == 0 {
        return Err(Error::InvalidArgument(
            "the number of starts must be at least 1".into(),
        ));
    }
    let ties = Ties::new(g);
    let w: i64 = if g.is_directed() { 2 } else { 1 };
    let arcs: usize = ties.rows.iter().map(|r| r.len()).sum();
    let c = Counts {
        pairs: (w * (n as i64) * (n as i64 - 1) / 2) as f64,
        ties: (arcs / 2) as f64,
        w,
    };
    if opts.fit == FactionFit::Correlation {
        if k == 1 || k == n {
            return Err(Error::InvalidArgument(
                "the correlation criterion needs more than one faction and fewer factions \
                 than nodes"
                    .into(),
            ));
        }
        if c.ties == 0.0 || c.ties == c.pairs {
            return Err(Error::InvalidArgument(
                "the correlation criterion is undefined for a network with no ties or with \
                 every possible tie"
                    .into(),
            ));
        }
    }
    let tenure = if opts.tenure == 0 {
        7.max(n / 4)
    } else {
        opts.tenure
    };
    // Starts differ in length, so threads claim them as they finish; runs
    // are sorted by start below.
    let threads = parallel::threads_for(
        opts.threads,
        opts.starts,
        (opts.starts as u64).saturating_mul((n * n) as u64 * 16),
    );
    let parts = parallel::for_each_claimed(opts.starts, threads, Vec::new, |acc, s| {
        acc.push((s, tabu_search(&ties, &c, &opts, tenure, s)));
    });
    let mut runs: Vec<(usize, (Score, Vec<usize>))> = parts.into_iter().flatten().collect();
    runs.sort_by_key(|r| r.0);
    let mut best: Option<(Score, Vec<usize>)> = None;
    for (_, run) in runs {
        if best.as_ref().map_or(true, |b| run.0 > b.0) {
            best = Some(run);
        }
    }
    let group = best.expect("at least one start").1;
    let membership = dense_labels(&group).0;
    let st = State::new(&ties, membership.clone(), k, w);
    Ok(Factions {
        membership,
        errors: c.errors(st.within, st.tied_within) as u64,
        correlation: c.correlation(st.within, st.tied_within),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    fn opts(k: usize, fit: FactionFit) -> FactionOptions {
        FactionOptions {
            k,
            fit,
            threads: 2,
            ..Default::default()
        }
    }

    #[test]
    fn two_cliques_and_a_bridge() {
        let g = und(6, &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)]);
        let f = factions(&g, opts(2, FactionFit::Errors)).unwrap();
        assert_eq!(f.membership, vec![0, 0, 0, 1, 1, 1]);
        assert_eq!(f.errors, 1);
        // P = 15, m = 7, W = 6, T = 6.
        let r = (15.0 * 6.0 - 7.0 * 6.0) / ((15.0 * 7.0 - 49.0) * (15.0 * 6.0 - 36.0f64)).sqrt();
        assert!((f.correlation - r).abs() < 1e-15);
        let f = factions(&g, opts(2, FactionFit::Correlation)).unwrap();
        assert_eq!(f.membership, vec![0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn directed_counts_ordered_pairs() {
        // 0 <-> 1 and 2 <-> 3, plus 1 -> 2 and a parallel arc and a loop.
        let d = Graph::from_edges(
            4,
            &[0, 1, 2, 3, 1, 1, 0],
            &[1, 0, 3, 2, 2, 2, 0],
            None,
            true,
        )
        .unwrap();
        let f = factions(&d, opts(2, FactionFit::Errors)).unwrap();
        assert_eq!(f.membership, vec![0, 0, 1, 1]);
        assert_eq!(f.errors, 1);
    }

    #[test]
    fn edge_cases() {
        let g = und(4, &[(0, 1)]);
        assert!(factions(&g, opts(0, FactionFit::Errors)).is_err());
        assert!(factions(&g, opts(5, FactionFit::Errors)).is_err());
        let f = factions(&g, opts(1, FactionFit::Errors)).unwrap();
        assert_eq!(f.membership, vec![0; 4]);
        assert_eq!(f.errors, 5);
        assert!(f.correlation.is_nan());
        assert!(factions(&g, opts(1, FactionFit::Correlation)).is_err());
        assert!(factions(&und(4, &[]), opts(2, FactionFit::Correlation)).is_err());
        let f = factions(&g, opts(4, FactionFit::Errors)).unwrap();
        assert_eq!(f.membership, vec![0, 1, 2, 3]);
        assert_eq!(f.errors, 1);
    }

    #[test]
    fn thread_count_does_not_matter() {
        let mut e = Vec::new();
        let mut rng = Stream::new(3, 0);
        for a in 0..30i64 {
            for b in (a + 1)..30 {
                let p = if a / 10 == b / 10 { 0.6 } else { 0.08 };
                if rng.uniform() < p {
                    e.push((a, b));
                }
            }
        }
        let g = und(30, &e);
        let mut o = opts(3, FactionFit::Errors);
        o.threads = 1;
        let a = factions(&g, o).unwrap();
        o.threads = 4;
        assert_eq!(a, factions(&g, o).unwrap());
    }
}
