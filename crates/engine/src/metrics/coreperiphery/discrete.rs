//! The discrete core-periphery model.
//!
//! Borgatti and Everett (1999, `borgatti1999models`) compare a network
//! with an idealized core-periphery pattern: a partition of the nodes into
//! a core and a periphery induces a pattern matrix delta with delta_ij = 1
//! when i and j are both in the core and 0 when both are in the periphery.
//! For pairs with one node in each, the paper gives three choices: 1 (its
//! equation 2, a core node is tied to everyone), 0 (equation 3), or no
//! value at all, so that those pairs are left out (equation 4, which the
//! authors recommend and UCINET uses by default). The fit of a partition is
//! the Pearson correlation between the adjacency matrix and delta over the
//! pairs that have a pattern value, and the core is the partition that
//! maximizes it.
//!
//! Search. The authors use a genetic algorithm and note that any
//! combinatorial optimization method will do. This module uses tabu search
//! (Glover 1989, `glover1989tabu`) over single-node moves between core and
//! periphery, from several random starts: each step makes the best move
//! that is not tabu (a tabu move is allowed only if it beats the best
//! partition found so far), a node that moves may not move back for
//! `tenure` steps, and a start ends after `max_stall` steps without
//! improvement. The best partition over all starts is returned, the
//! earliest start winning ties. The correlation is updated in constant
//! time per move. Like any heuristic this can miss the global optimum; the
//! validation suite compares the fit with exhaustive search on small
//! networks.
//!
//! Conventions:
//! - Pairs of distinct nodes only: the diagonal (self-loops) is ignored,
//!   as in the paper's analyses. Undirected networks use each unordered
//!   pair once (the upper half of the matrix); directed networks use both
//!   ordered pairs, and the pattern is the same in both directions.
//! - Binary by default: a pair has value 1 when any edge joins the two
//!   nodes in that direction (parallel edges count once). With weights,
//!   the value is the total weight of the edges, and weights must be
//!   non-negative.
//! - The correlation is undefined (and such partitions are never chosen)
//!   when the included values or the pattern values are all equal, for
//!   example when every node is in the core. A network whose values are
//!   all equal has no core-periphery fit and is rejected.
//! - Starts use streams (seed, start) of [`crate::rng::Stream`] and run in
//!   parallel, so results do not depend on the number of threads.

use super::PairRows;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;
use crate::rng::Stream;

/// The pattern value of pairs with one node in the core and one in the
/// periphery.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoreLinks {
    /// Left out of the fit (Borgatti and Everett's equation 4; UCINET's
    /// default).
    #[default]
    Missing,
    /// Expected to be tied (equation 2).
    Ones,
    /// Expected not to be tied (equation 3).
    Zeros,
}

impl CoreLinks {
    /// Decodes the C interface codes: 0 = missing, 1 = ones, 2 = zeros.
    pub fn from_code(code: i32) -> Result<CoreLinks> {
        match code {
            0 => Ok(CoreLinks::Missing),
            1 => Ok(CoreLinks::Ones),
            2 => Ok(CoreLinks::Zeros),
            _ => Err(Error::InvalidArgument(format!(
                "unknown code {code} for core-periphery ties"
            ))),
        }
    }
}

/// Options for [`core_periphery`].
#[derive(Debug, Clone, Copy)]
pub struct CorePeripheryOptions {
    /// Pattern value of core-periphery pairs.
    pub links: CoreLinks,
    /// Use edge weights as values.
    pub weighted: bool,
    /// Number of random starting partitions.
    pub starts: usize,
    /// Steps without improvement after which a start ends.
    pub max_stall: usize,
    /// Steps for which a node may not move back; 0 chooses max(7, n / 4).
    pub tenure: usize,
    /// Seed for the starting partitions.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for CorePeripheryOptions {
    fn default() -> Self {
        CorePeripheryOptions {
            links: CoreLinks::Missing,
            weighted: false,
            starts: 20,
            max_stall: 200,
            tenure: 0,
            seed: 1,
            threads: 0,
        }
    }
}

/// Result of [`core_periphery`].
#[derive(Debug, Clone, PartialEq)]
pub struct CorePeriphery {
    /// Whether each node is in the core.
    pub core: Vec<bool>,
    /// Correlation between the data and the ideal pattern.
    pub correlation: f64,
}

/// Sums over the pairs in the fit.
#[derive(Debug, Clone, Copy, Default)]
struct Sums {
    n: f64,
    a: f64,
    aa: f64,
    d: f64,
    ad: f64,
}

impl Sums {
    fn correlation(&self) -> f64 {
        let den =
            ((self.n * self.aa - self.a * self.a) * (self.n * self.d - self.d * self.d)).sqrt();
        if den > 0.0 && den.is_finite() {
            (self.n * self.ad - self.a * self.d) / den
        } else {
            f64::NAN
        }
    }
}

struct Model<'a> {
    pr: &'a PairRows,
    links: CoreLinks,
    total1: Vec<f64>,
    total2: Vec<f64>,
}

struct State {
    core: Vec<bool>,
    ncore: usize,
    /// Sums of s1 and s2 over each node's core neighbors.
    c1: Vec<f64>,
    c2: Vec<f64>,
    sums: Sums,
}

impl Model<'_> {
    fn state(&self, core: Vec<bool>) -> State {
        let n = core.len();
        let mut c1 = vec![0.0; n];
        let mut c2 = vec![0.0; n];
        for v in 0..n {
            for &(u, s1, s2) in &self.pr.rows[v] {
                if core[u as usize] {
                    c1[v] += s1;
                    c2[v] += s2;
                }
            }
        }
        let ncore = core.iter().filter(|&&c| c).count();
        let mut st = State {
            core,
            ncore,
            c1,
            c2,
            sums: Sums::default(),
        };
        st.sums = self.full_sums(&st);
        st
    }

    /// The sums computed from scratch.
    fn full_sums(&self, st: &State) -> Sums {
        let per = self.pr.per_node_pair as f64;
        let n = st.core.len() as f64;
        let nc = st.ncore as f64;
        let np = n - nc;
        // Totals over core-core, periphery-periphery, and mixed pairs.
        let (mut cc1, mut cc2, mut pp1, mut pp2) = (0.0, 0.0, 0.0, 0.0);
        for v in 0..st.core.len() {
            let p1 = self.total1[v] - st.c1[v];
            let p2 = self.total2[v] - st.c2[v];
            if st.core[v] {
                cc1 += st.c1[v];
                cc2 += st.c2[v];
            } else {
                pp1 += p1;
                pp2 += p2;
            }
        }
        // Every pair was counted from both ends.
        let (cc1, cc2, pp1, pp2) = (cc1 / 2.0, cc2 / 2.0, pp1 / 2.0, pp2 / 2.0);
        let ncc = per * nc * (nc - 1.0) / 2.0;
        let npp = per * np * (np - 1.0) / 2.0;
        let all = self.pr.pairs;
        match self.links {
            CoreLinks::Missing => Sums {
                n: ncc + npp,
                a: cc1 + pp1,
                aa: cc2 + pp2,
                d: ncc,
                ad: cc1,
            },
            CoreLinks::Zeros => Sums {
                n: all,
                a: self.pr.sum,
                aa: self.pr.sum_sq,
                d: ncc,
                ad: cc1,
            },
            CoreLinks::Ones => Sums {
                n: all,
                a: self.pr.sum,
                aa: self.pr.sum_sq,
                d: all - npp,
                ad: self.pr.sum - pp1,
            },
        }
    }

    /// The sums after node v changes side.
    fn flipped(&self, st: &State, v: usize) -> Sums {
        let per = self.pr.per_node_pair as f64;
        let n = st.core.len() as f64;
        let in_core = st.core[v];
        // Pairs between v and the other core nodes, and the other
        // periphery nodes.
        let nc = st.ncore as f64 - if in_core { 1.0 } else { 0.0 };
        let np = n - 1.0 - nc;
        let (kc, kp) = (per * nc, per * np);
        let (c1, c2) = (st.c1[v], st.c2[v]);
        let (p1, p2) = (self.total1[v] - c1, self.total2[v] - c2);
        let mut s = st.sums;
        // sign = +1 moving into the core, -1 moving out.
        let sg = if in_core { -1.0 } else { 1.0 };
        match self.links {
            CoreLinks::Missing => {
                // Pairs with the core switch between excluded and core-core;
                // pairs with the periphery between periphery-periphery and
                // excluded.
                s.n += sg * (kc - kp);
                s.a += sg * (c1 - p1);
                s.aa += sg * (c2 - p2);
                s.d += sg * kc;
                s.ad += sg * c1;
            }
            CoreLinks::Ones => {
                // Pairs with the periphery switch between 0 and 1.
                s.d += sg * kp;
                s.ad += sg * p1;
            }
            CoreLinks::Zeros => {
                // Pairs with the core switch between 0 and 1.
                s.d += sg * kc;
                s.ad += sg * c1;
            }
        }
        s
    }

    fn flip(&self, st: &mut State, v: usize) {
        st.sums = self.flipped(st, v);
        let into = !st.core[v];
        st.core[v] = into;
        if into {
            st.ncore += 1;
        } else {
            st.ncore -= 1;
        }
        let sg = if into { 1.0 } else { -1.0 };
        for &(u, s1, s2) in &self.pr.rows[v] {
            st.c1[u as usize] += sg * s1;
            st.c2[u as usize] += sg * s2;
        }
    }
}

fn value(s: &Sums) -> f64 {
    let r = s.correlation();
    if r.is_nan() {
        f64::NEG_INFINITY
    } else {
        r
    }
}

fn search(
    model: &Model,
    opts: &CorePeripheryOptions,
    tenure: usize,
    start: usize,
) -> (f64, Vec<bool>) {
    let n = model.total1.len();
    let mut rng = Stream::new(opts.seed, start as u64);
    let core: Vec<bool> = (0..n).map(|_| rng.next_u32() & 1 == 1).collect();
    let mut st = model.state(core);
    let mut best = value(&st.sums);
    let mut best_core = st.core.clone();
    let mut tabu = vec![0usize; n];
    let (mut stall, mut step) = (0usize, 0usize);
    let mut since_refresh = 0usize;
    while stall < opts.max_stall {
        step += 1;
        let mut chosen: Option<(f64, usize)> = None;
        for v in 0..n {
            let s = value(&model.flipped(&st, v));
            let allowed = tabu[v] < step || s > best;
            if allowed && chosen.map_or(true, |(b, _)| s > b) {
                chosen = Some((s, v));
            }
        }
        let Some((s, v)) = chosen else { break };
        model.flip(&mut st, v);
        tabu[v] = step + tenure;
        // Refresh the running sums now and then so rounding cannot build up.
        since_refresh += 1;
        if since_refresh == 64 {
            st.sums = model.full_sums(&st);
            since_refresh = 0;
        }
        if s > best {
            best = s;
            best_core.clone_from(&st.core);
            stall = 0;
        } else {
            stall += 1;
        }
    }
    // Report the exact fit of the best partition.
    let exact = model.state(best_core.clone()).sums.correlation();
    (exact, best_core)
}

/// Fits the discrete core-periphery model. See the module documentation.
pub fn core_periphery(g: &Graph, opts: CorePeripheryOptions) -> Result<CorePeriphery> {
    let n = g.node_count();
    if opts.starts == 0 {
        return Err(Error::InvalidArgument(
            "the number of starts must be at least 1".into(),
        ));
    }
    let pr = PairRows::new(g, opts.weighted)?;
    if n < 2 || pr.pairs * pr.sum_sq - pr.sum * pr.sum <= 0.0 {
        return Err(Error::InvalidArgument(
            "a core-periphery fit needs at least two nodes and pairs with different values \
             (some tied and some not)"
                .into(),
        ));
    }
    let total1: Vec<f64> = pr
        .rows
        .iter()
        .map(|r| r.iter().map(|p| p.1).sum())
        .collect();
    let total2: Vec<f64> = pr
        .rows
        .iter()
        .map(|r| r.iter().map(|p| p.2).sum())
        .collect();
    let model = Model {
        pr: &pr,
        links: opts.links,
        total1,
        total2,
    };
    let tenure = if opts.tenure == 0 {
        7.max(n / 4)
    } else {
        opts.tenure
    };
    let parts = parallel::for_each_item(opts.starts, opts.threads, Vec::new, |acc, s| {
        acc.push((s, search(&model, &opts, tenure, s)));
    });
    let mut runs: Vec<(usize, (f64, Vec<bool>))> = parts.into_iter().flatten().collect();
    runs.sort_by_key(|r| r.0);
    let mut best: Option<(f64, Vec<bool>)> = None;
    for (_, run) in runs {
        let better = match &best {
            None => true,
            Some(b) => run.0 > b.0 || (b.0.is_nan() && !run.0.is_nan()),
        };
        if better {
            best = Some(run);
        }
    }
    let (correlation, core) = best.expect("at least one start");
    Ok(CorePeriphery { core, correlation })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    /// Correlation of a partition, from the definition.
    fn direct(g: &Graph, core: &[bool], links: CoreLinks) -> f64 {
        let n = g.node_count();
        let mut a = vec![0.0; n * n];
        for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
            if s != t {
                a[s as usize * n + t as usize] = 1.0;
                if !g.is_directed() {
                    a[t as usize * n + s as usize] = 1.0;
                }
            }
        }
        let (mut x, mut y) = (Vec::new(), Vec::new());
        for i in 0..n {
            for j in 0..n {
                if i == j || (!g.is_directed() && j < i) {
                    continue;
                }
                let d = match (core[i], core[j], links) {
                    (true, true, _) => 1.0,
                    (false, false, _) => 0.0,
                    (_, _, CoreLinks::Missing) => continue,
                    (_, _, CoreLinks::Ones) => 1.0,
                    (_, _, CoreLinks::Zeros) => 0.0,
                };
                x.push(a[i * n + j]);
                y.push(d);
            }
        }
        let m = x.len() as f64;
        let (mx, my) = (x.iter().sum::<f64>() / m, y.iter().sum::<f64>() / m);
        let sxy: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
        let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
        let syy: f64 = y.iter().map(|b| (b - my) * (b - my)).sum();
        sxy / (sxx * syy).sqrt()
    }

    fn exhaustive(g: &Graph, links: CoreLinks) -> f64 {
        let n = g.node_count();
        (0..(1u32 << n))
            .map(|mask| {
                let core: Vec<bool> = (0..n).map(|v| mask >> v & 1 == 1).collect();
                direct(g, &core, links)
            })
            .filter(|r| !r.is_nan())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    #[test]
    fn star_has_its_center_as_core() {
        let g = und(6, &[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2)]);
        let r = core_periphery(&g, CorePeripheryOptions::default()).unwrap();
        assert!((r.correlation - exhaustive(&g, CoreLinks::Missing)).abs() < 1e-12);
        let core: Vec<bool> = r.core.clone();
        assert!((direct(&g, &core, CoreLinks::Missing) - r.correlation).abs() < 1e-12);
    }

    #[test]
    fn matches_exhaustive_search() {
        let mut rng = Stream::new(5, 0);
        for trial in 0..12 {
            let n = 9;
            let directed = trial % 2 == 1;
            let (mut s, mut t) = (Vec::new(), Vec::new());
            for a in 0..n as i64 {
                for b in 0..n as i64 {
                    if a == b || (!directed && b < a) {
                        continue;
                    }
                    let p = if a < 3 || b < 3 { 0.6 } else { 0.15 };
                    if rng.uniform() < p {
                        s.push(a);
                        t.push(b);
                    }
                }
            }
            let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
            for links in [CoreLinks::Missing, CoreLinks::Ones, CoreLinks::Zeros] {
                let opts = CorePeripheryOptions {
                    links,
                    threads: 2,
                    ..Default::default()
                };
                let r = core_periphery(&g, opts).unwrap();
                let want = exhaustive(&g, links);
                assert!(
                    (r.correlation - want).abs() < 1e-12,
                    "{links:?}: {} vs {want}",
                    r.correlation
                );
                assert!((direct(&g, &r.core, links) - r.correlation).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn rejects_constant_networks() {
        assert!(core_periphery(&und(4, &[]), CorePeripheryOptions::default()).is_err());
        let k3 = und(3, &[(0, 1), (1, 2), (0, 2)]);
        assert!(core_periphery(&k3, CorePeripheryOptions::default()).is_err());
    }
}
