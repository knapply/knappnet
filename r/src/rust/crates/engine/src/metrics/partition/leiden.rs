//! Leiden community detection.
//!
//! Traag, Waltman, and van Eck (2019, `traag2019louvain`) maximize, for an
//! undirected network with total edge weight m, node weights n_i, and
//! resolution gamma,
//!
//! ```text
//! Q = 1/(2m) sum_ij (A_ij - gamma n_i n_j) delta(c_i, c_j),
//! ```
//!
//! the form igraph's `cluster_leiden()` documents. With n_i the weighted
//! degree and gamma replaced by gamma / (2m) this is modularity
//! ([`crate::metrics::community::modularity`]); with n_i = 1 it is the
//! Constant Potts Model (CPM) of Traag, Van Dooren, and Nesterov (2011,
//! `traag2011narrow`), whose usual form sum_c (e_c - gamma n_c^2 / 2)
//! (NetworkX's `constant_potts_model()`) is m Q.
//!
//! One iteration follows Algorithm A.2 of the paper's supplement:
//! 1. Fast local moving: nodes wait in a queue (initially all, in random
//!    order); the node at its head moves to the neighboring or empty group
//!    that most increases Q, if any increases it, and its neighbors outside
//!    its new group join the queue if they are not in it.
//! 2. Refinement: within each group S, starting from singletons, every
//!    node v that is well connected to S (E(v, S - v) >= gamma n_v
//!    (N_S - n_v)) and still alone is merged into a well-connected refined
//!    group C of S with probability proportional to exp(dQ / beta) over
//!    the merges with dQ >= 0 (staying alone has dQ = 0), where dQ =
//!    E(v, C) - gamma n_v N_C is measured in edge weight. Refined groups
//!    are therefore connected.
//! 3. Aggregation: each refined group becomes a node, its group in step 1
//!    gives the initial partition of the aggregate network, and steps 1-3
//!    repeat until step 1 leaves every node of the aggregate network alone.
//!
//! Iterations repeat from the partition of the previous one. The paper
//! proves that after each iteration every group is gamma-connected, so
//! groups are connected whenever gamma > 0.
//!
//! Conventions:
//! - Undirected networks only, as igraph's `cluster_leiden()` requires.
//! - Weights (if used) are added as given and must be non-negative;
//!   parallel edges add up. A self-loop of weight w adds 2w to A_ii and to
//!   its node's degree, as in [`crate::metrics::community::modularity`],
//!   and never affects a move.
//! - Node weights default to the weighted degree for modularity and to 1
//!   for the CPM; given node weights must be finite and non-negative. For
//!   modularity, the resolution is divided by 2m, as igraph does.
//! - `n_iterations` iterations are run; a negative value runs until an
//!   iteration changes nothing (at most 10,000 iterations). `initial`
//!   gives a starting partition (any non-negative labels); otherwise each
//!   node starts alone.
//! - In the refinement, only refined groups adjacent to v are candidates;
//!   non-adjacent ones have dQ < 0 whenever gamma n_v N_C > 0 and would
//!   break connectivity when gamma = 0. Moves in step 1 need an increase
//!   of Q greater than zero; ties between equally good groups go to the
//!   first neighbor group in node order.
//! - Random numbers come from [`crate::rng::Stream`] with the seed
//!   (iteration i uses stream i), so a seed reproduces the result on any
//!   platform. The method runs on one thread.
//! - The quality returned is Q above on the input network (NaN when it has
//!   no edges). Groups are labelled 0, 1, ... in order of their first node.
//! - Partitions are not expected to match igraph's exactly (different
//!   random numbers), only in quality.

use std::collections::VecDeque;

use super::{additive_weights, relabel, require_undirected};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::community::dense_labels;
use crate::rng::Stream;

/// The quality function of [`leiden`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LeidenObjective {
    /// Modularity (node weights: weighted degree; resolution / 2m).
    Modularity,
    /// Constant Potts Model (node weights: 1).
    Cpm,
}

impl LeidenObjective {
    /// 0 = CPM (igraph's default), 1 = modularity.
    pub fn from_code(code: i32) -> Result<LeidenObjective> {
        match code {
            0 => Ok(LeidenObjective::Cpm),
            1 => Ok(LeidenObjective::Modularity),
            _ => Err(Error::InvalidArgument(format!(
                "unknown Leiden objective code {code}"
            ))),
        }
    }
}

/// Options for [`leiden`].
#[derive(Debug, Clone, Copy)]
pub struct LeidenOptions<'a> {
    /// Quality function.
    pub objective: LeidenObjective,
    /// Resolution gamma.
    pub resolution: f64,
    /// Randomness of the refinement.
    pub beta: f64,
    /// Number of iterations; negative runs until stable.
    pub n_iterations: i64,
    /// Use edge weights.
    pub weighted: bool,
    /// Node weights (default from the objective).
    pub node_weights: Option<&'a [f64]>,
    /// Starting partition (default: singletons).
    pub initial: Option<&'a [i64]>,
    /// Seed of the random stream.
    pub seed: u64,
}

/// Result of [`leiden`].
#[derive(Debug, Clone, PartialEq)]
pub struct Leiden {
    /// Group of each node.
    pub membership: Vec<usize>,
    /// Quality Q of the partition (see the module documentation).
    pub quality: f64,
    /// Iterations run.
    pub iterations: usize,
}

const MAX_ITERATIONS: usize = 10_000;

/// A level of the method: adjacency without loops (parallel edges merged),
/// node weights.
#[derive(Clone)]
struct Level {
    adj: Vec<Vec<(usize, f64)>>,
    nw: Vec<f64>,
}

impl Level {
    fn n(&self) -> usize {
        self.nw.len()
    }

    /// The network whose nodes are the groups of `part` (k groups).
    fn aggregate(&self, part: &[usize], k: usize) -> Level {
        let mut nw = vec![0.0; k];
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); k];
        for v in 0..self.n() {
            nw[part[v]] += self.nw[v];
            for &(u, x) in &self.adj[v] {
                if part[u] != part[v] {
                    rows[part[v]].push((part[u], x));
                }
            }
        }
        for row in rows.iter_mut() {
            row.sort_unstable_by_key(|&(u, _)| u);
            let mut merged: Vec<(usize, f64)> = Vec::with_capacity(row.len());
            for &(u, x) in row.iter() {
                match merged.last_mut() {
                    Some(last) if last.0 == u => last.1 += x,
                    _ => merged.push((u, x)),
                }
            }
            *row = merged;
        }
        Level { adj: rows, nw }
    }
}

/// Weights from v to each neighbouring group, gathered without clearing a
/// full array.
struct Gather {
    w: Vec<f64>,
    seen: Vec<bool>,
    list: Vec<usize>,
}

impl Gather {
    fn new(n: usize) -> Gather {
        Gather {
            w: vec![0.0; n],
            seen: vec![false; n],
            list: Vec::new(),
        }
    }

    fn collect(&mut self, adj: &[(usize, f64)], group: &[usize], keep: impl Fn(usize) -> bool) {
        for &(u, x) in adj {
            if !keep(u) {
                continue;
            }
            let c = group[u];
            if !self.seen[c] {
                self.seen[c] = true;
                self.list.push(c);
            }
            self.w[c] += x;
        }
    }

    fn clear(&mut self) {
        for &c in &self.list {
            self.w[c] = 0.0;
            self.seen[c] = false;
        }
        self.list.clear();
    }
}

/// Fast local moving (step 1). Returns whether any node moved.
fn move_nodes(lv: &Level, comm: &mut [usize], gamma: f64, rng: &mut Stream) -> bool {
    let n = lv.n();
    let mut size = vec![0.0f64; n];
    let mut count = vec![0usize; n];
    for v in 0..n {
        size[comm[v]] += lv.nw[v];
        count[comm[v]] += 1;
    }
    let mut empty: Vec<usize> = (0..n).rev().filter(|&c| count[c] == 0).collect();
    let mut order: Vec<usize> = (0..n).collect();
    rng.permutation(&mut order);
    let mut queue: VecDeque<usize> = order.into_iter().collect();
    let mut queued = vec![true; n];
    let mut g = Gather::new(n);
    let mut moved = false;
    while let Some(v) = queue.pop_front() {
        queued[v] = false;
        let old = comm[v];
        let nv = lv.nw[v];
        g.collect(&lv.adj[v], comm, |_| true);
        // Value of being in group c (v itself left out of its own group).
        let stay = g.w[old] - gamma * nv * (size[old] - nv);
        let mut best = old;
        let mut best_val = stay;
        for &c in &g.list {
            if c == old {
                continue;
            }
            let val = g.w[c] - gamma * nv * size[c];
            if val > best_val {
                best_val = val;
                best = c;
            }
        }
        // The empty group has value 0 (only useful when v is not alone).
        let mut to_empty = false;
        if 0.0 > best_val && count[old] > 1 {
            if let Some(&e) = empty.last() {
                best = e;
                to_empty = true;
            }
        }
        g.clear();
        if best != old {
            if to_empty {
                empty.pop();
            }
            size[old] -= nv;
            count[old] -= 1;
            if count[old] == 0 {
                empty.push(old);
            }
            size[best] += nv;
            count[best] += 1;
            comm[v] = best;
            moved = true;
            for &(u, _) in &lv.adj[v] {
                if !queued[u] && comm[u] != best {
                    queued[u] = true;
                    queue.push_back(u);
                }
            }
        }
    }
    moved
}

/// Refinement (step 2): a refined partition of `comm`, numbered densely.
fn refine(
    lv: &Level,
    comm: &[usize],
    gamma: f64,
    beta: f64,
    rng: &mut Stream,
) -> (Vec<usize>, usize) {
    let n = lv.n();
    let mut csize = vec![0.0f64; n];
    for v in 0..n {
        csize[comm[v]] += lv.nw[v];
    }
    let mut refined: Vec<usize> = (0..n).collect();
    let mut rsize: Vec<f64> = lv.nw.clone();
    let mut rcount = vec![1usize; n];
    // External weight of each refined group within its group.
    let mut ext: Vec<f64> = (0..n)
        .map(|v| {
            lv.adj[v]
                .iter()
                .filter(|&&(u, _)| comm[u] == comm[v])
                .map(|&(_, x)| x)
                .sum()
        })
        .collect();
    let well = |e: f64, s: f64, total: f64| e >= gamma * s * (total - s);
    let mut order: Vec<usize> = (0..n).collect();
    rng.permutation(&mut order);
    let mut g = Gather::new(n);
    let mut cand: Vec<(usize, f64)> = Vec::new();
    for &v in &order {
        let s_total = csize[comm[v]];
        if rcount[refined[v]] != 1 || !well(ext[v], lv.nw[v], s_total) {
            continue;
        }
        let own = refined[v];
        g.collect(&lv.adj[v], &refined, |u| comm[u] == comm[v]);
        cand.clear();
        cand.push((own, 0.0));
        for &r in &g.list {
            if r == own || !well(ext[r], rsize[r], s_total) {
                continue;
            }
            let dq = g.w[r] - gamma * lv.nw[v] * rsize[r];
            if dq >= 0.0 {
                cand.push((r, dq));
            }
        }
        let w_to = |r: usize, g: &Gather| if g.seen[r] { g.w[r] } else { 0.0 };
        let top = cand
            .iter()
            .map(|&(_, d)| d)
            .fold(f64::NEG_INFINITY, f64::max);
        let total: f64 = cand.iter().map(|&(_, d)| ((d - top) / beta).exp()).sum();
        let mut x = rng.uniform() * total;
        let mut chosen = cand[cand.len() - 1].0;
        for &(r, d) in &cand {
            x -= ((d - top) / beta).exp();
            if x < 0.0 {
                chosen = r;
                break;
            }
        }
        if chosen != own {
            let wvr = w_to(chosen, &g);
            ext[chosen] = ext[chosen] + ext[v] - 2.0 * wvr;
            rsize[chosen] += lv.nw[v];
            rcount[chosen] += 1;
            rcount[own] -= 1;
            refined[v] = chosen;
        }
        g.clear();
    }
    dense_labels(&refined)
}

/// Q of a partition of the input network (see the module documentation),
/// from edge weights, loop weights, and node weights.
fn quality(base: &Level, loops: &[f64], memb: &[usize], gamma: f64, two_m: f64) -> f64 {
    if two_m == 0.0 {
        return f64::NAN;
    }
    let k = memb.iter().copied().max().map_or(0, |x| x + 1);
    let mut inside = vec![0.0f64; k];
    let mut size = vec![0.0f64; k];
    for v in 0..base.n() {
        size[memb[v]] += base.nw[v];
        inside[memb[v]] += 2.0 * loops[v];
        for &(u, x) in &base.adj[v] {
            if memb[u] == memb[v] {
                inside[memb[v]] += x;
            }
        }
    }
    let q: f64 = (0..k).map(|c| inside[c] - gamma * size[c] * size[c]).sum();
    q / two_m
}

/// Leiden communities of an undirected network. See the module
/// documentation.
pub fn leiden(g: &Graph, opts: &LeidenOptions) -> Result<Leiden> {
    require_undirected(g, "The Leiden method")?;
    if !(opts.resolution.is_finite() && opts.resolution >= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "resolution must be finite and non-negative, not {}",
            opts.resolution
        )));
    }
    if !(opts.beta.is_finite() && opts.beta > 0.0) {
        return Err(Error::InvalidArgument(format!(
            "beta must be finite and positive, not {}",
            opts.beta
        )));
    }
    let n = g.node_count();
    let w = additive_weights(g, opts.weighted)?;
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut loops = vec![0.0f64; n];
    let mut degree = vec![0.0f64; n];
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        degree[u] += x;
        degree[v] += x;
        if u == v {
            loops[u] += x;
        } else {
            rows[u].push((v, x));
            rows[v].push((u, x));
        }
    }
    let two_m: f64 = degree.iter().sum();
    let nw: Vec<f64> = match opts.node_weights {
        Some(x) => {
            if x.len() != n {
                return Err(Error::LengthMismatch {
                    what: "node weights",
                    expected: n,
                    found: x.len(),
                });
            }
            if let Some(&bad) = x.iter().find(|v| !(v.is_finite() && **v >= 0.0)) {
                return Err(Error::InvalidArgument(format!(
                    "node weights must be finite and non-negative, found {bad}"
                )));
            }
            x.to_vec()
        }
        None => match opts.objective {
            LeidenObjective::Modularity => degree.clone(),
            LeidenObjective::Cpm => vec![1.0; n],
        },
    };
    let gamma = match opts.objective {
        LeidenObjective::Modularity if two_m > 0.0 => opts.resolution / two_m,
        LeidenObjective::Modularity => 0.0,
        LeidenObjective::Cpm => opts.resolution,
    };
    let base = Level { adj: rows, nw }.aggregate(&(0..n).collect::<Vec<_>>(), n);
    let mut memb: Vec<usize> = match opts.initial {
        None => (0..n).collect(),
        Some(init) => {
            if init.len() != n {
                return Err(Error::LengthMismatch {
                    what: "initial membership",
                    expected: n,
                    found: init.len(),
                });
            }
            if let Some(&bad) = init.iter().find(|&&c| c < 0) {
                return Err(Error::InvalidArgument(format!(
                    "initial membership labels must be non-negative, found {bad}"
                )));
            }
            dense_labels(init).0
        }
    };
    let limit = if opts.n_iterations < 0 {
        MAX_ITERATIONS
    } else {
        opts.n_iterations as usize
    };
    let mut iterations = 0;
    while iterations < limit {
        let mut rng = Stream::new(opts.seed, iterations as u64);
        iterations += 1;
        let before = memb.clone();
        // `map[v]`: the node of the current level holding input node v.
        let mut map: Vec<usize> = (0..n).collect();
        let mut lv = base.clone();
        let mut comm = memb.clone();
        loop {
            move_nodes(&lv, &mut comm, gamma, &mut rng);
            let (dense, k) = dense_labels(&comm);
            comm = dense;
            if k == lv.n() {
                break;
            }
            let (refined, kr) = refine(&lv, &comm, gamma, opts.beta, &mut rng);
            let mut next_comm = vec![0usize; kr];
            for v in 0..lv.n() {
                next_comm[refined[v]] = comm[v];
            }
            for m in map.iter_mut() {
                *m = refined[*m];
            }
            lv = lv.aggregate(&refined, kr);
            comm = next_comm;
        }
        memb = relabel(&map.iter().map(|&v| comm[v]).collect::<Vec<_>>());
        if opts.n_iterations < 0 && memb == relabel(&before) {
            break;
        }
    }
    let membership = relabel(&memb);
    let q = quality(&base, &loops, &membership, gamma, two_m);
    Ok(Leiden {
        membership,
        quality: q,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::community::modularity;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    fn opts(objective: LeidenObjective, resolution: f64, seed: u64) -> LeidenOptions<'static> {
        LeidenOptions {
            objective,
            resolution,
            beta: 0.01,
            n_iterations: -1,
            weighted: false,
            node_weights: None,
            initial: None,
            seed,
        }
    }

    fn ring_of_cliques() -> Graph {
        let mut e = Vec::new();
        for c in 0..5i64 {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    e.push((4 * c + i, 4 * c + j));
                }
            }
            e.push((4 * c, (4 * c + 4) % 20));
        }
        und(20, &e)
    }

    #[test]
    fn modularity_on_ring_of_cliques() {
        let g = ring_of_cliques();
        for seed in 0..5 {
            let r = leiden(&g, &opts(LeidenObjective::Modularity, 1.0, seed)).unwrap();
            for c in 0..5 {
                assert!((0..4).all(|i| r.membership[4 * c + i] == r.membership[4 * c]));
            }
            let memb: Vec<i64> = r.membership.iter().map(|&c| c as i64).collect();
            let q = modularity(&g, &memb, 1.0, false).unwrap();
            assert!((q - r.quality).abs() < 1e-12);
            assert_eq!(*r.membership.iter().max().unwrap(), 4);
        }
    }

    #[test]
    fn cpm_resolution_extremes() {
        let g = ring_of_cliques();
        // Density inside a clique is 1, between cliques 1/16: gamma 0.5
        // separates the cliques.
        let r = leiden(&g, &opts(LeidenObjective::Cpm, 0.5, 1)).unwrap();
        assert_eq!(*r.membership.iter().max().unwrap(), 4);
        // gamma above 1 leaves every node alone.
        let r = leiden(&g, &opts(LeidenObjective::Cpm, 1.5, 1)).unwrap();
        assert_eq!(r.membership, (0..20).collect::<Vec<_>>());
        // CPM quality times m is sum_c (e_c - gamma n_c^2 / 2).
        let r = leiden(&g, &opts(LeidenObjective::Cpm, 0.5, 1)).unwrap();
        let expect = 5.0 * (6.0 - 0.5 * 16.0 / 2.0) / 35.0;
        assert!((r.quality - expect).abs() < 1e-12);
    }

    #[test]
    fn initial_membership_and_iterations() {
        let g = ring_of_cliques();
        let init: Vec<i64> = (0..20).map(|v| (v / 4) as i64).collect();
        let mut o = opts(LeidenObjective::Modularity, 1.0, 2);
        o.initial = Some(&init);
        o.n_iterations = 1;
        let r = leiden(&g, &o).unwrap();
        assert_eq!(r.membership, dense_labels(&init).0);
        assert_eq!(r.iterations, 1);
        o.n_iterations = 0;
        assert_eq!(leiden(&g, &o).unwrap().membership, dense_labels(&init).0);
        // A disconnected starting group is split.
        let g = und(4, &[(0, 1), (2, 3)]);
        let init = [0, 0, 0, 0];
        let mut o = opts(LeidenObjective::Modularity, 1.0, 2);
        o.initial = Some(&init);
        assert_eq!(leiden(&g, &o).unwrap().membership, vec![0, 0, 1, 1]);
    }

    #[test]
    fn edge_cases() {
        let r = leiden(&und(3, &[]), &opts(LeidenObjective::Modularity, 1.0, 1)).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.quality.is_nan());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(leiden(&d, &opts(LeidenObjective::Cpm, 1.0, 1)).is_err());
        let mut o = opts(LeidenObjective::Cpm, 1.0, 1);
        o.beta = 0.0;
        assert!(leiden(&und(2, &[(0, 1)]), &o).is_err());
    }
}
