//! Infomap: communities that minimize the map equation.
//!
//! Rosvall and Bergstrom (2008, `rosvall2008maps`) describe a random walk
//! on the network with a two-level code: a module codebook for each group
//! and an index codebook for movements between groups. With p_a the
//! stationary visit rate of node a, q_i the rate at which the walk exits
//! group i, q = sum_i q_i, and plogp(x) = x log2 x, the description length
//! per step (the map equation, in bits) is
//!
//! ```text
//! L = plogp(q) - 2 sum_i plogp(q_i) - sum_a plogp(p_a)
//!     + sum_i plogp(q_i + sum_{a in i} p_a).
//! ```
//!
//! Flow:
//! - Undirected networks: p_a = s_a / (2m) (s_a the weighted degree, a
//!   self-loop counting twice) and q_i is the weight of the edges leaving
//!   group i over 2m. There is no teleportation, and node weights have no
//!   effect.
//! - Directed networks: the walk follows an out-link (in proportion to its
//!   weight) with probability 1 - tau and teleports with probability
//!   tau = 0.15 (always from a node without out-links), landing on node a
//!   with probability t_a proportional to its node weight (uniform by
//!   default). p is its stationary distribution (power iteration), and
//!   teleportation steps are part of the code, so
//!   q_i = sum_{a in i} p_a [tau_a (1 - T_i) + (1 - tau_a) w_a,out(i) / w_a],
//!   with T_i = sum_{a in i} t_a and w_a,out(i) the weight of a's out-links
//!   leaving i. Equation 9 of Rosvall and Bergstrom (2008) writes the
//!   teleportation term with uniform t as tau (n - n_i)/(n - 1), a jump
//!   that never lands on its own node; here teleportation is the same
//!   jump in the flow and in the code, as in PageRank.
//!
//! Search (Rosvall, Axelsson, and Bergstrom 2009, `rosvall2009map`): the
//! core algorithm visits nodes in random order and moves each to the
//! neighboring group that most decreases L, until no move decreases it,
//! then merges each group into one node and repeats on the merged network
//! until nothing moves. Two refinements then alternate until L stops
//! decreasing: single-node movements (the core algorithm on the original
//! nodes, starting from the current groups) and submodule movements (the
//! core algorithm within each group, on its own internal flow, proposes
//! submodules, which then move between groups as units under the full map
//! equation). The whole search is repeated `trials` times and the
//! partition of shortest description length kept (the first on ties).
//!
//! Conventions:
//! - Directed networks use edge directions (igraph's `cluster_infomap()`
//!   does too); undirected networks have none.
//! - igraph's `cluster_infomap()` gives undirected networks the directed
//!   flow, as if each edge ran in both directions, so it teleports and
//!   reports a different (larger) codelength than the undirected one
//!   here. Its results are reproduced by passing the network as directed
//!   with each edge listed both ways (a self-loop twice); the validation
//!   suite compares the two tools that way.
//! - Weights (if used) are added as given and must be non-negative;
//!   parallel edges add up. Node weights must be finite, non-negative, and
//!   not all zero.
//! - A move must decrease L by more than 1e-10 bits.
//! - Random numbers come from [`crate::rng::Stream`] with the seed (trial
//!   i uses stream i), so a seed reproduces the result on any platform.
//!   The method runs on one thread.
//! - An undirected network without edges has no random walk: every node
//!   is its own group and the codelength is NaN.
//! - Groups are labelled 0, 1, ... in order of their first node.
//!   [`map_equation`] is L of a given partition.

use super::relabel;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::community::dense_labels;
use crate::rng::Stream;

/// Teleportation probability of the directed random walk.
pub const TELEPORT: f64 = 0.15;
const MIN_GAIN: f64 = 1e-10;
const MAX_SWEEPS: usize = 1_000;
const MAX_ROUNDS: usize = 100;

/// Result of [`infomap`].
#[derive(Debug, Clone, PartialEq)]
pub struct Infomap {
    /// Group of each node.
    pub membership: Vec<usize>,
    /// Description length of the partition, in bits per step.
    pub codelength: f64,
}

fn plogp(x: f64) -> f64 {
    if x > 0.0 {
        x * x.log2()
    } else {
        0.0
    }
}

/// A network of nodes (or merged groups) with flows.
#[derive(Clone)]
struct Level {
    /// Visit rate.
    flow: Vec<f64>,
    /// Teleportation flow out of the node: sum of p_a tau_a.
    tele: Vec<f64>,
    /// Teleportation target mass: sum of t_a.
    target: Vec<f64>,
    /// Link flows to other nodes, and from other nodes.
    out: Vec<Vec<(usize, f64)>>,
    inn: Vec<Vec<(usize, f64)>>,
}

impl Level {
    fn n(&self) -> usize {
        self.flow.len()
    }

    /// Merges the groups of `part` (k groups) into nodes.
    fn aggregate(&self, part: &[usize], k: usize) -> Level {
        let mut flow = vec![0.0; k];
        let mut tele = vec![0.0; k];
        let mut target = vec![0.0; k];
        let mut out: Vec<Vec<(usize, f64)>> = vec![Vec::new(); k];
        for v in 0..self.n() {
            let c = part[v];
            flow[c] += self.flow[v];
            tele[c] += self.tele[v];
            target[c] += self.target[v];
            for &(u, f) in &self.out[v] {
                if part[u] != c {
                    out[c].push((part[u], f));
                }
            }
        }
        let out = merge_rows(out);
        let inn = transpose(&out, k);
        Level {
            flow,
            tele,
            target,
            out,
            inn,
        }
    }

    /// The nodes in `nodes` with only the links among them.
    fn subnetwork(&self, nodes: &[usize], index: &[usize]) -> Level {
        let k = nodes.len();
        let mut out: Vec<Vec<(usize, f64)>> = vec![Vec::new(); k];
        for (i, &v) in nodes.iter().enumerate() {
            for &(u, f) in &self.out[v] {
                if index[u] != usize::MAX {
                    out[i].push((index[u], f));
                }
            }
        }
        let inn = transpose(&out, k);
        Level {
            flow: nodes.iter().map(|&v| self.flow[v]).collect(),
            tele: nodes.iter().map(|&v| self.tele[v]).collect(),
            target: nodes.iter().map(|&v| self.target[v]).collect(),
            out,
            inn,
        }
    }
}

fn merge_rows(mut rows: Vec<Vec<(usize, f64)>>) -> Vec<Vec<(usize, f64)>> {
    for row in rows.iter_mut() {
        row.sort_unstable_by_key(|&(u, _)| u);
        let mut merged: Vec<(usize, f64)> = Vec::with_capacity(row.len());
        for &(u, f) in row.iter() {
            match merged.last_mut() {
                Some(last) if last.0 == u => last.1 += f,
                _ => merged.push((u, f)),
            }
        }
        *row = merged;
    }
    rows
}

fn transpose(out: &[Vec<(usize, f64)>], k: usize) -> Vec<Vec<(usize, f64)>> {
    let mut inn: Vec<Vec<(usize, f64)>> = vec![Vec::new(); k];
    for (v, row) in out.iter().enumerate() {
        for &(u, f) in row {
            inn[u].push((v, f));
        }
    }
    inn
}

/// The flow network of the input, and the node-level entropy term
/// sum_a plogp(p_a). `None` when an undirected network has no edges.
fn flow_network(g: &Graph, weighted: bool, node_weights: Option<&[f64]>) -> Result<Option<Level>> {
    let n = g.node_count();
    let w = crate::weights::additive_weights(g, weighted, "Infomap")?;
    let mut target = match node_weights {
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
            if n > 0 && x.iter().sum::<f64>() <= 0.0 {
                return Err(Error::InvalidArgument(
                    "node weights must not all be zero".into(),
                ));
            }
            x.to_vec()
        }
        None => vec![1.0; n],
    };
    let total_t: f64 = target.iter().sum();
    target.iter_mut().for_each(|t| *t /= total_t);
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut strength = vec![0.0f64; n];
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        strength[u] += x;
        rows[u].push((v, x));
        if !g.is_directed() {
            strength[v] += x;
            if u != v {
                rows[v].push((u, x));
            } else {
                // A loop counts twice toward an undirected degree.
                rows[u].push((u, x));
            }
        }
    }
    let rows = merge_rows(rows);
    let (flow, tele, out) = if g.is_directed() {
        let dangling: Vec<bool> = strength.iter().map(|&s| s <= 0.0).collect();
        let mut p = vec![1.0 / n as f64; n];
        let mut next = vec![0.0f64; n];
        for _ in 0..10_000 {
            next.iter_mut().for_each(|x| *x = 0.0);
            let mut jump = 0.0;
            for a in 0..n {
                if dangling[a] {
                    jump += p[a];
                } else {
                    jump += TELEPORT * p[a];
                    let s = (1.0 - TELEPORT) * p[a] / strength[a];
                    for &(b, x) in &rows[a] {
                        next[b] += s * x;
                    }
                }
            }
            for (b, x) in next.iter_mut().enumerate() {
                *x += jump * target[b];
            }
            let total: f64 = next.iter().sum();
            next.iter_mut().for_each(|x| *x /= total);
            let diff: f64 = p.iter().zip(&next).map(|(a, b)| (a - b).abs()).sum();
            std::mem::swap(&mut p, &mut next);
            if diff < 1e-15 {
                break;
            }
        }
        let tele: Vec<f64> = (0..n)
            .map(|a| if dangling[a] { p[a] } else { TELEPORT * p[a] })
            .collect();
        let out: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|a| {
                rows[a]
                    .iter()
                    .filter(|&&(b, _)| b != a)
                    .map(|&(b, x)| (b, (1.0 - TELEPORT) * p[a] * x / strength[a]))
                    .collect()
            })
            .collect();
        (p, tele, out)
    } else {
        let two_m: f64 = strength.iter().sum();
        if two_m <= 0.0 {
            return Ok(None);
        }
        let p: Vec<f64> = strength.iter().map(|s| s / two_m).collect();
        let out: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|a| {
                rows[a]
                    .iter()
                    .filter(|&&(b, _)| b != a)
                    .map(|&(b, x)| (b, x / two_m))
                    .collect()
            })
            .collect();
        (p, vec![0.0; n], out)
    };
    let inn = transpose(&out, n);
    Ok(Some(Level {
        flow,
        tele,
        target,
        out,
        inn,
    }))
}

/// Group totals of a partition of a level, and the terms of L that depend
/// on the partition.
struct State {
    comm: Vec<usize>,
    flow: Vec<f64>,
    tele: Vec<f64>,
    target: Vec<f64>,
    /// Link flow leaving each group.
    link: Vec<f64>,
    sum_exit: f64,
    sum_plogp_exit: f64,
    sum_plogp_exit_flow: f64,
}

impl State {
    fn exit(tele: f64, target: f64, link: f64) -> f64 {
        tele * (1.0 - target) + link
    }

    fn new(lv: &Level, comm: Vec<usize>) -> State {
        let n = lv.n();
        let mut s = State {
            comm,
            flow: vec![0.0; n],
            tele: vec![0.0; n],
            target: vec![0.0; n],
            link: vec![0.0; n],
            sum_exit: 0.0,
            sum_plogp_exit: 0.0,
            sum_plogp_exit_flow: 0.0,
        };
        for v in 0..n {
            let c = s.comm[v];
            s.flow[c] += lv.flow[v];
            s.tele[c] += lv.tele[v];
            s.target[c] += lv.target[v];
            for &(u, f) in &lv.out[v] {
                if s.comm[u] != c {
                    s.link[c] += f;
                }
            }
        }
        for c in 0..n {
            let e = State::exit(s.tele[c], s.target[c], s.link[c]);
            s.sum_exit += e;
            s.sum_plogp_exit += plogp(e);
            s.sum_plogp_exit_flow += plogp(e + s.flow[c]);
        }
        s
    }

    /// L without the constant node term.
    fn codelength(&self) -> f64 {
        plogp(self.sum_exit) - 2.0 * self.sum_plogp_exit + self.sum_plogp_exit_flow
    }
}

/// One node's flows to and from each neighbouring group.
struct Gather {
    to: Vec<f64>,
    from: Vec<f64>,
    seen: Vec<bool>,
    list: Vec<usize>,
}

impl Gather {
    fn new(n: usize) -> Gather {
        Gather {
            to: vec![0.0; n],
            from: vec![0.0; n],
            seen: vec![false; n],
            list: Vec::new(),
        }
    }

    fn touch(&mut self, c: usize) {
        if !self.seen[c] {
            self.seen[c] = true;
            self.list.push(c);
        }
    }

    fn clear(&mut self) {
        for &c in &self.list {
            self.to[c] = 0.0;
            self.from[c] = 0.0;
            self.seen[c] = false;
        }
        self.list.clear();
    }
}

/// Core local moves on one level until no move decreases L. Returns
/// whether any node moved.
fn local_moves(lv: &Level, st: &mut State, rng: &mut Stream) -> bool {
    let n = lv.n();
    let mut g = Gather::new(n);
    let mut order: Vec<usize> = (0..n).collect();
    let mut moved_any = false;
    for _ in 0..MAX_SWEEPS {
        rng.permutation(&mut order);
        let mut moved = false;
        for &v in &order {
            let a = st.comm[v];
            let mut out_v = 0.0;
            for &(u, f) in &lv.out[v] {
                out_v += f;
                let c = st.comm[u];
                g.touch(c);
                g.to[c] += f;
            }
            for &(u, f) in &lv.inn[v] {
                let c = st.comm[u];
                g.touch(c);
                g.from[c] += f;
            }
            let (pv, tv, sv) = (lv.flow[v], lv.tele[v], lv.target[v]);
            // Group a without v.
            let a_flow = st.flow[a] - pv;
            let a_tele = st.tele[a] - tv;
            let a_target = st.target[a] - sv;
            let a_link = st.link[a] - (out_v - g.to[a]) + g.from[a];
            let old_a = State::exit(st.tele[a], st.target[a], st.link[a]);
            let new_a = State::exit(a_tele, a_target, a_link);
            let mut best = a;
            let mut best_delta = -MIN_GAIN;
            let mut best_terms = (0.0, 0.0, 0.0);
            for &b in &g.list {
                if b == a {
                    continue;
                }
                let b_link = st.link[b] + (out_v - g.to[b]) - g.from[b];
                let old_b = State::exit(st.tele[b], st.target[b], st.link[b]);
                let new_b = State::exit(st.tele[b] + tv, st.target[b] + sv, b_link);
                let sum_exit = st.sum_exit - old_a - old_b + new_a + new_b;
                let sum_pe =
                    st.sum_plogp_exit - plogp(old_a) - plogp(old_b) + plogp(new_a) + plogp(new_b);
                let sum_pef =
                    st.sum_plogp_exit_flow - plogp(old_a + st.flow[a]) - plogp(old_b + st.flow[b])
                        + plogp(new_a + a_flow)
                        + plogp(new_b + st.flow[b] + pv);
                let delta = plogp(sum_exit) - 2.0 * sum_pe + sum_pef - st.codelength();
                if delta < best_delta {
                    best_delta = delta;
                    best = b;
                    best_terms = (sum_exit, sum_pe, sum_pef);
                }
            }
            if best != a {
                let b = best;
                st.link[b] += (out_v - g.to[b]) - g.from[b];
                st.link[a] = a_link;
                st.flow[a] = a_flow;
                st.tele[a] = a_tele;
                st.target[a] = a_target;
                st.flow[b] += pv;
                st.tele[b] += tv;
                st.target[b] += sv;
                st.sum_exit = best_terms.0;
                st.sum_plogp_exit = best_terms.1;
                st.sum_plogp_exit_flow = best_terms.2;
                st.comm[v] = b;
                moved = true;
                moved_any = true;
            }
            g.clear();
        }
        if !moved {
            break;
        }
    }
    moved_any
}

/// The core algorithm from a starting partition of `lv`: local moves, then
/// merging, until nothing moves. Returns the partition of `lv`'s nodes.
fn optimize(lv: &Level, init: Vec<usize>, rng: &mut Stream) -> Vec<usize> {
    let mut map: Vec<usize> = (0..lv.n()).collect();
    let mut level = lv.clone();
    let mut comm = dense_labels(&init).0;
    loop {
        let mut st = State::new(&level, comm);
        local_moves(&level, &mut st, rng);
        let (dense, k) = dense_labels(&st.comm);
        if k == level.n() {
            return map.iter().map(|&v| dense[v]).collect();
        }
        for m in map.iter_mut() {
            *m = dense[*m];
        }
        level = level.aggregate(&dense, k);
        comm = (0..k).collect();
    }
}

/// Description length of a partition of a level (without the node term).
fn partition_codelength(lv: &Level, comm: &[usize]) -> f64 {
    State::new(lv, dense_labels(comm).0).codelength()
}

/// One trial: core algorithm, then alternating single-node and submodule
/// movements while L decreases.
fn trial(lv: &Level, rng: &mut Stream) -> (Vec<usize>, f64) {
    let n = lv.n();
    let mut comm = optimize(lv, (0..n).collect(), rng);
    let mut len = partition_codelength(lv, &comm);
    for _ in 0..MAX_ROUNDS {
        let before = len;
        // Single-node movements.
        let fine = optimize(lv, comm.clone(), rng);
        let fine_len = partition_codelength(lv, &fine);
        if fine_len < len {
            comm = fine;
            len = fine_len;
        }
        // Submodule movements.
        let (dense, k) = dense_labels(&comm);
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); k];
        for v in 0..n {
            members[dense[v]].push(v);
        }
        let mut sub = vec![0usize; n];
        let mut index = vec![usize::MAX; n];
        let mut next = 0;
        for group in &members {
            for (i, &v) in group.iter().enumerate() {
                index[v] = i;
            }
            let part = optimize(
                &lv.subnetwork(group, &index),
                (0..group.len()).collect(),
                rng,
            );
            let (part, kp) = dense_labels(&part);
            for (i, &v) in group.iter().enumerate() {
                sub[v] = next + part[i];
                index[v] = usize::MAX;
            }
            next += kp;
        }
        let coarse_level = lv.aggregate(&sub, next);
        let mut start = vec![0usize; next];
        for v in 0..n {
            start[sub[v]] = dense[v];
        }
        let moved = optimize(&coarse_level, start, rng);
        let coarse: Vec<usize> = (0..n).map(|v| moved[sub[v]]).collect();
        let coarse_len = partition_codelength(lv, &coarse);
        if coarse_len < len {
            comm = coarse;
            len = coarse_len;
        }
        if before - len <= MIN_GAIN {
            break;
        }
    }
    (comm, len)
}

/// Infomap communities: the best of `trials` searches. See the module
/// documentation.
pub fn infomap(
    g: &Graph,
    weighted: bool,
    node_weights: Option<&[f64]>,
    trials: usize,
    seed: u64,
) -> Result<Infomap> {
    if trials == 0 {
        return Err(Error::InvalidArgument(
            "the number of trials must be at least 1".into(),
        ));
    }
    let n = g.node_count();
    let lv = match flow_network(g, weighted, node_weights)? {
        Some(lv) => lv,
        None => {
            return Ok(Infomap {
                membership: (0..n).collect(),
                codelength: f64::NAN,
            })
        }
    };
    let node_term: f64 = lv.flow.iter().map(|&p| plogp(p)).sum();
    let mut best: Option<(Vec<usize>, f64)> = None;
    for t in 0..trials {
        let mut rng = Stream::new(seed, t as u64);
        let (comm, len) = trial(&lv, &mut rng);
        if best.as_ref().map_or(true, |b| len < b.1) {
            best = Some((comm, len));
        }
    }
    let (comm, len) = best.unwrap();
    Ok(Infomap {
        membership: relabel(&comm),
        codelength: len - node_term,
    })
}

/// The map equation L (bits per step) of a partition (one label per node;
/// labels arbitrary non-negative integers). NaN for an undirected network
/// without edges.
pub fn map_equation(
    g: &Graph,
    membership: &[i64],
    weighted: bool,
    node_weights: Option<&[f64]>,
) -> Result<f64> {
    let n = g.node_count();
    if membership.len() != n {
        return Err(Error::LengthMismatch {
            what: "membership",
            expected: n,
            found: membership.len(),
        });
    }
    if let Some(&bad) = membership.iter().find(|&&c| c < 0) {
        return Err(Error::InvalidArgument(format!(
            "membership labels must be non-negative, found {bad}"
        )));
    }
    let lv = match flow_network(g, weighted, node_weights)? {
        Some(lv) => lv,
        None => return Ok(f64::NAN),
    };
    let node_term: f64 = lv.flow.iter().map(|&p| plogp(p)).sum();
    Ok(partition_codelength(&lv, &dense_labels(membership).0) - node_term)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn map_equation_by_hand() {
        // Two triangles joined by one edge, split in two: 2m = 14,
        // p = (2,2,3,3,2,2)/14, q_1 = q_2 = 1/14.
        let g = und(6, &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)]);
        let l = map_equation(&g, &[0, 0, 0, 1, 1, 1], false, None).unwrap();
        let p = [2.0, 2.0, 3.0, 3.0, 2.0, 2.0].map(|x| x / 14.0);
        let q = 1.0 / 14.0;
        let expect =
            plogp(2.0 * q) - 2.0 * 2.0 * plogp(q) - p.iter().map(|&x| plogp(x)).sum::<f64>()
                + 2.0 * plogp(q + 0.5);
        assert!((l - expect).abs() < 1e-12);
        // One group: the entropy of p.
        let one = map_equation(&g, &[0; 6], false, None).unwrap();
        let h = -p.iter().map(|&x| plogp(x)).sum::<f64>();
        assert!((one - h).abs() < 1e-12);
    }

    #[test]
    fn finds_ring_of_cliques() {
        let mut e = Vec::new();
        for c in 0..5i64 {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    e.push((4 * c + i, 4 * c + j));
                }
            }
            e.push((4 * c, (4 * c + 4) % 20));
        }
        let g = und(20, &e);
        let r = infomap(&g, false, None, 3, 1).unwrap();
        let expect: Vec<usize> = (0..20).map(|v| v / 4).collect();
        assert_eq!(r.membership, expect);
        let memb: Vec<i64> = expect.iter().map(|&c| c as i64).collect();
        let l = map_equation(&g, &memb, false, None).unwrap();
        assert!((r.codelength - l).abs() < 1e-12);
        assert_eq!(r, infomap(&g, false, None, 3, 1).unwrap());
    }

    #[test]
    fn directed_flow_sums_to_one() {
        // Two directed 3-cycles joined one way.
        let s = [0, 1, 2, 3, 4, 5, 2];
        let t = [1, 2, 0, 4, 5, 3, 3];
        let g = Graph::from_edges(6, &s, &t, None, true).unwrap();
        let lv = flow_network(&g, false, None).unwrap().unwrap();
        assert!((lv.flow.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let r = infomap(&g, false, None, 5, 2).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        let l = map_equation(&g, &[0, 0, 0, 1, 1, 1], false, None).unwrap();
        assert!((r.codelength - l).abs() < 1e-12);
    }

    #[test]
    fn edge_cases() {
        let r = infomap(&und(3, &[]), false, None, 1, 1).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2]);
        assert!(r.codelength.is_nan());
        assert!(infomap(&und(2, &[(0, 1)]), false, None, 0, 1).is_err());
        assert!(infomap(&und(2, &[(0, 1)]), false, Some(&[0.0, 0.0]), 1, 1).is_err());
    }
}
