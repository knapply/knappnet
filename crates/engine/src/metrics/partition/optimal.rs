//! Exact maximum modularity for small networks.
//!
//! Finding the partition of maximum modularity is NP-hard (Brandes,
//! Delling, Gaertler, Goerke, Hoefer, Nikoloski, and Wagner 2008,
//! `brandes2008modularity`), who also give the integer linear program that
//! igraph's `cluster_optimal()` solves with GLPK. Without a solver
//! dependency, this module solves the same problem exactly by branch and
//! bound over assignments of nodes to groups.
//!
//! With B_ij = A_ij + A_ji - gamma (k_i^out k_j^in + k_j^out k_i^in) / m'
//! for each pair i < j (undirected: m' = 2m and k^out = k^in = k with
//! A symmetric, so B_ij = 2 A_ij - gamma k_i k_j / m), modularity is a
//! constant plus (sum of B_ij over pairs in the same group) / m', so the
//! search maximizes F = sum of B_ij over pairs placed together. Facts used:
//! - Groups of an optimal partition are (weakly) connected, since B_ij <= 0
//!   for non-adjacent pairs: splitting a disconnected group never lowers
//!   F. Groups therefore never span components, F is a sum over
//!   components, and each component is solved on its own.
//! - Nodes are assigned in maximum-adjacency order (each next node the one
//!   most strongly tied to those already placed), each to an existing
//!   group or to a new one (groups opened in order, which removes the
//!   symmetry between group labels), best-scoring choices first.
//! - Bound: with the placed nodes fixed, each unplaced node can gain at
//!   most max(0, max over groups of its total B to the group's placed
//!   members), and each pair of unplaced nodes at most max(0, B_ij); a
//!   branch whose value plus these gains cannot beat the best partition
//!   found so far is cut. The starting incumbent is the better of
//!   singletons and one group.
//!
//! With whole-number weights and gamma = 1 every quantity is an integer
//! multiple of 1/m' and the search compares them exactly.
//!
//! Conventions:
//! - Directed and undirected networks, with modularity as
//!   [`crate::metrics::community::modularity`] defines it (igraph's
//!   `cluster_optimal()` also accepts both); weights (if used) must be
//!   non-negative, parallel edges add up, and self-loops only add a
//!   constant.
//! - Size limit: at most 60 nodes in any component, and at most 10^9
//!   branch-and-bound steps in all; beyond either the function returns an
//!   error rather than a partition that might not be optimal. Sparse
//!   networks of 30 to 40 nodes typically take well under a second;
//!   dense ones take much longer.
//! - Ties between optimal partitions go to the first found. Nodes of
//!   degree zero are groups of their own. A network without edges has
//!   modularity NaN and one group per node.
//! - Groups are labelled 0, 1, ... in order of their first node.

use super::relabel;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::community::modularity;

/// Largest component the search accepts.
pub const MAX_COMPONENT: usize = 60;
const MAX_STEPS: u64 = 1_000_000_000;

/// Result of [`optimal_modularity`].
#[derive(Debug, Clone, PartialEq)]
pub struct OptimalModularity {
    /// Group of each node in a partition of maximum modularity.
    pub membership: Vec<usize>,
    /// Its modularity.
    pub modularity: f64,
}

struct Search<'a> {
    /// Pair scores, dense k x k, in search order.
    b: &'a [f64],
    k: usize,
    /// suffix[i]: sum of positive B over pairs among nodes i.. (search order).
    suffix: Vec<f64>,
    group: Vec<usize>,
    /// gain[v * k + c]: total B from v to the placed members of group c.
    gain: Vec<f64>,
    groups: usize,
    best: f64,
    best_group: Vec<usize>,
    steps: u64,
}

impl Search<'_> {
    fn bound(&self, depth: usize) -> f64 {
        let mut ub = self.suffix[depth];
        for v in depth..self.k {
            let row = &self.gain[v * self.k..v * self.k + self.groups];
            ub += row.iter().copied().fold(0.0, f64::max);
        }
        ub
    }

    fn place(&mut self, v: usize, c: usize, sign: f64) {
        for u in (v + 1)..self.k {
            self.gain[u * self.k + c] += sign * self.b[v * self.k + u];
        }
    }

    fn run(&mut self, depth: usize, value: f64) -> Result<()> {
        self.steps += 1;
        if self.steps > MAX_STEPS {
            return Err(Error::InvalidArgument(format!(
                "the exact search for maximum modularity exceeded {MAX_STEPS} steps; the \
                 network is too large or too dense for an exact solution"
            )));
        }
        if depth == self.k {
            if value > self.best {
                self.best = value;
                self.best_group.clone_from(&self.group);
            }
            return Ok(());
        }
        if value + self.bound(depth) <= self.best {
            return Ok(());
        }
        let v = depth;
        let mut choices: Vec<(f64, usize)> = (0..self.groups)
            .map(|c| (self.gain[v * self.k + c], c))
            .collect();
        choices.push((0.0, self.groups));
        choices.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        for (g, c) in choices {
            let opened = c == self.groups;
            if opened {
                self.groups += 1;
            }
            self.group[v] = c;
            self.place(v, c, 1.0);
            let r = self.run(depth + 1, value + g);
            self.place(v, c, -1.0);
            if opened {
                self.groups -= 1;
            }
            r?;
        }
        Ok(())
    }
}

/// A partition of maximum modularity (exact). See the module documentation.
pub fn optimal_modularity(g: &Graph, resolution: f64, weighted: bool) -> Result<OptimalModularity> {
    if !(resolution.is_finite() && resolution >= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "resolution must be finite and non-negative, not {resolution}"
        )));
    }
    let n = g.node_count();
    let w = crate::weights::additive_weights(g, weighted, "modularity")?;
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut kout = vec![0.0f64; n];
    let mut kin = vec![0.0f64; n];
    let mut pair: std::collections::HashMap<(usize, usize), f64> = Default::default();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        kout[u] += x;
        kin[v] += x;
        if !g.is_directed() {
            kout[v] += x;
            kin[u] += x;
        }
        if u != v {
            *pair.entry((u.min(v), u.max(v))).or_insert(0.0) += x;
            adj[u].push(v);
            adj[v].push(u);
        }
    }
    // m' (2m undirected, m directed).
    let mp: f64 = kout.iter().sum();
    if mp == 0.0 {
        return Ok(OptimalModularity {
            membership: (0..n).collect(),
            modularity: f64::NAN,
        });
    }
    // Pair score, a positive multiple of the pair's contribution to
    // modularity: m' c (A_ij + A_ji) - gamma (k_i^out k_j^in + k_j^out k_i^in),
    // with c = 2 for undirected networks (where `pair` holds A_ij once)
    // and c = 1 for directed ones (where it holds A_ij + A_ji).
    let c = if g.is_directed() { 1.0 } else { 2.0 };
    let score = |i: usize, j: usize| -> f64 {
        let a = pair.get(&(i.min(j), i.max(j))).copied().unwrap_or(0.0);
        mp * c * a - resolution * (kout[i] * kin[j] + kout[j] * kin[i])
    };
    // Components of the underlying undirected network.
    let mut comp = vec![usize::MAX; n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for s in 0..n {
        if comp[s] != usize::MAX {
            continue;
        }
        let id = comps.len();
        comp[s] = id;
        let mut members = vec![s];
        let mut i = 0;
        while i < members.len() {
            let v = members[i];
            i += 1;
            for &u in &adj[v] {
                if comp[u] == usize::MAX {
                    comp[u] = id;
                    members.push(u);
                }
            }
        }
        comps.push(members);
    }
    if let Some(big) = comps.iter().map(|c| c.len()).max() {
        if big > MAX_COMPONENT {
            return Err(Error::InvalidArgument(format!(
                "exact maximum modularity is limited to components of at most \
                 {MAX_COMPONENT} nodes; this network has a component of {big}"
            )));
        }
    }
    let mut label = vec![0usize; n];
    let mut next = 0usize;
    let mut steps = 0u64;
    for members in &comps {
        if members.len() == 1 {
            label[members[0]] = next;
            next += 1;
            continue;
        }
        // Maximum-adjacency order within the component.
        let k = members.len();
        let full: Vec<Vec<f64>> = members
            .iter()
            .map(|&i| {
                members
                    .iter()
                    .map(|&j| if i == j { 0.0 } else { score(i, j) })
                    .collect()
            })
            .collect();
        let tie: Vec<Vec<f64>> = members
            .iter()
            .map(|&i| {
                members
                    .iter()
                    .map(|&j| pair.get(&(i.min(j), i.max(j))).copied().unwrap_or(0.0))
                    .collect()
            })
            .collect();
        let start = (0..k)
            .max_by(|&a, &b| {
                let (da, db) = (
                    kout[members[a]] + kin[members[a]],
                    kout[members[b]] + kin[members[b]],
                );
                da.total_cmp(&db).then(b.cmp(&a))
            })
            .unwrap();
        let mut order = vec![start];
        let mut placed = vec![false; k];
        placed[start] = true;
        let mut attach: Vec<f64> = tie[start].clone();
        while order.len() < k {
            let v = (0..k)
                .filter(|&v| !placed[v])
                .max_by(|&a, &b| attach[a].total_cmp(&attach[b]).then(b.cmp(&a)))
                .unwrap();
            placed[v] = true;
            order.push(v);
            for u in 0..k {
                attach[u] += tie[v][u];
            }
        }
        let b: Vec<f64> = order
            .iter()
            .flat_map(|&i| order.iter().map(move |&j| (i, j)))
            .map(|(i, j)| full[i][j])
            .collect();
        let mut suffix = vec![0.0f64; k + 1];
        for i in (0..k).rev() {
            let row: f64 = ((i + 1)..k).map(|j| b[i * k + j].max(0.0)).sum();
            suffix[i] = suffix[i + 1] + row;
        }
        // Incumbent: the better of one group and singletons.
        let together: f64 = (0..k)
            .flat_map(|i| ((i + 1)..k).map(move |j| (i, j)))
            .map(|(i, j)| b[i * k + j])
            .sum();
        let (best, best_group) = if together > 0.0 {
            (together, vec![0; k])
        } else {
            (0.0, (0..k).collect())
        };
        let mut s = Search {
            b: &b,
            k,
            suffix,
            group: vec![0; k],
            gain: vec![0.0; k * k],
            groups: 0,
            best,
            best_group,
            steps,
        };
        s.run(0, 0.0)?;
        steps = s.steps;
        let (groups, _) = crate::metrics::community::dense_labels(&s.best_group);
        let kg = groups.iter().copied().max().unwrap_or(0) + 1;
        for (pos, &i) in order.iter().enumerate() {
            label[members[i]] = next + groups[pos];
        }
        next += kg;
    }
    let membership = relabel(&label);
    let memb: Vec<i64> = membership.iter().map(|&c| c as i64).collect();
    let q = modularity(g, &memb, resolution, weighted)?;
    Ok(OptimalModularity {
        membership,
        modularity: q,
    })
}

/// Every partition of `n` nodes as restricted growth strings, for tests
/// and validation of tiny networks.
#[cfg(test)]
fn partitions(n: usize, f: &mut impl FnMut(&[usize])) {
    fn rec(a: &mut Vec<usize>, n: usize, max: usize, f: &mut impl FnMut(&[usize])) {
        if a.len() == n {
            f(a);
            return;
        }
        for c in 0..=max + 1 {
            a.push(c);
            rec(a, n, max.max(c), f);
            a.pop();
        }
    }
    if n == 0 {
        f(&[]);
        return;
    }
    let mut a = vec![0];
    rec(&mut a, n, 0, f);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    fn brute(g: &Graph, gamma: f64, weighted: bool) -> f64 {
        let mut best = f64::NEG_INFINITY;
        partitions(g.node_count(), &mut |p| {
            let m: Vec<i64> = p.iter().map(|&c| c as i64).collect();
            best = best.max(modularity(g, &m, gamma, weighted).unwrap());
        });
        best
    }

    #[test]
    fn matches_enumeration_on_random_networks() {
        for trial in 0..40u64 {
            let mut rng = Stream::new(99, trial);
            let n = 3 + rng.below(6) as usize;
            let directed = trial % 3 == 0;
            let (mut s, mut t, mut w) = (Vec::new(), Vec::new(), Vec::new());
            for i in 0..n {
                for j in 0..n {
                    if (directed || i < j) && i != j && rng.uniform() < 0.35 {
                        s.push(i as i64);
                        t.push(j as i64);
                        w.push(1.0 + rng.below(3) as f64);
                    }
                }
            }
            if s.is_empty() {
                continue;
            }
            let g = Graph::from_edges(n, &s, &t, Some(&w), directed).unwrap();
            for (gamma, weighted) in [(1.0, false), (1.0, true), (0.5, false), (1.7, true)] {
                let r = optimal_modularity(&g, gamma, weighted).unwrap();
                let b = brute(&g, gamma, weighted);
                assert!(
                    (r.modularity - b).abs() < 1e-12,
                    "{trial} {gamma} {} {b}",
                    r.modularity
                );
            }
        }
    }

    #[test]
    fn edge_cases() {
        let g = Graph::from_edges(3, &[], &[], None, false).unwrap();
        let r = optimal_modularity(&g, 1.0, false).unwrap();
        assert!(r.modularity.is_nan());
        assert_eq!(r.membership, vec![0, 1, 2]);
        let (s, t): (Vec<i64>, Vec<i64>) = (0..61).map(|i| (i, i + 1)).unzip();
        let path = Graph::from_edges(62, &s, &t, None, false).unwrap();
        assert!(optimal_modularity(&path, 1.0, false).is_err());
    }
}
