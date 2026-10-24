//! Node similarity from walks: SimRank and Panther.
//!
//! # SimRank
//!
//! SimRank (Jeh and Widom 2002, `jeh2002simrank`) scores two nodes as
//! similar when their in-neighbours are similar: s(u, u) = 1 and, for
//! u != v,
//!
//! s(u, v) = c / (|I(u)| |I(v)|) Σ_{a ∈ I(u)} Σ_{b ∈ I(v)} s(a, b),
//!
//! with I(u) the in-neighbours of u (the neighbours in an undirected
//! graph) and 0 when either set is empty. [`simrank`] iterates the matrix
//! form S ← c Wᵀ S W with the diagonal reset to 1, W the adjacency matrix
//! with each column divided by its sum (Jacobi iteration from S = I, the
//! "naive" method of the paper), with two stopping rules, those of
//! NetworkX 3.7's two implementations (`hagberg2008exploring`):
//!
//! - [`SimRankStop::Allclose`] (`simrank_similarity`, the public
//!   function): stop after the first iteration with |S_new − S_old| ≤
//!   tol + 1e-5 |S_new| for every entry (`numpy.allclose`). W holds the
//!   edge weights when `weighted` (NetworkX reads the `weight` attribute
//!   through `to_numpy_array`), else 1 per adjacent pair.
//! - [`SimRankStop::Relative`] (`_simrank_similarity_python`, the
//!   pairwise definition): stop after the first iteration with |S_new −
//!   S_old| ≤ tol (1 + |S_old|) for every entry; unweighted.
//!
//! As in NetworkX, reaching `max_iterations` iterations is an error, even
//! when the last one met the stopping rule.
//!
//! Conventions: direction is respected (in-neighbours); parallel edges
//! are merged without weights and their weights added with them; a
//! self-loop makes a node its own in-neighbour (weight counted once, as
//! the diagonal of NetworkX's adjacency matrix). Each entry is computed in
//! a fixed order, so results do not depend on the thread count.
//!
//! # Panther
//!
//! Panther (Zhang, Tang, Ma, Tong, Jing, and Li 2015, `zhang2015panther`)
//! estimates the similarity of u and v as the share of R random paths
//! that contain both: each path starts at a node drawn uniformly and takes
//! T steps of a random walk (each step to an out-neighbour with
//! probability proportional to the edge weight, or uniformly without
//! weights), with R = ⌊c / ε² (log₂ C(T, 2) + 1 + ln(1/δ))⌋ paths (the
//! paper's bound; ε defaults to 1/√m). [`panther_similarity`] returns the
//! k nodes most similar to a source, as NetworkX 3.7 does.
//!
//! Conventions, as in NetworkX's `panther_similarity`:
//! - Isolated nodes (no edge at all; a self-loop is an edge) are removed
//!   before sampling; a source that is isolated is an error. m counts the
//!   edges of the network without them (parallel edges once).
//! - Every node must have an out-neighbour (in directed networks a walk
//!   cannot continue from a node without one; NetworkX fails there too).
//! - Weights must be non-negative; a node's weights may not all be zero.
//! - Paths are sampled from the counter-based stream (`seed`, path index),
//!   so results are reproducible for a seed on every machine and thread
//!   count; they are not NetworkX's samples, which come from NumPy's
//!   generator. Scores agree with NetworkX's in distribution.
//! - Ties in the ranking (NetworkX's order is unspecified there) are
//!   broken by putting the source first, then by node index; the source
//!   is removed from the result, which therefore has k entries, or k + 1
//!   when more than k + 1 nodes tie with the source.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::matrix::Dense;
use crate::parallel;
use crate::rng::Stream;

// ---- Weighted adjacency -----------------------------------------------------

/// Distinct neighbours with total weights, rows sorted: out-rows and
/// in-rows (the same for undirected graphs).
struct Weighted {
    out: Vec<Vec<(u32, f64)>>,
    inc: Vec<Vec<(u32, f64)>>,
}

impl Weighted {
    fn new(g: &Graph, weighted: bool) -> Result<Weighted> {
        let n = g.node_count();
        let w = if weighted {
            Some(g.edge_weights().ok_or(Error::MissingWeights)?)
        } else {
            None
        };
        let mut out: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
        let mut inc: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
        for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
            let x = w.map_or(1.0, |w| w[e]);
            out[s as usize].push((t, x));
            inc[t as usize].push((s, x));
            if !g.is_directed() && s != t {
                out[t as usize].push((s, x));
                inc[s as usize].push((t, x));
            }
        }
        for rows in [&mut out, &mut inc] {
            for r in rows.iter_mut() {
                r.sort_by_key(|p| p.0);
                let mut merged: Vec<(u32, f64)> = Vec::with_capacity(r.len());
                for &(v, x) in r.iter() {
                    match merged.last_mut() {
                        Some(last) if last.0 == v => {
                            if weighted {
                                last.1 += x
                            }
                        }
                        _ => merged.push((v, x)),
                    }
                }
                *r = merged;
            }
        }
        Ok(Weighted { out, inc })
    }
}

// ---- SimRank -----------------------------------------------------------------

/// SimRank's stopping rule (see the module documentation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimRankStop {
    /// NetworkX's `simrank_similarity` (`numpy.allclose`).
    Allclose,
    /// NetworkX's pairwise implementation (relative change).
    Relative,
}

impl SimRankStop {
    /// Decodes the C interface codes: 0 = allclose, 1 = relative.
    pub fn from_code(code: i32) -> Result<SimRankStop> {
        match code {
            0 => Ok(SimRankStop::Allclose),
            1 => Ok(SimRankStop::Relative),
            _ => Err(Error::InvalidArgument(format!(
                "unknown SimRank stopping rule code {code}"
            ))),
        }
    }
}

/// Options for [`simrank`].
#[derive(Debug, Clone, Copy)]
pub struct SimRankOptions {
    /// Decay factor c in (0, 1) (NetworkX's `importance_factor`, 0.9).
    pub importance: f64,
    /// Iteration limit (1000); reaching it is an error.
    pub max_iterations: usize,
    /// Tolerance of the stopping rule (1e-4).
    pub tolerance: f64,
    /// Use edge weights in W (only with [`SimRankStop::Allclose`]).
    pub weighted: bool,
    /// Stopping rule.
    pub stop: SimRankStop,
    /// Threads (0 = all cores).
    pub threads: usize,
}

impl Default for SimRankOptions {
    fn default() -> Self {
        SimRankOptions {
            importance: 0.9,
            max_iterations: 1000,
            tolerance: 1e-4,
            weighted: false,
            stop: SimRankStop::Allclose,
            threads: 0,
        }
    }
}

/// The SimRank matrix (row-major, n by n). See the module documentation.
pub fn simrank(g: &Graph, opts: SimRankOptions) -> Result<Dense> {
    let n = g.node_count();
    if !(opts.importance > 0.0 && opts.importance < 1.0) {
        return Err(Error::InvalidArgument(
            "importance must lie strictly between 0 and 1".into(),
        ));
    }
    if opts.max_iterations == 0 {
        return Err(Error::InvalidArgument(
            "max_iterations must be at least 1".into(),
        ));
    }
    if opts.weighted && opts.stop == SimRankStop::Relative {
        return Err(Error::InvalidArgument(
            "the pairwise SimRank iteration is unweighted".into(),
        ));
    }
    let adj = Weighted::new(g, opts.weighted)?;
    // Column v of W: in-neighbours a of v with weight w_av / Σ_a w_av
    // (a zero column sum is replaced by 1, leaving the column zero).
    let cols: Vec<Vec<(u32, f64)>> = adj
        .inc
        .iter()
        .map(|r| {
            let s: f64 = r.iter().map(|p| p.1).sum();
            let s = if s == 0.0 { 1.0 } else { s };
            r.iter().map(|&(a, x)| (a, x / s)).collect()
        })
        .collect();
    let arcs: usize = cols.iter().map(Vec::len).sum();
    let work = (n as u64).saturating_mul(arcs as u64 + n as u64);
    let threads = parallel::threads_for(opts.threads, n, work);
    let mut s = vec![0.0; n * n];
    for i in 0..n {
        s[i * n + i] = 1.0;
    }
    let mut t = vec![0.0; n * n];
    let mut next = vec![0.0; n * n];
    let c = opts.importance;
    let chunk = n.max(1);
    let mut iterations = 0;
    loop {
        iterations += 1;
        // T = S W, row by row.
        parallel::fill_chunks(
            &mut t,
            chunk,
            threads,
            || (),
            |_, i, row| {
                let srow = &s[i * n..(i + 1) * n];
                for (v, x) in row.iter_mut().enumerate() {
                    *x = cols[v].iter().map(|&(b, w)| srow[b as usize] * w).sum();
                }
            },
        );
        // S' = c Wᵀ T with unit diagonal, row by row.
        parallel::fill_chunks(
            &mut next,
            chunk,
            threads,
            || (),
            |_, u, row| {
                row.iter_mut().for_each(|x| *x = 0.0);
                for &(a, w) in &cols[u] {
                    let trow = &t[a as usize * n..(a as usize + 1) * n];
                    for (x, &y) in row.iter_mut().zip(trow) {
                        *x += w * y;
                    }
                }
                for x in row.iter_mut() {
                    *x *= c;
                }
                row[u] = 1.0;
            },
        );
        let tol = opts.tolerance;
        let close = match opts.stop {
            SimRankStop::Allclose => s
                .iter()
                .zip(&next)
                .all(|(&old, &new)| (old - new).abs() <= tol + 1e-5 * new.abs()),
            SimRankStop::Relative => s
                .iter()
                .zip(&next)
                .all(|(&old, &new)| (new - old).abs() <= tol * (1.0 + old.abs())),
        };
        std::mem::swap(&mut s, &mut next);
        if close || iterations == opts.max_iterations {
            break;
        }
    }
    if iterations == opts.max_iterations {
        return Err(Error::InvalidArgument(format!(
            "SimRank did not converge after {} iterations",
            opts.max_iterations
        )));
    }
    Dense::from_row_major(n, n, s)
}

// ---- Panther -----------------------------------------------------------------

/// Options for [`panther_similarity`].
#[derive(Debug, Clone, Copy)]
pub struct PantherOptions {
    /// Number of similar nodes to return (NetworkX's default is 5).
    pub k: usize,
    /// Steps per random path, at least 2 (5).
    pub path_length: usize,
    /// Constant of the sample size bound (0.5).
    pub c: f64,
    /// Failure probability of the bound (0.1).
    pub delta: f64,
    /// Error bound ε; `None` uses 1/√m.
    pub eps: Option<f64>,
    /// Walk with probabilities proportional to edge weights.
    pub weighted: bool,
    /// Random seed.
    pub seed: u64,
    /// Threads (0 = all cores).
    pub threads: usize,
}

impl Default for PantherOptions {
    fn default() -> Self {
        PantherOptions {
            k: 5,
            path_length: 5,
            c: 0.5,
            delta: 0.1,
            eps: None,
            weighted: false,
            seed: 0,
            threads: 0,
        }
    }
}

/// The nodes most similar to a source, most similar first, with their
/// scores.
#[derive(Debug, Clone, PartialEq)]
pub struct Similar {
    /// Node indices of the input network.
    pub nodes: Vec<u32>,
    /// Their similarity scores.
    pub scores: Vec<f64>,
}

/// Random paths over the non-isolated nodes of a network: for every node,
/// the paths that contain it.
struct Paths {
    /// Kept node -> input node.
    keep: Vec<u32>,
    /// For every kept node, the increasing indices of the paths
    /// containing it.
    by_node: Vec<Vec<u32>>,
    /// Distinct kept nodes of every path.
    nodes_of: Vec<Vec<u32>>,
    samples: usize,
}

fn sample_paths(g: &Graph, source: usize, opts: PantherOptions) -> Result<(Paths, usize)> {
    let n = g.node_count();
    if source >= n {
        return Err(Error::InvalidArgument(format!(
            "source {source} is not a node of the network (0..{n})"
        )));
    }
    if opts.path_length < 2 {
        return Err(Error::InvalidArgument(
            "path_length must be at least 2".into(),
        ));
    }
    let adj = Weighted::new(g, opts.weighted)?;
    let isolated = |v: usize| adj.out[v].is_empty() && adj.inc[v].is_empty();
    if isolated(source) {
        return Err(Error::InvalidArgument(format!(
            "Panther similarity is not defined for the isolated source node {source}"
        )));
    }
    let keep: Vec<u32> = (0..n).filter(|&v| !isolated(v)).map(|v| v as u32).collect();
    let mut index = vec![u32::MAX; n];
    for (i, &v) in keep.iter().enumerate() {
        index[v as usize] = i as u32;
    }
    let kn = keep.len();
    // Edges of the kept network: distinct pairs (arcs), loops once.
    let edges: usize = keep
        .iter()
        .map(|&v| {
            let r = &adj.out[v as usize];
            if g.is_directed() {
                r.len()
            } else {
                r.iter().filter(|p| p.0 >= v).count()
            }
        })
        .sum();
    let eps = match opts.eps {
        Some(e) if e > 0.0 && e.is_finite() => e,
        Some(_) => return Err(Error::InvalidArgument("eps must be positive".into())),
        None => (1.0 / edges as f64).sqrt(),
    };
    if opts.c.is_nan()
        || opts.c <= 0.0
        || opts.delta.is_nan()
        || opts.delta <= 0.0
        || opts.delta >= 1.0
    {
        return Err(Error::InvalidArgument(
            "c must be positive and delta between 0 and 1".into(),
        ));
    }
    let t = opts.path_length as f64;
    let pairs = t * (t - 1.0) / 2.0;
    let bound = opts.c / (eps * eps) * (pairs.log2() + 1.0 + (1.0 / opts.delta).ln());
    if bound.is_nan() || bound >= 1e12 {
        return Err(Error::InvalidArgument(format!(
            "the Panther sample size {bound} is too large"
        )));
    }
    let samples = bound as usize;
    if samples == 0 {
        return Err(Error::InvalidArgument(
            "the Panther sample size is zero; decrease eps or increase c".into(),
        ));
    }
    // Cumulative transition weights of every kept node.
    let mut cum: Vec<Vec<(u32, f64)>> = Vec::with_capacity(kn);
    for &v in &keep {
        let mut total = 0.0;
        let mut row = Vec::with_capacity(adj.out[v as usize].len());
        for &(u, x) in &adj.out[v as usize] {
            if x.is_nan() || x < 0.0 {
                return Err(Error::InvalidArgument(
                    "Panther random walks need non-negative weights".into(),
                ));
            }
            total += x;
            row.push((index[u as usize], total));
        }
        if total.is_nan() || total <= 0.0 {
            return Err(Error::InvalidArgument(format!(
                "node {v} has no out-neighbour (or only zero weights), so random paths \
                 cannot continue from it"
            )));
        }
        cum.push(row);
    }
    let work = (samples as u64).saturating_mul(opts.path_length as u64 * 64);
    let threads = parallel::threads_for(opts.threads, samples, work);
    let nodes_of: Vec<Vec<u32>> = parallel::map_items(
        samples,
        threads,
        || (),
        |_, p| {
            let mut rng = Stream::new(opts.seed, p as u64);
            let mut v = rng.below(kn as u64) as usize;
            let mut path = Vec::with_capacity(opts.path_length + 1);
            path.push(v as u32);
            for _ in 0..opts.path_length {
                let row = &cum[v];
                let x = rng.uniform() * row[row.len() - 1].1;
                let i = row.partition_point(|q| q.1 <= x).min(row.len() - 1);
                v = row[i].0 as usize;
                path.push(v as u32);
            }
            path.sort_unstable();
            path.dedup();
            path
        },
    );
    let mut by_node: Vec<Vec<u32>> = vec![Vec::new(); kn];
    for (p, nodes) in nodes_of.iter().enumerate() {
        for &v in nodes {
            by_node[v as usize].push(p as u32);
        }
    }
    let s = index[source] as usize;
    Ok((
        Paths {
            keep,
            by_node,
            nodes_of,
            samples,
        },
        s,
    ))
}

/// Share of paths through kept node `s` that also pass through every kept
/// node.
fn similarity_row(p: &Paths, s: usize, counts: &mut [u32]) -> Vec<f64> {
    counts.iter_mut().for_each(|c| *c = 0);
    for &path in &p.by_node[s] {
        for &v in &p.nodes_of[path as usize] {
            counts[v as usize] += 1;
        }
    }
    let inv = 1.0 / p.samples as f64;
    counts.iter().map(|&c| c as f64 * inv).collect()
}

/// The `take` highest of `scores` (ties: `first`, then index).
fn top(scores: &[f64], take: usize, first: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| {
        scores[b]
            .total_cmp(&scores[a])
            .then((a != first).cmp(&(b != first)))
            .then(a.cmp(&b))
    });
    order.truncate(take);
    order
}

/// The `k` nodes most similar to `source` by Panther. See the module
/// documentation.
pub fn panther_similarity(g: &Graph, source: usize, opts: PantherOptions) -> Result<Similar> {
    let (p, s) = sample_paths(g, source, opts)?;
    let kn = p.keep.len();
    if kn < opts.k {
        return Err(Error::InvalidArgument(format!(
            "the number of requested nodes {} is greater than the number of nodes {kn}",
            opts.k
        )));
    }
    let mut counts = vec![0u32; kn];
    let row = similarity_row(&p, s, &mut counts);
    let best = top(&row, (opts.k + 1).min(kn), s);
    let mut out = Similar {
        nodes: Vec::new(),
        scores: Vec::new(),
    };
    for i in best {
        if i != s {
            out.nodes.push(p.keep[i]);
            out.scores.push(row[i]);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simrank_of_a_star() {
        // Leaves of a star share their only neighbour: s = c. The centre
        // and a leaf: in-neighbours {leaves} and {centre}, whose
        // similarities are 0, so s = 0.
        let g = Graph::from_edges(4, &[0, 0, 0], &[1, 2, 3], None, false).unwrap();
        let s = simrank(&g, SimRankOptions::default()).unwrap();
        assert!((s.get(1, 2) - 0.9).abs() < 1e-12);
        assert_eq!(s.get(0, 1), 0.0);
        assert_eq!(s.get(3, 3), 1.0);
    }

    #[test]
    fn panther_ranks_neighbours_of_a_path() {
        let g = Graph::from_edges(5, &[0, 1, 2, 3], &[1, 2, 3, 4], None, false).unwrap();
        let o = PantherOptions {
            k: 2,
            eps: Some(0.05),
            seed: 7,
            threads: 1,
            ..Default::default()
        };
        let r = panther_similarity(&g, 2, o).unwrap();
        assert_eq!(r.nodes.len(), 2);
        assert!(r.nodes.contains(&1) && r.nodes.contains(&3));
        let o4 = PantherOptions { threads: 4, ..o };
        assert_eq!(panther_similarity(&g, 2, o4).unwrap(), r);
    }
}
