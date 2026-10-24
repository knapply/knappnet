//! Measures of electrical current and random walks: current-flow
//! betweenness of nodes and edges (exact, for sources and targets, and
//! approximated from sampled pairs), current-flow closeness, resistance
//! distance, effective graph resistance, Kemeny's constant, and
//! second-order centrality.
//!
//! The network is read as an electrical circuit: every edge is a resistor
//! whose conductance is its weight (or 1), and a unit current enters at a
//! node s and leaves at a node t. With L the weighted Laplacian, the
//! potentials p solve L p = e_s - e_t. Writing C for the inverse of L with
//! one node of each component grounded (its row and column removed, then
//! padded with zeros), p = C (e_s - e_t), the current through edge
//! e = (u, v) is c_e (p_u - p_v), and the effective resistance between i
//! and j is R_ij = C_ii + C_jj - 2 C_ij (Klein and Randic 1993,
//! `klein1993resistance`). C comes from a Cholesky factorization of the
//! grounded Laplacian ([`crate::linalg::cholesky_inverse`]).
//!
//! Sources:
//! - **Current-flow betweenness** (Newman 2005, `newman2005measure`, as
//!   random-walk betweenness; Brandes and Fleischer 2005,
//!   `brandes2005centrality`): the throughput of v for the pair (s, t) is
//!   half the sum of the absolute currents on the edges at v, and v's
//!   betweenness sums it over unordered pairs of other nodes. Following
//!   Brandes and Fleischer, the current of edge e = (u, v) for every pair is
//!   F_e(s) - F_e(t) with F_e = c_e (C[u, .] - C[v, .]), so sorting F_e
//!   gives the sum over all pairs in O(n log n) per edge, O(mn log n + n^3)
//!   in all. The edge measure is half the sum over unordered pairs of the
//!   absolute current through the edge.
//! - **Approximation** (Brandes and Fleischer 2005, section 5): the
//!   average throughput over k random pairs, scaled by
//!   n(n - 1)/((n - 1)(n - 2)), with k = ceil((c / epsilon)^2 log n) and
//!   c = n(n - 1)/((n - 1)(n - 2)), NetworkX's choice; each pair needs one
//!   solve, by conjugate gradients (`hestenes1952methods`) on large
//!   networks.
//! - **Current-flow closeness** (Brandes and Fleischer 2005;
//!   `stephenson1989rethinking`): 1 / sum over w of R_vw.
//! - **Effective graph resistance** (Ellens et al. 2011,
//!   `ellens2011effective`): the sum of R_ij over unordered pairs, also
//!   called the Kirchhoff index.
//! - **Kemeny's constant** (Levene and Loizou 2002,
//!   `levene2002kemeny`): the expected number of steps of a random walk
//!   from any node to a node drawn from the stationary distribution,
//!   sum over i >= 2 of 1 / (1 - lambda_i) with lambda_i the eigenvalues of
//!   D^-1/2 A D^-1/2 other than the largest.
//! - **Second-order centrality** (Kermarrec, Le Merrer, Sericola, and
//!   Tredan 2011, `kermarrec2011second`): the standard deviation of the
//!   return times of a random walk to a node, on a network made regular by
//!   self-loops: with d_max the largest weighted degree, a node of degree
//!   d < d_max gets a self-loop of weight d_max - d, and
//!   sigma_v = sqrt(2 sum over u of M_uv - n(n + 1)), where M_uv is the
//!   mean first passage time from u to v (M_vv the mean return time). The
//!   passage times come from the fundamental matrix
//!   Z = (I - P + 1 pi')^-1: M_uv = (Z_vv - Z_uv) / pi_v and
//!   M_vv = 1 / pi_v (Kemeny and Snell's formulas), one inversion instead of
//!   NetworkX's n linear solves.
//!
//! Conventions:
//! - Direction is ignored and self-loops carry no current. Parallel edges
//!   (and arcs in opposite directions) are parallel resistors: their
//!   conductances add, and a pair's current is shared among them in
//!   proportion to conductance. NetworkX reads only simple undirected
//!   graphs for these measures and agrees there.
//! - Weights: `weights_as` strength (the default for current flow and
//!   Kemeny's constant) uses the weight as a conductance, as NetworkX's
//!   current-flow functions and `kemeny_constant` do; distance uses 1/w (a
//!   resistance), which is what NetworkX's `resistance_distance` and
//!   `effective_graph_resistance` do by default (`invert_weight = True`).
//!   Weights must be positive.
//! - Normalization: node and edge current-flow betweenness divide by
//!   (n - 1)(n - 2)/2, the number of unordered pairs of other nodes, as in
//!   NetworkX. With fewer than three nodes they are 0.
//! - Disconnected networks: current-flow betweenness and closeness, the
//!   approximation, Kemeny's constant, and second-order centrality are
//!   defined for connected networks only and are NaN otherwise (NetworkX
//!   raises an error). Resistance distances between components are
//!   infinite, and so is the effective graph resistance (as NetworkX's).
//!   A single node has current-flow closeness 1 and Kemeny's constant 0, as
//!   in NetworkX.
//! - Subsets ([`current_flow_betweenness_subset`]): the definition sums the
//!   throughput of v over ordered pairs s in S, t in T, s != t, with v not
//!   an endpoint, and normalized divides by (n - 1)(n - 2) (so that S and T
//!   equal to every node give the normalized measure). NetworkX's
//!   `current_flow_betweenness_centrality_subset` instead sums half the
//!   absolute current over all (s, t) in S x T, endpoints included, and adds
//!   the constant 1/(2 - n) to every node; with `networkx` set that formula
//!   is reproduced exactly (it agrees with the definition, normalized, only
//!   when S and T are every node).
//! - Second-order centrality: an existing self-loop counts in the degree
//!   and is replaced by the added loop, as in NetworkX's
//!   `second_order_centrality`; weights are used as given (unweighted
//!   networks use 1).
//! - Sums over edges and pairs are taken in fixed groups merged in order
//!   ([`crate::parallel::sum_in_order`]), so results do not depend on the
//!   number of threads.

use crate::error::{Error, Result};
use crate::graph::{Csr, Graph};
use crate::linalg;
use crate::matrix::Dense;
use crate::parallel;
use crate::rng::Stream;
use crate::union_find::UnionFind;
use crate::weights::WeightKind;

/// Options for the current-flow measures.
#[derive(Debug, Clone, Copy)]
pub struct FlowOptions {
    /// Use edge weights.
    pub weighted: bool,
    /// Weights as conductances (strength, the default) or resistances
    /// (distance).
    pub weights_as: WeightKind,
    /// Normalize (see the module documentation).
    pub normalized: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for FlowOptions {
    fn default() -> Self {
        FlowOptions {
            weighted: false,
            weights_as: WeightKind::Strength,
            normalized: true,
            threads: 0,
        }
    }
}

/// The circuit: conductances between distinct nodes, the conductance of
/// every input edge, and the components.
struct Circuit {
    n: usize,
    rows: Csr,
    degree: Vec<f64>,
    edge_cond: Vec<f64>,
    component: Vec<usize>,
    components: usize,
    threads: usize,
}

impl Circuit {
    #[allow(clippy::needless_range_loop)]
    fn new(g: &Graph, weighted: bool, kind: WeightKind, threads: usize) -> Result<Circuit> {
        let n = g.node_count();
        let m = g.edge_count();
        let (src, tgt) = (g.edge_sources(), g.edge_targets());
        let edge_cond: Vec<f64> = if weighted {
            let w = g.edge_weights().ok_or(Error::MissingWeights)?;
            if let Some(e) = w.iter().position(|&x| x <= 0.0) {
                return Err(Error::InvalidArgument(format!(
                    "edge {e} has weight {}; conductances and resistances must be strictly positive",
                    w[e]
                )));
            }
            match kind {
                WeightKind::Strength => w.to_vec(),
                WeightKind::Distance => w.iter().map(|x| 1.0 / x).collect(),
            }
        } else {
            vec![1.0; m]
        };
        let s: Vec<i64> = src.iter().map(|&v| v as i64).collect();
        let t: Vec<i64> = tgt.iter().map(|&v| v as i64).collect();
        let und = Graph::from_edges(n, &s, &t, Some(&edge_cond), false)?;
        let rows = und.simple_undirected(true);
        let degree = (0..n)
            .map(|v| rows.weights(v).map_or(0.0, |w| w.iter().sum()))
            .collect();
        let mut uf = UnionFind::new(n);
        for v in 0..n {
            for &u in rows.neighbors(v) {
                uf.union(v, u as usize);
            }
        }
        let mut label = vec![usize::MAX; n];
        let mut component = vec![0; n];
        let mut components = 0;
        for v in 0..n {
            let r = uf.find(v);
            if label[r] == usize::MAX {
                label[r] = components;
                components += 1;
            }
            component[v] = label[r];
        }
        let edge_cond = edge_cond
            .iter()
            .enumerate()
            .map(|(e, &c)| if src[e] == tgt[e] { 0.0 } else { c })
            .collect();
        Ok(Circuit {
            n,
            rows,
            degree,
            edge_cond,
            component,
            components,
            threads,
        })
    }

    fn connected(&self) -> bool {
        self.components <= 1
    }

    /// Conductance between u and v (summed over parallel edges).
    fn conductance(&self, u: usize, v: usize) -> f64 {
        let row = self.rows.neighbors(u);
        match row.binary_search(&(v as u32)) {
            Ok(i) => self.rows.weights(u).map_or(1.0, |w| w[i]),
            Err(_) => 0.0,
        }
    }

    /// The grounded inverse of the Laplacian restricted to `nodes` (a
    /// connected set in increasing order), grounding `nodes[0]`: a matrix
    /// over positions 1.. of `nodes`.
    fn grounded_inverse(&self, nodes: &[usize]) -> Result<Dense> {
        let k = nodes.len();
        if k <= 1 {
            return Ok(Dense::zeros(0, 0));
        }
        let mut pos = vec![usize::MAX; self.n];
        for (i, &v) in nodes.iter().enumerate() {
            pos[v] = i;
        }
        let mut a = Dense::zeros(k - 1, k - 1);
        for (i, &v) in nodes.iter().enumerate().skip(1) {
            a.set(i - 1, i - 1, self.degree[v]);
            let w = self.rows.weights(v);
            for (j, &u) in self.rows.neighbors(v).iter().enumerate() {
                let pu = pos[u as usize];
                if pu != usize::MAX && pu > 0 {
                    a.set(i - 1, pu - 1, -w.map_or(1.0, |w| w[j]));
                }
            }
        }
        linalg::spd_inverse_parallel(&a, self.threads).ok_or_else(|| {
            Error::InvalidArgument("the grounded Laplacian is not positive definite".into())
        })
    }
}

/// C(i, j) of a grounded inverse over nodes 0..n grounding node 0.
#[inline]
fn c_at(c: &Dense, i: usize, j: usize) -> f64 {
    if i == 0 || j == 0 {
        0.0
    } else {
        c.get(i - 1, j - 1)
    }
}

/// The sum of |x_s - x_t| over unordered pairs of entries.
fn pair_sum(x: &mut [f64]) -> f64 {
    x.sort_unstable_by(|a, b| a.total_cmp(b));
    let n = x.len() as f64;
    x.iter()
        .enumerate()
        .map(|(i, &v)| v * (2.0 * i as f64 - (n - 1.0)))
        .sum()
}

/// Node and edge values of the current-flow betweenness measures.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowBetweenness {
    /// One value per node.
    pub nodes: Vec<f64>,
    /// One value per edge, in input order.
    pub edges: Vec<f64>,
}

impl FlowBetweenness {
    fn undefined(n: usize, m: usize) -> FlowBetweenness {
        FlowBetweenness {
            nodes: vec![f64::NAN; n],
            edges: vec![f64::NAN; m],
        }
    }
}

/// The distinct pairs (u, v), u < v, of a circuit with their conductance.
fn pairs_of(c: &Circuit) -> Vec<(usize, usize, f64)> {
    let mut out = Vec::new();
    for u in 0..c.n {
        let w = c.rows.weights(u);
        for (j, &v) in c.rows.neighbors(u).iter().enumerate() {
            if (v as usize) > u {
                out.push((u, v as usize, w.map_or(1.0, |w| w[j])));
            }
        }
    }
    out
}

/// Spreads a value per adjacent pair over its parallel edges in proportion
/// to conductance.
fn per_edge(g: &Graph, c: &Circuit, pair_value: impl Fn(usize, usize) -> f64) -> Vec<f64> {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    (0..g.edge_count())
        .map(|e| {
            let (a, b) = (src[e] as usize, tgt[e] as usize);
            if a == b {
                return 0.0;
            }
            let (u, v) = (a.min(b), a.max(b));
            pair_value(u, v) * c.edge_cond[e] / c.conductance(u, v)
        })
        .collect()
}

/// Current-flow betweenness of every node and edge. See the module
/// documentation.
pub fn current_flow_betweenness(g: &Graph, opts: FlowOptions) -> Result<FlowBetweenness> {
    let n = g.node_count();
    let m = g.edge_count();
    let c = Circuit::new(g, opts.weighted, opts.weights_as, opts.threads)?;
    if n < 3 {
        return Ok(FlowBetweenness {
            nodes: vec![0.0; n],
            edges: vec![0.0; m],
        });
    }
    if !c.connected() {
        return Ok(FlowBetweenness::undefined(n, m));
    }
    let all: Vec<usize> = (0..n).collect();
    let inv = c.grounded_inverse(&all)?;
    let pairs = pairs_of(&c);
    let (node_raw, pair_raw) = sum_over_pairs(&pairs, n, opts.threads, |u, v, w, f| {
        for (i, x) in f.iter_mut().enumerate() {
            *x = w * (c_at(&inv, u, i) - c_at(&inv, v, i));
        }
        pair_sum(f)
    });
    let nf = n as f64;
    let scale = if opts.normalized {
        2.0 / ((nf - 1.0) * (nf - 2.0))
    } else {
        1.0
    };
    let nodes = node_raw
        .iter()
        .map(|&s| (0.5 * s - 0.5 * (nf - 1.0)) * scale)
        .collect();
    let edges = per_edge(g, &c, |u, v| {
        let k = pairs.partition_point(|p| (p.0, p.1) < (u, v));
        0.5 * pair_raw[k] * scale
    });
    Ok(FlowBetweenness { nodes, edges })
}

/// Evaluates `value(u, v, conductance, scratch)` for every adjacent pair in
/// parallel and returns, per node, the sum of the values of its pairs and,
/// per pair, the value. The scratch buffer holds n values.
fn sum_over_pairs<F>(
    pairs: &[(usize, usize, f64)],
    n: usize,
    threads: usize,
    value: F,
) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(usize, usize, f64, &mut [f64]) -> f64 + Sync,
{
    let threads = parallel::threads_for(
        threads,
        pairs.len(),
        (pairs.len() as u64).saturating_mul(n as u64 * 8),
    );
    let values = parallel::map_items(
        pairs.len(),
        threads,
        || vec![0.0; n],
        |f, k| {
            let (u, v, w) = pairs[k];
            value(u, v, w, f)
        },
    );
    let mut node = vec![0.0; n];
    for (k, &(u, v, _)) in pairs.iter().enumerate() {
        node[u] += values[k];
        node[v] += values[k];
    }
    (node, values)
}

/// Current-flow betweenness for sources `sources` and targets `targets`.
/// With `networkx`, NetworkX's formula (see the module documentation).
pub fn current_flow_betweenness_subset(
    g: &Graph,
    sources: &[u32],
    targets: &[u32],
    networkx: bool,
    opts: FlowOptions,
) -> Result<FlowBetweenness> {
    let n = g.node_count();
    let m = g.edge_count();
    for &v in sources.iter().chain(targets) {
        if v as usize >= n {
            return Err(Error::InvalidArgument(format!(
                "node index {v} is out of range for a network with {n} nodes"
            )));
        }
    }
    let c = Circuit::new(g, opts.weighted, opts.weights_as, opts.threads)?;
    if !c.connected() || n == 0 {
        return Ok(FlowBetweenness::undefined(n, m));
    }
    let all: Vec<usize> = (0..n).collect();
    let inv = c.grounded_inverse(&all)?;
    let pairs = pairs_of(&c);
    let nf = n as f64;
    // Ordered pairs of distinct nodes (sets) for the definition; the lists
    // as given for NetworkX.
    let (s_list, t_list): (Vec<usize>, Vec<usize>) = if networkx {
        (
            sources.iter().map(|&v| v as usize).collect(),
            targets.iter().map(|&v| v as usize).collect(),
        )
    } else {
        let mut s: Vec<usize> = sources.iter().map(|&v| v as usize).collect();
        let mut t: Vec<usize> = targets.iter().map(|&v| v as usize).collect();
        s.sort_unstable();
        s.dedup();
        t.sort_unstable();
        t.dedup();
        (s, t)
    };
    let (node_raw, pair_raw) = sum_over_pairs(&pairs, n, opts.threads, |u, v, w, f| {
        for (i, x) in f.iter_mut().enumerate() {
            *x = w * (c_at(&inv, u, i) - c_at(&inv, v, i));
        }
        let mut total = 0.0;
        for &s in &s_list {
            for &t in &t_list {
                if networkx || s != t {
                    total += 0.5 * (f[s] - f[t]).abs();
                }
            }
        }
        total
    });
    if networkx {
        let nb = if opts.normalized {
            (nf - 1.0) * (nf - 2.0)
        } else {
            2.0
        };
        let shift = 1.0 / (2.0 - nf);
        let nodes = node_raw.iter().map(|&b| b / nb + shift).collect();
        let edges = per_edge(g, &c, |u, v| {
            let k = pairs.partition_point(|p| (p.0, p.1) < (u, v));
            pair_raw[k] / nb
        });
        return Ok(FlowBetweenness { nodes, edges });
    }
    // Remove the throughput of each node as an endpoint: 1/2 per ordered
    // pair (s, t), s != t, in which it is s or t.
    let mut in_s = vec![false; n];
    let mut in_t = vec![false; n];
    for &s in &s_list {
        in_s[s] = true;
    }
    for &t in &t_list {
        in_t[t] = true;
    }
    let (ns, nt) = (s_list.len() as f64, t_list.len() as f64);
    let scale = if opts.normalized && n > 2 {
        1.0 / ((nf - 1.0) * (nf - 2.0))
    } else {
        1.0
    };
    let nodes = (0..n)
        .map(|v| {
            let as_source = if in_s[v] {
                nt - in_t[v] as u8 as f64
            } else {
                0.0
            };
            let as_target = if in_t[v] {
                ns - in_s[v] as u8 as f64
            } else {
                0.0
            };
            (node_raw[v] - 0.5 * (as_source + as_target)) * scale
        })
        .collect();
    let edges = per_edge(g, &c, |u, v| {
        let k = pairs.partition_point(|p| (p.0, p.1) < (u, v));
        pair_raw[k] * scale
    });
    Ok(FlowBetweenness { nodes, edges })
}

/// How [`current_flow_betweenness_approx`] chooses its pairs.
#[derive(Debug, Clone, Copy)]
pub enum Sample<'a> {
    /// These pairs (s, t), s != t.
    Pairs(&'a [(u32, u32)]),
    /// Random pairs, as many as NetworkX's formula gives for `epsilon` and
    /// `sample_weight` (an error if more than `kmax`), drawn from the
    /// counter-based stream of `seed`.
    Random {
        /// Absolute error tolerance.
        epsilon: f64,
        /// Largest number of pairs allowed.
        kmax: usize,
        /// Multiplier of the number of pairs.
        sample_weight: f64,
        /// Seed of the random stream.
        seed: u64,
    },
}

/// Current-flow betweenness of every node approximated from sampled pairs.
/// See the module documentation.
pub fn current_flow_betweenness_approx(
    g: &Graph,
    sample: Sample,
    opts: FlowOptions,
) -> Result<Vec<f64>> {
    approx_with(g, sample, opts, DENSE_LIMIT)
}

/// Largest network whose potentials come from the dense grounded inverse;
/// larger ones use conjugate gradients.
const DENSE_LIMIT: usize = 3000;

fn approx_with(
    g: &Graph,
    sample: Sample,
    opts: FlowOptions,
    dense_limit: usize,
) -> Result<Vec<f64>> {
    let n = g.node_count();
    let c = Circuit::new(g, opts.weighted, opts.weights_as, opts.threads)?;
    if !c.connected() {
        return Ok(vec![f64::NAN; n]);
    }
    if n < 3 {
        return Ok(vec![0.0; n]);
    }
    let nf = n as f64;
    let nb = (nf - 1.0) * (nf - 2.0);
    let cstar = nf * (nf - 1.0) / nb;
    let drawn: Vec<(u32, u32)>;
    let pairs: &[(u32, u32)] = match sample {
        Sample::Pairs(p) => {
            if let Some(&(s, t)) = p
                .iter()
                .find(|&&(s, t)| s == t || s as usize >= n || t as usize >= n)
            {
                return Err(Error::InvalidArgument(format!(
                    "sampled pair ({s}, {t}) must join two distinct nodes of the network"
                )));
            }
            p
        }
        Sample::Random {
            epsilon,
            kmax,
            sample_weight,
            seed,
        } => {
            if epsilon.is_nan() || epsilon <= 0.0 || sample_weight.is_nan() || sample_weight <= 0.0
            {
                return Err(Error::InvalidArgument(
                    "epsilon and the sample weight must be positive".into(),
                ));
            }
            let k = (sample_weight * ((cstar / epsilon).powi(2) * nf.ln()).ceil()).floor();
            if k > kmax as f64 {
                return Err(Error::InvalidArgument(format!(
                    "the approximation needs {k} pairs, more than kmax = {kmax}; \
                     increase kmax or epsilon"
                )));
            }
            drawn = (0..k as u64)
                .map(|i| {
                    let mut rng = Stream::new(seed, i);
                    let s = rng.below(n as u64) as u32;
                    let mut t = rng.below(n as u64 - 1) as u32;
                    if t >= s {
                        t += 1;
                    }
                    (s, t)
                })
                .collect();
            &drawn
        }
    };
    let k = pairs.len();
    if k == 0 {
        return Ok(vec![0.0; n]);
    }
    // Potentials: from the grounded inverse when its O(n^3) cost is below
    // that of one conjugate-gradient solve per pair (taken as 100 sweeps
    // over the network), otherwise by conjugate gradients.
    let nn = n as f64;
    let per_pair = 100.0 * (nn + c.rows.arc_count() as f64);
    let dense = n <= dense_limit && nn * nn * nn / 3.0 <= k as f64 * per_pair;
    let all: Vec<usize> = (0..n).collect();
    let inv = if dense {
        Some(c.grounded_inverse(&all)?)
    } else {
        None
    };
    let mut total = vec![0.0; n];
    let threads = parallel::threads_for(
        opts.threads,
        k,
        (k as u64).saturating_mul((n + c.rows.arc_count()) as u64 * 20),
    );
    parallel::sum_in_order(
        k,
        threads,
        || vec![0.0; n],
        || vec![0.0; n],
        |p, acc: &mut Vec<f64>, i| {
            let (s, t) = (pairs[i].0 as usize, pairs[i].1 as usize);
            match &inv {
                Some(inv) => {
                    for (v, x) in p.iter_mut().enumerate() {
                        *x = c_at(inv, v, s) - c_at(inv, v, t);
                    }
                }
                None => {
                    let mut b = vec![0.0; n];
                    b[s] = 1.0;
                    b[t] = -1.0;
                    let sol = linalg::conjugate_gradient(
                        |x, out| laplacian_apply(&c, x, out),
                        &c.degree,
                        &b,
                        1e-12,
                        20 * n,
                    );
                    p.copy_from_slice(&sol.x);
                }
            }
            for v in 0..n {
                if v == s || v == t {
                    continue;
                }
                let w = c.rows.weights(v);
                let mut tau = 0.0;
                for (j, &u) in c.rows.neighbors(v).iter().enumerate() {
                    tau += w.map_or(1.0, |w| w[j]) * (p[v] - p[u as usize]).abs();
                }
                acc[v] += 0.5 * tau;
            }
        },
        |acc| {
            for (x, y) in total.iter_mut().zip(acc.iter_mut()) {
                *x += *y;
                *y = 0.0;
            }
        },
    );
    let factor = if opts.normalized { 1.0 } else { nb / 2.0 };
    Ok(total
        .into_iter()
        .map(|x| x * cstar / k as f64 * factor)
        .collect())
}

/// out = L x for the circuit's Laplacian.
fn laplacian_apply(c: &Circuit, x: &[f64], out: &mut [f64]) {
    for v in 0..c.n {
        let w = c.rows.weights(v);
        let mut s = c.degree[v] * x[v];
        for (j, &u) in c.rows.neighbors(v).iter().enumerate() {
            s -= w.map_or(1.0, |w| w[j]) * x[u as usize];
        }
        out[v] = s;
    }
}

/// Resistance distances between every pair of nodes, as an n-by-n matrix
/// (infinite between components). See the module documentation.
pub fn resistance_distance(g: &Graph, opts: FlowOptions) -> Result<Dense> {
    let n = g.node_count();
    let c = Circuit::new(g, opts.weighted, opts.weights_as, opts.threads)?;
    let mut r = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            if c.component[i] != c.component[j] {
                r.set(i, j, f64::INFINITY);
            }
        }
    }
    for comp in 0..c.components {
        let nodes: Vec<usize> = (0..n).filter(|&v| c.component[v] == comp).collect();
        let inv = c.grounded_inverse(&nodes)?;
        for (a, &i) in nodes.iter().enumerate() {
            for (b, &j) in nodes.iter().enumerate() {
                let x = c_at(&inv, a, a) + c_at(&inv, b, b) - 2.0 * c_at(&inv, a, b);
                r.set(i, j, if i == j { 0.0 } else { x });
            }
        }
    }
    Ok(r)
}

/// Current-flow closeness of every node: 1 / the sum of its resistance
/// distances. See the module documentation.
pub fn current_flow_closeness(g: &Graph, opts: FlowOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    if n <= 1 {
        return Ok(vec![1.0; n]);
    }
    let r = resistance_distance(g, opts)?;
    Ok((0..n)
        .map(|v| {
            let s: f64 = r.row(v).iter().sum();
            if s.is_finite() {
                1.0 / s
            } else {
                f64::NAN
            }
        })
        .collect())
}

/// Effective graph resistance: the sum of resistance distances over
/// unordered pairs (infinite for a disconnected network, NaN without
/// nodes).
pub fn effective_graph_resistance(g: &Graph, opts: FlowOptions) -> Result<f64> {
    let n = g.node_count();
    if n == 0 {
        return Ok(f64::NAN);
    }
    let r = resistance_distance(g, opts)?;
    let mut total = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            total += r.get(i, j);
        }
    }
    Ok(total)
}

/// The weighted adjacency matrix as NetworkX builds it for a multigraph:
/// parallel edges add up, direction is ignored, and a self-loop's weight
/// is on the diagonal once.
fn adjacency(g: &Graph, weighted: bool) -> Result<Dense> {
    let n = g.node_count();
    let w = if weighted {
        let w = g.edge_weights().ok_or(Error::MissingWeights)?;
        if let Some(e) = w.iter().position(|&x| x < 0.0) {
            return Err(Error::InvalidArgument(format!(
                "edge {e} has negative weight {}; the weights must be non-negative",
                w[e]
            )));
        }
        Some(w)
    } else {
        None
    };
    let mut a = Dense::zeros(n, n);
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        let (s, t) = (s as usize, t as usize);
        let x = w.map_or(1.0, |w| w[e]);
        a.set(s, t, a.get(s, t) + x);
        if s != t {
            a.set(t, s, a.get(t, s) + x);
        }
    }
    Ok(a)
}

fn connected_dense(a: &Dense) -> bool {
    let n = a.rows();
    let mut uf = UnionFind::new(n);
    for i in 0..n {
        for j in 0..i {
            if a.get(i, j) != 0.0 {
                uf.union(i, j);
            }
        }
    }
    (0..n).all(|v| uf.same(v, 0))
}

/// Kemeny's constant of the random walk on the network (NaN if the network
/// is disconnected or has no nodes). Weights are used as given. See the
/// module documentation.
pub fn kemeny_constant(g: &Graph, weighted: bool) -> Result<f64> {
    let n = g.node_count();
    if n == 0 {
        return Ok(f64::NAN);
    }
    let a = adjacency(g, weighted)?;
    if !connected_dense(&a) {
        return Ok(f64::NAN);
    }
    let d: Vec<f64> = (0..n).map(|i| a.row(i).iter().sum()).collect();
    let mut h = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let x = a.get(i, j);
            if x != 0.0 && d[i] > 0.0 && d[j] > 0.0 {
                h.set(i, j, x / (d[i].sqrt() * d[j].sqrt()));
            }
        }
    }
    let (values, _) = linalg::symmetric_eigen(&h)?;
    Ok(values[..n - 1].iter().map(|&l| 1.0 / (1.0 - l)).sum())
}

/// Second-order centrality of every node (NaN for every node if the
/// network is disconnected). See the module documentation.
#[allow(clippy::needless_range_loop)]
pub fn second_order_centrality(g: &Graph, weighted: bool) -> Result<Vec<f64>> {
    let n = g.node_count();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut a = adjacency(g, weighted)?;
    if !connected_dense(&a) {
        return Ok(vec![f64::NAN; n]);
    }
    let deg: Vec<f64> = (0..n).map(|i| a.row(i).iter().sum()).collect();
    let dmax = deg.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    for i in 0..n {
        if deg[i] < dmax {
            a.set(i, i, dmax - deg[i]);
        }
    }
    let r: Vec<f64> = (0..n).map(|i| a.row(i).iter().sum()).collect();
    let total: f64 = r.iter().sum();
    let pi: Vec<f64> = r.iter().map(|x| x / total).collect();
    // Z = (I - P + 1 pi')^-1.
    let mut z = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let p = a.get(i, j) / r[i];
            let id = if i == j { 1.0 } else { 0.0 };
            z.set(i, j, id - p + pi[j]);
        }
    }
    let z = linalg::invert(&z, 1e-13).ok_or_else(|| {
        Error::InvalidArgument("the fundamental matrix of the random walk is singular".into())
    })?;
    let nf = n as f64;
    Ok((0..n)
        .map(|v| {
            let mut sum = 1.0 / pi[v];
            for u in 0..n {
                if u != v {
                    sum += (z.get(v, v) - z.get(u, v)) / pi[v];
                }
            }
            (2.0 * sum - nf * (nf + 1.0)).sqrt()
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(n: usize, s: &[i64], t: &[i64], w: Option<&[f64]>) -> Graph {
        Graph::from_edges(n, s, t, w, false).unwrap()
    }

    #[test]
    fn path_values() {
        // Path 0-1-2: all current between 0 and 2 passes node 1.
        let p = g(3, &[0, 1], &[1, 2], None);
        let raw = FlowOptions {
            normalized: false,
            ..Default::default()
        };
        let b = current_flow_betweenness(&p, raw).unwrap();
        assert!((b.nodes[1] - 1.0).abs() < 1e-12 && b.nodes[0].abs() < 1e-12);
        // Each edge carries current 1 for two of the three pairs.
        assert!((b.edges[0] - 1.0).abs() < 1e-12);
        let r = resistance_distance(&p, FlowOptions::default()).unwrap();
        assert!((r.get(0, 2) - 2.0).abs() < 1e-12);
        assert!(
            (effective_graph_resistance(&p, FlowOptions::default()).unwrap() - 4.0).abs() < 1e-12
        );
        let cl = current_flow_closeness(&p, FlowOptions::default()).unwrap();
        assert!((cl[1] - 0.5).abs() < 1e-12 && (cl[0] - 1.0 / 3.0).abs() < 1e-12);
        // Kemeny's constant of the path on 3 nodes: eigenvalues of the
        // normalized adjacency are -1, 0, 1, so K = 1/2 + 1.
        assert!((kemeny_constant(&p, false).unwrap() - 1.5).abs() < 1e-12);
    }

    #[test]
    fn subset_with_every_node_matches_the_full_measure() {
        let h = g(
            6,
            &[0, 0, 1, 2, 3, 4, 1, 2],
            &[1, 2, 3, 3, 4, 5, 2, 5],
            Some(&[1.0, 2.0, 0.5, 1.5, 3.0, 1.0, 1.0, 2.0]),
        );
        let opts = FlowOptions {
            weighted: true,
            ..Default::default()
        };
        let full = current_flow_betweenness(&h, opts).unwrap();
        let all: Vec<u32> = (0..6).collect();
        for networkx in [false, true] {
            let sub = current_flow_betweenness_subset(&h, &all, &all, networkx, opts).unwrap();
            for (a, b) in sub.nodes.iter().zip(&full.nodes) {
                assert!(
                    (a - b).abs() < 1e-12,
                    "{networkx} {:?} {:?}",
                    sub.nodes,
                    full.nodes
                );
            }
        }
        // The approximation over every ordered pair is exact.
        let pairs: Vec<(u32, u32)> = (0..6u32)
            .flat_map(|s| (0..6u32).filter(move |&t| t != s).map(move |t| (s, t)))
            .collect();
        let approx = current_flow_betweenness_approx(&h, Sample::Pairs(&pairs), opts).unwrap();
        for (a, b) in approx.iter().zip(&full.nodes) {
            assert!((a - b).abs() < 1e-12);
        }
        // Conjugate gradients give the same potentials as the inverse.
        let cg = approx_with(&h, Sample::Pairs(&pairs), opts, 0).unwrap();
        for (a, b) in cg.iter().zip(&approx) {
            assert!((a - b).abs() < 1e-9, "{cg:?} {approx:?}");
        }
    }

    /// The approximation with potentials from the dense inverse whatever
    /// the number of pairs.
    fn approx_dense(g: &Graph, sample: Sample, opts: FlowOptions) -> Vec<f64> {
        let n = g.node_count();
        let pairs: Vec<(u32, u32)> = match sample {
            Sample::Pairs(p) => p.to_vec(),
            Sample::Random {
                epsilon,
                sample_weight,
                seed,
                ..
            } => {
                let nf = n as f64;
                let cstar = nf * (nf - 1.0) / ((nf - 1.0) * (nf - 2.0));
                let k = (sample_weight * ((cstar / epsilon).powi(2) * nf.ln()).ceil()).floor();
                (0..k as u64)
                    .map(|i| {
                        let mut rng = Stream::new(seed, i);
                        let s = rng.below(n as u64) as u32;
                        let mut t = rng.below(n as u64 - 1) as u32;
                        if t >= s {
                            t += 1;
                        }
                        (s, t)
                    })
                    .collect()
            }
        };
        // Every ordered pair repeated makes the dense inverse the cheaper
        // choice without changing the mean.
        let mut many = Vec::new();
        for _ in 0..(n * n * n / 300 / pairs.len().max(1) + 1) {
            many.extend_from_slice(&pairs);
        }
        approx_with(g, Sample::Pairs(&many), opts, usize::MAX).unwrap()
    }

    #[test]
    fn conjugate_gradients_on_a_larger_random_network() {
        let mut rng = Stream::new(3, 0);
        let n = 300;
        let mut s: Vec<i64> = (0..n as i64 - 1).collect();
        let mut t: Vec<i64> = (1..n as i64).collect();
        for _ in 0..900 {
            s.push(rng.below(n as u64) as i64);
            t.push(rng.below(n as u64) as i64);
        }
        let w: Vec<f64> = (0..s.len()).map(|_| 0.5 + rng.uniform() * 3.0).collect();
        let h = Graph::from_edges(n, &s, &t, Some(&w), false).unwrap();
        let opts = FlowOptions {
            weighted: true,
            threads: 2,
            ..Default::default()
        };
        let sample = Sample::Random {
            epsilon: 2.0,
            kmax: 10_000,
            sample_weight: 1.0,
            seed: 11,
        };
        // Few pairs on a network this size choose conjugate gradients; force
        // the dense inverse for the comparison.
        let dense = approx_dense(&h, sample, opts);
        let cg = approx_with(&h, sample, opts, 0).unwrap();
        for (a, b) in dense.iter().zip(&cg) {
            assert!((a - b).abs() <= 1e-8 * (1.0 + a.abs()), "{a} {b}");
        }
    }
}
