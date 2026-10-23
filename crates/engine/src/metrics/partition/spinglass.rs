//! Spinglass community detection by simulated annealing.
//!
//! Reichardt and Bornholdt (2006, `reichardt2006statistical`) find
//! communities as the ground state of a Potts model with q spin states:
//! each node i carries a spin s_i in 1..q and the energy is
//!
//! ```text
//! H = - sum_{i<j} (A_ij - gamma p_ij) delta(s_i, s_j),
//! ```
//!
//! with p_ij = k_i k_j / (2m) for the configuration null model (where -H
//! is m times modularity plus a constant when gamma = 1) or p_ij = p = 2m / (n (n - 1))
//! for the simple (Erdos-Renyi) null model, m the total edge weight and
//! k_i the weighted degree. Traag and Bruggeman (2009,
//! `traag2009community`) extend it to signed networks, with A+ and A- the
//! positive weights and the magnitudes of the negative ones, each with
//! its own null model, degrees, and resolution:
//!
//! ```text
//! H = - sum_{i<j} (A+_ij - gamma p+_ij) delta(s_i, s_j)
//!     + sum_{i<j} (A-_ij - gamma_minus p-_ij) delta(s_i, s_j).
//! ```
//!
//! The ground state is approached by simulated annealing (Kirkpatrick,
//! Gelatt, and Vecchi 1983, `kirkpatrick1983optimization`) with
//! heat-bath (Gibbs) updates: a node is moved to spin s with probability
//! proportional to exp(-dH(s) / T), where dH(s) is the change in H. The
//! temperature starts at `start_temp` and is multiplied by `cool_fact`
//! after each sweep until it falls to `stop_temp` or below; a sweep visits
//! every node once in a new random order (`parallel = false`) or updates
//! every node at once from the previous configuration (`parallel =
//! true`). A final zero-temperature pass then moves nodes to their best
//! spin until no move lowers H, so the result is a local minimum of H.
//!
//! These are the options igraph's `cluster_spinglass()` documents
//! (`spins`, `parupdate`, `start.temp`, `stop.temp`, `cool.fact`,
//! `update.rule`, `gamma`, `implementation`, `gamma.minus`); igraph does
//! not document its schedule in more detail, so partitions are not
//! expected to match igraph's, only in quality.
//!
//! Conventions:
//! - Direction is ignored, as in igraph: the network is taken as its
//!   underlying simple undirected network, with the weights of parallel
//!   edges (and of reciprocal arcs) added up. Self-loops never change H
//!   and are ignored.
//! - Connected networks only, as igraph requires; disconnected ones are
//!   rejected (detect communities in each component separately).
//! - Weights (if used) are added as given. The original implementation
//!   requires non-negative weights; the signed implementation accepts any
//!   finite weights. Temperatures are in units of edge weight.
//! - `spins` is the largest number of groups; it is capped at n. The
//!   parallel update is not available for the signed implementation (as
//!   in igraph).
//! - Starting spins are drawn uniformly. Random numbers come from
//!   [`crate::rng::Stream`] with the seed (the start from stream 0, sweep
//!   s from stream s + 1), so a seed reproduces the result on any
//!   platform. The method runs on one thread.
//! - The modularity returned is [`crate::metrics::community::modularity`]
//!   (resolution 1) of the partition on the underlying undirected network,
//!   NaN when a weight is negative or there are no edges. Groups are
//!   labelled 0, 1, ... in order of their first node.

use super::relabel;
use crate::error::{Error, Result};
use crate::graph::{Csr, Graph};
use crate::rng::Stream;

/// Null model of [`spinglass`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinglassNull {
    /// p_ij = p = 2m / (n (n - 1)): same number of edges.
    Simple,
    /// p_ij = k_i k_j / (2m): same degrees.
    Config,
}

impl SpinglassNull {
    /// From the C ABI code: 0 = simple, 1 = config.
    pub fn from_code(code: i32) -> Result<SpinglassNull> {
        match code {
            0 => Ok(SpinglassNull::Simple),
            1 => Ok(SpinglassNull::Config),
            _ => Err(Error::InvalidArgument(format!(
                "unknown spinglass null model code {code} (0 = simple, 1 = config)"
            ))),
        }
    }
}

/// Options for [`spinglass`].
#[derive(Debug, Clone, Copy)]
pub struct SpinglassOptions {
    /// Number of spin states (largest number of groups).
    pub spins: usize,
    /// Update every node at once in each sweep.
    pub parallel: bool,
    /// Starting temperature (> 0).
    pub start_temp: f64,
    /// Stopping temperature (> 0).
    pub stop_temp: f64,
    /// Cooling factor, in (0, 1).
    pub cool_fact: f64,
    /// Null model.
    pub null: SpinglassNull,
    /// Resolution of the positive part.
    pub gamma: f64,
    /// Use the signed (Traag-Bruggeman) Hamiltonian, which allows negative
    /// weights.
    pub signed: bool,
    /// Resolution of the negative part (signed implementation only).
    pub gamma_minus: f64,
    /// Use edge weights.
    pub weighted: bool,
    /// Seed.
    pub seed: u64,
}

impl Default for SpinglassOptions {
    fn default() -> Self {
        SpinglassOptions {
            spins: 25,
            parallel: false,
            start_temp: 1.0,
            stop_temp: 0.01,
            cool_fact: 0.99,
            null: SpinglassNull::Config,
            gamma: 1.0,
            signed: false,
            gamma_minus: 1.0,
            weighted: false,
            seed: 0,
        }
    }
}

/// Result of [`spinglass`].
#[derive(Debug, Clone, PartialEq)]
pub struct Spinglass {
    /// Group of each node.
    pub membership: Vec<usize>,
    /// Modularity of the partition (see the module documentation).
    pub modularity: f64,
    /// Energy H of the final configuration.
    pub energy: f64,
}

struct State<'a> {
    adj: &'a Csr,
    kp: Vec<f64>,
    kn: Vec<f64>,
    /// Coefficients of the null models: config uses c * k_i * K_s, simple
    /// uses c * N_s.
    cp: f64,
    cn: f64,
    config: bool,
    spin: Vec<usize>,
    /// Totals per spin: positive degree, negative degree, node count.
    tp: Vec<f64>,
    tn: Vec<f64>,
    count: Vec<f64>,
    /// Scratch: weight from the node to each spin.
    lp: Vec<f64>,
    ln: Vec<f64>,
    touched: Vec<usize>,
}

impl State<'_> {
    fn add(&mut self, v: usize, s: usize, sign: f64) {
        self.tp[s] += sign * self.kp[v];
        self.tn[s] += sign * self.kn[v];
        self.count[s] += sign;
    }

    /// -dH(s) for every spin s with v removed from its spin, into `gain`.
    fn gains(&mut self, v: usize, gain: &mut [f64]) {
        let w = self.adj.weights(v);
        for (i, &u) in self.adj.neighbors(v).iter().enumerate() {
            let x = w.map_or(1.0, |w| w[i]);
            let s = self.spin[u as usize];
            if self.lp[s] == 0.0 && self.ln[s] == 0.0 {
                self.touched.push(s);
            }
            if x >= 0.0 {
                self.lp[s] += x;
            } else {
                self.ln[s] -= x;
            }
        }
        for (s, g) in gain.iter_mut().enumerate() {
            let (np, nn) = if self.config {
                (
                    self.cp * self.kp[v] * self.tp[s],
                    self.cn * self.kn[v] * self.tn[s],
                )
            } else {
                (self.cp * self.count[s], self.cn * self.count[s])
            };
            *g = (self.lp[s] - np) - (self.ln[s] - nn);
        }
        for &s in &self.touched {
            self.lp[s] = 0.0;
            self.ln[s] = 0.0;
        }
        self.touched.clear();
    }
}

/// Heat-bath choice: s with probability proportional to exp(gain[s] / t).
fn draw(gain: &[f64], t: f64, rng: &mut Stream) -> usize {
    let top = gain.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = gain.iter().map(|&g| ((g - top) / t).exp()).sum();
    let mut x = rng.uniform() * total;
    for (s, &g) in gain.iter().enumerate() {
        x -= ((g - top) / t).exp();
        if x < 0.0 {
            return s;
        }
    }
    // Rounding: fall back to the best spin.
    gain.iter()
        .enumerate()
        .fold(0, |b, (s, &g)| if g > gain[b] { s } else { b })
}

const MAX_QUENCH: usize = 10_000;

/// Spinglass communities of a connected network. See the module
/// documentation.
pub fn spinglass(g: &Graph, opts: &SpinglassOptions) -> Result<Spinglass> {
    let n = g.node_count();
    let bad = |what: &str| Err(Error::InvalidArgument(what.to_string()));
    if opts.spins == 0 {
        return bad("the number of spins must be positive");
    }
    if !(opts.start_temp.is_finite() && opts.start_temp > 0.0) {
        return bad("the starting temperature must be positive and finite");
    }
    if !(opts.stop_temp.is_finite() && opts.stop_temp > 0.0) {
        return bad("the stopping temperature must be positive and finite");
    }
    if !(opts.cool_fact > 0.0 && opts.cool_fact < 1.0) {
        return bad("the cooling factor must be strictly between 0 and 1");
    }
    if !(opts.gamma.is_finite() && opts.gamma >= 0.0) {
        return bad("gamma must be finite and non-negative");
    }
    if !(opts.gamma_minus.is_finite() && opts.gamma_minus >= 0.0) {
        return bad("gamma_minus must be finite and non-negative");
    }
    if opts.signed && opts.parallel {
        return bad("the parallel update is not available for the signed implementation");
    }
    if opts.weighted {
        let w = g.edge_weights().ok_or(Error::MissingWeights)?;
        if !opts.signed {
            if let Some(e) = w.iter().position(|&x| x < 0.0) {
                return Err(Error::InvalidArgument(format!(
                    "edge {e} has negative weight {}; the original spinglass implementation \
                     requires non-negative weights (use the signed implementation)",
                    w[e]
                )));
            }
        }
    }
    if n == 0 {
        return bad("spinglass communities need at least one node");
    }
    let adj = g.simple_undirected(opts.weighted);
    // Connectedness.
    let mut seen = vec![false; n];
    let mut stack = vec![0usize];
    seen[0] = true;
    while let Some(v) = stack.pop() {
        for &u in adj.neighbors(v) {
            if !seen[u as usize] {
                seen[u as usize] = true;
                stack.push(u as usize);
            }
        }
    }
    if seen.iter().any(|&s| !s) {
        return bad(
            "spinglass communities need a connected network; detect communities in each \
             component separately",
        );
    }
    let q = opts.spins.min(n);
    let mut kp = vec![0.0f64; n];
    let mut kn = vec![0.0f64; n];
    for v in 0..n {
        let w = adj.weights(v);
        for i in 0..adj.neighbors(v).len() {
            let x = w.map_or(1.0, |w| w[i]);
            if x >= 0.0 {
                kp[v] += x;
            } else {
                kn[v] -= x;
            }
        }
    }
    let mp: f64 = kp.iter().sum::<f64>() / 2.0;
    let mn: f64 = kn.iter().sum::<f64>() / 2.0;
    let config = opts.null == SpinglassNull::Config;
    let pairs = n as f64 * (n as f64 - 1.0);
    let coef = |gamma: f64, total: f64| -> f64 {
        if total == 0.0 {
            0.0
        } else if config {
            gamma / (2.0 * total)
        } else {
            gamma * 2.0 * total / pairs
        }
    };
    let mut st = State {
        adj: &adj,
        kp,
        kn,
        cp: coef(opts.gamma, mp),
        cn: coef(opts.gamma_minus, mn),
        config,
        spin: vec![0; n],
        tp: vec![0.0; q],
        tn: vec![0.0; q],
        count: vec![0.0; q],
        lp: vec![0.0; q],
        ln: vec![0.0; q],
        touched: Vec::new(),
    };
    let mut rng = Stream::new(opts.seed, 0);
    for v in 0..n {
        let s = rng.below(q as u64) as usize;
        st.spin[v] = s;
        st.add(v, s, 1.0);
    }
    let mut gain = vec![0.0f64; q];
    let mut order: Vec<usize> = (0..n).collect();
    let mut next = vec![0usize; n];
    let mut t = opts.start_temp;
    let mut sweep = 0u64;
    loop {
        sweep += 1;
        let mut rng = Stream::new(opts.seed, sweep);
        if opts.parallel {
            for (v, nx) in next.iter_mut().enumerate() {
                let s = st.spin[v];
                st.add(v, s, -1.0);
                st.gains(v, &mut gain);
                st.add(v, s, 1.0);
                *nx = draw(&gain, t, &mut rng);
            }
            for (v, &s) in next.iter().enumerate() {
                st.add(v, st.spin[v], -1.0);
                st.spin[v] = s;
                st.add(v, s, 1.0);
            }
        } else {
            rng.permutation(&mut order);
            for &v in &order {
                let s = st.spin[v];
                st.add(v, s, -1.0);
                st.gains(v, &mut gain);
                let s = draw(&gain, t, &mut rng);
                st.spin[v] = s;
                st.add(v, s, 1.0);
            }
        }
        if t <= opts.stop_temp {
            break;
        }
        t *= opts.cool_fact;
    }
    // Zero temperature: best moves in node order until none improves.
    for _ in 0..MAX_QUENCH {
        let mut changed = false;
        for v in 0..n {
            let s = st.spin[v];
            st.add(v, s, -1.0);
            st.gains(v, &mut gain);
            let mut best = s;
            for (c, &x) in gain.iter().enumerate() {
                if x > gain[best] {
                    best = c;
                }
            }
            if best != s {
                changed = true;
            }
            st.spin[v] = best;
            st.add(v, best, 1.0);
        }
        if !changed {
            break;
        }
    }
    // Energy of the final configuration, summed over pairs i < j.
    let mut energy = 0.0;
    for v in 0..n {
        let w = adj.weights(v);
        for (i, &u) in adj.neighbors(v).iter().enumerate() {
            if (u as usize) > v && st.spin[u as usize] == st.spin[v] {
                energy -= w.map_or(1.0, |w| w[i]);
            }
        }
    }
    for s in 0..q {
        let (np, nn) = if config {
            let sq_p: f64 = (0..n)
                .filter(|&v| st.spin[v] == s)
                .map(|v| st.kp[v] * st.kp[v])
                .sum();
            let sq_n: f64 = (0..n)
                .filter(|&v| st.spin[v] == s)
                .map(|v| st.kn[v] * st.kn[v])
                .sum();
            (
                st.cp * (st.tp[s] * st.tp[s] - sq_p) / 2.0,
                st.cn * (st.tn[s] * st.tn[s] - sq_n) / 2.0,
            )
        } else {
            let pairs = st.count[s] * (st.count[s] - 1.0) / 2.0;
            (st.cp * pairs, st.cn * pairs)
        };
        energy += np - nn;
    }
    let membership = relabel(&st.spin);
    let modularity =
        if g.edge_weights().is_some_and(|w| w.iter().any(|&x| x < 0.0)) && opts.weighted {
            f64::NAN
        } else {
            let memb: Vec<i64> = membership.iter().map(|&c| c as i64).collect();
            let und;
            let h = if g.is_directed() {
                und = Graph::from_edges(
                    n,
                    &g.edge_sources()
                        .iter()
                        .map(|&x| x as i64)
                        .collect::<Vec<_>>(),
                    &g.edge_targets()
                        .iter()
                        .map(|&x| x as i64)
                        .collect::<Vec<_>>(),
                    g.edge_weights(),
                    false,
                )?;
                &und
            } else {
                g
            };
            crate::metrics::community::modularity(h, &memb, 1.0, opts.weighted)?
        };
    Ok(Spinglass {
        membership,
        modularity,
        energy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::community::modularity;

    fn ring_of_cliques(k: usize, size: usize) -> Graph {
        let (mut s, mut t) = (Vec::new(), Vec::new());
        for c in 0..k {
            let base = c * size;
            for i in 0..size {
                for j in (i + 1)..size {
                    s.push((base + i) as i64);
                    t.push((base + j) as i64);
                }
            }
            s.push((base + size - 1) as i64);
            t.push((((c + 1) % k) * size) as i64);
        }
        Graph::from_edges(k * size, &s, &t, None, false).unwrap()
    }

    #[test]
    fn finds_cliques_and_energy_matches_modularity() {
        let g = ring_of_cliques(5, 5);
        for parallel in [false, true] {
            let opts = SpinglassOptions {
                parallel,
                seed: 3,
                ..Default::default()
            };
            let r = spinglass(&g, &opts).unwrap();
            let truth: Vec<usize> = (0..25).map(|v| v / 5).collect();
            assert_eq!(r.membership, truth, "parallel = {parallel}");
            // -H = m Q + sum_i k_i^2 / (4m) for the configuration model
            // with gamma = 1 (pairs i < j against all ordered pairs).
            let m = g.edge_count() as f64;
            let sq: f64 = (0..25)
                .map(|v| (g.out().neighbors(v).len() as f64).powi(2))
                .sum();
            assert!((-r.energy - m * r.modularity - sq / (4.0 * m)).abs() < 1e-9);
            let memb: Vec<i64> = r.membership.iter().map(|&c| c as i64).collect();
            assert!((modularity(&g, &memb, 1.0, false).unwrap() - r.modularity).abs() < 1e-12);
        }
    }

    #[test]
    fn signed_separates_negative_ties() {
        // Two positive triangles joined by negative edges.
        let s = [0, 1, 0, 3, 4, 3, 0, 1, 2];
        let t = [1, 2, 2, 4, 5, 5, 3, 4, 5];
        let w = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0];
        let g = Graph::from_edges(6, &s, &t, Some(&w), false).unwrap();
        let opts = SpinglassOptions {
            signed: true,
            weighted: true,
            seed: 1,
            ..Default::default()
        };
        let r = spinglass(&g, &opts).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 1, 1, 1]);
        assert!(r.modularity.is_nan());
        let orig = SpinglassOptions {
            weighted: true,
            ..Default::default()
        };
        assert!(spinglass(&g, &orig).is_err());
    }

    #[test]
    fn seeded_and_validated() {
        let g = ring_of_cliques(4, 4);
        let opts = SpinglassOptions {
            seed: 9,
            null: SpinglassNull::Simple,
            ..Default::default()
        };
        assert_eq!(spinglass(&g, &opts).unwrap(), spinglass(&g, &opts).unwrap());
        let two = Graph::from_edges(4, &[0, 2], &[1, 3], None, false).unwrap();
        assert!(spinglass(&two, &SpinglassOptions::default()).is_err());
        let bad = SpinglassOptions {
            cool_fact: 1.0,
            ..Default::default()
        };
        assert!(spinglass(&g, &bad).is_err());
        let one = Graph::from_edges(1, &[], &[], None, false).unwrap();
        assert_eq!(
            spinglass(&one, &SpinglassOptions::default())
                .unwrap()
                .membership,
            vec![0]
        );
    }
}
