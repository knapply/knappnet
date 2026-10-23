//! Markov chain Monte Carlo sampling of networks from a model.
//!
//! The sampler is the Metropolis-Hastings algorithm (Metropolis et al.
//! 1953, `metropolis1953equation`; Hastings 1970, `hastings1970monte`) on
//! single-dyad toggles, as introduced for these models by Snijders (2002,
//! `snijders2002markov`). A proposed toggle of dyad (i, j) is accepted with
//! probability min(1, exp(theta' delta) * q), where delta is the change in
//! the statistics and q the ratio of proposal probabilities.
//!
//! Two proposals are available:
//!
//! - **Random dyad**: a free dyad is drawn uniformly; q = 1.
//! - **Tie/no-tie (TNT)**: with probability 1/2 a tie is drawn uniformly
//!   from the current ties (proposing to remove it), otherwise a free dyad
//!   is drawn uniformly (Morris, Handcock, and Hunter 2008,
//!   `morris2008specification`). With E ties and D free dyads, removing a
//!   tie has q = E / (E + D) and adding one has q = (D + E + 1) / (E + 1);
//!   when the network is empty the dyad branch is taken with probability 1,
//!   so the first addition has q = (D + 1) / 2 and the removal of a lone tie
//!   q = 2 / (D + 1). TNT mixes far better on sparse networks and is the
//!   default, as in `ergm`.
//!
//! Chains run in parallel through [`crate::parallel`]; chain c draws from
//! stream c of the seed ([`crate::rng::Stream`]), so results are identical
//! for any thread count. Every chain starts from the observed network
//! unless asked to start from the empty network.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::parallel;
use crate::rng::Stream;

use super::model::Model;
use super::state::NetState;
use super::terms::Terms;

/// How toggles are proposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Proposal {
    /// Tie/no-tie.
    #[default]
    Tnt,
    /// A uniformly random free dyad.
    Random,
}

impl Proposal {
    /// Decodes the C interface codes: 0 = TNT, 1 = random dyad.
    pub fn from_code(code: i32) -> Result<Proposal> {
        match code {
            0 => Ok(Proposal::Tnt),
            1 => Ok(Proposal::Random),
            _ => Err(Error::InvalidArgument(format!(
                "unknown proposal code {code}"
            ))),
        }
    }
}

/// Settings for one simulation run.
#[derive(Debug, Clone)]
pub struct SampleOptions {
    /// Proposals discarded before the first sample.
    pub burnin: usize,
    /// Proposals between consecutive samples (at least 1).
    pub interval: usize,
    /// Samples per chain.
    pub nsim: usize,
    /// Independent chains.
    pub chains: usize,
    /// Proposal type.
    pub proposal: Proposal,
    /// Seed; chain c uses stream `stream_offset + c`.
    pub seed: u64,
    /// First stream index, so repeated runs within one estimation use
    /// fresh streams.
    pub stream_offset: u64,
    /// Threads (0 = all cores).
    pub threads: usize,
    /// Whether to keep each sampled network's ties.
    pub keep_networks: bool,
    /// Start from the empty network instead of the observed one.
    pub start_empty: bool,
}

impl Default for SampleOptions {
    fn default() -> Self {
        SampleOptions {
            burnin: 10_000,
            interval: 1_000,
            nsim: 1_000,
            chains: 1,
            proposal: Proposal::Tnt,
            seed: 0,
            stream_offset: 0,
            threads: 0,
            keep_networks: false,
            start_empty: false,
        }
    }
}

/// Networks and statistics sampled from a model.
#[derive(Debug, Clone)]
pub struct Sample {
    /// Number of statistics.
    pub width: usize,
    /// Samples per chain.
    pub nsim: usize,
    /// Number of chains.
    pub chains: usize,
    /// Statistics, chain by chain, sample by sample: `chains * nsim` rows
    /// of `width` values.
    pub stats: Vec<f64>,
    /// Ties of each sampled network, in the same order, if kept.
    pub networks: Option<Vec<Vec<(u32, u32)>>>,
    /// Proposals accepted, over all chains.
    pub accepted: u64,
    /// Proposals made, over all chains.
    pub proposed: u64,
    /// Number of nodes.
    pub n: usize,
    /// Whether ties are directed.
    pub directed: bool,
}

impl Sample {
    /// Row `r` of the statistics.
    #[inline]
    pub fn row(&self, r: usize) -> &[f64] {
        &self.stats[r * self.width..(r + 1) * self.width]
    }

    /// Number of rows.
    #[inline]
    pub fn rows(&self) -> usize {
        self.chains * self.nsim
    }

    /// Mean of each statistic over every row.
    pub fn mean(&self) -> Vec<f64> {
        let mut m = vec![0.0; self.width];
        for r in 0..self.rows() {
            for (k, v) in self.row(r).iter().enumerate() {
                m[k] += v;
            }
        }
        let rows = self.rows().max(1) as f64;
        m.iter_mut().for_each(|v| *v /= rows);
        m
    }
}

/// The free ties of a constrained network, kept separately so that the
/// tie branch of TNT never proposes a fixed tie.
struct FreeTies {
    ties: Vec<(u32, u32)>,
    pos: HashMap<u64, usize>,
}

impl FreeTies {
    fn new(state: &NetState, constraints: &super::state::Constraints) -> FreeTies {
        let mut f = FreeTies {
            ties: Vec::new(),
            pos: HashMap::new(),
        };
        for &(a, b) in state.ties() {
            if !constraints.is_fixed(a as usize, b as usize) {
                f.insert(a, b);
            }
        }
        f
    }

    fn insert(&mut self, a: u32, b: u32) {
        self.pos
            .insert((a as u64) << 32 | b as u64, self.ties.len());
        self.ties.push((a, b));
    }

    fn remove(&mut self, a: u32, b: u32) {
        if let Some(p) = self.pos.remove(&((a as u64) << 32 | b as u64)) {
            let last = self.ties.len() - 1;
            self.ties.swap(p, last);
            self.ties.pop();
            if p != last {
                let m = self.ties[p];
                self.pos.insert((m.0 as u64) << 32 | m.1 as u64, p);
            }
        }
    }
}

/// One Markov chain over a network state.
pub struct Chain<'a> {
    terms: &'a Terms,
    free: Option<&'a [(u32, u32)]>,
    free_ties: Option<FreeTies>,
    dyads: usize,
    /// Current network.
    pub state: NetState,
    /// Current statistics.
    pub stats: Vec<f64>,
    delta: Vec<f64>,
    rng: Stream,
    proposal: Proposal,
    accepted: u64,
    proposed: u64,
}

impl<'a> Chain<'a> {
    /// A chain over `state` with the model's terms, drawing from `rng`.
    pub fn new(model: &'a Model, state: NetState, rng: Stream, proposal: Proposal) -> Chain<'a> {
        let stats = model.terms().stats(&state);
        let free = model.constraints().map(|c| c.free_dyads());
        let free_ties = model.constraints().map(|c| FreeTies::new(&state, c));
        let dyads = model.free_dyad_count();
        Chain {
            terms: model.terms(),
            free,
            free_ties,
            dyads,
            state,
            stats,
            delta: vec![0.0; model.width()],
            rng,
            proposal,
            accepted: 0,
            proposed: 0,
        }
    }

    /// A uniformly random free dyad.
    fn random_dyad(&mut self) -> (usize, usize) {
        if let Some(free) = self.free {
            let d = free[self.rng.below(free.len() as u64) as usize];
            return (d.0 as usize, d.1 as usize);
        }
        let n = self.state.node_count() as u64;
        if self.state.is_directed() {
            let i = self.rng.below(n);
            let mut j = self.rng.below(n - 1);
            if j >= i {
                j += 1;
            }
            (i as usize, j as usize)
        } else {
            // An unordered pair: draw an ordered pair and sort it.
            let i = self.rng.below(n);
            let mut j = self.rng.below(n - 1);
            if j >= i {
                j += 1;
            }
            (i.min(j) as usize, i.max(j) as usize)
        }
    }

    /// Makes `steps` proposals at `theta`.
    pub fn run(&mut self, theta: &[f64], steps: usize) {
        if self.dyads == 0 {
            return;
        }
        for _ in 0..steps {
            self.proposed += 1;
            let d = self.dyads as f64;
            // Ties the proposal may remove: all of them, or the free ones.
            let e = match &self.free_ties {
                Some(f) => f.ties.len(),
                None => self.state.tie_count(),
            };
            let (i, j, log_q) = match self.proposal {
                Proposal::Random => {
                    let (i, j) = self.random_dyad();
                    (i, j, 0.0)
                }
                Proposal::Tnt => {
                    // Removing the last tie leaves an empty network, whose
                    // reverse move always takes the dyad branch.
                    let remove_q = if e == 1 {
                        (2.0 / (d + 1.0)).ln()
                    } else {
                        ((e as f64) / (e as f64 + d)).ln()
                    };
                    if e > 0 && self.rng.uniform() < 0.5 {
                        let k = self.rng.below(e as u64) as usize;
                        let (a, b) = match &self.free_ties {
                            Some(f) => f.ties[k],
                            None => self.state.ties()[k],
                        };
                        (a as usize, b as usize, remove_q)
                    } else {
                        let (i, j) = self.random_dyad();
                        if self.state.has(i, j) {
                            // A tie reached through the dyad branch has the
                            // same proposal probability as one drawn from
                            // the tie list.
                            (i, j, remove_q)
                        } else if e == 0 {
                            (i, j, ((d + 1.0) / 2.0).ln())
                        } else {
                            (i, j, ((d + e as f64 + 1.0) / (e as f64 + 1.0)).ln())
                        }
                    }
                }
            };
            let present = self.state.has(i, j);
            self.terms.change(&self.state, i, j, &mut self.delta);
            let sign = if present { -1.0 } else { 1.0 };
            let mut log_ratio = log_q;
            for (t, dl) in theta.iter().zip(&self.delta) {
                log_ratio += sign * t * dl;
            }
            if log_ratio >= 0.0 || self.rng.uniform().ln() < log_ratio {
                self.state.toggle(i, j);
                if let Some(f) = &mut self.free_ties {
                    let (a, b) = self.state.canonical(i, j);
                    if present {
                        f.remove(a, b);
                    } else {
                        f.insert(a, b);
                    }
                }
                for (s, dl) in self.stats.iter_mut().zip(&self.delta) {
                    *s += sign * dl;
                }
                self.accepted += 1;
            }
        }
    }
}

/// Draws a sample from `model` at coefficients `theta`.
pub fn simulate(model: &Model, theta: &[f64], opts: &SampleOptions) -> Result<Sample> {
    let p = model.width();
    if theta.len() != p {
        return Err(Error::LengthMismatch {
            what: "coefficients",
            expected: p,
            found: theta.len(),
        });
    }
    if let Some(k) = theta.iter().position(|t| !t.is_finite()) {
        return Err(Error::InvalidArgument(format!(
            "coefficient {k} is not finite"
        )));
    }
    if opts.interval == 0 {
        return Err(Error::InvalidArgument("interval must be at least 1".into()));
    }
    if opts.chains == 0 {
        return Err(Error::InvalidArgument(
            "there must be at least one chain".into(),
        ));
    }
    let n = model.state().node_count();
    let directed = model.state().is_directed();
    let start = if opts.start_empty {
        NetState::empty(n, directed)
    } else {
        model.state().clone()
    };
    struct ChainOut {
        chain: usize,
        stats: Vec<f64>,
        nets: Vec<Vec<(u32, u32)>>,
        accepted: u64,
        proposed: u64,
    }
    let outs =
        parallel::for_each_item(opts.chains, opts.threads, Vec::<ChainOut>::new, |acc, c| {
            let rng = Stream::new(opts.seed, opts.stream_offset + c as u64);
            let mut chain = Chain::new(model, start.clone(), rng, opts.proposal);
            let mut stats = Vec::with_capacity(opts.nsim * p);
            let mut nets = Vec::new();
            chain.run(theta, opts.burnin);
            for _ in 0..opts.nsim {
                chain.run(theta, opts.interval);
                stats.extend_from_slice(&chain.stats);
                if opts.keep_networks {
                    nets.push(chain.state.ties().to_vec());
                }
            }
            acc.push(ChainOut {
                chain: c,
                stats,
                nets,
                accepted: chain.accepted,
                proposed: chain.proposed,
            });
        });
    let mut all: Vec<ChainOut> = outs.into_iter().flatten().collect();
    all.sort_by_key(|o| o.chain);
    let mut stats = Vec::with_capacity(opts.chains * opts.nsim * p);
    let mut nets = Vec::new();
    let (mut accepted, mut proposed) = (0, 0);
    for o in all {
        stats.extend(o.stats);
        nets.extend(o.nets);
        accepted += o.accepted;
        proposed += o.proposed;
    }
    Ok(Sample {
        width: p,
        nsim: opts.nsim,
        chains: opts.chains,
        stats,
        networks: opts.keep_networks.then_some(nets),
        accepted,
        proposed,
        n,
        directed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ergm::terms::{Covariates, Term, Terms};

    fn edges_model(n: usize, directed: bool) -> Model {
        let state = NetState::from_edges(n, &[0, 1], &[1, 2], directed).unwrap();
        let terms = Terms::new(vec![Term::Edges], n, directed, Covariates::default()).unwrap();
        Model::new(state, terms, None).unwrap()
    }

    /// Under an edges-only model every dyad is an independent Bernoulli
    /// tie with probability logistic(theta), so the sampled tie count has a
    /// known mean, for both proposals.
    #[test]
    fn edges_model_has_bernoulli_mean() {
        for proposal in [Proposal::Tnt, Proposal::Random] {
            let model = edges_model(8, false);
            let theta = [-1.0];
            let opts = SampleOptions {
                burnin: 2000,
                interval: 20,
                nsim: 4000,
                chains: 2,
                proposal,
                seed: 3,
                ..Default::default()
            };
            let s = simulate(&model, &theta, &opts).unwrap();
            let p = 1.0 / (1.0 + 1f64.exp());
            let expect = 28.0 * p;
            let mean = s.mean()[0];
            // sd of the mean is well below 0.1 for this many draws.
            assert!(
                (mean - expect).abs() < 0.25,
                "{proposal:?}: {mean} vs {expect}"
            );
            assert_eq!(s.rows(), 8000);
        }
    }

    #[test]
    fn chains_do_not_depend_on_threads_and_networks_are_kept() {
        let model = edges_model(6, true);
        let base = SampleOptions {
            burnin: 100,
            interval: 5,
            nsim: 50,
            chains: 3,
            seed: 9,
            keep_networks: true,
            ..Default::default()
        };
        let a = simulate(&model, &[0.2], &base).unwrap();
        let b = simulate(
            &model,
            &[0.2],
            &SampleOptions {
                threads: 1,
                ..base.clone()
            },
        )
        .unwrap();
        assert_eq!(a.stats, b.stats);
        let nets = a.networks.as_ref().unwrap();
        assert_eq!(nets.len(), 150);
        for (r, net) in nets.iter().enumerate() {
            assert_eq!(net.len() as f64, a.row(r)[0]);
        }
    }
}
