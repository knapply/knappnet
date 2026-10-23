//! Ministep simulation of one period (Snijders 2001,
//! `snijders2001statistical`; RSiena manual section 6.1.3 and the
//! `Siena_algorithms` document, section "Outline of model dynamics").
//!
//! Starting from the network observed at the beginning of the period, the
//! continuous-time Markov chain is simulated step by step. With a constant
//! rate function every actor has rate rho, so the waiting time to the next
//! ministep is exponential with rate n rho and the actor is uniform. The
//! actor i evaluates every alter j: toggling the tie i -> j changes the
//! evaluation function by h_j = +/- sum_k beta_k Delta_kij(x) (+ when the
//! tie is created, - when dropped), and h_i = 0 for leaving the network
//! unchanged; option j is chosen with probability
//! exp(h_j) / sum_h exp(h_h).
//!
//! The simulation stops after one time unit (unconditional) or when the
//! number of tie variables that differ from the start of the period
//! reaches the observed number (conditional; Snijders 2001, section 4.2).
//! In the conditional case the rate does not matter, so the chain is run
//! with rho = 1 and the time at which it stops estimates the period's rate
//! parameter.
//!
//! Along the way the score function (derivative of the log-likelihood of
//! the simulated path with respect to the parameters; Schweinberger and
//! Snijders 2007, `schweinberger2007markov`) is accumulated: for each
//! choice, the chosen option's derivative minus its expectation under the
//! choice probabilities, and for a rate rho in the unconditional case
//! K / rho - n, where K is the number of ministeps in the period.

use super::data::Panel;
use super::effects::{change_statistics, ActorContext, Model};
use super::network::Network;
use crate::error::{Error, Result};
use crate::rng::Stream;

/// What happens in one period.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Stopping {
    /// Run for one time unit with the given rate per actor.
    Unconditional {
        /// Rate parameter rho_m of the period.
        rate: f64,
    },
    /// Run with rate 1 until this many tie variables differ from the
    /// start of the period.
    Conditional {
        /// The observed number of changed tie variables.
        distance: usize,
    },
}

/// The outcome of simulating one period.
#[derive(Debug, Clone)]
pub struct PeriodSimulation {
    /// The network at the end of the period.
    pub network: Network,
    /// Effect statistics of the final network (period covariates).
    pub statistics: Vec<f64>,
    /// Score for each effect parameter.
    pub scores: Vec<f64>,
    /// Score for the rate parameter (unconditional only; 0 otherwise).
    pub rate_score: f64,
    /// Number of tie variables that differ from the start of the period.
    pub distance: usize,
    /// Simulated time at the end (the stopping time under rate 1 in the
    /// conditional case, 1 in the unconditional case).
    pub time: f64,
    /// Number of ministeps.
    pub steps: usize,
}

/// Reusable buffers for simulating periods.
pub struct Simulator<'a> {
    panel: &'a Panel,
    model: &'a Model,
    ctx: ActorContext,
    delta: Vec<f64>,
    weights: Vec<f64>,
    /// Ministep limit per period for conditional simulation, after which
    /// the simulation is deemed unable to reach the observed distance.
    max_steps: usize,
}

impl<'a> Simulator<'a> {
    /// Buffers for `panel` and `model`.
    pub fn new(panel: &'a Panel, model: &'a Model) -> Simulator<'a> {
        let n = panel.n();
        Simulator {
            panel,
            model,
            ctx: ActorContext::new(n),
            delta: vec![0.0; n * model.len()],
            weights: vec![0.0; n],
            // RSiena stops a conditional simulation that runs far longer
            // than expected ("unlikely to terminate this epoch"). The
            // observed distance is reached in about n^2 steps even when
            // nearly every tie must flip, so this bound only catches
            // parameter values under which the chain gets stuck.
            max_steps: 1_000_000usize.max(200 * n * n),
        }
    }

    /// Simulates period `period` from `start` with effect parameters
    /// `beta`, drawing random numbers from `rng`.
    pub fn period(
        &mut self,
        period: usize,
        start: &Network,
        beta: &[f64],
        stopping: Stopping,
        rng: &mut Stream,
    ) -> Result<PeriodSimulation> {
        let n = self.panel.n();
        let k = self.model.len();
        debug_assert_eq!(beta.len(), k);
        let mut net = start.clone();
        let mut scores = vec![0.0; k];
        let mut time = 0.0;
        let mut steps = 0usize;
        let mut distance = 0usize;
        let (rate, target) = match stopping {
            Stopping::Unconditional { rate } => {
                if !(rate > 0.0 && rate.is_finite()) {
                    return Err(Error::InvalidArgument(format!(
                        "rate parameter must be positive and finite, not {rate}"
                    )));
                }
                (rate, None)
            }
            Stopping::Conditional { distance } => (1.0, Some(distance)),
        };
        let total_rate = rate * n as f64;
        loop {
            if let Some(t) = target {
                if distance >= t {
                    break;
                }
                if steps >= self.max_steps {
                    return Err(Error::InvalidArgument(format!(
                        "the conditional simulation of period {} did not reach the observed \
                         number of changes ({}) within {} ministeps; the parameters make the \
                         observed changes very unlikely. Try other initial values or \
                         unconditional estimation",
                        period + 1,
                        t,
                        self.max_steps
                    )));
                }
            }
            // Exponential waiting time with rate n rho; the uniform is in
            // [0, 1), so 1 - u is in (0, 1] and the logarithm is finite.
            let wait = -(1.0 - rng.uniform()).ln() / total_rate;
            if target.is_none() && time + wait > 1.0 {
                break;
            }
            time += wait;
            steps += 1;
            let i = rng.below(n as u64) as usize;
            self.ctx.prepare(self.model, &net, i);
            change_statistics(
                self.model,
                self.panel,
                period,
                &net,
                &self.ctx,
                i,
                &mut self.delta,
            );
            // h_j for every option; the no-change option j = i has h = 0.
            let row = net.row(i);
            let mut max_h = 0.0f64;
            for (j, &x) in row.iter().enumerate() {
                let h = if j == i {
                    0.0
                } else {
                    let sign = if x != 0 { -1.0 } else { 1.0 };
                    let d = &self.delta[j * k..(j + 1) * k];
                    sign * d.iter().zip(beta).map(|(a, b)| a * b).sum::<f64>()
                };
                self.weights[j] = h;
                if h > max_h {
                    max_h = h;
                }
            }
            let mut total = 0.0;
            for w in self.weights.iter_mut() {
                *w = (*w - max_h).exp();
                total += *w;
            }
            let u = rng.uniform() * total;
            let mut acc = 0.0;
            let mut chosen = n - 1;
            for (j, &w) in self.weights.iter().enumerate() {
                acc += w;
                if u < acc {
                    chosen = j;
                    break;
                }
            }
            // Score: derivative of the chosen option minus its expectation.
            for (j, &x) in row.iter().enumerate() {
                if j == i {
                    continue;
                }
                let p = self.weights[j] / total;
                let sign = if x != 0 { -1.0 } else { 1.0 };
                let d = &self.delta[j * k..(j + 1) * k];
                for kk in 0..k {
                    scores[kk] -= p * sign * d[kk];
                }
            }
            if chosen != i {
                let sign = if row[chosen] != 0 { -1.0 } else { 1.0 };
                let d = &self.delta[chosen * k..(chosen + 1) * k];
                for kk in 0..k {
                    scores[kk] += sign * d[kk];
                }
                let was_different = net.has(i, chosen) != start.has(i, chosen);
                net.toggle(i, chosen);
                if was_different {
                    distance -= 1;
                } else {
                    distance += 1;
                }
            }
        }
        let rate_score = match stopping {
            Stopping::Unconditional { rate } => steps as f64 / rate - n as f64,
            Stopping::Conditional { .. } => 0.0,
        };
        if target.is_none() {
            time = 1.0;
        }
        let statistics = self.model.statistics(self.panel, period, &net);
        Ok(PeriodSimulation {
            network: net,
            statistics,
            scores,
            rate_score,
            distance,
            time,
            steps,
        })
    }
}

/// A simulated panel: one network per period, each starting from the
/// observed network at the beginning of its period.
#[derive(Debug, Clone)]
pub struct PanelSimulation {
    /// The simulated network at the end of each period.
    pub networks: Vec<Network>,
    /// Effect statistics per period (period-major, `k` per period).
    pub statistics: Vec<f64>,
    /// Simulated time per period (see [`PeriodSimulation::time`]).
    pub times: Vec<f64>,
    /// Changed tie variables per period.
    pub distances: Vec<usize>,
}

/// Simulates every period of the panel once. `rates` has one value per
/// period and is used only when `conditional` is false. Stream `index` of
/// `seed` provides the random numbers.
pub fn simulate_panel(
    panel: &Panel,
    model: &Model,
    beta: &[f64],
    rates: &[f64],
    conditional: bool,
    seed: u64,
    index: u64,
) -> Result<PanelSimulation> {
    if beta.len() != model.len() {
        return Err(Error::LengthMismatch {
            what: "effect parameters",
            expected: model.len(),
            found: beta.len(),
        });
    }
    if !conditional && rates.len() != panel.period_count() {
        return Err(Error::LengthMismatch {
            what: "rate parameters",
            expected: panel.period_count(),
            found: rates.len(),
        });
    }
    let observed = panel.distances();
    let mut sim = Simulator::new(panel, model);
    let mut rng = Stream::new(seed, index);
    let mut out = PanelSimulation {
        networks: Vec::with_capacity(panel.period_count()),
        statistics: Vec::with_capacity(panel.period_count() * model.len()),
        times: Vec::with_capacity(panel.period_count()),
        distances: Vec::with_capacity(panel.period_count()),
    };
    for m in 0..panel.period_count() {
        let stopping = if conditional {
            Stopping::Conditional {
                distance: observed[m] as usize,
            }
        } else {
            Stopping::Unconditional { rate: rates[m] }
        };
        let r = sim.period(m, panel.wave(m), beta, stopping, &mut rng)?;
        out.statistics.extend(r.statistics);
        out.times.push(r.time);
        out.distances.push(r.distance);
        out.networks.push(r.network);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::saom::effects::{Effect, EffectKind};

    fn panel() -> Panel {
        let n = 12;
        let mut rng = Stream::new(3, 9);
        let mut waves = Vec::new();
        for _ in 0..2 {
            let mut a = vec![0u8; n * n];
            for i in 0..n {
                for j in 0..n {
                    if i != j && rng.uniform() < 0.25 {
                        a[i * n + j] = 1;
                    }
                }
            }
            waves.push(a);
        }
        Panel::from_adjacency(n, &waves).unwrap()
    }

    #[test]
    fn conditional_simulation_stops_at_the_observed_distance() {
        let p = panel();
        let m = Model::new(
            &p,
            &[
                Effect::structural(EffectKind::Density),
                Effect::structural(EffectKind::Recip),
            ],
        )
        .unwrap();
        let s = simulate_panel(&p, &m, &[-1.5, 1.0], &[], true, 1, 0).unwrap();
        assert_eq!(s.distances[0] as f64, p.distances()[0]);
        assert_eq!(s.networks[0].distance(p.wave(0)), s.distances[0]);
        assert!(s.times[0] > 0.0);
        let again = simulate_panel(&p, &m, &[-1.5, 1.0], &[], true, 1, 0).unwrap();
        assert_eq!(again.networks[0], s.networks[0]);
        assert_ne!(
            simulate_panel(&p, &m, &[-1.5, 1.0], &[], true, 1, 1)
                .unwrap()
                .networks[0],
            s.networks[0]
        );
    }

    #[test]
    fn unconditional_simulation_runs_one_time_unit() {
        let p = panel();
        let m = Model::new(&p, &[Effect::structural(EffectKind::Density)]).unwrap();
        let s = simulate_panel(&p, &m, &[-1.0], &[3.0], false, 5, 0).unwrap();
        assert_eq!(s.times[0], 1.0);
        assert!(simulate_panel(&p, &m, &[-1.0], &[0.0], false, 5, 0).is_err());
        assert!(simulate_panel(&p, &m, &[-1.0, 2.0], &[3.0], false, 5, 0).is_err());
    }

    #[test]
    fn choice_probabilities_follow_the_evaluation_function() {
        // Two actors: actor 0 toggles its single tie. With density beta,
        // the probability of creating the tie from the empty network is
        // exp(beta) / (1 + exp(beta)) at each ministep. Count first steps.
        let p = Panel::from_edge_lists(2, &[vec![], vec![(0, 1)]]).unwrap();
        let m = Model::new(&p, &[Effect::structural(EffectKind::Density)]).unwrap();
        let beta = 1.0f64;
        let (mut created, mut single) = (0u32, 0u32);
        let mut sim = Simulator::new(&p, &m);
        for t in 0..4000 {
            let mut rng = Stream::new(11, t);
            let r = sim
                .period(
                    0,
                    p.wave(0),
                    &[beta],
                    Stopping::Unconditional { rate: 0.5 },
                    &mut rng,
                )
                .unwrap();
            // Among runs with exactly one ministep (either actor), the
            // share that created a tie estimates exp(beta) / (1 + exp(beta)).
            if r.steps == 1 {
                single += 1;
                created += (r.network.tie_count() == 1) as u32;
            }
        }
        let share = created as f64 / single as f64;
        let expect = beta.exp() / (1.0 + beta.exp());
        assert!(single > 1000, "{single}");
        assert!((share - expect).abs() < 0.04, "{share} vs {expect}");
    }
}
