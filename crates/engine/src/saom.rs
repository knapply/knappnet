//! Stochastic actor-oriented models (SAOMs) for network panel data.
//!
//! A SAOM (Snijders 1996, `snijders1996stochastic`; Snijders 2001,
//! `snijders2001statistical`; Snijders, van de Bunt, and Steglich 2010,
//! `snijders2010introduction`; Snijders 2017, `snijders2017stochastic`)
//! describes how a directed network observed at two or more moments
//! changed between observations as a continuous-time Markov chain of
//! *ministeps*: actor i gets an opportunity to change at rate
//! lambda_i, and then adds or drops one outgoing tie, or leaves the
//! network as it is, choosing option x' with probability proportional to
//! exp(f_i(x')), where the evaluation function
//! f_i(x) = sum_k beta_k s_ik(x) is a weighted sum of *effects*. The
//! engine implements the model as defined in the RSiena manual (Ripley et
//! al. 2026, `ripley2026manual`, chapters 6 and 13, and the companion
//! document `Siena_algorithms`), which is the reference implementation
//! the results are validated against.
//!
//! Three layers are provided:
//!
//! - [`data`]: a [`Panel`] of M >= 2 binary directed networks on the same
//!   n actors, with constant and changing actor covariates and constant
//!   dyadic covariates.
//! - [`effects`]: the effects (RSiena short names in parentheses) with their
//!   statistics s_ik(x) and change statistics
//!   Delta_kij(x) = s_ik(x with i->j) - s_ik(x without i->j).
//! - [`simulate`] and [`mod@estimate`]: ministep simulation of a period, and
//!   estimation by the method of moments through the three-phase
//!   Robbins-Monro algorithm of Snijders (2001) as refined in RSiena, with
//!   derivatives by the score-function method of Schweinberger and Snijders
//!   (2007, `schweinberger2007markov`) and standard errors by the delta
//!   method. [`gof`] computes the auxiliary statistics used for goodness of
//!   fit (Lospinoso and Snijders 2019, `lospinoso2019goodness`).
//!
//! # Conventions
//!
//! - Networks are directed and binary. Self-loops are impossible (the
//!   diagonal is ignored) and parallel edges collapse to one tie, as in
//!   RSiena, which works on adjacency matrices. Weights are ignored.
//! - The waves are observations of the same actors in the same order; the
//!   time between consecutive waves is one unit. Period m runs from wave m
//!   to wave m + 1.
//! - Target statistics for period m are the effect statistics of wave
//!   m + 1 evaluated with the covariate values of wave m; a changing
//!   covariate is constant during a period at its value at the start of
//!   the period (RSiena manual, section 5.3).
//! - Actor covariates are centered by subtracting their mean over all
//!   actors (and, for changing covariates, over all waves supplied).
//!   Dyadic covariates are centered by the mean of their off-diagonal
//!   cells. Similarity sim_ij = 1 - |v_i - v_j| / range is centered by its
//!   mean over all ordered pairs of distinct actors (and waves). The
//!   `sameX` effect is not centered. These are RSiena's conventions
//!   (manual, section 5.3.2).
//! - Missing values, structural zeros and ones, and composition change are
//!   not supported: every tie and covariate value must be observed.
//! - Estimation is conditional by default: in each period the simulation
//!   runs until the number of ties that differ from the start of the
//!   period equals the number observed, so the rate parameters are not
//!   part of the Robbins-Monro iteration and are estimated afterwards from
//!   the simulated period lengths (Snijders 2001, section 4.2). With
//!   unconditional estimation the simulation runs for one time unit and
//!   each period's rate is a parameter whose statistic is the observed
//!   number of changed ties.
//! - Random numbers come from [`crate::rng::Stream`]; every simulation has
//!   its own stream, so results are reproducible for a given seed whatever
//!   the thread count (floating-point sums over simulations can differ in
//!   the last bits between thread counts).
//!
//! Not implemented (left for future work): behavior dynamics and
//! co-evolution, multiple networks, non-directed networks, creation and
//! endowment effects, rate effects beyond the constant rate, maximum
//! likelihood and Bayesian estimation, score-type tests, and time
//! heterogeneity tests.

pub mod data;
pub mod effects;
pub mod estimate;
pub mod gof;
pub mod linalg;
pub mod network;
pub mod simulate;

pub use data::{ActorCovariate, DyadicCovariate, Panel};
pub use effects::{Effect, EffectKind, Model};
pub use estimate::{estimate, EstimateOptions, Fit};
pub use simulate::{simulate_panel, PanelSimulation};
