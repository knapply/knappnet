//! Exponential random graph models (ERGMs).
//!
//! An ERGM gives the probability of a network y with n nodes as
//!
//! ```text
//!     P(Y = y) = exp(theta' g(y)) / c(theta),
//! ```
//!
//! where g(y) is a vector of network statistics, theta their coefficients,
//! and c(theta) the sum of exp(theta' g(y)) over every network on the same
//! nodes (Holland and Leinhardt 1981, `holland1981exponential`; Frank and
//! Strauss 1986, `frank1986markov`; Wasserman and Pattison 1996,
//! `wasserman1996logit`; Snijders et al. 2006, `snijders2006new`). The
//! statnet `ergm` package (Hunter et al. 2008, `hunter2008ergm`; Morris,
//! Handcock, and Hunter 2008, `morris2008specification`) is the reference
//! implementation; this module reproduces its term definitions and
//! parameterizations so that models and estimates are directly comparable.
//!
//! Components:
//!
//! - [`state`]: a network that can be toggled one dyad at a time, with
//!   fixed (structural) dyads.
//! - [`terms`]: model terms, their statistics, and their change
//!   statistics.
//! - [`mple`] and [`logit`]: maximum pseudolikelihood by logistic
//!   regression of dyad states on change statistics (Strauss and Ikeda
//!   1990, `strauss1990pseudolikelihood`), which is the exact maximum
//!   likelihood estimate for dyad-independent models.
//! - [`sampler`]: Metropolis-Hastings simulation with tie/no-tie or
//!   random-dyad proposals.
//! - [`mcmle`]: Monte Carlo maximum likelihood (Geyer and Thompson 1992,
//!   `geyer1992constrained`) with the stepping algorithm of Hummel, Hunter,
//!   and Handcock (2012, `hummel2012improving`), standard errors as in
//!   Hunter and Handcock (2006, `hunter2006inference`), and the
//!   log-likelihood by path sampling (Gelman and Meng 1998,
//!   `gelman1998simulating`; Hunter and Handcock 2006).
//! - [`exact`]: exact likelihood, expectations, and maximum likelihood by
//!   enumerating every network, for networks with few free dyads.
//! - [`gof`]: the goodness-of-fit statistics of Hunter, Goodreau, and
//!   Handcock (2008, `hunter2008goodness`).
//! - [`tergm`]: discrete-time temporal ERGMs by pooled pseudolikelihood with
//!   bootstrap confidence intervals (Hanneke, Fu, and Xing 2010,
//!   `hanneke2010discrete`; Desmarais and Cranmer 2012,
//!   `desmarais2012statistical`; Leifeld, Cranmer, and Desmarais 2018,
//!   `leifeld2018temporal`).
//!
//! Conventions, shared by every component and matching `ergm`:
//!
//! - Networks are binary and simple. Parallel edges collapse to one tie,
//!   self-loops are an error, and weights are ignored.
//! - Undirected ties are unordered pairs; directed ties are ordered.
//! - Node attributes are numeric or categorical; categorical levels are
//!   sorted, and `nodefactor` drops the first level by default.
//! - Dyads marked fixed keep their observed state: they are never toggled,
//!   never appear in the pseudolikelihood, and are excluded from the dyad
//!   count D that enters proposals and the null model (whose
//!   log-likelihood is -D log 2).
//! - Geometrically weighted terms use fixed decay parameters.
//! - Randomness comes from [`crate::rng::Stream`]: the sample for a given
//!   seed does not depend on the number of threads.

pub mod exact;
pub mod gof;
pub mod linalg;
pub mod logit;
pub mod mcmle;
pub mod model;
pub mod mple;
pub mod sampler;
pub mod state;
pub mod tergm;
pub mod terms;

pub use model::Model;
pub use sampler::{simulate, Proposal, Sample, SampleOptions};
pub use state::{Constraints, NetState};
pub use terms::{Covariates, DegreeKind, Endpoint, Partners, Term, Terms};
