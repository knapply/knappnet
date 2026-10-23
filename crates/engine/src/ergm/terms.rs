//! Model terms: sufficient statistics and their change statistics.
//!
//! A term contributes one or more statistics g_k(y) to an exponential
//! random graph model P(Y = y) = exp(theta' g(y)) / c(theta). The change
//! statistic of dyad (i, j) is
//!
//! ```text
//!     delta_ij(y) = g(y with y_ij = 1) - g(y with y_ij = 0),
//! ```
//!
//! the rest of the network held fixed (Wasserman and Pattison 1996,
//! `wasserman1996logit`; Hunter et al. 2008, `hunter2008ergm`). Every
//! estimation method here works from change statistics: the
//! pseudolikelihood regresses y_ij on delta_ij, and the sampler scores a
//! proposed toggle by theta' delta_ij.
//!
//! Terms, their names, and their parameterizations follow the statnet
//! `ergm` package exactly (Morris, Handcock, and Hunter 2008,
//! `morris2008specification`), so that `ergm_stats()` agrees with
//! `summary(formula)` there:
//!
//! | Term | Statistic |
//! |---|---|
//! | `edges` | number of ties |
//! | `mutual` | number of pairs tied in both directions (directed) |
//! | `nodecov(x)` | sum over ties of x_i + x_j (`nodeocov`: x_i, `nodeicov`: x_j, directed) |
//! | `nodefactor(a)` | for each kept level, number of tie endpoints at nodes with that level (`nodeofactor`, `nodeifactor`: tails or heads only) |
//! | `nodematch(a)` | ties whose endpoints share a level; `diff = TRUE`: one statistic per kept level |
//! | `absdiff(x, pow)` | sum over ties of abs(x_i - x_j)^pow |
//! | `edgecov(w)` | sum over ties of w_ij |
//! | `kstar(k)`, `istar(k)`, `ostar(k)` | number of k-stars (undirected; in- and out-stars directed) |
//! | `triangle` | triangles (undirected); transitive plus cyclic triples (directed) |
//! | `ttriple`, `ctriple` | transitive and cyclic triples (directed) |
//! | `gwesp(a)`, `gwdsp(a)`, `gwnsp(a)` | geometrically weighted edgewise, dyadwise, and non-edgewise shared partners with fixed decay |
//! | `gwdegree(a)`, `gwidegree(a)`, `gwodegree(a)` | geometrically weighted degree with fixed decay |
//! | `degree(k)`, `idegree(k)`, `odegree(k)` | nodes of degree exactly k |
//! | `isolates` | nodes with no ties |
//!
//! Geometrically weighted terms (Snijders et al. 2006, `snijders2006new`;
//! Hunter 2007, `hunter2007curved`) use statnet's `fixed = TRUE` form with
//! decay parameter alpha:
//!
//! ```text
//!     gwesp(alpha) = exp(alpha) * sum_{k >= 1} (1 - (1 - exp(-alpha))^k) EP_k,
//!     gwdegree(alpha) = exp(alpha) * sum_{k >= 1} (1 - (1 - exp(-alpha))^k) D_k,
//! ```
//!
//! where EP_k is the number of ties with exactly k shared partners and D_k
//! the number of nodes of degree k. Shared partners of a directed tie
//! i -> j are the outgoing two-path ("OTP") partners k with i -> k -> j,
//! statnet's default for directed networks; dyadwise statistics count
//! ordered pairs in directed networks and unordered pairs otherwise.
//!
//! `nodefactor` keeps every level but the first in sorted order (statnet's
//! `levels = -1` default); `nodematch(diff = TRUE)` keeps all levels. The
//! caller decides which levels to keep, so other choices are possible.
//!
//! Conventions: ties are binary, self-loops do not exist, and parallel
//! edges have collapsed (see [`super::state`]). Weights are ignored.

use crate::error::{Error, Result};
use crate::matrix::Dense;

use super::state::NetState;

/// Node attributes and dyadic covariates a model may refer to.
#[derive(Debug, Clone, Default)]
pub struct Covariates {
    /// Numeric node attributes, one vector of n values each.
    pub numeric: Vec<Vec<f64>>,
    /// Categorical node attributes as level codes `0..levels`.
    pub categorical: Vec<Vec<u32>>,
    /// Dyadic covariates, n-by-n.
    pub dyadic: Vec<Dense>,
}

impl Covariates {
    fn check(&self, n: usize) -> Result<()> {
        for (k, x) in self.numeric.iter().enumerate() {
            if x.len() != n {
                return Err(Error::LengthMismatch {
                    what: "numeric node attribute",
                    expected: n,
                    found: x.len(),
                });
            }
            if let Some(i) = x.iter().position(|v| !v.is_finite()) {
                return Err(Error::InvalidArgument(format!(
                    "numeric node attribute {k} has a missing or non-finite value at node {i}"
                )));
            }
        }
        for a in &self.categorical {
            if a.len() != n {
                return Err(Error::LengthMismatch {
                    what: "categorical node attribute",
                    expected: n,
                    found: a.len(),
                });
            }
        }
        for (k, w) in self.dyadic.iter().enumerate() {
            if w.rows() != n || w.cols() != n {
                return Err(Error::InvalidArgument(format!(
                    "dyadic covariate {k} is {}-by-{}, but the network has {n} nodes",
                    w.rows(),
                    w.cols()
                )));
            }
            if w.as_slice().iter().any(|v| !v.is_finite()) {
                return Err(Error::InvalidArgument(format!(
                    "dyadic covariate {k} has a missing or non-finite value"
                )));
            }
        }
        Ok(())
    }
}

/// Which endpoint a node-level term reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endpoint {
    /// Both endpoints (`nodecov`, `nodefactor`).
    Both,
    /// The tail of a directed tie (`nodeocov`, `nodeofactor`).
    Tail,
    /// The head of a directed tie (`nodeicov`, `nodeifactor`).
    Head,
}

/// Direction of a degree-based term.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeKind {
    /// Degree (undirected networks).
    Total,
    /// In-degree.
    In,
    /// Out-degree.
    Out,
}

/// Which shared-partner statistic a geometrically weighted term uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Partners {
    /// Edgewise: ties, by their shared partner count.
    Edgewise,
    /// Dyadwise: every dyad, by its shared partner count.
    Dyadwise,
    /// Non-edgewise: dyads without a tie.
    Nonedgewise,
}

/// One model term.
#[derive(Debug, Clone, PartialEq)]
pub enum Term {
    /// `edges`.
    Edges,
    /// `mutual` (directed only).
    Mutual,
    /// `nodecov`, `nodeicov`, `nodeocov`: numeric attribute `attr`.
    NodeCov {
        /// Index into [`Covariates::numeric`].
        attr: usize,
        /// Which endpoint(s) contribute.
        endpoint: Endpoint,
    },
    /// `nodefactor`, `nodeifactor`, `nodeofactor`: one statistic per level
    /// in `levels`.
    NodeFactor {
        /// Index into [`Covariates::categorical`].
        attr: usize,
        /// Which endpoint(s) contribute.
        endpoint: Endpoint,
        /// Level codes that get a statistic, in output order.
        levels: Vec<u32>,
    },
    /// `nodematch`: one statistic (`levels` empty) or one per level.
    NodeMatch {
        /// Index into [`Covariates::categorical`].
        attr: usize,
        /// Level codes that get a statistic (differential homophily), or
        /// empty for uniform homophily.
        levels: Vec<u32>,
    },
    /// `absdiff(attr, pow)`.
    AbsDiff {
        /// Index into [`Covariates::numeric`].
        attr: usize,
        /// Exponent.
        pow: f64,
    },
    /// `edgecov`: dyadic covariate `cov`.
    EdgeCov {
        /// Index into [`Covariates::dyadic`].
        cov: usize,
    },
    /// `kstar(k)`, `istar(k)`, `ostar(k)`.
    Star {
        /// Star size k >= 1.
        k: usize,
        /// Degree direction (`Total` for undirected `kstar`).
        kind: DegreeKind,
    },
    /// `triangle`.
    Triangle,
    /// `ttriple` (directed).
    TransitiveTriple,
    /// `ctriple` (directed).
    CyclicTriple,
    /// `gwesp`, `gwdsp`, `gwnsp` with fixed decay.
    GwPartners {
        /// Decay parameter alpha.
        decay: f64,
        /// Which dyads are counted.
        partners: Partners,
    },
    /// `gwdegree`, `gwidegree`, `gwodegree` with fixed decay.
    GwDegree {
        /// Decay parameter alpha.
        decay: f64,
        /// Degree direction.
        kind: DegreeKind,
    },
    /// `degree(k)`, `idegree(k)`, `odegree(k)`.
    Degree {
        /// The degree counted.
        k: usize,
        /// Degree direction.
        kind: DegreeKind,
    },
    /// `isolates`.
    Isolates,
}

impl Term {
    /// Number of statistics the term contributes.
    pub fn width(&self) -> usize {
        match self {
            Term::NodeFactor { levels, .. } => levels.len(),
            Term::NodeMatch { levels, .. } => levels.len().max(1),
            _ => 1,
        }
    }

    /// Whether the term's change statistic depends only on the dyad and
    /// the covariates, never on other ties (Strauss and Ikeda 1990,
    /// `strauss1990pseudolikelihood`: for such models the pseudolikelihood
    /// is the likelihood).
    pub fn dyad_independent(&self) -> bool {
        matches!(
            self,
            Term::Edges
                | Term::NodeCov { .. }
                | Term::NodeFactor { .. }
                | Term::NodeMatch { .. }
                | Term::AbsDiff { .. }
                | Term::EdgeCov { .. }
        )
    }

    fn check(&self, directed: bool, cov: &Covariates) -> Result<()> {
        let need_directed = |name: &str| {
            if directed {
                Ok(())
            } else {
                Err(Error::InvalidArgument(format!(
                    "term `{name}` is defined only for directed networks"
                )))
            }
        };
        let numeric = |attr: usize| {
            if attr < cov.numeric.len() {
                Ok(())
            } else {
                Err(Error::InvalidArgument(format!(
                    "numeric attribute index {attr} is out of range"
                )))
            }
        };
        // A level that no node has is allowed (its statistic is 0), since
        // factor levels can be declared without being used.
        let categorical = |attr: usize| {
            if attr >= cov.categorical.len() {
                return Err(Error::InvalidArgument(format!(
                    "categorical attribute index {attr} is out of range"
                )));
            }
            Ok(())
        };
        let kind_ok = |kind: DegreeKind, name: &str| match (kind, directed) {
            (DegreeKind::Total, true) => Err(Error::InvalidArgument(format!(
                "term `{name}` is defined only for undirected networks; use its in- or out- form"
            ))),
            (DegreeKind::In, false) | (DegreeKind::Out, false) => Err(Error::InvalidArgument(
                format!("term `{name}` is defined only for directed networks"),
            )),
            _ => Ok(()),
        };
        match self {
            Term::Edges | Term::Isolates | Term::Triangle => Ok(()),
            Term::Mutual => need_directed("mutual"),
            Term::TransitiveTriple => need_directed("ttriple"),
            Term::CyclicTriple => need_directed("ctriple"),
            Term::NodeCov { attr, endpoint } => {
                if *endpoint != Endpoint::Both {
                    need_directed("nodeicov/nodeocov")?;
                }
                numeric(*attr)
            }
            Term::NodeFactor {
                attr,
                endpoint,
                levels,
            } => {
                if *endpoint != Endpoint::Both {
                    need_directed("nodeifactor/nodeofactor")?;
                }
                if levels.is_empty() {
                    return Err(Error::InvalidArgument(
                        "nodefactor needs at least one level".into(),
                    ));
                }
                categorical(*attr)
            }
            Term::NodeMatch { attr, .. } => categorical(*attr),
            Term::AbsDiff { attr, pow } => {
                if !pow.is_finite() || *pow <= 0.0 {
                    return Err(Error::InvalidArgument(
                        "absdiff needs a positive, finite exponent".into(),
                    ));
                }
                numeric(*attr)
            }
            Term::EdgeCov { cov: c } => {
                if *c < cov.dyadic.len() {
                    Ok(())
                } else {
                    Err(Error::InvalidArgument(format!(
                        "dyadic covariate index {c} is out of range"
                    )))
                }
            }
            Term::Star { k, kind } => {
                if *k == 0 {
                    return Err(Error::InvalidArgument(
                        "star size must be at least 1".into(),
                    ));
                }
                kind_ok(*kind, "kstar/istar/ostar")
            }
            Term::GwPartners { decay, .. } => {
                if !decay.is_finite() {
                    return Err(Error::InvalidArgument(
                        "the decay parameter must be finite".into(),
                    ));
                }
                Ok(())
            }
            Term::GwDegree { decay, kind } => {
                if !decay.is_finite() {
                    return Err(Error::InvalidArgument(
                        "the decay parameter must be finite".into(),
                    ));
                }
                kind_ok(*kind, "gwdegree/gwidegree/gwodegree")
            }
            Term::Degree { kind, .. } => kind_ok(*kind, "degree/idegree/odegree"),
        }
    }
}

/// The geometric weighting function of statnet's fixed-decay terms:
/// g(k) = exp(a) (1 - (1 - exp(-a))^k), with g(0) = 0.
#[inline]
pub fn gw(decay: f64, k: usize) -> f64 {
    if k == 0 {
        0.0
    } else {
        decay.exp() * (1.0 - (1.0 - (-decay).exp()).powi(k as i32))
    }
}

/// The binomial coefficient C(n, k) as a double (0 when k > n).
#[inline]
pub fn choose(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut out = 1.0;
    for i in 0..k {
        out = out * (n - i) as f64 / (i + 1) as f64;
    }
    out
}

/// A model: a set of terms over a network with covariates.
#[derive(Debug, Clone)]
pub struct Terms {
    terms: Vec<Term>,
    offsets: Vec<usize>,
    width: usize,
    cov: Covariates,
}

impl Terms {
    /// Checks the terms against the network's direction and covariates.
    pub fn new(terms: Vec<Term>, n: usize, directed: bool, cov: Covariates) -> Result<Terms> {
        if terms.is_empty() {
            return Err(Error::InvalidArgument(
                "a model needs at least one term".into(),
            ));
        }
        cov.check(n)?;
        let mut offsets = Vec::with_capacity(terms.len() + 1);
        let mut width = 0;
        for t in &terms {
            t.check(directed, &cov)?;
            offsets.push(width);
            width += t.width();
        }
        offsets.push(width);
        Ok(Terms {
            terms,
            offsets,
            width,
            cov,
        })
    }

    /// Number of statistics.
    #[inline]
    pub fn width(&self) -> usize {
        self.width
    }

    /// The terms.
    pub fn terms(&self) -> &[Term] {
        &self.terms
    }

    /// The covariates.
    pub fn covariates(&self) -> &Covariates {
        &self.cov
    }

    /// Whether every term is dyad independent.
    pub fn dyad_independent(&self) -> bool {
        self.terms.iter().all(Term::dyad_independent)
    }

    /// The change statistics of dyad (i, j): g(y + ij) - g(y - ij), with
    /// the rest of `y` as it is. Writes `width()` values to `out`.
    pub fn change(&self, y: &NetState, i: usize, j: usize, out: &mut [f64]) {
        debug_assert_eq!(out.len(), self.width);
        let present = y.has(i, j);
        for (t, term) in self.terms.iter().enumerate() {
            let slot = &mut out[self.offsets[t]..self.offsets[t + 1]];
            self.change_term(term, y, i, j, present, slot);
        }
    }

    /// The statistics of the whole network, counted directly.
    pub fn stats(&self, y: &NetState) -> Vec<f64> {
        let mut out = vec![0.0; self.width];
        for (t, term) in self.terms.iter().enumerate() {
            let slot = &mut out[self.offsets[t]..self.offsets[t + 1]];
            self.stat_term(term, y, slot);
        }
        out
    }

    #[inline]
    fn level(&self, attr: usize, v: usize) -> u32 {
        self.cov.categorical[attr][v]
    }

    #[inline]
    fn value(&self, attr: usize, v: usize) -> f64 {
        self.cov.numeric[attr][v]
    }

    /// Degree of `v` in the network without the tie (i, j), for the
    /// direction a term asks about.
    #[inline]
    fn degree_without(
        y: &NetState,
        v: usize,
        kind: DegreeKind,
        i: usize,
        j: usize,
        present: bool,
    ) -> usize {
        let d = match kind {
            DegreeKind::Total => y.degree(v),
            DegreeKind::In => y.in_degree(v),
            DegreeKind::Out => y.out_degree(v),
        };
        if !present {
            return d;
        }
        // The tie i -> j counts toward i's out-degree and j's in-degree
        // (both nodes' degree in an undirected network).
        let counts = match kind {
            DegreeKind::Total => v == i || v == j,
            DegreeKind::Out => v == i,
            DegreeKind::In => v == j,
        };
        d - counts as usize
    }

    #[allow(clippy::too_many_arguments)]
    fn change_term(
        &self,
        term: &Term,
        y: &NetState,
        i: usize,
        j: usize,
        present: bool,
        out: &mut [f64],
    ) {
        let directed = y.is_directed();
        match term {
            Term::Edges => out[0] = 1.0,
            Term::Mutual => out[0] = y.has(j, i) as u8 as f64,
            Term::NodeCov { attr, endpoint } => {
                out[0] = match endpoint {
                    Endpoint::Both => self.value(*attr, i) + self.value(*attr, j),
                    Endpoint::Tail => self.value(*attr, i),
                    Endpoint::Head => self.value(*attr, j),
                }
            }
            Term::NodeFactor {
                attr,
                endpoint,
                levels,
            } => {
                let (li, lj) = (self.level(*attr, i), self.level(*attr, j));
                for (s, &l) in levels.iter().enumerate() {
                    out[s] = match endpoint {
                        Endpoint::Both => (li == l) as u8 as f64 + (lj == l) as u8 as f64,
                        Endpoint::Tail => (li == l) as u8 as f64,
                        Endpoint::Head => (lj == l) as u8 as f64,
                    };
                }
            }
            Term::NodeMatch { attr, levels } => {
                let (li, lj) = (self.level(*attr, i), self.level(*attr, j));
                if levels.is_empty() {
                    out[0] = (li == lj) as u8 as f64;
                } else {
                    for (s, &l) in levels.iter().enumerate() {
                        out[s] = (li == lj && li == l) as u8 as f64;
                    }
                }
            }
            Term::AbsDiff { attr, pow } => {
                out[0] = (self.value(*attr, i) - self.value(*attr, j))
                    .abs()
                    .powf(*pow)
            }
            Term::EdgeCov { cov } => {
                let (a, b) = y.canonical(i, j);
                out[0] = self.cov.dyadic[*cov].get(a as usize, b as usize);
            }
            Term::Star { k, kind } => {
                // A new tie completes a k-star at an endpoint with each
                // choice of k - 1 of that endpoint's other ties.
                out[0] = match kind {
                    DegreeKind::Total => {
                        choose(Self::degree_without(y, i, *kind, i, j, present), k - 1)
                            + choose(Self::degree_without(y, j, *kind, i, j, present), k - 1)
                    }
                    DegreeKind::Out => {
                        choose(Self::degree_without(y, i, *kind, i, j, present), k - 1)
                    }
                    DegreeKind::In => {
                        choose(Self::degree_without(y, j, *kind, i, j, present), k - 1)
                    }
                }
            }
            Term::Triangle => {
                out[0] = if directed {
                    transitive_change(y, i, j) + cyclic_change(y, i, j)
                } else {
                    y.shared_partners_undirected(i, j) as f64
                }
            }
            Term::TransitiveTriple => out[0] = transitive_change(y, i, j),
            Term::CyclicTriple => out[0] = cyclic_change(y, i, j),
            Term::GwPartners { decay, partners } => {
                out[0] = match partners {
                    Partners::Edgewise => gwesp_change(y, i, j, present, *decay),
                    Partners::Dyadwise => gwdsp_change(y, i, j, present, *decay),
                    Partners::Nonedgewise => {
                        gwdsp_change(y, i, j, present, *decay)
                            - gwesp_change(y, i, j, present, *decay)
                    }
                }
            }
            Term::GwDegree { decay, kind } => {
                let step = |v: usize| {
                    let d = Self::degree_without(y, v, *kind, i, j, present);
                    gw(*decay, d + 1) - gw(*decay, d)
                };
                out[0] = match kind {
                    DegreeKind::Total => step(i) + step(j),
                    DegreeKind::Out => step(i),
                    DegreeKind::In => step(j),
                }
            }
            Term::Degree { k, kind } => {
                let step = |v: usize| {
                    let d = Self::degree_without(y, v, *kind, i, j, present);
                    (d + 1 == *k) as u8 as f64 - (d == *k) as u8 as f64
                };
                out[0] = match kind {
                    DegreeKind::Total => step(i) + step(j),
                    DegreeKind::Out => step(i),
                    DegreeKind::In => step(j),
                }
            }
            Term::Isolates => {
                let was =
                    |v: usize| Self::degree_without(y, v, DegreeKind::Total, i, j, present) == 0;
                out[0] = -(was(i) as u8 as f64) - (was(j) as u8 as f64);
            }
        }
    }

    fn stat_term(&self, term: &Term, y: &NetState, out: &mut [f64]) {
        let n = y.node_count();
        let directed = y.is_directed();
        let ties = y.ties();
        match term {
            Term::Edges => out[0] = ties.len() as f64,
            Term::Mutual => {
                out[0] = ties
                    .iter()
                    .filter(|&&(a, b)| a < b && y.has(b as usize, a as usize))
                    .count() as f64
            }
            Term::NodeCov { attr, endpoint } => {
                out[0] = ties
                    .iter()
                    .map(|&(a, b)| match endpoint {
                        Endpoint::Both => {
                            self.value(*attr, a as usize) + self.value(*attr, b as usize)
                        }
                        Endpoint::Tail => self.value(*attr, a as usize),
                        Endpoint::Head => self.value(*attr, b as usize),
                    })
                    .sum()
            }
            Term::NodeFactor {
                attr,
                endpoint,
                levels,
            } => {
                for (s, &l) in levels.iter().enumerate() {
                    out[s] = ties
                        .iter()
                        .map(|&(a, b)| {
                            let ta = (self.level(*attr, a as usize) == l) as u8 as f64;
                            let hb = (self.level(*attr, b as usize) == l) as u8 as f64;
                            match endpoint {
                                Endpoint::Both => ta + hb,
                                Endpoint::Tail => ta,
                                Endpoint::Head => hb,
                            }
                        })
                        .sum();
                }
            }
            Term::NodeMatch { attr, levels } => {
                if levels.is_empty() {
                    out[0] = ties
                        .iter()
                        .filter(|&&(a, b)| {
                            self.level(*attr, a as usize) == self.level(*attr, b as usize)
                        })
                        .count() as f64;
                } else {
                    for (s, &l) in levels.iter().enumerate() {
                        out[s] = ties
                            .iter()
                            .filter(|&&(a, b)| {
                                self.level(*attr, a as usize) == l
                                    && self.level(*attr, b as usize) == l
                            })
                            .count() as f64;
                    }
                }
            }
            Term::AbsDiff { attr, pow } => {
                out[0] = ties
                    .iter()
                    .map(|&(a, b)| {
                        (self.value(*attr, a as usize) - self.value(*attr, b as usize))
                            .abs()
                            .powf(*pow)
                    })
                    .sum()
            }
            Term::EdgeCov { cov } => {
                out[0] = ties
                    .iter()
                    .map(|&(a, b)| self.cov.dyadic[*cov].get(a as usize, b as usize))
                    .sum()
            }
            Term::Star { k, kind } => {
                out[0] = (0..n)
                    .map(|v| {
                        let d = match kind {
                            DegreeKind::Total => y.degree(v),
                            DegreeKind::In => y.in_degree(v),
                            DegreeKind::Out => y.out_degree(v),
                        };
                        choose(d, *k)
                    })
                    .sum()
            }
            Term::Triangle => {
                out[0] = if directed {
                    count_transitive(y) + count_cyclic(y)
                } else {
                    // Each triangle is seen from each of its three edges.
                    ties.iter()
                        .map(|&(a, b)| y.shared_partners_undirected(a as usize, b as usize))
                        .sum::<usize>() as f64
                        / 3.0
                }
            }
            Term::TransitiveTriple => out[0] = count_transitive(y),
            Term::CyclicTriple => out[0] = count_cyclic(y),
            Term::GwPartners { decay, partners } => {
                let sp = |a: usize, b: usize| {
                    if directed {
                        y.two_paths(a, b)
                    } else {
                        y.shared_partners_undirected(a, b)
                    }
                };
                let edgewise: f64 = ties
                    .iter()
                    .map(|&(a, b)| gw(*decay, sp(a as usize, b as usize)))
                    .sum();
                out[0] = match partners {
                    Partners::Edgewise => edgewise,
                    Partners::Dyadwise | Partners::Nonedgewise => {
                        let dyadwise: f64 = y
                            .dyads()
                            .iter()
                            .map(|&(a, b)| gw(*decay, sp(a as usize, b as usize)))
                            .sum();
                        if *partners == Partners::Dyadwise {
                            dyadwise
                        } else {
                            dyadwise - edgewise
                        }
                    }
                }
            }
            Term::GwDegree { decay, kind } => {
                out[0] = (0..n)
                    .map(|v| {
                        gw(
                            *decay,
                            match kind {
                                DegreeKind::Total => y.degree(v),
                                DegreeKind::In => y.in_degree(v),
                                DegreeKind::Out => y.out_degree(v),
                            },
                        )
                    })
                    .sum()
            }
            Term::Degree { k, kind } => {
                out[0] = (0..n)
                    .filter(|&v| {
                        *k == match kind {
                            DegreeKind::Total => y.degree(v),
                            DegreeKind::In => y.in_degree(v),
                            DegreeKind::Out => y.out_degree(v),
                        }
                    })
                    .count() as f64
            }
            Term::Isolates => out[0] = (0..n).filter(|&v| y.degree(v) == 0).count() as f64,
        }
    }
}

/// Transitive triples completed by adding i -> j: the new arc can be the
/// first, second, or shortcut arc of a triple (a -> b, b -> c, a -> c).
fn transitive_change(y: &NetState, i: usize, j: usize) -> f64 {
    let mut count = 0usize;
    // i -> j as a -> b: j -> k and i -> k.
    // i -> j as a -> c: i -> k and k -> j.
    for &k in y.out_neighbors(i) {
        let k = k as usize;
        if k == j {
            continue;
        }
        count += y.has(j, k) as usize + y.has(k, j) as usize;
    }
    // i -> j as b -> c: k -> i and k -> j.
    for &k in y.in_neighbors(i) {
        let k = k as usize;
        if k != j && y.has(k, j) {
            count += 1;
        }
    }
    count as f64
}

/// Cyclic triples completed by adding i -> j: j -> k -> i.
fn cyclic_change(y: &NetState, i: usize, j: usize) -> f64 {
    y.out_neighbors(j)
        .iter()
        .filter(|&&k| k as usize != i && y.has(k as usize, i))
        .count() as f64
}

fn count_transitive(y: &NetState) -> f64 {
    // Each transitive triple a -> b, b -> c, a -> c is counted once at its
    // shortcut arc a -> c, as the number of two-paths a -> b -> c.
    y.ties()
        .iter()
        .map(|&(a, c)| y.two_paths(a as usize, c as usize))
        .sum::<usize>() as f64
}

fn count_cyclic(y: &NetState) -> f64 {
    // Each cycle a -> b -> c -> a is seen from each of its three arcs.
    y.ties()
        .iter()
        .map(|&(a, b)| y.two_paths(b as usize, a as usize))
        .sum::<usize>() as f64
        / 3.0
}

/// Shared partners of the (possibly ordered) pair (a, b) in the network
/// without the tie (i, j).
#[inline]
fn partners_without(y: &NetState, a: usize, b: usize, i: usize, j: usize, present: bool) -> usize {
    let sp = if y.is_directed() {
        y.two_paths(a, b)
    } else {
        y.shared_partners_undirected(a, b)
    };
    if !present {
        return sp;
    }
    // Does the tie (i, j) act as a shared-partner path for (a, b)?
    let uses = if y.is_directed() {
        // Two-path a -> i -> j with j == b, or a -> k -> b with (a, k) = (i, j).
        (a == i && j != b && y.has(j, b)) || (b == j && i != a && y.has(a, i))
    } else {
        // Partner j of (a, b) via a == i, or partner i via a or b == j.
        (a == i && y.has(j, b) && j != b)
            || (b == i && y.has(j, a) && j != a)
            || (a == j && y.has(i, b) && i != b)
            || (b == j && y.has(i, a) && i != a)
    };
    sp - uses as usize
}

/// Change in gwesp from adding i -> j: the new tie's own shared partners,
/// plus every tie that gains i or j as a shared partner.
fn gwesp_change(y: &NetState, i: usize, j: usize, present: bool, decay: f64) -> f64 {
    let own = partners_without(y, i, j, i, j, present);
    let mut delta = gw(decay, own);
    let bump = |a: usize, b: usize| {
        let s = partners_without(y, a, b, i, j, present);
        gw(decay, s + 1) - gw(decay, s)
    };
    if y.is_directed() {
        // Ties i -> b with j -> b gain partner j; ties a -> j with a -> i
        // gain partner i.
        for &b in y.out_neighbors(i) {
            let b = b as usize;
            if b != j && y.has(j, b) {
                delta += bump(i, b);
            }
        }
        for &a in y.in_neighbors(j) {
            let a = a as usize;
            if a != i && y.has(a, i) {
                delta += bump(a, j);
            }
        }
    } else {
        // Every common neighbor u: ties (i, u) and (j, u) each gain one.
        for &u in y.out_neighbors(i) {
            let u = u as usize;
            if u != j && y.has(u, j) {
                delta += bump(i, u) + bump(j, u);
            }
        }
    }
    delta
}

/// Change in gwdsp from adding i -> j: every dyad that gains a two-path
/// through the new tie.
fn gwdsp_change(y: &NetState, i: usize, j: usize, present: bool, decay: f64) -> f64 {
    let mut delta = 0.0;
    let bump = |a: usize, b: usize| {
        let s = partners_without(y, a, b, i, j, present);
        gw(decay, s + 1) - gw(decay, s)
    };
    if y.is_directed() {
        // a -> i -> j gives the ordered pair (a, j) a partner; i -> j -> b
        // gives (i, b) one.
        for &a in y.in_neighbors(i) {
            let a = a as usize;
            if a != j {
                delta += bump(a, j);
            }
        }
        for &b in y.out_neighbors(j) {
            let b = b as usize;
            if b != i {
                delta += bump(i, b);
            }
        }
    } else {
        for &u in y.out_neighbors(i) {
            let u = u as usize;
            if u != j {
                delta += bump(u, j);
            }
        }
        for &u in y.out_neighbors(j) {
            let u = u as usize;
            if u != i {
                delta += bump(u, i);
            }
        }
    }
    delta
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gw_and_choose() {
        assert_eq!(gw(0.5, 0), 0.0);
        assert!((gw(0.0, 3) - 1.0).abs() < 1e-12);
        assert!((gw(0.5, 1) - 1.0).abs() < 1e-12);
        assert_eq!(choose(5, 2), 10.0);
        assert_eq!(choose(2, 3), 0.0);
        assert_eq!(choose(4, 0), 1.0);
    }

    #[test]
    fn widths_and_checks() {
        let cov = Covariates {
            categorical: vec![vec![0, 1, 1, 2]],
            ..Default::default()
        };
        let t = Terms::new(
            vec![
                Term::Edges,
                Term::NodeFactor {
                    attr: 0,
                    endpoint: Endpoint::Both,
                    levels: vec![1, 2],
                },
                Term::NodeMatch {
                    attr: 0,
                    levels: vec![],
                },
            ],
            4,
            false,
            cov.clone(),
        )
        .unwrap();
        assert_eq!(t.width(), 4);
        assert!(t.dyad_independent());
        assert!(Terms::new(vec![Term::Mutual], 4, false, cov.clone()).is_err());
        assert!(Terms::new(
            vec![Term::NodeFactor {
                attr: 3,
                endpoint: Endpoint::Both,
                levels: vec![1]
            }],
            4,
            false,
            cov
        )
        .is_err());
    }
}
