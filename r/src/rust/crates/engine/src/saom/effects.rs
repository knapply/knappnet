//! Effects of the network evaluation function.
//!
//! Every effect is defined by its actor statistic s_ik(x), taken verbatim
//! from chapter 13 of the RSiena manual (Ripley et al. 2026,
//! `ripley2026manual`); the short name RSiena uses is given in
//! parentheses. x_ij = 1 if i has a tie to j, x_i+ is i's out-degree, x_+i
//! its in-degree, v_i an actor covariate (centered), w_ij a dyadic
//! covariate (centered), and sums run over all actors other than the ones
//! already named (x_jj = 0 throughout).
//!
//! | Effect | s_ik(x) |
//! |---|---|
//! | out-degree (`density`) | x_i+ |
//! | reciprocity (`recip`) | sum_j x_ij x_ji |
//! | transitive triplets (`transTrip`) | sum_{j,h} x_ij x_ih x_hj |
//! | transitive ties (`transTies`) | sum_j x_ij max_h (x_ih x_hj) |
//! | transitive reciprocated triplets (`transRecTrip`) | sum_{j,h} x_ij x_ji x_ih x_hj |
//! | 3-cycles (`cycle3`) | sum_{j,h} x_ij x_jh x_hi (the network statistic, the number of 3-cycles, is the sum over actors divided by 3, as RSiena reports it) |
//! | GWESP (`gwespFF`, `gwespFB`, `gwespBF`, `gwespBB`) | sum_j x_ij e^a {1 - (1 - e^-a)^esp_ij} |
//! | in-degree popularity (`inPop`) | sum_j x_ij x_+j |
//! | in-degree popularity sqrt (`inPopSqrt`) | sum_j x_ij sqrt(x_+j) |
//! | out-degree popularity (`outPop`) | sum_j x_ij x_j+ |
//! | out-degree popularity sqrt (`outPopSqrt`) | sum_j x_ij sqrt(x_j+) |
//! | in-degree activity (`inAct`) | x_i+ x_+i |
//! | in-degree activity sqrt (`inActSqrt`) | x_i+ sqrt(x_+i) |
//! | out-degree activity (`outAct`) | x_i+^2 |
//! | out-degree activity sqrt (`outActSqrt`) | x_i+^1.5 |
//! | covariate ego (`egoX`) | v_i x_i+ |
//! | covariate alter (`altX`) | sum_j x_ij v_j |
//! | covariate similarity (`simX`) | sum_j x_ij (sim_ij - mean sim) |
//! | same covariate (`sameX`) | sum_j x_ij 1{v_i = v_j} |
//! | covariate ego x alter (`egoXaltX`) | v_i sum_j x_ij v_j |
//! | dyadic covariate (`X`) | sum_j x_ij w_ij |
//!
//! For the GWESP effects (Hunter 2007, `hunter2007curved`, as adapted to
//! actor-oriented models in the manual) esp_ij counts the actors h with a
//! two-path i -> h -> j (FF), i -> h <- j (FB), i <- h -> j (BF), or
//! i <- h <- j (BB), and a = parameter / 100 (RSiena's default parameter
//! is 69, a = log 2). RSiena implements them as *elementary* effects: the
//! change statistic for the tie i -> j is the weight of that tie alone,
//! e^a {1 - (1 - e^-a)^esp_ij}, and does not include the change the tie
//! makes to the shared-partner counts of i's other ties. The engine does
//! the same, so estimates match RSiena's; the statistic reported as a
//! target is still sum_i s_ik(x).
//!
//! The change statistic Delta_kij(x) = s_ik(x^+ij) - s_ik(x^-ij), with
//! x^+ij and x^-ij the network with and without the tie i -> j, is what a
//! ministep needs (manual, section 6.1.3). It is computed incrementally
//! from per-actor two-path counts, and a randomized test checks it against
//! the difference of the statistics.

use super::data::Panel;
use super::network::Network;
use crate::error::{Error, Result};

/// The kind of an effect, by its RSiena short name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EffectKind {
    /// `density`
    Density,
    /// `recip`
    Recip,
    /// `transTrip`
    TransTrip,
    /// `transTies`
    TransTies,
    /// `transRecTrip`
    TransRecTrip,
    /// `cycle3`
    Cycle3,
    /// `gwespFF`
    GwespFF,
    /// `gwespFB`
    GwespFB,
    /// `gwespBF`
    GwespBF,
    /// `gwespBB`
    GwespBB,
    /// `inPop`
    InPop,
    /// `inPopSqrt`
    InPopSqrt,
    /// `outPop`
    OutPop,
    /// `outPopSqrt`
    OutPopSqrt,
    /// `inAct`
    InAct,
    /// `inActSqrt`
    InActSqrt,
    /// `outAct`
    OutAct,
    /// `outActSqrt`
    OutActSqrt,
    /// `egoX` (actor covariate)
    EgoX,
    /// `altX` (actor covariate)
    AltX,
    /// `simX` (actor covariate)
    SimX,
    /// `sameX` (actor covariate)
    SameX,
    /// `egoXaltX` (actor covariate)
    EgoXAltX,
    /// `X` (dyadic covariate)
    DyadX,
}

impl EffectKind {
    /// All kinds, in the order of their C interface codes.
    pub const ALL: [EffectKind; 24] = [
        EffectKind::Density,
        EffectKind::Recip,
        EffectKind::TransTrip,
        EffectKind::TransTies,
        EffectKind::TransRecTrip,
        EffectKind::Cycle3,
        EffectKind::GwespFF,
        EffectKind::GwespFB,
        EffectKind::GwespBF,
        EffectKind::GwespBB,
        EffectKind::InPop,
        EffectKind::InPopSqrt,
        EffectKind::OutPop,
        EffectKind::OutPopSqrt,
        EffectKind::InAct,
        EffectKind::InActSqrt,
        EffectKind::OutAct,
        EffectKind::OutActSqrt,
        EffectKind::EgoX,
        EffectKind::AltX,
        EffectKind::SimX,
        EffectKind::SameX,
        EffectKind::EgoXAltX,
        EffectKind::DyadX,
    ];

    /// Decodes a C interface code (the position in [`EffectKind::ALL`]).
    pub fn from_code(code: i32) -> Result<EffectKind> {
        usize::try_from(code)
            .ok()
            .and_then(|c| EffectKind::ALL.get(c).copied())
            .ok_or_else(|| Error::InvalidArgument(format!("unknown effect code {code}")))
    }

    /// Decodes an RSiena short name.
    pub fn from_name(name: &str) -> Result<EffectKind> {
        EffectKind::ALL
            .iter()
            .copied()
            .find(|k| k.name() == name)
            .ok_or_else(|| Error::InvalidArgument(format!("unknown effect `{name}`")))
    }

    /// The RSiena short name.
    pub fn name(self) -> &'static str {
        match self {
            EffectKind::Density => "density",
            EffectKind::Recip => "recip",
            EffectKind::TransTrip => "transTrip",
            EffectKind::TransTies => "transTies",
            EffectKind::TransRecTrip => "transRecTrip",
            EffectKind::Cycle3 => "cycle3",
            EffectKind::GwespFF => "gwespFF",
            EffectKind::GwespFB => "gwespFB",
            EffectKind::GwespBF => "gwespBF",
            EffectKind::GwespBB => "gwespBB",
            EffectKind::InPop => "inPop",
            EffectKind::InPopSqrt => "inPopSqrt",
            EffectKind::OutPop => "outPop",
            EffectKind::OutPopSqrt => "outPopSqrt",
            EffectKind::InAct => "inAct",
            EffectKind::InActSqrt => "inActSqrt",
            EffectKind::OutAct => "outAct",
            EffectKind::OutActSqrt => "outActSqrt",
            EffectKind::EgoX => "egoX",
            EffectKind::AltX => "altX",
            EffectKind::SimX => "simX",
            EffectKind::SameX => "sameX",
            EffectKind::EgoXAltX => "egoXaltX",
            EffectKind::DyadX => "X",
        }
    }

    /// Whether the effect needs an actor covariate.
    pub fn needs_actor_covariate(self) -> bool {
        matches!(
            self,
            EffectKind::EgoX
                | EffectKind::AltX
                | EffectKind::SimX
                | EffectKind::SameX
                | EffectKind::EgoXAltX
        )
    }

    /// Whether the effect needs a dyadic covariate.
    pub fn needs_dyadic_covariate(self) -> bool {
        matches!(self, EffectKind::DyadX)
    }

    /// Whether the effect has an internal parameter (the GWESP effects,
    /// a = parameter / 100).
    pub fn has_parameter(self) -> bool {
        matches!(
            self,
            EffectKind::GwespFF | EffectKind::GwespFB | EffectKind::GwespBF | EffectKind::GwespBB
        )
    }
}

/// One effect of the evaluation function.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Effect {
    /// Which effect.
    pub kind: EffectKind,
    /// Index of the actor or dyadic covariate it uses, if any.
    pub covariate: Option<usize>,
    /// Internal effect parameter (GWESP: 100 a; RSiena's default is 69).
    pub parameter: f64,
}

impl Effect {
    /// A structural effect.
    pub fn structural(kind: EffectKind) -> Effect {
        Effect {
            kind,
            covariate: None,
            parameter: if kind.has_parameter() { 69.0 } else { 0.0 },
        }
    }

    /// A covariate effect.
    pub fn with_covariate(kind: EffectKind, covariate: usize) -> Effect {
        Effect {
            kind,
            covariate: Some(covariate),
            parameter: 0.0,
        }
    }

    /// A GWESP effect with internal parameter 100 a.
    pub fn with_parameter(kind: EffectKind, parameter: f64) -> Effect {
        Effect {
            kind,
            covariate: None,
            parameter,
        }
    }

    fn alpha(&self) -> f64 {
        self.parameter / 100.0
    }
}

/// The evaluation function's effects, checked against a panel.
#[derive(Debug, Clone)]
pub struct Model {
    effects: Vec<Effect>,
    /// Which two-path counts the effects need.
    needs: Needs,
}

#[derive(Debug, Clone, Copy, Default)]
struct Needs {
    /// sum_h x_ih x_hj (two-paths i -> h -> j).
    forward: bool,
    /// sum_h x_ih x_jh (out-stars).
    outstar: bool,
    /// sum_h x_hi x_hj (in-stars).
    instar: bool,
    /// sum_h x_hi x_jh (two-paths j -> h -> i).
    backward: bool,
}

impl Model {
    /// Checks the effects against the panel's covariates.
    pub fn new(panel: &Panel, effects: &[Effect]) -> Result<Model> {
        if effects.is_empty() {
            return Err(Error::InvalidArgument(
                "the model needs at least one effect".into(),
            ));
        }
        let mut needs = Needs::default();
        for (k, e) in effects.iter().enumerate() {
            if e.kind.needs_actor_covariate() {
                match e.covariate {
                    Some(c) if c < panel.actor_covariates().len() => {
                        if e.kind == EffectKind::SimX && panel.actor_covariates()[c].range() <= 0.0
                        {
                            return Err(Error::InvalidArgument(format!(
                                "effect {} ({}) needs a covariate that varies between actors",
                                k + 1,
                                e.kind.name()
                            )));
                        }
                    }
                    _ => {
                        return Err(Error::InvalidArgument(format!(
                            "effect {} ({}) needs an actor covariate",
                            k + 1,
                            e.kind.name()
                        )))
                    }
                }
            } else if e.kind.needs_dyadic_covariate() {
                match e.covariate {
                    Some(c) if c < panel.dyadic_covariates().len() => {}
                    _ => {
                        return Err(Error::InvalidArgument(format!(
                            "effect {} ({}) needs a dyadic covariate",
                            k + 1,
                            e.kind.name()
                        )))
                    }
                }
            }
            if e.kind.has_parameter() && !(e.parameter.is_finite() && e.parameter >= 0.0) {
                return Err(Error::InvalidArgument(format!(
                    "effect {} ({}) needs a non-negative internal parameter",
                    k + 1,
                    e.kind.name()
                )));
            }
            match e.kind {
                EffectKind::TransTrip => {
                    needs.forward = true;
                    needs.outstar = true;
                }
                EffectKind::TransTies | EffectKind::TransRecTrip | EffectKind::GwespFF => {
                    needs.forward = true
                }
                EffectKind::GwespFB => needs.outstar = true,
                EffectKind::GwespBF => needs.instar = true,
                EffectKind::Cycle3 | EffectKind::GwespBB => needs.backward = true,
                _ => {}
            }
        }
        Ok(Model {
            effects: effects.to_vec(),
            needs,
        })
    }

    /// The effects.
    pub fn effects(&self) -> &[Effect] {
        &self.effects
    }

    /// Number of effects.
    pub fn len(&self) -> usize {
        self.effects.len()
    }

    /// Whether there are no effects (never, for a constructed model).
    pub fn is_empty(&self) -> bool {
        self.effects.is_empty()
    }

    /// Statistics sum_i s_ik(x) of every effect on `net`, with the
    /// covariates in force during `period`. For `cycle3` the sum counts
    /// every 3-cycle once per actor on it, and RSiena reports the number of
    /// 3-cycles, so the sum is divided by 3 (a constant factor, which leaves
    /// estimates and standard errors unchanged).
    pub fn statistics(&self, panel: &Panel, period: usize, net: &Network) -> Vec<f64> {
        let mut out = vec![0.0; self.effects.len()];
        let mut ctx = ActorContext::new(net.n());
        for i in 0..net.n() {
            ctx.prepare(self, net, i);
            for (k, e) in self.effects.iter().enumerate() {
                out[k] += actor_statistic(e, panel, period, net, &ctx, i);
            }
        }
        for (k, e) in self.effects.iter().enumerate() {
            if e.kind == EffectKind::Cycle3 {
                out[k] /= 3.0;
            }
        }
        out
    }

    /// Statistics s_ik(x) of every effect for one actor.
    pub fn actor_statistics(
        &self,
        panel: &Panel,
        period: usize,
        net: &Network,
        i: usize,
    ) -> Vec<f64> {
        let mut ctx = ActorContext::new(net.n());
        ctx.prepare(self, net, i);
        self.effects
            .iter()
            .map(|e| actor_statistic(e, panel, period, net, &ctx, i))
            .collect()
    }
}

/// Per-actor two-path counts, reused across ministeps.
pub(crate) struct ActorContext {
    actor: usize,
    forward: Vec<u32>,
    outstar: Vec<u32>,
    instar: Vec<u32>,
    backward: Vec<u32>,
}

impl ActorContext {
    pub(crate) fn new(n: usize) -> ActorContext {
        ActorContext {
            actor: usize::MAX,
            forward: vec![0; n],
            outstar: vec![0; n],
            instar: vec![0; n],
            backward: vec![0; n],
        }
    }

    /// Computes the two-path counts from and to actor `i` on `net`.
    pub(crate) fn prepare(&mut self, model: &Model, net: &Network, i: usize) {
        self.actor = i;
        let n = net.n();
        let needs = model.needs;
        if needs.forward {
            self.forward.iter_mut().for_each(|c| *c = 0);
        }
        if needs.outstar {
            self.outstar.iter_mut().for_each(|c| *c = 0);
        }
        if needs.instar {
            self.instar.iter_mut().for_each(|c| *c = 0);
        }
        if needs.backward {
            self.backward.iter_mut().for_each(|c| *c = 0);
        }
        if needs.forward || needs.outstar {
            for h in net.out_neighbors(i) {
                if needs.forward {
                    for (c, &x) in self.forward.iter_mut().zip(net.row(h)) {
                        *c += x as u32;
                    }
                }
                if needs.outstar {
                    for (c, &x) in self.outstar.iter_mut().zip(net.col(h)) {
                        *c += x as u32;
                    }
                }
            }
        }
        if needs.instar || needs.backward {
            for h in net.in_neighbors(i) {
                if needs.instar {
                    for (c, &x) in self.instar.iter_mut().zip(net.row(h)) {
                        *c += x as u32;
                    }
                }
                if needs.backward {
                    for (c, &x) in self.backward.iter_mut().zip(net.col(h)) {
                        *c += x as u32;
                    }
                }
            }
        }
        debug_assert_eq!(self.forward.len(), n);
    }
}

#[inline]
fn gwesp_weight(alpha: f64, esp: u32) -> f64 {
    if esp == 0 {
        return 0.0;
    }
    alpha.exp() * (1.0 - (1.0 - (-alpha).exp()).powi(esp as i32))
}

/// Actor i's statistic s_ik(x) for one effect, given i's two-path counts.
fn actor_statistic(
    e: &Effect,
    panel: &Panel,
    period: usize,
    net: &Network,
    ctx: &ActorContext,
    i: usize,
) -> f64 {
    let row = net.row(i);
    let out_i = net.outdeg(i) as f64;
    let in_i = net.indeg(i) as f64;
    let sum_over_alters = |f: &dyn Fn(usize) -> f64| -> f64 {
        row.iter()
            .enumerate()
            .filter(|(_, &x)| x != 0)
            .map(|(j, _)| f(j))
            .sum()
    };
    match e.kind {
        EffectKind::Density => out_i,
        EffectKind::Recip => sum_over_alters(&|j| net.has(j, i) as u8 as f64),
        EffectKind::TransTrip => sum_over_alters(&|j| ctx.forward[j] as f64),
        EffectKind::TransTies => sum_over_alters(&|j| (ctx.forward[j] > 0) as u8 as f64),
        EffectKind::TransRecTrip => sum_over_alters(&|j| {
            if net.has(j, i) {
                ctx.forward[j] as f64
            } else {
                0.0
            }
        }),
        EffectKind::Cycle3 => sum_over_alters(&|j| ctx.backward[j] as f64),
        EffectKind::GwespFF => sum_over_alters(&|j| gwesp_weight(e.alpha(), ctx.forward[j])),
        EffectKind::GwespFB => sum_over_alters(&|j| gwesp_weight(e.alpha(), ctx.outstar[j])),
        EffectKind::GwespBF => sum_over_alters(&|j| gwesp_weight(e.alpha(), ctx.instar[j])),
        EffectKind::GwespBB => sum_over_alters(&|j| gwesp_weight(e.alpha(), ctx.backward[j])),
        EffectKind::InPop => sum_over_alters(&|j| net.indeg(j) as f64),
        EffectKind::InPopSqrt => sum_over_alters(&|j| (net.indeg(j) as f64).sqrt()),
        EffectKind::OutPop => sum_over_alters(&|j| net.outdeg(j) as f64),
        EffectKind::OutPopSqrt => sum_over_alters(&|j| (net.outdeg(j) as f64).sqrt()),
        EffectKind::InAct => out_i * in_i,
        EffectKind::InActSqrt => out_i * in_i.sqrt(),
        EffectKind::OutAct => out_i * out_i,
        EffectKind::OutActSqrt => out_i * out_i.sqrt(),
        EffectKind::EgoX => {
            let v = panel.actor_covariates()[e.covariate.unwrap()].centered(period);
            v[i] * out_i
        }
        EffectKind::AltX => {
            let v = panel.actor_covariates()[e.covariate.unwrap()].centered(period);
            sum_over_alters(&|j| v[j])
        }
        EffectKind::SimX => {
            let c = &panel.actor_covariates()[e.covariate.unwrap()];
            sum_over_alters(&|j| c.similarity(period, i, j))
        }
        EffectKind::SameX => {
            let v = panel.actor_covariates()[e.covariate.unwrap()].raw(period);
            sum_over_alters(&|j| (v[i] == v[j]) as u8 as f64)
        }
        EffectKind::EgoXAltX => {
            let v = panel.actor_covariates()[e.covariate.unwrap()].centered(period);
            v[i] * sum_over_alters(&|j| v[j])
        }
        EffectKind::DyadX => {
            let w = &panel.dyadic_covariates()[e.covariate.unwrap()];
            sum_over_alters(&|j| w.centered(i, j))
        }
    }
}

/// Change statistics Delta_kij(x) for actor i and every alter j, written
/// to `out[j * k + kk]` (effect kk of alter j); `out[i * k ..]` is 0.
///
/// `ctx` must have been prepared for actor i on `net`. The values do not
/// depend on whether the tie i -> j is present in `net`.
pub(crate) fn change_statistics(
    model: &Model,
    panel: &Panel,
    period: usize,
    net: &Network,
    ctx: &ActorContext,
    i: usize,
    out: &mut [f64],
) {
    debug_assert_eq!(ctx.actor, i);
    let n = net.n();
    let k = model.len();
    debug_assert_eq!(out.len(), n * k);
    let row = net.row(i);
    let col = net.col(i);
    for (kk, e) in model.effects().iter().enumerate() {
        match e.kind {
            EffectKind::Density => {
                for j in 0..n {
                    out[j * k + kk] = 1.0;
                }
            }
            EffectKind::Recip => {
                for j in 0..n {
                    out[j * k + kk] = col[j] as f64;
                }
            }
            EffectKind::TransTrip => {
                // The tie closes forward[j] two-paths and forms an out-star
                // with each of the outstar[j] alters j also points to.
                for j in 0..n {
                    out[j * k + kk] = (ctx.forward[j] + ctx.outstar[j]) as f64;
                }
            }
            EffectKind::TransTies => {
                // Besides making j itself transitively reached, the tie
                // adds the two-path i -> j -> j' for every alter j' of j;
                // only alters j' of i with no other two-path gain from it.
                for j in 0..n {
                    let mut d = (ctx.forward[j] > 0) as u32;
                    let row_j = net.row(j);
                    for jp in 0..n {
                        if row[jp] != 0 && row_j[jp] != 0 && jp != j {
                            // The two-path count without the tie i -> j.
                            let without = ctx.forward[jp] - (row[j] as u32);
                            if without == 0 {
                                d += 1;
                            }
                        }
                    }
                    out[j * k + kk] = d as f64;
                }
            }
            EffectKind::TransRecTrip => {
                // x_ji forward[j] for the closing tie, plus one for every
                // reciprocated alter j' of i that j points to (the tie
                // i -> j then is the i -> h leg of a triplet with h = j).
                for j in 0..n {
                    let mut d = if col[j] != 0 {
                        ctx.forward[j] as f64
                    } else {
                        0.0
                    };
                    let row_j = net.row(j);
                    for jp in 0..n {
                        if jp != j && row[jp] != 0 && col[jp] != 0 && row_j[jp] != 0 {
                            d += 1.0;
                        }
                    }
                    out[j * k + kk] = d;
                }
            }
            EffectKind::Cycle3 => {
                for j in 0..n {
                    out[j * k + kk] = ctx.backward[j] as f64;
                }
            }
            EffectKind::GwespFF => {
                for j in 0..n {
                    out[j * k + kk] = gwesp_weight(e.alpha(), ctx.forward[j]);
                }
            }
            EffectKind::GwespFB => {
                for j in 0..n {
                    out[j * k + kk] = gwesp_weight(e.alpha(), ctx.outstar[j]);
                }
            }
            EffectKind::GwespBF => {
                for j in 0..n {
                    out[j * k + kk] = gwesp_weight(e.alpha(), ctx.instar[j]);
                }
            }
            EffectKind::GwespBB => {
                for j in 0..n {
                    out[j * k + kk] = gwesp_weight(e.alpha(), ctx.backward[j]);
                }
            }
            EffectKind::InPop => {
                // j's in-degree with the tie present.
                for j in 0..n {
                    out[j * k + kk] = (net.indeg(j) - row[j] as u32 + 1) as f64;
                }
            }
            EffectKind::InPopSqrt => {
                for j in 0..n {
                    out[j * k + kk] = ((net.indeg(j) - row[j] as u32 + 1) as f64).sqrt();
                }
            }
            EffectKind::OutPop => {
                for j in 0..n {
                    out[j * k + kk] = net.outdeg(j) as f64;
                }
            }
            EffectKind::OutPopSqrt => {
                for j in 0..n {
                    out[j * k + kk] = (net.outdeg(j) as f64).sqrt();
                }
            }
            EffectKind::InAct => {
                let d = net.indeg(i) as f64;
                for j in 0..n {
                    out[j * k + kk] = d;
                }
            }
            EffectKind::InActSqrt => {
                let d = (net.indeg(i) as f64).sqrt();
                for j in 0..n {
                    out[j * k + kk] = d;
                }
            }
            EffectKind::OutAct => {
                // (o + 1)^2 - o^2 with o the out-degree without the tie.
                for j in 0..n {
                    let o = (net.outdeg(i) - row[j] as u32) as f64;
                    out[j * k + kk] = 2.0 * o + 1.0;
                }
            }
            EffectKind::OutActSqrt => {
                for j in 0..n {
                    let o = (net.outdeg(i) - row[j] as u32) as f64;
                    out[j * k + kk] = (o + 1.0).powf(1.5) - o.powf(1.5);
                }
            }
            EffectKind::EgoX => {
                let v = panel.actor_covariates()[e.covariate.unwrap()].centered(period);
                for j in 0..n {
                    out[j * k + kk] = v[i];
                }
            }
            EffectKind::AltX => {
                let v = panel.actor_covariates()[e.covariate.unwrap()].centered(period);
                for j in 0..n {
                    out[j * k + kk] = v[j];
                }
            }
            EffectKind::SimX => {
                let c = &panel.actor_covariates()[e.covariate.unwrap()];
                for j in 0..n {
                    out[j * k + kk] = c.similarity(period, i, j);
                }
            }
            EffectKind::SameX => {
                let v = panel.actor_covariates()[e.covariate.unwrap()].raw(period);
                for j in 0..n {
                    out[j * k + kk] = (v[i] == v[j]) as u8 as f64;
                }
            }
            EffectKind::EgoXAltX => {
                let v = panel.actor_covariates()[e.covariate.unwrap()].centered(period);
                for j in 0..n {
                    out[j * k + kk] = v[i] * v[j];
                }
            }
            EffectKind::DyadX => {
                let w = &panel.dyadic_covariates()[e.covariate.unwrap()];
                for j in 0..n {
                    out[j * k + kk] = w.centered(i, j);
                }
            }
        }
        out[i * k + kk] = 0.0;
    }
}

/// Change statistics of every effect for actor `i` on wave `wave` of the
/// panel, with the covariates of period `period`: `n * k` values, alter
/// major.
pub fn change_statistics_on_wave(
    model: &Model,
    panel: &Panel,
    period: usize,
    wave: usize,
    i: usize,
) -> Vec<f64> {
    let net = panel.wave(wave);
    let mut ctx = ActorContext::new(net.n());
    ctx.prepare(model, net, i);
    let mut out = vec![0.0; net.n() * model.len()];
    change_statistics(model, panel, period, net, &ctx, i, &mut out);
    out
}

/// Target statistics: for each period m, the statistics of wave m + 1
/// with the covariates of wave m. Returned period-major (`k` values per
/// period).
pub fn targets(model: &Model, panel: &Panel) -> Vec<f64> {
    let mut out = Vec::with_capacity(panel.period_count() * model.len());
    for m in 0..panel.period_count() {
        out.extend(model.statistics(panel, m, panel.wave(m + 1)));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    fn random_panel(seed: u64, n: usize, density: f64) -> Panel {
        let mut rng = Stream::new(seed, 0);
        let mut waves = Vec::new();
        for _ in 0..3 {
            let mut a = vec![0u8; n * n];
            for i in 0..n {
                for j in 0..n {
                    if i != j && rng.uniform() < density {
                        a[i * n + j] = 1;
                    }
                }
            }
            waves.push(a);
        }
        let mut p = Panel::from_adjacency(n, &waves).unwrap();
        let v: Vec<f64> = (0..n).map(|i| (i % 4) as f64).collect();
        p.add_actor_covariate(&v).unwrap();
        let changing: Vec<f64> = (0..3 * n).map(|i| ((i * 7) % 5) as f64).collect();
        p.add_actor_covariate(&changing).unwrap();
        let w: Vec<f64> = (0..n * n).map(|i| ((i * 3) % 7) as f64 / 2.0).collect();
        p.add_dyadic_covariate(&w).unwrap();
        p
    }

    fn all_effects() -> Vec<Effect> {
        let mut es = Vec::new();
        for kind in EffectKind::ALL {
            if kind.needs_actor_covariate() {
                es.push(Effect::with_covariate(kind, 0));
                es.push(Effect::with_covariate(kind, 1));
            } else if kind.needs_dyadic_covariate() {
                es.push(Effect::with_covariate(kind, 0));
            } else if kind.has_parameter() {
                es.push(Effect::with_parameter(kind, 69.0));
                es.push(Effect::with_parameter(kind, 120.0));
            } else {
                es.push(Effect::structural(kind));
            }
        }
        es
    }

    #[test]
    fn hand_checked_statistics() {
        // 0 -> 1, 1 -> 2, 0 -> 2, 2 -> 0: one transitive triplet (0,2 via 1),
        // one 3-cycle (0 -> 1 -> 2 -> 0), one reciprocated pair (0, 2).
        let p = Panel::from_edge_lists(3, &[vec![], vec![(0, 1), (1, 2), (0, 2), (2, 0)]]).unwrap();
        let m = Model::new(
            &p,
            &[
                Effect::structural(EffectKind::Density),
                Effect::structural(EffectKind::Recip),
                Effect::structural(EffectKind::TransTrip),
                Effect::structural(EffectKind::TransTies),
                Effect::structural(EffectKind::Cycle3),
                Effect::structural(EffectKind::InPop),
                Effect::structural(EffectKind::OutAct),
                Effect::with_parameter(EffectKind::GwespFF, 69.0),
            ],
        )
        .unwrap();
        let s = targets(&m, &p);
        // cycle3: each actor's s_i counts the cycle through it once, and the
        // statistic is the number of cycles: 1.
        // inPop: sum_j indeg_j^2 = 1 + 1 + 4 = 6. outAct: 4 + 1 + 1 = 6.
        // gwespFF with a = log 2: tie 0 -> 2 has esp 1, weight 1.
        assert_eq!(s, vec![4.0, 2.0, 1.0, 1.0, 1.0, 6.0, 6.0, 1.0]);
    }

    #[test]
    fn change_statistics_equal_differences_of_statistics() {
        for seed in 0..6 {
            let n = 7 + seed as usize;
            let p = random_panel(seed, n, 0.3);
            let m = Model::new(&p, &all_effects()).unwrap();
            for wave in 0..2 {
                let base = p.wave(wave).clone();
                for i in 0..n {
                    let delta = change_statistics_on_wave(&m, &p, wave, wave, i);
                    for j in 0..n {
                        if j == i {
                            continue;
                        }
                        let mut with = base.clone();
                        with.set(i, j, true);
                        let mut without = base.clone();
                        without.set(i, j, false);
                        let a = m.actor_statistics(&p, wave, &with, i);
                        let b = m.actor_statistics(&p, wave, &without, i);
                        for k in 0..m.len() {
                            if m.effects()[k].kind.has_parameter() {
                                // GWESP is elementary: the change statistic
                                // is the weight of the tie i -> j alone.
                                let tie_weight = {
                                    let esp = match m.effects()[k].kind {
                                        EffectKind::GwespFF => with
                                            .out_neighbors(i)
                                            .filter(|&h| with.has(h, j))
                                            .count(),
                                        EffectKind::GwespFB => with
                                            .out_neighbors(i)
                                            .filter(|&h| with.has(j, h))
                                            .count(),
                                        EffectKind::GwespBF => {
                                            with.in_neighbors(i).filter(|&h| with.has(h, j)).count()
                                        }
                                        _ => {
                                            with.in_neighbors(i).filter(|&h| with.has(j, h)).count()
                                        }
                                    };
                                    gwesp_weight(m.effects()[k].alpha(), esp as u32)
                                };
                                assert!(
                                    (delta[j * m.len() + k] - tie_weight).abs() < 1e-9,
                                    "{} i {i} j {j}",
                                    m.effects()[k].kind.name()
                                );
                                continue;
                            }
                            let expect = a[k] - b[k];
                            assert!(
                                (delta[j * m.len() + k] - expect).abs() < 1e-9,
                                "seed {seed} wave {wave} effect {} i {i} j {j}: {} vs {}",
                                m.effects()[k].kind.name(),
                                delta[j * m.len() + k],
                                expect
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn model_checks_covariates() {
        let p = random_panel(1, 5, 0.3);
        assert!(Model::new(&p, &[Effect::structural(EffectKind::EgoX)]).is_err());
        assert!(Model::new(&p, &[Effect::with_covariate(EffectKind::EgoX, 7)]).is_err());
        assert!(Model::new(&p, &[Effect::with_covariate(EffectKind::DyadX, 1)]).is_err());
        assert!(Model::new(&p, &[Effect::with_parameter(EffectKind::GwespFF, -1.0)]).is_err());
        assert!(Model::new(&p, &[]).is_err());
        assert_eq!(
            EffectKind::from_name("gwespFF").unwrap(),
            EffectKind::GwespFF
        );
        assert_eq!(EffectKind::from_code(23).unwrap(), EffectKind::DyadX);
        assert!(EffectKind::from_code(24).is_err());
    }
}
