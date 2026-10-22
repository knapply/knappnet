//! Mixing patterns: mixing matrices, assortativity coefficients, the E-I
//! index, and within-group tie shares.
//!
//! Sources:
//! - Assortativity (Newman 2002, `newman2002assortative`; Newman 2003,
//!   `newman2003mixing`). With e_ij the fraction of arcs from category i to
//!   category j, a_i = sum_j e_ij, and b_j = sum_i e_ij, the nominal
//!   coefficient is r = (sum_i e_ii - sum_i a_i b_i) / (1 - sum_i a_i b_i)
//!   (Newman 2003, eq. 2); its numerator alone is the modularity of the
//!   partition. The coefficient for numeric values is the Pearson
//!   correlation of the values at the two ends of every arc (eq. 21 and
//!   26); for directed networks the source and target may carry different
//!   values. Degree assortativity uses out-degree at the source and
//!   in-degree at the target (igraph's `assortativity_degree()`, NetworkX's
//!   `degree_assortativity_coefficient()` defaults); Newman's "excess"
//!   degrees (degree minus one) give the same correlation.
//! - The E-I index (Krackhardt and Stern 1988, `krackhardt1988informal`):
//!   (E - I) / (E + I), where E counts external ties (between groups) and I
//!   internal ties (within a group), at the level of the network, of each
//!   group (ties of its members), and of each node (its own ties).
//! - Mixing matrices count ties by the categories of their endpoints; the
//!   layout (folded upper triangle for undirected networks, generalized
//!   matrices crossing two different node attributes, `direction`) follows
//!   the R package homophily (Knapp 2019, `knapp2019homophily`).
//! - Within-group tie shares compare the share of a group's ties that stay
//!   inside it with the share expected if ties ignored groups, the
//!   "baseline homophily" of McPherson, Smith-Lovin, and Cook (2001,
//!   `mcpherson2001birds`).
//!
//! # Conventions
//!
//! Arcs: a directed network contributes each arc once (or, ignoring
//! direction, in both directions); an undirected network contributes each
//! edge in both directions. Parallel edges each count. Weights are not used
//! except by weighted degree assortativity, whose values are strengths.
//!
//! | | igraph | NetworkX | homophily |
//! |---|---|---|---|
//! | Undirected self-loop in assortativity | two arcs (each edge is a reciprocal pair) | one arc | two arcs |
//! | Directed nominal assortativity | arcs as directed (`directed = TRUE`) | arcs as directed | direction ignored (symmetrized mixing matrix): igraph's `directed = FALSE` |
//! | Unnormalized assortativity | modularity (nominal), covariance (values) | (none) | (none) |
//! | Undirected mixing matrix | (none) | full symmetric matrix of arc counts | folded: diagonal 2 x internal edges, upper triangle 2 x edges between |
//!
//! The E-I index is defined on binary ties, so weights are ignored, and
//! parallel edges each count as a tie (the R function warns about them, as
//! homophily does). Groups ([`EiCounting::Standard`], the default):
//! I_g is the number of ties among members of g and E_g the number of ties
//! between a member and a non-member (with `direction` `Out` or `In`,
//! only ties sent or received by members). This is the group index of
//! Krackhardt and Stern. [`EiCounting::Homophily`] reproduces homophily's
//! `ei_index(scope = "group")`, which takes I_g from the diagonal of its
//! folded mixing matrix and E_g from the row of its upper triangle, so a
//! tie between two groups counts only for the group listed first. Network
//! and node indices are the same under both, except that homophily counts
//! an undirected self-loop twice at its node. A self-loop is an internal
//! tie. A group or node without ties has an undefined index (NaN).
//!
//! UCINET reports the E-I index at the same three levels; its help page
//! (Network > Cohesion > E-I Index) could not be retrieved, so agreement
//! with UCINET's group and node counting is not claimed.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::{degree, strength, DegreeMode};
use crate::metrics::eigenvector::LoopWeight;

/// Options for the assortativity coefficients.
#[derive(Debug, Clone, Copy)]
pub struct AssortOptions {
    /// Treat arcs of a directed network as directed.
    pub directed: bool,
    /// How many arcs an undirected self-loop contributes: two (igraph,
    /// homophily), one (NetworkX), or none.
    pub loops: LoopWeight,
    /// Return the coefficient (true) or its numerator: modularity for
    /// categories, covariance for values.
    pub normalized: bool,
}

impl Default for AssortOptions {
    fn default() -> Self {
        AssortOptions {
            directed: true,
            loops: LoopWeight::Twice,
            normalized: true,
        }
    }
}

/// Calls `f(source, target)` for every arc as the options read the graph.
fn for_each_arc(g: &Graph, opts: &AssortOptions, mut f: impl FnMut(usize, usize)) {
    let directed = g.is_directed() && opts.directed;
    let loop_copies = match opts.loops {
        LoopWeight::Twice => 2,
        LoopWeight::Once => 1,
        LoopWeight::Ignore => 0,
    };
    for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
        let (s, t) = (s as usize, t as usize);
        if s == t {
            // A directed loop is one arc; ignoring direction it is an
            // undirected loop.
            let copies = if directed { 1 } else { loop_copies };
            for _ in 0..copies {
                f(s, s);
            }
        } else if directed {
            f(s, t);
        } else {
            f(s, t);
            f(t, s);
        }
    }
}

fn check_len<T>(v: &[T], n: usize, what: &'static str) -> Result<()> {
    if v.len() != n {
        return Err(Error::LengthMismatch {
            what,
            expected: n,
            found: v.len(),
        });
    }
    Ok(())
}

/// Checks category codes 0..k (or -1 for missing when `missing` is true).
fn check_codes(codes: &[i64], k: usize, missing: bool) -> Result<()> {
    if let Some(&c) = codes
        .iter()
        .find(|&&c| !(c >= 0 && (c as u64) < k as u64) && !(missing && c == -1))
    {
        return Err(Error::InvalidArgument(format!(
            "category code {c} is outside 0..{k}"
        )));
    }
    Ok(())
}

/// Nominal assortativity of categories `types` (codes 0..k). See the
/// module documentation.
pub fn assortativity_nominal(
    g: &Graph,
    types: &[i64],
    k: usize,
    opts: AssortOptions,
) -> Result<f64> {
    let n = g.node_count();
    check_len(types, n, "types")?;
    check_codes(types, k, false)?;
    let mut e = vec![0.0f64; k * k];
    let mut total = 0.0;
    for_each_arc(g, &opts, |s, t| {
        e[types[s] as usize * k + types[t] as usize] += 1.0;
        total += 1.0;
    });
    if total == 0.0 {
        return Ok(f64::NAN);
    }
    let (mut trace, mut ab) = (0.0, 0.0);
    for i in 0..k {
        trace += e[i * k + i] / total;
        let a: f64 = (0..k).map(|j| e[i * k + j]).sum::<f64>() / total;
        let b: f64 = (0..k).map(|j| e[j * k + i]).sum::<f64>() / total;
        ab += a * b;
    }
    if !opts.normalized {
        return Ok(trace - ab);
    }
    if ab >= 1.0 {
        return Ok(f64::NAN);
    }
    Ok((trace - ab) / (1.0 - ab))
}

/// Assortativity of numeric node values: the Pearson correlation of
/// `values` at arc sources and `values_in` (default `values`) at arc
/// targets. See the module documentation.
pub fn assortativity_values(
    g: &Graph,
    values: &[f64],
    values_in: Option<&[f64]>,
    opts: AssortOptions,
) -> Result<f64> {
    let n = g.node_count();
    check_len(values, n, "values")?;
    let target = values_in.unwrap_or(values);
    check_len(target, n, "values_in")?;
    if values.iter().chain(target).any(|v| !v.is_finite()) {
        return Err(Error::InvalidArgument(
            "assortativity needs finite values for every node".into(),
        ));
    }
    // Two passes (means, then centered moments) for accuracy.
    let (mut count, mut sx, mut sy) = (0.0f64, 0.0f64, 0.0f64);
    for_each_arc(g, &opts, |s, t| {
        count += 1.0;
        sx += values[s];
        sy += target[t];
    });
    if count == 0.0 {
        return Ok(f64::NAN);
    }
    let (mx, my) = (sx / count, sy / count);
    let (mut cxy, mut cxx, mut cyy) = (0.0f64, 0.0f64, 0.0f64);
    for_each_arc(g, &opts, |s, t| {
        let dx = values[s] - mx;
        let dy = target[t] - my;
        cxy += dx * dy;
        cxx += dx * dx;
        cyy += dy * dy;
    });
    if !opts.normalized {
        return Ok(cxy / count);
    }
    if cxx == 0.0 || cyy == 0.0 {
        return Ok(f64::NAN);
    }
    Ok(cxy / (cxx.sqrt() * cyy.sqrt()))
}

/// Degree assortativity: out-degree at the source against in-degree at
/// the target, or total degree when direction is ignored. With `weighted`
/// the values are strengths. Degrees count every edge end, so an
/// undirected self-loop adds 2 to its node's degree, whatever
/// `opts.loops` says about the arcs it contributes.
pub fn assortativity_degree(g: &Graph, weighted: bool, opts: AssortOptions) -> Result<f64> {
    let value = |mode: DegreeMode| -> Result<Vec<f64>> {
        if weighted {
            strength(g, mode, true)
        } else {
            Ok(degree(g, mode, true))
        }
    };
    if g.is_directed() && opts.directed {
        let out = value(DegreeMode::Out)?;
        let inc = value(DegreeMode::In)?;
        assortativity_values(g, &out, Some(&inc), opts)
    } else {
        let all = value(DegreeMode::All)?;
        assortativity_values(g, &all, None, opts)
    }
}

/// A mixing matrix: counts of arcs by the categories of their endpoints.
#[derive(Debug, Clone, PartialEq)]
pub struct Mixing {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Counts, row-major.
    pub counts: Vec<f64>,
    /// Whether rows hold the second dimension and columns the first
    /// (homophily's layout for `direction = "in"` in directed networks).
    pub transposed: bool,
}

/// The mixing matrix of node categories `dim1` (codes 0..k1) by `dim2`
/// (codes 0..k2), with -1 marking a missing category (arcs touching it are
/// not counted). `same` says that the two dimensions are one and the same
/// attribute. Layout, following homophily's `as_mixing_matrix()`:
/// - Undirected: every edge counts in both directions (a loop twice); when
///   `same`, the matrix is folded, entry (i, j) for i < j holding the
///   count for (i, j) plus that for (j, i), with zeros below the diagonal.
/// - Directed, `same`: rows are source categories and columns target
///   categories; `direction` `In` transposes the result.
/// - Directed, different dimensions: `Out` counts dim1 of the source by
///   dim2 of the target; `In` counts dim2 of the source (rows) by dim1 of
///   the target (columns); `All` counts both arc directions, dim1 of
///   either endpoint by dim2 of the other.
pub fn mixing_matrix(
    g: &Graph,
    dim1: &[i64],
    k1: usize,
    dim2: &[i64],
    k2: usize,
    same: bool,
    direction: DegreeMode,
) -> Result<Mixing> {
    let n = g.node_count();
    check_len(dim1, n, "dim1")?;
    check_len(dim2, n, "dim2")?;
    check_codes(dim1, k1, true)?;
    check_codes(dim2, k2, true)?;
    if same && (k1 != k2 || dim1 != dim2) {
        return Err(Error::InvalidArgument(
            "`same` requires identical dimensions".into(),
        ));
    }
    let mut counts = vec![0.0f64; k1 * k2];
    let mut add = |a: usize, b: usize| {
        let (r, c) = (dim1[a], dim2[b]);
        if r >= 0 && c >= 0 {
            counts[r as usize * k2 + c as usize] += 1.0;
        }
    };
    let both = !g.is_directed() || (!same && direction == DegreeMode::All);
    let reverse = g.is_directed() && !same && direction == DegreeMode::In;
    for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
        let (s, t) = (s as usize, t as usize);
        if both {
            add(s, t);
            add(t, s);
        } else if reverse {
            add(t, s);
        } else {
            add(s, t);
        }
    }
    if !g.is_directed() && same {
        for i in 0..k1 {
            for j in (i + 1)..k1 {
                counts[i * k1 + j] += counts[j * k1 + i];
                counts[j * k1 + i] = 0.0;
            }
        }
    }
    if g.is_directed() && direction == DegreeMode::In {
        let mut t = vec![0.0; k1 * k2];
        for i in 0..k1 {
            for j in 0..k2 {
                t[j * k1 + i] = counts[i * k2 + j];
            }
        }
        return Ok(Mixing {
            rows: k2,
            cols: k1,
            counts: t,
            transposed: true,
        });
    }
    Ok(Mixing {
        rows: k1,
        cols: k2,
        counts,
        transposed: false,
    })
}

/// The level at which the E-I index is reported.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EiScope {
    /// One value for the network.
    #[default]
    Global,
    /// One value per group.
    Group,
    /// One value per node.
    Node,
}

impl EiScope {
    /// Decodes the C interface codes: 0 = global, 1 = group, 2 = node.
    pub fn from_code(code: i32) -> Result<EiScope> {
        match code {
            0 => Ok(EiScope::Global),
            1 => Ok(EiScope::Group),
            2 => Ok(EiScope::Node),
            _ => Err(Error::InvalidArgument(format!(
                "unknown E-I scope code {code}"
            ))),
        }
    }
}

/// How ties are counted for the E-I index of groups (and of nodes with
/// undirected self-loops). See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EiCounting {
    /// Krackhardt and Stern (1988).
    #[default]
    Standard,
    /// The R package homophily.
    Homophily,
}

impl EiCounting {
    /// Decodes the C interface codes: 0 = standard, 1 = homophily.
    pub fn from_code(code: i32) -> Result<EiCounting> {
        match code {
            0 => Ok(EiCounting::Standard),
            1 => Ok(EiCounting::Homophily),
            _ => Err(Error::InvalidArgument(format!(
                "unknown E-I counting code {code}"
            ))),
        }
    }
}

fn ratio(e: f64, i: f64) -> f64 {
    if e + i == 0.0 {
        f64::NAN
    } else {
        (e - i) / (e + i)
    }
}

/// External and internal tie counts behind the E-I index, one pair per
/// group (scope `Group`), per node (`Node`), or for the network
/// (`Global`). `groups` holds codes 0..k. `direction` selects ties sent
/// (`Out`), received (`In`), or both (`All`) for groups and nodes of
/// directed networks.
pub fn ei_counts(
    g: &Graph,
    groups: &[i64],
    k: usize,
    scope: EiScope,
    direction: DegreeMode,
    counting: EiCounting,
) -> Result<Vec<(f64, f64)>> {
    let n = g.node_count();
    check_len(groups, n, "groups")?;
    check_codes(groups, k, false)?;
    let directed = g.is_directed();
    let sent = !directed || direction != DegreeMode::In;
    let received = !directed || direction != DegreeMode::Out;
    let edges = g.edge_sources().iter().zip(g.edge_targets());
    let mut out = match scope {
        EiScope::Global => vec![(0.0, 0.0)],
        EiScope::Group => vec![(0.0, 0.0); k],
        EiScope::Node => vec![(0.0, 0.0); n],
    };
    if scope == EiScope::Group && counting == EiCounting::Homophily {
        return Ok(homophily_groups(g, groups, k, direction));
    }
    for (&s, &t) in edges {
        let (s, t) = (s as usize, t as usize);
        let (gs, gt) = (groups[s] as usize, groups[t] as usize);
        let internal = gs == gt;
        match scope {
            EiScope::Global => {
                if internal {
                    out[0].1 += 1.0;
                } else {
                    out[0].0 += 1.0;
                }
            }
            EiScope::Group => {
                if internal {
                    out[gs].1 += 1.0;
                } else {
                    if sent {
                        out[gs].0 += 1.0;
                    }
                    if received {
                        out[gt].0 += 1.0;
                    }
                }
            }
            EiScope::Node => {
                let bump = |x: &mut (f64, f64), by: f64| {
                    if internal {
                        x.1 += by;
                    } else {
                        x.0 += by;
                    }
                };
                if s == t {
                    // One tie; homophily counts an undirected loop (and,
                    // with direction "all", a directed one) twice.
                    let twice = counting == EiCounting::Homophily && sent && received;
                    bump(&mut out[s], if twice { 2.0 } else { 1.0 });
                    continue;
                }
                if sent {
                    bump(&mut out[s], 1.0);
                }
                if received {
                    bump(&mut out[t], 1.0);
                }
            }
        }
    }
    Ok(out)
}

/// homophily's group counts: I_g from the diagonal and E_g from row g of
/// the upper triangle of its mixing matrix.
fn homophily_groups(g: &Graph, groups: &[i64], k: usize, direction: DegreeMode) -> Vec<(f64, f64)> {
    let m = mixing_matrix(g, groups, k, groups, k, true, direction)
        .expect("codes were checked by the caller");
    (0..k)
        .map(|i| {
            let e: f64 = ((i + 1)..k).map(|j| m.counts[i * k + j]).sum();
            (e, m.counts[i * k + i])
        })
        .collect()
}

/// The E-I index (E - I) / (E + I) at the given scope; NaN where there are
/// no ties. See [`ei_counts`] and the module documentation.
pub fn ei_index(
    g: &Graph,
    groups: &[i64],
    k: usize,
    scope: EiScope,
    direction: DegreeMode,
    counting: EiCounting,
) -> Result<Vec<f64>> {
    Ok(ei_counts(g, groups, k, scope, direction, counting)?
        .into_iter()
        .map(|(e, i)| ratio(e, i))
        .collect())
}

/// Within-group tie shares of one group.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroupShare {
    /// Number of members.
    pub size: usize,
    /// Ties from members to members.
    pub internal: f64,
    /// Ties from members to non-members.
    pub external: f64,
    /// internal / (internal + external); NaN without ties.
    pub observed: f64,
    /// (size - 1) / (n - 1): the share of a member's possible alters that
    /// belong to the group, which is the expected share if ties ignored
    /// groups. NaN when n < 2.
    pub expected: f64,
}

/// Share of each group's ties that stay within the group, against the
/// share expected under random mixing. Ties are counted from the members'
/// side: each tie between two members counts once for each (so `internal`
/// is twice the number of internal ties in an undirected network), and
/// with `direction` `Out` or `In` only ties sent or received by members
/// count. Self-loops are ignored, since a node is not its own alter.
pub fn group_shares(
    g: &Graph,
    groups: &[i64],
    k: usize,
    direction: DegreeMode,
) -> Result<Vec<GroupShare>> {
    let n = g.node_count();
    check_len(groups, n, "groups")?;
    check_codes(groups, k, false)?;
    let directed = g.is_directed();
    let sent = !directed || direction != DegreeMode::In;
    let received = !directed || direction != DegreeMode::Out;
    let mut size = vec![0usize; k];
    for &c in groups {
        size[c as usize] += 1;
    }
    let mut internal = vec![0.0f64; k];
    let mut external = vec![0.0f64; k];
    for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
        let (s, t) = (s as usize, t as usize);
        if s == t {
            continue;
        }
        let (gs, gt) = (groups[s] as usize, groups[t] as usize);
        let tally = |x: &mut Vec<f64>, y: &mut Vec<f64>| {
            if sent {
                if gs == gt {
                    x[gs] += 1.0;
                } else {
                    y[gs] += 1.0;
                }
            }
            if received {
                if gs == gt {
                    x[gt] += 1.0;
                } else {
                    y[gt] += 1.0;
                }
            }
        };
        tally(&mut internal, &mut external);
    }
    Ok((0..k)
        .map(|c| {
            let total = internal[c] + external[c];
            GroupShare {
                size: size[c],
                internal: internal[c],
                external: external[c],
                observed: if total > 0.0 {
                    internal[c] / total
                } else {
                    f64::NAN
                },
                expected: if n >= 2 {
                    (size[c] as f64 - 1.0) / (n as f64 - 1.0)
                } else {
                    f64::NAN
                },
            }
        })
        .collect())
}
