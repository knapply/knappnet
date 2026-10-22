//! Two-mode (bipartite, affiliation) networks: projections, density, and
//! centrality normalized for two-mode data.
//!
//! A two-mode network has two disjoint node sets (for example people and
//! the events they attend) and ties only between the sets. Here the sets
//! are given by a node type: `false` for the first mode and `true` for the
//! second, igraph's convention. Every edge must join the two modes.
//!
//! Sources:
//! - Projections. The one-mode projection onto a mode joins two of its
//!   nodes when they share a partner in the other mode, with weight the
//!   number of shared partners: the off-diagonal of A A^T for the
//!   incidence matrix A (Breiger 1974, `breiger1974duality`). Newman's
//!   collaboration weighting gives each shared partner e the weight
//!   1/(d_e - 1), where d_e is the number of partners of e (Newman 2001,
//!   `newman2001scientific`). For valued data the product weighting sums
//!   w_ie w_je over shared partners, the off-diagonal of W W^T. The
//!   diagonal (each node's number of partners, or sum of squared weights)
//!   can be kept as self-loops, which together with both projections is
//!   the dual-projection approach of Everett and Borgatti (2013,
//!   `everett2013dual`).
//! - Two-mode density is the number of ties divided by n1 n2, the number of
//!   possible ties between the modes (Borgatti and Everett 1997,
//!   `borgatti1997network`; Wasserman and Faust 1994, `wasserman1994social`,
//!   chapter 8).
//! - Centrality normalized for two-mode data (Borgatti and Everett 1997;
//!   Borgatti and Halgin 2011, `borgatti2011analyzing`). Scores are the
//!   usual one-mode scores of the bipartite graph; only the normalization
//!   changes, because a node can be adjacent only to the other mode. For a
//!   node whose own mode has n_o nodes, with n_i nodes in the other mode:
//!   - degree: d / n_i;
//!   - closeness: (n_i + 2(n_o - 1)) / S, where S is the sum of distances
//!     to the other nodes (the smallest possible sum is distance 1 to the
//!     other mode and 2 to the node's own mode);
//!   - betweenness: b / b_max with
//!     b_max = [n_i^2 (s+1)^2 + n_i (s+1)(2t - s - 1) - t(2s - t + 3)] / 2,
//!     s = (n_o - 1) div n_i, t = (n_o - 1) mod n_i. This is the maximum of
//!     Borgatti and Everett (1997, section 5.3), in the single expression
//!     NetworkX uses; for n_o <= n_i it reduces to
//!     n_i(n_i - 1)/2 + (n_o - 1)(n_o - 2)/2 + (n_o - 1)(n_i - 1). The
//!     property tests confirm it is the largest betweenness any bipartite
//!     graph of those sizes attains, by exhaustive search over small sizes;
//!   - eigenvector: the principal eigenvector of the bipartite adjacency
//!     matrix, which is the pair of leading singular vectors of the
//!     incidence matrix (Borgatti and Everett 1997; Borgatti and Halgin
//!     2011). Its two halves have equal length, so dividing the unit
//!     eigenvector by sqrt(1/2), the normalization of Borgatti and Everett
//!     (1997, section 5.4), gives each mode Euclidean length 1. That is
//!     what is returned. Borgatti and Everett's further two-mode
//!     normalization of eigenvector scores has no closed form (they divide
//!     by scores of a conjectured extremal graph) and is not implemented.
//!
//! Conventions:
//! - Every function reads the network as a binary or valued affiliation
//!   matrix and ignores direction: an arc either way is a tie.
//! - Parallel edges: for ties, degree, density, distances, betweenness,
//!   eigenvector scores, and count and Newman projections they collapse
//!   into one tie (shared partners are distinct partners, as in NetworkX's
//!   projections and igraph's `bipartite_projection()`). The product
//!   projection sums their weights, as `as.matrix()` does.
//! - Self-loops cannot occur: a loop joins a node to its own mode.
//! - Weights are ignored by everything except the product projection.
//!   Centralities are defined for binary ties, and distances count steps.
//! - Closeness in networks that are not connected: S sums the distances to
//!   reachable nodes only. With [`Unreachable::WassermanFaust`] the score is
//!   then multiplied by r / (n - 1), the share of the other n - 1 nodes the
//!   node reaches, and a node that reaches nobody scores 0; this is
//!   NetworkX's `bipartite.closeness_centrality(normalized = True)`. With
//!   [`Unreachable::Ignore`] the score is not rescaled (NetworkX with
//!   `normalized = False`), and a node that reaches nobody has NaN
//!   (NetworkX reports 0).
//! - Normalizations whose denominator is 0 or undefined (an empty other
//!   mode, or b_max = 0 when the other mode has a single node) give NaN, as
//!   does density when either mode is empty, and every eigenvector score
//!   when there are no ties.
//! - Projections contain one edge per pair of nodes with at least one
//!   shared partner, listed by increasing first and then second node index,
//!   each pair once with the smaller index first.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::betweenness::{betweenness, BetweennessOptions};
use crate::metrics::closeness::{distance_profiles, ClosenessOptions, Unreachable};
use crate::metrics::degree::{degree, DegreeMode};
use crate::metrics::eigenvector::{
    eigenvector_centrality, EigenScale, EigenvectorOptions, LoopWeight,
};
use crate::metrics::topography::{simplify, ParallelWeight};
use crate::parallel;

/// Checks that `types` has one entry per node and that every edge joins
/// the two modes. Returns the sizes of the first and second modes.
pub fn check_types(g: &Graph, types: &[bool]) -> Result<(usize, usize)> {
    if types.len() != g.node_count() {
        return Err(Error::LengthMismatch {
            what: "node types",
            expected: g.node_count(),
            found: types.len(),
        });
    }
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        if types[s as usize] == types[t as usize] {
            return Err(Error::InvalidArgument(format!(
                "edge {e} joins nodes {s} and {t}, which are in the same mode; \
                 in a two-mode network every edge joins the two modes"
            )));
        }
    }
    let second = types.iter().filter(|&&x| x).count();
    Ok((types.len() - second, second))
}

/// The binary affiliation graph: undirected, parallel edges collapsed.
fn affiliation(g: &Graph) -> Result<Graph> {
    simplify(g, true, ParallelWeight::Max)
}

/// Sizes of each node's own mode and of the other mode.
fn sizes(types: &[bool], first: usize, second: usize, v: usize) -> (f64, f64) {
    if types[v] {
        (second as f64, first as f64)
    } else {
        (first as f64, second as f64)
    }
}

/// Number of ties divided by the number of possible ties between the
/// modes, n1 n2. Parallel edges count once. NaN if a mode is empty.
pub fn twomode_density(g: &Graph, types: &[bool]) -> Result<f64> {
    let (first, second) = check_types(g, types)?;
    let possible = first as f64 * second as f64;
    if possible == 0.0 {
        return Ok(f64::NAN);
    }
    Ok(affiliation(g)?.edge_count() as f64 / possible)
}

/// Number of distinct partners of every node, divided by the size of the
/// other mode when `normalized`.
pub fn twomode_degree(g: &Graph, types: &[bool], normalized: bool) -> Result<Vec<f64>> {
    let (first, second) = check_types(g, types)?;
    let a = affiliation(g)?;
    let mut d = degree(&a, DegreeMode::All, false);
    if normalized {
        for (v, x) in d.iter_mut().enumerate() {
            let (_, other) = sizes(types, first, second, v);
            *x /= other;
        }
    }
    Ok(d)
}

/// Options for [`twomode_closeness`].
#[derive(Debug, Clone, Copy)]
pub struct TwoModeClosenessOptions {
    /// Use Borgatti and Everett's normalization; otherwise 1 / S.
    pub normalized: bool,
    /// `Ignore` or `WassermanFaust`; see the module documentation.
    pub unreachable: Unreachable,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for TwoModeClosenessOptions {
    fn default() -> Self {
        TwoModeClosenessOptions {
            normalized: true,
            unreachable: Unreachable::WassermanFaust,
            threads: 0,
        }
    }
}

/// Closeness normalized for two-mode data. See the module documentation.
pub fn twomode_closeness(
    g: &Graph,
    types: &[bool],
    opts: TwoModeClosenessOptions,
) -> Result<Vec<f64>> {
    let (first, second) = check_types(g, types)?;
    match opts.unreachable {
        Unreachable::Ignore => {}
        Unreachable::WassermanFaust if opts.normalized => {}
        Unreachable::WassermanFaust => {
            return Err(Error::InvalidArgument(
                "the Wasserman-Faust rule applies to normalized closeness only".into(),
            ))
        }
        _ => {
            return Err(Error::InvalidArgument(
                "two-mode closeness supports the ignore and Wasserman-Faust rules \
                 for unreachable nodes"
                    .into(),
            ))
        }
    }
    let a = affiliation(g)?;
    let n = a.node_count();
    let copts = ClosenessOptions {
        mode: DegreeMode::All,
        threads: opts.threads,
        ..Default::default()
    };
    let profiles = distance_profiles(&a, &copts)?;
    Ok(profiles
        .iter()
        .enumerate()
        .map(|(v, p)| {
            if p.reached == 0 {
                return match opts.unreachable {
                    Unreachable::WassermanFaust => 0.0,
                    _ => f64::NAN,
                };
            }
            let (own, other) = sizes(types, first, second, v);
            let mut c = if opts.normalized {
                (other + 2.0 * (own - 1.0)) / p.total
            } else {
                1.0 / p.total
            };
            if opts.unreachable == Unreachable::WassermanFaust {
                c *= p.reached as f64 / (n - 1) as f64;
            }
            c
        })
        .collect())
}

/// Largest betweenness a node can have in a bipartite graph whose own mode
/// has `own` nodes and whose other mode has `other` nodes (Borgatti and
/// Everett 1997). NaN if `other` is 0.
pub fn betweenness_max(own: usize, other: usize) -> f64 {
    if other == 0 || own == 0 {
        return f64::NAN;
    }
    // Exact integer arithmetic; the bracket is always even.
    let (m, s, t) = (
        other as i128,
        ((own - 1) / other) as i128,
        ((own - 1) % other) as i128,
    );
    let twice = m * m * (s + 1) * (s + 1) + m * (s + 1) * (2 * t - s - 1) - t * (2 * s - t + 3);
    twice as f64 / 2.0
}

/// Betweenness normalized for two-mode data. See the module documentation.
pub fn twomode_betweenness(
    g: &Graph,
    types: &[bool],
    normalized: bool,
    threads: usize,
) -> Result<Vec<f64>> {
    let (first, second) = check_types(g, types)?;
    let a = affiliation(g)?;
    let mut b = betweenness(
        &a,
        BetweennessOptions {
            threads,
            ..Default::default()
        },
    )?;
    if normalized {
        let max_first = betweenness_max(first, second);
        let max_second = betweenness_max(second, first);
        for (v, x) in b.iter_mut().enumerate() {
            let max = if types[v] { max_second } else { max_first };
            *x = if max > 0.0 { *x / max } else { f64::NAN };
        }
    }
    Ok(b)
}

/// Result of [`twomode_eigenvector`].
#[derive(Debug, Clone, PartialEq)]
pub struct TwoModeEigenvector {
    /// Scores in node order; each mode has Euclidean length 1.
    pub values: Vec<f64>,
    /// The principal eigenvalue of the bipartite adjacency matrix, which is
    /// the largest singular value of the incidence matrix.
    pub eigenvalue: f64,
    /// Whether the network is connected (the scores are unique only then).
    pub connected: bool,
}

/// Eigenvector centrality for two-mode data: the leading singular vectors
/// of the binary incidence matrix, each mode scaled to length 1. See the
/// module documentation.
pub fn twomode_eigenvector(
    g: &Graph,
    types: &[bool],
    tol: f64,
    max_iter: usize,
    threads: usize,
) -> Result<TwoModeEigenvector> {
    check_types(g, types)?;
    let a = affiliation(g)?;
    let e = eigenvector_centrality(
        &a,
        EigenvectorOptions {
            mode: DegreeMode::All,
            weighted: false,
            scale: EigenScale::Unit,
            loops: LoopWeight::Ignore,
            tol,
            max_iter,
            threads,
        },
    )?;
    let mut values = e.values;
    if a.edge_count() == 0 {
        values.iter_mut().for_each(|x| *x = f64::NAN);
        return Ok(TwoModeEigenvector {
            values,
            eigenvalue: 0.0,
            connected: e.connected,
        });
    }
    // Each mode's half of the principal eigenvector has the same length, so
    // this is the division by sqrt(1/2); computing each norm directly keeps
    // rounding from leaking between modes.
    for mode in [false, true] {
        let norm = values
            .iter()
            .zip(types)
            .filter(|(_, &t)| t == mode)
            .map(|(x, _)| x * x)
            .sum::<f64>()
            .sqrt();
        if norm > 0.0 {
            for (x, &t) in values.iter_mut().zip(types) {
                if t == mode {
                    *x /= norm;
                }
            }
        }
    }
    Ok(TwoModeEigenvector {
        values,
        eigenvalue: e.eigenvalue,
        connected: e.connected,
    })
}

/// How projected ties are weighted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ProjectionWeight {
    /// Number of shared partners.
    #[default]
    Count,
    /// Newman's collaboration weight: sum of 1/(d_e - 1) over shared
    /// partners e.
    Newman,
    /// Sum of w_ie w_je over shared partners (requires weights).
    Product,
}

impl ProjectionWeight {
    /// Decodes the C interface codes: 0 = count, 1 = Newman, 2 = product.
    pub fn from_code(code: i32) -> Result<ProjectionWeight> {
        match code {
            0 => Ok(ProjectionWeight::Count),
            1 => Ok(ProjectionWeight::Newman),
            2 => Ok(ProjectionWeight::Product),
            _ => Err(Error::InvalidArgument(format!(
                "unknown projection weight code {code}"
            ))),
        }
    }
}

/// Options for [`project`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ProjectionOptions {
    /// Project onto the second mode (`true`) or the first (`false`).
    pub second: bool,
    /// How ties are weighted.
    pub weight: ProjectionWeight,
    /// Keep the diagonal as self-loops (not defined for Newman weights).
    pub loops: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// A one-mode projection as an edge list over the original node indices.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Projection {
    /// Smaller endpoint of each edge.
    pub from: Vec<u32>,
    /// Larger endpoint (equal to `from` for a diagonal loop).
    pub to: Vec<u32>,
    /// Edge weights.
    pub weights: Vec<f64>,
}

/// Distinct partners of every node in CSR form, with the summed weight of
/// the edges to each partner.
struct Partners {
    offsets: Vec<usize>,
    nodes: Vec<u32>,
    weights: Vec<f64>,
}

impl Partners {
    fn new(g: &Graph) -> Partners {
        let n = g.node_count();
        let w = g.edge_weights();
        let mut pairs: Vec<(u32, u32, f64)> = Vec::with_capacity(2 * g.edge_count());
        for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
            let x = w.map_or(1.0, |w| w[e]);
            pairs.push((s, t, x));
            pairs.push((t, s, x));
        }
        // Sorting by (node, partner, edge order) keeps weight sums in input
        // order, so they do not depend on the sort.
        pairs.sort_by_key(|p| (p.0, p.1));
        let mut offsets = vec![0usize; n + 1];
        let mut nodes = Vec::with_capacity(pairs.len());
        let mut weights = Vec::with_capacity(pairs.len());
        let mut last: Option<(u32, u32)> = None;
        for (s, t, x) in pairs {
            if last == Some((s, t)) {
                *weights.last_mut().expect("a partner was pushed") += x;
            } else {
                nodes.push(t);
                weights.push(x);
                offsets[s as usize + 1] += 1;
                last = Some((s, t));
            }
        }
        for v in 0..n {
            offsets[v + 1] += offsets[v];
        }
        Partners {
            offsets,
            nodes,
            weights,
        }
    }

    #[inline]
    fn of(&self, v: usize) -> (&[u32], &[f64]) {
        let r = self.offsets[v]..self.offsets[v + 1];
        (&self.nodes[r.clone()], &self.weights[r])
    }
}

/// Per-thread state for [`project`]: a sparse accumulator.
struct Accumulator {
    sums: Vec<f64>,
    seen: Vec<bool>,
    touched: Vec<u32>,
    rows: Vec<(u32, Vec<(u32, f64)>)>,
}

/// One-mode projection of a two-mode network. See the module
/// documentation.
pub fn project(g: &Graph, types: &[bool], opts: ProjectionOptions) -> Result<Projection> {
    check_types(g, types)?;
    if opts.weight == ProjectionWeight::Product && !g.is_weighted() {
        return Err(Error::MissingWeights);
    }
    if opts.weight == ProjectionWeight::Newman && opts.loops {
        return Err(Error::InvalidArgument(
            "Newman's collaboration weights have no diagonal; use loops = false".into(),
        ));
    }
    let n = g.node_count();
    let partners = Partners::new(g);
    let target = opts.second;
    let partials = parallel::for_each_item(
        n,
        opts.threads,
        || Accumulator {
            sums: vec![0.0; n],
            seen: vec![false; n],
            touched: Vec::new(),
            rows: Vec::new(),
        },
        |acc, i| {
            if types[i] != target {
                return;
            }
            let (events, wi) = partners.of(i);
            for (&e, &w_ie) in events.iter().zip(wi) {
                let (others, we) = partners.of(e as usize);
                let share = match opts.weight {
                    ProjectionWeight::Newman if others.len() > 1 => 1.0 / (others.len() - 1) as f64,
                    _ => 1.0,
                };
                // Partners are sorted, so the pairs j >= i are a suffix.
                let start = others.partition_point(|&j| (j as usize) < i);
                for (&j, &w_je) in others[start..].iter().zip(&we[start..]) {
                    let j = j as usize;
                    if j == i && !opts.loops {
                        continue;
                    }
                    let add = match opts.weight {
                        ProjectionWeight::Count => 1.0,
                        ProjectionWeight::Newman => share,
                        ProjectionWeight::Product => w_ie * w_je,
                    };
                    if !acc.seen[j] {
                        acc.seen[j] = true;
                        acc.touched.push(j as u32);
                    }
                    acc.sums[j] += add;
                }
            }
            if acc.touched.is_empty() {
                return;
            }
            acc.touched.sort_unstable();
            let row: Vec<(u32, f64)> = acc
                .touched
                .iter()
                .map(|&j| (j, acc.sums[j as usize]))
                .collect();
            for &j in &acc.touched {
                acc.sums[j as usize] = 0.0;
                acc.seen[j as usize] = false;
            }
            acc.touched.clear();
            acc.rows.push((i as u32, row));
        },
    );
    let mut rows: Vec<(u32, Vec<(u32, f64)>)> =
        partials.into_iter().flat_map(|acc| acc.rows).collect();
    rows.sort_unstable_by_key(|r| r.0);
    let total = rows.iter().map(|r| r.1.len()).sum();
    let mut out = Projection {
        from: Vec::with_capacity(total),
        to: Vec::with_capacity(total),
        weights: Vec::with_capacity(total),
    };
    for (i, row) in rows {
        for (j, w) in row {
            out.from.push(i);
            out.to.push(j);
            out.weights.push(w);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// People 0, 1, 2 and events 3, 4: 0-3, 1-3, 1-4, 2-4, plus a repeat
    /// 1-3 with weights 1, 2, 3, 4, 5.
    fn small() -> (Graph, Vec<bool>) {
        let g = Graph::from_edges(
            5,
            &[0, 1, 1, 4, 3],
            &[3, 3, 4, 2, 1],
            Some(&[1.0, 2.0, 3.0, 4.0, 5.0]),
            false,
        )
        .unwrap();
        (g, vec![false, false, false, true, true])
    }

    #[test]
    fn rejects_edges_within_a_mode() {
        let g = Graph::from_edges(3, &[0], &[1], None, false).unwrap();
        assert!(check_types(&g, &[false, false, true]).is_err());
        assert!(check_types(&g, &[false, true]).is_err());
        assert_eq!(check_types(&g, &[false, true, true]).unwrap(), (1, 2));
    }

    #[test]
    fn projections_count_newman_and_product() {
        let (g, types) = small();
        let opts = |weight, second, loops| ProjectionOptions {
            second,
            weight,
            loops,
            threads: 2,
        };
        let p = project(&g, &types, opts(ProjectionWeight::Count, false, false)).unwrap();
        assert_eq!(
            (p.from, p.to, p.weights),
            (vec![0, 1], vec![1, 2], vec![1.0, 1.0])
        );
        let p = project(&g, &types, opts(ProjectionWeight::Count, true, true)).unwrap();
        assert_eq!(
            (p.from, p.to, p.weights),
            (vec![3, 3, 4], vec![3, 4, 4], vec![2.0, 1.0, 2.0])
        );
        let p = project(&g, &types, opts(ProjectionWeight::Newman, false, false)).unwrap();
        assert_eq!(p.weights, vec![1.0, 1.0]);
        // Edge 1-3 has summed weight 2 + 5 = 7.
        let p = project(&g, &types, opts(ProjectionWeight::Product, false, true)).unwrap();
        assert_eq!(p.from, vec![0, 0, 1, 1, 2]);
        assert_eq!(p.to, vec![0, 1, 1, 2, 2]);
        assert_eq!(p.weights, vec![1.0, 7.0, 49.0 + 9.0, 12.0, 16.0]);
        assert!(project(&g, &types, opts(ProjectionWeight::Newman, false, true)).is_err());
    }

    #[test]
    fn normalizations() {
        let (g, types) = small();
        assert_eq!(twomode_density(&g, &types).unwrap(), 4.0 / 6.0);
        assert_eq!(
            twomode_degree(&g, &types, true).unwrap(),
            vec![0.5, 1.0, 0.5, 2.0 / 3.0, 2.0 / 3.0]
        );
        // Node 1 (first mode, n_o = 3, n_i = 2) reaches 3 and 4 at 1 and
        // 0 and 2 at 2: S = 6, normalized (2 + 4) / 6.
        let c = twomode_closeness(&g, &types, TwoModeClosenessOptions::default()).unwrap();
        assert_eq!(c[1], 1.0);
        let b = twomode_betweenness(&g, &types, false, 1).unwrap();
        assert_eq!(b, vec![0.0, 4.0, 0.0, 3.0, 3.0]);
        // b_max for n_o = 3, n_i = 2: s = 1, t = 0, (16 + 0 - 0) / 2 = 4.
        assert_eq!(betweenness_max(3, 2), 4.0);
        let nb = twomode_betweenness(&g, &types, true, 1).unwrap();
        assert_eq!(nb[1], 1.0);
        assert_eq!(nb[3], 3.0 / betweenness_max(2, 3));
        assert!(betweenness_max(3, 1) == 0.0 && betweenness_max(1, 0).is_nan());
    }

    #[test]
    fn eigenvector_modes_have_unit_length() {
        let (g, types) = small();
        let e = twomode_eigenvector(&g, &types, 1e-14, 100_000, 1).unwrap();
        for mode in [false, true] {
            let s: f64 = e
                .values
                .iter()
                .zip(&types)
                .filter(|(_, &t)| t == mode)
                .map(|(x, _)| x * x)
                .sum();
            assert!((s - 1.0).abs() < 1e-12);
        }
        // A path of five nodes: the largest eigenvalue is 2 cos(pi / 6).
        assert!((e.eigenvalue - 3f64.sqrt()).abs() < 1e-9);
        let empty = Graph::from_edges(2, &[], &[], None, false).unwrap();
        let e = twomode_eigenvector(&empty, &[false, true], 1e-12, 100, 1).unwrap();
        assert!(e.values.iter().all(|x| x.is_nan()));
    }

    #[test]
    fn disconnected_closeness() {
        // 0-2 and an isolated first-mode node 1.
        let g = Graph::from_edges(3, &[0], &[2], None, false).unwrap();
        let types = [false, false, true];
        let wf = twomode_closeness(&g, &types, TwoModeClosenessOptions::default()).unwrap();
        // Node 0: S = 1, (1 + 2) / 1 times r / (n - 1) = 1/2. Node 2:
        // (2 + 0) / 1 times 1/2.
        assert_eq!(wf, vec![1.5, 0.0, 1.0]);
        let ig = TwoModeClosenessOptions {
            unreachable: Unreachable::Ignore,
            ..Default::default()
        };
        let c = twomode_closeness(&g, &types, ig).unwrap();
        assert_eq!(c[0], 3.0);
        assert!(c[1].is_nan());
    }
}
