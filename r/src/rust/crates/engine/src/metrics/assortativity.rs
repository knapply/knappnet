//! Degree and attribute mixing as NetworkX 3.7 computes it: the (x, y)
//! pairs of node values along arcs, mixing dicts and matrices, the
//! assortativity coefficients over a subset of nodes, average neighbour
//! degree, and average degree connectivity.
//!
//! Sources: assortativity and mixing matrices (Newman 2002,
//! `newman2002assortative`; Newman 2003, `newman2003mixing`); average
//! neighbour degree (Barrat et al. 2004, `barrat2004architecture`);
//! average degree connectivity k_nn(k) (Pastor-Satorras, Vázquez, and
//! Vespignani 2001, `pastorsatorras2001dynamical`; Barrat et al. 2004 for
//! the weighted form). The conventions below are those of NetworkX's
//! `assortativity` package (`hagberg2008exploring`), so a backend can
//! return its results exactly; [`crate::metrics::mixing`] holds the
//! igraph-style coefficients over all nodes.
//!
//! Pairs ([`degree_pairs`], NetworkX's `node_degree_xy`): for every node u
//! of the subset and every arc u → v leaving it (every incident edge in an
//! undirected graph, a self-loop once) with v in the subset, the pair
//! (x-degree of u, y-degree of v). Directed graphs take x and y among in-,
//! out-, and total degree (NetworkX offers in and out); undirected graphs
//! use the degree. Degrees count parallel edges and an undirected loop
//! twice; with `weighted` they are strengths. [`attribute_pairs`]
//! (`node_attribute_xy`) pairs the labels of u and v for every arc leaving
//! a node u of the subset, whether or not v is in the subset, as NetworkX
//! does. Parallel edges give one pair each.
//!
//! Coefficients ([`degree_assortativity`], [`attribute_assortativity`],
//! [`numeric_assortativity`]): the nominal coefficient is (Σ e_ii − Σ a_i b_i) / (1 − Σ
//! a_i b_i) over the normalized mixing matrix e of the attribute pairs;
//! the numeric coefficients (degree and numeric attributes) are the
//! Pearson correlation of the pairs, which is what NetworkX's `_numeric_ac`
//! computes from the mixing matrix. With a constant variable the result is
//! NaN; without pairs the numeric coefficients are NaN and the nominal one
//! is 0, as NetworkX returns them (its nominal formula on an empty matrix
//! gives 0).
//!
//! [`average_neighbor_degree`] (NetworkX's `average_neighbor_degree`): for
//! a node with source degree k > 0, the sum over its arcs (successors for
//! source "out", predecessors for "in", both for "in+out"; neighbours in
//! undirected graphs) of the target degree of the other end, each term
//! multiplied by the arc weight when weighted, divided by k (the weighted
//! source degree when weighted); 0 when k = 0. Target degrees are never
//! weighted.
//!
//! [`average_degree_connectivity`] (`average_degree_connectivity`): for
//! every source degree k (unweighted) of the nodes of the subset, the sum
//! of the target degrees of their neighbours (weighted by the arc when
//! weighted) divided by the sum of their source degrees (weighted when
//! weighted). As in NetworkX, the neighbours for source "in+out" of a
//! directed graph are the successors only, and for source "in" the
//! predecessors, with the arc from the predecessor weighting the term.
//!
//! Self-loops count as NetworkX counts them in a `Graph` or `DiGraph`;
//! parallel edges each count (NetworkX's multigraph handling of these two
//! functions reads weights from edge-key dicts and is not reproduced).

use crate::error::{Error, Result};
use crate::graph::{Csr, Graph};

/// Which degree of a node: NetworkX's "in", "out", or "in+out".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeKind {
    /// In-degree.
    In,
    /// Out-degree.
    Out,
    /// In plus out (the degree of an undirected graph).
    Total,
}

impl DegreeKind {
    /// Decodes the C interface codes: 1 = out, 2 = in, 3 = in+out (as the
    /// shared `mode` codes).
    pub fn from_code(code: i32) -> Result<DegreeKind> {
        match code {
            1 => Ok(DegreeKind::Out),
            2 => Ok(DegreeKind::In),
            3 => Ok(DegreeKind::Total),
            _ => Err(Error::InvalidArgument(format!(
                "unknown degree kind code {code}"
            ))),
        }
    }
}

/// Degree of every node of the requested kind (strength if weighted);
/// undirected graphs ignore the kind and count a loop twice.
fn degrees(g: &Graph, kind: DegreeKind, weighted: bool) -> Result<Vec<f64>> {
    let n = g.node_count();
    if weighted && !g.is_weighted() {
        return Err(Error::MissingWeights);
    }
    let row = |c: &Csr, v: usize| -> f64 {
        if weighted {
            c.weights(v).map_or(0.0, |w| w.iter().sum())
        } else {
            c.degree(v) as f64
        }
    };
    let loops = |v: usize| -> f64 {
        let nb = g.out().neighbors(v);
        let w = g.out().weights(v);
        nb.iter()
            .enumerate()
            .filter(|(_, &u)| u as usize == v)
            .map(|(i, _)| {
                if weighted {
                    w.map_or(1.0, |w| w[i])
                } else {
                    1.0
                }
            })
            .sum()
    };
    Ok((0..n)
        .map(|v| {
            if !g.is_directed() {
                return row(g.out(), v) + loops(v);
            }
            match kind {
                DegreeKind::Out => row(g.out(), v),
                DegreeKind::In => row(g.inc(), v),
                DegreeKind::Total => row(g.out(), v) + row(g.inc(), v),
            }
        })
        .collect())
}

fn membership(n: usize, nodes: Option<&[u32]>) -> Result<Vec<bool>> {
    match nodes {
        None => Ok(vec![true; n]),
        Some(list) => {
            let mut m = vec![false; n];
            for &v in list {
                if v as usize >= n {
                    return Err(Error::InvalidArgument(format!(
                        "node {v} is outside 0..{n}"
                    )));
                }
                m[v as usize] = true;
            }
            Ok(m)
        }
    }
}

/// The (x-degree, y-degree) pairs along arcs, as NetworkX's
/// `node_degree_xy`, in order of the source node and then of its arcs.
pub fn degree_pairs(
    g: &Graph,
    x: DegreeKind,
    y: DegreeKind,
    weighted: bool,
    nodes: Option<&[u32]>,
) -> Result<Vec<(f64, f64)>> {
    let n = g.node_count();
    let inside = membership(n, nodes)?;
    let dx = degrees(g, x, weighted)?;
    let dy = degrees(g, y, weighted)?;
    let mut out = Vec::new();
    for u in 0..n {
        if !inside[u] {
            continue;
        }
        for &v in g.out().neighbors(u) {
            if inside[v as usize] {
                out.push((dx[u], dy[v as usize]));
            }
        }
    }
    Ok(out)
}

/// The (label of u, label of v) pairs for every arc u → v leaving a node
/// of the subset, as NetworkX's `node_attribute_xy`.
pub fn attribute_pairs(
    g: &Graph,
    labels: &[f64],
    nodes: Option<&[u32]>,
) -> Result<Vec<(f64, f64)>> {
    let n = g.node_count();
    if labels.len() != n {
        return Err(Error::LengthMismatch {
            what: "node values",
            expected: n,
            found: labels.len(),
        });
    }
    let inside = membership(n, nodes)?;
    let mut out = Vec::new();
    for u in 0..n {
        if inside[u] {
            for &v in g.out().neighbors(u) {
                out.push((labels[u], labels[v as usize]));
            }
        }
    }
    Ok(out)
}

/// Mixing counts: the distinct values met in the pairs (sorted) and the
/// number of pairs for every (x value, y value), row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct MixingCounts {
    /// Distinct values of either end, increasing.
    pub values: Vec<f64>,
    /// `counts[i * k + j]`: pairs with x = `values[i]` and y = `values[j]`.
    pub counts: Vec<f64>,
}

/// Counts pairs by value (NetworkX's `mixing_dict`, as a matrix over the
/// sorted values).
pub fn mixing_counts(pairs: &[(f64, f64)]) -> Result<MixingCounts> {
    if pairs.iter().any(|p| p.0.is_nan() || p.1.is_nan()) {
        return Err(Error::InvalidArgument("node values must not be NaN".into()));
    }
    let mut values: Vec<f64> = pairs.iter().flat_map(|&(a, b)| [a, b]).collect();
    values.sort_by(|a, b| a.total_cmp(b));
    values.dedup();
    let k = values.len();
    let index = |x: f64| {
        values
            .binary_search_by(|v| v.total_cmp(&x))
            .expect("a value")
    };
    let mut counts = vec![0.0; k * k];
    for &(a, b) in pairs {
        counts[index(a) * k + index(b)] += 1.0;
    }
    Ok(MixingCounts { values, counts })
}

/// Pearson correlation of the pairs (NaN without pairs or variance).
pub fn pair_correlation(pairs: &[(f64, f64)]) -> f64 {
    let m = pairs.len() as f64;
    if pairs.is_empty() {
        return f64::NAN;
    }
    let (sx, sy) = pairs
        .iter()
        .fold((0.0, 0.0), |(a, b), &(x, y)| (a + x, b + y));
    let (mx, my) = (sx / m, sy / m);
    let (mut cxy, mut cxx, mut cyy) = (0.0, 0.0, 0.0);
    for &(x, y) in pairs {
        let (dx, dy) = (x - mx, y - my);
        cxy += dx * dy;
        cxx += dx * dx;
        cyy += dy * dy;
    }
    if cxx == 0.0 || cyy == 0.0 {
        return f64::NAN;
    }
    cxy / (cxx.sqrt() * cyy.sqrt())
}

/// Nominal assortativity of the pairs' labels: (Σ e_ii − Σ a_i b_i) / (1 −
/// Σ a_i b_i) over the normalized mixing matrix (NaN when the denominator
/// is 0; 0 without pairs, the value NetworkX computes from an empty
/// matrix).
pub fn pair_nominal(pairs: &[(f64, f64)]) -> Result<f64> {
    let m = mixing_counts(pairs)?;
    let k = m.values.len();
    let total = pairs.len() as f64;
    if total == 0.0 {
        return Ok(0.0);
    }
    let (mut trace, mut ab) = (0.0, 0.0);
    for i in 0..k {
        trace += m.counts[i * k + i] / total;
        let a: f64 = (0..k).map(|j| m.counts[j * k + i]).sum::<f64>() / total;
        let b: f64 = (0..k).map(|j| m.counts[i * k + j]).sum::<f64>() / total;
        ab += a * b;
    }
    if ab == 1.0 {
        return Ok(f64::NAN);
    }
    Ok((trace - ab) / (1.0 - ab))
}

/// Degree assortativity as NetworkX's `degree_assortativity_coefficient`
/// (and `degree_pearson_correlation_coefficient`): the correlation of
/// [`degree_pairs`].
pub fn degree_assortativity(
    g: &Graph,
    x: DegreeKind,
    y: DegreeKind,
    weighted: bool,
    nodes: Option<&[u32]>,
) -> Result<f64> {
    Ok(pair_correlation(&degree_pairs(g, x, y, weighted, nodes)?))
}

/// Nominal assortativity as NetworkX's
/// `attribute_assortativity_coefficient`, over [`attribute_pairs`].
pub fn attribute_assortativity(g: &Graph, labels: &[f64], nodes: Option<&[u32]>) -> Result<f64> {
    pair_nominal(&attribute_pairs(g, labels, nodes)?)
}

/// Numeric assortativity as NetworkX's
/// `numeric_assortativity_coefficient`: the correlation of
/// [`attribute_pairs`], keeping only the pairs whose two values are both
/// values of nodes of the subset (NetworkX builds its mixing matrix over
/// those values and drops the other pairs).
pub fn numeric_assortativity(g: &Graph, values: &[f64], nodes: Option<&[u32]>) -> Result<f64> {
    let pairs = attribute_pairs(g, values, nodes)?;
    let inside = membership(g.node_count(), nodes)?;
    let mut seen: Vec<f64> = (0..g.node_count())
        .filter(|&v| inside[v])
        .map(|v| values[v])
        .collect();
    seen.sort_by(|a, b| a.total_cmp(b));
    seen.dedup();
    let known = |x: f64| seen.binary_search_by(|v| v.total_cmp(&x)).is_ok();
    let kept: Vec<(f64, f64)> = pairs
        .into_iter()
        .filter(|p| known(p.0) && known(p.1))
        .collect();
    Ok(pair_correlation(&kept))
}

/// Average neighbour degree of every node of `nodes` (all nodes when
/// `None`), in that order. See the module documentation.
pub fn average_neighbor_degree(
    g: &Graph,
    source: DegreeKind,
    target: DegreeKind,
    weighted: bool,
    nodes: Option<&[u32]>,
) -> Result<Vec<f64>> {
    let n = g.node_count();
    if weighted && !g.is_weighted() {
        return Err(Error::MissingWeights);
    }
    let (source, target) = if g.is_directed() {
        (source, target)
    } else {
        (DegreeKind::Total, DegreeKind::Total)
    };
    let t_deg = degrees(g, target, false)?;
    let s_deg = degrees(g, source, weighted)?;
    let all: Vec<u32> = (0..n as u32).collect();
    let list = nodes.unwrap_or(&all);
    membership(n, Some(list))?;
    let sum_over = |c: &Csr, v: usize| -> f64 {
        let w = c.weights(v);
        c.neighbors(v)
            .iter()
            .enumerate()
            .map(|(i, &u)| {
                let x = if weighted {
                    w.map_or(1.0, |w| w[i])
                } else {
                    1.0
                };
                x * t_deg[u as usize]
            })
            .sum()
    };
    Ok(list
        .iter()
        .map(|&v| {
            let v = v as usize;
            let k = s_deg[v];
            if k == 0.0 {
                return 0.0;
            }
            let total = if !g.is_directed() {
                sum_over(g.out(), v)
            } else {
                match source {
                    DegreeKind::Out => sum_over(g.out(), v),
                    DegreeKind::In => sum_over(g.inc(), v),
                    DegreeKind::Total => sum_over(g.out(), v) + sum_over(g.inc(), v),
                }
            };
            total / k
        })
        .collect())
}

/// Average degree connectivity: the distinct source degrees k of the
/// nodes of the subset (increasing) and k_nn(k) for each. See the module
/// documentation.
pub fn average_degree_connectivity(
    g: &Graph,
    source: DegreeKind,
    target: DegreeKind,
    weighted: bool,
    nodes: Option<&[u32]>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let n = g.node_count();
    if weighted && !g.is_weighted() {
        return Err(Error::MissingWeights);
    }
    let (source, target) = if g.is_directed() {
        (source, target)
    } else {
        (DegreeKind::Total, DegreeKind::Total)
    };
    let inside = membership(n, nodes)?;
    let k_of = degrees(g, source, false)?;
    let kw = degrees(g, source, weighted)?;
    let t_deg = degrees(g, target, false)?;
    // Neighbours: predecessors for "in", successors otherwise (NetworkX's
    // DiGraph.neighbors for "in+out").
    let csr = if g.is_directed() && source == DegreeKind::In {
        g.inc()
    } else {
        g.out()
    };
    let mut acc: Vec<(f64, f64, f64)> = Vec::new();
    for v in 0..n {
        if !inside[v] {
            continue;
        }
        let w = csr.weights(v);
        let s: f64 = csr
            .neighbors(v)
            .iter()
            .enumerate()
            .map(|(i, &u)| {
                let x = if weighted {
                    w.map_or(1.0, |w| w[i])
                } else {
                    1.0
                };
                x * t_deg[u as usize]
            })
            .sum();
        acc.push((k_of[v], s, kw[v]));
    }
    acc.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (mut keys, mut vals) = (Vec::new(), Vec::new());
    let mut i = 0;
    while i < acc.len() {
        let k = acc[i].0;
        let (mut dsum, mut dnorm) = (0.0, 0.0);
        while i < acc.len() && acc[i].0 == k {
            dsum += acc[i].1;
            dnorm += acc[i].2;
            i += 1;
        }
        keys.push(k);
        vals.push(if dnorm == 0.0 { dsum } else { dsum / dnorm });
    }
    Ok((keys, vals))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neighbour_degrees_of_a_star_with_a_tail() {
        // Star centre 0 with leaves 1, 2, 3; 3 also tied to 4.
        let g = Graph::from_edges(5, &[0, 0, 0, 3], &[1, 2, 3, 4], None, false).unwrap();
        let t = DegreeKind::Total;
        let a = average_neighbor_degree(&g, t, t, false, None).unwrap();
        assert_eq!(a, vec![4.0 / 3.0, 3.0, 3.0, 2.0, 2.0]);
        let (k, v) = average_degree_connectivity(&g, t, t, false, None).unwrap();
        assert_eq!(k, vec![1.0, 2.0, 3.0]);
        // k = 1: nodes 1, 2, 4 with neighbour degrees 3, 3, 2 over 3.
        assert_eq!(v, vec![8.0 / 3.0, 2.0, 4.0 / 3.0]);
        let p = degree_pairs(&g, t, t, false, None).unwrap();
        assert_eq!(p.len(), 8);
        assert!(pair_correlation(&p) < 0.0);
    }
}
