//! Trophic levels and coherence, and local and global reaching centrality.
//!
//! Sources:
//! - **Trophic levels** (Levine 1980, `levine1980several`): a basal node
//!   (no incoming arcs) has level 1, and every other node one more than the
//!   weighted average level of the nodes with arcs into it,
//!   s_i = 1 + (1 / k_i) sum_j w_ji s_j, k_i the weighted in-degree. The
//!   levels solve a linear system over the non-basal nodes.
//! - **Trophic coherence** (Johnson, Dominguez-Garcia, Donetti, and Munoz
//!   2014, `johnson2014trophic`): the trophic difference of an arc (u, v) is
//!   s_v - s_u, and the incoherence parameter q is the standard deviation
//!   of the differences over the arcs.
//! - **Reaching centrality** (Mones, Vicsek, and Vicsek 2012,
//!   `mones2012hierarchy`): the local reaching centrality of v is the share
//!   of the other nodes it reaches (directed, unweighted); for weighted
//!   networks it averages, over the nodes v reaches, the mean weight of the
//!   edges along the shortest path to them (weights read as strengths, so
//!   paths use lengths W/w with W the total weight), divided by n - 1. The
//!   global reaching centrality is the mean of max(C_R) - C_R(v) over the
//!   nodes, divided by n - 1 rather than n.
//!
//! Conventions (matching NetworkX):
//! - Trophic measures need a directed network; weights are used as given
//!   (positive), parallel arcs add their weights, and a self-loop counts in
//!   the in-degree and the average. The levels are NaN when there is no
//!   basal node, when some node cannot be reached from a basal node, or
//!   when the system is singular (NetworkX raises an error). The
//!   differences are given per arc, in input order. The incoherence
//!   parameter is the population standard deviation over distinct ordered
//!   pairs of adjacent nodes (NetworkX keys its differences by node pair);
//!   unless `cannibalism`, self-loops are removed first, and the levels
//!   recomputed without them. It is NaN without arcs.
//! - Reaching centrality follows NetworkX's `local_reaching_centrality`
//!   and `global_reaching_centrality`, which use three definitions: on
//!   unweighted directed networks the share of nodes reached; on unweighted
//!   undirected networks the sum of 1/d over the other nodes reached,
//!   divided by n - 1; and on weighted networks the mean-edge-weight
//!   average above, divided by the mean weight W/m when normalized. Among
//!   tied shortest paths the engine takes the one its search settles first
//!   (see [`crate::paths`]), which can average other weights than the path
//!   NetworkX finds. Weights must be positive. The values are NaN for a
//!   network with fewer than two nodes or no edges (NetworkX raises).

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::linalg;
use crate::matrix::Dense;
use crate::metrics::degree::DegreeMode;
use crate::paths::search::{SearchOptions, Searcher};
use crate::paths::{ArcLengths, PathGraph, PathGraphOptions, NO_NODE};

fn require_directed(g: &Graph, what: &str) -> Result<()> {
    if !g.is_directed() {
        return Err(Error::InvalidArgument(format!(
            "{what} are defined for directed networks"
        )));
    }
    Ok(())
}

fn positive_weights(g: &Graph, weighted: bool) -> Result<Option<&[f64]>> {
    if !weighted {
        return Ok(None);
    }
    let w = g.edge_weights().ok_or(Error::MissingWeights)?;
    if let Some(e) = w.iter().position(|&x| x <= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "edge {e} has weight {}; this measure needs positive weights",
            w[e]
        )));
    }
    Ok(Some(w))
}

/// Trophic level of every node (NaN for every node when undefined). See
/// the module documentation.
pub fn trophic_levels(g: &Graph, weighted: bool) -> Result<Vec<f64>> {
    require_directed(g, "trophic levels")?;
    let n = g.node_count();
    let w = positive_weights(g, weighted)?;
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut indeg = vec![0usize; n];
    for &t in tgt {
        indeg[t as usize] += 1;
    }
    let basal: Vec<u32> = (0..n as u32).filter(|&v| indeg[v as usize] == 0).collect();
    if basal.is_empty() {
        return Ok(vec![f64::NAN; n]);
    }
    // Every node must be reachable from a basal node.
    let mut seen = vec![false; n];
    let mut queue = basal.clone();
    for &b in &basal {
        seen[b as usize] = true;
    }
    let mut head = 0;
    while head < queue.len() {
        let v = queue[head] as usize;
        head += 1;
        for &u in g.out().neighbors(v) {
            if !seen[u as usize] {
                seen[u as usize] = true;
                queue.push(u);
            }
        }
    }
    if seen.iter().any(|&s| !s) {
        return Ok(vec![f64::NAN; n]);
    }
    // (I - P) y = 1 over non-basal nodes, P_ij = w_ji / k_i; then
    // s = y + 1 (the basal terms are folded in, as NetworkX does).
    let others: Vec<usize> = (0..n).filter(|&v| indeg[v] > 0).collect();
    let mut pos = vec![usize::MAX; n];
    for (i, &v) in others.iter().enumerate() {
        pos[v] = i;
    }
    let k = others.len();
    let mut strength = vec![0.0; n];
    for e in 0..g.edge_count() {
        strength[tgt[e] as usize] += w.map_or(1.0, |w| w[e]);
    }
    let mut a = Dense::identity(k);
    for e in 0..g.edge_count() {
        let (s, t) = (src[e] as usize, tgt[e] as usize);
        if pos[s] == usize::MAX {
            continue;
        }
        let x = w.map_or(1.0, |w| w[e]) / strength[t];
        let (i, j) = (pos[t], pos[s]);
        a.set(i, j, a.get(i, j) - x);
    }
    let Some(y) = linalg::solve(&a, &vec![1.0; k], 1e-13) else {
        return Ok(vec![f64::NAN; n]);
    };
    let mut levels = vec![1.0; n];
    for (i, &v) in others.iter().enumerate() {
        levels[v] = y[i] + 1.0;
    }
    Ok(levels)
}

/// Trophic difference s_v - s_u of every edge (u, v), in input order.
pub fn trophic_differences(g: &Graph, weighted: bool) -> Result<Vec<f64>> {
    let levels = trophic_levels(g, weighted)?;
    Ok(g.edge_sources()
        .iter()
        .zip(g.edge_targets())
        .map(|(&s, &t)| levels[t as usize] - levels[s as usize])
        .collect())
}

/// Trophic incoherence parameter. See the module documentation.
pub fn trophic_incoherence(g: &Graph, weighted: bool, cannibalism: bool) -> Result<f64> {
    require_directed(g, "trophic levels")?;
    let h;
    let g = if cannibalism {
        g
    } else {
        let keep: Vec<usize> = (0..g.edge_count())
            .filter(|&e| g.edge_sources()[e] != g.edge_targets()[e])
            .collect();
        let s: Vec<i64> = keep.iter().map(|&e| g.edge_sources()[e] as i64).collect();
        let t: Vec<i64> = keep.iter().map(|&e| g.edge_targets()[e] as i64).collect();
        let w: Option<Vec<f64>> = g
            .edge_weights()
            .map(|w| keep.iter().map(|&e| w[e]).collect());
        h = Graph::from_edges(g.node_count(), &s, &t, w.as_deref(), true)?;
        &h
    };
    let diffs = trophic_differences(g, weighted)?;
    let mut pairs: Vec<(u32, u32, f64)> = g
        .edge_sources()
        .iter()
        .zip(g.edge_targets())
        .zip(&diffs)
        .map(|((&s, &t), &d)| (s, t, d))
        .collect();
    pairs.sort_by_key(|p| (p.0, p.1));
    pairs.dedup_by_key(|p| (p.0, p.1));
    if pairs.is_empty() {
        return Ok(f64::NAN);
    }
    let k = pairs.len() as f64;
    let mean = pairs.iter().map(|p| p.2).sum::<f64>() / k;
    let var = pairs
        .iter()
        .map(|p| (p.2 - mean) * (p.2 - mean))
        .sum::<f64>()
        / k;
    Ok(var.sqrt())
}

/// Options for [`local_reaching_centrality`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ReachingOptions {
    /// Use edge weights (as strengths; see the module documentation).
    pub weighted: bool,
    /// Divide by the mean edge weight (weighted only).
    pub normalized: bool,
}

/// Local reaching centrality of every node. See the module documentation.
pub fn local_reaching_centrality(g: &Graph, opts: ReachingOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    let w = positive_weights(g, opts.weighted)?;
    let m = g.edge_count();
    let total_weight: f64 = w.map_or(m as f64, |w| w.iter().sum());
    if n <= 1 || total_weight <= 0.0 {
        return Ok(vec![f64::NAN; n]);
    }
    let lengths: Vec<f64> = match w {
        Some(w) => w.iter().map(|x| total_weight / x).collect(),
        None => vec![1.0; m],
    };
    let pg = PathGraph::new(
        g,
        PathGraphOptions {
            mode: DegreeMode::Out,
            lengths: ArcLengths::PerEdge(&lengths),
            collapse_parallel: false,
        },
    )?;
    let norm = if opts.normalized && w.is_some() {
        total_weight / m as f64
    } else {
        1.0
    };
    let nm1 = (n - 1) as f64;
    let mut ws = Searcher::new(n);
    let mut sum_w = vec![0.0; n];
    let mut hops = vec![0.0; n];
    let mut out = Vec::with_capacity(n);
    for v in 0..n {
        ws.run(&pg, &[v as u32], SearchOptions::default())?;
        let reached = ws.order().len();
        if w.is_none() && g.is_directed() {
            out.push((reached - 1) as f64 / nm1);
            continue;
        }
        let mut total = 0.0;
        for &u in ws.order() {
            let u = u as usize;
            let p = ws.parent(u);
            if p == NO_NODE {
                sum_w[u] = 0.0;
                hops[u] = 0.0;
                continue;
            }
            let e = pg.edge(ws.parent_arc(u) as usize) as usize;
            sum_w[u] = sum_w[p as usize] + w.map_or(1.0, |w| w[e]);
            hops[u] = hops[p as usize] + 1.0;
            total += match w {
                Some(_) => sum_w[u] / hops[u],
                None => 1.0 / hops[u],
            };
        }
        out.push(total / norm / nm1);
    }
    Ok(out)
}

/// Global reaching centrality. See the module documentation.
pub fn global_reaching_centrality(g: &Graph, opts: ReachingOptions) -> Result<f64> {
    let lrc = local_reaching_centrality(g, opts)?;
    let n = lrc.len();
    if n <= 1 || lrc.iter().any(|x| x.is_nan()) {
        return Ok(f64::NAN);
    }
    let max = lrc.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok(lrc.iter().map(|c| max - c).sum::<f64>() / (n - 1) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn food_chain_levels() {
        // 0 -> 1 -> 2 and 0 -> 2: level(2) = 1 + (1 + 2)/2 = 2.5.
        let g = Graph::from_edges(3, &[0, 1, 0], &[1, 2, 2], None, true).unwrap();
        let l = trophic_levels(&g, false).unwrap();
        assert!((l[0] - 1.0).abs() < 1e-12 && (l[1] - 2.0).abs() < 1e-12);
        assert!((l[2] - 2.5).abs() < 1e-12);
        let d = trophic_differences(&g, false).unwrap();
        assert!((d[2] - 1.5).abs() < 1e-12);
        // Differences 1, 0.5, 1.5: standard deviation sqrt(1/6).
        let q = trophic_incoherence(&g, false, false).unwrap();
        assert!((q - (1.0f64 / 6.0).sqrt()).abs() < 1e-12);
        let cyc = Graph::from_edges(2, &[0, 1], &[1, 0], None, true).unwrap();
        assert!(trophic_levels(&cyc, false).unwrap()[0].is_nan());
    }

    #[test]
    fn reaching_on_a_directed_star() {
        // 0 -> 1, 0 -> 2, 0 -> 3: the hub reaches everyone.
        let g = Graph::from_edges(4, &[0, 0, 0], &[1, 2, 3], None, true).unwrap();
        let l = local_reaching_centrality(&g, ReachingOptions::default()).unwrap();
        assert_eq!(l, vec![1.0, 0.0, 0.0, 0.0]);
        let grc = global_reaching_centrality(&g, ReachingOptions::default()).unwrap();
        assert!((grc - 1.0).abs() < 1e-12);
    }
}
