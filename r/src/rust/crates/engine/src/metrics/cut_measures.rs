//! Measures of a cut between node sets S and T: cut size, volume,
//! normalized cut, conductance, and edge, mixing, node, and boundary
//! expansion.
//!
//! Sources: normalized cut (Shi and Malik 2000, `shi2000normalized`);
//! conductance (Kannan, Vempala, and Vetta 2004, `kannan2004clusterings`);
//! edge, node, and boundary expansion (Hoory, Linial, and Wigderson 2006,
//! `hoory2006expander`). Every formula and convention, mixing expansion
//! and volume included, is that of NetworkX 3.7's `algorithms.cuts`
//! (`hagberg2008exploring`):
//!
//! - cut size: the edges joining S and T (T defaults to the complement of
//!   S); in a directed graph the arcs from S to T plus those from T to S.
//!   Written out for sets that may overlap, as NetworkX evaluates it: in an
//!   undirected graph an edge {a, b} counts once if (a ∈ S and b ∈ T) or
//!   (b ∈ S and a ∈ T) (with the complement: if exactly one end is in S);
//!   in a directed graph each arc a → b counts once for a ∈ S with b ∈ T or
//!   (b ∈ S and a ∈ T), and once more for a ∈ T with b ∈ S or (b ∈ T and a
//!   ∈ S). With `weighted` the weights are added instead.
//! - volume of S: the sum of the degrees of its nodes (out-degrees in a
//!   directed graph; strengths when weighted; an undirected loop counts
//!   twice).
//! - normalized cut: cut (1/vol(S) + 1/vol(T)); conductance: cut /
//!   min(vol(S), vol(T)); edge expansion: cut / min(|S|, |T|); mixing
//!   expansion: cut / (2 m), m the number of edges (parallel edges and
//!   loops included, never weighted).
//! - node expansion: |N(S)| / |S|, N(S) every neighbour (successor) of a
//!   node of S, S itself included; boundary expansion: |N(S) \ S| / |S|.
//!
//! Conventions: parallel edges each count (as in a NetworkX multigraph);
//! sets are sets (repeated nodes count once). Where NetworkX divides by
//! zero and raises `ZeroDivisionError` the result is the IEEE value (inf,
//! or NaN for 0/0). NetworkX's `cut_size` of a directed graph without T
//! raises a `TypeError`; here T is the complement there too.

use crate::error::{Error, Result};
use crate::graph::Graph;

/// Every measure of one cut. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutMeasures {
    /// Cut size (total weight when weighted).
    pub cut_size: f64,
    /// Volume of S.
    pub volume_s: f64,
    /// Volume of T.
    pub volume_t: f64,
    /// Normalized cut size.
    pub normalized_cut: f64,
    /// Conductance.
    pub conductance: f64,
    /// Edge expansion.
    pub edge_expansion: f64,
    /// Mixing expansion.
    pub mixing_expansion: f64,
    /// Node expansion (of S).
    pub node_expansion: f64,
    /// Boundary expansion (of S).
    pub boundary_expansion: f64,
}

fn set_of(n: usize, nodes: &[u32], what: &str) -> Result<(Vec<bool>, usize)> {
    let mut m = vec![false; n];
    let mut count = 0;
    for &v in nodes {
        if v as usize >= n {
            return Err(Error::InvalidArgument(format!(
                "{what} contains node {v}, outside 0..{n}"
            )));
        }
        if !m[v as usize] {
            m[v as usize] = true;
            count += 1;
        }
    }
    Ok((m, count))
}

/// The measures of the cut between `s` and `t` (the complement of `s`
/// when `None`). See the module documentation.
pub fn cut_measures(
    g: &Graph,
    s: &[u32],
    t: Option<&[u32]>,
    weighted: bool,
) -> Result<CutMeasures> {
    let n = g.node_count();
    if weighted && !g.is_weighted() {
        return Err(Error::MissingWeights);
    }
    let (in_s, ns) = set_of(n, s, "S")?;
    let (in_t, nt) = match t {
        Some(t) => set_of(n, t, "T")?,
        None => {
            let m: Vec<bool> = in_s.iter().map(|&x| !x).collect();
            (m, n - ns)
        }
    };
    let w = g.edge_weights();
    let mut cut = 0.0;
    for (e, (&a, &b)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        let (a, b) = (a as usize, b as usize);
        let x = if weighted {
            w.map_or(1.0, |w| w[e])
        } else {
            1.0
        };
        if g.is_directed() {
            if in_s[a] && (in_t[b] || (in_s[b] && in_t[a])) {
                cut += x;
            }
            if in_t[a] && (in_s[b] || (in_t[b] && in_s[a])) {
                cut += x;
            }
        } else if (in_s[a] && in_t[b]) || (in_s[b] && in_t[a]) {
            cut += x;
        }
    }
    let volume = |inside: &[bool]| -> f64 {
        let mut v = 0.0;
        for (e, (&a, &b)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
            let x = if weighted {
                w.map_or(1.0, |w| w[e])
            } else {
                1.0
            };
            if inside[a as usize] {
                v += x;
            }
            if !g.is_directed() && inside[b as usize] {
                v += x;
            }
        }
        v
    };
    let (vs, vt) = (volume(&in_s), volume(&in_t));
    let mut reach = vec![false; n];
    for v in (0..n).filter(|&v| in_s[v]) {
        for &u in g.out().neighbors(v) {
            reach[u as usize] = true;
        }
    }
    let nbhd = reach.iter().filter(|&&x| x).count() as f64;
    let boundary = (0..n).filter(|&v| reach[v] && !in_s[v]).count() as f64;
    let ns = ns as f64;
    Ok(CutMeasures {
        cut_size: cut,
        volume_s: vs,
        volume_t: vt,
        normalized_cut: cut * (1.0 / vs + 1.0 / vt),
        conductance: cut / vs.min(vt),
        edge_expansion: cut / ns.min(nt as f64),
        mixing_expansion: cut / (2.0 * g.edge_count() as f64),
        node_expansion: nbhd / ns,
        boundary_expansion: boundary / ns,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_path_cut_in_the_middle() {
        // 0-1-2-3, S = {0, 1}: one edge cut; volumes 3 and 3.
        let g = Graph::from_edges(4, &[0, 1, 2], &[1, 2, 3], None, false).unwrap();
        let c = cut_measures(&g, &[0, 1], None, false).unwrap();
        assert_eq!(c.cut_size, 1.0);
        assert_eq!((c.volume_s, c.volume_t), (3.0, 3.0));
        assert_eq!(c.normalized_cut, 2.0 / 3.0);
        assert_eq!(c.conductance, 1.0 / 3.0);
        assert_eq!(c.edge_expansion, 0.5);
        assert_eq!(c.mixing_expansion, 1.0 / 6.0);
        // N(S) = {0, 1, 2}; boundary {2}.
        assert_eq!(c.node_expansion, 1.5);
        assert_eq!(c.boundary_expansion, 0.5);
        // Directed cycle 0 -> 1 -> 2 -> 0 with S = {0}: one arc each way.
        let d = Graph::from_edges(3, &[0, 1, 2], &[1, 2, 0], None, true).unwrap();
        assert_eq!(cut_measures(&d, &[0], None, false).unwrap().cut_size, 2.0);
    }
}
