//! Parallel edges under the role-relative rule of the module
//! documentation of [`super`]: edge multiplicity, the simplicity test, and
//! `simplify`, which groups parallel edges and combines their weights,
//! returning an index mapping so bindings can aggregate other attributes.
//!
//! Sources: igraph's `simplify()` (`csardi2006igraph`) defines multiple
//! edges by endpoints alone and combines attributes with `edge.attr.comb`;
//! Pajek's "Remove Multiple Lines" offers sum, count, min, max, and single
//! (`denooy2018exploratory`). The choices here are the same; the key adds
//! the relation and the time spell so that a multiplex or temporal network
//! with one edge per relation and spell is simple.

use super::validate::{edge_groups, Groups, KeyColumns, Resolved};
use crate::error::{Error, Result};

/// How the weights of a group of parallel edges are combined.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Combine {
    /// Their sum (the usual choice for strengths).
    #[default]
    Sum,
    /// Their mean.
    Mean,
    /// The smallest (the usual choice for distances).
    Min,
    /// The largest.
    Max,
    /// The first edge's weight, in row order.
    First,
    /// The number of edges in the group, ignoring weights.
    Count,
}

impl Combine {
    /// Decodes the C interface codes: 0 = sum, 1 = mean, 2 = min, 3 = max,
    /// 4 = first, 5 = count.
    pub fn from_code(code: i32) -> Result<Combine> {
        match code {
            0 => Ok(Combine::Sum),
            1 => Ok(Combine::Mean),
            2 => Ok(Combine::Min),
            3 => Ok(Combine::Max),
            4 => Ok(Combine::First),
            5 => Ok(Combine::Count),
            _ => Err(Error::InvalidArgument(format!(
                "unknown combine code {code}"
            ))),
        }
    }
}

/// The result of [`simplify`].
#[derive(Debug, Clone, PartialEq)]
pub struct Simplified {
    /// For every original edge, the index of the edge it became in the
    /// simplified network, or `None` for a removed loop. Kept edges are
    /// numbered in order of first appearance.
    pub group: Vec<Option<u32>>,
    /// The original row of each kept edge (its group's first row), so
    /// bindings can copy that row's attributes.
    pub keep: Vec<usize>,
    /// Combined weight of each kept edge, when the network is weighted or
    /// `Combine::Count` was asked for.
    pub weight: Option<Vec<f64>>,
}

/// Number of edges sharing each edge's key (1 for an edge with no
/// parallel). Needs the groups from validation.
pub fn multiplicity(groups: &Groups) -> Vec<u32> {
    groups
        .group
        .iter()
        .map(|&g| groups.sizes[g as usize])
        .collect()
}

/// Whether the network has neither loops nor parallel edges.
pub fn is_simple(resolved: &Resolved) -> bool {
    resolved.loops == 0
        && resolved
            .groups
            .as_ref()
            .map_or(true, |g| g.sizes.len() == g.group.len())
}

/// Groups parallel edges (same endpoints with their layers, relation, and
/// spell) and, with `drop_loops`, removes self-loops (a coupling edge is
/// not one). `combine` says how the weights of a group are combined;
/// without weights only `Combine::Count` produces a weight column. Rows
/// flagged as unobserved pairs are grouped like any other; the binding
/// keeps or drops them by their flag.
pub fn simplify(resolved: &Resolved, combine: Combine, drop_loops: bool) -> Simplified {
    let groups = match &resolved.groups {
        Some(g) => g.clone(),
        None => edge_groups(
            resolved.n,
            &resolved.src,
            &resolved.tgt,
            resolved.directed,
            KeyColumns::of(resolved),
        ),
    };
    let m = groups.group.len();
    let k = groups.sizes.len();
    // Renumber the kept groups in order of first appearance.
    let mut new_index = vec![None; k];
    let mut keep = Vec::with_capacity(k);
    for (g, &first) in groups.first.iter().enumerate() {
        if drop_loops && resolved.is_loop(first) {
            continue;
        }
        new_index[g] = Some(keep.len() as u32);
        keep.push(first);
    }
    let group: Vec<Option<u32>> = groups
        .group
        .iter()
        .map(|&g| new_index[g as usize])
        .collect();
    let weight = match (&resolved.weights, combine) {
        (_, Combine::Count) => Some(
            keep.iter()
                .map(|&first| groups.sizes[groups.group[first] as usize] as f64)
                .collect(),
        ),
        (None, _) => None,
        (Some(w), combine) => {
            let mut acc = vec![0.0; keep.len()];
            let mut seen = vec![false; keep.len()];
            for e in 0..m {
                let Some(j) = group[e] else { continue };
                let j = j as usize;
                let x = w[e];
                if !seen[j] {
                    seen[j] = true;
                    acc[j] = x;
                    continue;
                }
                match combine {
                    Combine::Sum | Combine::Mean => acc[j] += x,
                    Combine::Min => acc[j] = acc[j].min(x),
                    Combine::Max => acc[j] = acc[j].max(x),
                    Combine::First => {}
                    Combine::Count => unreachable!(),
                }
            }
            if combine == Combine::Mean {
                for (j, &first) in keep.iter().enumerate() {
                    acc[j] /= groups.sizes[groups.group[first] as usize] as f64;
                }
            }
            Some(acc)
        }
    };
    Simplified {
        group,
        keep,
        weight,
    }
}

#[cfg(test)]
mod tests {
    use super::super::columns::Column;
    use super::super::roles::{EdgeRole, WeightType};
    use super::super::Frame;
    use super::*;

    /// Undirected: 1-2, 2-1, 1-1, 2-3, 2-3, 1-1 with weights 1..6.
    fn frame() -> Frame {
        let mut f = Frame::new(false);
        f.add_edge_column("from", Column::Int(vec![1, 2, 1, 2, 2, 1], None))
            .unwrap();
        f.add_edge_column("to", Column::Int(vec![2, 1, 1, 3, 3, 1], None))
            .unwrap();
        f.add_edge_column("w", Column::Float(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Strength));
        f
    }

    #[test]
    fn multiplicity_and_simplicity() {
        let f = frame();
        let v = f.validate(1);
        let r = v.resolved.unwrap();
        assert!(!is_simple(&r));
        assert_eq!(
            multiplicity(r.groups.as_ref().unwrap()),
            vec![2, 2, 2, 2, 2, 2]
        );
        let mut f = Frame::new(true);
        f.add_edge_column("from", Column::Int(vec![1, 2], None))
            .unwrap();
        f.add_edge_column("to", Column::Int(vec![2, 1], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let r = f.validate(1).resolved.unwrap();
        assert!(is_simple(&r));
        assert_eq!(multiplicity(r.groups.as_ref().unwrap()), vec![1, 1]);
    }

    #[test]
    fn simplify_combines_weights() {
        let f = frame();
        let r = f.validate(1).resolved.unwrap();
        let s = simplify(&r, Combine::Sum, false);
        assert_eq!(s.keep, vec![0, 2, 3]);
        assert_eq!(
            s.group,
            vec![Some(0), Some(0), Some(1), Some(2), Some(2), Some(1)]
        );
        assert_eq!(s.weight.unwrap(), vec![3.0, 9.0, 9.0]);
        let s = simplify(&r, Combine::Mean, true);
        assert_eq!(s.keep, vec![0, 3]);
        assert_eq!(
            s.group,
            vec![Some(0), Some(0), None, Some(1), Some(1), None]
        );
        assert_eq!(s.weight.unwrap(), vec![1.5, 4.5]);
        assert_eq!(
            simplify(&r, Combine::Min, false).weight.unwrap(),
            vec![1.0, 3.0, 4.0]
        );
        assert_eq!(
            simplify(&r, Combine::Max, false).weight.unwrap(),
            vec![2.0, 6.0, 5.0]
        );
        assert_eq!(
            simplify(&r, Combine::First, false).weight.unwrap(),
            vec![1.0, 3.0, 4.0]
        );
        assert_eq!(
            simplify(&r, Combine::Count, true).weight.unwrap(),
            vec![2.0, 2.0]
        );
        // Without weights only a count gives a weight column; groups are
        // recomputed when validation did not keep them.
        let mut unweighted = r.clone();
        unweighted.weights = None;
        unweighted.groups = None;
        assert!(simplify(&unweighted, Combine::Sum, false).weight.is_none());
        let s = simplify(&unweighted, Combine::Count, false);
        assert_eq!(s.weight.unwrap(), vec![2.0, 2.0, 2.0]);
        // Directed: 1 -> 2 and 2 -> 1 stay apart.
        let mut f = frame();
        f.declarations_mut().directed = true;
        let r = f.validate(1).resolved.unwrap();
        let s = simplify(&r, Combine::Sum, false);
        assert_eq!(s.keep, vec![0, 1, 2, 3]);
    }

    #[test]
    fn combine_codes() {
        assert_eq!(Combine::from_code(5).unwrap(), Combine::Count);
        assert!(Combine::from_code(6).is_err());
    }
}
