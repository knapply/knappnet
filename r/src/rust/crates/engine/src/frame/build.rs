//! Building the CSR [`Graph`] from a validated frame.
//!
//! Node `i` of the graph is row `i` of the node table, followed by the
//! nodes added from the edges in order of first appearance; edge `e` is
//! row `e` of the edge table. Weights are passed when the frame has a
//! valid weight role. Nothing is reordered, so results aligned with the
//! graph are aligned with the tables.

use super::columns::Column;
use super::roles::EdgeRole;
use super::validate::Resolved;
use super::Frame;
use crate::error::{Error, Result};
use crate::graph::Graph;

/// How graph node indices relate to the tables.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeMap {
    /// Rows of the node table (nodes `0..table_rows` are those rows).
    pub table_rows: usize,
    /// Nodes added from the edges, in node order after the table's rows:
    /// the edge row and endpoint (`true` for target) the identifier first
    /// appeared in.
    pub added: Vec<(usize, bool)>,
    /// The identifiers of the added nodes, as a column of the endpoint
    /// type, so bindings can append them to their node table.
    pub added_ids: Column,
}

/// A graph built from a frame, with the mapping and the resolved columns.
#[derive(Debug, Clone)]
pub struct Built {
    /// The graph.
    pub graph: Graph,
    /// Node index mapping.
    pub map: NodeMap,
    /// What validation resolved (indices, modes, relations, times).
    pub resolved: Resolved,
    /// The notes validation produced (there are no errors).
    pub report: super::Report,
}

/// Validates `frame` and builds the graph, or returns the first error's
/// message (every error is in [`super::Validation`]; callers that want
/// them all validate first).
pub fn build(frame: &Frame, threads: usize) -> Result<Built> {
    let validation = frame.validate(threads);
    build_validated(frame, validation)
}

/// Builds the graph from a validation already run on `frame`.
pub fn build_validated(frame: &Frame, validation: super::Validation) -> Result<Built> {
    if !validation.report.is_valid() {
        return Err(Error::InvalidArgument(validation.report.error_text()));
    }
    let resolved = validation
        .resolved
        .ok_or_else(|| Error::InvalidArgument("the frame's endpoints did not resolve".into()))?;
    let src: Vec<i64> = resolved.src.iter().map(|&i| i as i64).collect();
    let tgt: Vec<i64> = resolved.tgt.iter().map(|&i| i as i64).collect();
    let graph = Graph::from_edges(
        resolved.n,
        &src,
        &tgt,
        resolved.weights.as_deref(),
        frame.is_directed(),
    )?;
    let map = NodeMap {
        table_rows: frame.node_rows(),
        added: resolved.added.clone(),
        added_ids: added_ids(frame, &resolved.added),
    };
    Ok(Built {
        graph,
        map,
        resolved,
        report: validation.report,
    })
}

/// The identifiers of the added nodes, taken from the endpoint columns.
fn added_ids(frame: &Frame, added: &[(usize, bool)]) -> Column {
    let (_, src) = frame.edge_role(EdgeRole::Source).expect("validated");
    let (_, tgt) = frame.edge_role(EdgeRole::Target).expect("validated");
    let from_src: Vec<usize> = added.iter().filter(|a| !a.1).map(|a| a.0).collect();
    let from_tgt: Vec<usize> = added.iter().filter(|a| a.1).map(|a| a.0).collect();
    // Keep node order: take rows from each column, then interleave.
    let s = src.take(&from_src);
    let t = tgt.take(&from_tgt);
    let (mut i, mut j) = (0, 0);
    let mut order = Vec::with_capacity(added.len());
    for a in added {
        if a.1 {
            order.push((true, j));
            j += 1;
        } else {
            order.push((false, i));
            i += 1;
        }
    }
    match (s, t) {
        (Column::Int(sv, sm), Column::Int(tv, tm)) => {
            let pick = |(from_t, k): &(bool, usize)| if *from_t { tv[*k] } else { sv[*k] };
            let values: Vec<i64> = order.iter().map(pick).collect();
            let any_missing = sm.is_some() || tm.is_some();
            let missing = any_missing.then(|| {
                order
                    .iter()
                    .map(|(from_t, k)| {
                        let m = if *from_t { &tm } else { &sm };
                        m.as_ref().is_some_and(|m| m[*k])
                    })
                    .collect()
            });
            Column::Int(values, missing)
        }
        (Column::Float(sv), Column::Float(tv)) => Column::Float(
            order
                .iter()
                .map(|(from_t, k)| if *from_t { tv[*k] } else { sv[*k] })
                .collect(),
        ),
        (Column::Str(s), Column::Str(t)) => {
            let values: Vec<Option<&str>> = order
                .iter()
                .map(|(from_t, k)| {
                    if *from_t {
                        t.get_str(*k)
                    } else {
                        s.get_str(*k)
                    }
                })
                .collect();
            Column::Str(super::columns::Strings::from_strs(&values))
        }
        // Integer and double endpoints were matched as numbers; the
        // identifiers are reported as doubles.
        (s, t) => Column::Float(
            order
                .iter()
                .map(|(from_t, k)| {
                    let c = if *from_t { &t } else { &s };
                    match c {
                        Column::Int(v, _) => v[*k] as f64,
                        Column::Float(v) => v[*k],
                        _ => f64::NAN,
                    }
                })
                .collect(),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::super::columns::Strings;
    use super::super::roles::{NodeRole, WeightType};
    use super::*;

    fn strs(v: &[&str]) -> Column {
        Column::Str(Strings::from_strs(
            &v.iter().map(|s| Some(*s)).collect::<Vec<_>>(),
        ))
    }

    #[test]
    fn builds_the_same_graph_as_index_construction() {
        let mut f = Frame::new(true);
        f.add_node_column("id", strs(&["a", "b"])).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column("from", strs(&["a", "b", "c"])).unwrap();
        f.add_edge_column("to", strs(&["b", "c", "a"])).unwrap();
        f.add_edge_column("w", Column::Float(vec![1.0, 2.0, 3.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Distance));
        assert!(build(&f, 1).is_err());
        f.declarations_mut().add_nodes = true;
        let b = build(&f, 1).unwrap();
        let expected =
            Graph::from_edges(3, &[0, 1, 2], &[1, 2, 0], Some(&[1.0, 2.0, 3.0]), true).unwrap();
        assert_eq!(b.graph.edge_sources(), expected.edge_sources());
        assert_eq!(b.graph.edge_targets(), expected.edge_targets());
        assert_eq!(b.graph.edge_weights(), expected.edge_weights());
        assert_eq!(b.graph.out().neighbors(2), &[0]);
        assert_eq!(b.map.table_rows, 2);
        assert_eq!(b.map.added, vec![(1, true)]);
        assert_eq!(b.map.added_ids.render(0), "c");
        assert_eq!(b.report.notes().count(), 1);
    }

    #[test]
    fn added_ids_keep_node_order_and_type() {
        let mut f = Frame::new(false);
        f.add_edge_column("from", Column::Int(vec![5, 7, 9], None))
            .unwrap();
        f.add_edge_column("to", Column::Int(vec![6, 5, 7], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let b = build(&f, 1).unwrap();
        assert_eq!(b.map.added_ids, Column::Int(vec![5, 6, 7, 9], None));
        assert_eq!(b.graph.node_count(), 4);
        // Mixed integer and double endpoints come back as doubles.
        let mut f = Frame::new(false);
        f.add_edge_column("from", Column::Int(vec![1, 2], None))
            .unwrap();
        f.add_edge_column("to", Column::Float(vec![2.0, 3.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let b = build(&f, 1).unwrap();
        assert_eq!(b.map.added_ids, Column::Float(vec![1.0, 2.0, 3.0]));
        let mut f = Frame::new(false);
        f.add_edge_column("from", Column::Float(vec![0.5])).unwrap();
        f.add_edge_column("to", Column::Float(vec![1.5])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let b = build(&f, 1).unwrap();
        assert_eq!(b.map.added_ids, Column::Float(vec![0.5, 1.5]));
    }
}
