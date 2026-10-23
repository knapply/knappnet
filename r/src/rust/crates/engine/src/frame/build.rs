//! Building the CSR [`Graph`] from a validated frame, and subgraphs of
//! one layer or one ego from a built frame.
//!
//! Node `i` of the graph is row `i` of the node table, followed by the
//! nodes added from the edges in order of first appearance. The graph's
//! edges are the observed rows of the edge table (rows not flagged as
//! unobserved pairs) in row order, so without an `edge_missing` role
//! edge `e` is row `e`; [`EdgeMap`] gives the row of every graph edge
//! and [`missing_pairs`] the pairs left out. Weights are passed when the
//! frame has a valid weight role. Nothing is reordered, so results
//! aligned with the graph are aligned with the tables through the maps.

use super::columns::Column;
use super::roles::EdgeRole;
use super::validate::{EdgeClass, Resolved};
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

/// How graph edge indices relate to the edge table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeMap {
    /// The edge table row of every graph edge, in graph edge order;
    /// `None` when the two coincide (no row is an unobserved pair).
    pub rows: Option<Vec<usize>>,
    /// Number of graph edges.
    pub count: usize,
}

impl EdgeMap {
    /// The edge table row of graph edge `e`.
    #[inline]
    pub fn row(&self, e: usize) -> usize {
        self.rows.as_ref().map_or(e, |r| r[e])
    }
}

/// A graph built from a frame, with the mappings and the resolved columns.
#[derive(Debug, Clone)]
pub struct Built {
    /// The graph of the observed ties.
    pub graph: Graph,
    /// Node index mapping.
    pub map: NodeMap,
    /// Edge index mapping.
    pub edges: EdgeMap,
    /// What validation resolved (indices, modes, relations, times,
    /// layers, egos, missing flags).
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
    let rows: Option<Vec<usize>> = resolved
        .missing
        .as_ref()
        .filter(|m| m.iter().any(|&x| x))
        .map(|m| (0..m.len()).filter(|&e| !m[e]).collect());
    let (graph, count) = graph_of(&resolved, rows.as_deref(), frame.is_directed())?;
    let map = NodeMap {
        table_rows: frame.node_rows(),
        added: resolved.added.clone(),
        added_ids: added_ids(frame, &resolved.added),
    };
    Ok(Built {
        graph,
        map,
        edges: EdgeMap { rows, count },
        resolved,
        report: validation.report,
    })
}

/// The graph of the edge rows `rows` (every row when `None`) over all
/// `resolved.n` nodes, with the number of its edges.
fn graph_of(resolved: &Resolved, rows: Option<&[usize]>, directed: bool) -> Result<(Graph, usize)> {
    let pick = |v: &[u32]| -> Vec<i64> {
        match rows {
            Some(rows) => rows.iter().map(|&e| v[e] as i64).collect(),
            None => v.iter().map(|&i| i as i64).collect(),
        }
    };
    let src = pick(&resolved.src);
    let tgt = pick(&resolved.tgt);
    let weights: Option<Vec<f64>> = resolved.weights.as_ref().map(|w| match rows {
        Some(rows) => rows.iter().map(|&e| w[e]).collect(),
        None => w.clone(),
    });
    let graph = Graph::from_edges(resolved.n, &src, &tgt, weights.as_deref(), directed)?;
    Ok((graph, src.len()))
}

/// The unobserved pairs of a built frame: the node indices of each row
/// flagged as missing, with the row, in row order.
pub fn missing_pairs(built: &Built) -> Vec<(u32, u32, usize)> {
    let r = &built.resolved;
    match &r.missing {
        Some(m) => (0..m.len())
            .filter(|&e| m[e])
            .map(|e| (r.src[e], r.tgt[e], e))
            .collect(),
        None => Vec::new(),
    }
}

/// A subgraph of a built frame: the graph of some nodes and some of the
/// observed edge rows, with the mappings back to the tables.
#[derive(Debug, Clone)]
pub struct Subgraph {
    /// The graph; its node `i` is node `nodes[i]` of the built graph.
    pub graph: Graph,
    /// For every subgraph node, the built graph's node index (a node
    /// table row or an added node), ascending.
    pub nodes: Vec<u32>,
    /// For every subgraph edge, the edge table row, ascending.
    pub edges: Vec<usize>,
}

/// Builds the subgraph of the nodes for which `keep_node` holds and the
/// observed edge rows for which `keep_edge` holds (an edge whose endpoint
/// is dropped is dropped too). Nodes and edges keep their relative
/// order.
pub fn subgraph(
    built: &Built,
    keep_node: impl Fn(usize) -> bool,
    keep_edge: impl Fn(usize) -> bool,
) -> Result<Subgraph> {
    let r = &built.resolved;
    let mut index = vec![u32::MAX; r.n];
    let mut nodes = Vec::new();
    for (i, slot) in index.iter_mut().enumerate() {
        if keep_node(i) {
            *slot = nodes.len() as u32;
            nodes.push(i as u32);
        }
    }
    let edges: Vec<usize> = (0..r.src.len())
        .filter(|&e| {
            r.observed(e)
                && keep_edge(e)
                && index[r.src[e] as usize] != u32::MAX
                && index[r.tgt[e] as usize] != u32::MAX
        })
        .collect();
    let src: Vec<i64> = edges
        .iter()
        .map(|&e| index[r.src[e] as usize] as i64)
        .collect();
    let tgt: Vec<i64> = edges
        .iter()
        .map(|&e| index[r.tgt[e] as usize] as i64)
        .collect();
    let weights: Option<Vec<f64>> = r
        .weights
        .as_ref()
        .map(|w| edges.iter().map(|&e| w[e]).collect());
    let graph = Graph::from_edges(nodes.len(), &src, &tgt, weights.as_deref(), r.directed)?;
    Ok(Subgraph {
        graph,
        nodes,
        edges,
    })
}

/// The graph of one layer: every actor, and the observed intralayer
/// edges whose layer is `layer` (an index into the resolved layer
/// labels). Coupling and interlayer edges are left out.
pub fn layer_subgraph(built: &Built, layer: usize) -> Result<Subgraph> {
    let r = &built.resolved;
    let layers = r
        .layers
        .as_ref()
        .ok_or_else(|| Error::InvalidArgument("the frame has no layer roles".into()))?;
    if layer >= layers.labels.len() {
        return Err(Error::InvalidArgument(format!(
            "layer index {layer} is out of range; the network has {} layer(s)",
            layers.labels.len()
        )));
    }
    let layer = layer as u32;
    subgraph(
        built,
        |_| true,
        |e| layers.src[e] == layer && layers.class(e, &r.src, &r.tgt) == EdgeClass::Intralayer,
    )
}

/// The graph of one ego's network: the nodes and the observed edges
/// whose ego is `ego` (an index into the resolved ego labels).
pub fn ego_subgraph(built: &Built, ego: usize) -> Result<Subgraph> {
    let r = &built.resolved;
    let egos = r
        .egos
        .as_ref()
        .ok_or_else(|| Error::InvalidArgument("the frame has no ego role".into()))?;
    if ego >= egos.labels.len() {
        return Err(Error::InvalidArgument(format!(
            "ego index {ego} is out of range; the network has {} ego(s)",
            egos.labels.len()
        )));
    }
    let ego = ego as u32;
    subgraph(built, |i| egos.node[i] == ego, |e| egos.edge[e] == ego)
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
    fn unobserved_pairs_are_left_out_of_the_graph() {
        let mut f = Frame::new(false);
        f.add_edge_column("from", Column::Int(vec![1, 1, 2, 3], None))
            .unwrap();
        f.add_edge_column("to", Column::Int(vec![2, 3, 3, 4], None))
            .unwrap();
        f.add_edge_column("w", Column::Float(vec![1.0, 2.0, 3.0, 4.0]))
            .unwrap();
        f.add_edge_column("na", Column::Bool(vec![false, true, false, true], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        f.set_weight_type(Some(WeightType::Strength));
        let b = build(&f, 1).unwrap();
        assert_eq!(b.graph.node_count(), 4);
        assert_eq!(b.graph.edge_count(), 2);
        assert_eq!(b.graph.edge_sources(), &[0, 1]);
        assert_eq!(b.graph.edge_targets(), &[1, 2]);
        assert_eq!(b.graph.edge_weights(), Some(&[1.0, 3.0][..]));
        assert_eq!(b.edges.rows, Some(vec![0, 2]));
        assert_eq!((b.edges.count, b.edges.row(1)), (2, 2));
        assert_eq!(missing_pairs(&b), vec![(0, 2, 1), (2, 3, 3)]);
        // Without flagged rows the mapping is the identity.
        f.clear_edge_role(EdgeRole::EdgeMissing);
        let b = build(&f, 1).unwrap();
        assert_eq!(b.edges.rows, None);
        assert_eq!((b.edges.count, b.edges.row(3)), (4, 3));
        assert!(missing_pairs(&b).is_empty());
    }

    #[test]
    fn layer_and_ego_subgraphs() {
        // Actors a, b, c in layers x and y; a-b in both, b-c in y, a
        // coupled across layers, a-c interlayer, and an unobserved b-c
        // in x.
        let mut f = Frame::new(false);
        f.add_node_column("id", strs(&["a", "b", "c"])).unwrap();
        f.add_node_column("ego", strs(&["a", "a", "a"])).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.set_node_role(NodeRole::Ego, "ego").unwrap();
        f.add_edge_column("from", strs(&["a", "a", "b", "a", "a", "b"]))
            .unwrap();
        f.add_edge_column("to", strs(&["b", "b", "c", "a", "c", "c"]))
            .unwrap();
        f.add_edge_column("ls", strs(&["x", "y", "y", "x", "x", "x"]))
            .unwrap();
        f.add_edge_column("lt", strs(&["x", "y", "y", "y", "y", "x"]))
            .unwrap();
        f.add_edge_column(
            "na",
            Column::Bool(vec![false, false, false, false, false, true], None),
        )
        .unwrap();
        f.add_edge_column("w", Column::Float(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "ls").unwrap();
        f.set_edge_role(EdgeRole::TargetLayer, "lt").unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Strength));
        let b = build(&f, 1).unwrap();
        assert_eq!(b.graph.edge_count(), 5);
        let x = layer_subgraph(&b, 0).unwrap();
        assert_eq!(x.nodes, vec![0, 1, 2]);
        assert_eq!(x.edges, vec![0]);
        assert_eq!(x.graph.edge_weights(), Some(&[1.0][..]));
        let y = layer_subgraph(&b, 1).unwrap();
        assert_eq!(y.edges, vec![1, 2]);
        assert_eq!(y.graph.edge_sources(), &[0, 1]);
        assert!(layer_subgraph(&b, 2).is_err());
        let e = ego_subgraph(&b, 0).unwrap();
        assert_eq!(e.nodes, vec![0, 1, 2]);
        assert_eq!(e.edges, vec![0, 1, 2, 3, 4]);
        assert!(ego_subgraph(&b, 1).is_err());
        // A subgraph of two nodes renumbers them and drops the edges
        // that leave.
        let s = subgraph(&b, |i| i != 1, |_| true).unwrap();
        assert_eq!(s.nodes, vec![0, 2]);
        assert_eq!(s.edges, vec![3, 4]);
        assert_eq!(s.graph.edge_sources(), &[0, 0]);
        assert_eq!(s.graph.edge_targets(), &[0, 1]);
        // Ego subgraphs on a frame of two egos.
        let mut g = Frame::new(true);
        g.add_node_column("id", strs(&["p", "p1", "q", "q1"]))
            .unwrap();
        g.add_node_column("ego", strs(&["p", "p", "q", "q"]))
            .unwrap();
        g.set_node_role(NodeRole::Id, "id").unwrap();
        g.set_node_role(NodeRole::Ego, "ego").unwrap();
        g.add_edge_column("from", strs(&["p", "q", "q1"])).unwrap();
        g.add_edge_column("to", strs(&["p1", "q1", "q"])).unwrap();
        g.set_edge_role(EdgeRole::Source, "from").unwrap();
        g.set_edge_role(EdgeRole::Target, "to").unwrap();
        let b = build(&g, 1).unwrap();
        let q = ego_subgraph(&b, 1).unwrap();
        assert_eq!(q.nodes, vec![2, 3]);
        assert_eq!(q.edges, vec![1, 2]);
        assert_eq!(q.graph.edge_sources(), &[0, 1]);
        assert!(ego_subgraph(&b, 2).is_err());
        let mut h = Frame::new(true);
        h.add_edge_column("from", strs(&["p"])).unwrap();
        h.add_edge_column("to", strs(&["q"])).unwrap();
        h.set_edge_role(EdgeRole::Source, "from").unwrap();
        h.set_edge_role(EdgeRole::Target, "to").unwrap();
        let b = build(&h, 1).unwrap();
        assert!(ego_subgraph(&b, 0).is_err() && layer_subgraph(&b, 0).is_err());
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
