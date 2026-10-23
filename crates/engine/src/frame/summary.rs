//! Derived properties for printing a network: everything the summary card
//! of a binding shows, computed in one call from a built frame.
//!
//! Counts, directedness, modes and their sizes, the weight column with its
//! meaning and range, relations and their sizes, layers with their
//! intralayer sizes and the coupling and interlayer totals, the temporal
//! kind and time range, node spells and changes, geometry flags, egos
//! with their sizes, missing-data counts, loops, parallel edges,
//! isolates, weak components (`crate::metrics::topography::components`),
//! and density (`crate::metrics::topography::density`, loops excluded
//! from the denominator; NaN where undefined). Isolates, components, and
//! density are computed on the graph of observed ties.

use super::build::Built;
use super::roles::{CoordinateKind, EdgeRole, NodeRole, WeightType};
use super::validate::EdgeClass;
use super::Frame;
use crate::metrics::topography::{components, density, ComponentKind};

/// How a temporal network records time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeKind {
    /// No edge times.
    Static,
    /// A start time per edge (a contact sequence).
    Contacts,
    /// A start and an end per edge (a link stream).
    EdgeSpells,
    /// Discrete waves, wave `k` being `[k, k + 1)`.
    Waves,
}

impl TimeKind {
    /// The C interface code: 0 = static, 1 = contacts, 2 = edge spells,
    /// 3 = waves.
    pub fn code(self) -> i32 {
        match self {
            TimeKind::Static => 0,
            TimeKind::Contacts => 1,
            TimeKind::EdgeSpells => 2,
            TimeKind::Waves => 3,
        }
    }
}

/// The derived properties of a network.
#[derive(Debug, Clone, PartialEq)]
pub struct Summary {
    /// Number of nodes (table rows plus added nodes).
    pub nodes: usize,
    /// Number of edge rows (observed ties and unobserved pairs).
    pub edges: usize,
    /// Number of observed ties (the built graph's edges).
    pub graph_edges: usize,
    /// Whether edges are directed.
    pub directed: bool,
    /// Mode labels in order of first appearance (empty without a mode
    /// role).
    pub mode_labels: Vec<String>,
    /// Number of nodes in each mode.
    pub mode_counts: Vec<usize>,
    /// Whether edges were declared to join different modes only.
    pub between_modes_only: Option<bool>,
    /// The weight column, when the network is weighted.
    pub weight_column: Option<String>,
    /// The weights' declared meaning.
    pub weight_type: Option<WeightType>,
    /// Smallest weight (NaN when unweighted or without edges).
    pub weight_min: f64,
    /// Largest weight (NaN when unweighted or without edges).
    pub weight_max: f64,
    /// Relation labels in order of first appearance (empty without a
    /// relation role).
    pub relation_labels: Vec<String>,
    /// Number of edges in each relation.
    pub relation_counts: Vec<usize>,
    /// Whether a sign role is set.
    pub signed: bool,
    /// Number of negative edges in a signed network.
    pub negative_edges: usize,
    /// Layer labels (empty without layer roles).
    pub layer_labels: Vec<String>,
    /// Number of intralayer edges in each layer.
    pub layer_counts: Vec<usize>,
    /// Whether a layer table was given.
    pub layer_table: bool,
    /// Whether interlayer edges were declared to be couplings only.
    pub couplings_only: Option<bool>,
    /// Number of intralayer edges.
    pub intralayer: usize,
    /// Number of coupling edges (one actor, two layers).
    pub coupling: usize,
    /// Number of interlayer edges between different actors.
    pub interlayer: usize,
    /// How time is recorded.
    pub time: TimeKind,
    /// Earliest start time (NaN when static or without edges).
    pub time_start: f64,
    /// Latest end time (NaN when static or without edges).
    pub time_end: f64,
    /// Whether nodes have presence spells from the `node_start` and
    /// `node_end` roles.
    pub node_spells: bool,
    /// Whether nodes have presence spells from a node spells table.
    pub node_spell_table: bool,
    /// Rows of the node spells table.
    pub node_spell_rows: usize,
    /// Rows of the node changes table.
    pub node_change_rows: usize,
    /// Whether the node table carries geometry.
    pub node_geometry: bool,
    /// Whether the edge table carries geometry.
    pub edge_geometry: bool,
    /// The coordinate reference system, when geometry is present.
    pub crs: Option<String>,
    /// Whether representative coordinates are set, and what kind.
    pub coordinates: Option<CoordinateKind>,
    /// Ego labels (empty without an ego role).
    pub ego_labels: Vec<String>,
    /// Number of nodes in each ego's network.
    pub ego_node_counts: Vec<usize>,
    /// Number of edge rows in each ego's network.
    pub ego_edge_counts: Vec<usize>,
    /// Number of egos that are nodes.
    pub egos_as_nodes: usize,
    /// Whether egos were declared to be nodes.
    pub egos_are_nodes: Option<bool>,
    /// Number of self-loops.
    pub loops: usize,
    /// Number of parallel edges (edges repeating an earlier edge's key).
    pub parallel: usize,
    /// Number of nodes without edges.
    pub isolates: usize,
    /// Number of weak components.
    pub components: usize,
    /// Density with loops excluded from the denominator; NaN when
    /// undefined.
    pub density: f64,
    /// Nodes marked as unobserved.
    pub missing_nodes: usize,
    /// Edges marked as unobserved.
    pub missing_edges: usize,
}

/// Computes the derived properties of a built frame.
pub fn summarize(frame: &Frame, built: &Built) -> Summary {
    let r = &built.resolved;
    let g = &built.graph;
    let m = r.src.len();

    let counts = |codes: &[u32], k: usize| {
        let mut c = vec![0usize; k];
        for &x in codes {
            if x != u32::MAX {
                c[x as usize] += 1;
            }
        }
        c
    };
    let (mode_labels, mode_counts) = match &r.modes {
        Some(d) => {
            let mut c = counts(&d.codes, d.labels.len());
            // Added nodes have no mode; they are not counted anywhere.
            c.truncate(d.labels.len());
            (d.labels.clone(), c)
        }
        None => (Vec::new(), Vec::new()),
    };
    let (relation_labels, relation_counts) = match &r.relations {
        Some(d) => (d.labels.clone(), counts(&d.codes, d.labels.len())),
        None => (Vec::new(), Vec::new()),
    };
    let (weight_min, weight_max) = match &r.weights {
        Some(w) if !w.is_empty() => w
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
                (lo.min(x), hi.max(x))
            }),
        _ => (f64::NAN, f64::NAN),
    };
    let negative_edges = match frame.edge_role(EdgeRole::Sign) {
        Some((_, super::Column::Int(v, _))) => v.iter().filter(|&&x| x < 0).count(),
        Some((_, super::Column::Bool(v, _))) => v.iter().filter(|&&x| !x).count(),
        _ => 0,
    };
    let time = if frame.has_edge_role(EdgeRole::Wave) {
        TimeKind::Waves
    } else if frame.has_edge_role(EdgeRole::End) {
        TimeKind::EdgeSpells
    } else if frame.has_edge_role(EdgeRole::Start) {
        TimeKind::Contacts
    } else {
        TimeKind::Static
    };
    let (time_start, time_end) = match (&r.start, m) {
        (Some(s), m) if m > 0 => {
            let lo = s.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = r
                .end
                .as_ref()
                .unwrap_or(s)
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            (lo, hi)
        }
        _ => (f64::NAN, f64::NAN),
    };
    let geo = frame.geometry();
    let mut has_edge = vec![false; g.node_count()];
    for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
        has_edge[s as usize] = true;
        has_edge[t as usize] = true;
    }
    let flagged = |col: Option<(&str, &super::Column)>| match col {
        Some((_, super::Column::Bool(v, _))) => v.iter().filter(|&&x| x).count(),
        _ => 0,
    };
    let (mut layer_labels, mut layer_counts) = (Vec::new(), Vec::new());
    let (mut intralayer, mut coupling, mut interlayer) = (0, 0, 0);
    if let Some(l) = &r.layers {
        layer_labels = l.labels.clone();
        layer_counts = vec![0; l.labels.len()];
        for e in 0..m {
            match l.class(e, &r.src, &r.tgt) {
                EdgeClass::Intralayer => {
                    intralayer += 1;
                    layer_counts[l.src[e] as usize] += 1;
                }
                EdgeClass::Coupling => coupling += 1,
                EdgeClass::Interlayer => interlayer += 1,
            }
        }
    }
    let (ego_labels, ego_node_counts, ego_edge_counts, egos_as_nodes) = match &r.egos {
        Some(eg) => (
            eg.labels.clone(),
            counts(&eg.node, eg.labels.len()),
            counts(&eg.edge, eg.labels.len()),
            eg.node_of.iter().filter(|x| x.is_some()).count(),
        ),
        None => (Vec::new(), Vec::new(), Vec::new(), 0),
    };
    Summary {
        nodes: g.node_count(),
        edges: m,
        graph_edges: g.edge_count(),
        directed: frame.is_directed(),
        mode_labels,
        mode_counts,
        between_modes_only: frame.declarations().between_modes_only,
        weight_column: frame
            .edge_role(EdgeRole::Weight)
            .filter(|_| r.weights.is_some())
            .map(|(n, _)| n.to_string()),
        weight_type: frame.weight_type().filter(|_| r.weights.is_some()),
        weight_min,
        weight_max,
        relation_labels,
        relation_counts,
        signed: frame.has_edge_role(EdgeRole::Sign),
        negative_edges,
        layer_labels,
        layer_counts,
        layer_table: r.layers.as_ref().is_some_and(|l| l.from_table),
        couplings_only: frame.declarations().couplings_only,
        intralayer,
        coupling,
        interlayer,
        time,
        time_start,
        time_end,
        node_spells: r.node_spells.is_some(),
        node_spell_table: r.spell_table.is_some(),
        node_spell_rows: r.spell_table.as_ref().map_or(0, |t| t.rows.len()),
        node_change_rows: frame.rows(super::Which::NodeChanges),
        node_geometry: geo.nodes,
        edge_geometry: geo.edges,
        crs: geo.node_crs.clone().or_else(|| geo.edge_crs.clone()),
        coordinates: frame.has_node_role(NodeRole::X).then_some(geo.coordinates),
        ego_labels,
        ego_node_counts,
        ego_edge_counts,
        egos_as_nodes,
        egos_are_nodes: frame.declarations().egos_are_nodes,
        loops: r.loops,
        parallel: r.groups.as_ref().map_or(0, |gr| m - gr.sizes.len()),
        isolates: has_edge.iter().filter(|&&x| !x).count(),
        components: components(g, ComponentKind::Weak).count(),
        density: density(g, false),
        missing_nodes: flagged(frame.node_role(NodeRole::NodeMissing)),
        missing_edges: flagged(frame.edge_role(EdgeRole::EdgeMissing)),
    }
}

#[cfg(test)]
mod tests {
    use super::super::columns::{Column, Strings};
    use super::super::roles::Declarations;
    use super::*;

    fn strs(v: &[&str]) -> Column {
        Column::Str(Strings::from_strs(
            &v.iter().map(|s| Some(*s)).collect::<Vec<_>>(),
        ))
    }

    #[test]
    fn summary_of_a_rich_network() {
        let mut f = Frame::new(false);
        f.add_node_column("id", strs(&["a", "b", "c", "d", "e"]))
            .unwrap();
        f.add_node_column("mode", strs(&["p", "p", "e", "e", "p"]))
            .unwrap();
        f.add_node_column(
            "na",
            Column::Bool(vec![false, false, false, false, true], None),
        )
        .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.set_node_role(NodeRole::Mode, "mode").unwrap();
        f.set_node_role(NodeRole::NodeMissing, "na").unwrap();
        f.add_edge_column("from", strs(&["a", "a", "b", "b"]))
            .unwrap();
        f.add_edge_column("to", strs(&["c", "c", "d", "d"]))
            .unwrap();
        f.add_edge_column("w", Column::Float(vec![1.0, 2.0, 3.0, 4.0]))
            .unwrap();
        f.add_edge_column("rel", strs(&["x", "y", "x", "x"]))
            .unwrap();
        f.add_edge_column("s", Column::Float(vec![1.0, 2.0, 3.0, 4.0]))
            .unwrap();
        f.add_edge_column("e", Column::Float(vec![5.0, 2.0, 3.0, 9.0]))
            .unwrap();
        f.add_edge_column("sign", Column::Int(vec![1, -1, 1, -1], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_edge_role(EdgeRole::Relation, "rel").unwrap();
        f.set_edge_role(EdgeRole::Start, "s").unwrap();
        f.set_edge_role(EdgeRole::End, "e").unwrap();
        f.set_edge_role(EdgeRole::Sign, "sign").unwrap();
        f.set_weight_type(Some(WeightType::Strength));
        *f.declarations_mut() = Declarations {
            directed: false,
            between_modes_only: Some(true),
            ..Declarations::default()
        };
        f.geometry_mut().nodes = true;
        f.geometry_mut().node_crs = Some("EPSG:4326".into());
        let b = f.build(1).unwrap();
        let s = summarize(&f, &b);
        assert_eq!((s.nodes, s.edges, s.directed), (5, 4, false));
        assert_eq!(s.mode_labels, vec!["p", "e"]);
        assert_eq!(s.mode_counts, vec![3, 2]);
        assert_eq!(s.between_modes_only, Some(true));
        assert_eq!(s.weight_column.as_deref(), Some("w"));
        assert_eq!(s.weight_type, Some(WeightType::Strength));
        assert_eq!((s.weight_min, s.weight_max), (1.0, 4.0));
        assert_eq!(s.relation_labels, vec!["x", "y"]);
        assert_eq!(s.relation_counts, vec![3, 1]);
        assert!(s.signed);
        assert_eq!(s.negative_edges, 2);
        assert_eq!(s.time, TimeKind::EdgeSpells);
        assert_eq!((s.time_start, s.time_end), (1.0, 9.0));
        assert!(!s.node_spells);
        assert!(s.node_geometry && !s.edge_geometry);
        assert_eq!(s.crs.as_deref(), Some("EPSG:4326"));
        assert_eq!(s.coordinates, None);
        // b-d at [3, 3) and b-d at [4, 9) are different spells: not parallel.
        assert_eq!((s.loops, s.parallel, s.isolates), (0, 0, 1));
        assert_eq!(s.components, 3);
        assert!((s.density - 0.4).abs() < 1e-12);
        assert_eq!((s.missing_nodes, s.missing_edges), (1, 0));
    }

    #[test]
    fn summary_of_layers_egos_spells_and_missing() {
        use super::super::roles::{ChangeRole, SpellRole};
        let mut f = Frame::new(false);
        f.add_node_column("id", strs(&["a", "b", "c"])).unwrap();
        f.add_node_column("ego", strs(&["a", "a", "a"])).unwrap();
        f.add_node_column("na", Column::Bool(vec![false, true, false], None))
            .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.set_node_role(NodeRole::Ego, "ego").unwrap();
        f.set_node_role(NodeRole::NodeMissing, "na").unwrap();
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
        f.add_edge_column("s", Column::Float(vec![1.0; 6])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "ls").unwrap();
        f.set_edge_role(EdgeRole::TargetLayer, "lt").unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        f.set_edge_role(EdgeRole::Start, "s").unwrap();
        f.add_column(super::super::Which::Layers, "layer", strs(&["y", "x", "z"]))
            .unwrap();
        f.set_layer_role(super::super::roles::LayerRole::Layer, "layer")
            .unwrap();
        f.add_column(super::super::Which::NodeSpells, "id", strs(&["a", "a"]))
            .unwrap();
        f.add_column(
            super::super::Which::NodeSpells,
            "s",
            Column::Float(vec![0.0, 5.0]),
        )
        .unwrap();
        f.add_column(
            super::super::Which::NodeSpells,
            "e",
            Column::Float(vec![2.0, 9.0]),
        )
        .unwrap();
        f.set_spell_role(SpellRole::Id, "id").unwrap();
        f.set_spell_role(SpellRole::Start, "s").unwrap();
        f.set_spell_role(SpellRole::End, "e").unwrap();
        f.add_column(super::super::Which::NodeChanges, "id", strs(&["b"]))
            .unwrap();
        f.add_column(
            super::super::Which::NodeChanges,
            "t",
            Column::Float(vec![1.0]),
        )
        .unwrap();
        f.set_change_role(ChangeRole::Id, "id").unwrap();
        f.set_change_role(ChangeRole::Time, "t").unwrap();
        f.declarations_mut().couplings_only = Some(false);
        f.declarations_mut().egos_are_nodes = Some(true);
        let b = f.build(1).unwrap();
        let s = summarize(&f, &b);
        assert_eq!((s.nodes, s.edges, s.graph_edges), (3, 6, 5));
        assert_eq!(s.layer_labels, vec!["y", "x", "z"]);
        assert_eq!(s.layer_counts, vec![2, 2, 0]);
        assert!(s.layer_table);
        assert_eq!(s.couplings_only, Some(false));
        assert_eq!((s.intralayer, s.coupling, s.interlayer), (4, 1, 1));
        assert_eq!(s.ego_labels, vec!["a"]);
        assert_eq!((s.ego_node_counts, s.ego_edge_counts), (vec![3], vec![6]));
        assert_eq!((s.egos_as_nodes, s.egos_are_nodes), (1, Some(true)));
        assert_eq!(s.time, TimeKind::Contacts);
        assert!(!s.node_spells && s.node_spell_table);
        assert_eq!((s.node_spell_rows, s.node_change_rows), (2, 1));
        assert_eq!((s.missing_nodes, s.missing_edges), (1, 1));
        // The loop-free coupling is not a loop; isolates, components, and
        // density come from the observed graph (b-c in x is unobserved,
        // but b-c in y is a tie).
        assert_eq!((s.loops, s.parallel, s.isolates), (0, 0, 0));
        assert_eq!(s.components, 1);
        assert!((s.density - 5.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn summary_of_a_plain_network() {
        let mut f = Frame::new(true);
        f.add_edge_column("from", Column::Int(vec![1, 1, 2], None))
            .unwrap();
        f.add_edge_column("to", Column::Int(vec![2, 2, 2], None))
            .unwrap();
        f.add_edge_column("wave", Column::Int(vec![1, 1, 2], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::Wave, "wave").unwrap();
        let b = f.build(1).unwrap();
        let s = summarize(&f, &b);
        assert_eq!((s.nodes, s.edges, s.loops, s.parallel), (2, 3, 1, 1));
        assert!(s.mode_labels.is_empty() && s.relation_labels.is_empty());
        assert!(s.weight_column.is_none() && s.weight_min.is_nan());
        assert_eq!(s.time, TimeKind::Waves);
        assert_eq!((s.time_start, s.time_end), (1.0, 3.0));
        assert_eq!(s.components, 1);
        assert!((s.density - 1.5).abs() < 1e-12);
        assert_eq!(TimeKind::Contacts.code(), 1);
        let mut empty = Frame::new(true);
        empty
            .add_edge_column("from", Column::Int(vec![], None))
            .unwrap();
        empty
            .add_edge_column("to", Column::Int(vec![], None))
            .unwrap();
        empty.add_edge_column("t", Column::Float(vec![])).unwrap();
        empty.set_edge_role(EdgeRole::Source, "from").unwrap();
        empty.set_edge_role(EdgeRole::Target, "to").unwrap();
        empty.set_edge_role(EdgeRole::Start, "t").unwrap();
        let b = empty.build(1).unwrap();
        let s = summarize(&empty, &b);
        assert_eq!(s.time, TimeKind::Contacts);
        assert!(s.time_start.is_nan() && s.density.is_nan());
        assert_eq!(s.components, 0);
    }
}
