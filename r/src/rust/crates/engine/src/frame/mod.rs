//! The data model: two typed tables, the roles that give their columns
//! structural meaning, three optional long tables, declared properties,
//! validation with structured diagnostics, and construction of a
//! [`Graph`](crate::Graph) from them. Every language binding builds its
//! network object through this module, so the checks and their messages
//! are the same in R and Python. The specification this module implements
//! is `docs/DATA_MODEL.md`.
//!
//! # The model
//!
//! A network is a pair of tables (`wasserman1994social`, chapters 2 and 3,
//! distinguish the structural variables, the ties, from the composition
//! variables, the attributes; the two-table form is the edge list with a
//! separate vertex table that `butts2008network` and igraph's
//! `graph_from_data_frame()` take):
//!
//! - **nodes**: one row per node; an identifier column (role `id`) and
//!   any attributes.
//! - **edges**: one row per edge; two endpoint columns (roles `source` and
//!   `target`) holding node identifiers, and any attributes.
//!
//! and three optional long tables attached to the network:
//!
//! - **layers**: one row per layer of a multilayer network (role `layer`
//!   is the identifier; attributes stay with the binding);
//! - **node_spells**: one row per presence spell of a node (roles `id`,
//!   `start`, `end`; several rows per node);
//! - **node_changes**: one row per change of a node attribute (roles `id`
//!   and `time`; the attribute columns stay with the binding).
//!
//! The engine receives only the columns that carry a role, as typed
//! columns ([`Column`]: integer, double, text, or logical, each with
//! missing values); attributes stay with the binding. A frame without a
//! node table takes its nodes from the edges.
//!
//! # Roles and their definitions
//!
//! Node roles ([`NodeRole`]) and edge roles ([`EdgeRole`]) are:
//!
//! - `id`, `source`, `target`: identifiers. Integer, double with whole
//!   values, or text; the endpoint columns must share one type and the
//!   identifier column must match it. Missing identifiers and empty
//!   strings are errors; so are endpoints absent from the node table
//!   unless the frame declares `add_nodes` (or has no node table), in
//!   which case they are appended in order of first appearance and
//!   reported.
//! - `weight` with a [`WeightType`]: a network is weighted when one numeric
//!   edge column is declared the weight, and its meaning must be declared
//!   once, `strength` (larger brings nodes closer: frequency, intensity,
//!   capacity; `barrat2004architecture`) or `distance` (larger keeps them
//!   apart: cost, length, latency; `newman2001scientific`,
//!   `opsahl2010node`). There is no default. Weights must be finite;
//!   distances must be positive.
//! - `mode`: a declared, disjoint set of nodes of one kind. Any discrete
//!   node column (text, integer, or logical; igraph's logical `type` is
//!   accepted as is). Two modes with `between_modes_only` is a two-mode
//!   (affiliation) network in the sense of `wasserman1994social` chapter 8
//!   and `borgatti1997network`: every edge joins the two modes; k modes
//!   with the constraint is a k-partite network; modes without it allow
//!   ties within a mode (a multimode or multilevel network). Being
//!   bipartite as a graph (2-colourable) is a separate, structural
//!   property that needs no declaration and is not checked here. Coupling
//!   edges of a multilayer network (below) are exempt from the rule.
//! - `relation`: the relation an edge belongs to, in a multiplex network
//!   (several relations on one node set, each edge in exactly one:
//!   `wasserman1994social` chapter 3; the edge-coloured multigraph of
//!   `kivela2014multilayer`). Any discrete edge column without missing
//!   values.
//! - `sign`: `-1`/`1` as an integer column, or a logical column where
//!   `true` is positive (`cartwright1956structural`). The sign is
//!   separate from the weight, which then holds the magnitude and must
//!   not be negative, so a negative tie is never read as a negative
//!   length.
//! - `source_layer`, `target_layer`: the layers of the two endpoints in a
//!   multilayer network (`kivela2014multilayer`,
//!   `dedomenico2013mathematical`; the actor-by-layer vocabulary of
//!   multinet, `magnani2021analysis`, and pymnet, `nurmi2024pymnet`).
//!   Actors stay one row each in the node table; the two discrete edge
//!   columns are set together and have no missing values. An edge whose
//!   layers are equal is intralayer; one whose layers differ is a
//!   coupling when it joins an actor to itself and an interlayer edge
//!   otherwise. The optional layer table fixes the layer identifiers and
//!   their order; every edge layer must then be one of them. The
//!   declaration `couplings_only` makes every interlayer edge a coupling
//!   (the multiplex special case).
//! - `start`, `end`: times of a temporal network, as doubles in the
//!   binding's unit (seconds, days). Spells are half-open `[start, end)`,
//!   the convention of networkDynamic (`butts2026networkdynamic`), so
//!   consecutive spells and window boundaries are unambiguous; an edge
//!   with `start` only, or with `end == start`, is an instantaneous
//!   contact active exactly at `start`. In the vocabulary of
//!   `latapy2018stream`, `start` alone gives a contact sequence, `start`
//!   and `end` a link stream, and node spells as well a stream graph
//!   (`holme2012temporal` call the second an interval graph; that name
//!   is avoided because graph theory uses it for something else). Edge
//!   times must be finite; an edge must not end before it starts.
//! - `wave`: an integer edge column for panel data observed at discrete
//!   waves (`snijders2010introduction`); wave `k` is the spell
//!   `[k, k + 1)`, so every temporal operation applies unchanged.
//! - `node_start`, `node_end`: a node's presence spell, missing meaning
//!   unbounded. An edge may only be active while both endpoints are
//!   present (Pajek's consistency rule for time intervals,
//!   `denooy2018exploratory`): `node_start <= start` and
//!   `end <= node_end` for both endpoints, where a contact is active at
//!   its `start` if `start < node_end` or the node's spell is itself the
//!   instant `start`. A node with several presence spells uses the
//!   `node_spells` table instead (roles [`SpellRole`]): rows of `id`,
//!   `start`, `end`, sorted per node, which must not overlap unless
//!   `overlapping_spells` is declared; an edge must then lie inside one
//!   spell of each endpoint. A node absent from the table is present
//!   throughout. Time-varying node attributes live in the `node_changes`
//!   table (roles [`ChangeRole`]): rows of `id` and the finite `time` a
//!   change takes effect; the values stay with the binding.
//! - `x`, `y`: representative coordinates of the nodes in a spatial
//!   network (`barthelemy2011spatial`: nodes embedded in a metric space;
//!   adjacency still comes from the edge table, never from geometry).
//!   Geometry itself stays with the binding; the frame records whether
//!   each table carries geometry and in which coordinate reference
//!   system ([`Geometry`]), and the two must agree.
//! - `ego` (nodes, and optionally edges): the ego whose network a row
//!   belongs to, so that one frame holds a collection of ego networks
//!   (`marsden2002egocentric`, `everett2005ego`, `perry2018egocentric`)
//!   and measures run per ego. Text or integer, no missing values. An
//!   edge's ego is its own column when set, else its source's, and both
//!   endpoints must belong to it. A node whose identifier is an ego
//!   identifier is that ego and must belong to its own network; the
//!   declaration `egos_are_nodes` requires every ego to be a node.
//! - `node_missing`, `edge_missing`: logical columns marking
//!   non-respondents and unobserved pairs, as distinct from absent ties
//!   (`kossinets2006effects`; network's `na` attribute,
//!   `butts2008network`). A row with `edge_missing` true is kept and
//!   counted, carries its relation, layers, and spell like any row, but
//!   is not a tie of the built graph; a pair that is both observed and
//!   unobserved (one key, both flags) is an error. A non-respondent's
//!   unreported ties are implied by rule and never materialized. Measures
//!   that refuse unresolved missingness read the counts in the summary.
//!
//! # Declared properties
//!
//! [`Declarations`] state whether the network is directed, and optionally
//! whether loops and parallel edges are allowed, whether edges must join
//! different modes, whether egos are nodes, whether interlayer edges are
//! couplings only, and whether node spells may overlap. A declaration is
//! enforced: violating rows are errors, as in `as.network.data.frame()`
//! of `butts2008network`, and contradictory declarations are errors
//! rather than silent overrides. A property left undeclared is observed
//! and reported as a note.
//!
//! # Self-loops and parallel edges
//!
//! A self-loop is an edge whose endpoints are the same node in the same
//! layer; loops are kept (`csardi2006igraph` and NetworkX,
//! `hagberg2008exploring`, allow them, and `newman2018networks`,
//! chapter 6, defines multigraphs with self-edges and multiedges;
//! `wasserman1994social` define graphs without them), counted, and each
//! measure documents how it treats them. A coupling edge is not a loop.
//! Two edges are parallel when they have the same endpoints with the
//! same layers (ordered in a directed network, unordered otherwise)
//! **and** the same values of every role that distinguishes edges: the
//! relation and the time spell. A friendship tie and an advice tie
//! between the same two nodes are not parallel, and neither are two
//! spells of one tie or one tie in two layers. A network with neither
//! loops nor parallel edges is simple. [`multiplicity`] implements the
//! rule, with `simplify` returning an index mapping so bindings can
//! aggregate attributes.
//!
//! # Construction
//!
//! [`build()`] validates and then builds the CSR [`Graph`](crate::Graph):
//! node `i` of the graph is row `i` of the node table, followed by nodes
//! added from the edges in order of first appearance; the graph's edges
//! are the observed rows of the edge table (those not flagged
//! `edge_missing`) in row order, so without missing rows edge `e` is row
//! `e`. Both mappings are returned so bindings can align results with
//! their tables; [`build::missing_pairs`] lists the unobserved pairs, and
//! [`build::layer_subgraph`] and [`build::ego_subgraph`] build the graph
//! of one layer or one ego with their row mappings.
//!
//! # Diagnostics
//!
//! Validation never stops at the first problem: it returns every error and
//! note it can establish ([`Report`]), each with a stable code, a message,
//! and the offending rows (the first [`diagnostics::ROW_CAP`], plus the
//! count), so bindings phrase row lists their own way. The codes are
//! listed on [`Code`].
//!
//! # Performance
//!
//! Validation and construction take O(n + m) time (plus O(s log s) to
//! sort the s rows of the node spells table): identifiers are hashed
//! once, parallel edges are found by sorting the edge key, and column
//! scans run in parallel chunks through [`crate::parallel`] with results
//! merged in row order, so the diagnostics do not depend on the thread
//! count.

pub mod build;
pub mod columns;
pub mod diagnostics;
pub mod multiplicity;
pub mod roles;
pub mod summary;
pub mod validate;

pub use build::{
    build, build_validated, ego_subgraph, layer_subgraph, missing_pairs, subgraph, Built, EdgeMap,
    NodeMap, Subgraph,
};
pub use columns::{Column, Kind, Mask, Strings, Table, Which};
pub use diagnostics::{Code, Diagnostic, Report, Severity};
pub use multiplicity::{simplify, Combine, Simplified};
pub use roles::{
    ChangeRole, CoordinateKind, Declaration, Declarations, EdgeRole, Geometry, LayerRole, NodeRole,
    SpellRole, WeightType,
};
pub use summary::{summarize, Summary, TimeKind};
pub use validate::{validate, EdgeClass, Egos, Layers, NodeSpells, Resolved, Validation};

use crate::error::{Error, Result};

/// The tables, their roles, and the declared properties of the network.
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    nodes: Option<Table>,
    edges: Table,
    layers: Option<Table>,
    node_spells: Option<Table>,
    node_changes: Option<Table>,
    node_roles: [Option<usize>; NodeRole::ALL.len()],
    edge_roles: [Option<usize>; EdgeRole::ALL.len()],
    layer_roles: [Option<usize>; LayerRole::ALL.len()],
    spell_roles: [Option<usize>; SpellRole::ALL.len()],
    change_roles: [Option<usize>; ChangeRole::ALL.len()],
    weight_type: Option<WeightType>,
    declarations: Declarations,
    geometry: Geometry,
}

impl Frame {
    /// An empty frame without a node table.
    pub fn new(directed: bool) -> Frame {
        Frame {
            nodes: None,
            edges: Table::new(),
            layers: None,
            node_spells: None,
            node_changes: None,
            node_roles: [None; NodeRole::ALL.len()],
            edge_roles: [None; EdgeRole::ALL.len()],
            layer_roles: [None; LayerRole::ALL.len()],
            spell_roles: [None; SpellRole::ALL.len()],
            change_roles: [None; ChangeRole::ALL.len()],
            weight_type: None,
            declarations: Declarations {
                directed,
                ..Declarations::default()
            },
            geometry: Geometry::default(),
        }
    }

    /// Gives the frame a node table with `nrow` rows and no columns yet
    /// (the row count is checked against every node column added later).
    pub fn set_node_rows(&mut self, nrow: usize) -> Result<()> {
        match &self.nodes {
            Some(t) if t.ncol() > 0 && t.nrow() != nrow => Err(Error::LengthMismatch {
                what: "node table",
                expected: t.nrow(),
                found: nrow,
            }),
            _ => {
                self.nodes = Some(Table::with_rows(nrow));
                Ok(())
            }
        }
    }

    /// Adds a column to `which`, creating an optional table if needed.
    pub fn add_column(&mut self, which: Which, name: &str, column: Column) -> Result<usize> {
        let table = match which {
            Which::Nodes => self.nodes.get_or_insert_with(Table::new),
            Which::Edges => &mut self.edges,
            Which::Layers => self.layers.get_or_insert_with(Table::new),
            Which::NodeSpells => self.node_spells.get_or_insert_with(Table::new),
            Which::NodeChanges => self.node_changes.get_or_insert_with(Table::new),
        };
        table.push(which, name, column)
    }

    /// Adds a node column, creating the node table if needed.
    pub fn add_node_column(&mut self, name: &str, column: Column) -> Result<usize> {
        self.add_column(Which::Nodes, name, column)
    }

    /// Adds an edge column.
    pub fn add_edge_column(&mut self, name: &str, column: Column) -> Result<usize> {
        self.add_column(Which::Edges, name, column)
    }

    /// The table `which`, if the frame has it (the edge table always
    /// exists).
    pub fn table(&self, which: Which) -> Option<&Table> {
        match which {
            Which::Nodes => self.nodes.as_ref(),
            Which::Edges => Some(&self.edges),
            Which::Layers => self.layers.as_ref(),
            Which::NodeSpells => self.node_spells.as_ref(),
            Which::NodeChanges => self.node_changes.as_ref(),
        }
    }

    fn position(&self, which: Which, name: &str) -> Result<usize> {
        self.table(which)
            .and_then(|t| t.position(name))
            .ok_or_else(|| Error::InvalidArgument(format!("`{which}` has no column `{name}`")))
    }

    /// Assigns a node role to the column called `name`.
    pub fn set_node_role(&mut self, role: NodeRole, name: &str) -> Result<()> {
        let pos = self.position(Which::Nodes, name)?;
        self.node_roles[role_index(&NodeRole::ALL, role)] = Some(pos);
        Ok(())
    }

    /// Assigns an edge role to the column called `name`.
    pub fn set_edge_role(&mut self, role: EdgeRole, name: &str) -> Result<()> {
        let pos = self.position(Which::Edges, name)?;
        self.edge_roles[role_index(&EdgeRole::ALL, role)] = Some(pos);
        Ok(())
    }

    /// Assigns a layer table role to the column called `name`.
    pub fn set_layer_role(&mut self, role: LayerRole, name: &str) -> Result<()> {
        let pos = self.position(Which::Layers, name)?;
        self.layer_roles[role_index(&LayerRole::ALL, role)] = Some(pos);
        Ok(())
    }

    /// Assigns a node spells table role to the column called `name`.
    pub fn set_spell_role(&mut self, role: SpellRole, name: &str) -> Result<()> {
        let pos = self.position(Which::NodeSpells, name)?;
        self.spell_roles[role_index(&SpellRole::ALL, role)] = Some(pos);
        Ok(())
    }

    /// Assigns a node changes table role to the column called `name`.
    pub fn set_change_role(&mut self, role: ChangeRole, name: &str) -> Result<()> {
        let pos = self.position(Which::NodeChanges, name)?;
        self.change_roles[role_index(&ChangeRole::ALL, role)] = Some(pos);
        Ok(())
    }

    /// Removes a node role.
    pub fn clear_node_role(&mut self, role: NodeRole) {
        self.node_roles[role_index(&NodeRole::ALL, role)] = None;
    }

    /// Removes an edge role.
    pub fn clear_edge_role(&mut self, role: EdgeRole) {
        self.edge_roles[role_index(&EdgeRole::ALL, role)] = None;
    }

    /// Removes a layer table role.
    pub fn clear_layer_role(&mut self, role: LayerRole) {
        self.layer_roles[role_index(&LayerRole::ALL, role)] = None;
    }

    /// Removes a node spells table role.
    pub fn clear_spell_role(&mut self, role: SpellRole) {
        self.spell_roles[role_index(&SpellRole::ALL, role)] = None;
    }

    /// Removes a node changes table role.
    pub fn clear_change_role(&mut self, role: ChangeRole) {
        self.change_roles[role_index(&ChangeRole::ALL, role)] = None;
    }

    /// Declares the meaning of the weights.
    pub fn set_weight_type(&mut self, weight_type: Option<WeightType>) {
        self.weight_type = weight_type;
    }

    /// The declared meaning of the weights.
    pub fn weight_type(&self) -> Option<WeightType> {
        self.weight_type
    }

    /// The declared properties.
    pub fn declarations(&self) -> &Declarations {
        &self.declarations
    }

    /// The declared properties, for editing.
    pub fn declarations_mut(&mut self) -> &mut Declarations {
        &mut self.declarations
    }

    /// Whether edges are directed.
    pub fn is_directed(&self) -> bool {
        self.declarations.directed
    }

    /// The geometry flags.
    pub fn geometry(&self) -> &Geometry {
        &self.geometry
    }

    /// The geometry flags, for editing.
    pub fn geometry_mut(&mut self) -> &mut Geometry {
        &mut self.geometry
    }

    /// The node table, if the frame has one.
    pub fn nodes(&self) -> Option<&Table> {
        self.nodes.as_ref()
    }

    /// The edge table.
    pub fn edges(&self) -> &Table {
        &self.edges
    }

    /// The layer table, if the frame has one.
    pub fn layers(&self) -> Option<&Table> {
        self.layers.as_ref()
    }

    /// The node spells table, if the frame has one.
    pub fn node_spells(&self) -> Option<&Table> {
        self.node_spells.as_ref()
    }

    /// The node changes table, if the frame has one.
    pub fn node_changes(&self) -> Option<&Table> {
        self.node_changes.as_ref()
    }

    /// Rows of the node table (0 without one).
    pub fn node_rows(&self) -> usize {
        self.nodes.as_ref().map_or(0, |t| t.nrow())
    }

    /// Rows of the edge table.
    pub fn edge_rows(&self) -> usize {
        self.edges.nrow()
    }

    /// Rows of the table `which` (0 for an absent table).
    pub fn rows(&self, which: Which) -> usize {
        self.table(which).map_or(0, |t| t.nrow())
    }

    /// The column holding a node role, with its name.
    pub fn node_role(&self, role: NodeRole) -> Option<(&str, &Column)> {
        let t = self.nodes.as_ref()?;
        let i = self.node_roles[role_index(&NodeRole::ALL, role)]?;
        Some((&t.names()[i], t.column(i)))
    }

    /// The column holding an edge role, with its name.
    pub fn edge_role(&self, role: EdgeRole) -> Option<(&str, &Column)> {
        let i = self.edge_roles[role_index(&EdgeRole::ALL, role)]?;
        Some((&self.edges.names()[i], self.edges.column(i)))
    }

    /// The column holding a layer table role, with its name.
    pub fn layer_role(&self, role: LayerRole) -> Option<(&str, &Column)> {
        let t = self.layers.as_ref()?;
        let i = self.layer_roles[role_index(&LayerRole::ALL, role)]?;
        Some((&t.names()[i], t.column(i)))
    }

    /// The column holding a node spells table role, with its name.
    pub fn spell_role(&self, role: SpellRole) -> Option<(&str, &Column)> {
        let t = self.node_spells.as_ref()?;
        let i = self.spell_roles[role_index(&SpellRole::ALL, role)]?;
        Some((&t.names()[i], t.column(i)))
    }

    /// The column holding a node changes table role, with its name.
    pub fn change_role(&self, role: ChangeRole) -> Option<(&str, &Column)> {
        let t = self.node_changes.as_ref()?;
        let i = self.change_roles[role_index(&ChangeRole::ALL, role)]?;
        Some((&t.names()[i], t.column(i)))
    }

    /// Whether a node role is assigned.
    pub fn has_node_role(&self, role: NodeRole) -> bool {
        self.node_role(role).is_some()
    }

    /// Whether an edge role is assigned.
    pub fn has_edge_role(&self, role: EdgeRole) -> bool {
        self.edge_role(role).is_some()
    }

    /// Validates the frame. See [`validate()`].
    pub fn validate(&self, threads: usize) -> Validation {
        validate(self, threads)
    }

    /// Validates the frame and builds the graph. See [`build()`].
    pub fn build(&self, threads: usize) -> Result<Built> {
        build(self, threads)
    }
}

fn role_index<T: PartialEq + Copy>(all: &[T], role: T) -> usize {
    all.iter()
        .position(|r| *r == role)
        .expect("every role is listed")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roles_point_at_columns() {
        let mut f = Frame::new(false);
        f.add_edge_column("from", Column::Int(vec![1, 2], None))
            .unwrap();
        f.add_edge_column("to", Column::Int(vec![2, 3], None))
            .unwrap();
        assert!(f.set_edge_role(EdgeRole::Source, "nope").is_err());
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        assert_eq!(f.edge_role(EdgeRole::Source).unwrap().0, "from");
        assert!(!f.has_edge_role(EdgeRole::Weight));
        f.clear_edge_role(EdgeRole::Source);
        assert!(!f.has_edge_role(EdgeRole::Source));
        assert!(f.nodes().is_none());
        assert!(f.set_node_role(NodeRole::Id, "id").is_err());
        f.add_node_column("id", Column::Int(vec![1, 2, 3], None))
            .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        assert_eq!(f.node_rows(), 3);
        assert!(f.set_node_rows(4).is_err());
        f.clear_node_role(NodeRole::Id);
        assert!(!f.has_node_role(NodeRole::Id));
        let mut g = Frame::new(true);
        g.set_node_rows(5).unwrap();
        assert_eq!(g.node_rows(), 5);
        assert!(g.add_node_column("id", Column::Int(vec![1], None)).is_err());
        g.set_node_rows(1).unwrap();
        g.add_node_column("id", Column::Int(vec![1], None)).unwrap();
        g.set_weight_type(Some(WeightType::Distance));
        assert_eq!(g.weight_type(), Some(WeightType::Distance));
        g.declarations_mut().loops = Some(false);
        assert_eq!(g.declarations().loops, Some(false));
        g.geometry_mut().nodes = true;
        assert!(g.geometry().nodes);
        assert!(g.is_directed());
    }

    #[test]
    fn long_tables_and_their_roles() {
        let mut f = Frame::new(false);
        assert!(f.layers().is_none() && f.table(Which::NodeSpells).is_none());
        assert_eq!(f.rows(Which::Layers), 0);
        assert!(f.set_layer_role(LayerRole::Layer, "layer").is_err());
        f.add_column(Which::Layers, "layer", Column::Int(vec![1, 2], None))
            .unwrap();
        f.set_layer_role(LayerRole::Layer, "layer").unwrap();
        assert_eq!(f.layer_role(LayerRole::Layer).unwrap().0, "layer");
        assert_eq!(f.rows(Which::Layers), 2);
        f.clear_layer_role(LayerRole::Layer);
        assert!(f.layer_role(LayerRole::Layer).is_none());
        f.add_column(Which::NodeSpells, "id", Column::Int(vec![1], None))
            .unwrap();
        f.add_column(Which::NodeSpells, "s", Column::Float(vec![0.0]))
            .unwrap();
        f.set_spell_role(SpellRole::Id, "id").unwrap();
        f.set_spell_role(SpellRole::Start, "s").unwrap();
        assert!(f.set_spell_role(SpellRole::End, "e").is_err());
        assert_eq!(f.spell_role(SpellRole::Start).unwrap().0, "s");
        f.clear_spell_role(SpellRole::Start);
        assert!(f.spell_role(SpellRole::Start).is_none());
        assert!(f.node_spells().is_some());
        f.add_column(Which::NodeChanges, "id", Column::Int(vec![1], None))
            .unwrap();
        f.add_column(Which::NodeChanges, "t", Column::Float(vec![2.0]))
            .unwrap();
        f.set_change_role(ChangeRole::Id, "id").unwrap();
        f.set_change_role(ChangeRole::Time, "t").unwrap();
        assert_eq!(f.change_role(ChangeRole::Time).unwrap().0, "t");
        f.clear_change_role(ChangeRole::Id);
        assert!(f.change_role(ChangeRole::Id).is_none());
        assert_eq!(f.node_changes().unwrap().nrow(), 1);
        assert!(f
            .add_column(Which::NodeChanges, "t", Column::Float(vec![2.0]))
            .is_err());
    }
}
