//! Roles, which give columns their structural meaning, and the declared
//! properties of a network. The module documentation of [`super`] states
//! what each role means and cites its source.

use std::fmt;

use crate::error::{Error, Result};

/// A role a node column can play.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NodeRole {
    /// The node identifier (integer, double with whole values, or text).
    Id,
    /// The node's mode: any discrete column (text, integer, or logical).
    Mode,
    /// The start of the node's presence spell (double).
    NodeStart,
    /// The end of the node's presence spell (double).
    NodeEnd,
    /// Representative x coordinate or longitude (double).
    X,
    /// Representative y coordinate or latitude (double).
    Y,
    /// Logical: `true` marks a node whose ties were not observed (a
    /// non-respondent).
    NodeMissing,
    /// The ego whose network the node belongs to (text or integer).
    Ego,
}

impl NodeRole {
    /// Every node role, in code order.
    pub const ALL: [NodeRole; 8] = [
        NodeRole::Id,
        NodeRole::Mode,
        NodeRole::NodeStart,
        NodeRole::NodeEnd,
        NodeRole::X,
        NodeRole::Y,
        NodeRole::NodeMissing,
        NodeRole::Ego,
    ];

    /// The role's name, as bindings and messages use it.
    pub fn name(self) -> &'static str {
        match self {
            NodeRole::Id => "id",
            NodeRole::Mode => "mode",
            NodeRole::NodeStart => "node_start",
            NodeRole::NodeEnd => "node_end",
            NodeRole::X => "x",
            NodeRole::Y => "y",
            NodeRole::NodeMissing => "node_missing",
            NodeRole::Ego => "ego",
        }
    }

    /// Decodes the C interface codes, which are positions in [`Self::ALL`].
    pub fn from_code(code: i32) -> Result<NodeRole> {
        usize::try_from(code)
            .ok()
            .and_then(|i| NodeRole::ALL.get(i).copied())
            .ok_or_else(|| Error::InvalidArgument(format!("unknown node role code {code}")))
    }
}

/// A role an edge column can play.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EdgeRole {
    /// The edge's source endpoint (same type as the node identifier).
    Source,
    /// The edge's target endpoint.
    Target,
    /// The edge weight (double or integer); its meaning is [`WeightType`].
    Weight,
    /// The relation the edge belongs to: any discrete column.
    Relation,
    /// The edge's sign: an integer column holding -1 and 1, or a logical
    /// column where `true` is positive. Provisional.
    Sign,
    /// Start time of the edge's spell, or the time of a contact (double).
    Start,
    /// End time of the edge's spell (double). Without it every timed edge
    /// is a contact, active exactly at its start.
    End,
    /// Observation wave (integer); wave `k` is the spell `[k, k + 1)`.
    Wave,
    /// Logical: `true` marks a pair whose status was not observed; the
    /// row is kept and counted but is not a tie of the built graph.
    EdgeMissing,
    /// The layer of the source endpoint in a multilayer network: any
    /// discrete column, set together with [`EdgeRole::TargetLayer`].
    SourceLayer,
    /// The layer of the target endpoint.
    TargetLayer,
    /// The ego whose network the edge belongs to; optional alongside the
    /// node role, and must agree with both endpoints.
    Ego,
}

impl EdgeRole {
    /// Every edge role, in code order.
    pub const ALL: [EdgeRole; 12] = [
        EdgeRole::Source,
        EdgeRole::Target,
        EdgeRole::Weight,
        EdgeRole::Relation,
        EdgeRole::Sign,
        EdgeRole::Start,
        EdgeRole::End,
        EdgeRole::Wave,
        EdgeRole::EdgeMissing,
        EdgeRole::SourceLayer,
        EdgeRole::TargetLayer,
        EdgeRole::Ego,
    ];

    /// The role's name, as bindings and messages use it.
    pub fn name(self) -> &'static str {
        match self {
            EdgeRole::Source => "source",
            EdgeRole::Target => "target",
            EdgeRole::Weight => "weight",
            EdgeRole::Relation => "relation",
            EdgeRole::Sign => "sign",
            EdgeRole::Start => "start",
            EdgeRole::End => "end",
            EdgeRole::Wave => "wave",
            EdgeRole::EdgeMissing => "edge_missing",
            EdgeRole::SourceLayer => "source_layer",
            EdgeRole::TargetLayer => "target_layer",
            EdgeRole::Ego => "ego",
        }
    }

    /// Decodes the C interface codes, which are positions in [`Self::ALL`].
    pub fn from_code(code: i32) -> Result<EdgeRole> {
        usize::try_from(code)
            .ok()
            .and_then(|i| EdgeRole::ALL.get(i).copied())
            .ok_or_else(|| Error::InvalidArgument(format!("unknown edge role code {code}")))
    }
}

/// The role of a column of the optional layer table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LayerRole {
    /// The layer identifier: any discrete column, unique and not missing.
    Layer,
}

impl LayerRole {
    /// Every layer role, in code order.
    pub const ALL: [LayerRole; 1] = [LayerRole::Layer];

    /// The role's name.
    pub fn name(self) -> &'static str {
        "layer"
    }

    /// Decodes the C interface codes, which are positions in [`Self::ALL`].
    pub fn from_code(code: i32) -> Result<LayerRole> {
        usize::try_from(code)
            .ok()
            .and_then(|i| LayerRole::ALL.get(i).copied())
            .ok_or_else(|| Error::InvalidArgument(format!("unknown layer role code {code}")))
    }
}

/// The roles of the optional node spells table (one row per presence
/// spell of a node; several rows per node).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpellRole {
    /// The node identifier (same type as the node table's).
    Id,
    /// Start of the spell (double); missing means unbounded.
    Start,
    /// End of the spell (double); missing means unbounded.
    End,
}

impl SpellRole {
    /// Every spell role, in code order.
    pub const ALL: [SpellRole; 3] = [SpellRole::Id, SpellRole::Start, SpellRole::End];

    /// The role's name.
    pub fn name(self) -> &'static str {
        match self {
            SpellRole::Id => "id",
            SpellRole::Start => "start",
            SpellRole::End => "end",
        }
    }

    /// Decodes the C interface codes, which are positions in [`Self::ALL`].
    pub fn from_code(code: i32) -> Result<SpellRole> {
        usize::try_from(code)
            .ok()
            .and_then(|i| SpellRole::ALL.get(i).copied())
            .ok_or_else(|| Error::InvalidArgument(format!("unknown node spell role code {code}")))
    }
}

/// The roles of the optional node changes table (one row per change of
/// a node attribute; the attribute columns stay with the binding).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChangeRole {
    /// The node identifier (same type as the node table's).
    Id,
    /// The time the change takes effect (double, finite).
    Time,
}

impl ChangeRole {
    /// Every change role, in code order.
    pub const ALL: [ChangeRole; 2] = [ChangeRole::Id, ChangeRole::Time];

    /// The role's name.
    pub fn name(self) -> &'static str {
        match self {
            ChangeRole::Id => "id",
            ChangeRole::Time => "time",
        }
    }

    /// Decodes the C interface codes, which are positions in [`Self::ALL`].
    pub fn from_code(code: i32) -> Result<ChangeRole> {
        usize::try_from(code)
            .ok()
            .and_then(|i| ChangeRole::ALL.get(i).copied())
            .ok_or_else(|| Error::InvalidArgument(format!("unknown node change role code {code}")))
    }
}

/// The declared meaning of edge weights (`newman2001scientific`,
/// `opsahl2010node`). There is no default: a weighted network must say
/// which it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeightType {
    /// Larger values bring nodes closer (frequency, intensity, capacity).
    Strength,
    /// Larger values keep nodes apart (cost, length, latency); must be
    /// positive.
    Distance,
}

impl WeightType {
    /// The name used in messages and summaries.
    pub fn name(self) -> &'static str {
        match self {
            WeightType::Strength => "strength",
            WeightType::Distance => "distance",
        }
    }

    /// Decodes the C interface codes: 0 = distance, 1 = strength (the codes
    /// of `weights_as` elsewhere in the interface).
    pub fn from_code(code: i32) -> Result<WeightType> {
        match code {
            0 => Ok(WeightType::Distance),
            1 => Ok(WeightType::Strength),
            _ => Err(Error::InvalidArgument(format!(
                "unknown weight type code {code}"
            ))),
        }
    }
}

impl fmt::Display for WeightType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// What kind of coordinates the representative points are in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoordinateKind {
    /// Not stated.
    #[default]
    Unknown,
    /// Longitude and latitude in degrees; latitude must lie in [-90, 90].
    LonLat,
    /// Planar coordinates in any unit.
    Planar,
}

impl CoordinateKind {
    /// Decodes the C interface codes: 0 = unknown, 1 = lonlat, 2 = planar.
    pub fn from_code(code: i32) -> Result<CoordinateKind> {
        match code {
            0 => Ok(CoordinateKind::Unknown),
            1 => Ok(CoordinateKind::LonLat),
            2 => Ok(CoordinateKind::Planar),
            _ => Err(Error::InvalidArgument(format!(
                "unknown coordinate kind code {code}"
            ))),
        }
    }
}

/// Geometry stays with the binding (R's sf, Python's shapely); the engine
/// records only whether each table carries geometry and in which
/// coordinate reference system, so that it can check the two agree and
/// report them.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Geometry {
    /// Whether the node table carries a geometry column.
    pub nodes: bool,
    /// Whether the edge table carries a geometry column.
    pub edges: bool,
    /// The node geometry's coordinate reference system, as the binding
    /// names it (an EPSG code, WKT, or proj string); `None` when unset.
    pub node_crs: Option<String>,
    /// The edge geometry's coordinate reference system.
    pub edge_crs: Option<String>,
    /// What the representative coordinates (roles `x` and `y`) are.
    pub coordinates: CoordinateKind,
}

/// Properties the user declares about a network. `None` means "not
/// declared": the property is then observed and reported as a note rather
/// than enforced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Declarations {
    /// Whether edges are ordered pairs.
    pub directed: bool,
    /// Whether self-loops are allowed. `Some(false)` makes every loop an
    /// error (network's `loops = FALSE`, `butts2008network`).
    pub loops: Option<bool>,
    /// Whether parallel edges are allowed (network's `multiple`).
    pub multiple: Option<bool>,
    /// Whether every edge must join two different modes (a two-mode or
    /// k-partite network). Requires a mode role.
    pub between_modes_only: Option<bool>,
    /// Whether endpoints absent from the node table are added to it, in
    /// order of first appearance, rather than being errors. Implied when
    /// the frame has no node table.
    pub add_nodes: bool,
    /// Whether every ego identifier is a node identifier (`Some(true)`:
    /// an ego that is not a node is an error). Requires an ego role.
    pub egos_are_nodes: Option<bool>,
    /// Whether every interlayer edge joins one actor to itself
    /// (`Some(true)`: the multiplex special case; an interlayer edge
    /// between two actors is an error). Requires the layer roles.
    pub couplings_only: Option<bool>,
    /// Whether the presence spells of one node may overlap
    /// (`Some(true)`); otherwise overlapping rows of the node spells
    /// table are errors.
    pub overlapping_spells: Option<bool>,
}

impl Default for Declarations {
    fn default() -> Declarations {
        Declarations {
            directed: true,
            loops: None,
            multiple: None,
            between_modes_only: None,
            add_nodes: false,
            egos_are_nodes: None,
            couplings_only: None,
            overlapping_spells: None,
        }
    }
}

/// Which declaration a C interface call sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Declaration {
    /// [`Declarations::directed`].
    Directed,
    /// [`Declarations::loops`].
    Loops,
    /// [`Declarations::multiple`].
    Multiple,
    /// [`Declarations::between_modes_only`].
    BetweenModesOnly,
    /// [`Declarations::add_nodes`].
    AddNodes,
    /// [`Declarations::egos_are_nodes`].
    EgosAreNodes,
    /// [`Declarations::couplings_only`].
    CouplingsOnly,
    /// [`Declarations::overlapping_spells`].
    OverlappingSpells,
}

impl Declaration {
    /// Decodes the C interface codes: 0 = directed, 1 = loops,
    /// 2 = multiple, 3 = between_modes_only, 4 = add_nodes,
    /// 5 = egos_are_nodes, 6 = couplings_only, 7 = overlapping_spells.
    pub fn from_code(code: i32) -> Result<Declaration> {
        match code {
            0 => Ok(Declaration::Directed),
            1 => Ok(Declaration::Loops),
            2 => Ok(Declaration::Multiple),
            3 => Ok(Declaration::BetweenModesOnly),
            4 => Ok(Declaration::AddNodes),
            5 => Ok(Declaration::EgosAreNodes),
            6 => Ok(Declaration::CouplingsOnly),
            7 => Ok(Declaration::OverlappingSpells),
            _ => Err(Error::InvalidArgument(format!(
                "unknown declaration code {code}"
            ))),
        }
    }
}

impl Declarations {
    /// Sets one declaration from a C interface value: for the optional
    /// ones, a negative value clears the declaration.
    pub fn set(&mut self, which: Declaration, value: i32) {
        let flag = if value < 0 { None } else { Some(value != 0) };
        match which {
            Declaration::Directed => self.directed = value > 0,
            Declaration::Loops => self.loops = flag,
            Declaration::Multiple => self.multiple = flag,
            Declaration::BetweenModesOnly => self.between_modes_only = flag,
            Declaration::AddNodes => self.add_nodes = value > 0,
            Declaration::EgosAreNodes => self.egos_are_nodes = flag,
            Declaration::CouplingsOnly => self.couplings_only = flag,
            Declaration::OverlappingSpells => self.overlapping_spells = flag,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_round_trip() {
        for (i, r) in NodeRole::ALL.iter().enumerate() {
            assert_eq!(NodeRole::from_code(i as i32).unwrap(), *r);
        }
        for (i, r) in EdgeRole::ALL.iter().enumerate() {
            assert_eq!(EdgeRole::from_code(i as i32).unwrap(), *r);
        }
        assert!(NodeRole::from_code(-1).is_err());
        assert!(EdgeRole::from_code(99).is_err());
        assert_eq!(WeightType::from_code(1).unwrap(), WeightType::Strength);
        assert!(WeightType::from_code(2).is_err());
        assert_eq!(
            CoordinateKind::from_code(1).unwrap(),
            CoordinateKind::LonLat
        );
        assert!(CoordinateKind::from_code(3).is_err());
        assert!(Declaration::from_code(8).is_err());
        assert_eq!(LayerRole::from_code(0).unwrap(), LayerRole::Layer);
        assert!(LayerRole::from_code(1).is_err());
        assert_eq!(SpellRole::from_code(2).unwrap(), SpellRole::End);
        assert!(SpellRole::from_code(3).is_err());
        assert_eq!(ChangeRole::from_code(1).unwrap(), ChangeRole::Time);
        assert!(ChangeRole::from_code(-1).is_err());
        assert_eq!(
            (
                LayerRole::Layer.name(),
                SpellRole::Id.name(),
                ChangeRole::Time.name()
            ),
            ("layer", "id", "time")
        );
    }

    #[test]
    fn declarations_set_and_clear() {
        let mut d = Declarations::default();
        d.set(Declaration::Loops, 0);
        assert_eq!(d.loops, Some(false));
        d.set(Declaration::Loops, -1);
        assert_eq!(d.loops, None);
        d.set(Declaration::Directed, 0);
        assert!(!d.directed);
        d.set(Declaration::AddNodes, 1);
        assert!(d.add_nodes);
        d.set(Declaration::Multiple, 1);
        d.set(Declaration::BetweenModesOnly, 1);
        assert_eq!((d.multiple, d.between_modes_only), (Some(true), Some(true)));
        d.set(Declaration::EgosAreNodes, 1);
        d.set(Declaration::CouplingsOnly, 0);
        d.set(Declaration::OverlappingSpells, -1);
        assert_eq!(
            (d.egos_are_nodes, d.couplings_only, d.overlapping_spells),
            (Some(true), Some(false), None)
        );
    }
}
