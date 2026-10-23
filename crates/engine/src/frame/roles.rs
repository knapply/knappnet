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
    /// non-respondent). Provisional; see the module documentation.
    NodeMissing,
}

impl NodeRole {
    /// Every node role, in code order.
    pub const ALL: [NodeRole; 7] = [
        NodeRole::Id,
        NodeRole::Mode,
        NodeRole::NodeStart,
        NodeRole::NodeEnd,
        NodeRole::X,
        NodeRole::Y,
        NodeRole::NodeMissing,
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
    /// Logical: `true` marks a dyad whose status was not observed.
    /// Provisional.
    EdgeMissing,
}

impl EdgeRole {
    /// Every edge role, in code order.
    pub const ALL: [EdgeRole; 9] = [
        EdgeRole::Source,
        EdgeRole::Target,
        EdgeRole::Weight,
        EdgeRole::Relation,
        EdgeRole::Sign,
        EdgeRole::Start,
        EdgeRole::End,
        EdgeRole::Wave,
        EdgeRole::EdgeMissing,
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
}

impl Default for Declarations {
    fn default() -> Declarations {
        Declarations {
            directed: true,
            loops: None,
            multiple: None,
            between_modes_only: None,
            add_nodes: false,
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
}

impl Declaration {
    /// Decodes the C interface codes: 0 = directed, 1 = loops,
    /// 2 = multiple, 3 = between_modes_only, 4 = add_nodes.
    pub fn from_code(code: i32) -> Result<Declaration> {
        match code {
            0 => Ok(Declaration::Directed),
            1 => Ok(Declaration::Loops),
            2 => Ok(Declaration::Multiple),
            3 => Ok(Declaration::BetweenModesOnly),
            4 => Ok(Declaration::AddNodes),
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
        assert!(Declaration::from_code(5).is_err());
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
    }
}
