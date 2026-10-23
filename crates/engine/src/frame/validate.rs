//! Validation of a [`Frame`]: every check listed in the module
//! documentation of [`super`], in O(n + m) time, returning a [`Report`]
//! and, when the endpoints resolve, the node index mapping and the other
//! resolved columns that [`super::build()`] and [`super::summary`] use.
//!
//! The checks come in four groups, after `as.network.data.frame()` in
//! `butts2008network` (whose validation of endpoint types, missing
//! endpoints, unknown and duplicate vertices, and the loops, multiple,
//! and bipartite flags is the model) with the temporal, spatial, mode,
//! relation, and weight checks added:
//!
//! 1. structure: endpoint and identifier columns, their types, missing
//!    values, unknown and duplicate identifiers;
//! 2. declared properties: loops, parallel edges, modes;
//! 3. attribute sanity: weights, modes, relations, signs, times,
//!    coordinates, missing-data roles;
//! 4. advisory notes: nodes added, isolates, 0-based identifiers, coerced
//!    identifiers, loops and parallel edges when nothing was declared.
//!
//! Every check that scans a column runs in parallel chunks and merges its
//! rows in row order, so a report is the same for any thread count.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hash, Hasher};

use super::columns::{Column, Kind, Strings, Which};
use super::diagnostics::{Code, Diagnostic, Report, ROW_CAP};
use super::roles::{
    ChangeRole, CoordinateKind, EdgeRole, LayerRole, NodeRole, SpellRole, WeightType,
};
use super::Frame;
use crate::parallel::for_each_item;

/// A discrete column resolved to integer codes: `codes[i]` indexes
/// `labels`, in order of first appearance, or is `u32::MAX` for a
/// missing value.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Discrete {
    /// One code per row.
    pub codes: Vec<u32>,
    /// The distinct values, rendered as text, in order of first appearance.
    pub labels: Vec<String>,
}

/// Parallel-edge groups: edges with the same key (endpoints, relation,
/// spell) share a group.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Groups {
    /// Group of every edge, numbered from 0 in order of first appearance.
    pub group: Vec<u32>,
    /// Number of edges in each group.
    pub sizes: Vec<u32>,
    /// The first edge row of each group.
    pub first: Vec<usize>,
}

/// The layers of a multilayer network, resolved to codes.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Layers {
    /// The layer labels: the layer table's identifiers in row order when
    /// a table is given, else the edge layers in order of first
    /// appearance (source before target within a row).
    pub labels: Vec<String>,
    /// Layer code of every edge's source endpoint.
    pub src: Vec<u32>,
    /// Layer code of every edge's target endpoint.
    pub tgt: Vec<u32>,
    /// Whether the labels come from a layer table.
    pub from_table: bool,
}

/// How an edge of a multilayer network sits among the layers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClass {
    /// Both endpoints in one layer.
    Intralayer,
    /// One actor joined to itself in another layer.
    Coupling,
    /// Two actors in different layers.
    Interlayer,
}

impl EdgeClass {
    /// The C interface code: 0 = intralayer, 1 = coupling, 2 = interlayer.
    pub fn code(self) -> i32 {
        match self {
            EdgeClass::Intralayer => 0,
            EdgeClass::Coupling => 1,
            EdgeClass::Interlayer => 2,
        }
    }
}

impl Layers {
    /// The class of edge `e`, given the resolved endpoints.
    pub fn class(&self, e: usize, src: &[u32], tgt: &[u32]) -> EdgeClass {
        if self.src[e] == self.tgt[e] {
            EdgeClass::Intralayer
        } else if src[e] == tgt[e] {
            EdgeClass::Coupling
        } else {
            EdgeClass::Interlayer
        }
    }
}

/// The egos of a collection of ego networks, resolved to codes.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Egos {
    /// The ego labels, in order of first appearance in the node table.
    pub labels: Vec<String>,
    /// Ego code of every node (`u32::MAX` for a node added from the
    /// edges, which validation reports as an error).
    pub node: Vec<u32>,
    /// Ego code of every edge: its own ego column when set, else its
    /// source's.
    pub edge: Vec<u32>,
    /// For every ego, the node index of the ego itself when the ego
    /// identifier is a node identifier.
    pub node_of: Vec<Option<u32>>,
}

/// Presence spells per node from the node spells table, sorted by start
/// within each node; a node without rows has no entries and is present
/// throughout.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct NodeSpells {
    /// `offsets[i]..offsets[i + 1]` are node `i`'s spells (`n + 1`
    /// values).
    pub offsets: Vec<usize>,
    /// Start of each spell (`-inf` when unbounded).
    pub start: Vec<f64>,
    /// End of each spell (`inf` when unbounded).
    pub end: Vec<f64>,
    /// The table row of each spell.
    pub rows: Vec<usize>,
}

impl NodeSpells {
    /// Whether the spell `[s, e)` (a contact when `e == s`) lies inside
    /// one of node `i`'s spells. Spells are searched by start; when they
    /// may overlap (declared), every spell of the node is tried.
    pub fn contains(&self, i: usize, s: f64, e: f64, overlapping: bool) -> bool {
        let (lo, hi) = (self.offsets[i], self.offsets[i + 1]);
        if lo == hi {
            return true;
        }
        let inside = |a: f64, b: f64| a <= s && e <= b && (e > s || s < b || a == b);
        if overlapping {
            return (lo..hi).any(|k| inside(self.start[k], self.end[k]));
        }
        // The last spell starting at or before s.
        let starts = &self.start[lo..hi];
        let k = starts.partition_point(|&a| a <= s);
        k > 0 && inside(self.start[lo + k - 1], self.end[lo + k - 1])
    }
}

/// What validation resolved, available when the endpoints map to nodes.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    /// Number of nodes: the node table's rows plus the added nodes.
    pub n: usize,
    /// Source node index of every edge.
    pub src: Vec<u32>,
    /// Target node index of every edge.
    pub tgt: Vec<u32>,
    /// Nodes added from the edges, in node order after the table's rows:
    /// the edge row the identifier first appeared in and whether it was
    /// the target (`true`) or the source.
    pub added: Vec<(usize, bool)>,
    /// The identifier type after coercion.
    pub id_kind: Kind,
    /// The mode of every node, when a mode role is set and valid.
    pub modes: Option<Discrete>,
    /// The relation of every edge, when a relation role is set and valid.
    pub relations: Option<Discrete>,
    /// Edge weights as doubles, when the weight role is set and valid.
    pub weights: Option<Vec<f64>>,
    /// Start times, when the network is temporal and the times are valid
    /// (waves give `k`).
    pub start: Option<Vec<f64>>,
    /// End times: the end role, `k + 1` for waves, or the start for
    /// contacts.
    pub end: Option<Vec<f64>>,
    /// Node presence spells `(start, end)` with infinities for unbounded
    /// ends, when node spell roles are set and valid.
    pub node_spells: Option<Vec<(f64, f64)>>,
    /// Presence spells from the node spells table, when the table is
    /// given and valid.
    pub spell_table: Option<NodeSpells>,
    /// The layers, when the layer roles are set and valid.
    pub layers: Option<Layers>,
    /// The egos, when the ego role is set and valid.
    pub egos: Option<Egos>,
    /// Whether each edge row is an unobserved pair, when the
    /// `edge_missing` role is set and valid.
    pub missing: Option<Vec<bool>>,
    /// Parallel-edge groups, when every distinguishing role resolved.
    pub groups: Option<Groups>,
    /// Number of self-loops (same node, same layer).
    pub loops: usize,
    /// Whether edges are directed (the declaration, kept here so that
    /// later operations on the resolved edges need no frame).
    pub directed: bool,
}

impl Resolved {
    /// Whether edge `e` is an observed tie (not flagged missing).
    #[inline]
    pub fn observed(&self, e: usize) -> bool {
        !self.missing.as_ref().is_some_and(|m| m[e])
    }

    /// Number of observed ties.
    pub fn observed_count(&self) -> usize {
        match &self.missing {
            Some(m) => m.iter().filter(|&&x| !x).count(),
            None => self.src.len(),
        }
    }

    /// Whether edge `e` is a self-loop: the same node in the same layer.
    #[inline]
    pub fn is_loop(&self, e: usize) -> bool {
        self.src[e] == self.tgt[e] && self.layers.as_ref().map_or(true, |l| l.src[e] == l.tgt[e])
    }
}

/// The outcome of validation.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Validation {
    /// Every error and note.
    pub report: Report,
    /// The resolved columns, when the endpoints could be mapped to nodes.
    /// It can be present alongside errors in other columns.
    pub resolved: Option<Resolved>,
}

// ---- Hashing -------------------------------------------------------------

/// A multiplicative hasher for identifiers and edge keys: the standard
/// library's hasher is designed to resist collision attacks, which an
/// identifier column does not need, and this one is several times faster
/// on the short keys hashed here. Words are folded with a rotate, xor,
/// and multiply, and the final avalanche is the SplitMix64 mixing
/// function (`steele2014fast`).
#[derive(Default, Clone, Copy)]
struct Fast(u64);

const MULT: u64 = 0x9E37_79B9_7F4A_7C15;

impl Fast {
    #[inline]
    fn add(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(MULT);
    }
}

impl Hasher for Fast {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            let mut w = [0u8; 8];
            w.copy_from_slice(c);
            self.add(u64::from_le_bytes(w));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut w = [0u8; 8];
            w[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(w) ^ ((rest.len() as u64) << 56));
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }

    #[inline]
    fn write_i64(&mut self, i: i64) {
        self.add(i as u64);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        // So that the low bits, which the table indexes by, depend on every
        // input word.
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<Fast>>;

fn fast_map<K, V>(capacity: usize) -> FastMap<K, V> {
    HashMap::with_capacity_and_hasher(capacity, BuildHasherDefault::default())
}

// ---- Row scans -----------------------------------------------------------

const CHUNK: usize = 1 << 16;

/// The rows in `0..len` for which `pred` holds: the first [`ROW_CAP`] in
/// ascending order, and the count. Scans chunks in parallel and merges
/// them in row order.
fn scan_rows(
    len: usize,
    threads: usize,
    pred: impl Fn(usize) -> bool + Sync,
) -> (Vec<usize>, usize) {
    let chunks = (len + CHUNK - 1) / CHUNK;
    if chunks <= 1 || threads == 1 {
        let mut rows = Vec::new();
        let mut count = 0;
        for i in 0..len {
            if pred(i) {
                count += 1;
                if rows.len() < ROW_CAP {
                    rows.push(i);
                }
            }
        }
        return (rows, count);
    }
    let states = for_each_item(
        chunks,
        threads,
        Vec::<(usize, Vec<usize>, usize)>::new,
        |state, chunk| {
            let lo = chunk * CHUNK;
            let hi = (lo + CHUNK).min(len);
            let mut rows = Vec::new();
            let mut count = 0;
            for i in lo..hi {
                if pred(i) {
                    count += 1;
                    if rows.len() < ROW_CAP {
                        rows.push(i);
                    }
                }
            }
            if count > 0 {
                state.push((chunk, rows, count));
            }
        },
    );
    let mut parts: Vec<(usize, Vec<usize>, usize)> = states.into_iter().flatten().collect();
    parts.sort_unstable_by_key(|p| p.0);
    let mut rows = Vec::new();
    let mut count = 0;
    for (_, part_rows, part_count) in parts {
        count += part_count;
        for r in part_rows {
            if rows.len() < ROW_CAP {
                rows.push(r);
            }
        }
    }
    (rows, count)
}

/// Reports `code` for the rows where `pred` holds, if any.
fn report_rows(
    report: &mut Report,
    code: Code,
    table: Which,
    message: impl Into<String>,
    len: usize,
    threads: usize,
    pred: impl Fn(usize) -> bool + Sync,
) -> bool {
    let (rows, count) = scan_rows(len, threads, pred);
    if count > 0 {
        report.push(Diagnostic::capped(code, table, message, rows, count));
        true
    } else {
        false
    }
}

// ---- Identifiers ---------------------------------------------------------

/// An identifier column seen as keys of one kind.
#[derive(Clone, Copy)]
enum IdView<'a> {
    Int(&'a [i64], Option<&'a [bool]>),
    Float(&'a [f64]),
    Str(&'a Strings),
}

impl<'a> IdView<'a> {
    fn of(c: &'a Column) -> Option<IdView<'a>> {
        match c {
            Column::Int(v, m) => Some(IdView::Int(v, m.as_deref())),
            Column::Float(v) => Some(IdView::Float(v)),
            Column::Str(s) => Some(IdView::Str(s)),
            Column::Bool(..) => None,
        }
    }

    fn kind(self) -> Kind {
        match self {
            IdView::Int(..) => Kind::Int,
            IdView::Float(_) => Kind::Float,
            IdView::Str(_) => Kind::Str,
        }
    }

    #[inline]
    fn int(self, i: usize) -> Option<i64> {
        match self {
            IdView::Int(v, m) => (!m.is_some_and(|m| m[i])).then(|| v[i]),
            _ => None,
        }
    }

    #[inline]
    fn float(self, i: usize) -> Option<u64> {
        let x = match self {
            IdView::Int(v, m) => (!m.is_some_and(|m| m[i])).then(|| v[i] as f64)?,
            IdView::Float(v) => v[i],
            IdView::Str(_) => return None,
        };
        // -0.0 and 0.0 are the same identifier.
        (!x.is_nan()).then(|| (x + 0.0).to_bits())
    }

    #[inline]
    fn str(self, i: usize) -> Option<&'a [u8]> {
        match self {
            IdView::Str(s) => s.get(i),
            _ => None,
        }
    }
}

/// A text identifier as a hash key: short values are stored inline so a
/// lookup compares bytes in the table itself rather than following a
/// pointer into the column (the second cache miss that dominated lookups
/// on large networks); longer values are compared through the pointer.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum StrKey<'a> {
    Short([u8; 16], u8),
    Long(&'a [u8]),
}

impl<'a> StrKey<'a> {
    #[inline]
    fn new(bytes: &'a [u8]) -> StrKey<'a> {
        if bytes.len() <= 16 {
            let mut k = [0u8; 16];
            k[..bytes.len()].copy_from_slice(bytes);
            StrKey::Short(k, bytes.len() as u8)
        } else {
            StrKey::Long(bytes)
        }
    }

    #[inline]
    fn is_empty(self) -> bool {
        matches!(self, StrKey::Short(_, 0))
    }
}

/// The result of mapping identifiers to indices.
struct Mapped {
    n: usize,
    src: Vec<u32>,
    tgt: Vec<u32>,
    added: Vec<(usize, bool)>,
    ok: bool,
}

/// The mapping together with the identifier index it was built with.
struct Mapping<'a> {
    mapped: Mapped,
    index: IdIndex<'a>,
}

/// The identifier table after mapping, for looking up the identifiers of
/// the long tables.
enum IdIndex<'a> {
    Int(FastMap<i64, u32>),
    Float(FastMap<u64, u32>),
    Str(FastMap<StrKey<'a>, u32>),
}

/// What a lookup of one identifier found.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Lookup {
    Missing,
    Unknown,
    Node(u32),
}

impl<'a> IdIndex<'a> {
    /// Whether identifiers of `kind` can be looked up in this index
    /// (numbers in a numeric index, text in a text index).
    fn accepts(&self, kind: Kind) -> bool {
        match self {
            IdIndex::Str(_) => kind == Kind::Str,
            _ => kind == Kind::Int || kind == Kind::Float,
        }
    }

    /// Looks up row `i` of `view`.
    fn find(&self, view: IdView<'a>, i: usize) -> Lookup {
        let found = |x: Option<&u32>| x.map_or(Lookup::Unknown, |&n| Lookup::Node(n));
        match self {
            IdIndex::Int(map) => match view {
                IdView::Int(..) => match view.int(i) {
                    None => Lookup::Missing,
                    Some(k) => found(map.get(&k)),
                },
                IdView::Float(v) => {
                    if v[i].is_nan() {
                        Lookup::Missing
                    } else if v[i].fract() == 0.0 && v[i].abs() < 9.0e15 {
                        found(map.get(&(v[i] as i64)))
                    } else {
                        Lookup::Unknown
                    }
                }
                IdView::Str(_) => Lookup::Unknown,
            },
            IdIndex::Float(map) => match view.float(i) {
                None if matches!(view, IdView::Str(_)) => Lookup::Unknown,
                None => Lookup::Missing,
                Some(k) => found(map.get(&k)),
            },
            IdIndex::Str(map) => match view.str(i) {
                None if matches!(view, IdView::Str(_)) => Lookup::Missing,
                None => Lookup::Unknown,
                Some(b) => found(map.get(&StrKey::new(b))),
            },
        }
    }
}

/// Endpoint states after the parallel lookup pass.
const UNKNOWN: u32 = u32::MAX;
const MISSING: u32 = u32::MAX - 1;
const EMPTY: u32 = u32::MAX - 2;

/// Maps node identifiers and edge endpoints to indices with one hash
/// table, reporting missing, duplicate, and unknown identifiers. Node
/// identifiers are inserted in row order; endpoints are looked up in
/// parallel chunks, and the ones not found are resolved in row order
/// afterwards, so nodes added from the edges are numbered by first
/// appearance whatever the thread count.
fn map_ids<K: Hash + Eq + Copy + Sync>(
    frame: &Frame,
    node_key: impl Fn(usize) -> Option<K>,
    src_key: impl Fn(usize) -> Option<K> + Sync,
    tgt_key: impl Fn(usize) -> Option<K> + Sync,
    empty: impl Fn(K) -> bool + Sync,
    threads: usize,
    report: &mut Report,
) -> (Mapped, FastMap<K, u32>) {
    let has_nodes = frame.nodes().is_some();
    let n0 = frame.node_rows();
    let m = frame.edge_rows();
    let add_nodes = !has_nodes || frame.declarations().add_nodes;
    let mut index: FastMap<K, u32> = fast_map(n0 + if add_nodes { m.min(1 << 20) } else { 0 });
    let mut ok = true;

    let mut missing = Vec::new();
    let mut empties = Vec::new();
    let mut dups = Vec::new();
    let (mut n_missing, mut n_empty, mut n_dup) = (0, 0, 0);
    for i in 0..n0 {
        match node_key(i) {
            None => {
                n_missing += 1;
                push_capped(&mut missing, i);
            }
            Some(k) if empty(k) => {
                n_empty += 1;
                push_capped(&mut empties, i);
            }
            Some(k) => {
                if index.insert(k, i as u32).is_some() {
                    n_dup += 1;
                    push_capped(&mut dups, i);
                }
            }
        }
    }
    let id_name = frame
        .node_role(NodeRole::Id)
        .map(|(name, _)| name.to_string())
        .unwrap_or_default();
    if n_missing > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_MISSING_ID,
            Which::Nodes,
            format!("node identifiers (column `{id_name}` of `nodes`) must not be missing"),
            missing,
            n_missing,
        ));
    }
    if n_empty > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_EMPTY_ID,
            Which::Nodes,
            format!("node identifiers (column `{id_name}` of `nodes`) must not be empty strings"),
            empties,
            n_empty,
        ));
    }
    if n_dup > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_DUPLICATE_ID,
            Which::Nodes,
            format!(
                "node identifiers (column `{id_name}` of `nodes`) must be unique; the listed rows repeat an earlier identifier"
            ),
            dups,
            n_dup,
        ));
    }

    // Parallel lookup of every endpoint against the node table.
    let lookup = |key: Option<K>| match key {
        None => MISSING,
        Some(k) if empty(k) => EMPTY,
        Some(k) => index.get(&k).copied().unwrap_or(UNKNOWN),
    };
    let chunks = (m + CHUNK - 1) / CHUNK;
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    if chunks <= 1 || threads == 1 {
        for e in 0..m {
            src.push(lookup(src_key(e)));
            tgt.push(lookup(tgt_key(e)));
        }
    } else {
        let mut parts: Vec<(usize, Vec<u32>, Vec<u32>)> = for_each_item(
            chunks,
            threads,
            Vec::<(usize, Vec<u32>, Vec<u32>)>::new,
            |state, chunk| {
                let lo = chunk * CHUNK;
                let hi = (lo + CHUNK).min(m);
                let s: Vec<u32> = (lo..hi).map(|e| lookup(src_key(e))).collect();
                let t: Vec<u32> = (lo..hi).map(|e| lookup(tgt_key(e))).collect();
                state.push((chunk, s, t));
            },
        )
        .into_iter()
        .flatten()
        .collect();
        parts.sort_unstable_by_key(|p| p.0);
        for (_, s, t) in parts {
            src.extend_from_slice(&s);
            tgt.extend_from_slice(&t);
        }
    }

    // Row-order pass over what the lookup could not settle.
    let mut added = Vec::new();
    let mut missing = Vec::new();
    let mut empties = Vec::new();
    let mut unknown = Vec::new();
    let (mut n_missing, mut n_empty, mut n_unknown) = (0, 0, 0);
    for e in 0..m {
        for is_target in [false, true] {
            let slot = if is_target { &mut tgt[e] } else { &mut src[e] };
            match *slot {
                MISSING => {
                    n_missing += 1;
                    push_capped(&mut missing, e);
                    *slot = 0;
                }
                EMPTY => {
                    n_empty += 1;
                    push_capped(&mut empties, e);
                    *slot = 0;
                }
                UNKNOWN => {
                    let k = if is_target { tgt_key(e) } else { src_key(e) }.expect("looked up");
                    // An earlier row may have added this identifier.
                    if let Some(&i) = index.get(&k) {
                        *slot = i;
                    } else if add_nodes {
                        let i = (n0 + added.len()) as u32;
                        index.insert(k, i);
                        added.push((e, is_target));
                        *slot = i;
                    } else {
                        n_unknown += 1;
                        push_capped(&mut unknown, e);
                        *slot = 0;
                    }
                }
                _ => {}
            }
        }
    }
    let (src_name, tgt_name) = endpoint_names(frame);
    if n_missing > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_MISSING_ENDPOINT,
            Which::Edges,
            format!(
                "edge endpoints (columns `{src_name}` and `{tgt_name}` of `edges`) must not be missing"
            ),
            missing,
            n_missing,
        ));
    }
    if n_empty > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_EMPTY_ID,
            Which::Edges,
            format!(
                "edge endpoints (columns `{src_name}` and `{tgt_name}` of `edges`) must not be empty strings"
            ),
            empties,
            n_empty,
        ));
    }
    if n_unknown > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_UNKNOWN_NODE,
            Which::Edges,
            format!(
                "edge endpoints (columns `{src_name}` and `{tgt_name}` of `edges`) must be identifiers in column `{id_name}` of `nodes`; add the missing nodes, or declare add_nodes to take them from the edges"
            ),
            unknown,
            n_unknown,
        ));
    }
    if has_nodes && !added.is_empty() {
        let rows: Vec<usize> = added.iter().map(|(e, _)| *e).collect();
        report.push(Diagnostic::rows(
            Code::N_NODES_ADDED,
            Which::Edges,
            format!(
                "{} node(s) absent from `nodes` were added from the edge endpoints, in order of first appearance",
                added.len()
            ),
            &rows,
        ));
    }
    let n = n0 + added.len();
    (
        Mapped {
            n,
            src,
            tgt,
            added,
            ok,
        },
        index,
    )
}

#[inline]
fn push_capped(rows: &mut Vec<usize>, i: usize) {
    if rows.len() < ROW_CAP {
        rows.push(i);
    }
}

fn endpoint_names(frame: &Frame) -> (String, String) {
    let name = |r| {
        frame
            .edge_role(r)
            .map(|(n, _)| n.to_string())
            .unwrap_or_default()
    };
    (name(EdgeRole::Source), name(EdgeRole::Target))
}

/// Resolves the identifier types of the endpoint and identifier columns
/// (`None` when they are unusable, after reporting) and whether numeric
/// coercion happened.
fn id_kind(frame: &Frame, report: &mut Report) -> Option<(Kind, bool)> {
    let (src_name, src) = frame.edge_role(EdgeRole::Source)?;
    let (tgt_name, tgt) = frame.edge_role(EdgeRole::Target)?;
    let (sv, tv) = match (IdView::of(src), IdView::of(tgt)) {
        (Some(s), Some(t)) => (s, t),
        _ => {
            report.push(Diagnostic::network(
                Code::E_ENDPOINT_TYPE,
                format!(
                    "edge endpoints (columns `{src_name}` and `{tgt_name}` of `edges`) must be integer, double, or text, not logical"
                ),
            ));
            return None;
        }
    };
    let mut kinds = vec![sv.kind(), tv.kind()];
    let text = kinds.iter().filter(|k| **k == Kind::Str).count();
    if text == 1 {
        report.push(Diagnostic::network(
            Code::E_ENDPOINT_TYPE,
            format!(
                "edge endpoints (columns `{src_name}` and `{tgt_name}` of `edges`) hold different types ({} and {}); convert them to the same type first",
                sv.kind().name(),
                tv.kind().name()
            ),
        ));
        return None;
    }
    if let Some((id_name, id)) = frame.node_role(NodeRole::Id) {
        let iv = match IdView::of(id) {
            Some(v) => v,
            None => {
                report.push(Diagnostic::network(
                    Code::E_ID_TYPE,
                    format!(
                        "node identifiers (column `{id_name}` of `nodes`) must be integer, double, or text, not logical"
                    ),
                ));
                return None;
            }
        };
        if (iv.kind() == Kind::Str) != (text == 2) {
            report.push(Diagnostic::network(
                Code::E_ID_TYPE,
                format!(
                    "node identifiers (column `{id_name}` of `nodes`, {}) must have the same type as the edge endpoints (columns `{src_name}` and `{tgt_name}` of `edges`, {})",
                    iv.kind().name(),
                    sv.kind().name()
                ),
            ));
            return None;
        }
        kinds.push(iv.kind());
    }
    if text == 2 {
        return Some((Kind::Str, false));
    }
    let coerced = kinds.contains(&Kind::Int) && kinds.contains(&Kind::Float);
    if coerced {
        report.push(Diagnostic::network(
            Code::N_IDS_COERCED,
            "identifiers of integer and double type were matched as numbers; give the identifier and endpoint columns one type to avoid the conversion",
        ));
        return Some((Kind::Float, true));
    }
    Some((kinds[0], false))
}

/// Maps identifiers to node indices by their resolved kind.
fn resolve_endpoints<'a>(
    frame: &'a Frame,
    kind: Kind,
    threads: usize,
    report: &mut Report,
) -> Mapping<'a> {
    let src = IdView::of(frame.edge_role(EdgeRole::Source).expect("checked").1).expect("checked");
    let tgt = IdView::of(frame.edge_role(EdgeRole::Target).expect("checked").1).expect("checked");
    let id = frame
        .node_role(NodeRole::Id)
        .and_then(|(_, c)| IdView::of(c));
    let no_nodes = frame.node_rows() == 0 || id.is_none();
    match kind {
        Kind::Int => {
            let (mapped, index) = map_ids(
                frame,
                |i| if no_nodes { None } else { id.unwrap().int(i) },
                |e| src.int(e),
                |e| tgt.int(e),
                |_| false,
                threads,
                report,
            );
            Mapping {
                mapped,
                index: IdIndex::Int(index),
            }
        }
        Kind::Float => {
            let (mapped, index) = map_ids(
                frame,
                |i| if no_nodes { None } else { id.unwrap().float(i) },
                |e| src.float(e),
                |e| tgt.float(e),
                |_| false,
                threads,
                report,
            );
            Mapping {
                mapped,
                index: IdIndex::Float(index),
            }
        }
        _ => {
            let (mapped, index) = map_ids(
                frame,
                |i| {
                    if no_nodes {
                        None
                    } else {
                        id.unwrap().str(i).map(StrKey::new)
                    }
                },
                |e| src.str(e).map(StrKey::new),
                |e| tgt.str(e).map(StrKey::new),
                |k: StrKey| k.is_empty(),
                threads,
                report,
            );
            Mapping {
                mapped,
                index: IdIndex::Str(index),
            }
        }
    }
}

/// Looks every row of an identifier column of a long table up in the
/// node index, reporting missing and unknown identifiers with the
/// table's name; `None` when the column's type cannot be matched.
fn lookup_ids(
    index: &IdIndex<'_>,
    which: Which,
    name: &str,
    column: &Column,
    id_name: &str,
    report: &mut Report,
) -> Option<Vec<u32>> {
    let view = match IdView::of(column) {
        Some(v) if index.accepts(v.kind()) => v,
        _ => {
            report.push(Diagnostic::network(
                Code::E_ID_TYPE,
                format!(
                    "identifiers in column `{name}` of `{which}` ({}) must have the same type as the node identifiers (column `{id_name}` of `nodes`)",
                    column.kind().name()
                ),
            ));
            return None;
        }
    };
    let n = column.len();
    let mut out = vec![0u32; n];
    let mut missing = Vec::new();
    let mut unknown = Vec::new();
    let (mut n_missing, mut n_unknown) = (0, 0);
    for (i, slot) in out.iter_mut().enumerate() {
        match index.find(view, i) {
            Lookup::Node(k) => *slot = k,
            Lookup::Missing => {
                n_missing += 1;
                push_capped(&mut missing, i);
            }
            Lookup::Unknown => {
                n_unknown += 1;
                push_capped(&mut unknown, i);
            }
        }
    }
    let mut ok = true;
    if n_missing > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_MISSING_ENDPOINT,
            which,
            format!("identifiers (column `{name}` of `{which}`) must not be missing"),
            missing,
            n_missing,
        ));
    }
    if n_unknown > 0 {
        ok = false;
        report.push(Diagnostic::capped(
            Code::E_UNKNOWN_NODE,
            which,
            format!(
                "identifiers (column `{name}` of `{which}`) must be identifiers in column `{id_name}` of `nodes`"
            ),
            unknown,
            n_unknown,
        ));
    }
    ok.then_some(out)
}

// ---- Discrete columns ----------------------------------------------------

/// Codes a discrete column (text, integer, or logical) by first
/// appearance. `None` for a double column.
fn discrete(column: &Column) -> Option<(Discrete, Vec<usize>, usize)> {
    let n = column.len();
    let mut codes = vec![u32::MAX; n];
    let mut labels = Vec::new();
    let mut missing = Vec::new();
    let mut n_missing = 0;
    match column {
        Column::Float(_) => return None,
        Column::Bool(v, m) => {
            let mut seen = [u32::MAX; 2];
            for i in 0..n {
                if m.as_ref().is_some_and(|m| m[i]) {
                    n_missing += 1;
                    push_capped(&mut missing, i);
                    continue;
                }
                let slot = v[i] as usize;
                if seen[slot] == u32::MAX {
                    seen[slot] = labels.len() as u32;
                    labels.push(if v[i] { "TRUE" } else { "FALSE" }.to_string());
                }
                codes[i] = seen[slot];
            }
        }
        Column::Int(v, m) => {
            let mut seen: FastMap<i64, u32> = fast_map(16);
            for i in 0..n {
                if m.as_ref().is_some_and(|m| m[i]) {
                    n_missing += 1;
                    push_capped(&mut missing, i);
                    continue;
                }
                let next = labels.len() as u32;
                let code = *seen.entry(v[i]).or_insert(next);
                if code == next {
                    labels.push(v[i].to_string());
                }
                codes[i] = code;
            }
        }
        Column::Str(s) => {
            let mut seen: FastMap<&[u8], u32> = fast_map(16);
            for (i, slot) in codes.iter_mut().enumerate() {
                let bytes = match s.get(i) {
                    Some(b) => b,
                    None => {
                        n_missing += 1;
                        push_capped(&mut missing, i);
                        continue;
                    }
                };
                let next = labels.len() as u32;
                let code = *seen.entry(bytes).or_insert(next);
                if code == next {
                    labels.push(String::from_utf8_lossy(bytes).into_owned());
                }
                *slot = code;
            }
        }
    }
    Some((Discrete { codes, labels }, missing, n_missing))
}

/// Codes the values of several discrete columns of one kind by first
/// appearance (the layers of the two edge columns and the layer table;
/// the egos of the node and edge columns), with lookup-only passes for
/// columns whose values must already be known.
struct Coder<'a> {
    ints: FastMap<i64, u32>,
    strs: FastMap<&'a [u8], u32>,
    bools: [u32; 2],
    labels: Vec<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Coded {
    Missing,
    Unknown,
    Code(u32),
}

impl<'a> Coder<'a> {
    fn new() -> Coder<'a> {
        Coder {
            ints: fast_map(16),
            strs: fast_map(16),
            bools: [u32::MAX; 2],
            labels: Vec::new(),
        }
    }

    /// Codes row `i` of `column`, adding a new label when `insert`.
    fn code(&mut self, column: &'a Column, i: usize, insert: bool) -> Coded {
        let next = self.labels.len() as u32;
        match column {
            Column::Int(v, m) => {
                if m.as_ref().is_some_and(|m| m[i]) {
                    return Coded::Missing;
                }
                match (self.ints.get(&v[i]), insert) {
                    (Some(&c), _) => Coded::Code(c),
                    (None, false) => Coded::Unknown,
                    (None, true) => {
                        self.ints.insert(v[i], next);
                        self.labels.push(v[i].to_string());
                        Coded::Code(next)
                    }
                }
            }
            Column::Str(s) => match s.get(i) {
                None => Coded::Missing,
                Some(b) => match (self.strs.get(b), insert) {
                    (Some(&c), _) => Coded::Code(c),
                    (None, false) => Coded::Unknown,
                    (None, true) => {
                        self.strs.insert(b, next);
                        self.labels.push(String::from_utf8_lossy(b).into_owned());
                        Coded::Code(next)
                    }
                },
            },
            Column::Bool(v, m) => {
                if m.as_ref().is_some_and(|m| m[i]) {
                    return Coded::Missing;
                }
                let slot = v[i] as usize;
                if self.bools[slot] != u32::MAX {
                    Coded::Code(self.bools[slot])
                } else if insert {
                    self.bools[slot] = next;
                    self.labels
                        .push(if v[i] { "TRUE" } else { "FALSE" }.to_string());
                    Coded::Code(next)
                } else {
                    Coded::Unknown
                }
            }
            Column::Float(_) => Coded::Unknown,
        }
    }
}

/// A numeric column as doubles (`None` for text or logical).
fn numeric(column: &Column) -> Option<Vec<f64>> {
    match column {
        Column::Int(v, m) => Some(
            v.iter()
                .enumerate()
                .map(|(i, &x)| {
                    if m.as_ref().is_some_and(|m| m[i]) {
                        f64::NAN
                    } else {
                        x as f64
                    }
                })
                .collect(),
        ),
        Column::Float(v) => Some(v.clone()),
        _ => None,
    }
}

// ---- The checks ----------------------------------------------------------

/// Validates `frame`. `threads` is 0 for every core.
pub fn validate(frame: &Frame, threads: usize) -> Validation {
    let mut report = Report::default();
    let decl = frame.declarations().clone();
    let m = frame.edge_rows();
    let n0 = frame.node_rows();

    // 1. Structure.
    if !frame.has_edge_role(EdgeRole::Source) || !frame.has_edge_role(EdgeRole::Target) {
        report.push(Diagnostic::network(
            Code::E_ENDPOINTS_UNDECLARED,
            "the edge table needs a source and a target column",
        ));
    }
    if frame.nodes().is_some() && !frame.has_node_role(NodeRole::Id) {
        report.push(Diagnostic::network(
            Code::E_ID_UNDECLARED,
            "the node table needs an identifier column",
        ));
    }
    let mut mapping = None;
    let mut id_kind_resolved = Kind::Int;
    if !report.has(Code::E_ENDPOINTS_UNDECLARED) && !report.has(Code::E_ID_UNDECLARED) {
        if let Some((kind, _)) = id_kind(frame, &mut report) {
            id_kind_resolved = kind;
            mapping = Some(resolve_endpoints(frame, kind, threads, &mut report));
        }
    }
    // The long tables need their roles.
    let layer_table_ok = frame.layers().is_none() || frame.layer_role(LayerRole::Layer).is_some();
    if !layer_table_ok {
        report.push(Diagnostic::network(
            Code::E_TABLE_ROLE_UNDECLARED,
            "the layer table needs a layer identifier column",
        ));
    }
    let spell_table_ok = frame.node_spells().is_none()
        || SpellRole::ALL
            .iter()
            .all(|&r| frame.spell_role(r).is_some());
    if !spell_table_ok {
        report.push(Diagnostic::network(
            Code::E_TABLE_ROLE_UNDECLARED,
            "the node spells table needs identifier, start, and end columns",
        ));
    }
    let change_table_ok = frame.node_changes().is_none()
        || ChangeRole::ALL
            .iter()
            .all(|&r| frame.change_role(r).is_some());
    if !change_table_ok {
        report.push(Diagnostic::network(
            Code::E_TABLE_ROLE_UNDECLARED,
            "the node changes table needs identifier and time columns",
        ));
    }

    // Declarations that contradict each other, whatever the data.
    let mode_role = frame.node_role(NodeRole::Mode);
    let has_source_layer = frame.has_edge_role(EdgeRole::SourceLayer);
    let has_target_layer = frame.has_edge_role(EdgeRole::TargetLayer);
    if has_source_layer != has_target_layer {
        report.push(Diagnostic::network(
            Code::E_LAYER_ROLE_PAIR,
            "a multilayer network needs both a source_layer and a target_layer column",
        ));
    }
    if decl.couplings_only.is_some() && !(has_source_layer && has_target_layer) {
        report.push(Diagnostic::network(
            Code::E_DECLARATION_CONFLICT,
            "interlayer edges were declared (couplings_only), but no layer columns are set",
        ));
    }
    let has_node_ego = frame.has_node_role(NodeRole::Ego);
    if decl.egos_are_nodes.is_some() && !has_node_ego {
        report.push(Diagnostic::network(
            Code::E_DECLARATION_CONFLICT,
            "egos were declared (egos_are_nodes), but no ego column is set on `nodes`",
        ));
    }
    if frame.has_edge_role(EdgeRole::Ego) && !has_node_ego {
        report.push(Diagnostic::network(
            Code::E_DECLARATION_CONFLICT,
            "an ego column on `edges` needs an ego column on `nodes`",
        ));
    }
    let has_spell_roles =
        frame.has_node_role(NodeRole::NodeStart) || frame.has_node_role(NodeRole::NodeEnd);
    if has_spell_roles && frame.node_spells().is_some() {
        report.push(Diagnostic::network(
            Code::E_DECLARATION_CONFLICT,
            "node presence is given twice: by node_start or node_end columns and by a node spells table; use one",
        ));
    }
    if decl.between_modes_only == Some(true) {
        if mode_role.is_none() {
            report.push(Diagnostic::network(
                Code::E_MODE_UNDECLARED,
                "edges were declared to join different modes only (between_modes_only), but no mode column is set",
            ));
        }
        if decl.loops == Some(true) {
            report.push(Diagnostic::network(
                Code::E_DECLARATION_CONFLICT,
                "loops cannot be allowed when every edge must join different modes: a self-loop stays within one mode",
            ));
        }
    }
    let has_start = frame.has_edge_role(EdgeRole::Start);
    let has_end = frame.has_edge_role(EdgeRole::End);
    let has_wave = frame.has_edge_role(EdgeRole::Wave);
    if has_end && !has_start {
        report.push(Diagnostic::network(
            Code::E_DECLARATION_CONFLICT,
            "an end time column needs a start time column",
        ));
    }
    if has_wave && (has_start || has_end) {
        report.push(Diagnostic::network(
            Code::E_DECLARATION_CONFLICT,
            "a wave column and start or end time columns cannot both be set; waves are spells [k, k + 1)",
        ));
    }

    // 3a. Modes (needed for the within-mode check below).
    let mut modes = None;
    if let Some((name, col)) = mode_role {
        match discrete(col) {
            None => report.push(Diagnostic::network(
                Code::E_MODE_TYPE,
                format!("the mode column `{name}` of `nodes` must be text, integer, or logical, not double"),
            )),
            Some((d, missing, count)) => {
                if count > 0 {
                    report.push(Diagnostic::capped(
                        Code::E_MODE_MISSING,
                        Which::Nodes,
                        format!("every node needs a mode; column `{name}` of `nodes` has missing values"),
                        missing,
                        count,
                    ));
                } else {
                    modes = Some(d);
                }
            }
        }
    }

    // 3b. Relations.
    let mut relations = None;
    if let Some((name, col)) = frame.edge_role(EdgeRole::Relation) {
        match discrete(col) {
            None => report.push(Diagnostic::network(
                Code::E_RELATION_TYPE,
                format!("the relation column `{name}` of `edges` must be text, integer, or logical, not double"),
            )),
            Some((d, missing, count)) => {
                if count > 0 {
                    report.push(Diagnostic::capped(
                        Code::E_RELATION_MISSING,
                        Which::Edges,
                        format!("every edge needs a relation; column `{name}` of `edges` has missing values"),
                        missing,
                        count,
                    ));
                } else {
                    relations = Some(d);
                }
            }
        }
    }

    // 3c. Signs.
    let signed = frame.has_edge_role(EdgeRole::Sign);
    if let Some((name, col)) = frame.edge_role(EdgeRole::Sign) {
        match col {
            Column::Int(v, mask) => {
                report_rows(
                    &mut report,
                    Code::E_SIGN_VALUE,
                    Which::Edges,
                    format!("signs (column `{name}` of `edges`) must be -1 or 1 and not missing"),
                    m,
                    threads,
                    |e| mask.as_ref().is_some_and(|m| m[e]) || (v[e] != 1 && v[e] != -1),
                );
            }
            Column::Bool(_, mask) => {
                if let Some(mask) = mask {
                    report_rows(
                        &mut report,
                        Code::E_SIGN_VALUE,
                        Which::Edges,
                        format!("signs (column `{name}` of `edges`) must not be missing"),
                        m,
                        threads,
                        |e| mask[e],
                    );
                }
            }
            _ => report.push(Diagnostic::network(
                Code::E_SIGN_TYPE,
                format!(
                    "the sign column `{name}` of `edges` must be integer (-1 and 1) or logical, not {}",
                    col.kind().name()
                ),
            )),
        }
    }

    // 3d. Weights.
    let mut weights = None;
    if let Some((name, col)) = frame.edge_role(EdgeRole::Weight) {
        match numeric(col) {
            None => report.push(Diagnostic::network(
                Code::E_WEIGHT_TYPE,
                format!(
                    "the weight column `{name}` of `edges` must be numeric, not {}",
                    col.kind().name()
                ),
            )),
            Some(w) => {
                let mut ok = true;
                let weight_type = match frame.weight_type() {
                    Some(t) => Some(t),
                    None => {
                        report.push(Diagnostic::network(
                            Code::E_WEIGHT_TYPE_MISSING,
                            format!(
                                "weights (column `{name}` of `edges`) need a declared meaning: strength (larger brings nodes closer) or distance (larger keeps them apart)"
                            ),
                        ));
                        None
                    }
                };
                if report_rows(
                    &mut report,
                    Code::E_WEIGHT_NOT_FINITE,
                    Which::Edges,
                    format!("weights (column `{name}` of `edges`) must not be missing, NaN, or infinite"),
                    m,
                    threads,
                    |e| !w[e].is_finite(),
                ) {
                    ok = false;
                }
                if weight_type == Some(WeightType::Distance)
                    && report_rows(
                        &mut report,
                        Code::E_DISTANCE_NOT_POSITIVE,
                        Which::Edges,
                        format!("weights (column `{name}` of `edges`) are distances and must be positive"),
                        m,
                        threads,
                        |e| w[e].is_finite() && w[e] <= 0.0,
                    )
                {
                    ok = false;
                }
                if signed
                    && report_rows(
                        &mut report,
                        Code::E_WEIGHT_NEGATIVE,
                        Which::Edges,
                        format!(
                            "weights (column `{name}` of `edges`) of a signed network are magnitudes and must not be negative; the sign column carries the sign"
                        ),
                        m,
                        threads,
                        |e| w[e].is_finite() && w[e] < 0.0,
                    )
                {
                    ok = false;
                }
                if ok && weight_type.is_some() {
                    weights = Some(w);
                }
            }
        }
    }

    // 3e. Times.
    let (mut start, mut end) = (None, None);
    let mut times_ok = true;
    if has_start && !has_wave {
        let (sname, scol) = frame.edge_role(EdgeRole::Start).expect("checked");
        match numeric(scol) {
            None => {
                times_ok = false;
                report.push(Diagnostic::network(
                    Code::E_TIME_TYPE,
                    format!(
                        "start times (column `{sname}` of `edges`) must be numeric, not {}",
                        scol.kind().name()
                    ),
                ));
            }
            Some(s) => {
                if report_rows(
                    &mut report,
                    Code::E_TIME_MISSING,
                    Which::Edges,
                    format!("start times (column `{sname}` of `edges`) must not be missing"),
                    m,
                    threads,
                    |e| s[e].is_nan(),
                ) {
                    times_ok = false;
                }
                if report_rows(
                    &mut report,
                    Code::E_TIME_NOT_FINITE,
                    Which::Edges,
                    format!("start times (column `{sname}` of `edges`) must be finite"),
                    m,
                    threads,
                    |e| s[e].is_infinite(),
                ) {
                    times_ok = false;
                }
                start = Some(s);
            }
        }
        if let Some((ename, ecol)) = frame.edge_role(EdgeRole::End) {
            match numeric(ecol) {
                None => {
                    times_ok = false;
                    report.push(Diagnostic::network(
                        Code::E_TIME_TYPE,
                        format!(
                            "end times (column `{ename}` of `edges`) must be numeric, not {}",
                            ecol.kind().name()
                        ),
                    ));
                }
                Some(e_) => {
                    if report_rows(
                        &mut report,
                        Code::E_TIME_MISSING,
                        Which::Edges,
                        format!("end times (column `{ename}` of `edges`) must not be missing"),
                        m,
                        threads,
                        |e| e_[e].is_nan(),
                    ) {
                        times_ok = false;
                    }
                    if report_rows(
                        &mut report,
                        Code::E_TIME_NOT_FINITE,
                        Which::Edges,
                        format!("end times (column `{ename}` of `edges`) must be finite"),
                        m,
                        threads,
                        |e| e_[e].is_infinite(),
                    ) {
                        times_ok = false;
                    }
                    if let Some(s) = &start {
                        if report_rows(
                            &mut report,
                            Code::E_SPELL_ORDER,
                            Which::Edges,
                            format!(
                                "edges must not end before they start (columns `{sname}` and `{ename}` of `edges`)"
                            ),
                            m,
                            threads,
                            |e| s[e].is_finite() && e_[e].is_finite() && e_[e] < s[e],
                        ) {
                            times_ok = false;
                        }
                    }
                    end = Some(e_);
                }
            }
        }
    }
    if has_wave && !has_start && !has_end {
        let (wname, wcol) = frame.edge_role(EdgeRole::Wave).expect("checked");
        let whole = match wcol {
            Column::Int(..) => true,
            Column::Float(v) => v
                .iter()
                .all(|x| x.is_nan() || (x.is_finite() && x.fract() == 0.0)),
            _ => false,
        };
        match (whole, numeric(wcol)) {
            (true, Some(w)) => {
                if report_rows(
                    &mut report,
                    Code::E_WAVE_MISSING,
                    Which::Edges,
                    format!(
                        "every edge needs a wave; column `{wname}` of `edges` has missing values"
                    ),
                    m,
                    threads,
                    |e| w[e].is_nan(),
                ) {
                    times_ok = false;
                }
                end = Some(w.iter().map(|k| k + 1.0).collect());
                start = Some(w);
            }
            _ => {
                times_ok = false;
                report.push(Diagnostic::network(
                    Code::E_WAVE_TYPE,
                    format!("the wave column `{wname}` of `edges` must hold whole numbers"),
                ));
            }
        }
    }
    if !times_ok {
        start = None;
        end = None;
    }

    // Node spells.
    let mut node_spells = None;
    let ns_role = frame.node_role(NodeRole::NodeStart);
    let ne_role = frame.node_role(NodeRole::NodeEnd);
    if ns_role.is_some() || ne_role.is_some() {
        let mut spells_ok = true;
        let mut bound = |role: Option<(&str, &Column)>, what: &str, fill: f64| -> Vec<f64> {
            match role {
                None => vec![fill; n0],
                Some((name, col)) => match numeric(col) {
                    Some(v) => v
                        .into_iter()
                        .map(|x| if x.is_nan() { fill } else { x })
                        .collect(),
                    None => {
                        spells_ok = false;
                        report.push(Diagnostic::network(
                            Code::E_TIME_TYPE,
                            format!("node {what} times (column `{name}` of `nodes`) must be numeric, not {}", col.kind().name()),
                        ));
                        vec![fill; n0]
                    }
                },
            }
        };
        let a = bound(ns_role, "start", f64::NEG_INFINITY);
        let b = bound(ne_role, "end", f64::INFINITY);
        if spells_ok
            && report_rows(
                &mut report,
                Code::E_NODE_SPELL_ORDER,
                Which::Nodes,
                "nodes must not leave before they arrive (node end earlier than node start)",
                n0,
                threads,
                |i| b[i] < a[i],
            )
        {
            spells_ok = false;
        }
        if spells_ok {
            node_spells = Some(a.into_iter().zip(b).collect::<Vec<_>>());
        }
    }

    // 3f. Coordinates and geometry.
    let x_role = frame.node_role(NodeRole::X);
    let y_role = frame.node_role(NodeRole::Y);
    if x_role.is_some() != y_role.is_some() {
        report.push(Diagnostic::network(
            Code::E_COORD_ROLE_PAIR,
            "coordinates need both an x (longitude) and a y (latitude) column",
        ));
    } else if let (Some((xn, xc)), Some((yn, yc))) = (x_role, y_role) {
        for (name, col) in [(xn, xc), (yn, yc)] {
            match numeric(col) {
                None => report.push(Diagnostic::network(
                    Code::E_COORD_TYPE,
                    format!(
                        "the coordinate column `{name}` of `nodes` must be numeric, not {}",
                        col.kind().name()
                    ),
                )),
                Some(v) => {
                    report_rows(
                        &mut report,
                        Code::E_COORD_NOT_FINITE,
                        Which::Nodes,
                        format!("the coordinate column `{name}` of `nodes` has missing or non-finite values"),
                        n0,
                        threads,
                        |i| !v[i].is_finite(),
                    );
                    if name == yn && frame.geometry().coordinates == CoordinateKind::LonLat {
                        report_rows(
                            &mut report,
                            Code::E_LATITUDE_RANGE,
                            Which::Nodes,
                            format!("the latitude column `{name}` of `nodes` has values outside [-90, 90]"),
                            n0,
                            threads,
                            |i| v[i].is_finite() && !(-90.0..=90.0).contains(&v[i]),
                        );
                    }
                }
            }
        }
    }
    let geo = frame.geometry();
    if let (true, true, Some(a), Some(b)) = (geo.nodes, geo.edges, &geo.node_crs, &geo.edge_crs) {
        if a != b {
            report.push(Diagnostic::network(
                Code::E_CRS_MISMATCH,
                format!(
                    "`nodes` and `edges` have different coordinate reference systems ({a} and {b}); transform one table before building the network"
                ),
            ));
        }
    }

    // 3g. Missing-data roles.
    let mut flagged = |role_col: Option<(&str, &Column)>, table: Which, code: Code, what: &str| {
        if let Some((name, col)) = role_col {
            match col {
                Column::Bool(v, mask) => {
                    if let Some(mask) = mask {
                        if mask.iter().any(|&x| x) {
                            report.push(Diagnostic::network(
                                Code::E_MISSING_ROLE_TYPE,
                                format!("the {what} column `{name}` of `{table}` must not have missing values; use FALSE for observed rows"),
                            ));
                            return;
                        }
                    }
                    report_rows(
                        &mut report,
                        code,
                        table,
                        format!(
                            "{} {what}(s) are marked as unobserved in column `{name}` of `{table}`",
                            v.iter().filter(|&&x| x).count()
                        ),
                        v.len(),
                        threads,
                        |i| v[i],
                    );
                }
                _ => report.push(Diagnostic::network(
                    Code::E_MISSING_ROLE_TYPE,
                    format!(
                        "the {what} column `{name}` of `{table}` must be logical, not {}",
                        col.kind().name()
                    ),
                )),
            }
        }
    };
    flagged(
        frame.node_role(NodeRole::NodeMissing),
        Which::Nodes,
        Code::N_MISSING_NODES,
        "node_missing",
    );
    flagged(
        frame.edge_role(EdgeRole::EdgeMissing),
        Which::Edges,
        Code::N_MISSING_EDGES,
        "edge_missing",
    );
    let missing: Option<Vec<bool>> = match frame.edge_role(EdgeRole::EdgeMissing) {
        Some((_, Column::Bool(v, None))) => Some(v.clone()),
        _ => None,
    };

    // 3h. Layers: the two edge columns and the layer table share one
    // coding.
    let mut layers = None;
    if let (Some((sn, sc)), Some((tn, tc))) = (
        frame.edge_role(EdgeRole::SourceLayer),
        frame.edge_role(EdgeRole::TargetLayer),
    ) {
        let table = frame.layer_role(LayerRole::Layer);
        let kind_ok = sc.kind() == tc.kind()
            && sc.kind() != Kind::Float
            && table.map_or(true, |(_, lc)| lc.kind() == sc.kind());
        if !kind_ok {
            report.push(Diagnostic::network(
                Code::E_LAYER_TYPE,
                format!(
                    "layer columns (`{sn}` and `{tn}` of `edges`{}) must be text, integer, or logical, and of one type",
                    table.map_or(String::new(), |(ln, _)| format!(", `{ln}` of `layers`"))
                ),
            ));
        } else {
            let mut coder = Coder::new();
            let mut ok = true;
            if let Some((ln, lc)) = table {
                let mut missing_ids = Vec::new();
                let mut dups = Vec::new();
                let (mut n_missing, mut n_dup) = (0, 0);
                for i in 0..lc.len() {
                    let before = coder.labels.len();
                    match coder.code(lc, i, true) {
                        Coded::Missing => {
                            n_missing += 1;
                            push_capped(&mut missing_ids, i);
                        }
                        Coded::Code(c) if (c as usize) < before => {
                            n_dup += 1;
                            push_capped(&mut dups, i);
                        }
                        _ => {}
                    }
                }
                if n_missing > 0 {
                    ok = false;
                    report.push(Diagnostic::capped(
                        Code::E_LAYER_ID_MISSING,
                        Which::Layers,
                        format!(
                            "layer identifiers (column `{ln}` of `layers`) must not be missing"
                        ),
                        missing_ids,
                        n_missing,
                    ));
                }
                if n_dup > 0 {
                    ok = false;
                    report.push(Diagnostic::capped(
                        Code::E_LAYER_ID_DUPLICATE,
                        Which::Layers,
                        format!("layer identifiers (column `{ln}` of `layers`) must be unique; the listed rows repeat an earlier identifier"),
                        dups,
                        n_dup,
                    ));
                }
            }
            let insert = table.is_none();
            let mut lsrc = vec![0u32; m];
            let mut ltgt = vec![0u32; m];
            let mut missing_rows = Vec::new();
            let mut unknown_rows = Vec::new();
            let (mut n_missing, mut n_unknown) = (0, 0);
            for e in 0..m {
                let a = coder.code(sc, e, insert);
                let b = coder.code(tc, e, insert);
                if a == Coded::Missing || b == Coded::Missing {
                    n_missing += 1;
                    push_capped(&mut missing_rows, e);
                } else if a == Coded::Unknown || b == Coded::Unknown {
                    n_unknown += 1;
                    push_capped(&mut unknown_rows, e);
                }
                if let Coded::Code(c) = a {
                    lsrc[e] = c;
                }
                if let Coded::Code(c) = b {
                    ltgt[e] = c;
                }
            }
            if n_missing > 0 {
                ok = false;
                report.push(Diagnostic::capped(
                    Code::E_LAYER_MISSING,
                    Which::Edges,
                    format!("every edge needs a layer at each endpoint; columns `{sn}` and `{tn}` of `edges` have missing values"),
                    missing_rows,
                    n_missing,
                ));
            }
            if n_unknown > 0 {
                ok = false;
                report.push(Diagnostic::capped(
                    Code::E_LAYER_UNKNOWN,
                    Which::Edges,
                    format!(
                        "edge layers (columns `{sn}` and `{tn}` of `edges`) must be identifiers in column `{}` of `layers`; add the missing layers to the table",
                        table.map(|(ln, _)| ln).unwrap_or_default()
                    ),
                    unknown_rows,
                    n_unknown,
                ));
            }
            if ok {
                layers = Some(Layers {
                    labels: coder.labels,
                    src: lsrc,
                    tgt: ltgt,
                    from_table: table.is_some(),
                });
            }
        }
    }

    // 2. Structure-dependent checks and declared properties.
    let mapping = match mapping {
        Some(map) if map.mapped.ok => map,
        _ => {
            return Validation {
                report,
                resolved: None,
            }
        }
    };
    let Mapping { mapped, index } = mapping;
    let Mapped {
        n, src, tgt, added, ..
    } = mapped;
    let (src_name, tgt_name) = endpoint_names(frame);
    let id_name = frame
        .node_role(NodeRole::Id)
        .map(|(name, _)| name.to_string())
        .unwrap_or_default();

    // Loops: the same node in the same layer. A coupling edge (the same
    // actor in two layers) is not a loop.
    let layer_codes = layers.as_ref().map(|l| (&l.src[..], &l.tgt[..]));
    let (loop_rows, loops) = scan_rows(m, threads, |e| {
        src[e] == tgt[e] && layer_codes.map_or(true, |(ls, lt)| ls[e] == lt[e])
    });
    match decl.loops {
        Some(false) if loops > 0 => report.push(Diagnostic::capped(
            Code::E_LOOP_DECLARED,
            Which::Edges,
            "self-loops were declared absent (loops = FALSE), but the listed edges join a node to itself",
            loop_rows,
            loops,
        )),
        None if loops > 0 => report.push(Diagnostic::capped(
            Code::N_LOOPS,
            Which::Edges,
            format!("the network has {loops} self-loop(s); each measure's documentation says how they count"),
            loop_rows,
            loops,
        )),
        _ => {}
    }

    // Nodes added from the edges have no mode and no presence spell.
    if modes.is_some() && !added.is_empty() {
        let rows: Vec<usize> = added.iter().map(|(e, _)| *e).collect();
        report.push(Diagnostic::rows(
            Code::E_MODE_MISSING,
            Which::Edges,
            format!(
                "every node needs a mode, but the listed edges introduce {} node(s) absent from `nodes`; add them to `nodes` with a mode",
                added.len()
            ),
            &rows,
        ));
        modes = None;
    }
    if let Some(ns) = &mut node_spells {
        ns.resize(n, (f64::NEG_INFINITY, f64::INFINITY));
    }

    // Edges within a mode. A coupling edge joins one actor to itself
    // across layers and is exempt.
    if let (Some(true), Some(md)) = (decl.between_modes_only, &modes) {
        let codes = &md.codes;
        report_rows(
            &mut report,
            Code::E_MODE_WITHIN,
            Which::Edges,
            format!(
                "edges were declared to join different modes only, but the listed edges join two nodes of the same mode (column `{}` of `nodes`)",
                mode_role.map(|(n, _)| n).unwrap_or_default()
            ),
            m,
            threads,
            |e| {
                codes[src[e] as usize] == codes[tgt[e] as usize]
                    && !(src[e] == tgt[e] && layer_codes.is_some_and(|(ls, lt)| ls[e] != lt[e]))
            },
        );
    }

    // Interlayer edges between different actors, when only couplings
    // were declared.
    if let (Some(true), Some((ls, lt))) = (decl.couplings_only, layer_codes) {
        report_rows(
            &mut report,
            Code::E_INTERLAYER_DECLARED,
            Which::Edges,
            "interlayer edges were declared to be couplings only (couplings_only), but the listed edges join two different actors in different layers",
            m,
            threads,
            |e| ls[e] != lt[e] && src[e] != tgt[e],
        );
    }

    // Egos.
    let mut egos = None;
    if let Some((en, ec)) = frame.node_role(NodeRole::Ego) {
        let edge_ego = frame.edge_role(EdgeRole::Ego);
        let kind_ok = matches!(ec.kind(), Kind::Int | Kind::Str)
            && edge_ego.map_or(true, |(_, c)| c.kind() == ec.kind());
        if !kind_ok {
            report.push(Diagnostic::network(
                Code::E_EGO_TYPE,
                format!(
                    "ego columns (`{en}` of `nodes`{}) must be text or integer, and of one type",
                    edge_ego.map_or(String::new(), |(n, _)| format!(", `{n}` of `edges`"))
                ),
            ));
        } else {
            let mut ok = true;
            let mut coder = Coder::new();
            let mut node_codes = vec![u32::MAX; n];
            let (mut rows, mut count) = (Vec::new(), 0);
            for (i, slot) in node_codes.iter_mut().enumerate().take(n0) {
                match coder.code(ec, i, true) {
                    Coded::Code(c) => *slot = c,
                    _ => {
                        count += 1;
                        push_capped(&mut rows, i);
                    }
                }
            }
            if count > 0 {
                ok = false;
                report.push(Diagnostic::capped(
                    Code::E_EGO_MISSING,
                    Which::Nodes,
                    format!("every node needs an ego; column `{en}` of `nodes` has missing values"),
                    rows,
                    count,
                ));
            }
            if !added.is_empty() {
                ok = false;
                let rows: Vec<usize> = added.iter().map(|(e, _)| *e).collect();
                report.push(Diagnostic::rows(
                    Code::E_EGO_MISSING,
                    Which::Edges,
                    format!(
                        "every node needs an ego, but the listed edges introduce {} node(s) absent from `nodes`; add them to `nodes` with an ego",
                        added.len()
                    ),
                    &rows,
                ));
            }
            // Which egos are nodes: every node identifier is looked up
            // among the ego labels.
            let mut node_of = vec![None; coder.labels.len()];
            let mut ego_of_node = vec![None; n0];
            if ok {
                if let Some((_, idc)) = frame.node_role(NodeRole::Id) {
                    for i in 0..n0 {
                        let ego_of_id = match (idc, ec.kind()) {
                            (Column::Int(..), Kind::Int) | (Column::Str(_), Kind::Str) => {
                                match coder.code(idc, i, false) {
                                    Coded::Code(c) => Some(c),
                                    _ => None,
                                }
                            }
                            (Column::Float(v), Kind::Int) => {
                                let x = v[i];
                                if x.is_finite() && x.fract() == 0.0 && x.abs() < 9.0e15 {
                                    coder.ints.get(&(x as i64)).copied()
                                } else {
                                    None
                                }
                            }
                            _ => None,
                        };
                        if let Some(c) = ego_of_id {
                            node_of[c as usize] = Some(i as u32);
                            ego_of_node[i] = Some(c);
                        }
                    }
                }
                let self_rows = report_rows(
                    &mut report,
                    Code::E_EGO_SELF,
                    Which::Nodes,
                    format!(
                        "a node that is an ego must belong to its own ego network, but the listed nodes (column `{en}` of `nodes`) belong to another ego"
                    ),
                    n0,
                    threads,
                    |i| ego_of_node[i].is_some_and(|c| c != node_codes[i]),
                );
                if self_rows {
                    ok = false;
                }
                if decl.egos_are_nodes == Some(true)
                    && report_rows(
                        &mut report,
                        Code::E_EGO_UNKNOWN,
                        Which::Nodes,
                        format!(
                            "egos were declared to be nodes (egos_are_nodes), but the ego of the listed nodes (column `{en}` of `nodes`) is not a node identifier (column `{id_name}`)"
                        ),
                        n0,
                        threads,
                        |i| node_of[node_codes[i] as usize].is_none(),
                    )
                {
                    ok = false;
                }
            }
            // Edge egos: the edge's own column, else the source's; both
            // endpoints must agree.
            let mut edge_codes = vec![u32::MAX; m];
            if ok {
                let (mut rows, mut count) = (Vec::new(), 0);
                for (e, slot) in edge_codes.iter_mut().enumerate() {
                    *slot = match edge_ego {
                        Some((_, c)) => match coder.code(c, e, true) {
                            Coded::Code(k) => k,
                            _ => {
                                count += 1;
                                push_capped(&mut rows, e);
                                u32::MAX
                            }
                        },
                        None => node_codes[src[e] as usize],
                    };
                }
                if count > 0 {
                    ok = false;
                    report.push(Diagnostic::capped(
                        Code::E_EGO_MISSING,
                        Which::Edges,
                        format!(
                            "every edge needs an ego; column `{}` of `edges` has missing values",
                            edge_ego.map(|(n, _)| n).unwrap_or_default()
                        ),
                        rows,
                        count,
                    ));
                }
            }
            if ok
                && report_rows(
                    &mut report,
                    Code::E_EGO_MISMATCH,
                    Which::Edges,
                    "an edge must join two nodes of one ego and carry that ego; the listed edges join nodes of different egos, or carry an ego that differs from their endpoints'",
                    m,
                    threads,
                    |e| {
                        let c = edge_codes[e];
                        node_codes[src[e] as usize] != c || node_codes[tgt[e] as usize] != c
                    },
                )
            {
                ok = false;
            }
            if ok {
                egos = Some(Egos {
                    labels: coder.labels,
                    node: node_codes,
                    edge: edge_codes,
                    node_of,
                });
            }
        }
    }

    // The node spells table: identifiers looked up, spells ordered, then
    // sorted per node and checked for overlap.
    let mut spell_table = None;
    if spell_table_ok && frame.node_spells().is_some() {
        let (idn, idc) = frame.spell_role(SpellRole::Id).expect("checked");
        let rows_n = idc.len();
        let ids = lookup_ids(&index, Which::NodeSpells, idn, idc, &id_name, &mut report);
        let mut ok = true;
        let mut bound = |role: SpellRole, fill: f64| -> Vec<f64> {
            let (name, col) = frame.spell_role(role).expect("checked");
            match numeric(col) {
                Some(v) => v
                    .into_iter()
                    .map(|x| if x.is_nan() { fill } else { x })
                    .collect(),
                None => {
                    ok = false;
                    report.push(Diagnostic::network(
                        Code::E_TIME_TYPE,
                        format!(
                            "spell {} times (column `{name}` of `node_spells`) must be numeric, not {}",
                            role.name(),
                            col.kind().name()
                        ),
                    ));
                    vec![fill; rows_n]
                }
            }
        };
        let a = bound(SpellRole::Start, f64::NEG_INFINITY);
        let b = bound(SpellRole::End, f64::INFINITY);
        if ok
            && report_rows(
                &mut report,
                Code::E_NODE_SPELL_ORDER,
                Which::NodeSpells,
                "nodes must not leave before they arrive (spell end earlier than spell start)",
                rows_n,
                threads,
                |i| b[i] < a[i],
            )
        {
            ok = false;
        }
        if let (true, Some(ids)) = (ok, ids) {
            // Counting sort by node, then sort each node's spells by
            // start and end.
            let mut offsets = vec![0usize; n + 1];
            for &i in &ids {
                offsets[i as usize + 1] += 1;
            }
            for i in 0..n {
                offsets[i + 1] += offsets[i];
            }
            let mut rows = vec![0usize; rows_n];
            let mut fill = offsets.clone();
            for (r, &i) in ids.iter().enumerate() {
                rows[fill[i as usize]] = r;
                fill[i as usize] += 1;
            }
            for i in 0..n {
                rows[offsets[i]..offsets[i + 1]].sort_by(|&x, &y| {
                    (a[x], b[x], x)
                        .partial_cmp(&(a[y], b[y], y))
                        .expect("no NaN")
                });
            }
            let start: Vec<f64> = rows.iter().map(|&r| a[r]).collect();
            let end: Vec<f64> = rows.iter().map(|&r| b[r]).collect();
            if decl.overlapping_spells != Some(true) {
                // A spell overlaps when it starts before an earlier spell
                // of its node ends.
                let mut overlapping = vec![false; rows_n];
                for i in 0..n {
                    let mut latest_end = f64::NEG_INFINITY;
                    for k in offsets[i]..offsets[i + 1] {
                        if start[k] < latest_end {
                            overlapping[rows[k]] = true;
                        }
                        latest_end = latest_end.max(end[k]);
                    }
                }
                if report_rows(
                    &mut report,
                    Code::E_NODE_SPELL_OVERLAP,
                    Which::NodeSpells,
                    "the presence spells of one node must not overlap (declare overlapping_spells to allow it); the listed spells start before an earlier spell of the same node ends",
                    rows_n,
                    threads,
                    |r| overlapping[r],
                ) {
                    ok = false;
                }
            }
            if ok {
                spell_table = Some(NodeSpells {
                    offsets,
                    start,
                    end,
                    rows,
                });
            }
        }
    }

    // The node changes table: identifiers looked up, times finite.
    if change_table_ok && frame.node_changes().is_some() {
        let (idn, idc) = frame.change_role(ChangeRole::Id).expect("checked");
        lookup_ids(&index, Which::NodeChanges, idn, idc, &id_name, &mut report);
        let (tn, tc) = frame.change_role(ChangeRole::Time).expect("checked");
        match numeric(tc) {
            None => report.push(Diagnostic::network(
                Code::E_TIME_TYPE,
                format!(
                    "change times (column `{tn}` of `node_changes`) must be numeric, not {}",
                    tc.kind().name()
                ),
            )),
            Some(t) => {
                report_rows(
                    &mut report,
                    Code::E_TIME_MISSING,
                    Which::NodeChanges,
                    format!("change times (column `{tn}` of `node_changes`) must not be missing"),
                    t.len(),
                    threads,
                    |i| t[i].is_nan(),
                );
                report_rows(
                    &mut report,
                    Code::E_TIME_NOT_FINITE,
                    Which::NodeChanges,
                    format!("change times (column `{tn}` of `node_changes`) must be finite"),
                    t.len(),
                    threads,
                    |i| t[i].is_infinite(),
                );
            }
        }
    }

    // Edge spells within node spells, from the roles or the table.
    if let Some(s) = &start {
        let e_ = end.as_deref();
        let message = "an edge can only be active while both of its endpoints are present; the listed edges fall outside a node's presence spell";
        if let Some(ns) = &node_spells {
            let inside = |node: u32, s: f64, e: f64| {
                let (a, b) = ns[node as usize];
                a <= s && e <= b && (e > s || s < b || a == b)
            };
            report_rows(
                &mut report,
                Code::E_SPELL_OUTSIDE_NODE,
                Which::Edges,
                message,
                m,
                threads,
                |i| {
                    let e = e_.map_or(s[i], |e| e[i]);
                    !(inside(src[i], s[i], e) && inside(tgt[i], s[i], e))
                },
            );
        } else if let Some(table) = &spell_table {
            let overlapping = decl.overlapping_spells == Some(true);
            report_rows(
                &mut report,
                Code::E_SPELL_OUTSIDE_NODE,
                Which::Edges,
                message,
                m,
                threads,
                |i| {
                    let e = e_.map_or(s[i], |e| e[i]);
                    !(table.contains(src[i] as usize, s[i], e, overlapping)
                        && table.contains(tgt[i] as usize, s[i], e, overlapping))
                },
            );
        }
    }

    // Parallel edges, keyed on endpoints with their layers, relation, and
    // spell.
    let groups_possible = (frame.edge_role(EdgeRole::Relation).is_none() || relations.is_some())
        && (!(has_start || has_wave) || start.is_some())
        && (!(has_source_layer && has_target_layer) || layers.is_some());
    let mut groups = None;
    if groups_possible {
        let g = edge_groups(
            n,
            &src,
            &tgt,
            decl.directed,
            KeyColumns {
                relation: relations.as_ref().map(|r| &r.codes[..]),
                start: start.as_deref(),
                end: end.as_deref(),
                layers: layer_codes,
            },
        );
        let parallel_rows: Vec<usize> = (0..m)
            .filter(|&e| g.first[g.group[e] as usize] != e)
            .take(ROW_CAP)
            .collect();
        let parallel = m - g.sizes.len();
        let what = {
            let mut parts = vec![format!("columns `{src_name}` and `{tgt_name}`")];
            if layers.is_some() {
                parts.push("layers".to_string());
            }
            if let Some((rn, _)) = frame.edge_role(EdgeRole::Relation) {
                parts.push(format!("relation `{rn}`"));
            }
            if has_start || has_wave {
                parts.push("time".to_string());
            }
            parts.join(", ")
        };
        // A pair cannot be both an observed tie and an unobserved pair.
        if let Some(flags) = &missing {
            let k = g.sizes.len();
            let mut observed = vec![false; k];
            let mut unobserved = vec![false; k];
            for (&grp, &flag) in g.group.iter().zip(flags) {
                if flag {
                    unobserved[grp as usize] = true;
                } else {
                    observed[grp as usize] = true;
                }
            }
            report_rows(
                &mut report,
                Code::E_MISSING_CONFLICT,
                Which::Edges,
                format!(
                    "a pair cannot be both an observed tie and an unobserved pair; the listed rows are marked unobserved (column `{}` of `edges`) but another row with the same {what} is an observed tie",
                    frame.edge_role(EdgeRole::EdgeMissing).map(|(n, _)| n).unwrap_or_default()
                ),
                m,
                threads,
                |e| {
                    let grp = g.group[e] as usize;
                    flags[e] && observed[grp] && unobserved[grp]
                },
            );
        }
        match decl.multiple {
            Some(false) if parallel > 0 => report.push(Diagnostic::capped(
                Code::E_PARALLEL_DECLARED,
                Which::Edges,
                format!(
                    "parallel edges were declared absent (multiple = FALSE), but the listed edges repeat an earlier edge with the same {what}"
                ),
                parallel_rows,
                parallel,
            )),
            None if parallel > 0 => report.push(Diagnostic::capped(
                Code::N_PARALLEL,
                Which::Edges,
                format!(
                    "the network has {parallel} parallel edge(s) (same {what}); each measure's documentation says how they count"
                ),
                parallel_rows,
                parallel,
            )),
            _ => {}
        }
        groups = Some(g);
    }

    // 4. Notes: isolates and 0-based identifiers.
    let mut has_edge = vec![false; n];
    for e in 0..m {
        has_edge[src[e] as usize] = true;
        has_edge[tgt[e] as usize] = true;
    }
    report_rows(
        &mut report,
        Code::N_ISOLATES,
        Which::Nodes,
        format!(
            "{} node(s) have no edges",
            has_edge.iter().filter(|&&x| !x).count()
        ),
        n0,
        threads,
        |i| !has_edge[i],
    );
    if id_kind_resolved == Kind::Int && n > 1 {
        let (mut lo, mut hi) = (i64::MAX, i64::MIN);
        let mut visit = |k: Option<i64>| {
            if let Some(k) = k {
                lo = lo.min(k);
                hi = hi.max(k);
            }
        };
        if let Some((_, id)) = frame.node_role(NodeRole::Id) {
            let v = IdView::of(id).expect("checked");
            for i in 0..n0 {
                visit(v.int(i));
            }
        }
        let sv =
            IdView::of(frame.edge_role(EdgeRole::Source).expect("checked").1).expect("checked");
        let tv =
            IdView::of(frame.edge_role(EdgeRole::Target).expect("checked").1).expect("checked");
        for &(e, is_target) in &added {
            visit(if is_target { tv.int(e) } else { sv.int(e) });
        }
        if lo == 0 && hi == (n as i64) - 1 {
            report.push(Diagnostic::network(
                Code::N_ZERO_BASED_IDS,
                format!(
                    "the integer identifiers run from 0 to {}; if they are positions in a 1-based table, add 1 to them",
                    n - 1
                ),
            ));
        }
    }

    Validation {
        report,
        resolved: Some(Resolved {
            n,
            src,
            tgt,
            added,
            id_kind: id_kind_resolved,
            modes,
            relations,
            weights,
            start,
            end,
            node_spells,
            spell_table,
            layers,
            egos,
            missing,
            groups,
            loops,
            directed: decl.directed,
        }),
    }
}

/// The resolved columns that distinguish edges with the same endpoints.
#[derive(Clone, Copy, Default)]
pub(super) struct KeyColumns<'a> {
    pub relation: Option<&'a [u32]>,
    pub start: Option<&'a [f64]>,
    pub end: Option<&'a [f64]>,
    pub layers: Option<(&'a [u32], &'a [u32])>,
}

impl<'a> KeyColumns<'a> {
    /// The key columns of a resolved frame.
    pub fn of(resolved: &'a Resolved) -> KeyColumns<'a> {
        KeyColumns {
            relation: resolved.relations.as_ref().map(|r| &r.codes[..]),
            start: resolved.start.as_deref(),
            end: resolved.end.as_deref(),
            layers: resolved.layers.as_ref().map(|l| (&l.src[..], &l.tgt[..])),
        }
    }
}

/// Groups edges by their key: endpoints with their layers (ordered when
/// directed), relation code, and spell, as the module documentation
/// defines parallel edges. Groups are numbered in order of first
/// appearance.
///
/// Edges are ordered by endpoint pair with two stable counting sorts (the
/// same O(n + m) passes that build the CSR rows), then each run of equal
/// pairs is split by layers, relation, and spell. This touches memory
/// sequentially, where a hash table keyed on the whole tuple did not.
pub(super) fn edge_groups(
    n: usize,
    src: &[u32],
    tgt: &[u32],
    directed: bool,
    key: KeyColumns<'_>,
) -> Groups {
    let KeyColumns {
        relation,
        start,
        end,
        layers,
    } = key;
    let m = src.len();
    let mut a = Vec::with_capacity(m);
    let mut b = Vec::with_capacity(m);
    // The layers of the oriented endpoints: in an undirected network the
    // node-layer pairs are ordered, so (a, la) <= (b, lb).
    let mut la = Vec::with_capacity(if layers.is_some() { m } else { 0 });
    let mut lb = Vec::with_capacity(la.capacity());
    for e in 0..m {
        let (ls, lt) = layers.map_or((0, 0), |(s, t)| (s[e], t[e]));
        let swap = !directed && (src[e] > tgt[e] || (src[e] == tgt[e] && ls > lt));
        if swap {
            a.push(tgt[e]);
            b.push(src[e]);
        } else {
            a.push(src[e]);
            b.push(tgt[e]);
        }
        if layers.is_some() {
            la.push(if swap { lt } else { ls });
            lb.push(if swap { ls } else { lt });
        }
    }
    // Stable counting sort by b, then by a: equal pairs end up adjacent,
    // in row order.
    let mut order = vec![0usize; m];
    let mut count = vec![0usize; n + 1];
    for &x in &b {
        count[x as usize + 1] += 1;
    }
    for i in 0..n {
        count[i + 1] += count[i];
    }
    for e in 0..m {
        order[count[b[e] as usize]] = e;
        count[b[e] as usize] += 1;
    }
    let mut count = vec![0usize; n + 1];
    for &x in &a {
        count[x as usize + 1] += 1;
    }
    for i in 0..n {
        count[i + 1] += count[i];
    }
    let mut sorted = vec![0usize; m];
    for &e in &order {
        sorted[count[a[e] as usize]] = e;
        count[a[e] as usize] += 1;
    }
    drop(order);

    // Split each run of equal pairs by (layers, relation, start, end).
    let key = |e: usize| {
        let l = if layers.is_some() {
            ((la[e] as u64) << 32) | lb[e] as u64
        } else {
            0
        };
        let r = relation.map_or(0, |r| r[e]);
        let s = start.map_or(0, |s| (s[e] + 0.0).to_bits());
        let t = end.map_or(s, |t| (t[e] + 0.0).to_bits());
        (l, r, s, t)
    };
    let mut tentative = vec![0u32; m];
    let mut first_row: Vec<usize> = Vec::new();
    let mut run: Vec<((u64, u32, u64, u64), usize)> = Vec::new();
    let mut i = 0;
    while i < m {
        let mut j = i + 1;
        let (pa, pb) = (a[sorted[i]], b[sorted[i]]);
        while j < m && a[sorted[j]] == pa && b[sorted[j]] == pb {
            j += 1;
        }
        if j - i == 1 {
            tentative[sorted[i]] = first_row.len() as u32;
            first_row.push(sorted[i]);
        } else {
            run.clear();
            run.extend(sorted[i..j].iter().map(|&e| (key(e), e)));
            run.sort_unstable();
            let mut prev = None;
            for &(k, e) in &run {
                if prev != Some(k) {
                    prev = Some(k);
                    first_row.push(e);
                }
                tentative[e] = first_row.len() as u32 - 1;
            }
        }
        i = j;
    }

    // Renumber by first appearance.
    let mut relabel = vec![u32::MAX; first_row.len()];
    let mut group = Vec::with_capacity(m);
    let mut sizes: Vec<u32> = Vec::new();
    let mut first = Vec::new();
    for &t in &tentative {
        let t = t as usize;
        if relabel[t] == u32::MAX {
            relabel[t] = sizes.len() as u32;
            sizes.push(0);
            first.push(first_row[t]);
        }
        let g = relabel[t];
        sizes[g as usize] += 1;
        group.push(g);
    }
    Groups {
        group,
        sizes,
        first,
    }
}

#[cfg(test)]
mod tests {
    use super::super::columns::Strings;
    use super::*;

    fn ints(v: &[i64]) -> Column {
        Column::Int(v.to_vec(), None)
    }

    fn strs(v: &[&str]) -> Column {
        Column::Str(Strings::from_strs(
            &v.iter().map(|s| Some(*s)).collect::<Vec<_>>(),
        ))
    }

    /// nodes a, b, c; edges a-b, b-c, c-a.
    fn triangle() -> Frame {
        let mut f = Frame::new(false);
        f.add_node_column("id", strs(&["a", "b", "c"])).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column("from", strs(&["a", "b", "c"])).unwrap();
        f.add_edge_column("to", strs(&["b", "c", "a"])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f
    }

    fn codes(v: &Validation) -> Vec<Code> {
        v.report.diagnostics.iter().map(|d| d.code).collect()
    }

    #[test]
    fn valid_frame_has_no_diagnostics() {
        let v = triangle().validate(1);
        assert!(codes(&v).is_empty(), "{:?}", v.report);
        let r = v.resolved.unwrap();
        assert_eq!((r.n, r.src, r.tgt), (3, vec![0, 1, 2], vec![1, 2, 0]));
        assert_eq!(r.id_kind, Kind::Str);
        assert!(r.groups.unwrap().sizes.iter().all(|&s| s == 1));
    }

    #[test]
    fn missing_roles() {
        let mut f = Frame::new(true);
        f.add_edge_column("from", ints(&[1])).unwrap();
        f.add_node_column("id", ints(&[1])).unwrap();
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![Code::E_ENDPOINTS_UNDECLARED, Code::E_ID_UNDECLARED]
        );
        assert!(v.resolved.is_none());
    }

    #[test]
    fn endpoint_and_id_types() {
        let mut f = triangle();
        f.add_edge_column("num", ints(&[1, 2, 3])).unwrap();
        f.set_edge_role(EdgeRole::Target, "num").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_ENDPOINT_TYPE]);
        let mut f = triangle();
        f.add_edge_column("flag", Column::Bool(vec![true; 3], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Target, "flag").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_ENDPOINT_TYPE]);
        let mut f = triangle();
        f.add_node_column("num", ints(&[1, 2, 3])).unwrap();
        f.set_node_role(NodeRole::Id, "num").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_ID_TYPE]);
        let mut f = triangle();
        f.add_node_column("flag", Column::Bool(vec![true; 3], None))
            .unwrap();
        f.set_node_role(NodeRole::Id, "flag").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_ID_TYPE]);
    }

    #[test]
    fn ids_coerced_between_integer_and_double() {
        let mut f = Frame::new(true);
        f.add_node_column("id", Column::Float(vec![1.0, 2.0, 3.0]))
            .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column("from", ints(&[1, 2])).unwrap();
        f.add_edge_column("to", ints(&[2, 3])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::N_IDS_COERCED]);
        let r = v.resolved.unwrap();
        assert_eq!(
            (r.src, r.tgt, r.id_kind),
            (vec![0, 1], vec![1, 2], Kind::Float)
        );
        // Doubles with fractions are legal identifiers in their own right.
        let mut f = Frame::new(true);
        f.add_edge_column("from", Column::Float(vec![0.5, -0.0]))
            .unwrap();
        f.add_edge_column("to", Column::Float(vec![0.0, 0.5]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty(), "{:?}", v.report);
        assert_eq!(v.resolved.unwrap().n, 2);
    }

    #[test]
    fn missing_empty_duplicate_unknown() {
        let mut f = Frame::new(false);
        f.add_node_column(
            "id",
            Column::Str(Strings::from_strs(&[Some("a"), None, Some(""), Some("a")])),
        )
        .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column(
            "from",
            Column::Str(Strings::from_strs(&[Some("a"), None, Some("z")])),
        )
        .unwrap();
        f.add_edge_column(
            "to",
            Column::Str(Strings::from_strs(&[Some("a"), Some(""), Some("a")])),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        let r = &v.report;
        assert_eq!(r.find(Code::E_MISSING_ID).unwrap().rows, vec![1]);
        assert_eq!(r.find(Code::E_EMPTY_ID).unwrap().rows, vec![2]);
        assert_eq!(r.find(Code::E_DUPLICATE_ID).unwrap().rows, vec![3]);
        assert_eq!(r.find(Code::E_MISSING_ENDPOINT).unwrap().rows, vec![1]);
        assert_eq!(r.find(Code::E_UNKNOWN_NODE).unwrap().rows, vec![2]);
        let empties: Vec<&Diagnostic> = r
            .diagnostics
            .iter()
            .filter(|d| d.code == Code::E_EMPTY_ID)
            .collect();
        assert_eq!(empties[1].table, Some(Which::Edges));
        assert_eq!(empties[1].rows, vec![1]);
        assert!(v.resolved.is_none());
        assert!(!r.is_valid());
    }

    #[test]
    fn nodes_added_from_edges() {
        let mut f = triangle();
        f.add_edge_column("to2", strs(&["b", "d", "d"])).unwrap();
        f.set_edge_role(EdgeRole::Target, "to2").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_UNKNOWN_NODE]);
        f.declarations_mut().add_nodes = true;
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::N_NODES_ADDED]);
        assert_eq!(v.report.find(Code::N_NODES_ADDED).unwrap().rows, vec![1]);
        let r = v.resolved.unwrap();
        assert_eq!((r.n, r.added), (4, vec![(1, true)]));
        assert_eq!(r.tgt, vec![1, 3, 3]);
        // Without a node table every node comes from the edges and nothing
        // is noted.
        let mut f = Frame::new(true);
        f.add_edge_column("from", strs(&["x", "y"])).unwrap();
        f.add_edge_column("to", strs(&["y", "z"])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty());
        let r = v.resolved.unwrap();
        assert_eq!(r.added, vec![(0, false), (0, true), (1, true)]);
        assert_eq!((r.src, r.tgt), (vec![0, 1], vec![1, 2]));
    }

    #[test]
    fn loops_and_parallel_declared_or_noted() {
        let mut f = Frame::new(false);
        f.add_edge_column("from", ints(&[1, 1, 2, 2])).unwrap();
        f.add_edge_column("to", ints(&[1, 2, 1, 3])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::N_LOOPS, Code::N_PARALLEL]);
        assert_eq!(v.report.find(Code::N_LOOPS).unwrap().rows, vec![0]);
        assert_eq!(v.report.find(Code::N_PARALLEL).unwrap().rows, vec![2]);
        f.declarations_mut().loops = Some(false);
        f.declarations_mut().multiple = Some(false);
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![Code::E_LOOP_DECLARED, Code::E_PARALLEL_DECLARED]
        );
        f.declarations_mut().loops = Some(true);
        f.declarations_mut().multiple = Some(true);
        assert!(codes(&f.validate(1)).is_empty());
        // Directed: 2 -> 1 is not parallel to 1 -> 2.
        f.declarations_mut().directed = true;
        f.declarations_mut().multiple = Some(false);
        assert!(codes(&f.validate(1)).is_empty());
    }

    #[test]
    fn parallel_edges_are_relative_to_relation_and_time() {
        let mut f = Frame::new(true);
        f.add_edge_column("from", ints(&[1, 1, 1])).unwrap();
        f.add_edge_column("to", ints(&[2, 2, 2])).unwrap();
        f.add_edge_column("rel", strs(&["friend", "advice", "friend"]))
            .unwrap();
        f.add_edge_column("t", Column::Float(vec![1.0, 1.0, 2.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::N_PARALLEL).unwrap().count, 2);
        f.set_edge_role(EdgeRole::Relation, "rel").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::N_PARALLEL).unwrap().rows, vec![2]);
        assert_eq!(
            v.resolved
                .as_ref()
                .unwrap()
                .relations
                .as_ref()
                .unwrap()
                .labels,
            vec!["friend", "advice"]
        );
        f.set_edge_role(EdgeRole::Start, "t").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty(), "{:?}", v.report);
        let g = v.resolved.unwrap().groups.unwrap();
        assert_eq!(g.sizes, vec![1, 1, 1]);
    }

    #[test]
    fn modes() {
        let mut f = triangle();
        f.add_node_column("mode", strs(&["person", "event", "person"]))
            .unwrap();
        f.set_node_role(NodeRole::Mode, "mode").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty());
        assert_eq!(
            v.resolved.unwrap().modes.unwrap().labels,
            vec!["person", "event"]
        );
        f.declarations_mut().between_modes_only = Some(true);
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_MODE_WITHIN]);
        assert_eq!(v.report.find(Code::E_MODE_WITHIN).unwrap().rows, vec![2]);
        f.declarations_mut().loops = Some(true);
        assert!(f.validate(1).report.has(Code::E_DECLARATION_CONFLICT));
        f.declarations_mut().loops = None;
        f.clear_node_role(NodeRole::Mode);
        assert_eq!(codes(&f.validate(1)), vec![Code::E_MODE_UNDECLARED]);
        let mut f = triangle();
        f.add_node_column("mode", Column::Float(vec![1.0; 3]))
            .unwrap();
        f.set_node_role(NodeRole::Mode, "mode").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_MODE_TYPE]);
        let mut f = triangle();
        f.add_node_column(
            "mode",
            Column::Bool(vec![true; 3], Some(vec![false, true, false])),
        )
        .unwrap();
        f.set_node_role(NodeRole::Mode, "mode").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_MODE_MISSING]);
        assert_eq!(v.report.find(Code::E_MODE_MISSING).unwrap().rows, vec![1]);
    }

    #[test]
    fn relations_and_signs() {
        let mut f = triangle();
        f.add_edge_column("rel", Column::Float(vec![1.0; 3]))
            .unwrap();
        f.set_edge_role(EdgeRole::Relation, "rel").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_RELATION_TYPE]);
        let mut f = triangle();
        f.add_edge_column(
            "rel",
            Column::Int(vec![1; 3], Some(vec![false, false, true])),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::Relation, "rel").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_RELATION_MISSING]);
        let mut f = triangle();
        f.add_edge_column("sign", ints(&[1, -1, 0])).unwrap();
        f.set_edge_role(EdgeRole::Sign, "sign").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_SIGN_VALUE]);
        assert_eq!(v.report.find(Code::E_SIGN_VALUE).unwrap().rows, vec![2]);
        let mut f = triangle();
        f.add_edge_column(
            "sign",
            Column::Bool(vec![true; 3], Some(vec![true, false, false])),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::Sign, "sign").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_SIGN_VALUE]);
        let mut f = triangle();
        f.add_edge_column("sign", strs(&["+", "-", "+"])).unwrap();
        f.set_edge_role(EdgeRole::Sign, "sign").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_SIGN_TYPE]);
    }

    #[test]
    fn weights() {
        let mut f = triangle();
        f.add_edge_column("w", strs(&["1", "2", "3"])).unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_WEIGHT_TYPE]);
        let mut f = triangle();
        f.add_edge_column("w", Column::Float(vec![1.0, f64::NAN, f64::INFINITY]))
            .unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![Code::E_WEIGHT_TYPE_MISSING, Code::E_WEIGHT_NOT_FINITE]
        );
        assert_eq!(
            v.report.find(Code::E_WEIGHT_NOT_FINITE).unwrap().rows,
            vec![1, 2]
        );
        let mut f = triangle();
        f.add_edge_column("w", Column::Int(vec![1, 0, -1], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Distance));
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_DISTANCE_NOT_POSITIVE]);
        assert_eq!(
            v.report.find(Code::E_DISTANCE_NOT_POSITIVE).unwrap().rows,
            vec![1, 2]
        );
        f.set_weight_type(Some(WeightType::Strength));
        let v = f.validate(1);
        assert!(codes(&v).is_empty());
        assert_eq!(v.resolved.unwrap().weights.unwrap(), vec![1.0, 0.0, -1.0]);
        f.add_edge_column("sign", ints(&[1, 1, -1])).unwrap();
        f.set_edge_role(EdgeRole::Sign, "sign").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_WEIGHT_NEGATIVE]);
        assert_eq!(
            v.report.find(Code::E_WEIGHT_NEGATIVE).unwrap().rows,
            vec![2]
        );
    }

    #[test]
    fn times() {
        let mut f = triangle();
        f.add_edge_column("end", Column::Float(vec![1.0; 3]))
            .unwrap();
        f.set_edge_role(EdgeRole::End, "end").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_DECLARATION_CONFLICT]);
        f.add_edge_column("start", strs(&["a", "b", "c"])).unwrap();
        f.set_edge_role(EdgeRole::Start, "start").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_TIME_TYPE]);
        let mut f = triangle();
        f.add_edge_column(
            "start",
            Column::Float(vec![1.0, f64::NAN, f64::NEG_INFINITY]),
        )
        .unwrap();
        f.add_edge_column("end", Column::Float(vec![0.5, 2.0, 3.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Start, "start").unwrap();
        f.set_edge_role(EdgeRole::End, "end").unwrap();
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![
                Code::E_TIME_MISSING,
                Code::E_TIME_NOT_FINITE,
                Code::E_SPELL_ORDER
            ]
        );
        assert_eq!(v.report.find(Code::E_SPELL_ORDER).unwrap().rows, vec![0]);
        assert!(v.resolved.unwrap().start.is_none());
        // Contacts and waves.
        let mut f = triangle();
        f.add_edge_column("t", ints(&[1, 2, 3])).unwrap();
        f.set_edge_role(EdgeRole::Start, "t").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty());
        let r = v.resolved.unwrap();
        assert_eq!(r.start.unwrap(), vec![1.0, 2.0, 3.0]);
        assert!(r.end.is_none());
        let mut f = triangle();
        f.add_edge_column("wave", Column::Float(vec![1.0, 1.0, 2.5]))
            .unwrap();
        f.set_edge_role(EdgeRole::Wave, "wave").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_WAVE_TYPE]);
        let mut f = triangle();
        f.add_edge_column(
            "wave",
            Column::Int(vec![1, 1, 2], Some(vec![false, false, true])),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::Wave, "wave").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_WAVE_MISSING]);
        let mut f = triangle();
        f.add_edge_column("wave", ints(&[1, 1, 2])).unwrap();
        f.add_edge_column("t", ints(&[1, 2, 3])).unwrap();
        f.set_edge_role(EdgeRole::Wave, "wave").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty());
        let r = v.resolved.unwrap();
        assert_eq!(r.end.unwrap(), vec![2.0, 2.0, 3.0]);
        f.set_edge_role(EdgeRole::Start, "t").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_DECLARATION_CONFLICT]);
    }

    #[test]
    fn node_spells() {
        let mut f = triangle();
        f.add_node_column("arrive", Column::Float(vec![0.0, 5.0, f64::NAN]))
            .unwrap();
        f.add_node_column("leave", Column::Float(vec![10.0, 4.0, 8.0]))
            .unwrap();
        f.set_node_role(NodeRole::NodeStart, "arrive").unwrap();
        f.set_node_role(NodeRole::NodeEnd, "leave").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_NODE_SPELL_ORDER]);
        assert_eq!(
            v.report.find(Code::E_NODE_SPELL_ORDER).unwrap().rows,
            vec![1]
        );
        let mut f = triangle();
        f.add_node_column("arrive", Column::Float(vec![0.0, 5.0, f64::NAN]))
            .unwrap();
        f.add_node_column("leave", Column::Float(vec![10.0, 9.0, 8.0]))
            .unwrap();
        f.set_node_role(NodeRole::NodeStart, "arrive").unwrap();
        f.set_node_role(NodeRole::NodeEnd, "leave").unwrap();
        // a-b at [1, 2) starts before b arrives; b-c at 6 is fine; c-a at
        // [8, 9) runs past c's leaving at 8.
        f.add_edge_column("s", Column::Float(vec![1.0, 6.0, 8.0]))
            .unwrap();
        f.add_edge_column("e", Column::Float(vec![2.0, 6.0, 9.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Start, "s").unwrap();
        f.set_edge_role(EdgeRole::End, "e").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::E_SPELL_OUTSIDE_NODE]);
        assert_eq!(
            v.report.find(Code::E_SPELL_OUTSIDE_NODE).unwrap().rows,
            vec![0, 2]
        );
        let r = v.resolved.unwrap();
        assert_eq!(r.node_spells.unwrap()[2], (f64::NEG_INFINITY, 8.0));
        let mut f = triangle();
        f.add_node_column("arrive", strs(&["x", "y", "z"])).unwrap();
        f.set_node_role(NodeRole::NodeStart, "arrive").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_TIME_TYPE]);
    }

    #[test]
    fn coordinates_and_crs() {
        let mut f = triangle();
        f.add_node_column("x", Column::Float(vec![0.0, 1.0, f64::NAN]))
            .unwrap();
        f.set_node_role(NodeRole::X, "x").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_COORD_ROLE_PAIR]);
        f.add_node_column("y", Column::Float(vec![0.0, 91.0, 2.0]))
            .unwrap();
        f.set_node_role(NodeRole::Y, "y").unwrap();
        f.geometry_mut().coordinates = CoordinateKind::LonLat;
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![Code::E_COORD_NOT_FINITE, Code::E_LATITUDE_RANGE]
        );
        assert_eq!(v.report.find(Code::E_LATITUDE_RANGE).unwrap().rows, vec![1]);
        f.geometry_mut().coordinates = CoordinateKind::Planar;
        assert_eq!(codes(&f.validate(1)), vec![Code::E_COORD_NOT_FINITE]);
        let mut f = triangle();
        f.add_node_column("x", strs(&["0", "1", "2"])).unwrap();
        f.add_node_column("y", Column::Float(vec![0.0; 3])).unwrap();
        f.set_node_role(NodeRole::X, "x").unwrap();
        f.set_node_role(NodeRole::Y, "y").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_COORD_TYPE]);
        let mut f = triangle();
        let g = f.geometry_mut();
        g.nodes = true;
        g.edges = true;
        g.node_crs = Some("EPSG:4326".into());
        g.edge_crs = Some("EPSG:3857".into());
        assert_eq!(codes(&f.validate(1)), vec![Code::E_CRS_MISMATCH]);
        f.geometry_mut().edge_crs = Some("EPSG:4326".into());
        assert!(codes(&f.validate(1)).is_empty());
    }

    #[test]
    fn missing_data_roles() {
        let mut f = triangle();
        f.add_node_column("na", Column::Bool(vec![false, true, false], None))
            .unwrap();
        f.set_node_role(NodeRole::NodeMissing, "na").unwrap();
        f.add_edge_column("na", Column::Bool(vec![true, true, false], None))
            .unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![Code::N_MISSING_NODES, Code::N_MISSING_EDGES]
        );
        assert_eq!(v.report.find(Code::N_MISSING_EDGES).unwrap().count, 2);
        let mut f = triangle();
        f.add_node_column("na", ints(&[0, 1, 0])).unwrap();
        f.set_node_role(NodeRole::NodeMissing, "na").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_MISSING_ROLE_TYPE]);
        let mut f = triangle();
        f.add_edge_column(
            "na",
            Column::Bool(vec![false; 3], Some(vec![true, false, false])),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        assert_eq!(codes(&f.validate(1)), vec![Code::E_MISSING_ROLE_TYPE]);
    }

    /// Nodes a, b; edges a-b in layer 1, a-b in layer 2, a-a across
    /// layers 1 and 2 (a coupling), b-b in layer 2 (a loop), a-b from
    /// layer 1 to layer 2 (interlayer).
    fn multilayer() -> Frame {
        let mut f = Frame::new(false);
        f.add_node_column("id", strs(&["a", "b"])).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column("from", strs(&["a", "a", "a", "b", "a"]))
            .unwrap();
        f.add_edge_column("to", strs(&["b", "b", "a", "b", "b"]))
            .unwrap();
        f.add_edge_column("ls", ints(&[1, 2, 1, 2, 1])).unwrap();
        f.add_edge_column("lt", ints(&[1, 2, 2, 2, 2])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "ls").unwrap();
        f.set_edge_role(EdgeRole::TargetLayer, "lt").unwrap();
        f
    }

    #[test]
    fn layers_are_classified_and_keyed() {
        let f = multilayer();
        let v = f.validate(1);
        // The loop b-b is the only loop; the coupling is not one, and
        // the two a-b edges in different layers are not parallel.
        assert_eq!(codes(&v), vec![Code::N_LOOPS]);
        assert_eq!(v.report.find(Code::N_LOOPS).unwrap().rows, vec![3]);
        let r = v.resolved.unwrap();
        let l = r.layers.as_ref().unwrap();
        assert_eq!(l.labels, vec!["1", "2"]);
        assert!(!l.from_table);
        let classes: Vec<EdgeClass> = (0..5).map(|e| l.class(e, &r.src, &r.tgt)).collect();
        assert_eq!(
            classes,
            vec![
                EdgeClass::Intralayer,
                EdgeClass::Intralayer,
                EdgeClass::Coupling,
                EdgeClass::Intralayer,
                EdgeClass::Interlayer
            ]
        );
        assert!(r.is_loop(3) && !r.is_loop(2));
        assert_eq!(r.groups.unwrap().sizes, vec![1, 1, 1, 1, 1]);
        // Undirected: a coupling written the other way round is parallel
        // to the first, and so is b-a in layer 2 to a-b in layer 2.
        let mut f = multilayer();
        f.add_edge_column("from2", strs(&["a", "b", "a", "b", "a"]))
            .unwrap();
        f.add_edge_column("lt2", ints(&[1, 2, 1, 2, 2])).unwrap();
        f.add_edge_column("ls2", ints(&[1, 2, 2, 2, 1])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from2").unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "ls2").unwrap();
        f.set_edge_role(EdgeRole::TargetLayer, "lt2").unwrap();
        let mut g = Frame::new(false);
        g.add_edge_column("from", strs(&["a", "b", "a", "a"]))
            .unwrap();
        g.add_edge_column("to", strs(&["b", "a", "a", "a"]))
            .unwrap();
        g.add_edge_column("ls", ints(&[2, 2, 1, 2])).unwrap();
        g.add_edge_column("lt", ints(&[2, 2, 2, 1])).unwrap();
        g.set_edge_role(EdgeRole::Source, "from").unwrap();
        g.set_edge_role(EdgeRole::Target, "to").unwrap();
        g.set_edge_role(EdgeRole::SourceLayer, "ls").unwrap();
        g.set_edge_role(EdgeRole::TargetLayer, "lt").unwrap();
        let v = g.validate(1);
        assert_eq!(v.report.find(Code::N_PARALLEL).unwrap().rows, vec![1, 3]);
        g.declarations_mut().directed = true;
        let v = g.validate(1);
        assert!(!v.report.has(Code::N_PARALLEL), "{:?}", v.report);
        // Declarations.
        let mut f = multilayer();
        f.declarations_mut().couplings_only = Some(true);
        let v = f.validate(1);
        assert_eq!(
            v.report.find(Code::E_INTERLAYER_DECLARED).unwrap().rows,
            vec![4]
        );
        f.declarations_mut().loops = Some(false);
        assert_eq!(v.report.find(Code::N_LOOPS).unwrap().rows, vec![3]);
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_LOOP_DECLARED).unwrap().rows, vec![3]);
        // Modes: the coupling is exempt from the between-modes rule.
        let mut f = multilayer();
        f.add_node_column("mode", strs(&["p", "e"])).unwrap();
        f.set_node_role(NodeRole::Mode, "mode").unwrap();
        f.declarations_mut().between_modes_only = Some(true);
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_MODE_WITHIN).unwrap().rows, vec![3]);
    }

    #[test]
    fn layer_roles_types_and_table() {
        let mut f = multilayer();
        f.clear_edge_role(EdgeRole::TargetLayer);
        let v = f.validate(1);
        assert!(v.report.has(Code::E_LAYER_ROLE_PAIR));
        assert!(v.resolved.unwrap().layers.is_none());
        let mut f = multilayer();
        f.declarations_mut().couplings_only = Some(false);
        f.clear_edge_role(EdgeRole::SourceLayer);
        let v = f.validate(1);
        assert!(v.report.has(Code::E_DECLARATION_CONFLICT));
        let mut f = multilayer();
        f.add_edge_column("lf", Column::Float(vec![1.0; 5]))
            .unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "lf").unwrap();
        assert!(f.validate(1).report.has(Code::E_LAYER_TYPE));
        let mut f = multilayer();
        f.add_edge_column("lsx", strs(&["1", "2", "1", "2", "1"]))
            .unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "lsx").unwrap();
        assert!(f.validate(1).report.has(Code::E_LAYER_TYPE));
        let mut f = multilayer();
        f.add_edge_column(
            "lsm",
            Column::Int(
                vec![1, 2, 1, 2, 1],
                Some(vec![false, true, false, false, false]),
            ),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "lsm").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_LAYER_MISSING).unwrap().rows, vec![1]);
        // A layer table fixes the labels and their order; unknown layers
        // are errors.
        let mut f = multilayer();
        f.add_column(Which::Layers, "layer", ints(&[2, 1, 3]))
            .unwrap();
        f.add_column(Which::Layers, "medium", strs(&["x", "y", "z"]))
            .unwrap();
        let v = f.validate(1);
        assert!(v.report.has(Code::E_TABLE_ROLE_UNDECLARED));
        f.set_layer_role(super::super::roles::LayerRole::Layer, "layer")
            .unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::N_LOOPS]);
        let l = v.resolved.unwrap().layers.unwrap();
        assert_eq!(l.labels, vec!["2", "1", "3"]);
        assert!(l.from_table);
        assert_eq!(l.src, vec![1, 0, 1, 0, 1]);
        assert_eq!(l.tgt, vec![1, 0, 0, 0, 0]);
        f.add_column(Which::Layers, "layer2", ints(&[2, 1, 1]))
            .unwrap();
        f.set_layer_role(super::super::roles::LayerRole::Layer, "layer2")
            .unwrap();
        let v = f.validate(1);
        assert_eq!(
            v.report.find(Code::E_LAYER_ID_DUPLICATE).unwrap().rows,
            vec![2]
        );
        f.add_column(
            Which::Layers,
            "layer3",
            Column::Int(vec![1, 0, 0], Some(vec![false, true, true])),
        )
        .unwrap();
        f.set_layer_role(super::super::roles::LayerRole::Layer, "layer3")
            .unwrap();
        let v = f.validate(1);
        assert_eq!(
            v.report.find(Code::E_LAYER_ID_MISSING).unwrap().rows,
            vec![1, 2]
        );
        assert_eq!(
            v.report.find(Code::E_LAYER_UNKNOWN).unwrap().rows,
            vec![1, 2, 3, 4]
        );
        f.add_column(Which::Layers, "layer4", strs(&["1", "2", "3"]))
            .unwrap();
        f.set_layer_role(super::super::roles::LayerRole::Layer, "layer4")
            .unwrap();
        assert!(f.validate(1).report.has(Code::E_LAYER_TYPE));
    }

    /// Egos ann and bob as nodes, each with alters; edges within egos.
    fn ego_frame() -> Frame {
        let mut f = Frame::new(false);
        f.add_node_column("id", strs(&["ann", "a1", "a2", "bob", "b1"]))
            .unwrap();
        f.add_node_column("ego", strs(&["ann", "ann", "ann", "bob", "bob"]))
            .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.set_node_role(NodeRole::Ego, "ego").unwrap();
        f.add_edge_column("from", strs(&["ann", "ann", "a1", "bob"]))
            .unwrap();
        f.add_edge_column("to", strs(&["a1", "a2", "a2", "b1"]))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f
    }

    #[test]
    fn egos() {
        let mut f = ego_frame();
        f.declarations_mut().egos_are_nodes = Some(true);
        let v = f.validate(1);
        assert!(codes(&v).is_empty(), "{:?}", v.report);
        let eg = v.resolved.unwrap().egos.unwrap();
        assert_eq!(eg.labels, vec!["ann", "bob"]);
        assert_eq!(eg.node, vec![0, 0, 0, 1, 1]);
        assert_eq!(eg.edge, vec![0, 0, 0, 1]);
        assert_eq!(eg.node_of, vec![Some(0), Some(3)]);
        // An explicit edge ego column must agree.
        f.add_edge_column("ego", strs(&["ann", "ann", "bob", "bob"]))
            .unwrap();
        f.set_edge_role(EdgeRole::Ego, "ego").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_EGO_MISMATCH).unwrap().rows, vec![2]);
        // An edge across egos.
        let mut f = ego_frame();
        f.add_edge_column("to2", strs(&["a1", "a2", "b1", "b1"]))
            .unwrap();
        f.set_edge_role(EdgeRole::Target, "to2").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_EGO_MISMATCH).unwrap().rows, vec![2]);
        // An ego that belongs to another ego's network.
        let mut f = ego_frame();
        f.add_node_column("ego2", strs(&["ann", "ann", "ann", "ann", "bob"]))
            .unwrap();
        f.set_node_role(NodeRole::Ego, "ego2").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_EGO_SELF).unwrap().rows, vec![3]);
        // Egos that are not nodes: fine unless declared.
        let mut f = ego_frame();
        f.add_node_column("ego3", strs(&["e1", "e1", "e1", "e2", "e2"]))
            .unwrap();
        f.set_node_role(NodeRole::Ego, "ego3").unwrap();
        let v = f.validate(1);
        assert!(codes(&v).is_empty(), "{:?}", v.report);
        assert_eq!(v.resolved.unwrap().egos.unwrap().node_of, vec![None, None]);
        f.declarations_mut().egos_are_nodes = Some(true);
        let v = f.validate(1);
        assert_eq!(
            v.report.find(Code::E_EGO_UNKNOWN).unwrap().rows,
            vec![0, 1, 2, 3, 4]
        );
        // Types and missing values.
        let mut f = ego_frame();
        f.add_node_column("egof", Column::Float(vec![1.0; 5]))
            .unwrap();
        f.set_node_role(NodeRole::Ego, "egof").unwrap();
        assert!(f.validate(1).report.has(Code::E_EGO_TYPE));
        let mut f = ego_frame();
        f.add_edge_column("egoi", ints(&[1, 1, 1, 2])).unwrap();
        f.set_edge_role(EdgeRole::Ego, "egoi").unwrap();
        assert!(f.validate(1).report.has(Code::E_EGO_TYPE));
        let mut f = ego_frame();
        f.add_node_column(
            "egom",
            Column::Str(Strings::from_strs(&[
                Some("ann"),
                None,
                Some("ann"),
                Some("bob"),
                Some("bob"),
            ])),
        )
        .unwrap();
        f.set_node_role(NodeRole::Ego, "egom").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_EGO_MISSING).unwrap().rows, vec![1]);
        let mut f = ego_frame();
        f.add_edge_column(
            "egom",
            Column::Str(Strings::from_strs(&[
                Some("ann"),
                None,
                Some("ann"),
                Some("bob"),
            ])),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::Ego, "egom").unwrap();
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_EGO_MISSING).unwrap().rows, vec![1]);
        assert_eq!(
            v.report.find(Code::E_EGO_MISSING).unwrap().table,
            Some(Which::Edges)
        );
        // Added nodes have no ego; an edge ego column needs the node one.
        let mut f = ego_frame();
        f.add_edge_column("to3", strs(&["a1", "a2", "a2", "zed"]))
            .unwrap();
        f.set_edge_role(EdgeRole::Target, "to3").unwrap();
        f.declarations_mut().add_nodes = true;
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_EGO_MISSING).unwrap().rows, vec![3]);
        let mut f = ego_frame();
        f.clear_node_role(NodeRole::Ego);
        f.add_edge_column("ego", strs(&["ann", "ann", "ann", "bob"]))
            .unwrap();
        f.set_edge_role(EdgeRole::Ego, "ego").unwrap();
        assert!(f.validate(1).report.has(Code::E_DECLARATION_CONFLICT));
        f.clear_edge_role(EdgeRole::Ego);
        f.declarations_mut().egos_are_nodes = Some(false);
        assert!(f.validate(1).report.has(Code::E_DECLARATION_CONFLICT));
        // Integer egos against double identifiers match as numbers.
        let mut f = Frame::new(true);
        f.add_node_column("id", Column::Float(vec![1.0, 2.0, 3.0]))
            .unwrap();
        f.add_node_column("ego", ints(&[1, 1, 3])).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.set_node_role(NodeRole::Ego, "ego").unwrap();
        f.add_edge_column("from", ints(&[1])).unwrap();
        f.add_edge_column("to", ints(&[2])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert_eq!(
            v.resolved.unwrap().egos.unwrap().node_of,
            vec![Some(0), Some(2)]
        );
    }

    fn with_spell_table(f: &mut Frame, ids: &[&str], s: &[f64], e: &[f64]) {
        use super::super::roles::SpellRole;
        f.add_column(Which::NodeSpells, "id", strs(ids)).unwrap();
        f.add_column(Which::NodeSpells, "s", Column::Float(s.to_vec()))
            .unwrap();
        f.add_column(Which::NodeSpells, "e", Column::Float(e.to_vec()))
            .unwrap();
        f.set_spell_role(SpellRole::Id, "id").unwrap();
        f.set_spell_role(SpellRole::Start, "s").unwrap();
        f.set_spell_role(SpellRole::End, "e").unwrap();
    }

    #[test]
    fn node_spell_table() {
        use super::super::roles::SpellRole;
        let mut f = triangle();
        // a: [8, 20) then [0, 5); b: [0, 20); c: unbounded.
        with_spell_table(
            &mut f,
            &["a", "a", "b"],
            &[8.0, 0.0, 0.0],
            &[20.0, 5.0, 20.0],
        );
        f.add_edge_column("s", Column::Float(vec![1.0, 6.0, 6.0]))
            .unwrap();
        f.add_edge_column("e", Column::Float(vec![10.0, 6.0, 7.0]))
            .unwrap();
        f.set_edge_role(EdgeRole::Start, "s").unwrap();
        f.set_edge_role(EdgeRole::End, "e").unwrap();
        let v = f.validate(1);
        // a-b at [1, 10) spans a's absence; c-a at [6, 7) falls in it.
        assert_eq!(codes(&v), vec![Code::E_SPELL_OUTSIDE_NODE]);
        assert_eq!(
            v.report.find(Code::E_SPELL_OUTSIDE_NODE).unwrap().rows,
            vec![0, 2]
        );
        let t = v.resolved.unwrap().spell_table.unwrap();
        assert_eq!(t.offsets, vec![0, 2, 3, 3]);
        assert_eq!(t.start, vec![0.0, 8.0, 0.0]);
        assert_eq!(t.end, vec![5.0, 20.0, 20.0]);
        assert_eq!(t.rows, vec![1, 0, 2]);
        assert!(t.contains(2, 100.0, 200.0, false));
        assert!(t.contains(0, 4.0, 5.0, false) && !t.contains(0, 4.0, 6.0, false));
        assert!(!t.contains(0, 5.0, 5.0, false));
        assert!(t.contains(0, 8.0, 8.0, false));
        // Overlap unless declared; the spell that starts inside an
        // earlier one is listed.
        let mut f = triangle();
        with_spell_table(&mut f, &["a", "a", "a"], &[0.0, 5.0, 3.0], &[5.0, 9.0, 4.0]);
        let v = f.validate(1);
        assert_eq!(
            v.report.find(Code::E_NODE_SPELL_OVERLAP).unwrap().rows,
            vec![2]
        );
        let mut g = triangle();
        with_spell_table(
            &mut g,
            &["a", "a", "a"],
            &[0.0, 5.0, 3.0],
            &[10.0, 9.0, 4.0],
        );
        let v = g.validate(1);
        assert_eq!(
            v.report.find(Code::E_NODE_SPELL_OVERLAP).unwrap().rows,
            vec![1, 2]
        );
        f.declarations_mut().overlapping_spells = Some(true);
        let v = f.validate(1);
        assert!(codes(&v).is_empty(), "{:?}", v.report);
        let t = v.resolved.unwrap().spell_table.unwrap();
        assert_eq!(t.rows, vec![0, 2, 1]);
        assert!(t.contains(0, 3.5, 4.0, true) && t.contains(0, 6.0, 8.0, true));
        assert!(!t.contains(0, 3.5, 8.0, true));
        // Missing times are unbounded; order and identifiers are checked.
        let mut f = triangle();
        with_spell_table(
            &mut f,
            &["a", "z", "b"],
            &[f64::NAN, 1.0, 5.0],
            &[3.0, 2.0, 4.0],
        );
        let v = f.validate(1);
        assert_eq!(v.report.find(Code::E_UNKNOWN_NODE).unwrap().rows, vec![1]);
        assert_eq!(
            v.report.find(Code::E_UNKNOWN_NODE).unwrap().table,
            Some(Which::NodeSpells)
        );
        assert_eq!(
            v.report.find(Code::E_NODE_SPELL_ORDER).unwrap().rows,
            vec![2]
        );
        let mut f = triangle();
        with_spell_table(&mut f, &["a"], &[f64::NAN], &[3.0]);
        let v = f.validate(1);
        assert_eq!(
            v.resolved.unwrap().spell_table.unwrap().start,
            vec![f64::NEG_INFINITY]
        );
        // Roles and conflicts.
        let mut f = triangle();
        f.add_column(Which::NodeSpells, "id", strs(&["a"])).unwrap();
        f.set_spell_role(SpellRole::Id, "id").unwrap();
        assert!(f.validate(1).report.has(Code::E_TABLE_ROLE_UNDECLARED));
        let mut f = triangle();
        with_spell_table(&mut f, &["a"], &[0.0], &[1.0]);
        f.add_node_column("arrive", Column::Float(vec![0.0; 3]))
            .unwrap();
        f.set_node_role(NodeRole::NodeStart, "arrive").unwrap();
        assert!(f.validate(1).report.has(Code::E_DECLARATION_CONFLICT));
        let mut f = triangle();
        f.add_column(Which::NodeSpells, "id", ints(&[1])).unwrap();
        f.add_column(Which::NodeSpells, "s", strs(&["x"])).unwrap();
        f.add_column(Which::NodeSpells, "e", Column::Float(vec![1.0]))
            .unwrap();
        f.set_spell_role(SpellRole::Id, "id").unwrap();
        f.set_spell_role(SpellRole::Start, "s").unwrap();
        f.set_spell_role(SpellRole::End, "e").unwrap();
        let v = f.validate(1);
        assert!(v.report.has(Code::E_ID_TYPE) && v.report.has(Code::E_TIME_TYPE));
        let mut f = triangle();
        with_spell_table(&mut f, &["a"], &[0.0], &[1.0]);
        f.add_column(
            Which::NodeSpells,
            "idm",
            Column::Str(Strings::from_strs(&[None::<&str>])),
        )
        .unwrap();
        f.set_spell_role(SpellRole::Id, "idm").unwrap();
        let v = f.validate(1);
        assert_eq!(
            v.report.find(Code::E_MISSING_ENDPOINT).unwrap().table,
            Some(Which::NodeSpells)
        );
    }

    #[test]
    fn node_change_table() {
        use super::super::roles::ChangeRole;
        let mut f = triangle();
        f.add_column(Which::NodeChanges, "id", strs(&["a", "z", "b", "c"]))
            .unwrap();
        f.add_column(
            Which::NodeChanges,
            "t",
            Column::Float(vec![1.0, 2.0, f64::NAN, f64::INFINITY]),
        )
        .unwrap();
        let v = f.validate(1);
        assert!(v.report.has(Code::E_TABLE_ROLE_UNDECLARED));
        f.set_change_role(ChangeRole::Id, "id").unwrap();
        f.set_change_role(ChangeRole::Time, "t").unwrap();
        let v = f.validate(1);
        assert_eq!(
            codes(&v),
            vec![
                Code::E_UNKNOWN_NODE,
                Code::E_TIME_MISSING,
                Code::E_TIME_NOT_FINITE
            ]
        );
        assert_eq!(v.report.find(Code::E_UNKNOWN_NODE).unwrap().rows, vec![1]);
        assert_eq!(v.report.find(Code::E_TIME_MISSING).unwrap().rows, vec![2]);
        assert_eq!(
            v.report.find(Code::E_TIME_NOT_FINITE).unwrap().table,
            Some(Which::NodeChanges)
        );
        f.add_column(Which::NodeChanges, "ts", strs(&["x"; 4]))
            .unwrap();
        f.set_change_role(ChangeRole::Time, "ts").unwrap();
        assert!(f.validate(1).report.has(Code::E_TIME_TYPE));
        let mut f = triangle();
        f.add_column(Which::NodeChanges, "id", strs(&["a", "a"]))
            .unwrap();
        f.add_column(Which::NodeChanges, "t", ints(&[1, 2]))
            .unwrap();
        f.set_change_role(ChangeRole::Id, "id").unwrap();
        f.set_change_role(ChangeRole::Time, "t").unwrap();
        assert!(codes(&f.validate(1)).is_empty());
    }

    #[test]
    fn unobserved_pairs_are_kept_but_not_ties() {
        let mut f = Frame::new(false);
        f.add_edge_column("from", ints(&[1, 1, 2, 2])).unwrap();
        f.add_edge_column("to", ints(&[2, 3, 1, 3])).unwrap();
        f.add_edge_column("na", Column::Bool(vec![false, true, true, true], None))
            .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        let v = f.validate(1);
        // 2-1 is both observed (row 0) and unobserved (row 2).
        assert_eq!(
            codes(&v),
            vec![
                Code::N_MISSING_EDGES,
                Code::E_MISSING_CONFLICT,
                Code::N_PARALLEL
            ]
        );
        assert_eq!(
            v.report.find(Code::E_MISSING_CONFLICT).unwrap().rows,
            vec![2]
        );
        let r = v.resolved.unwrap();
        assert_eq!(r.missing, Some(vec![false, true, true, true]));
        assert!(r.observed(0) && !r.observed(1));
        assert_eq!(r.observed_count(), 1);
        f.declarations_mut().directed = true;
        let v = f.validate(1);
        assert!(v.report.is_valid(), "{:?}", v.report);
    }

    #[test]
    fn isolates_and_zero_based_ids() {
        let mut f = Frame::new(true);
        f.add_node_column("id", ints(&[0, 1, 2])).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column("from", ints(&[0])).unwrap();
        f.add_edge_column("to", ints(&[1])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = f.validate(1);
        assert_eq!(codes(&v), vec![Code::N_ISOLATES, Code::N_ZERO_BASED_IDS]);
        assert_eq!(v.report.find(Code::N_ISOLATES).unwrap().rows, vec![2]);
        // 1-based or sparse identifiers are not flagged.
        let mut f = Frame::new(true);
        f.add_edge_column("from", ints(&[1, 2])).unwrap();
        f.add_edge_column("to", ints(&[2, 3])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        assert!(codes(&f.validate(1)).is_empty());
        let mut f = Frame::new(true);
        f.add_edge_column("from", ints(&[0, 5])).unwrap();
        f.add_edge_column("to", ints(&[5, 0])).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        assert!(codes(&f.validate(1)).is_empty());
    }

    #[test]
    fn scans_agree_across_threads_and_chunks() {
        let len = 3 * CHUNK + 17;
        let pred = |i: usize| i % 1000 == 7;
        let (r1, c1) = scan_rows(len, 1, pred);
        let (r4, c4) = scan_rows(len, 4, pred);
        assert_eq!((r1.len(), c1), (ROW_CAP, (len + 992) / 1000));
        assert_eq!((r1, c1), (r4, c4));
        assert_eq!(scan_rows(0, 4, pred), (vec![], 0));
        let (rows, count) = scan_rows(len, 3, |i| i == len - 1);
        assert_eq!((rows, count), (vec![len - 1], 1));
    }

    #[test]
    fn large_frames_give_the_same_report_on_every_thread_count() {
        // More rows than a chunk, with faults spread across chunks and
        // nodes added from the edges, so the parallel lookup, the
        // row-order pass, and the chunked scans are all exercised.
        let m = 3 * CHUNK + 1234;
        let n = 50_000;
        let mut src = Vec::with_capacity(m);
        let mut tgt = Vec::with_capacity(m);
        let mut w = Vec::with_capacity(m);
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        for e in 0..m {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let a = (x % n) as i64;
            let b = ((x >> 20) % (n + 100)) as i64;
            src.push(Some(format!("v{a}")));
            tgt.push(if e % 70_001 == 5 {
                None
            } else {
                Some(format!("v{b}"))
            });
            w.push(if e % 100_003 == 7 {
                f64::NAN
            } else {
                1.0 + (e % 5) as f64
            });
        }
        let ids: Vec<Option<String>> = (0..n).map(|i| Some(format!("v{i}"))).collect();
        let mut f = Frame::new(true);
        f.add_node_column("id", Column::Str(Strings::from_strs(&ids)))
            .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_edge_column("from", Column::Str(Strings::from_strs(&src)))
            .unwrap();
        f.add_edge_column("to", Column::Str(Strings::from_strs(&tgt)))
            .unwrap();
        f.add_edge_column("w", Column::Float(w)).unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Strength));
        f.declarations_mut().add_nodes = true;
        let one = f.validate(1);
        let three = f.validate(3);
        assert_eq!(one.report, three.report);
        assert!(one.report.has(Code::E_MISSING_ENDPOINT));
        assert!(one.report.has(Code::E_WEIGHT_NOT_FINITE));
        assert!(one.resolved.is_none());
        // Without the missing endpoints, the added nodes come in first
        // appearance order on every thread count.
        let fixed: Vec<Option<String>> = tgt
            .iter()
            .map(|t| t.clone().or_else(|| Some("v1".into())))
            .collect();
        let mut g = f.clone();
        g.add_edge_column("to2", Column::Str(Strings::from_strs(&fixed)))
            .unwrap();
        g.set_edge_role(EdgeRole::Target, "to2").unwrap();
        let one = g.validate(1);
        let three = g.validate(3);
        assert_eq!(one.report, three.report);
        let (r1, r3) = (one.resolved.unwrap(), three.resolved.unwrap());
        assert_eq!(r1.added, r3.added);
        assert_eq!((r1.src, r1.tgt), (r3.src, r3.tgt));
        assert!(r1.n > n as usize);
        assert!(three.report.has(Code::N_NODES_ADDED));
    }

    #[test]
    fn hasher_separates_short_keys() {
        let mut h1 = Fast::default();
        h1.write(b"ab");
        let mut h2 = Fast::default();
        h2.write(b"ab\0");
        assert_ne!(h1.finish(), h2.finish());
        let mut h3 = Fast::default();
        h3.write(b"abcdefghij");
        let mut h4 = Fast::default();
        h4.write(b"abcdefghik");
        assert_ne!(h3.finish(), h4.finish());
    }
}
