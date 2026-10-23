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
use super::roles::{CoordinateKind, EdgeRole, NodeRole, WeightType};
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
    /// Parallel-edge groups, when every distinguishing role resolved.
    pub groups: Option<Groups>,
    /// Number of self-loops.
    pub loops: usize,
    /// Whether edges are directed (the declaration, kept here so that
    /// later operations on the resolved edges need no frame).
    pub directed: bool,
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

/// The result of mapping identifiers to indices.
struct Mapping {
    n: usize,
    src: Vec<u32>,
    tgt: Vec<u32>,
    added: Vec<(usize, bool)>,
    ok: bool,
}

/// Maps node identifiers and edge endpoints to indices with one hash
/// table, reporting missing, duplicate, and unknown identifiers.
fn map_ids<K: Hash + Eq + Copy>(
    frame: &Frame,
    node_key: impl Fn(usize) -> Option<K>,
    src_key: impl Fn(usize) -> Option<K>,
    tgt_key: impl Fn(usize) -> Option<K>,
    empty: impl Fn(K) -> bool,
    report: &mut Report,
) -> Mapping {
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

    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    let mut added = Vec::new();
    let mut missing = Vec::new();
    let mut empties = Vec::new();
    let mut unknown = Vec::new();
    let (mut n_missing, mut n_empty, mut n_unknown) = (0, 0, 0);
    for e in 0..m {
        for (is_target, key) in [(false, src_key(e)), (true, tgt_key(e))] {
            let idx = match key {
                None => {
                    n_missing += 1;
                    push_capped(&mut missing, e);
                    0
                }
                Some(k) if empty(k) => {
                    n_empty += 1;
                    push_capped(&mut empties, e);
                    0
                }
                Some(k) => match index.get(&k) {
                    Some(&i) => i,
                    None if add_nodes => {
                        let i = (n0 + added.len()) as u32;
                        index.insert(k, i);
                        added.push((e, is_target));
                        i
                    }
                    None => {
                        n_unknown += 1;
                        push_capped(&mut unknown, e);
                        0
                    }
                },
            };
            if is_target {
                tgt.push(idx);
            } else {
                src.push(idx);
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
    Mapping {
        n,
        src,
        tgt,
        added,
        ok,
    }
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
fn resolve_endpoints(frame: &Frame, kind: Kind, report: &mut Report) -> Mapping {
    let src = IdView::of(frame.edge_role(EdgeRole::Source).expect("checked").1).expect("checked");
    let tgt = IdView::of(frame.edge_role(EdgeRole::Target).expect("checked").1).expect("checked");
    let id = frame
        .node_role(NodeRole::Id)
        .and_then(|(_, c)| IdView::of(c));
    let no_nodes = frame.node_rows() == 0 || id.is_none();
    match kind {
        Kind::Int => map_ids(
            frame,
            |i| if no_nodes { None } else { id.unwrap().int(i) },
            |e| src.int(e),
            |e| tgt.int(e),
            |_| false,
            report,
        ),
        Kind::Float => map_ids(
            frame,
            |i| if no_nodes { None } else { id.unwrap().float(i) },
            |e| src.float(e),
            |e| tgt.float(e),
            |_| false,
            report,
        ),
        _ => map_ids(
            frame,
            |i| if no_nodes { None } else { id.unwrap().str(i) },
            |e| src.str(e),
            |e| tgt.str(e),
            |k| k.is_empty(),
            report,
        ),
    }
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
            mapping = Some(resolve_endpoints(frame, kind, &mut report));
        }
    }

    // Declarations that contradict each other, whatever the data.
    let mode_role = frame.node_role(NodeRole::Mode);
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

    // 2. Structure-dependent checks and declared properties.
    let mapping = match mapping {
        Some(map) if map.ok => map,
        _ => {
            return Validation {
                report,
                resolved: None,
            }
        }
    };
    let Mapping {
        n, src, tgt, added, ..
    } = mapping;
    let (src_name, tgt_name) = endpoint_names(frame);

    // Loops.
    let (loop_rows, loops) = scan_rows(m, threads, |e| src[e] == tgt[e]);
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

    // Edges within a mode.
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
            |e| codes[src[e] as usize] == codes[tgt[e] as usize],
        );
    }

    // Edge spells within node spells.
    if let (Some(s), Some(ns)) = (&start, &node_spells) {
        let e_ = end.as_deref();
        let inside = |node: u32, s: f64, e: f64| {
            let (a, b) = ns[node as usize];
            a <= s && e <= b && (e > s || s < b || a == b)
        };
        report_rows(
            &mut report,
            Code::E_SPELL_OUTSIDE_NODE,
            Which::Edges,
            "an edge can only be active while both of its endpoints are present; the listed edges fall outside a node's presence spell",
            m,
            threads,
            |i| {
                let e = e_.map_or(s[i], |e| e[i]);
                !(inside(src[i], s[i], e) && inside(tgt[i], s[i], e))
            },
        );
    }

    // Parallel edges, keyed on endpoints, relation, and spell.
    let groups_possible = (frame.edge_role(EdgeRole::Relation).is_none() || relations.is_some())
        && (!(has_start || has_wave) || start.is_some());
    let mut groups = None;
    if groups_possible {
        let g = edge_groups(
            &src,
            &tgt,
            decl.directed,
            relations.as_ref().map(|r| &r.codes[..]),
            start.as_deref(),
            end.as_deref(),
        );
        let parallel_rows: Vec<usize> = (0..m)
            .filter(|&e| g.first[g.group[e] as usize] != e)
            .take(ROW_CAP)
            .collect();
        let parallel = m - g.sizes.len();
        let what = {
            let mut parts = vec![format!("columns `{src_name}` and `{tgt_name}`")];
            if let Some((rn, _)) = frame.edge_role(EdgeRole::Relation) {
                parts.push(format!("relation `{rn}`"));
            }
            if has_start || has_wave {
                parts.push("time".to_string());
            }
            parts.join(", ")
        };
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
            groups,
            loops,
            directed: decl.directed,
        }),
    }
}

/// Groups edges by their key: endpoints (ordered when directed), relation
/// code, and spell, as the module documentation defines parallel edges.
pub(super) fn edge_groups(
    src: &[u32],
    tgt: &[u32],
    directed: bool,
    relation: Option<&[u32]>,
    start: Option<&[f64]>,
    end: Option<&[f64]>,
) -> Groups {
    let m = src.len();
    let mut map: FastMap<(u32, u32, u32, u64, u64), u32> = fast_map(m);
    let mut group = Vec::with_capacity(m);
    let mut sizes: Vec<u32> = Vec::new();
    let mut first = Vec::new();
    for e in 0..m {
        let (a, b) = if directed || src[e] <= tgt[e] {
            (src[e], tgt[e])
        } else {
            (tgt[e], src[e])
        };
        let r = relation.map_or(0, |r| r[e]);
        let s = start.map_or(0, |s| (s[e] + 0.0).to_bits());
        let t = end.map_or(s, |t| (t[e] + 0.0).to_bits());
        let next = sizes.len() as u32;
        let g = *map.entry((a, b, r, s, t)).or_insert(next);
        if g == next {
            sizes.push(0);
            first.push(e);
        }
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
