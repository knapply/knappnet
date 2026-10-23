//! C ABI for the data model: building a frame from typed column buffers,
//! assigning roles and declarations, validating it, reading the
//! diagnostics, building the graph, and reading the derived properties.
//!
//! A frame handle is filled in steps: `net_frame_new`, then columns
//! (`net_frame_add_*`), roles (`net_frame_set_role`), the weight type,
//! declarations, and geometry flags. `net_frame_validate` runs every
//! check and keeps the diagnostics on the handle; `net_frame_build`
//! validates (if not done yet) and returns a graph handle. Any change to
//! the frame after validation discards the diagnostics and the build.
//!
//! Integer codes used by these functions:
//! - `table`: 0 = nodes, 1 = edges, 2 = layers, 3 = node_spells,
//!   4 = node_changes.
//! - node `role`: 0 = id, 1 = mode, 2 = node_start, 3 = node_end, 4 = x,
//!   5 = y, 6 = node_missing, 7 = ego. Edge `role`: 0 = source,
//!   1 = target, 2 = weight, 3 = relation, 4 = sign, 5 = start, 6 = end,
//!   7 = wave, 8 = edge_missing, 9 = source_layer, 10 = target_layer,
//!   11 = ego. Layer table `role`: 0 = layer. Node spells `role`: 0 = id,
//!   1 = start, 2 = end. Node changes `role`: 0 = id, 1 = time.
//! - `weight_type`: 0 = distance, 1 = strength, -1 = none.
//! - declaration `which`: 0 = directed, 1 = loops, 2 = multiple,
//!   3 = between_modes_only, 4 = add_nodes, 5 = egos_are_nodes,
//!   6 = couplings_only, 7 = overlapping_spells; `value` 0 = false,
//!   1 = true, -1 = not declared (for 1 to 3 and 5 to 7).
//! - coordinate `kind`: 0 = unknown, 1 = longitude and latitude,
//!   2 = planar.
//! - diagnostic `severity`: 0 = error, 1 = note; diagnostic `table`:
//!   a table code, or -1 for the network as a whole.
//! - `combine`: 0 = sum, 1 = mean, 2 = min, 3 = max, 4 = first,
//!   5 = count.
//! - edge `class` (`net_frame_edge_layers`): 0 = intralayer,
//!   1 = coupling, 2 = interlayer.
//! - subgraph `which` (`net_frame_subgraph`): 0 = one layer, 1 = one ego.
//! - summary numbers (`net_frame_summary`, [`SUMMARY_LEN`] doubles):
//!   see the header comment of `net_frame_summary`.
//!
//! Missing masks are byte arrays (non-zero = missing) and may be NULL
//! when nothing is missing. Logical values are `int32_t` (non-zero =
//! true). Text columns are UTF-8 bytes with `n + 1` offsets.

use std::ffi::{c_char, CStr, CString};

use engine::frame::multiplicity::{is_simple, multiplicity};
use engine::frame::{
    build_validated, ego_subgraph, layer_subgraph, missing_pairs, simplify, summarize, Built,
    ChangeRole, Column, Combine, CoordinateKind, Declaration, EdgeRole, Frame, LayerRole, NodeRole,
    Severity, SpellRole, Strings, Subgraph, Summary, Validation, WeightType, Which,
};
use engine::Error;

use crate::{count, deref, guard, slice, slice_mut, thread_count, GraphHandle};

/// Number of doubles `net_frame_summary` writes.
pub const SUMMARY_LEN: usize = 37;

/// Opaque frame handle.
pub struct FrameHandle {
    frame: Frame,
    validation: Option<Validation>,
    /// C strings of each diagnostic's code and message, aligned with the
    /// report.
    strings: Vec<(CString, CString)>,
    built: Option<Built>,
    summary: Option<Summary>,
    labels: Vec<Vec<CString>>,
}

impl FrameHandle {
    fn touch(&mut self) {
        self.validation = None;
        self.strings.clear();
        self.built = None;
        self.summary = None;
        self.labels.clear();
    }

    fn validate(&mut self, threads: usize) {
        if self.validation.is_none() {
            let v = self.frame.validate(threads);
            self.strings = v
                .report
                .diagnostics
                .iter()
                .map(|d| (cstring(d.code.as_str()), cstring(&d.message)))
                .collect();
            self.validation = Some(v);
        }
    }

    fn built(&mut self, threads: usize) -> Result<&Built, Error> {
        if self.built.is_none() {
            self.validate(threads);
            let v = self.validation.clone().expect("validated");
            self.built = Some(build_validated(&self.frame, v)?);
        }
        Ok(self.built.as_ref().expect("built"))
    }

    fn summary(&mut self, threads: usize) -> Result<&Summary, Error> {
        if self.summary.is_none() {
            let built = self.built(threads)?.clone();
            let s = summarize(&self.frame, &built);
            let strs = |v: &[String]| v.iter().map(|s| cstring(s)).collect::<Vec<_>>();
            self.labels = vec![
                strs(&s.mode_labels),
                strs(&s.relation_labels),
                s.weight_column.iter().map(|s| cstring(s)).collect(),
                s.crs.iter().map(|s| cstring(s)).collect(),
                strs(&s.layer_labels),
                strs(&s.ego_labels),
            ];
            self.summary = Some(s);
        }
        Ok(self.summary.as_ref().expect("summarized"))
    }
}

fn cstring(s: &str) -> CString {
    CString::new(s.replace('\0', "")).expect("NUL bytes removed")
}

unsafe fn frame_mut<'a>(h: *mut FrameHandle) -> Result<&'a mut FrameHandle, Error> {
    h.as_mut()
        .ok_or_else(|| Error::InvalidArgument("frame handle is null".into()))
}

unsafe fn frame_ref<'a>(h: *const FrameHandle) -> Result<&'a FrameHandle, Error> {
    deref(h, "frame")
}

unsafe fn name_arg<'a>(p: *const c_char, what: &str) -> Result<&'a str, Error> {
    if p.is_null() {
        return Err(Error::InvalidArgument(format!("{what} is null")));
    }
    CStr::from_ptr(p)
        .to_str()
        .map_err(|_| Error::InvalidArgument(format!("{what} is not valid UTF-8")))
}

fn which(table: i32) -> Result<Which, Error> {
    Which::from_code(table)
}

fn table_code(which: Option<Which>) -> i32 {
    match which {
        Some(w) => Which::ALL.iter().position(|x| *x == w).expect("listed") as i32,
        None => -1,
    }
}

unsafe fn mask(missing: *const u8, n: usize) -> Option<Vec<bool>> {
    if missing.is_null() {
        None
    } else {
        let m: Vec<bool> = slice(missing, n).iter().map(|&x| x != 0).collect();
        m.iter().any(|&x| x).then_some(m)
    }
}

unsafe fn add_column(
    h: *mut FrameHandle,
    table: i32,
    name: *const c_char,
    column: Column,
) -> Result<(), Error> {
    let f = frame_mut(h)?;
    let name = name_arg(name, "column name")?;
    f.touch();
    f.frame.add_column(which(table)?, name, column)?;
    Ok(())
}

/// A new, empty frame; directed != 0 for a directed network. Release
/// with `net_frame_free`.
#[no_mangle]
pub extern "C" fn net_frame_new(directed: i32) -> *mut FrameHandle {
    Box::into_raw(Box::new(FrameHandle {
        frame: Frame::new(directed != 0),
        validation: None,
        strings: Vec::new(),
        built: None,
        summary: None,
        labels: Vec::new(),
    }))
}

/// Releases a frame. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_frame_free(h: *mut FrameHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Gives the frame a node table of `n` rows (needed only when no node
/// column will be added; every node column must then have `n` rows).
///
/// # Safety
/// `h` must be a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_node_rows(
    h: *mut FrameHandle,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.touch();
        f.frame.set_node_rows(count(n, "node count")?)
    })
}

/// Adds an integer column of `n` values; `missing` may be null.
///
/// # Safety
/// `values` must hold `n` elements and `missing` null or `n` bytes.
#[no_mangle]
pub unsafe extern "C" fn net_frame_add_int(
    h: *mut FrameHandle,
    table: i32,
    name: *const c_char,
    values: *const i64,
    missing: *const u8,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "row count")?;
        let column = Column::Int(slice(values, n).to_vec(), mask(missing, n));
        add_column(h, table, name, column)
    })
}

/// Adds a double column of `n` values; NaN is missing.
///
/// # Safety
/// `values` must hold `n` elements.
#[no_mangle]
pub unsafe extern "C" fn net_frame_add_double(
    h: *mut FrameHandle,
    table: i32,
    name: *const c_char,
    values: *const f64,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "row count")?;
        add_column(h, table, name, Column::Float(slice(values, n).to_vec()))
    })
}

/// Adds a logical column of `n` values (non-zero = true); `missing` may
/// be null.
///
/// # Safety
/// `values` must hold `n` elements and `missing` null or `n` bytes.
#[no_mangle]
pub unsafe extern "C" fn net_frame_add_logical(
    h: *mut FrameHandle,
    table: i32,
    name: *const c_char,
    values: *const i32,
    missing: *const u8,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "row count")?;
        let v: Vec<bool> = slice(values, n).iter().map(|&x| x != 0).collect();
        add_column(h, table, name, Column::Bool(v, mask(missing, n)))
    })
}

/// Adds a text column of `n` values: value `i` is bytes
/// `offsets[i]..offsets[i + 1]` of `bytes` (UTF-8, no terminators);
/// `missing` may be null.
///
/// # Safety
/// `offsets` must hold `n + 1` elements, `bytes` the last offset's
/// bytes, and `missing` null or `n` bytes.
#[no_mangle]
pub unsafe extern "C" fn net_frame_add_string(
    h: *mut FrameHandle,
    table: i32,
    name: *const c_char,
    offsets: *const i64,
    bytes: *const c_char,
    missing: *const u8,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = count(n, "row count")?;
        let offs = slice(offsets, n + 1);
        let mut converted = Vec::with_capacity(n + 1);
        for &o in offs {
            converted.push(count(o, "text offset")?);
        }
        let total = if n == 0 && offs.is_empty() {
            converted.push(0);
            0
        } else {
            converted[n]
        };
        let data = slice(bytes as *const u8, total).to_vec();
        let strings = Strings::new(converted, data, mask(missing, n))?;
        add_column(h, table, name, Column::Str(strings))
    })
}

/// Assigns `role` of `table` to the column called `name`, or clears the
/// role when `name` is null. The role codes depend on the table (see
/// the module codes).
///
/// # Safety
/// `h` must be a live handle; `name` null or NUL-terminated.
#[no_mangle]
pub unsafe extern "C" fn net_frame_set_role(
    h: *mut FrameHandle,
    table: i32,
    role: i32,
    name: *const c_char,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.touch();
        let name = if name.is_null() {
            None
        } else {
            Some(name_arg(name, "column name")?)
        };
        match which(table)? {
            Which::Nodes => {
                let role = NodeRole::from_code(role)?;
                match name {
                    Some(name) => f.frame.set_node_role(role, name),
                    None => {
                        f.frame.clear_node_role(role);
                        Ok(())
                    }
                }
            }
            Which::Edges => {
                let role = EdgeRole::from_code(role)?;
                match name {
                    Some(name) => f.frame.set_edge_role(role, name),
                    None => {
                        f.frame.clear_edge_role(role);
                        Ok(())
                    }
                }
            }
            Which::Layers => {
                let role = LayerRole::from_code(role)?;
                match name {
                    Some(name) => f.frame.set_layer_role(role, name),
                    None => {
                        f.frame.clear_layer_role(role);
                        Ok(())
                    }
                }
            }
            Which::NodeSpells => {
                let role = SpellRole::from_code(role)?;
                match name {
                    Some(name) => f.frame.set_spell_role(role, name),
                    None => {
                        f.frame.clear_spell_role(role);
                        Ok(())
                    }
                }
            }
            Which::NodeChanges => {
                let role = ChangeRole::from_code(role)?;
                match name {
                    Some(name) => f.frame.set_change_role(role, name),
                    None => {
                        f.frame.clear_change_role(role);
                        Ok(())
                    }
                }
            }
        }
    })
}

/// Declares the meaning of the weights: 0 = distance, 1 = strength,
/// -1 = none.
///
/// # Safety
/// `h` must be a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_set_weight_type(
    h: *mut FrameHandle,
    weight_type: i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.touch();
        let t = if weight_type < 0 {
            None
        } else {
            Some(WeightType::from_code(weight_type)?)
        };
        f.frame.set_weight_type(t);
        Ok(())
    })
}

/// Sets a declaration (see the module codes).
///
/// # Safety
/// `h` must be a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_declare(
    h: *mut FrameHandle,
    which: i32,
    value: i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.touch();
        f.frame
            .declarations_mut()
            .set(Declaration::from_code(which)?, value);
        Ok(())
    })
}

/// Records whether `table` carries geometry and, if so, its coordinate
/// reference system (`crs` may be null when unknown).
///
/// # Safety
/// `h` must be a live handle; `crs` null or NUL-terminated.
#[no_mangle]
pub unsafe extern "C" fn net_frame_set_geometry(
    h: *mut FrameHandle,
    table: i32,
    present: i32,
    crs: *const c_char,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.touch();
        let crs = if crs.is_null() {
            None
        } else {
            Some(name_arg(crs, "crs")?.to_string())
        };
        let g = f.frame.geometry_mut();
        match which(table)? {
            Which::Nodes => {
                g.nodes = present != 0;
                g.node_crs = crs;
            }
            Which::Edges => {
                g.edges = present != 0;
                g.edge_crs = crs;
            }
            other => {
                return Err(Error::InvalidArgument(format!(
                    "only `nodes` and `edges` carry geometry, not `{other}`"
                )))
            }
        }
        Ok(())
    })
}

/// Says what the representative coordinates (roles x and y) are:
/// 0 = unknown, 1 = longitude and latitude, 2 = planar.
///
/// # Safety
/// `h` must be a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_set_coordinates(
    h: *mut FrameHandle,
    kind: i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.touch();
        f.frame.geometry_mut().coordinates = CoordinateKind::from_code(kind)?;
        Ok(())
    })
}

/// Runs every check. Writes the number of errors and notes; the
/// diagnostics are then read with `net_frame_diagnostic`. Returns 0 even
/// when the frame is invalid (the diagnostics say so); non-zero only
/// when the call itself fails.
///
/// # Safety
/// `h` must be a live handle; `errors` and `notes` writable or null.
#[no_mangle]
pub unsafe extern "C" fn net_frame_validate(
    h: *mut FrameHandle,
    threads: i32,
    errors: *mut i64,
    notes: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        f.validate(thread_count(threads));
        let report = &f.validation.as_ref().expect("validated").report;
        if !errors.is_null() {
            *errors = report.errors().count() as i64;
        }
        if !notes.is_null() {
            *notes = report.notes().count() as i64;
        }
        Ok(())
    })
}

/// Number of diagnostics, or -1 for a null or unvalidated handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_diagnostic_count(h: *const FrameHandle) -> i64 {
    match frame_ref(h) {
        Ok(f) => f
            .validation
            .as_ref()
            .map_or(-1, |v| v.report.diagnostics.len() as i64),
        Err(_) => -1,
    }
}

/// Reads diagnostic `i`: its code and message (valid while the handle
/// lives and the frame is unchanged), severity (0 = error, 1 = note),
/// table (a table code, or -1 for the network), the number of affected
/// rows in all, and how many of them `net_frame_diagnostic_rows` lists.
/// Any output pointer may be null.
///
/// # Safety
/// `h` must be a live, validated handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_diagnostic(
    h: *const FrameHandle,
    i: i64,
    code: *mut *const c_char,
    message: *mut *const c_char,
    severity: *mut i32,
    table: *mut i32,
    count_out: *mut i64,
    n_rows: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_ref(h)?;
        let v = f
            .validation
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("the frame has not been validated".into()))?;
        let i = count(i, "diagnostic index")?;
        let d = v
            .report
            .diagnostics
            .get(i)
            .ok_or_else(|| Error::InvalidArgument(format!("no diagnostic {i}")))?;
        if !code.is_null() {
            *code = f.strings[i].0.as_ptr();
        }
        if !message.is_null() {
            *message = f.strings[i].1.as_ptr();
        }
        if !severity.is_null() {
            *severity = match d.severity() {
                Severity::Error => 0,
                Severity::Note => 1,
            };
        }
        if !table.is_null() {
            *table = table_code(d.table);
        }
        if !count_out.is_null() {
            *count_out = d.count as i64;
        }
        if !n_rows.is_null() {
            *n_rows = d.rows.len() as i64;
        }
        Ok(())
    })
}

/// Copies the listed 0-based rows of diagnostic `i` into `rows`.
///
/// # Safety
/// `rows` must hold as many values as `net_frame_diagnostic` reported.
#[no_mangle]
pub unsafe extern "C" fn net_frame_diagnostic_rows(
    h: *const FrameHandle,
    i: i64,
    rows: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_ref(h)?;
        let v = f
            .validation
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("the frame has not been validated".into()))?;
        let i = count(i, "diagnostic index")?;
        let d = v
            .report
            .diagnostics
            .get(i)
            .ok_or_else(|| Error::InvalidArgument(format!("no diagnostic {i}")))?;
        let out = slice_mut(rows, d.rows.len());
        for (o, &r) in out.iter_mut().zip(&d.rows) {
            *o = r as i64;
        }
        Ok(())
    })
}

/// Validates (if not done yet) and builds the graph. Fails, with the
/// errors as the message, when the frame is invalid. Release `*out` with
/// `net_graph_free`; the frame keeps its own copy for the accessors
/// below.
///
/// # Safety
/// `h` must be a live handle; `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_frame_build(
    h: *mut FrameHandle,
    threads: i32,
    out: *mut *mut GraphHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        let built = f.built(thread_count(threads))?;
        *out = Box::into_raw(Box::new(GraphHandle(built.graph.clone())));
        Ok(())
    })
}

/// Number of nodes of the built graph (table rows plus added nodes), or
/// -1 before a build.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_node_count(h: *const FrameHandle) -> i64 {
    match frame_ref(h) {
        Ok(f) => f.built.as_ref().map_or(-1, |b| b.graph.node_count() as i64),
        Err(_) => -1,
    }
}

/// Number of nodes added from the edges, or -1 before a build.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_added_count(h: *const FrameHandle) -> i64 {
    match frame_ref(h) {
        Ok(f) => f.built.as_ref().map_or(-1, |b| b.map.added.len() as i64),
        Err(_) => -1,
    }
}

/// For each added node, in node order after the table's rows, the edge
/// row its identifier first appeared in and whether that was the target
/// (1) or the source (0) endpoint, so the caller can take the identifier
/// from its own endpoint columns.
///
/// # Safety
/// `rows` and `is_target` must hold `net_frame_added_count` values.
#[no_mangle]
pub unsafe extern "C" fn net_frame_added(
    h: *const FrameHandle,
    rows: *mut i64,
    is_target: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let k = b.map.added.len();
        let rows = slice_mut(rows, k);
        let targets = slice_mut(is_target, k);
        for (i, &(row, t)) in b.map.added.iter().enumerate() {
            rows[i] = row as i64;
            targets[i] = t as i32;
        }
        Ok(())
    })
}

unsafe fn built_ref<'a>(h: *const FrameHandle) -> Result<&'a Built, Error> {
    frame_ref(h)?
        .built
        .as_ref()
        .ok_or_else(|| Error::InvalidArgument("the frame has not been built".into()))
}

/// The 0-based node index of every edge's source and target.
///
/// # Safety
/// `src` and `tgt` must hold m values (the edge table's rows).
#[no_mangle]
pub unsafe extern "C" fn net_frame_edge_index(
    h: *const FrameHandle,
    src: *mut i64,
    tgt: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let m = b.resolved.src.len();
        for (o, &i) in slice_mut(src, m).iter_mut().zip(&b.resolved.src) {
            *o = i as i64;
        }
        for (o, &i) in slice_mut(tgt, m).iter_mut().zip(&b.resolved.tgt) {
            *o = i as i64;
        }
        Ok(())
    })
}

/// The integer codes of the modes (`which == 0`, n values: one per node
/// of the node table, in order of first appearance of each mode label),
/// relations (`which == 1`, m values), node egos (`which == 2`, one per
/// node), or edge egos (`which == 3`, m values). Fails when the role is
/// not set.
///
/// # Safety
/// `out` must hold as many values as the selector says.
#[no_mangle]
pub unsafe extern "C" fn net_frame_codes(
    h: *const FrameHandle,
    which: i32,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let (codes, role): (Option<&[u32]>, &str) = match which {
            0 => (b.resolved.modes.as_ref().map(|d| &d.codes[..]), "mode"),
            1 => (
                b.resolved.relations.as_ref().map(|d| &d.codes[..]),
                "relation",
            ),
            2 => (b.resolved.egos.as_ref().map(|e| &e.node[..]), "ego"),
            3 => (b.resolved.egos.as_ref().map(|e| &e.edge[..]), "ego"),
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown codes selector {which}"
                )))
            }
        };
        let codes =
            codes.ok_or_else(|| Error::InvalidArgument(format!("the frame has no {role} role")))?;
        for (o, &c) in slice_mut(out, codes.len()).iter_mut().zip(codes) {
            *o = c as i64;
        }
        Ok(())
    })
}

/// The layer codes of every edge's endpoints (indices into the summary's
/// layer labels) and each edge's class (0 = intralayer, 1 = coupling,
/// 2 = interlayer); any output may be null. Fails without layer roles.
///
/// # Safety
/// Each non-null output must hold m values.
#[no_mangle]
pub unsafe extern "C" fn net_frame_edge_layers(
    h: *const FrameHandle,
    src_layer: *mut i64,
    tgt_layer: *mut i64,
    class: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let r = &b.resolved;
        let l = r
            .layers
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("the frame has no layer roles".into()))?;
        let m = l.src.len();
        if !src_layer.is_null() {
            for (o, &c) in slice_mut(src_layer, m).iter_mut().zip(&l.src) {
                *o = c as i64;
            }
        }
        if !tgt_layer.is_null() {
            for (o, &c) in slice_mut(tgt_layer, m).iter_mut().zip(&l.tgt) {
                *o = c as i64;
            }
        }
        if !class.is_null() {
            for (e, o) in slice_mut(class, m).iter_mut().enumerate() {
                *o = l.class(e, &r.src, &r.tgt).code();
            }
        }
        Ok(())
    })
}

/// Number of edges of the built graph (the observed ties), or -1 before
/// a build.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_graph_edge_count(h: *const FrameHandle) -> i64 {
    match frame_ref(h) {
        Ok(f) => f.built.as_ref().map_or(-1, |b| b.edges.count as i64),
        Err(_) => -1,
    }
}

/// The edge table row of every graph edge, in graph edge order (the
/// identity unless rows are flagged as unobserved pairs).
///
/// # Safety
/// `out` must hold `net_frame_graph_edge_count` values.
#[no_mangle]
pub unsafe extern "C" fn net_frame_edge_rows(
    h: *const FrameHandle,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        for (e, o) in slice_mut(out, b.edges.count).iter_mut().enumerate() {
            *o = b.edges.row(e) as i64;
        }
        Ok(())
    })
}

/// The unobserved pairs: the node indices and the edge table row of
/// every row flagged as missing, in row order. Any output may be null.
/// The count is the summary's number of unobserved edges.
///
/// # Safety
/// Each non-null output must hold as many values as there are flagged
/// rows.
#[no_mangle]
pub unsafe extern "C" fn net_frame_missing_pairs(
    h: *const FrameHandle,
    src: *mut i64,
    tgt: *mut i64,
    rows: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let pairs = missing_pairs(b);
        let k = pairs.len();
        if !src.is_null() {
            for (o, p) in slice_mut(src, k).iter_mut().zip(&pairs) {
                *o = p.0 as i64;
            }
        }
        if !tgt.is_null() {
            for (o, p) in slice_mut(tgt, k).iter_mut().zip(&pairs) {
                *o = p.1 as i64;
            }
        }
        if !rows.is_null() {
            for (o, p) in slice_mut(rows, k).iter_mut().zip(&pairs) {
                *o = p.2 as i64;
            }
        }
        Ok(())
    })
}

fn subgraph_of(b: &Built, which: i32, i: i64) -> Result<Subgraph, Error> {
    let i = count(i, "subgraph index")?;
    match which {
        0 => layer_subgraph(b, i),
        1 => ego_subgraph(b, i),
        _ => Err(Error::InvalidArgument(format!(
            "unknown subgraph selector {which}; 0 = layer, 1 = ego"
        ))),
    }
}

/// The size of the subgraph of layer `i` (`which == 0`: every node and
/// the observed intralayer edges of that layer) or of ego `i`
/// (`which == 1`: the nodes and observed edges of that ego): writes its
/// node and edge counts.
///
/// # Safety
/// `nodes` and `edges` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_frame_subgraph_size(
    h: *const FrameHandle,
    which: i32,
    i: i64,
    nodes: *mut i64,
    edges: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let s = subgraph_of(built_ref(h)?, which, i)?;
        *nodes = s.nodes.len() as i64;
        *edges = s.edges.len() as i64;
        Ok(())
    })
}

/// Builds the subgraph described by `net_frame_subgraph_size` into
/// `*out` (release with `net_graph_free`) and writes, for each of its
/// nodes, the built graph's node index, and for each of its edges, the
/// edge table row; either array may be null.
///
/// # Safety
/// `node_index` and `edge_rows` must be null or hold the counts
/// `net_frame_subgraph_size` reported; `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_frame_subgraph(
    h: *const FrameHandle,
    which: i32,
    i: i64,
    out: *mut *mut GraphHandle,
    node_index: *mut i64,
    edge_rows: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let s = subgraph_of(built_ref(h)?, which, i)?;
        if !node_index.is_null() {
            for (o, &k) in slice_mut(node_index, s.nodes.len())
                .iter_mut()
                .zip(&s.nodes)
            {
                *o = k as i64;
            }
        }
        if !edge_rows.is_null() {
            for (o, &k) in slice_mut(edge_rows, s.edges.len()).iter_mut().zip(&s.edges) {
                *o = k as i64;
            }
        }
        *out = Box::into_raw(Box::new(GraphHandle(s.graph)));
        Ok(())
    })
}

/// Number of presence spells the built frame resolved: one per node from
/// the `node_start` and `node_end` roles, the sorted rows of the node
/// spells table, or 0 without either. -1 before a build.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_node_spell_count(h: *const FrameHandle) -> i64 {
    match frame_ref(h) {
        Ok(f) => f.built.as_ref().map_or(-1, |b| {
            let r = &b.resolved;
            if let Some(t) = &r.spell_table {
                t.rows.len() as i64
            } else {
                r.node_spells.as_ref().map_or(0, |s| s.len() as i64)
            }
        }),
        Err(_) => -1,
    }
}

/// The resolved presence spells, sorted by node and start: `offsets`
/// (n + 1 values) delimit each node's spells in `start` and `end`
/// (`net_frame_node_spell_count` values each, infinities for unbounded
/// ends), and `rows` gives each spell's row of the node spells table
/// (or the node's row when the spells come from the node roles). Any
/// output may be null.
///
/// # Safety
/// Outputs must be null or hold the sizes above.
#[no_mangle]
pub unsafe extern "C" fn net_frame_node_spells(
    h: *const FrameHandle,
    offsets: *mut i64,
    start: *mut f64,
    end: *mut f64,
    rows: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let r = &b.resolved;
        let n = r.n;
        let (offs, s, e, rw): (Vec<usize>, Vec<f64>, Vec<f64>, Vec<usize>) =
            if let Some(t) = &r.spell_table {
                (
                    t.offsets.clone(),
                    t.start.clone(),
                    t.end.clone(),
                    t.rows.clone(),
                )
            } else if let Some(sp) = &r.node_spells {
                (
                    (0..=n).collect(),
                    sp.iter().map(|x| x.0).collect(),
                    sp.iter().map(|x| x.1).collect(),
                    (0..n).collect(),
                )
            } else {
                (vec![0; n + 1], Vec::new(), Vec::new(), Vec::new())
            };
        if !offsets.is_null() {
            for (o, &k) in slice_mut(offsets, n + 1).iter_mut().zip(&offs) {
                *o = k as i64;
            }
        }
        if !start.is_null() {
            slice_mut(start, s.len()).copy_from_slice(&s);
        }
        if !end.is_null() {
            slice_mut(end, e.len()).copy_from_slice(&e);
        }
        if !rows.is_null() {
            for (o, &k) in slice_mut(rows, rw.len()).iter_mut().zip(&rw) {
                *o = k as i64;
            }
        }
        Ok(())
    })
}

/// Number of edges sharing each edge's key (endpoints, relation, spell):
/// 1 for an edge without parallels. `out` holds m values.
///
/// # Safety
/// `out` must hold m values.
#[no_mangle]
pub unsafe extern "C" fn net_frame_multiplicity(
    h: *const FrameHandle,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let g = b
            .resolved
            .groups
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("parallel edges were not resolved".into()))?;
        let mult = multiplicity(g);
        for (o, &k) in slice_mut(out, mult.len()).iter_mut().zip(&mult) {
            *o = k as i64;
        }
        Ok(())
    })
}

/// Whether the network has neither self-loops nor parallel edges (1) or
/// not (0).
///
/// # Safety
/// `out` must be writable.
#[no_mangle]
pub unsafe extern "C" fn net_frame_is_simple(
    h: *const FrameHandle,
    out: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        *out = is_simple(&b.resolved) as i32;
        Ok(())
    })
}

/// Groups parallel edges and, with `drop_loops != 0`, removes
/// self-loops. `group[e]` (m values) receives the 0-based index of the
/// kept edge that edge `e` became, numbered in order of first
/// appearance, or -1 for a removed loop; the kept edge is the first row
/// of its group. `weight` (m values, or NULL) receives at every row the
/// combined weight of that row's group (`combine`: 0 = sum, 1 = mean,
/// 2 = min, 3 = max, 4 = first, 5 = count); without a weight role only
/// count writes it, and the rest leave it untouched. `kept` receives the
/// number of kept edges.
///
/// # Safety
/// `group` must hold m values, `weight` null or m values, `kept`
/// writable.
#[no_mangle]
pub unsafe extern "C" fn net_frame_simplify(
    h: *const FrameHandle,
    combine: i32,
    drop_loops: i32,
    group: *mut i64,
    weight: *mut f64,
    kept: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let b = built_ref(h)?;
        let s = simplify(&b.resolved, Combine::from_code(combine)?, drop_loops != 0);
        let m = s.group.len();
        let out = slice_mut(group, m);
        for (e, g) in s.group.iter().enumerate() {
            out[e] = g.map_or(-1, |g| g as i64);
        }
        if let (false, Some(w)) = (weight.is_null(), &s.weight) {
            let out = slice_mut(weight, m);
            for (e, g) in s.group.iter().enumerate() {
                if let Some(g) = g {
                    out[e] = w[*g as usize];
                }
            }
        }
        *kept = s.keep.len() as i64;
        Ok(())
    })
}

/// The derived properties, as `SUMMARY_LEN` doubles:
/// 0 nodes, 1 edges, 2 directed, 3 number of modes, 4 between_modes_only
/// (-1 undeclared), 5 weighted, 6 weight type (0 distance, 1 strength,
/// -1 none), 7 smallest weight, 8 largest weight, 9 number of relations,
/// 10 signed, 11 negative edges, 12 time kind (0 static, 1 contacts,
/// 2 edge spells, 3 waves), 13 earliest start, 14 latest end, 15 node
/// spells, 16 node geometry, 17 edge geometry, 18 coordinate kind (-1
/// none, 0 unknown, 1 lonlat, 2 planar), 19 loops, 20 parallel edges,
/// 21 isolates, 22 weak components, 23 density, 24 unobserved nodes,
/// 25 unobserved edges, 26 number of layers, 27 layer table given,
/// 28 intralayer edges, 29 coupling edges, 30 interlayer edges,
/// 31 couplings_only (-1 undeclared), 32 number of egos, 33 egos that
/// are nodes, 34 egos_are_nodes (-1 undeclared), 35 node spell rows
/// (table), 36 node change rows. Field 15 is 1 for spells from the node
/// roles and 2 for a node spells table. Field 1 counts edge rows; the
/// built graph has field 1 minus field 25 edges. Undefined values are
/// NaN. Builds the graph if needed.
///
/// # Safety
/// `out` must hold `SUMMARY_LEN` doubles.
#[no_mangle]
pub unsafe extern "C" fn net_frame_summary(
    h: *mut FrameHandle,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_mut(h)?;
        let s = f.summary(thread_count(threads))?;
        let flag = |b: bool| b as i32 as f64;
        let opt = |o: Option<bool>| o.map_or(-1.0, flag);
        let values = [
            s.nodes as f64,
            s.edges as f64,
            flag(s.directed),
            s.mode_labels.len() as f64,
            opt(s.between_modes_only),
            flag(s.weight_column.is_some()),
            match s.weight_type {
                Some(WeightType::Distance) => 0.0,
                Some(WeightType::Strength) => 1.0,
                None => -1.0,
            },
            s.weight_min,
            s.weight_max,
            s.relation_labels.len() as f64,
            flag(s.signed),
            s.negative_edges as f64,
            s.time.code() as f64,
            s.time_start,
            s.time_end,
            if s.node_spell_table {
                2.0
            } else {
                flag(s.node_spells)
            },
            flag(s.node_geometry),
            flag(s.edge_geometry),
            match s.coordinates {
                None => -1.0,
                Some(CoordinateKind::Unknown) => 0.0,
                Some(CoordinateKind::LonLat) => 1.0,
                Some(CoordinateKind::Planar) => 2.0,
            },
            s.loops as f64,
            s.parallel as f64,
            s.isolates as f64,
            s.components as f64,
            s.density,
            s.missing_nodes as f64,
            s.missing_edges as f64,
            s.layer_labels.len() as f64,
            flag(s.layer_table),
            s.intralayer as f64,
            s.coupling as f64,
            s.interlayer as f64,
            opt(s.couplings_only),
            s.ego_labels.len() as f64,
            s.egos_as_nodes as f64,
            opt(s.egos_are_nodes),
            s.node_spell_rows as f64,
            s.node_change_rows as f64,
        ];
        slice_mut(out, SUMMARY_LEN).copy_from_slice(&values);
        Ok(())
    })
}

/// Label `i` of the summary's modes (`which == 0`), relations
/// (`which == 1`), layers (`which == 4`), or egos (`which == 5`), the
/// weight column name (`which == 2`, `i == 0`), or the coordinate
/// reference system (`which == 3`, `i == 0`); NULL when out of range or
/// before `net_frame_summary`. Valid while the handle lives and the
/// frame is unchanged.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_frame_summary_label(
    h: *const FrameHandle,
    which: i32,
    i: i64,
) -> *const c_char {
    let f = match frame_ref(h) {
        Ok(f) => f,
        Err(_) => return std::ptr::null(),
    };
    let i = match usize::try_from(i) {
        Ok(i) => i,
        Err(_) => return std::ptr::null(),
    };
    usize::try_from(which)
        .ok()
        .and_then(|w| f.labels.get(w))
        .and_then(|v| v.get(i))
        .map_or(std::ptr::null(), |s| s.as_ptr())
}

/// The number of nodes in each mode (`which == 0`, as many values as
/// modes), edges in each relation (`which == 1`), intralayer edges in
/// each layer (`which == 2`), nodes in each ego's network (`which == 3`),
/// or edge rows in each ego's network (`which == 4`).
///
/// # Safety
/// `out` must hold as many values as the summary reports modes,
/// relations, layers, or egos.
#[no_mangle]
pub unsafe extern "C" fn net_frame_summary_counts(
    h: *const FrameHandle,
    which: i32,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let f = frame_ref(h)?;
        let s = f
            .summary
            .as_ref()
            .ok_or_else(|| Error::InvalidArgument("call net_frame_summary first".into()))?;
        let counts = match which {
            0 => &s.mode_counts,
            1 => &s.relation_counts,
            2 => &s.layer_counts,
            3 => &s.ego_node_counts,
            4 => &s.ego_edge_counts,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown counts selector {which}"
                )))
            }
        };
        for (o, &c) in slice_mut(out, counts.len()).iter_mut().zip(counts) {
            *o = c as i64;
        }
        Ok(())
    })
}
