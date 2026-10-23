//! C ABI for reading and writing network files.
//!
//! A network handle holds three tables, addressed by number: 0 = graph
//! attributes (one row), 1 = nodes (first column `id`), 2 = edges (first
//! columns `from` and `to`). Column kinds are 0 = logical, 1 = integer,
//! 2 = double, 3 = text.
//!
//! Reading: `net_read` fills a handle; callers then query each table's
//! columns. Writing: `net_network_new` creates an empty handle; callers add
//! columns and call `net_write`.

use std::ffi::{c_char, CStr, CString};
use std::path::Path;

use engine::io::{self, Column, Format, Kind, Network, ReadOptions, Table, Value, WriteOptions};
use engine::Error;

use crate::{deref, guard};

/// Opaque network handle.
pub struct NetworkHandle {
    directed: bool,
    tables: [Table; 3],
    names: [Vec<CString>; 3],
    strings: [Vec<Option<Vec<Option<CString>>>>; 3],
}

impl NetworkHandle {
    fn new(directed: bool, tables: [Table; 3]) -> NetworkHandle {
        let names = tables
            .clone()
            .map(|t| t.names.iter().map(|n| cstring(n)).collect());
        let strings = tables.clone().map(|t| vec![None; t.columns.len()]);
        NetworkHandle {
            directed,
            tables,
            names,
            strings,
        }
    }

    fn from_network(net: Network) -> NetworkHandle {
        let mut graph = Table::default();
        for (name, value) in &net.graph {
            let col = Column::from_values(vec![value.clone()], None)
                .unwrap_or(Column::Str(vec![value.to_text()]));
            // Duplicate graph attribute names keep the first value.
            let _ = graph.push(name.clone(), col);
        }
        NetworkHandle::new(net.directed, [graph, net.nodes, net.edges])
    }

    fn to_network(&self) -> Result<Network, Error> {
        let g = &self.tables[0];
        let graph = (0..g.columns.len())
            .filter(|_| g.nrow() > 0)
            .map(|i| (g.names[i].clone(), g.columns[i].get(0)))
            .collect();
        Ok(Network {
            directed: self.directed,
            graph,
            nodes: self.tables[1].clone(),
            edges: self.tables[2].clone(),
        })
    }
}

fn cstring(s: &str) -> CString {
    CString::new(s.replace('\0', "")).expect("NUL bytes removed")
}

unsafe fn opt_str<'a>(p: *const c_char) -> Option<&'a str> {
    if p.is_null() {
        None
    } else {
        CStr::from_ptr(p).to_str().ok()
    }
}

unsafe fn handle<'a>(h: *const NetworkHandle) -> Result<&'a NetworkHandle, Error> {
    deref(h, "network")
}

fn table_index(which: i32) -> Result<usize, Error> {
    match which {
        0..=2 => Ok(which as usize),
        _ => Err(Error::InvalidArgument(format!("no table {which}"))),
    }
}

fn kind_code(k: Kind) -> i32 {
    match k {
        Kind::Bool => 0,
        Kind::Int => 1,
        Kind::Float => 2,
        Kind::Str => 3,
    }
}

/// Reads a network file. `format` may be null to infer it from the
/// extension (for `.xml`, GraphML or DyNetML from the root element).
/// `directed` is -1 (not specified), 0, or 1; it only affects
/// formats that do not record direction.
///
/// # Safety
/// `path` must be a NUL-terminated string; `format` null or NUL-terminated;
/// `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_read(
    path: *const c_char,
    format: *const c_char,
    directed: i32,
    out: *mut *mut NetworkHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let path = opt_str(path).ok_or_else(|| Error::InvalidArgument("invalid path".into()))?;
        let format = match opt_str(format) {
            Some(name) => Some(
                Format::from_name(name)
                    .ok_or_else(|| Error::InvalidArgument(format!("unknown format {name:?}")))?,
            ),
            None => None,
        };
        let opts = ReadOptions {
            directed: if directed < 0 {
                None
            } else {
                Some(directed != 0)
            },
        };
        let net = io::read_file(Path::new(path), format, &opts)?;
        *out = Box::into_raw(Box::new(NetworkHandle::from_network(net)));
        Ok(())
    })
}

/// Creates an empty network for writing.
#[no_mangle]
pub extern "C" fn net_network_new(directed: i32) -> *mut NetworkHandle {
    let empty = [Table::default(), Table::default(), Table::default()];
    Box::into_raw(Box::new(NetworkHandle::new(directed != 0, empty)))
}

/// Releases a network handle. Null is ignored.
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_network_free(h: *mut NetworkHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// 1 if directed, 0 if not, -1 for a null handle.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_network_directed(h: *const NetworkHandle) -> i32 {
    handle(h).map(|n| n.directed as i32).unwrap_or(-1)
}

/// Number of rows in a table, or -1 on error.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_table_nrow(h: *const NetworkHandle, which: i32) -> i64 {
    match (handle(h), table_index(which)) {
        (Ok(n), Ok(t)) => n.tables[t].nrow() as i64,
        _ => -1,
    }
}

/// Number of columns in a table, or -1 on error.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_table_ncol(h: *const NetworkHandle, which: i32) -> i32 {
    match (handle(h), table_index(which)) {
        (Ok(n), Ok(t)) => n.tables[t].columns.len() as i32,
        _ => -1,
    }
}

/// Name of a column (valid while the handle lives), or null on error.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_table_name(
    h: *const NetworkHandle,
    which: i32,
    col: i32,
) -> *const c_char {
    match (handle(h), table_index(which)) {
        (Ok(n), Ok(t)) => n.names[t]
            .get(col as usize)
            .map(|c| c.as_ptr())
            .unwrap_or(std::ptr::null()),
        _ => std::ptr::null(),
    }
}

/// Kind of a column (0 logical, 1 integer, 2 double, 3 text), or -1.
///
/// # Safety
/// `h` must be null or a live handle.
#[no_mangle]
pub unsafe extern "C" fn net_table_kind(h: *const NetworkHandle, which: i32, col: i32) -> i32 {
    match (handle(h), table_index(which)) {
        (Ok(n), Ok(t)) => n.tables[t]
            .columns
            .get(col as usize)
            .map(|c| kind_code(c.kind()))
            .unwrap_or(-1),
        _ => -1,
    }
}

/// Copies a logical, integer, or double column into `out` as doubles, with
/// `na[i]` set to 1 for missing values. Integers beyond 2^53 lose precision.
///
/// # Safety
/// `h` must be a live handle; `out` and `na` must hold `nrow` elements.
#[no_mangle]
pub unsafe extern "C" fn net_table_numeric(
    h: *const NetworkHandle,
    which: i32,
    col: i32,
    out: *mut f64,
    na: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let n = handle(h)?;
        let t = &n.tables[table_index(which)?];
        let c = t
            .columns
            .get(col as usize)
            .ok_or_else(|| Error::InvalidArgument(format!("no column {col}")))?;
        let rows = t.nrow();
        let out = crate::slice_mut(out, rows);
        let na = crate::slice_mut(na, rows);
        for i in 0..rows {
            let v = match c.get(i) {
                Value::Bool(b) => Some(b as i64 as f64),
                Value::Int(x) => Some(x as f64),
                Value::Float(x) => Some(x),
                Value::Null => None,
                Value::Str(_) => return Err(Error::InvalidArgument("column is text".into())),
            };
            out[i] = v.unwrap_or(0.0);
            na[i] = v.is_none() as i32;
        }
        Ok(())
    })
}

/// Text value of a cell (valid while the handle lives), or null if the value
/// is missing or the column is not text.
///
/// # Safety
/// `h` must be null or a live handle created by this library.
#[no_mangle]
pub unsafe extern "C" fn net_table_string(
    h: *mut NetworkHandle,
    which: i32,
    col: i32,
    row: i64,
) -> *const c_char {
    let n = match h.as_mut() {
        Some(n) => n,
        None => return std::ptr::null(),
    };
    let t = match table_index(which) {
        Ok(t) => t,
        Err(_) => return std::ptr::null(),
    };
    let c = col as usize;
    let column = match n.tables[t].columns.get(c) {
        Some(Column::Str(v)) => v,
        _ => return std::ptr::null(),
    };
    let cache = n.strings[t][c]
        .get_or_insert_with(|| column.iter().map(|s| s.as_deref().map(cstring)).collect());
    match cache.get(row as usize) {
        Some(Some(s)) => s.as_ptr(),
        _ => std::ptr::null(),
    }
}

fn push_column(n: &mut NetworkHandle, which: i32, name: &str, column: Column) -> Result<(), Error> {
    let t = table_index(which)?;
    n.tables[t].push(name.to_string(), column)?;
    n.names[t].push(cstring(name));
    n.strings[t].push(None);
    Ok(())
}

/// Adds a logical (kind 0), integer (1), or double (2) column. `na` may be
/// null when nothing is missing.
///
/// # Safety
/// `h` must be a live handle; `name` NUL-terminated; `vals` (and `na` if
/// non-null) must hold `n` elements.
#[no_mangle]
pub unsafe extern "C" fn net_table_add_numeric(
    h: *mut NetworkHandle,
    which: i32,
    name: *const c_char,
    kind: i32,
    vals: *const f64,
    na: *const i32,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = h
            .as_mut()
            .ok_or_else(|| Error::InvalidArgument("network handle is null".into()))?;
        let name =
            opt_str(name).ok_or_else(|| Error::InvalidArgument("invalid column name".into()))?;
        let rows = n.max(0) as usize;
        let vals = crate::slice(vals, rows);
        let missing = |i: usize| !na.is_null() && *na.add(i) != 0;
        let column = match kind {
            0 => Column::Bool(
                (0..rows)
                    .map(|i| (!missing(i)).then(|| vals[i] != 0.0))
                    .collect(),
            ),
            1 => Column::Int(
                (0..rows)
                    .map(|i| (!missing(i)).then(|| vals[i] as i64))
                    .collect(),
            ),
            2 => Column::Float((0..rows).map(|i| (!missing(i)).then(|| vals[i])).collect()),
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "invalid column kind {kind}"
                )))
            }
        };
        push_column(net, which, name, column)
    })
}

/// Adds a text column. Null entries are missing values; invalid UTF-8 is
/// replaced.
///
/// # Safety
/// `h` must be a live handle; `name` NUL-terminated; `vals` must hold `n`
/// pointers, each null or NUL-terminated.
#[no_mangle]
pub unsafe extern "C" fn net_table_add_string(
    h: *mut NetworkHandle,
    which: i32,
    name: *const c_char,
    vals: *const *const c_char,
    n: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = h
            .as_mut()
            .ok_or_else(|| Error::InvalidArgument("network handle is null".into()))?;
        let name =
            opt_str(name).ok_or_else(|| Error::InvalidArgument("invalid column name".into()))?;
        let ptrs = crate::slice(vals, n.max(0) as usize);
        let column = Column::Str(
            ptrs.iter()
                .map(|&p| (!p.is_null()).then(|| CStr::from_ptr(p).to_string_lossy().into_owned()))
                .collect(),
        );
        push_column(net, which, name, column)
    })
}

/// Writes a network file. `format` may be null to infer it from the
/// extension; `weight` may be null to use the `weight` column.
///
/// # Safety
/// `h` must be a live handle; strings null or NUL-terminated.
#[no_mangle]
pub unsafe extern "C" fn net_write(
    h: *const NetworkHandle,
    path: *const c_char,
    format: *const c_char,
    weight: *const c_char,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let net = handle(h)?.to_network()?;
        let path = opt_str(path).ok_or_else(|| Error::InvalidArgument("invalid path".into()))?;
        let format = match opt_str(format) {
            Some(name) => Some(
                Format::from_name(name)
                    .ok_or_else(|| Error::InvalidArgument(format!("unknown format {name:?}")))?,
            ),
            None => None,
        };
        let opts = WriteOptions {
            weight: opt_str(weight).map(str::to_string),
        };
        io::write_file(Path::new(path), &net, format, &opts)
    })
}
