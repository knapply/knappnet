//! A small typed, columnar table used to exchange networks with files and
//! language bindings.

use std::collections::HashMap;

use crate::error::{Error, Result};

/// A single attribute value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Missing.
    Null,
    /// Boolean.
    Bool(bool),
    /// 64-bit integer.
    Int(i64),
    /// 64-bit float.
    Float(f64),
    /// Text.
    Str(String),
}

impl Value {
    /// Renders the value as text, or `None` when it is missing.
    pub fn to_text(&self) -> Option<String> {
        match self {
            Value::Null => None,
            Value::Bool(b) => Some(if *b { "true" } else { "false" }.to_string()),
            Value::Int(i) => Some(i.to_string()),
            Value::Float(f) => Some(format_float(*f)),
            Value::Str(s) => Some(s.clone()),
        }
    }
}

/// Formats a float so it always reads back as a float (with a decimal point
/// or exponent) and round-trips exactly.
pub fn format_float(f: f64) -> String {
    if f.is_nan() {
        return "NaN".into();
    }
    if f.is_infinite() {
        return if f > 0.0 { "Infinity" } else { "-Infinity" }.into();
    }
    let a = f.abs();
    let mut s = if a != 0.0 && !(1e-5..1e16).contains(&a) {
        format!("{f:e}")
    } else {
        format!("{f}")
    };
    if !s.contains(['.', 'e', 'E']) {
        s.push_str(".0");
    } else if let Some(pos) = s.find('e') {
        if !s[..pos].contains('.') {
            s.insert_str(pos, ".0");
        }
    }
    s
}

/// The type of a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Booleans.
    Bool,
    /// 64-bit integers.
    Int,
    /// 64-bit floats.
    Float,
    /// Text.
    Str,
}

/// A typed column with missing values.
#[derive(Debug, Clone, PartialEq)]
pub enum Column {
    /// Booleans.
    Bool(Vec<Option<bool>>),
    /// 64-bit integers.
    Int(Vec<Option<i64>>),
    /// 64-bit floats.
    Float(Vec<Option<f64>>),
    /// Text.
    Str(Vec<Option<String>>),
}

impl Column {
    /// Number of rows.
    pub fn len(&self) -> usize {
        match self {
            Column::Bool(v) => v.len(),
            Column::Int(v) => v.len(),
            Column::Float(v) => v.len(),
            Column::Str(v) => v.len(),
        }
    }

    /// Whether the column has no rows.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The column's type.
    pub fn kind(&self) -> Kind {
        match self {
            Column::Bool(_) => Kind::Bool,
            Column::Int(_) => Kind::Int,
            Column::Float(_) => Kind::Float,
            Column::Str(_) => Kind::Str,
        }
    }

    /// The value at row `i`.
    pub fn get(&self, i: usize) -> Value {
        let missing = Value::Null;
        match self {
            Column::Bool(v) => v[i].map(Value::Bool).unwrap_or(missing),
            Column::Int(v) => v[i].map(Value::Int).unwrap_or(missing),
            Column::Float(v) => v[i].map(Value::Float).unwrap_or(missing),
            Column::Str(v) => v[i].clone().map(Value::Str).unwrap_or(missing),
        }
    }

    /// Converts values to `kind`, or infers the narrowest kind that holds
    /// every value when `kind` is `None` (Bool, then Int, then Float, then
    /// Str). Text that does not parse as the requested kind is an error.
    pub fn from_values(values: Vec<Value>, kind: Option<Kind>) -> Result<Column> {
        let kind = kind.unwrap_or_else(|| infer_kind(&values));
        let bad =
            |v: &Value| Error::InvalidArgument(format!("value {v:?} cannot be stored as {kind:?}"));
        Ok(match kind {
            Kind::Bool => Column::Bool(
                values
                    .iter()
                    .map(|v| match v {
                        Value::Null => Ok(None),
                        Value::Bool(b) => Ok(Some(*b)),
                        Value::Int(i) if *i == 0 || *i == 1 => Ok(Some(*i == 1)),
                        Value::Str(s) => parse_bool(s).map(Some).ok_or_else(|| bad(v)),
                        _ => Err(bad(v)),
                    })
                    .collect::<Result<_>>()?,
            ),
            Kind::Int => Column::Int(
                values
                    .iter()
                    .map(|v| match v {
                        Value::Null => Ok(None),
                        Value::Int(i) => Ok(Some(*i)),
                        Value::Bool(b) => Ok(Some(*b as i64)),
                        Value::Float(f) if f.fract() == 0.0 && f.abs() < 9.0e15 => {
                            Ok(Some(*f as i64))
                        }
                        Value::Str(s) => s.trim().parse().map(Some).map_err(|_| bad(v)),
                        _ => Err(bad(v)),
                    })
                    .collect::<Result<_>>()?,
            ),
            Kind::Float => Column::Float(
                values
                    .iter()
                    .map(|v| match v {
                        Value::Null => Ok(None),
                        Value::Int(i) => Ok(Some(*i as f64)),
                        Value::Float(f) => Ok(Some(*f)),
                        Value::Bool(b) => Ok(Some(*b as i64 as f64)),
                        Value::Str(s) => parse_float(s).map(Some).ok_or_else(|| bad(v)),
                    })
                    .collect::<Result<_>>()?,
            ),
            Kind::Str => Column::Str(values.iter().map(Value::to_text).collect()),
        })
    }
}

fn infer_kind(values: &[Value]) -> Kind {
    let mut kind: Option<Kind> = None;
    for v in values {
        let k = match v {
            Value::Null => continue,
            Value::Bool(_) => Kind::Bool,
            Value::Int(_) => Kind::Int,
            Value::Float(_) => Kind::Float,
            Value::Str(_) => return Kind::Str,
        };
        kind = Some(match (kind, k) {
            (None, k) => k,
            (Some(a), b) if a == b => a,
            (Some(Kind::Int), Kind::Float) | (Some(Kind::Float), Kind::Int) => Kind::Float,
            _ => return Kind::Str,
        });
    }
    kind.unwrap_or(Kind::Str)
}

/// Parses common boolean spellings.
pub fn parse_bool(s: &str) -> Option<bool> {
    match s.trim().to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" | "t" => Some(true),
        "false" | "0" | "no" | "f" => Some(false),
        _ => None,
    }
}

/// Parses a float, accepting the spellings used by common network formats.
pub fn parse_float(s: &str) -> Option<f64> {
    let t = s.trim();
    match t.to_ascii_lowercase().as_str() {
        "nan" => Some(f64::NAN),
        "inf" | "infinity" | "+inf" | "+infinity" => Some(f64::INFINITY),
        "-inf" | "-infinity" => Some(f64::NEG_INFINITY),
        _ => t.parse().ok(),
    }
}

/// A table of named, equal-length columns.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table {
    /// Column names, unique.
    pub names: Vec<String>,
    /// Columns, aligned with `names`.
    pub columns: Vec<Column>,
}

impl Table {
    /// Number of rows (0 for a table with no columns).
    pub fn nrow(&self) -> usize {
        self.columns.first().map(Column::len).unwrap_or(0)
    }

    /// The column called `name`, if any.
    pub fn column(&self, name: &str) -> Option<&Column> {
        self.names
            .iter()
            .position(|n| n == name)
            .map(|i| &self.columns[i])
    }

    /// Adds a column, checking its length and the uniqueness of its name.
    pub fn push(&mut self, name: impl Into<String>, column: Column) -> Result<()> {
        let name = name.into();
        if self.names.contains(&name) {
            return Err(Error::InvalidArgument(format!(
                "duplicate column name {name:?}"
            )));
        }
        if !self.columns.is_empty() && column.len() != self.nrow() {
            return Err(Error::LengthMismatch {
                what: "column",
                expected: self.nrow(),
                found: column.len(),
            });
        }
        self.names.push(name);
        self.columns.push(column);
        Ok(())
    }

    /// String values of a column (the identifier and endpoint columns).
    pub fn strings(&self, name: &str) -> Result<Vec<String>> {
        let col = self
            .column(name)
            .ok_or_else(|| Error::InvalidArgument(format!("missing column {name:?}")))?;
        (0..col.len())
            .map(|i| {
                col.get(i).to_text().ok_or_else(|| {
                    Error::InvalidArgument(format!("column {name:?} has a missing value"))
                })
            })
            .collect()
    }
}

/// Accumulates rows whose columns appear over time, then builds a typed
/// table. Columns that first appear after some rows are back-filled with
/// missing values.
#[derive(Debug, Default)]
pub struct TableBuilder {
    names: Vec<String>,
    index: HashMap<String, usize>,
    values: Vec<Vec<Value>>,
    kinds: Vec<Option<Kind>>,
    rows: usize,
}

impl TableBuilder {
    /// Creates a builder whose first columns are `names` (typically the
    /// identifier columns), all text.
    pub fn with_columns(names: &[&str]) -> TableBuilder {
        let mut b = TableBuilder::default();
        for n in names {
            b.column(n, Some(Kind::Str));
        }
        b
    }

    /// Declares a column (if new) and returns its position. A declared
    /// kind is used when the table is built; otherwise the kind is inferred.
    pub fn column(&mut self, name: &str, kind: Option<Kind>) -> usize {
        if let Some(&i) = self.index.get(name) {
            if kind.is_some() {
                self.kinds[i] = kind;
            }
            return i;
        }
        let i = self.names.len();
        self.names.push(name.to_string());
        self.index.insert(name.to_string(), i);
        self.values.push(vec![Value::Null; self.rows]);
        self.kinds.push(kind);
        i
    }

    /// Appends a row. Unmentioned columns get missing values.
    pub fn push_row(&mut self, row: Vec<(String, Value)>) {
        for (name, value) in row {
            let i = self.column(&name, None);
            let col = &mut self.values[i];
            col.resize(self.rows, Value::Null);
            col.push(value);
        }
        self.rows += 1;
        for col in &mut self.values {
            col.resize(self.rows, Value::Null);
        }
    }

    /// Number of rows so far.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Builds the table.
    pub fn finish(self) -> Result<Table> {
        let mut t = Table::default();
        for ((name, values), kind) in self.names.into_iter().zip(self.values).zip(self.kinds) {
            let col = Column::from_values(values, kind)?;
            t.names.push(name);
            t.columns.push(col);
        }
        Ok(t)
    }
}

/// A network as two tables plus graph-level attributes.
///
/// `nodes` starts with a text column `id`; `edges` starts with text columns
/// `from` and `to` holding node identifiers. All other columns are
/// attributes.
#[derive(Debug, Clone, PartialEq)]
pub struct Network {
    /// Whether edges are directed.
    pub directed: bool,
    /// Graph-level attributes.
    pub graph: Vec<(String, Value)>,
    /// One row per node.
    pub nodes: Table,
    /// One row per edge.
    pub edges: Table,
}

impl Network {
    /// Checks the structural columns and that every endpoint is a node.
    pub fn validate(&self) -> Result<()> {
        let ids = self.nodes.strings("id")?;
        let mut seen = std::collections::HashSet::with_capacity(ids.len());
        for id in &ids {
            if !seen.insert(id.as_str()) {
                return Err(Error::InvalidArgument(format!(
                    "node identifier {id:?} appears more than once"
                )));
            }
        }
        for col in ["from", "to"] {
            for v in self.edges.strings(col)? {
                if !seen.contains(v.as_str()) {
                    return Err(Error::InvalidArgument(format!(
                        "edge endpoint {v:?} is not a node"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Attribute columns of the node table (everything except `id`).
    pub fn node_attributes(&self) -> impl Iterator<Item = (&String, &Column)> {
        self.nodes
            .names
            .iter()
            .zip(&self.nodes.columns)
            .filter(|(n, _)| n.as_str() != "id")
    }

    /// Attribute columns of the edge table (everything except `from`, `to`).
    pub fn edge_attributes(&self) -> impl Iterator<Item = (&String, &Column)> {
        self.edges
            .names
            .iter()
            .zip(&self.edges.columns)
            .filter(|(n, _)| n.as_str() != "from" && n.as_str() != "to")
    }

    /// Uses a text column as node identifiers when every node has a value
    /// and the values are unique: the column replaces `id` and edge
    /// endpoints are translated. Returns whether the relabeling happened.
    pub fn use_as_ids(&mut self, name: &str) -> Result<bool> {
        let pos = match self.nodes.names.iter().position(|n| n == name) {
            Some(p) => p,
            None => return Ok(false),
        };
        let labels: Option<Vec<String>> = match &self.nodes.columns[pos] {
            Column::Str(v) => v.iter().cloned().collect(),
            _ => None,
        };
        let labels = match labels {
            Some(l) => l,
            None => return Ok(false),
        };
        let unique: std::collections::HashSet<&String> = labels.iter().collect();
        if unique.len() != labels.len() {
            return Ok(false);
        }
        let old = self.node_index()?;
        let id_pos = self
            .nodes
            .names
            .iter()
            .position(|n| n == "id")
            .expect("id column");
        for end in ["from", "to"] {
            let e = self
                .edges
                .names
                .iter()
                .position(|n| n == end)
                .expect("endpoint column");
            let mapped = self
                .edges
                .strings(end)?
                .into_iter()
                .map(|v| Some(labels[old[&v]].clone()))
                .collect();
            self.edges.columns[e] = Column::Str(mapped);
        }
        self.nodes.columns[id_pos] = Column::Str(labels.into_iter().map(Some).collect());
        self.nodes.names.remove(pos);
        self.nodes.columns.remove(pos);
        Ok(true)
    }

    /// Node identifier to row index.
    pub fn node_index(&self) -> Result<HashMap<String, usize>> {
        Ok(self
            .nodes
            .strings("id")?
            .into_iter()
            .enumerate()
            .map(|(i, s)| (s, i))
            .collect())
    }
}

/// Builds a [`Network`] while reading a file: adds nodes on first mention
/// and keeps attribute names clear of the structural columns.
#[derive(Debug)]
pub struct NetworkBuilder {
    /// Whether edges are directed.
    pub directed: bool,
    /// Graph-level attributes.
    pub graph: Vec<(String, Value)>,
    nodes: TableBuilder,
    edges: TableBuilder,
    node_rows: HashMap<String, usize>,
}

impl NetworkBuilder {
    /// Creates an empty builder.
    pub fn new(directed: bool) -> NetworkBuilder {
        NetworkBuilder {
            directed,
            graph: Vec::new(),
            nodes: TableBuilder::with_columns(&["id"]),
            edges: TableBuilder::with_columns(&["from", "to"]),
            node_rows: HashMap::new(),
        }
    }

    /// Declares a node attribute's type.
    pub fn node_kind(&mut self, name: &str, kind: Kind) {
        self.nodes.column(&attr_name(name, &["id"]), Some(kind));
    }

    /// Declares an edge attribute's type.
    pub fn edge_kind(&mut self, name: &str, kind: Kind) {
        self.edges
            .column(&attr_name(name, &["from", "to"]), Some(kind));
    }

    /// Adds a node. Returns false if the identifier already exists.
    pub fn add_node(&mut self, id: &str, attrs: Vec<(String, Value)>) -> bool {
        if self.node_rows.contains_key(id) {
            return false;
        }
        let mut row = vec![("id".to_string(), Value::Str(id.to_string()))];
        row.extend(attrs.into_iter().map(|(k, v)| (attr_name(&k, &["id"]), v)));
        self.node_rows.insert(id.to_string(), self.nodes.rows());
        self.nodes.push_row(row);
        true
    }

    /// Whether a node exists.
    pub fn has_node(&self, id: &str) -> bool {
        self.node_rows.contains_key(id)
    }

    /// Adds an edge, creating missing endpoint nodes.
    pub fn add_edge(&mut self, from: &str, to: &str, attrs: Vec<(String, Value)>) {
        for id in [from, to] {
            if !self.has_node(id) {
                self.add_node(id, Vec::new());
            }
        }
        let mut row = vec![
            ("from".to_string(), Value::Str(from.to_string())),
            ("to".to_string(), Value::Str(to.to_string())),
        ];
        row.extend(
            attrs
                .into_iter()
                .map(|(k, v)| (attr_name(&k, &["from", "to"]), v)),
        );
        self.edges.push_row(row);
    }

    /// Builds the network.
    pub fn finish(self) -> Result<Network> {
        Ok(Network {
            directed: self.directed,
            graph: self.graph,
            nodes: self.nodes.finish()?,
            edges: self.edges.finish()?,
        })
    }
}

/// Renames attributes that would collide with structural columns.
fn attr_name(name: &str, reserved: &[&str]) -> String {
    if reserved.contains(&name) {
        format!("{name}_attr")
    } else {
        name.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floats_round_trip_and_look_like_floats() {
        for f in [
            0.0,
            1.0,
            -2.5,
            1e-7,
            3.0e20,
            0.1 + 0.2,
            f64::MAX,
            f64::MIN_POSITIVE,
        ] {
            let s = format_float(f);
            assert!(s.contains(['.', 'e']), "{s}");
            assert_eq!(s.parse::<f64>().unwrap(), f, "{s}");
        }
        assert_eq!(format_float(1.0), "1.0");
        assert_eq!(format_float(1e300), "1.0e300");
    }

    #[test]
    fn kinds_are_inferred() {
        let c =
            Column::from_values(vec![Value::Int(1), Value::Null, Value::Float(2.5)], None).unwrap();
        assert_eq!(c, Column::Float(vec![Some(1.0), None, Some(2.5)]));
        let c = Column::from_values(vec![Value::Int(1), Value::Str("x".into())], None).unwrap();
        assert_eq!(c.kind(), Kind::Str);
        let c = Column::from_values(vec![Value::Str("true".into())], Some(Kind::Bool)).unwrap();
        assert_eq!(c, Column::Bool(vec![Some(true)]));
    }

    #[test]
    fn builder_backfills_and_reserves_names() {
        let mut b = NetworkBuilder::new(true);
        b.add_node("a", vec![]);
        b.add_node("b", vec![("color".into(), Value::Str("red".into()))]);
        b.add_edge("a", "c", vec![("from".into(), Value::Int(1))]);
        let n = b.finish().unwrap();
        assert_eq!(n.nodes.names, vec!["id", "color"]);
        assert_eq!(n.nodes.nrow(), 3);
        assert_eq!(
            n.nodes.column("color").unwrap(),
            &Column::Str(vec![None, Some("red".into()), None])
        );
        assert_eq!(n.edges.names, vec!["from", "to", "from_attr"]);
        n.validate().unwrap();
    }
}
