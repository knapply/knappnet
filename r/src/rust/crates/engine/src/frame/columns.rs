//! Typed columns with missing values, and the two tables of a frame.
//!
//! Columns are stored the way bindings hold them: integers and logicals
//! as a value array plus an optional missing mask, doubles with NaN as the
//! missing value (R's `NA_real_` is a NaN, and every check here treats NaN
//! and missing alike), and text as one byte buffer with offsets, so a
//! column of ten million identifiers is two allocations rather than ten
//! million. Nothing here is specific to networks; [`super::validate()`] gives
//! the columns their meaning through roles.

use std::fmt;

use crate::error::{Error, Result};

/// The type of a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// 64-bit integers.
    Int,
    /// 64-bit floating point numbers.
    Float,
    /// UTF-8 text.
    Str,
    /// Logicals.
    Bool,
}

impl Kind {
    /// The name used in messages.
    pub fn name(self) -> &'static str {
        match self {
            Kind::Int => "integer",
            Kind::Float => "double",
            Kind::Str => "text",
            Kind::Bool => "logical",
        }
    }
}

/// A missing-value mask: `None` means no value is missing.
pub type Mask = Option<Vec<bool>>;

fn mask_get(mask: &Mask, i: usize) -> bool {
    mask.as_ref().is_some_and(|m| m[i])
}

fn check_mask(mask: &Mask, len: usize, what: &str) -> Result<()> {
    match mask {
        Some(m) if m.len() != len => Err(Error::LengthMismatch {
            what: if what == "nodes" {
                "node column mask"
            } else {
                "edge column mask"
            },
            expected: len,
            found: m.len(),
        }),
        _ => Ok(()),
    }
}

/// Text stored as one byte buffer with `n + 1` offsets.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Strings {
    offsets: Vec<usize>,
    bytes: Vec<u8>,
    missing: Mask,
}

impl Strings {
    /// Builds text storage from offsets (`n + 1` values, non-decreasing,
    /// ending at `bytes.len()`) and a missing mask.
    pub fn new(offsets: Vec<usize>, bytes: Vec<u8>, missing: Mask) -> Result<Strings> {
        let bad = |msg: &str| Error::InvalidArgument(format!("text column: {msg}"));
        if offsets.is_empty() {
            return Err(bad("offsets must hold n + 1 values"));
        }
        if offsets[0] != 0 || offsets[offsets.len() - 1] != bytes.len() {
            return Err(bad("offsets must start at 0 and end at the byte count"));
        }
        if offsets.windows(2).any(|w| w[1] < w[0]) {
            return Err(bad("offsets must not decrease"));
        }
        if std::str::from_utf8(&bytes).is_err() {
            return Err(bad("values must be valid UTF-8"));
        }
        let n = offsets.len() - 1;
        if let Some(m) = &missing {
            if m.len() != n {
                return Err(bad("the missing mask must have one entry per value"));
            }
        }
        Ok(Strings {
            offsets,
            bytes,
            missing,
        })
    }

    /// Builds text storage from string slices (`None` is missing).
    pub fn from_strs<S: AsRef<str>>(values: &[Option<S>]) -> Strings {
        let mut offsets = Vec::with_capacity(values.len() + 1);
        let mut bytes = Vec::new();
        let mut missing = vec![false; values.len()];
        let mut any_missing = false;
        offsets.push(0);
        for (i, v) in values.iter().enumerate() {
            match v {
                Some(s) => bytes.extend_from_slice(s.as_ref().as_bytes()),
                None => {
                    missing[i] = true;
                    any_missing = true;
                }
            }
            offsets.push(bytes.len());
        }
        Strings {
            offsets,
            bytes,
            missing: any_missing.then_some(missing),
        }
    }

    /// Number of values.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Whether there are no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The bytes of value `i`, or `None` when it is missing.
    #[inline]
    pub fn get(&self, i: usize) -> Option<&[u8]> {
        if mask_get(&self.missing, i) {
            None
        } else {
            Some(&self.bytes[self.offsets[i]..self.offsets[i + 1]])
        }
    }

    /// Value `i` as text, or `None` when it is missing.
    pub fn get_str(&self, i: usize) -> Option<&str> {
        // The constructor checked the whole buffer; offsets may still split
        // a multi-byte character, in which case the value is reported as
        // its lossy rendering rather than a panic.
        self.get(i).map(|b| std::str::from_utf8(b).unwrap_or(""))
    }

    /// The offsets (`n + 1` values).
    pub fn offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// The byte buffer.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

/// A typed column with missing values.
#[derive(Debug, Clone, PartialEq)]
pub enum Column {
    /// Integers with a missing mask.
    Int(Vec<i64>, Mask),
    /// Doubles; NaN is missing.
    Float(Vec<f64>),
    /// Text.
    Str(Strings),
    /// Logicals with a missing mask.
    Bool(Vec<bool>, Mask),
}

impl Column {
    /// Number of rows.
    pub fn len(&self) -> usize {
        match self {
            Column::Int(v, _) => v.len(),
            Column::Float(v) => v.len(),
            Column::Str(s) => s.len(),
            Column::Bool(v, _) => v.len(),
        }
    }

    /// Whether the column has no rows.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The column's type.
    pub fn kind(&self) -> Kind {
        match self {
            Column::Int(..) => Kind::Int,
            Column::Float(_) => Kind::Float,
            Column::Str(_) => Kind::Str,
            Column::Bool(..) => Kind::Bool,
        }
    }

    /// Whether row `i` is missing (NaN counts as missing in a double
    /// column).
    #[inline]
    pub fn is_missing(&self, i: usize) -> bool {
        match self {
            Column::Int(_, m) => mask_get(m, i),
            Column::Float(v) => v[i].is_nan(),
            Column::Str(s) => s.get(i).is_none(),
            Column::Bool(_, m) => mask_get(m, i),
        }
    }

    /// Row `i` rendered for a message, or `NA` when missing.
    pub fn render(&self, i: usize) -> String {
        if self.is_missing(i) {
            return "NA".into();
        }
        match self {
            Column::Int(v, _) => v[i].to_string(),
            Column::Float(v) => format!("{}", v[i]),
            Column::Str(s) => s.get_str(i).unwrap_or("").to_string(),
            Column::Bool(v, _) => if v[i] { "TRUE" } else { "FALSE" }.to_string(),
        }
    }

    /// Checks that a mask, if any, has one entry per value.
    fn check(&self, table: &str) -> Result<()> {
        match self {
            Column::Int(v, m) => check_mask(m, v.len(), table),
            Column::Bool(v, m) => check_mask(m, v.len(), table),
            _ => Ok(()),
        }
    }

    /// The rows `rows` of this column as a new column.
    pub fn take(&self, rows: &[usize]) -> Column {
        match self {
            Column::Int(v, m) => Column::Int(
                rows.iter().map(|&i| v[i]).collect(),
                m.as_ref().map(|m| rows.iter().map(|&i| m[i]).collect()),
            ),
            Column::Float(v) => Column::Float(rows.iter().map(|&i| v[i]).collect()),
            Column::Str(s) => {
                let values: Vec<Option<&str>> = rows.iter().map(|&i| s.get_str(i)).collect();
                Column::Str(Strings::from_strs(&values))
            }
            Column::Bool(v, m) => Column::Bool(
                rows.iter().map(|&i| v[i]).collect(),
                m.as_ref().map(|m| rows.iter().map(|&i| m[i]).collect()),
            ),
        }
    }
}

/// One of the two tables of a frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Which {
    /// The node table.
    Nodes,
    /// The edge table.
    Edges,
}

impl Which {
    /// The table's name in messages.
    pub fn name(self) -> &'static str {
        match self {
            Which::Nodes => "nodes",
            Which::Edges => "edges",
        }
    }
}

impl fmt::Display for Which {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// A table: named columns of equal length. Only the columns a role refers
/// to need to be present; attributes stay with the binding.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table {
    nrow: Option<usize>,
    names: Vec<String>,
    columns: Vec<Column>,
}

impl Table {
    /// A table with no columns and no fixed row count.
    pub fn new() -> Table {
        Table::default()
    }

    /// Fixes the row count before any column is added, so that a table
    /// with no role columns (a node table given only for its row count)
    /// still has a size.
    pub fn with_rows(nrow: usize) -> Table {
        Table {
            nrow: Some(nrow),
            ..Table::default()
        }
    }

    /// Number of rows (0 for a table without columns or a fixed count).
    pub fn nrow(&self) -> usize {
        self.nrow.unwrap_or(0)
    }

    /// Whether the row count has been fixed by a column or `with_rows`.
    pub fn has_rows(&self) -> bool {
        self.nrow.is_some()
    }

    /// Column names in order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Number of columns.
    pub fn ncol(&self) -> usize {
        self.columns.len()
    }

    /// Adds a column. Its length must match the table's rows, and its name
    /// must be new.
    pub fn push(&mut self, which: Which, name: &str, column: Column) -> Result<usize> {
        column.check(which.name())?;
        if self.names.iter().any(|n| n == name) {
            return Err(Error::InvalidArgument(format!(
                "`{which}` already has a column named `{name}`"
            )));
        }
        match self.nrow {
            Some(n) if n != column.len() => {
                return Err(Error::LengthMismatch {
                    what: match which {
                        Which::Nodes => "node column",
                        Which::Edges => "edge column",
                    },
                    expected: n,
                    found: column.len(),
                })
            }
            Some(_) => {}
            None => self.nrow = Some(column.len()),
        }
        self.names.push(name.to_string());
        self.columns.push(column);
        Ok(self.columns.len() - 1)
    }

    /// Position of the column called `name`.
    pub fn position(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// The column at position `i`.
    pub fn column(&self, i: usize) -> &Column {
        &self.columns[i]
    }

    /// The column called `name`.
    pub fn get(&self, name: &str) -> Option<&Column> {
        self.position(name).map(|i| &self.columns[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strings_round_trip() {
        let s = Strings::from_strs(&[Some("ab"), None, Some(""), Some("é")]);
        assert_eq!(s.len(), 4);
        assert_eq!(s.get_str(0), Some("ab"));
        assert_eq!(s.get(1), None);
        assert_eq!(s.get(2), Some(&b""[..]));
        assert_eq!(s.get_str(3), Some("é"));
        let rebuilt = Strings::new(s.offsets().to_vec(), s.bytes().to_vec(), None).unwrap();
        assert_eq!(rebuilt.get_str(1), Some(""));
        assert!(Strings::new(vec![0, 3], vec![1, 2], None).is_err());
        assert!(Strings::new(vec![0, 2, 1], vec![1, 2], None).is_err());
        assert!(Strings::new(vec![0, 1], vec![0xff], None).is_err());
        assert!(Strings::new(vec![], vec![], None).is_err());
    }

    #[test]
    fn table_checks_lengths_and_names() {
        let mut t = Table::new();
        t.push(Which::Edges, "a", Column::Float(vec![1.0, 2.0]))
            .unwrap();
        assert_eq!(t.nrow(), 2);
        assert!(t
            .push(Which::Edges, "b", Column::Int(vec![1], None))
            .is_err());
        assert!(t
            .push(Which::Edges, "a", Column::Int(vec![1, 2], None))
            .is_err());
        assert!(t
            .push(Which::Edges, "c", Column::Int(vec![1, 2], Some(vec![true])))
            .is_err());
        let fixed = Table::with_rows(3);
        assert_eq!(fixed.nrow(), 3);
        assert!(fixed.has_rows());
    }

    #[test]
    fn missing_values_and_rendering() {
        let c = Column::Float(vec![1.5, f64::NAN]);
        assert!(!c.is_missing(0) && c.is_missing(1));
        assert_eq!(c.render(1), "NA");
        let b = Column::Bool(vec![true, false], Some(vec![false, true]));
        assert_eq!(b.render(0), "TRUE");
        assert_eq!(b.render(1), "NA");
        let i = Column::Int(vec![7, 8], None);
        assert_eq!(i.take(&[1]).render(0), "8");
        let s = Column::Str(Strings::from_strs(&[Some("x"), None]));
        assert_eq!(s.take(&[1, 0]).render(1), "x");
    }
}
