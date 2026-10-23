//! Structured diagnostics: what validation found, with stable codes that
//! bindings turn into their own messages.
//!
//! Every diagnostic has a [`Code`] (a stable string such as
//! `E_UNKNOWN_NODE`), a severity (errors stop construction; notes are
//! advisory and can be silenced), the table it concerns, a message in
//! plain English that a binding may show as is, the first
//! [`ROW_CAP`] offending 0-based row indices, and the total count. The
//! message never lists rows, so bindings phrase them idiomatically
//! ("rows 3, 7, and 12 (and 4 more) of `edges`").

use std::fmt;

use super::columns::Which;

/// How many offending rows a diagnostic keeps. The count is always exact.
pub const ROW_CAP: usize = 50;

/// Whether a diagnostic stops construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// The frame cannot become a network until this is fixed.
    Error,
    /// Advisory; the frame is still valid.
    Note,
}

macro_rules! codes {
    ($( $variant:ident => $severity:ident, $text:literal; )*) => {
        /// The stable code of a diagnostic. See each variant's
        /// documentation for the condition it reports. Variants are named
        /// exactly as the codes bindings see, so the two never drift.
        #[allow(non_camel_case_types)]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub enum Code {
            $(
                #[doc = $text]
                $variant,
            )*
        }

        impl Code {
            /// Every code, in a fixed order.
            pub const ALL: &'static [Code] = &[$(Code::$variant),*];

            /// The stable string form, `E_...` for errors and `N_...` for
            /// notes.
            pub fn as_str(self) -> &'static str {
                match self {
                    $(Code::$variant => stringify!($variant),)*
                }
            }

            /// The condition the code reports, as a message template.
            pub fn template(self) -> &'static str {
                match self {
                    $(Code::$variant => $text,)*
                }
            }

            /// Whether the code is an error or a note.
            pub fn severity(self) -> Severity {
                match self {
                    $(Code::$variant => Severity::$severity,)*
                }
            }
        }
    };
}

codes! {
    // Structure.
    E_ENDPOINTS_UNDECLARED => Error, "the edge table has no source and target roles";
    E_ID_UNDECLARED => Error, "the node table has no identifier role";
    E_TABLE_ROLE_UNDECLARED => Error, "a layers, node_spells, or node_changes table lacks one of its required roles";
    E_ENDPOINT_TYPE => Error, "the source and target columns must be integer, double, or text and of the same type";
    E_ID_TYPE => Error, "the node identifier column must have the same type as the edge endpoints";
    E_MISSING_ENDPOINT => Error, "edge endpoints must not be missing";
    E_EMPTY_ID => Error, "identifiers must not be empty strings";
    E_MISSING_ID => Error, "node identifiers must not be missing";
    E_DUPLICATE_ID => Error, "node identifiers must be unique";
    E_UNKNOWN_NODE => Error, "edge endpoints must be identifiers in the node table";
    // Declarations.
    E_LOOP_DECLARED => Error, "self-loops are present but were declared absent";
    E_PARALLEL_DECLARED => Error, "parallel edges are present but were declared absent";
    E_MODE_WITHIN => Error, "edges within a mode are present but edges were declared to join different modes only";
    E_MODE_UNDECLARED => Error, "edges were declared to join different modes, but no mode role is set";
    E_INTERLAYER_DECLARED => Error, "interlayer edges between different actors are present but were declared absent (couplings_only)";
    E_EGO_UNKNOWN => Error, "egos were declared to be nodes, but some ego identifiers are not node identifiers";
    E_DECLARATION_CONFLICT => Error, "the declarations contradict each other";
    // Weights.
    E_WEIGHT_TYPE => Error, "the weight column must be numeric";
    E_WEIGHT_TYPE_MISSING => Error, "a weighted network must declare whether weights are strengths or distances";
    E_WEIGHT_NOT_FINITE => Error, "weights must not be missing, NaN, or infinite";
    E_DISTANCE_NOT_POSITIVE => Error, "weights declared as distances must be positive";
    E_WEIGHT_NEGATIVE => Error, "weights of a signed network are magnitudes and must not be negative";
    // Modes, relations, signs.
    E_MODE_TYPE => Error, "the mode column must be discrete (text, integer, or logical)";
    E_MODE_MISSING => Error, "every node must have a mode";
    E_RELATION_TYPE => Error, "the relation column must be discrete (text, integer, or logical)";
    E_RELATION_MISSING => Error, "every edge must belong to a relation";
    E_SIGN_TYPE => Error, "the sign column must be integer or logical";
    E_SIGN_VALUE => Error, "signs must be -1 or 1 (or FALSE and TRUE) and not missing";
    // Layers.
    E_LAYER_ROLE_PAIR => Error, "the source_layer and target_layer roles must be set together";
    E_LAYER_TYPE => Error, "layer columns must be discrete (text, integer, or logical) and of one type";
    E_LAYER_MISSING => Error, "every edge needs a layer at each endpoint";
    E_LAYER_UNKNOWN => Error, "edge layers must be identifiers in the layer table";
    E_LAYER_ID_MISSING => Error, "layer identifiers must not be missing";
    E_LAYER_ID_DUPLICATE => Error, "layer identifiers must be unique";
    // Egos.
    E_EGO_TYPE => Error, "ego columns must be text or integer and of one type";
    E_EGO_MISSING => Error, "every node and edge needs an ego";
    E_EGO_SELF => Error, "a node that is an ego must belong to its own ego network";
    E_EGO_MISMATCH => Error, "an edge must join two nodes of one ego, and carry that ego";
    // Time.
    E_TIME_TYPE => Error, "time columns must be numeric";
    E_TIME_MISSING => Error, "edge times must not be missing";
    E_TIME_NOT_FINITE => Error, "edge times must be finite";
    E_SPELL_ORDER => Error, "an edge must not end before it starts";
    E_NODE_SPELL_ORDER => Error, "a node must not leave before it arrives";
    E_NODE_SPELL_OVERLAP => Error, "the presence spells of one node must not overlap";
    E_SPELL_OUTSIDE_NODE => Error, "an edge can only be active while both of its endpoints are present";
    E_WAVE_TYPE => Error, "the wave column must be integer";
    E_WAVE_MISSING => Error, "every edge must belong to a wave";
    // Space.
    E_COORD_TYPE => Error, "coordinate columns must be numeric";
    E_COORD_ROLE_PAIR => Error, "the x and y roles must be set together";
    E_COORD_NOT_FINITE => Error, "coordinates must not be missing or infinite";
    E_LATITUDE_RANGE => Error, "latitudes must lie between -90 and 90";
    E_CRS_MISMATCH => Error, "the node and edge geometries must share one coordinate reference system";
    // Missing-data roles.
    E_MISSING_ROLE_TYPE => Error, "missing-data roles must be logical columns without missing values";
    E_MISSING_CONFLICT => Error, "a pair cannot be both an observed tie and an unobserved pair";
    // Notes.
    N_NODES_ADDED => Note, "nodes were added from edge endpoints absent from the node table";
    N_ISOLATES => Note, "some nodes have no edges";
    N_ZERO_BASED_IDS => Note, "the integer identifiers are 0 to n - 1; check that they are not positions in a 1-based table";
    N_IDS_COERCED => Note, "identifiers of integer and double type were matched as numbers";
    N_LOOPS => Note, "the network has self-loops; each measure's documentation says how they count";
    N_PARALLEL => Note, "the network has parallel edges (same endpoints, relation, and time); each measure's documentation says how they count";
    N_MISSING_NODES => Note, "some nodes are marked as unobserved";
    N_MISSING_EDGES => Note, "some edges are marked as unobserved";
}

impl fmt::Display for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl Code {
    /// Looks a code up by its string form.
    pub fn parse(s: &str) -> Option<Code> {
        Code::ALL.iter().copied().find(|c| c.as_str() == s)
    }
}

/// One finding of validation.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    /// The stable code.
    pub code: Code,
    /// The table whose rows are listed; `None` for the network as a whole.
    pub table: Option<Which>,
    /// A message in plain English, with specifics (column names, values)
    /// but no row list.
    pub message: String,
    /// The first [`ROW_CAP`] offending 0-based rows, ascending.
    pub rows: Vec<usize>,
    /// How many rows are affected in all (at least `rows.len()`).
    pub count: usize,
}

impl Diagnostic {
    /// A diagnostic about the network as a whole.
    pub fn network(code: Code, message: impl Into<String>) -> Diagnostic {
        Diagnostic {
            code,
            table: None,
            message: message.into(),
            rows: Vec::new(),
            count: 0,
        }
    }

    /// A diagnostic about rows of a table. `rows` holds the affected rows
    /// (any number); only the first [`ROW_CAP`] are kept.
    pub fn rows(
        code: Code,
        table: Which,
        message: impl Into<String>,
        rows: &[usize],
    ) -> Diagnostic {
        Diagnostic {
            code,
            table: Some(table),
            message: message.into(),
            rows: rows.iter().copied().take(ROW_CAP).collect(),
            count: rows.len(),
        }
    }

    /// A diagnostic about rows of a table whose list was already capped
    /// while scanning.
    pub fn capped(
        code: Code,
        table: Which,
        message: impl Into<String>,
        rows: Vec<usize>,
        count: usize,
    ) -> Diagnostic {
        debug_assert!(rows.len() <= ROW_CAP && rows.len() <= count);
        Diagnostic {
            code,
            table: Some(table),
            message: message.into(),
            rows,
            count,
        }
    }

    /// The severity, taken from the code.
    pub fn severity(&self) -> Severity {
        self.code.severity()
    }
}

/// Everything validation found.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Report {
    /// Diagnostics in the order they were found: structural errors first,
    /// then declarations, attributes, and notes.
    pub diagnostics: Vec<Diagnostic>,
}

impl Report {
    /// Adds a diagnostic.
    pub fn push(&mut self, d: Diagnostic) {
        self.diagnostics.push(d);
    }

    /// Whether no error was found.
    pub fn is_valid(&self) -> bool {
        self.errors().next().is_none()
    }

    /// The errors.
    pub fn errors(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics
            .iter()
            .filter(|d| d.severity() == Severity::Error)
    }

    /// The notes.
    pub fn notes(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics
            .iter()
            .filter(|d| d.severity() == Severity::Note)
    }

    /// The first diagnostic with `code`.
    pub fn find(&self, code: Code) -> Option<&Diagnostic> {
        self.diagnostics.iter().find(|d| d.code == code)
    }

    /// Whether a diagnostic with `code` was found.
    pub fn has(&self, code: Code) -> bool {
        self.find(code).is_some()
    }

    /// A one-paragraph rendering of the errors, for callers that want a
    /// single message (the C interface's error string, for instance).
    pub fn error_text(&self) -> String {
        let mut out = String::new();
        for d in self.errors() {
            if !out.is_empty() {
                out.push_str("; ");
            }
            out.push_str(&d.message);
            if let (Some(table), true) = (d.table, d.count > 0) {
                out.push_str(&format!(" ({} of `{table}`)", rows_text(&d.rows, d.count)));
            }
        }
        out
    }
}

/// "row 3", "rows 3 and 7", or "rows 3, 7, 9, 12, 15, ... (40 in all)",
/// with 1-based positions.
pub fn rows_text(rows: &[usize], count: usize) -> String {
    let shown: Vec<String> = rows.iter().take(5).map(|r| (r + 1).to_string()).collect();
    match count {
        0 => String::new(),
        1 => format!("row {}", shown[0]),
        2 if shown.len() == 2 => format!("rows {} and {}", shown[0], shown[1]),
        k if k <= 5 && shown.len() == k => {
            format!("rows {}, and {}", shown[..k - 1].join(", "), shown[k - 1])
        }
        k => format!("rows {}, ... ({k} in all)", shown.join(", ")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_have_stable_strings_and_severities() {
        assert_eq!(Code::E_UNKNOWN_NODE.as_str(), "E_UNKNOWN_NODE");
        assert_eq!(Code::parse("N_LOOPS"), Some(Code::N_LOOPS));
        assert_eq!(Code::parse("nope"), None);
        for c in Code::ALL {
            let s = c.as_str();
            assert!(s.starts_with("E_") || s.starts_with("N_"), "{s}");
            assert_eq!(
                c.severity(),
                if s.starts_with("E_") {
                    Severity::Error
                } else {
                    Severity::Note
                }
            );
            assert!(!c.template().is_empty());
            assert_eq!(Code::parse(s), Some(*c));
        }
        let mut sorted: Vec<&str> = Code::ALL.iter().map(|c| c.as_str()).collect();
        sorted.dedup();
        assert_eq!(sorted.len(), Code::ALL.len());
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), Code::ALL.len());
    }

    #[test]
    fn report_separates_errors_and_notes() {
        let mut r = Report::default();
        r.push(Diagnostic::network(Code::N_LOOPS, "two loops"));
        assert!(r.is_valid());
        let rows: Vec<usize> = (0..70).collect();
        r.push(Diagnostic::rows(
            Code::E_UNKNOWN_NODE,
            Which::Edges,
            "bad",
            &rows,
        ));
        assert!(!r.is_valid());
        let d = r.find(Code::E_UNKNOWN_NODE).unwrap();
        assert_eq!(d.rows.len(), ROW_CAP);
        assert_eq!(d.count, 70);
        assert_eq!(r.errors().count(), 1);
        assert_eq!(r.notes().count(), 1);
        assert!(r
            .error_text()
            .contains("bad (rows 1, 2, 3, 4, 5, ... (70 in all) of `edges`)"));
    }

    #[test]
    fn rows_text_forms() {
        assert_eq!(rows_text(&[2], 1), "row 3");
        assert_eq!(rows_text(&[0, 4], 2), "rows 1 and 5");
        assert_eq!(rows_text(&[0, 1, 2], 3), "rows 1, 2, and 3");
        assert_eq!(
            rows_text(&[0, 1, 2, 3, 4, 5], 6),
            "rows 1, 2, 3, 4, 5, ... (6 in all)"
        );
        assert_eq!(rows_text(&[0, 1], 9), "rows 1, 2, ... (9 in all)");
        assert_eq!(rows_text(&[], 0), "");
    }
}
