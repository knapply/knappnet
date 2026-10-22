//! Error types.

use std::fmt;

/// Errors produced while building a graph or running an algorithm.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Source and target arrays have different lengths.
    LengthMismatch {
        /// Which input is the wrong length.
        what: &'static str,
        /// The required length.
        expected: usize,
        /// The length supplied.
        found: usize,
    },
    /// An edge endpoint is negative or not less than the node count.
    NodeOutOfRange {
        /// Index of the offending edge.
        edge: usize,
        /// The invalid endpoint.
        node: i64,
        /// The node count.
        n: usize,
    },
    /// A weight is NaN or infinite.
    NonFiniteWeight {
        /// Index of the offending edge.
        edge: usize,
    },
    /// An algorithm requires strictly positive weights.
    NonPositiveWeight {
        /// Index of the offending edge.
        edge: usize,
        /// Its weight.
        weight: f64,
    },
    /// The graph is too large for the index type.
    TooLarge {
        /// Which quantity is too large.
        what: &'static str,
        /// Its value.
        value: usize,
    },
    /// The algorithm needs edge weights but the graph has none.
    MissingWeights,
    /// An argument is outside its allowed set of values.
    InvalidArgument(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::LengthMismatch {
                what,
                expected,
                found,
            } => write!(f, "{what} has length {found}, but {expected} was expected"),
            Error::NodeOutOfRange { edge, node, n } => write!(
                f,
                "edge {edge} refers to node index {node}, but valid indices are 0..{n}"
            ),
            Error::NonFiniteWeight { edge } => {
                write!(f, "edge {edge} has a missing, NaN, or infinite weight")
            }
            Error::NonPositiveWeight { edge, weight } => write!(
                f,
                "edge {edge} has weight {weight}; this algorithm treats weights as \
                 distances and requires them to be strictly positive"
            ),
            Error::TooLarge { what, value } => {
                write!(f, "{what} ({value}) exceeds the supported maximum")
            }
            Error::MissingWeights => write!(f, "the graph has no edge weights"),
            Error::InvalidArgument(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for Error {}

/// Result type used throughout the engine.
pub type Result<T> = std::result::Result<T, Error>;
