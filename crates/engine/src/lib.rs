//! Core graph storage and network analysis algorithms.
//!
//! The engine is language-neutral: it works on 0-based integer node indices
//! and plain numeric arrays. Language bindings live in separate crates and
//! packages. Every algorithm documents its conventions and cites its source
//! in `docs/references.bib`.

pub mod error;
pub mod graph;
pub mod metrics;
pub mod parallel;

pub use error::{Error, Result};
pub use graph::Graph;
