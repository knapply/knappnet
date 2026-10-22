#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
//!
//! ## Organization
//!
//! The engine is language-neutral: it works on 0-based integer node indices
//! and plain numeric arrays. Language bindings live in separate crates and
//! packages. Every algorithm documents its conventions and cites its source
//! in the project's `docs/references.bib`.

pub mod error;
pub mod graph;
pub mod io;
pub mod metrics;
pub mod parallel;

pub use error::{Error, Result};
pub use graph::Graph;
pub use metrics::betweenness::{betweenness, BetweennessOptions};
pub use metrics::degree::{degree, strength, DegreeMode};
pub use metrics::triad::{triad_census, TRIAD_NAMES};
