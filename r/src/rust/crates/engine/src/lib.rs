#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
//!
//! ## Organization
//!
//! The engine is language-neutral: it works on 0-based integer node indices
//! and plain numeric arrays. Language bindings live in separate crates and
//! packages. Every algorithm documents its conventions and cites its source
//! in the project's `docs/references.bib`.

pub mod bitset;
pub mod conventions;
pub mod dynamics;
pub mod ergm;
pub mod error;
pub mod generate;
pub mod graph;
pub mod io;
pub mod layout;
pub mod logs;
pub mod matrix;
pub mod metrics;
pub mod parallel;
pub mod qap;
pub mod rng;
pub mod saom;
pub mod spatial;
pub mod temporal;
pub mod weights;

pub use error::{Error, Result};
pub use graph::Graph;
pub use metrics::betweenness::{betweenness, BetweennessOptions};
pub use metrics::degree::{degree, strength, DegreeMode};
pub use metrics::triad::{triad_census, TRIAD_NAMES};
pub use weights::{PathAlgorithm, WeightKind};
