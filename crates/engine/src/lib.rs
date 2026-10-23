#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
//!
//! ## Organization
//!
//! The engine is language-neutral: it works on 0-based integer node indices
//! and plain numeric arrays. Language bindings live in separate crates and
//! packages. Every algorithm documents its conventions and cites its source
//! in the project's `docs/references.bib`.
//!
//! | Area | Modules |
//! |---|---|
//! | Storage and shared machinery | [`graph`], [`bitset`], [`matrix`], [`linalg`], [`union_find`], [`rng`], [`parallel`], [`error`], [`conventions`] |
//! | Paths and weights | [`weights`] (how weights are read and which shortest-path algorithm applies) |
//! | Network files | [`io`] |
//! | Measures (centrality, cohesion, positions, structure) | [`metrics`] |
//! | Statistical models | [`qap`], [`ergm`], [`saom`] |
//! | Dynamics | [`temporal`] |
//! | Space | [`spatial`] |
//! | Drawing | [`layout`] |
//! | Random networks | [`generate`] |

// Storage and shared machinery.
pub mod bitset;
pub mod conventions;
pub mod error;
pub mod graph;
pub mod linalg;
pub mod matrix;
pub mod parallel;
pub mod rng;
pub mod union_find;
pub mod weights;

// Files.
pub mod io;

// Measures.
pub mod metrics;

// Statistical models.
pub mod ergm;
pub mod qap;
pub mod saom;

// Dynamics, space, drawing, and random networks.
pub mod generate;
pub mod layout;
pub mod spatial;
pub mod temporal;

pub use error::{Error, Result};
pub use graph::Graph;
pub use metrics::betweenness::{betweenness, BetweennessOptions};
pub use metrics::degree::{degree, strength, DegreeMode};
pub use metrics::triad::{triad_census, TRIAD_NAMES};
pub use weights::{PathAlgorithm, WeightKind};
