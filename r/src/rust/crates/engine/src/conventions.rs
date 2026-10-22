//! How reference tools normalize centrality and treat disconnected networks.
//!
//! The same measure name means different numbers in different programs.
//! This module is the single place where those differences are recorded,
//! and [`Convention`] is the switch that reproduces each tool. Measures take
//! explicit options (for example [`crate::metrics::closeness::Unreachable`]);
//! the functions here say which options each tool's convention implies.
//!
//! Sources, read for this module (versions current in October 2026):
//! - igraph (`csardi2006igraph`): C library 0.10 sources
//!   `src/centrality/closeness.c`, `centralization.c`, `eigenvector.c`, and
//!   `src/properties/basic_properties.c`; R package 2.x wrappers
//!   (`centr_degree(loops = TRUE)`, `centr_clo(mode = "out")`).
//! - NetworkX (`hagberg2008exploring`) 3.7: `closeness_centrality`,
//!   `harmonic_centrality`, `eigenvector_centrality_numpy`, `density`.
//! - sna (`butts2008sna`): `closeness`, `evcent`, `degree`, `betweenness`,
//!   `centralization`, `gden` (R sources `nli.R`, `gli.R`).
//! - UCINET (`borgatti2002ucinet`): the help pages for Closeness,
//!   Eigenvector, and Freeman betweenness (`ucinet_help_closeness`,
//!   `ucinet_help_eigenvector`, `ucinet_help_betweenness`) and the published
//!   port of UCINET's routines by its authors (`borgatti2026xucinet`), which
//!   states each formula from UCINET's Delphi sources and checks it against
//!   UCINET 6.849 output. UCINET itself cannot be run here, so UCINET
//!   conventions are validated only once golden files from UCINET are
//!   committed (see `validation/manifest.d/topography.yaml`).
//!
//! # Normalization and disconnected networks
//!
//! n is the number of nodes, d(i, j) a shortest-path length, r_i the number
//! of other nodes i reaches, and S_i the sum of distances to them.
//!
//! | | igraph | NetworkX | sna | UCINET |
//! |---|---|---|---|---|
//! | Ties | edges as given (parallel edges count, loops count) | as igraph | adjacency matrix: parallel edges collapse, diagonal ignored | as sna |
//! | Parallel edges on shortest paths | distinct paths | one path, through the lightest edge | one path | one path |
//! | Normalized degree | d / (n - 1) | d / (n - 1) | (no option) | d / (n - 1); valued data d / (max value (n - 1)) |
//! | Normalized betweenness | 2b / ((n-1)(n-2)) undirected, b / ((n-1)(n-2)) directed | as igraph | (no option; `rescale` divides by the sum) | as igraph, as a percentage |
//! | Closeness | 1 / S_i, or r_i / S_i normalized; only reachable nodes count; NaN if r_i = 0 | (r_i / S_i)(r_i / (n - 1)) (Wasserman and Faust 1994, p. 201); 0 if r_i = 0; distances *to* the node | (n - 1) / S_i, 0 if any node is unreachable | (n - 1) / S_i with every unreachable distance set to the largest observed distance plus 1 (default; also n, or ignore) |
//! | Closeness direction (default) | out | in | out | out and in |
//! | Harmonic | sum 1/d, optionally / (n - 1) | sum 1/d, distances to the node, never normalized | `suminvdir`: sum 1/d / (n - 1) | "reciprocal": sum 1/d / (n - 1) |
//! | Eigenvector scale | maximum = 1 | unit Euclidean norm | unit Euclidean norm | unit Euclidean norm |
//! | Eigenvector, directed | in-ties (left eigenvector); R default `directed = FALSE` ignores direction | in-ties | out-ties (right eigenvector) | symmetrized by maximum |
//! | Eigenvector, disconnected | computed; zero outside the dominant component | refused | computed | lambda-squared scores within components when more than one component has ties |
//! | Centralization | Freeman 1978/79: sum(max - c_i) / theoretical maximum | none | as igraph, own maxima | as igraph, percentages, own maxima |
//!
//! # Freeman centralization maxima
//!
//! Centralization is sum_i (max_j c_j - c_i) divided by the largest value
//! that sum can take in a network of the same size (Freeman 1978,
//! `freeman1978centrality`), which a star attains. Tools differ in what the
//! star may contain (loops) and in how closeness treats a star's
//! unreachable pairs:
//!
//! | Measure | igraph | sna | UCINET |
//! |---|---|---|---|
//! | Degree, undirected | (n-1)(n-2); (n-1)n with loops | (n-1)(n-2); (n-1)(n-1) with `diag` | (n-1)(n-2) |
//! | Degree in or out, directed | (n-1)^2; (n-1)n with loops | (n-1)^2; (n-1)n with `diag` | (n-1)^2 |
//! | Degree all, directed | 2(n-1)(n-2); 2(n-1)^2 with loops | (n-1)(2n-4); (n-1)(2n-3) with `diag` | (not reported) |
//! | Betweenness | (n-1)^2 (n-2), halved if undirected | as igraph | as igraph |
//! | Closeness, undirected or mode all | (n-1)(n-2) / (2n-3) | as igraph | as igraph, and also for directed data |
//! | Closeness in or out, directed | (n-1)^2 / n | as igraph | (n-1)(n-2) / (2n-3) |
//!
//! Closeness centralization uses each tool's normalized closeness. igraph's
//! is NaN when some node reaches no other; UCINET reports none unless every
//! node reaches every other.
//!
//! UCINET's maxima are those of version 6.849 (2026): its log for the
//! `campnet` data reports an in-degree centralization of 0.1246, which is
//! (18 * 5 - 54) / 17^2. Logs from builds copyrighted 2002-12 report
//! directed degree centralization over (n-1)(n-2) instead (for example
//! 53.542% for a 43-node network with maximum out-degree 42 and 884 ties,
//! which is 100 (43 * 42 - 884) / (42 * 41)); results from those versions
//! will differ for directed data.
//!
//! UCINET reports normalized scores and centralization as percentages in
//! some routines and proportions in others; this engine always uses
//! proportions.

use crate::error::{Error, Result};
use crate::metrics::closeness::Unreachable;
use crate::metrics::eigenvector::EigenScale;

/// A reference tool whose conventions a measure should reproduce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Convention {
    /// igraph (the default).
    #[default]
    Igraph,
    /// NetworkX.
    Networkx,
    /// The R package sna.
    Sna,
    /// UCINET.
    Ucinet,
}

impl Convention {
    /// Decodes the C interface codes: 0 = igraph, 1 = NetworkX, 2 = sna,
    /// 3 = UCINET.
    pub fn from_code(code: i32) -> Result<Convention> {
        match code {
            0 => Ok(Convention::Igraph),
            1 => Ok(Convention::Networkx),
            2 => Ok(Convention::Sna),
            3 => Ok(Convention::Ucinet),
            _ => Err(Error::InvalidArgument(format!(
                "unknown convention code {code}"
            ))),
        }
    }

    /// The tool's name, for messages.
    pub fn name(self) -> &'static str {
        match self {
            Convention::Igraph => "igraph",
            Convention::Networkx => "NetworkX",
            Convention::Sna => "sna",
            Convention::Ucinet => "UCINET",
        }
    }

    /// Whether the tool works on an adjacency matrix, so that parallel
    /// edges collapse into one tie and self-loops are ignored.
    pub fn uses_adjacency_matrix(self) -> bool {
        matches!(self, Convention::Sna | Convention::Ucinet)
    }

    /// How closeness treats pairs with no path between them.
    pub fn closeness_unreachable(self) -> Unreachable {
        match self {
            Convention::Igraph => Unreachable::Ignore,
            Convention::Networkx => Unreachable::WassermanFaust,
            Convention::Sna => Unreachable::Zero,
            Convention::Ucinet => Unreachable::MaxPlusOne,
        }
    }

    /// Whether the tool reports normalized closeness by default. igraph
    /// reports 1 / S_i unless asked; the others report only normalized
    /// scores.
    pub fn closeness_normalized(self) -> bool {
        !matches!(self, Convention::Igraph)
    }

    /// How eigenvector centrality is scaled.
    pub fn eigenvector_scale(self) -> EigenScale {
        match self {
            Convention::Igraph => EigenScale::Max,
            _ => EigenScale::Unit,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_round_trip() {
        for (code, c) in [
            (0, Convention::Igraph),
            (1, Convention::Networkx),
            (2, Convention::Sna),
            (3, Convention::Ucinet),
        ] {
            assert_eq!(Convention::from_code(code).unwrap(), c);
        }
        assert!(Convention::from_code(4).is_err());
    }
}
