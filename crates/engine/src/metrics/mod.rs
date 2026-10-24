//! Network measures. Each module documents its conventions and sources.
//!
//! | Group | Modules |
//! |---|---|
//! | Centrality | [`degree`], [`betweenness`], [`closeness`], [`eigenvector`], [`spectral`] (PageRank, Katz, power, HITS), [`reach`] (reach, edge betweenness, information, bridging), [`brokerage`] |
//! | Structure and description | [`descriptive`] (reciprocity, transitivity, distances), [`topography`] (density, components, centralization), [`triad`], [`mixing`], [`ego`], [`bipartite`], [`linkpred`] (link prediction) |
//! | Cohesion and groups | [`cores`], [`cliques`], [`factions`] (bicomponents, n-cliques, k-plexes, lambda sets), [`community`] (modularity, Louvain), [`partition`] (Girvan-Newman, leading eigenvector, fast greedy), [`resilience`] (key players, flow betweenness) |
//! | Positions and roles | [`equivalence`] (structural equivalence, CONCOR), [`regular`], [`blockmodel`], [`coreperiphery`], [`scaling`], [`similarity`] (SimRank, Panther), [`edit_distance`], [`wl`] (colour refinement, Weisfeiler-Lehman hashes) |

// Centrality.
pub mod betweenness;
pub mod brokerage;
pub mod closeness;
pub mod degree;
pub mod eigenvector;
pub mod reach;
pub mod spectral;

// Structure and description.
pub mod bipartite;
pub mod descriptive;
pub mod ego;
pub mod linkpred;
pub mod mixing;
pub mod topography;
pub mod triad;

// Cohesion and groups.
pub mod cliques;
pub mod community;
pub mod cores;
pub mod factions;
pub mod partition;
pub mod resilience;

// Positions and roles.
pub mod blockmodel;
pub mod coreperiphery;
pub mod edit_distance;
pub mod equivalence;
pub mod regular;
pub mod scaling;
pub mod similarity;
pub mod wl;
