//! Network measures. Each module documents its conventions and sources.
//!
//! | Group | Modules |
//! |---|---|
//! | Centrality | [`degree`], [`betweenness`], [`closeness`], [`eigenvector`], [`spectral`] (PageRank, Katz, power, HITS), [`reach`] (reach, edge betweenness, information, bridging), [`brokerage`] |
//! | Structure and description | [`descriptive`] (reciprocity, transitivity, distances), [`topography`] (density, components, centralization), [`triad`], [`mixing`], [`ego`], [`bipartite`] |
//! | Connectivity | [`connectivity`] (node and edge connectivity, minimum cuts, disjoint paths, k-components, k-edge-components, condensation, local bridges) |
//! | Trees and DAGs | [`trees`] (spanning trees, branchings, arborescences), [`dag`] (topological orders, transitive closure and reduction, longest paths, dominators, lowest common ancestors) |
//! | Cohesion and groups | [`cores`], [`cliques`], [`factions`] (bicomponents, n-cliques, k-plexes, lambda sets), [`community`] (modularity, Louvain), [`partition`] (Girvan-Newman, leading eigenvector, fast greedy), [`resilience`] (key players, flow betweenness) |
//! | Positions and roles | [`equivalence`] (structural equivalence, CONCOR), [`regular`], [`blockmodel`], [`coreperiphery`], [`scaling`] |

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
pub mod mixing;
pub mod topography;
pub mod triad;

// Connectivity.
pub mod connectivity;

// Trees and DAGs.
pub mod dag;
pub mod trees;

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
pub mod equivalence;
pub mod regular;
pub mod scaling;
