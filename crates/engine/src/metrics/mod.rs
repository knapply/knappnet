//! Network measures. Each module documents its conventions and sources.
//!
//! | Group | Modules |
//! |---|---|
//! | Centrality | [`degree`], [`betweenness`], [`closeness`], [`eigenvector`], [`spectral`] (PageRank, Katz, power, HITS), [`reach`] (reach, edge betweenness, information, bridging), [`brokerage`], [`path_centrality`] (subset betweenness, percolation, load), [`communicability`], [`current_flow`] (current-flow betweenness and closeness, resistance distance, Kemeny's constant, second-order centrality), [`spreading`] (VoteRank, Laplacian centrality, dispersion), [`hierarchy`] (trophic levels, reaching centrality) |
//! | Structure and description | [`descriptive`] (reciprocity, transitivity, distances), [`distance_measures`] (center, periphery, barycenter, Wiener indices, efficiency, closeness vitality), [`topography`] (density, components, centralization), [`triad`], [`mixing`], [`ego`], [`bipartite`] |
//! | Connectivity | [`connectivity`] (node and edge connectivity, minimum cuts, disjoint paths, k-components, k-edge-components, condensation, local bridges) |
//! | Trees and DAGs | [`trees`] (spanning trees, branchings, arborescences), [`dag`] (topological orders, transitive closure and reduction, longest paths, dominators, lowest common ancestors) |
//! | Cycles | [`cycles`] (cycle bases, minimum cycle bases, simple and chordless cycles, girth, finding a cycle) |
//! | Cohesion and groups | [`cores`], [`cliques`], [`factions`] (bicomponents, n-cliques, k-plexes, lambda sets), [`community`] (modularity, Louvain), [`partition`] (Girvan-Newman, leading eigenvector, fast greedy), [`resilience`] (key players, flow betweenness) |
//! | Positions and roles | [`equivalence`] (structural equivalence, CONCOR), [`regular`], [`blockmodel`], [`coreperiphery`], [`scaling`] |

// Centrality.
pub mod betweenness;
pub mod brokerage;
pub mod closeness;
pub mod communicability;
pub mod current_flow;
pub mod degree;
pub mod eigenvector;
pub mod hierarchy;
pub mod path_centrality;
pub mod reach;
pub mod spectral;
pub mod spreading;

// Structure and description.
pub mod bipartite;
pub mod descriptive;
pub mod distance_measures;
pub mod ego;
pub mod mixing;
pub mod topography;
pub mod triad;

// Connectivity.
pub mod connectivity;

// Trees and DAGs.
pub mod dag;
pub mod trees;

// Cycles.
pub mod cycles;

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
