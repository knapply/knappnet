//! Network measures. Each module documents its conventions and sources.
//!
//! | Group | Modules |
//! |---|---|
//! | Centrality | [`degree`], [`betweenness`], [`closeness`], [`eigenvector`], [`spectral`] (PageRank, Katz, power, HITS), [`reach`] (reach, edge betweenness, information, bridging), [`brokerage`], [`path_centrality`] (subset betweenness, percolation, load), [`communicability`], [`current_flow`] (current-flow betweenness and closeness, resistance distance, Kemeny's constant, second-order centrality), [`spreading`] (VoteRank, Laplacian centrality, dispersion), [`hierarchy`] (trophic levels, reaching centrality) |
//! | Structure and description | [`descriptive`] (reciprocity, transitivity, distances), [`squares`] (square clustering, generalized degree), [`distance_measures`] (center, periphery, barycenter, Wiener indices, efficiency, closeness vitality), [`topography`] (density, components, centralization), [`triad`], [`mixing`], [`assortativity`] (NetworkX mixing), [`cut_measures`] (conductance, expansion), [`ego`], [`bipartite`], [`linkpred`] (link prediction) |
//! | Connectivity | [`connectivity`] (node and edge connectivity, minimum cuts, disjoint paths, k-components, k-edge-components, condensation, local bridges) |
//! | Trees and DAGs | [`trees`] (spanning trees, branchings, arborescences), [`dag`] (topological orders, transitive closure and reduction, longest paths, dominators, lowest common ancestors) |
//! | Cycles | [`cycles`] (cycle bases, minimum cycle bases, simple and chordless cycles, girth, finding a cycle) |
//! | Cohesion and groups | [`cores`] (cores, shells, trusses, onion layers), [`cliques`], [`matching`], [`coloring`] (greedy colouring, independent sets), [`factions`] (bicomponents, n-cliques, k-plexes, lambda sets), [`community`] (modularity, Louvain), [`partition`] (Girvan-Newman, leading eigenvector, fast greedy), [`resilience`] (key players, flow betweenness), [`rewiring`] (edge swaps, reference graphs, small-world coefficients, rich club) |
//! | Positions and roles | [`equivalence`] (structural equivalence, CONCOR), [`regular`], [`blockmodel`], [`coreperiphery`], [`scaling`], [`similarity`] (SimRank, Panther), [`edit_distance`], [`isomorphism`] (VF2++, tree isomorphism), [`wl`] (colour refinement, Weisfeiler-Lehman hashes) |

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
pub mod assortativity;
pub mod bipartite;
pub mod cut_measures;
pub mod descriptive;
pub mod distance_measures;
pub mod ego;
pub mod linkpred;
pub mod mixing;
pub mod squares;
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
pub mod coloring;
pub mod community;
pub mod cores;
pub mod factions;
pub mod matching;
pub mod partition;
pub mod resilience;
pub mod rewiring;

// Positions and roles.
pub mod blockmodel;
pub mod coreperiphery;
pub mod edit_distance;
pub mod equivalence;
pub mod isomorphism;
pub mod regular;
pub mod scaling;
pub mod similarity;
pub mod wl;
