/*
 * Registration of the .Call entry points, grouped by the C file that
 * defines them and sorted by name within each group.
 */
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* glue.c: Graphs. */
extern SEXP C_betweenness(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_degree(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_graph_alive(SEXP);
extern SEXP C_graph_new(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_triad_census(SEXP, SEXP);

/* frame.c: The data model. */
extern SEXP C_frame_build(SEXP, SEXP);
extern SEXP C_frame_codes(SEXP, SEXP, SEXP);
extern SEXP C_frame_index(SEXP, SEXP);
extern SEXP C_frame_is_simple(SEXP);
extern SEXP C_frame_multiplicity(SEXP, SEXP);
extern SEXP C_frame_new(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_frame_simplify(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_frame_summary(SEXP, SEXP);
extern SEXP C_frame_validate(SEXP, SEXP);

/* files.c: Network files. */
extern SEXP C_read_network(SEXP, SEXP, SEXP);
extern SEXP C_write_network(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* descriptive.c: Descriptive and structural measures. */
extern SEXP C_clustering(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cut_structure(SEXP);
extern SEXP C_distances(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_eccentricity(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_edge_counts(SEXP);
extern SEXP C_path_summary(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_reachability(SEXP, SEXP, SEXP);
extern SEXP C_reciprocity(SEXP, SEXP);
extern SEXP C_transitivity(SEXP, SEXP, SEXP);

/* topography.c: Topography and centrality. */
extern SEXP C_centralization(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_closeness(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_components(SEXP, SEXP);
extern SEXP C_density(SEXP, SEXP);
extern SEXP C_eigenvector(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_fragmentation(SEXP, SEXP);
extern SEXP C_graph_simplify(SEXP, SEXP, SEXP);

/* spectral.c: Spectral centralities and mixing patterns. */
extern SEXP C_assortativity_degree(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_assortativity_nominal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_assortativity_values(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ei_index(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_group_shares(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_hits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_katz(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mixing_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pagerank(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_power_centrality(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* roles.c: Reach, equivalence, blockmodels, scaling, and space. */
extern SEXP C_automorphic_equivalence(SEXP, SEXP);
extern SEXP C_block_model(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_bridging_coefficient(SEXP);
extern SEXP C_classical_mds(SEXP, SEXP, SEXP);
extern SEXP C_edge_betweenness(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_edge_lengths(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_group_geometry(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_information_centrality(SEXP, SEXP);
extern SEXP C_mds_stress(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_proximity(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_reach(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_regular_equivalence(SEXP, SEXP, SEXP);

/* positions.c: Brokerage and positions. */
extern SEXP C_brokerage(SEXP, SEXP, SEXP);
extern SEXP C_concor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_structural_equivalence(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_structural_holes(SEXP, SEXP, SEXP, SEXP);

/* twomode.c: Ego networks and two-mode networks. */
extern SEXP C_ego_members(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ego_sizes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ego_summary(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_project(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_twomode_centrality(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_twomode_density(SEXP, SEXP);
extern SEXP C_twomode_eigenvector(SEXP, SEXP, SEXP, SEXP, SEXP);

/* subgroups.c: Cores, cliques, and communities. */
extern SEXP C_clique_number(SEXP, SEXP);
extern SEXP C_coreness(SEXP, SEXP);
extern SEXP C_louvain(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_max_cliques(SEXP, SEXP, SEXP);
extern SEXP C_modularity(SEXP, SEXP, SEXP, SEXP);

/* groups.c: Cohesive subsets, partitions, core-periphery, and random graphs. */
extern SEXP C_bicomponents(SEXP);
extern SEXP C_compactness(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_core_periphery(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_core_periphery_continuous(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_factions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_fast_greedy(SEXP, SEXP, SEXP);
extern SEXP C_girvan_newman(SEXP, SEXP);
extern SEXP C_k_plexes(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_lambda_connectivity(SEXP);
extern SEXP C_leading_eigen(SEXP, SEXP);
extern SEXP C_n_cliques(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_configuration(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_gnm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_gnp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_preferential(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_sbm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_small_world(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_random_uman(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_small_world(SEXP, SEXP, SEXP, SEXP, SEXP);

/* cyber.c: Key players, group centrality, and flow betweenness. */
extern SEXP C_flow_betweenness(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_group_centrality(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_key_players(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_kpp_score(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* qap.c: Quadratic assignment procedure. */
extern SEXP C_mrqap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP);
extern SEXP C_qap_cor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* temporal.c: Temporal networks. */
extern SEXP C_temporal_active_at(SEXP, SEXP);
extern SEXP C_temporal_aggregate(SEXP, SEXP);
extern SEXP C_temporal_burstiness(SEXP, SEXP, SEXP);
extern SEXP C_temporal_closeness(SEXP, SEXP, SEXP);
extern SEXP C_temporal_inter_event(SEXP, SEXP, SEXP);
extern SEXP C_temporal_new(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_temporal_paths(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_temporal_reach(SEXP, SEXP, SEXP);
extern SEXP C_temporal_windows(SEXP, SEXP, SEXP);

/* layout.c: Layouts. */
extern SEXP C_layout_fr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_layout_simple(SEXP, SEXP, SEXP);
extern SEXP C_layout_stress(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_layout_stress_value(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_layout_tree(SEXP, SEXP);
extern SEXP C_layout_twomode(SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    /* glue.c */
    {"C_betweenness", (DL_FUNC)&C_betweenness, 5},
    {"C_degree", (DL_FUNC)&C_degree, 4},
    {"C_graph_alive", (DL_FUNC)&C_graph_alive, 1},
    {"C_graph_new", (DL_FUNC)&C_graph_new, 5},
    {"C_triad_census", (DL_FUNC)&C_triad_census, 2},
    /* frame.c */
    {"C_frame_build", (DL_FUNC)&C_frame_build, 2},
    {"C_frame_codes", (DL_FUNC)&C_frame_codes, 3},
    {"C_frame_index", (DL_FUNC)&C_frame_index, 2},
    {"C_frame_is_simple", (DL_FUNC)&C_frame_is_simple, 1},
    {"C_frame_multiplicity", (DL_FUNC)&C_frame_multiplicity, 2},
    {"C_frame_new", (DL_FUNC)&C_frame_new, 9},
    {"C_frame_simplify", (DL_FUNC)&C_frame_simplify, 5},
    {"C_frame_summary", (DL_FUNC)&C_frame_summary, 2},
    {"C_frame_validate", (DL_FUNC)&C_frame_validate, 2},
    /* files.c */
    {"C_read_network", (DL_FUNC)&C_read_network, 3},
    {"C_write_network", (DL_FUNC)&C_write_network, 7},
    /* descriptive.c */
    {"C_clustering", (DL_FUNC)&C_clustering, 5},
    {"C_cut_structure", (DL_FUNC)&C_cut_structure, 1},
    {"C_distances", (DL_FUNC)&C_distances, 6},
    {"C_eccentricity", (DL_FUNC)&C_eccentricity, 5},
    {"C_edge_counts", (DL_FUNC)&C_edge_counts, 1},
    {"C_path_summary", (DL_FUNC)&C_path_summary, 7},
    {"C_reachability", (DL_FUNC)&C_reachability, 3},
    {"C_reciprocity", (DL_FUNC)&C_reciprocity, 2},
    {"C_transitivity", (DL_FUNC)&C_transitivity, 3},
    /* topography.c */
    {"C_centralization", (DL_FUNC)&C_centralization, 7},
    {"C_closeness", (DL_FUNC)&C_closeness, 8},
    {"C_components", (DL_FUNC)&C_components, 2},
    {"C_density", (DL_FUNC)&C_density, 2},
    {"C_eigenvector", (DL_FUNC)&C_eigenvector, 8},
    {"C_fragmentation", (DL_FUNC)&C_fragmentation, 2},
    {"C_graph_simplify", (DL_FUNC)&C_graph_simplify, 3},
    /* spectral.c */
    {"C_assortativity_degree", (DL_FUNC)&C_assortativity_degree, 5},
    {"C_assortativity_nominal", (DL_FUNC)&C_assortativity_nominal, 6},
    {"C_assortativity_values", (DL_FUNC)&C_assortativity_values, 6},
    {"C_ei_index", (DL_FUNC)&C_ei_index, 6},
    {"C_group_shares", (DL_FUNC)&C_group_shares, 4},
    {"C_hits", (DL_FUNC)&C_hits, 7},
    {"C_katz", (DL_FUNC)&C_katz, 10},
    {"C_mixing_matrix", (DL_FUNC)&C_mixing_matrix, 7},
    {"C_pagerank", (DL_FUNC)&C_pagerank, 10},
    {"C_power_centrality", (DL_FUNC)&C_power_centrality, 9},
    /* roles.c */
    {"C_automorphic_equivalence", (DL_FUNC)&C_automorphic_equivalence, 2},
    {"C_block_model", (DL_FUNC)&C_block_model, 5},
    {"C_bridging_coefficient", (DL_FUNC)&C_bridging_coefficient, 1},
    {"C_classical_mds", (DL_FUNC)&C_classical_mds, 3},
    {"C_edge_betweenness", (DL_FUNC)&C_edge_betweenness, 5},
    {"C_edge_lengths", (DL_FUNC)&C_edge_lengths, 5},
    {"C_group_geometry", (DL_FUNC)&C_group_geometry, 6},
    {"C_information_centrality", (DL_FUNC)&C_information_centrality, 2},
    {"C_mds_stress", (DL_FUNC)&C_mds_stress, 4},
    {"C_proximity", (DL_FUNC)&C_proximity, 7},
    {"C_reach", (DL_FUNC)&C_reach, 7},
    {"C_regular_equivalence", (DL_FUNC)&C_regular_equivalence, 3},
    /* positions.c */
    {"C_brokerage", (DL_FUNC)&C_brokerage, 3},
    {"C_concor", (DL_FUNC)&C_concor, 6},
    {"C_structural_equivalence", (DL_FUNC)&C_structural_equivalence, 4},
    {"C_structural_holes", (DL_FUNC)&C_structural_holes, 4},
    /* twomode.c */
    {"C_ego_members", (DL_FUNC)&C_ego_members, 6},
    {"C_ego_sizes", (DL_FUNC)&C_ego_sizes, 6},
    {"C_ego_summary", (DL_FUNC)&C_ego_summary, 5},
    {"C_project", (DL_FUNC)&C_project, 6},
    {"C_twomode_centrality", (DL_FUNC)&C_twomode_centrality, 6},
    {"C_twomode_density", (DL_FUNC)&C_twomode_density, 2},
    {"C_twomode_eigenvector", (DL_FUNC)&C_twomode_eigenvector, 5},
    /* subgroups.c */
    {"C_clique_number", (DL_FUNC)&C_clique_number, 2},
    {"C_coreness", (DL_FUNC)&C_coreness, 2},
    {"C_louvain", (DL_FUNC)&C_louvain, 4},
    {"C_max_cliques", (DL_FUNC)&C_max_cliques, 3},
    {"C_modularity", (DL_FUNC)&C_modularity, 4},
    /* groups.c */
    {"C_bicomponents", (DL_FUNC)&C_bicomponents, 1},
    {"C_compactness", (DL_FUNC)&C_compactness, 5},
    {"C_core_periphery", (DL_FUNC)&C_core_periphery, 8},
    {"C_core_periphery_continuous", (DL_FUNC)&C_core_periphery_continuous, 8},
    {"C_factions", (DL_FUNC)&C_factions, 8},
    {"C_fast_greedy", (DL_FUNC)&C_fast_greedy, 3},
    {"C_girvan_newman", (DL_FUNC)&C_girvan_newman, 2},
    {"C_k_plexes", (DL_FUNC)&C_k_plexes, 5},
    {"C_lambda_connectivity", (DL_FUNC)&C_lambda_connectivity, 1},
    {"C_leading_eigen", (DL_FUNC)&C_leading_eigen, 2},
    {"C_n_cliques", (DL_FUNC)&C_n_cliques, 5},
    {"C_random_configuration", (DL_FUNC)&C_random_configuration, 7},
    {"C_random_gnm", (DL_FUNC)&C_random_gnm, 7},
    {"C_random_gnp", (DL_FUNC)&C_random_gnp, 7},
    {"C_random_preferential", (DL_FUNC)&C_random_preferential, 5},
    {"C_random_sbm", (DL_FUNC)&C_random_sbm, 7},
    {"C_random_small_world", (DL_FUNC)&C_random_small_world, 6},
    {"C_random_uman", (DL_FUNC)&C_random_uman, 6},
    {"C_small_world", (DL_FUNC)&C_small_world, 5},
    /* cyber.c */
    {"C_flow_betweenness", (DL_FUNC)&C_flow_betweenness, 4},
    {"C_group_centrality", (DL_FUNC)&C_group_centrality, 8},
    {"C_key_players", (DL_FUNC)&C_key_players, 11},
    {"C_kpp_score", (DL_FUNC)&C_kpp_score, 7},
    /* qap.c */
    {"C_mrqap", (DL_FUNC)&C_mrqap, 15},
    {"C_qap_cor", (DL_FUNC)&C_qap_cor, 10},
    /* temporal.c */
    {"C_temporal_active_at", (DL_FUNC)&C_temporal_active_at, 2},
    {"C_temporal_aggregate", (DL_FUNC)&C_temporal_aggregate, 2},
    {"C_temporal_burstiness", (DL_FUNC)&C_temporal_burstiness, 3},
    {"C_temporal_closeness", (DL_FUNC)&C_temporal_closeness, 3},
    {"C_temporal_inter_event", (DL_FUNC)&C_temporal_inter_event, 3},
    {"C_temporal_new", (DL_FUNC)&C_temporal_new, 6},
    {"C_temporal_paths", (DL_FUNC)&C_temporal_paths, 5},
    {"C_temporal_reach", (DL_FUNC)&C_temporal_reach, 3},
    {"C_temporal_windows", (DL_FUNC)&C_temporal_windows, 3},
    /* layout.c */
    {"C_layout_fr", (DL_FUNC)&C_layout_fr, 8},
    {"C_layout_simple", (DL_FUNC)&C_layout_simple, 3},
    {"C_layout_stress", (DL_FUNC)&C_layout_stress, 9},
    {"C_layout_stress_value", (DL_FUNC)&C_layout_stress_value, 5},
    {"C_layout_tree", (DL_FUNC)&C_layout_tree, 2},
    {"C_layout_twomode", (DL_FUNC)&C_layout_twomode, 3},
    {NULL, NULL, 0}};

/* The symbol name must match the package name; tools/rename.py updates it. */
void R_init_knappnet(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
