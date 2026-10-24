# Function names

Every exported function of the R package starts with `net_`, so that
completion lists the whole package after three keystrokes and no name
collides with base R, igraph, sna, network, tidygraph, or sf (decision 1
in `docs/DECISIONS.md`). The package is unreleased, so the old names were
removed without deprecated aliases. This file records the rules and the
complete map from the names used before October 2026.

## Rules

1. **`net_` plus the measure or action as an analyst says it**:
   `net_degree()`, `net_betweenness()`, `net_density()`,
   `net_subnetwork()`. Words that only repeated the prefix are dropped
   (`graph_density` is `net_density`, `ego_network` is `net_ego`,
   `layout_network` is `net_layout`).
2. **Standard short forms for long names**, where the short form is the
   one the literature and other tools use: `net_eigen` (eigenvector
   centrality), `net_harmonic`, `net_bonacich` (power centrality),
   `net_hub` and `net_authority`, `net_struct_equiv`,
   `net_regular_equiv`, `net_auto_equiv`, `net_mds`, `net_hclust`,
   `net_assort_degree`, `net_degree_dist`, `net_dw_reach` and
   `net_dw_fragmentation` (distance-weighted).
3. **Family words are kept where they help discovery and stay short**:
   `net_random_gnp`, `net_random_sbm`, `net_plot_degree`,
   `net_plot_partition`, `net_temporal_reach`, `net_twomode_degree`,
   `net_ego_size`, `net_key_players`. Community detection is named by
   method: `net_louvain`, `net_girvan_newman`, `net_leading_eigen`,
   `net_fast_greedy`.
4. **Conversions**: `as_net_frame()` is the S3 generic that converts
   other objects into a net_frame (decision 16; it replaced the planned
   `net_from()`), and `net_to_igraph()`, `net_to_network()`,
   `net_to_tbl_graph()`, `net_to_sf()`, `net_to_sfnetwork()`, and
   `net_to_egor()` convert out. File input and output are `net_read()` and
   `net_write()`; the constructor is `net_frame()`.
5. **Predicates** start `net_is_`: `net_is_twomode()`, `net_is_bipartite()`,
   `net_is_simple()`, `net_is_temporal()`, `net_is_spatial()`,
   `net_is_multilayer()`.
6. **Length**: at most 20 characters where clarity allows. The longest
   names are `net_twomode_betweenness` (23), `net_temporal_distances`,
   `net_temporal_closeness`, and `net_random_small_world` (22), whose
   family word and measure both carry meaning.
7. **S3 methods keep their generics' names**: `print()`, `summary()`,
   `plot()`, `as.matrix()`, and the registrations for suggested packages
   (`igraph::as.igraph()`, `network::as.network()`,
   `tidygraph::as_tbl_graph()`, `sf::st_as_sf()`,
   `sfnetworks::as_sfnetwork()`).

Classes are `net_frame`, `temporal_net_frame`, `spatial_net_frame`,
`spatiotemporal_net_frame`, and `multilayer_net_frame` (decision 17). The
section "Whole-network structure" replaces "Topography" (decision 20).

## Map of old names

| Old name | New name |
|---|---|
| `aggregate_network` | `net_aggregate` |
| `as_netframe` | `as_net_frame` |
| `as_twomode` | `net_as_twomode` |
| `assortativity` | `net_assortativity` |
| `assortativity_degree` | `net_assort_degree` |
| `assortativity_nominal` | `net_assort_nominal` |
| `authority_score` | `net_authority` |
| `automorphic_equivalence` | `net_auto_equiv` |
| `average_clustering` | `net_avg_clustering` |
| `average_degree` | `net_avg_degree` |
| `betweenness` | `net_betweenness` |
| `bicomponents` | `net_bicomponents` |
| `block_summary` | `net_block_summary` |
| `blockmodel` | `net_blockmodel` |
| `bridges` | `net_bridges` |
| `bridging_centrality` | `net_bridging` |
| `bridging_coefficient` | `net_bridging_coef` |
| `brokerage` | `net_brokerage` |
| `burstiness` | `net_burstiness` |
| `centralization` | `net_centralization` |
| `centralize` | `net_centralize` |
| `check_geometry` | `net_check_geometry` |
| `classical_mds` | `net_mds` |
| `clique_number` | `net_clique_number` |
| `closeness` | `net_closeness` |
| `cluster_fast_greedy` | `net_fast_greedy` |
| `cluster_girvan_newman` | `net_girvan_newman` |
| `cluster_hierarchical` | `net_hclust` |
| `cluster_leading_eigen` | `net_leading_eigen` |
| `cluster_louvain` | `net_louvain` |
| `clustering` | `net_clustering` |
| `combine_relations` | `net_combine_relations` |
| `compactness` | `net_compactness` |
| `components` | `net_components` |
| `concor` | `net_concor` |
| `constraint` | `net_constraint` |
| `core_periphery` | `net_core_periphery` |
| `coreness` | `net_coreness` |
| `cug_test` | `net_cug_test` |
| `cut_clusters` | `net_cut_clusters` |
| `cut_points` | `net_cut_points` |
| `degree` | `net_degree` |
| `degree_distribution` | `net_degree_dist` |
| `diameter` | `net_diameter` |
| `dichotomize` | `net_dichotomize` |
| `distance_decay` | `net_distance_decay` |
| `distance_weighted_fragmentation` | `net_dw_fragmentation` |
| `distance_weighted_reach` | `net_dw_reach` |
| `distances` | `net_distances` |
| `eccentricity` | `net_eccentricity` |
| `edge_betweenness` | `net_edge_betweenness` |
| `edge_lengths` | `net_edge_lengths` |
| `effective_size` | `net_effective_size` |
| `efficiency` | `net_efficiency` |
| `ego_betweenness` | `net_ego_betweenness` |
| `ego_density` | `net_ego_density` |
| `ego_network` | `net_ego` |
| `ego_networks` | `net_egos` |
| `ego_size` | `net_ego_size` |
| `ego_summary` | `net_ego_summary` |
| `ei_index` | `net_ei_index` |
| `eigenvector_centrality` | `net_eigen` |
| `example_network` | `net_example` |
| `example_networks` | `net_examples` |
| `factions` | `net_factions` |
| `filter_nodes` | `net_filter_nodes` |
| `flow_betweenness` | `net_flow_betweenness` |
| `fragmentation` | `net_fragmentation` |
| `from_lines` | `net_from_lines` |
| `from_polygons` | `net_from_polygons` |
| `graph_density` | `net_density` |
| `group_centrality` | `net_group_centrality` |
| `group_geometry` | `net_group_geometry` |
| `harmonic_centrality` | `net_harmonic` |
| `hierarchy` | `net_hierarchy` |
| `homophily` | `net_homophily` |
| `hub_score` | `net_hub` |
| `incidence_matrix` | `net_incidence` |
| `information_centrality` | `net_information` |
| `inter_event_times` | `net_inter_event` |
| `is_temporal` | `net_is_temporal` |
| `is_twomode` | `net_is_twomode` |
| `join_attributes` | `net_join` |
| `k_plexes` | `net_k_plexes` |
| `katz_centrality` | `net_katz` |
| `key_player_score` | `net_key_player_score` |
| `key_players` | `net_key_players` |
| `lambda_sets` | `net_lambda_sets` |
| `largest_component` | `net_largest_component` |
| `layout_network` | `net_layout` |
| `layout_stress` | `net_layout_stress` |
| `max_cliques` | `net_max_cliques` |
| `mds_distances` | `net_mds_distances` |
| `mds_stress` | `net_mds_stress` |
| `mean_distance` | `net_mean_distance` |
| `mixing_matrix` | `net_mixing` |
| `modularity` | `net_modularity` |
| `mrqap` | `net_mrqap` |
| `n_clans` | `net_n_clans` |
| `n_cliques` | `net_n_cliques` |
| `netframe` | `net_frame` |
| `network_at` | `net_at` |
| `network_window` | `net_window` |
| `pagerank` | `net_pagerank` |
| `plot_centrality` | `net_plot_centrality` |
| `plot_component_sizes` | `net_plot_components` |
| `plot_coreness` | `net_plot_coreness` |
| `plot_degree_distribution` | `net_plot_degree` |
| `plot_edge_lengths` | `net_plot_edge_lengths` |
| `plot_ego` | `net_plot_ego` |
| `plot_matrix` | `net_plot_matrix` |
| `plot_mixing_matrix` | `net_plot_mixing` |
| `plot_partition` | `net_plot_partition` |
| `plot_robustness` | `net_plot_robustness` |
| `plot_rolling` | `net_plot_rolling` |
| `plot_scores` | `net_plot_scores` |
| `plot_spatial` | `net_plot_spatial` |
| `plot_timeline` | `net_plot_timeline` |
| `plot_triad_census` | `net_plot_triads` |
| `power_centrality` | `net_bonacich` |
| `project` | `net_project` |
| `proximity_network` | `net_proximity` |
| `proximity_prestige` | `net_prox_prestige` |
| `qap_cor` | `net_qap_cor` |
| `radius` | `net_radius` |
| `random_configuration` | `net_random_config` |
| `random_gnm` | `net_random_gnm` |
| `random_gnp` | `net_random_gnp` |
| `random_preferential` | `net_random_pref` |
| `random_sbm` | `net_random_sbm` |
| `random_small_world` | `net_random_small_world` |
| `random_uman` | `net_random_uman` |
| `reach` | `net_reach` |
| `reachability` | `net_reachability` |
| `reachability_ratio` | `net_reach_ratio` |
| `read_network` | `net_read` |
| `reciprocity` | `net_reciprocity` |
| `recode_attribute` | `net_recode` |
| `regular_equivalence` | `net_regular_equiv` |
| `remove_isolates` | `net_remove_isolates` |
| `remove_pendants` | `net_remove_pendants` |
| `rolling` | `net_rolling` |
| `small_world` | `net_small_world` |
| `snapshots` | `net_snapshots` |
| `spatial_filter` | `net_spatial_filter` |
| `spatial_join` | `net_spatial_join` |
| `spatial_network` | `net_spatial` |
| `strength` | `net_strength` |
| `structural_equivalence` | `net_struct_equiv` |
| `subnetwork` | `net_subnetwork` |
| `symmetrize` | `net_symmetrize` |
| `temporal_closeness` | `net_temporal_closeness` |
| `temporal_distances` | `net_temporal_distances` |
| `temporal_network` | `net_temporal` |
| `temporal_reach` | `net_temporal_reach` |
| `to_igraph` | `net_to_igraph` |
| `to_network` | `net_to_network` |
| `to_sf` | `net_to_sf` |
| `to_sfnetwork` | `net_to_sfnetwork` |
| `to_tbl_graph` | `net_to_tbl_graph` |
| `transform_network` | `net_transform` |
| `transitivity` | `net_transitivity` |
| `transpose_network` | `net_transpose` |
| `triad_census` | `net_triad_census` |
| `twomode` | `net_twomode` |
| `twomode_betweenness` | `net_twomode_betweenness` |
| `twomode_closeness` | `net_twomode_closeness` |
| `twomode_degree` | `net_twomode_degree` |
| `twomode_density` | `net_twomode_density` |
| `twomode_eigenvector` | `net_twomode_eigen` |
| `write_network` | `net_write` |

## Functions added with the new data model

These had no earlier name.

| Function | Purpose |
|---|---|
| `as_net_frame()` | convert other objects into a net_frame (generic) |
| `net_diagnostics()` | the validation errors and notes as a data frame |
| `net_simplify()` | collapse parallel edges, optionally drop loops |
| `net_is_simple()` | test for loops and parallel edges |
| `net_multiplicity()` | size of each edge's group of parallel edges |
| `net_is_bipartite()` | test 2-colourability (a graph property) |
| `net_is_spatial()` | test for the spatial class |
| `net_split_relations()` | one network per relation of a multiplex network |
| `net_panel()` | a panel (waves) from a list of networks |
| `net_waves()` | split a panel back into one network per wave |
| `net_node_spells()` | node presence spells |
| `net_nodes_at()` | node attributes in force at a time |
| `net_as_signed()` | split a signed value column into sign and magnitude |
| `net_resolve_missing()` | resolve unobserved pairs and non-respondents explicitly |
| `net_as_multilayer()` | lift a multiplex network to layer form |
| `net_is_multilayer()` | test for the multilayer class |
| `net_edge_layers()` | classify edges as intralayer, coupling, or interlayer |
| `net_split_layers()` | one network per layer |
| `net_ego_frame()` | ego networks of a whole network in one net_frame |
| `net_split_egos()` | one network per ego |
| `net_to_egor()` | convert ego networks to egor's three tables |
| `net_plot_layers()` | one panel per layer |
| `net_plot_snapshots()` | one panel per time window or wave |
| `net_plot_egos()` | one panel per ego |
