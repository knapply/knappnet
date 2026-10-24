/* Compiles a C program against netcapi.h and the static library, checking
 * that the header matches the exported ABI. Run by tools/check_c_abi.sh. */
#include <math.h>
#include <stdio.h>
#include "netcapi.h"

/* Flows and cuts on the directed network 0 -> 1 -> 3, 0 -> 2 -> 3. */
static int flow_smoke(void) {
  char err[256];
  int64_t src[] = {0, 0, 1, 2};
  int64_t tgt[] = {1, 2, 3, 3};
  double cap[] = {2.0, 1.0, 1.0, INFINITY};
  GraphHandle *g = NULL;
  if (net_graph_new(4, src, tgt, NULL, 4, 1, &g, err, sizeof err) != 0) return 1;
  int64_t s[] = {0}, t[] = {3};
  double value = 0.0, flow[4], cost = 0.0, weight[4];
  int32_t side[4], cut[4];
  if (net_max_flow(g, s, 1, t, 1, cap, &value, flow, err, sizeof err) != 0 || value != 2.0)
    return 1;
  if (net_min_cut(g, s, 1, t, 1, cap, 1, &value, side, cut, err, sizeof err) != 0 ||
      side[1] != 1 || cut[2] != 1)
    return 1;
  if (net_min_node_cut_capacitated(g, s, 1, t, 1, NULL, NULL, 0, &value, side, NULL, err,
                                   sizeof err) != 0 || value != 2.0)
    return 1;
  double demand[] = {-2.0, 0.0, 0.0, 2.0};
  if (net_min_cost_flow(g, demand, cap, NULL, &cost, flow, err, sizeof err) != 0) return 1;
  if (net_max_flow_min_cost(g, 0, 3, cap, NULL, &value, &cost, NULL, err, sizeof err) != 0)
    return 1;
  int64_t parent[4];
  /* Gomory-Hu trees need an undirected network: this must fail. */
  if (net_gomory_hu_tree(g, NULL, parent, weight, err, sizeof err) == 0) return 1;
  net_graph_free(g);
  return 0;
}

/* Connectivity on two triangles 0-1-2 and 3-4-5 joined by 2-3. */
static int connectivity_smoke(void) {
  char err[256];
  int64_t src[] = {0, 1, 2, 3, 4, 5, 2};
  int64_t tgt[] = {1, 2, 0, 4, 5, 3, 3};
  GraphHandle *g = NULL;
  if (net_graph_new(6, src, tgt, NULL, 7, 0, &g, err, sizeof err) != 0) return 1;
  int64_t value = 0, size = 0, from[7], to[7], nodes[6];
  if (net_local_connectivity(g, 0, 5, 0, &value, err, sizeof err) != 0 || value != 1) return 1;
  if (net_connectivity(g, 1, 1, &value, err, sizeof err) != 0 || value != 1) return 1;
  double pairs[36], average = 0.0;
  if (net_all_pairs_node_connectivity(g, 1, pairs, &average, err, sizeof err) != 0) return 1;
  if (net_min_edge_cut(g, -1, -1, 1, from, to, &size, err, sizeof err) != 0 || size != 1)
    return 1;
  if (net_min_node_cut(g, 0, 5, 1, 1, nodes, &size, err, sizeof err) != 0 || size != 1)
    return 1;
  ListsHandle *paths = NULL;
  if (net_disjoint_paths(g, 0, 5, 1, &paths, err, sizeof err) != 0) return 1;
  if (net_lists_count(paths) != 1) return 1;
  int64_t offsets[2], items[8];
  double values[1];
  if (net_lists_total(paths) > 8) return 1;
  if (net_lists_copy(paths, offsets, items, values, err, sizeof err) != 0) return 1;
  net_lists_free(paths);
  ListsHandle *kc = NULL;
  if (net_k_components(g, 1, &kc, err, sizeof err) != 0 || net_lists_count(kc) != 3) return 1;
  net_lists_free(kc);
  int32_t labels[7], flags[7];
  double span[7];
  if (net_k_edge_components(g, 2, 1, labels, err, sizeof err) != 0 || labels[3] != 1) return 1;
  if (net_local_bridges(g, flags, span, err, sizeof err) != 0 || flags[6] != 1) return 1;
  int32_t semi = 0;
  /* Directed-only functions refuse an undirected network. */
  if (net_is_semiconnected(g, &semi, err, sizeof err) == 0) return 1;
  if (net_condensation(g, 0, labels, &value, err, sizeof err) == 0) return 1;
  if (net_attracting_components(g, 1, labels, err, sizeof err) == 0) return 1;
  net_graph_free(g);
  return 0;
}

/* Trees and DAGs: a weighted square with a diagonal, and a diamond DAG. */
static int trees_smoke(void) {
  char err[256];
  int64_t src[] = {0, 1, 2, 3, 0};
  int64_t tgt[] = {1, 2, 3, 0, 2};
  double w[] = {1.0, 2.0, 1.0, 2.0, 1.0};
  GraphHandle *u = NULL, *d = NULL;
  if (net_graph_new(4, src, tgt, w, 5, 0, &u, err, sizeof err) != 0) return 1;
  int64_t edges[4], size = 0, nodes[4], gen[4], path[4], dom[4];
  double weight = 0.0, length = 0.0;
  if (net_spanning_tree(u, 1, 0, 1, edges, &size, &weight, err, sizeof err) != 0 ||
      size != 3 || weight != 3.0)
    return 1;
  int64_t ds[] = {0, 0, 1, 2};
  int64_t dt[] = {1, 2, 3, 3};
  if (net_graph_new(4, ds, dt, NULL, 4, 1, &d, err, sizeof err) != 0) return 1;
  if (net_branching(d, 0, 3, -1, edges, &size, &weight, err, sizeof err) != 0 || size != 3)
    return 1;
  if (net_topological_sort(d, 1, nodes, gen, err, sizeof err) != 0 || gen[3] != 2) return 1;
  if (net_lexicographical_topological_sort(d, nodes, err, sizeof err) != 0) return 1;
  ListsHandle *l = NULL;
  if (net_all_topological_sorts(d, 1, 0, &l, err, sizeof err) != 0 || net_lists_count(l) != 2)
    return 1;
  net_lists_free(l);
  int32_t dag = 0;
  if (net_is_dag(d, &dag, err, sizeof err) != 0 || dag != 1) return 1;
  if (net_transitive_closure(d, 0, &l, err, sizeof err) != 0 || net_lists_total(l) != 5)
    return 1;
  net_lists_free(l);
  if (net_transitive_reduction(d, &l, err, sizeof err) != 0) return 1;
  net_lists_free(l);
  if (net_dag_longest_path(d, 0, 0, path, &size, &length, err, sizeof err) != 0 ||
      length != 2.0)
    return 1;
  if (net_immediate_dominators(d, 0, dom, err, sizeof err) != 0 || dom[3] != 0) return 1;
  if (net_dominance_frontiers(d, 0, &l, err, sizeof err) != 0) return 1;
  net_lists_free(l);
  int64_t a[] = {1}, b[] = {2}, anc[1];
  if (net_lowest_common_ancestors(d, a, b, 1, -2, 0, anc, err, sizeof err) != 0 || anc[0] != 0)
    return 1;
  net_graph_free(u);
  net_graph_free(d);
  return 0;
}

/* Cycles of a 4-cycle with a chord. */
static int cycles_smoke(void) {
  char err[256];
  int64_t src[] = {0, 1, 2, 3, 0};
  int64_t tgt[] = {1, 2, 3, 0, 2};
  GraphHandle *g = NULL;
  if (net_graph_new(4, src, tgt, NULL, 5, 0, &g, err, sizeof err) != 0) return 1;
  ListsHandle *l = NULL;
  if (net_cycle_basis(g, 1, -1, &l, err, sizeof err) != 0 || net_lists_count(l) != 2) return 1;
  net_lists_free(l);
  if (net_minimum_cycle_basis(g, 0, 1, &l, err, sizeof err) != 0 || net_lists_count(l) != 2)
    return 1;
  net_lists_free(l);
  if (net_simple_cycles(g, -1, 0, &l, err, sizeof err) != 0 || net_lists_count(l) != 3) return 1;
  net_lists_free(l);
  if (net_chordless_cycles(g, -1, 0, &l, err, sizeof err) != 0 || net_lists_count(l) != 2)
    return 1;
  net_lists_free(l);
  double gi = 0.0;
  if (net_girth(g, &gi, err, sizeof err) != 0 || gi != 3.0) return 1;
  int64_t edges[5], size = 0;
  int32_t rev[5];
  if (net_find_cycle(g, -1, 0, 1, edges, rev, &size, err, sizeof err) != 0 || size != 4)
    return 1;
  net_graph_free(g);
  return 0;
}

/* Link prediction, Weisfeiler-Lehman hashes, SimRank, Panther, and edit
 * distance on a path 0-1-2-3. */
static int similarity_checks(void) {
  int64_t src[] = {0, 1, 2};
  int64_t tgt[] = {1, 2, 3};
  char err[256];
  GraphHandle *g = NULL;
  if (net_graph_new(4, src, tgt, NULL, 3, 0, &g, err, sizeof err) != 0) return 1;
  int64_t us[] = {0}, vs[] = {2}, count = 0;
  double score = 0;
  if (net_link_prediction(g, 0, 0.8, 0.001, NULL, us, vs, 1, 1, &score, err, sizeof err) != 0 ||
      score != 1.0) return 1;
  if (net_non_edge_count(g, &count, err, sizeof err) != 0 || count != 3) return 1;
  int64_t a[3], b[3];
  double sc[3];
  if (net_link_prediction_non_edges(g, 1, 0.8, 0.001, NULL, 1, a, b, sc, 3, err,
                                    sizeof err) != 0 || a[0] != 0 || b[0] != 2) return 1;
  char hash[32];
  if (net_wl_graph_hash(g, NULL, NULL, 3, 16, 1, hash, sizeof hash, err, sizeof err) != 0)
    return 1;
  char sub[4 * 3 * 32];
  if (net_wl_subgraph_hashes(g, NULL, NULL, 3, 16, 0, 1, sub, sizeof sub, err,
                             sizeof err) != 0) return 1;
  double s[16];
  if (net_simrank(g, 0.9, 1000, 1e-4, 0, 0, 1, s, err, sizeof err) != 0 || s[0] != 1.0)
    return 1;
  int64_t nodes[4], found = 0;
  double ps[4];
  if (net_panther_similarity(g, 1, 2, 5, 0.5, 0.1, 0.2, 0, 1, 1, nodes, ps, 4, &found, err,
                             sizeof err) != 0 || found < 2) return 1;
  double d = -1;
  int32_t ok = 0;
  if (net_graph_edit_distance(g, g, NULL, NULL, NULL, NULL, NULL, -1, -1, INFINITY, 0, &d,
                              &ok, err, sizeof err) != 0 || !ok || d != 0.0) return 1;
  net_graph_free(g);
  return 0;
}

/* Isomorphism, trusses, and onion layers on a 4-cycle. */
static int structure_checks(void) {
  int64_t src[] = {0, 1, 2, 3};
  int64_t tgt[] = {1, 2, 3, 0};
  char err[256];
  GraphHandle *g = NULL;
  if (net_graph_new(4, src, tgt, NULL, 4, 0, &g, err, sizeof err) != 0) return 1;
  int64_t count = 0;
  if (net_vf2pp_count(g, g, NULL, NULL, 0, 0, 100, &count, err, sizeof err) != 0 ||
      count != 8) return 1;
  MatcherHandle *m = NULL;
  if (net_vf2pp_new(g, g, NULL, NULL, 0, 0, &m, err, sizeof err) != 0) return 1;
  int64_t maps[4 * 8], found = 0;
  if (net_vf2pp_next(m, 8, maps, &found, err, sizeof err) != 0 || found != 8) return 1;
  net_vf2pp_free(m);
  net_vf2pp_free(NULL);
  int32_t ok = 0;
  if (net_could_be_isomorphic(g, g, 7, 1, &ok, err, sizeof err) != 0 || !ok) return 1;
  int64_t tsrc[] = {0, 1, 2}, ttgt[] = {1, 2, 3}, map[4];
  GraphHandle *t = NULL;
  if (net_graph_new(4, tsrc, ttgt, NULL, 3, 0, &t, err, sizeof err) != 0) return 1;
  if (net_tree_isomorphism(t, t, -1, -1, map, &ok, err, sizeof err) != 0 || !ok) return 1;
  if (net_tree_isomorphism(g, g, -1, -1, map, &ok, err, sizeof err) == 0) return 1;
  int64_t nodes[4], truss[4], edges[4], ne = 0, nn = 0, layers[4];
  if (net_core_subset(g, 1, 0, 0, nodes, &count, err, sizeof err) != 0 || count != 4) return 1;
  if (net_truss_numbers(g, truss, err, sizeof err) != 0 || truss[0] != 2) return 1;
  if (net_k_truss(g, 3, edges, &ne, nodes, &nn, err, sizeof err) != 0 || ne != 0) return 1;
  if (net_onion_layers(g, layers, err, sizeof err) != 0 || layers[0] != 1) return 1;
  net_graph_free(t);
  net_graph_free(g);
  return 0;
}

/* Mixing, cuts, matchings, colourings, and squares on a path 0-1-2-3. */
static int mixing_checks(void) {
  int64_t src[] = {0, 1, 2};
  int64_t tgt[] = {1, 2, 3};
  char err[256];
  GraphHandle *g = NULL;
  if (net_graph_new(4, src, tgt, NULL, 3, 0, &g, err, sizeof err) != 0) return 1;
  double vals[4], counts[16], r = 0, avg[4], ks[4], kv[4], cut[9], sq[4];
  int64_t k = 0, nk = 0;
  if (net_degree_mixing(g, 1, 2, 0, NULL, 0, 4, vals, counts, &k, err, sizeof err) != 0 ||
      k != 2) return 1;
  double labels[] = {0, 1, 0, 1};
  if (net_attribute_mixing(g, labels, NULL, 0, 4, vals, counts, &k, err, sizeof err) != 0)
    return 1;
  if (net_assortativity_subset(g, 0, 1, 2, 0, NULL, NULL, 0, &r, err, sizeof err) != 0)
    return 1;
  if (net_average_neighbor_degree(g, 1, 1, 0, NULL, 0, avg, err, sizeof err) != 0 ||
      avg[0] != 2.0) return 1;
  if (net_average_degree_connectivity(g, 3, 3, 0, NULL, 0, ks, kv, &nk, err, sizeof err) != 0 ||
      nk != 2) return 1;
  int64_t s[] = {0, 1};
  if (net_cut_measures(g, s, 2, NULL, 0, 0, cut, err, sizeof err) != 0 || cut[0] != 1.0)
    return 1;
  int64_t edges[3], found = 0, colours[4], nodes[4];
  if (net_max_weight_matching(g, 0, 0, edges, &found, err, sizeof err) != 0 || found != 2)
    return 1;
  if (net_min_weight_matching(g, 0, edges, &found, err, sizeof err) != 0 || found != 2)
    return 1;
  if (net_maximal_matching(g, edges, &found, err, sizeof err) != 0 || found != 2) return 1;
  int64_t us[] = {0, 2}, vs[] = {1, 3};
  int32_t chk[3];
  if (net_check_matching(g, us, vs, 2, chk, err, sizeof err) != 0 || !chk[2]) return 1;
  if (net_greedy_color(g, 6, 0, colours, err, sizeof err) != 0) return 1;
  if (net_maximal_independent_set(g, NULL, 0, 1, nodes, &found, err, sizeof err) != 0)
    return 1;
  if (net_square_clustering(g, 1, sq, err, sizeof err) != 0) return 1;
  int64_t off[5], gk[6], gc[6], total = 0;
  if (net_generalized_degree(g, 1, off, gk, gc, 6, &total, err, sizeof err) != 0 ||
      total != 4) return 1;
  net_graph_free(g);
  return 0;
}

/* Swaps, reference graphs, small-world coefficients, and the rich club
 * on a 6-cycle with two chords. */
static int rewiring_checks(void) {
  int64_t src[] = {0, 1, 2, 3, 4, 5, 0, 1};
  int64_t tgt[] = {1, 2, 3, 4, 5, 0, 3, 4};
  char err[256];
  GraphHandle *g = NULL, *d = NULL;
  if (net_graph_new(6, src, tgt, NULL, 8, 0, &g, err, sizeof err) != 0) return 1;
  if (net_graph_new(6, src, tgt, NULL, 8, 1, &d, err, sizeof err) != 0) return 1;
  int64_t from[8], to[8], m = 0, swaps = 0;
  if (net_double_edge_swap(g, 2, 1000, 1, from, to, &m, err, sizeof err) != 0 || m != 8)
    return 1;
  if (net_connected_double_edge_swap(g, 5, 1, from, to, &m, &swaps, err, sizeof err) != 0)
    return 1;
  if (net_directed_edge_swap(d, 1, 1000, 1, from, to, &m, err, sizeof err) != 0 || m != 8)
    return 1;
  if (net_random_reference(g, 1, 1, 1, from, to, &m, err, sizeof err) != 0) return 1;
  if (net_lattice_reference(g, 1, NULL, 1, 1, from, to, &m, err, sizeof err) != 0) return 1;
  double s = 0, w = 0, rc[6];
  int64_t k = 0;
  if (net_small_world_sigma(g, 1, 2, 1, 1, &s, err, sizeof err) != 0) return 1;
  if (net_small_world_omega(g, 1, 2, 1, 1, &w, err, sizeof err) != 0) return 1;
  if (net_rich_club(g, 0, 10, 1, 1, 1, rc, &k, err, sizeof err) != 0 || k < 1) return 1;
  net_graph_free(d);
  net_graph_free(g);
  return 0;
}

int main(void) {
  int64_t src[] = {0, 1, 2};
  int64_t tgt[] = {1, 2, 3};
  char err[256];
  GraphHandle *g = NULL;
  if (net_graph_new(4, src, tgt, NULL, 3, 0, &g, err, sizeof err) != 0) {
    fprintf(stderr, "graph_new failed: %s\n", err);
    return 1;
  }
  double bc[4], deg[4], census[16];
  if (net_betweenness(g, 0, 0, 0, 1, bc, err, sizeof err) != 0) return 1;
  if (net_degree(g, 3, 1, 0, deg, err, sizeof err) != 0) return 1;
  if (net_triad_census(g, 1, census, err, sizeof err) != 0) return 1;
  if (bc[1] != 2.0 || bc[2] != 2.0 || deg[0] != 1.0 || deg[1] != 2.0) return 1;
  /* Errors must come back as messages, not crashes. */
  int64_t bad[] = {9};
  GraphHandle *h = NULL;
  if (net_graph_new(2, bad, bad, NULL, 1, 1, &h, err, sizeof err) == 0) return 1;
  if (net_degree(g, 7, 1, 0, deg, err, sizeof err) == 0) return 1;
  net_graph_free(g);
  net_graph_free(NULL);

  /* Files: build a network, write it, read it back. */
  NetworkHandle *n = net_network_new(1);
  const char *ids[] = {"a", "b", "c"};
  const char *from[] = {"a", "b"};
  const char *to[] = {"b", "c"};
  double w[] = {1.5, 2.0};
  if (net_table_add_string(n, 1, "id", ids, 3, err, sizeof err) != 0) return 1;
  if (net_table_add_string(n, 2, "from", from, 2, err, sizeof err) != 0) return 1;
  if (net_table_add_string(n, 2, "to", to, 2, err, sizeof err) != 0) return 1;
  if (net_table_add_numeric(n, 2, "weight", 2, w, NULL, 2, err, sizeof err) != 0) return 1;
  const char *path = "netcapi-smoke.graphml";
  if (net_write(n, path, NULL, NULL, err, sizeof err) != 0) {
    fprintf(stderr, "write failed: %s\n", err);
    return 1;
  }
  net_network_free(n);
  NetworkHandle *r = NULL;
  if (net_read(path, NULL, -1, &r, err, sizeof err) != 0) {
    fprintf(stderr, "read failed: %s\n", err);
    return 1;
  }
  remove(path);
  if (net_table_nrow(r, 1) != 3 || net_table_nrow(r, 2) != 2) return 1;
  if (net_table_kind(r, 2, 2) != 2) return 1;
  double back[2];
  int32_t na[2];
  if (net_table_numeric(r, 2, 2, back, na, err, sizeof err) != 0 || back[1] != 2.0) return 1;
  const char *s = net_table_string(r, 1, 0, 2);
  if (s == NULL || s[0] != 'c') return 1;
  net_network_free(r);
  if (net_read("missing.graphml", NULL, -1, &r, err, sizeof err) == 0) return 1;

  /* Shortest paths: a weighted path graph and a k-shortest-paths query. */
  int64_t ps[] = {0, 1, 0};
  int64_t pt[] = {1, 2, 2};
  double pw[] = {1.0, 1.0, 5.0};
  GraphHandle *pgr = NULL;
  PathGraphHandle *pgh = NULL;
  if (net_graph_new(3, ps, pt, pw, 3, 1, &pgr, err, sizeof err) != 0) return 1;
  if (net_path_graph_new(pgr, 1, 1, 0, NULL, -1, 0, &pgh, err, sizeof err) != 0) return 1;
  int64_t src0[] = {0};
  double pdist[3];
  int64_t ppar[3];
  if (net_shortest_paths(pgh, src0, 1, INFINITY, -1, 0, pdist, ppar, NULL, NULL, err,
                         sizeof err) != 0) return 1;
  if (pdist[2] != 2.0 || ppar[2] != 1 || ppar[0] != -1) return 1;
  PathListsHandle *kp = NULL;
  if (net_k_shortest_paths(pgh, 0, 2, 5, &kp, err, sizeof err) != 0) return 1;
  if (net_path_lists_count(kp) != 2 || net_path_lists_total(kp) != 5) return 1;
  int64_t koff[3], knodes[5];
  double klen[2];
  if (net_path_lists_copy(kp, koff, knodes, klen, err, sizeof err) != 0) return 1;
  if (klen[0] != 2.0 || klen[1] != 5.0 || koff[1] != 3 || knodes[4] != 2) return 1;
  net_path_lists_free(kp);
  net_path_graph_free(pgh);
  net_path_graph_free(NULL);
  net_graph_free(pgr);
  if (flow_smoke() != 0) {
    fprintf(stderr, "flow functions failed\n");
    return 1;
  }
  if (connectivity_smoke() != 0) {
    fprintf(stderr, "connectivity functions failed\n");
    return 1;
  }
  if (trees_smoke() != 0) {
    fprintf(stderr, "tree and DAG functions failed\n");
    return 1;
  }
  if (cycles_smoke() != 0) {
    fprintf(stderr, "cycle functions failed\n");
    return 1;
  }
  if (rewiring_checks() != 0) {
    fprintf(stderr, "rewiring checks failed\n");
    return 1;
  }
  if (mixing_checks() != 0) {
    fprintf(stderr, "mixing checks failed\n");
    return 1;
  }
  if (structure_checks() != 0) {
    fprintf(stderr, "structure checks failed\n");
    return 1;
  }
  if (similarity_checks() != 0) {
    fprintf(stderr, "similarity checks failed\n");
    return 1;
  }
  puts("C ABI smoke test passed");
  return 0;
}
