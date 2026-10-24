/* Compiles a C program against netcapi.h and the static library, checking
 * that the header matches the exported ABI. Run by tools/check_c_abi.sh. */
#include <math.h>
#include <stdio.h>
#include "netcapi.h"

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
