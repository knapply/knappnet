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
  if (similarity_checks() != 0) {
    fprintf(stderr, "similarity checks failed\n");
    return 1;
  }
  puts("C ABI smoke test passed");
  return 0;
}
