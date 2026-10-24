/* Compiles a C program against netcapi.h and the static library, checking
 * that the header matches the exported ABI. Run by tools/check_c_abi.sh. */
#include <math.h>
#include <stdio.h>
#include "netcapi.h"

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
  puts("C ABI smoke test passed");
  return 0;
}
