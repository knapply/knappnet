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
  if (flow_smoke() != 0) {
    fprintf(stderr, "flow functions failed\n");
    return 1;
  }
  puts("C ABI smoke test passed");
  return 0;
}
