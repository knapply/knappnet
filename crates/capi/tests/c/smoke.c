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
  if (connectivity_smoke() != 0) {
    fprintf(stderr, "connectivity functions failed\n");
    return 1;
  }
  puts("C ABI smoke test passed");
  return 0;
}
