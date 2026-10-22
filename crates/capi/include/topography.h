/* ---- Topography and centrality ----------------------------------------
 * Codes: mode 1 = out, 2 = in, 3 = all; weights_as 0 = distance,
 * 1 = strength; kind 0 = weak, 1 = strong; unreachable 0 = ignore (igraph),
 * 1 = Wasserman-Faust (NetworkX), 2 = zero (sna), 3 = longest length plus 1
 * (UCINET); measure 0 = degree, 1 = betweenness, 2 = closeness;
 * convention 0 = igraph, 1 = NetworkX, 2 = sna, 3 = UCINET; scale
 * 0 = largest score 1, 1 = unit Euclidean norm; eigenvector loops
 * 0 = twice, 1 = once, 2 = ignore. threads == 0 uses all cores. */

/* The adjacency-matrix reading of a graph (sna, UCINET): parallel edges
 * collapse, keeping the largest (keep == 0) or smallest (keep == 1)
 * weight, and self-loops are dropped. undirected != 0 also symmetrizes a
 * directed graph. Release *out with net_graph_free. */
int32_t net_graph_simplify(const GraphHandle *h, int32_t undirected,
                           int32_t keep, GraphHandle **out, char *err,
                           size_t err_len);

/* Edges divided by possible edges; loops != 0 counts possible self-loops.
 * Writes one double. */
int32_t net_density(const GraphHandle *h, int32_t loops, double *out,
                    char *err, size_t err_len);

/* membership: n labels, numbered from 0 in order of each component's
 * lowest-numbered node. count: one value. */
int32_t net_components(const GraphHandle *h, int32_t kind, int32_t *membership,
                       int64_t *count, char *err, size_t err_len);

/* Share of ordered pairs of distinct nodes in different components. */
int32_t net_fragmentation(const GraphHandle *h, int32_t kind, double *out,
                          char *err, size_t err_len);

/* Freeman centralization; writes one double. convention 1 (NetworkX) is
 * an error, since NetworkX has no centralization index. */
int32_t net_centralization(const GraphHandle *h, int32_t measure, int32_t mode,
                           int32_t loops, int32_t normalized,
                           int32_t convention, int32_t threads, double *out,
                           char *err, size_t err_len);

/* Closeness (harmonic == 0) or harmonic centrality (harmonic != 0); out
 * holds n doubles. unreachable only matters for closeness. */
int32_t net_closeness(const GraphHandle *h, int32_t harmonic, int32_t mode,
                      int32_t weighted, int32_t weights_as, int32_t normalized,
                      int32_t unreachable, int32_t threads, double *out,
                      char *err, size_t err_len);

/* Eigenvector centrality by power iteration; out holds n doubles.
 * tol <= 0 and max_iter <= 0 select the defaults (1e-12, 10000).
 * *connected is 1 if the graph is (strongly) connected, else 0. */
int32_t net_eigenvector(const GraphHandle *h, int32_t mode, int32_t weighted,
                        int32_t scale, int32_t loops, double tol,
                        int64_t max_iter, int32_t threads, double *out,
                        double *eigenvalue, int32_t *connected, char *err,
                        size_t err_len);
