/* ---- Descriptive and structural measures ------------------------------
 * Codes: mode 1 = out, 2 = in, 3 = all; weights_as 0 = distance,
 * 1 = strength; reciprocity mode 0 = default (reciprocated arcs over
 * arcs), 1 = ratio (mutual over non-null dyads), 2 = dyadic (symmetric
 * over all dyads), 3 = correlation (Garlaschelli-Loffredo rho);
 * weighting 0 = binary, 1 = Barrat, 2 = Onnela; undefined 0 = NaN
 * (igraph), 1 = zero (NetworkX). sources are 0-based node indices; a NULL
 * pointer or n_sources < 0 means every node. threads == 0 uses all
 * cores. */

/* Edge count (parallel edges and loops included), or -1 for NULL. */
int64_t net_graph_edge_count(const GraphHandle *h);

/* Reciprocity on the adjacency matrix (parallel arcs once, loops
 * ignored); writes one double, NaN when the denominator is 0. */
int32_t net_reciprocity(const GraphHandle *h, int32_t mode, double *out,
                        char *err, size_t err_len);

/* Global transitivity of the simple undirected reading: closed over
 * connected triples. Writes one double. */
int32_t net_transitivity(const GraphHandle *h, int32_t undefined,
                         int32_t threads, double *out, char *err,
                         size_t err_len);

/* Local clustering coefficients (average == 0; out holds n doubles) or
 * their mean (average != 0; one double). With undefined == 0 the mean
 * leaves out nodes with fewer than two neighbors. Weighted coefficients
 * need positive weights. */
int32_t net_clustering(const GraphHandle *h, int32_t average,
                       int32_t weighting, int32_t undefined, int32_t threads,
                       double *out, char *err, size_t err_len);

/* Shortest-path distances as a column-major matrix with one row per
 * source and one column per node; unreachable nodes are at distance
 * infinity. out holds (number of sources) * n doubles. */
int32_t net_distances(const GraphHandle *h, const int64_t *sources,
                      int64_t n_sources, int32_t mode, int32_t weighted,
                      int32_t weights_as, int32_t threads, double *out,
                      char *err, size_t err_len);

/* Totals over the searches from the sources; out holds 7 doubles: sources,
 * ordered pairs joined by a path, sum of their distances, longest finite
 * distance, sources not reaching every other node, mean distance, and
 * diameter (both infinite when unconnected == 0 and some source does not
 * reach every other node). */
int32_t net_path_summary(const GraphHandle *h, const int64_t *sources,
                         int64_t n_sources, int32_t mode, int32_t weighted,
                         int32_t weights_as, int32_t unconnected,
                         int32_t threads, double *out, char *err,
                         size_t err_len);

/* Largest distance from each node to the nodes it reaches (0 if none);
 * out holds n doubles. */
int32_t net_eccentricity(const GraphHandle *h, int32_t mode, int32_t weighted,
                         int32_t weights_as, int32_t threads, double *out,
                         char *err, size_t err_len);

/* Number of other nodes each node reaches (mode 1), is reached from (2),
 * or is joined to ignoring direction (3); out holds n doubles. */
int32_t net_reachability(const GraphHandle *h, int32_t mode, int32_t threads,
                         double *out, char *err, size_t err_len);

/* Cut points (cut_point: n flags) and bridges (bridge: one flag per edge,
 * in input order) of the undirected reading in which parallel edges and
 * opposite arcs stay separate edges and loops are ignored. */
int32_t net_cut_structure(const GraphHandle *h, int32_t *cut_point,
                          int32_t *bridge, char *err, size_t err_len);

/* Writes 3 doubles: self-loops, edges repeating an earlier edge between the
 * same (ordered, if directed) pair of nodes, and nodes without an edge to
 * another node. */
int32_t net_edge_counts(const GraphHandle *h, double *out, char *err,
                        size_t err_len);
