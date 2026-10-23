/*
 * C interface to the network analysis engine.
 *
 * All functions return 0 on success. On failure they return non-zero and
 * write a NUL-terminated message into `err` (at most `err_len` bytes).
 * Node indices are 0-based. Output buffers are allocated by the caller.
 */
#ifndef NETCAPI_H
#define NETCAPI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GraphHandle GraphHandle;

int32_t net_graph_new(int64_t n, const int64_t *src, const int64_t *tgt,
                      const double *weights, int64_t m, int32_t directed,
                      GraphHandle **out, char *err, size_t err_len);

void net_graph_free(GraphHandle *h);

int64_t net_graph_node_count(const GraphHandle *h);

/* mode: 1 = out, 2 = in, 3 = all. weighted != 0 returns strength. */
int32_t net_degree(const GraphHandle *h, int32_t mode, int32_t loops,
                   int32_t weighted, double *out, char *err, size_t err_len);

/* weights_as: 0 = distance, 1 = strength (lengths 1/w); only used when
 * weighted != 0. threads == 0 uses all available cores. */
int32_t net_betweenness(const GraphHandle *h, int32_t weighted,
                        int32_t weights_as, int32_t normalized,
                        int32_t threads, double *out, char *err,
                        size_t err_len);

/* Writes 16 counts in MAN order: 003 012 102 021D 021U 021C 111D 111U
 * 030T 030C 201 120D 120U 120C 210 300. */
int32_t net_triad_census(const GraphHandle *h, int32_t threads, double *out,
                         char *err, size_t err_len);

/* ---- Files ----------------------------------------------------------
 * A network handle holds three tables: 0 = graph attributes (one row),
 * 1 = nodes (first column "id"), 2 = edges (first columns "from", "to").
 * Column kinds: 0 = logical, 1 = integer, 2 = double, 3 = text.
 * Formats: graphml, gexf, gml, pajek, dl, json, dynetml (ORA), dot (write
 * only). Reading a .xml file without a format chooses GraphML or DyNetML
 * from the root element; writing one produces GraphML. */

typedef struct NetworkHandle NetworkHandle;

/* format may be NULL to infer from the extension. directed is -1 (not
 * specified), 0, or 1; it only matters for formats without direction. */
int32_t net_read(const char *path, const char *format, int32_t directed,
                 NetworkHandle **out, char *err, size_t err_len);

NetworkHandle *net_network_new(int32_t directed);
void net_network_free(NetworkHandle *h);
int32_t net_network_directed(const NetworkHandle *h);

int64_t net_table_nrow(const NetworkHandle *h, int32_t which);
int32_t net_table_ncol(const NetworkHandle *h, int32_t which);
const char *net_table_name(const NetworkHandle *h, int32_t which, int32_t col);
int32_t net_table_kind(const NetworkHandle *h, int32_t which, int32_t col);

/* Logical, integer, and double columns as doubles; na[i] = 1 if missing. */
int32_t net_table_numeric(const NetworkHandle *h, int32_t which, int32_t col,
                          double *out, int32_t *na, char *err, size_t err_len);

/* NULL if missing or not a text column. Valid while the handle lives. */
const char *net_table_string(NetworkHandle *h, int32_t which, int32_t col,
                             int64_t row);

/* kind: 0 logical, 1 integer, 2 double. na may be NULL. */
int32_t net_table_add_numeric(NetworkHandle *h, int32_t which, const char *name,
                              int32_t kind, const double *vals,
                              const int32_t *na, int64_t n, char *err,
                              size_t err_len);

/* NULL entries are missing values. */
int32_t net_table_add_string(NetworkHandle *h, int32_t which, const char *name,
                             const char *const *vals, int64_t n, char *err,
                             size_t err_len);

/* format and weight may be NULL. */
int32_t net_write(const NetworkHandle *h, const char *path, const char *format,
                  const char *weight, char *err, size_t err_len);

/* ---- Cohesive subgroups ----------------------------------------------
 * Cores, cliques, and communities. Declarations to be merged into
 * netcapi.h. */

/* Coreness (k-core number) of every node, written to out[0..n). mode:
 * 1 = out, 2 = in, 3 = all; ignored for undirected graphs. */
int32_t net_coreness(const GraphHandle *h, int32_t mode, double *out,
                     char *err, size_t err_len);

/* A collection of node groups, read in two steps: query the number of
 * groups and members, allocate, then copy. Members of group i are
 * members[offsets[i] .. offsets[i + 1]), 0-based and increasing. */
typedef struct GroupsHandle GroupsHandle;

int64_t net_groups_count(const GroupsHandle *h);
int64_t net_groups_total(const GroupsHandle *h);
/* offsets must hold count + 1 values and members total values. */
int32_t net_groups_copy(const GroupsHandle *h, int64_t *offsets,
                        int64_t *members, char *err, size_t err_len);
void net_groups_free(GroupsHandle *h);

/* Maximal cliques of the undirected simple graph underlying the network
 * (directions, loops, and edge multiplicity ignored) with at least
 * min_size members, in lexicographic order. Release the result with
 * net_groups_free. threads == 0 uses all available cores. */
int32_t net_max_cliques(const GraphHandle *h, int64_t min_size,
                        int32_t threads, GroupsHandle **out, char *err,
                        size_t err_len);

/* Size of the largest clique (0 for a graph with no nodes). */
int32_t net_clique_number(const GraphHandle *h, int32_t threads,
                          int64_t *out, char *err, size_t err_len);

/* Modularity of a partition: membership holds n non-negative labels.
 * Weights (weighted != 0) must be non-negative. Writes NaN when the
 * network has no edges or total weight zero. */
int32_t net_modularity(const GraphHandle *h, const int64_t *membership,
                       int64_t n, double resolution, int32_t weighted,
                       double *out, char *err, size_t err_len);

/* Louvain community detection on an undirected graph. Writes 0-based
 * group labels, numbered by first appearance, to membership[0..n) and the
 * partition's modularity to *q. The same seed gives the same result. */
int32_t net_louvain(const GraphHandle *h, double resolution, int32_t weighted,
                    uint64_t seed, int64_t *membership, double *q, char *err,
                    size_t err_len);

/*
 * Brokerage and positions: Burt's structural holes, Gould-Fernandez
 * brokerage, structural equivalence, and CONCOR. To be merged into
 * netcapi.h.
 */

/* measure: 0 = constraint, 1 = effective size, 2 = efficiency. With
 * weighted != 0, weights are tie strengths (strictly positive). Nodes
 * without contacts get NaN. threads == 0 uses all available cores. */
int32_t net_structural_holes(const GraphHandle *h, int32_t measure,
                             int32_t weighted, int32_t threads, double *out,
                             char *err, size_t err_len);

/* groups holds one code per node (n_groups must equal n). out receives an
 * n x 6 column-major matrix: coordinator, itinerant, representative,
 * gatekeeper, liaison, total. */
int32_t net_brokerage(const GraphHandle *h, const int64_t *groups,
                      int64_t n_groups, int32_t threads, double *out,
                      char *err, size_t err_len);

/* method: 0 = Euclidean, 1 = correlation, 2 = Hamming. out receives an
 * n x n column-major matrix. */
int32_t net_structural_equivalence(const GraphHandle *h, int32_t method,
                                   int32_t weighted, int32_t threads,
                                   double *out, char *err, size_t err_len);

/* out receives n block numbers starting at 1. diagnostics may be NULL; it
 * receives two counts: splits that hit max_iter without converging, and
 * splits in which a correlation with the reference node came within 1e-6
 * of zero (so rounding may have decided a node's block). */
int32_t net_concor(const GraphHandle *h, int32_t depth, int32_t weighted,
                   double tolerance, int32_t max_iter, int32_t threads,
                   int32_t *out, int32_t *diagnostics, char *err,
                   size_t err_len);

/* ---- Quadratic assignment procedure (QAP) ----------------------------
 * Matrices are n-by-n doubles in column-major order (as R stores them);
 * NaN marks a missing value and drops the dyad (listwise deletion).
 * mode: 0 = digraph (every ordered pair), 1 = graph (each unordered pair,
 * read from the lower triangle). diag != 0 also uses the diagonal.
 * With exhaustive != 0 all n! relabelings are evaluated (n <= 10) and reps
 * and seed are ignored; otherwise reps random relabelings are drawn,
 * draw i from stream i of seed, so results do not depend on threads.
 * threads == 0 uses all available cores. Counts are returned as doubles. */

/* Correlation test. summary receives 7 values: observed r, dyads used,
 * draws, then the number of draws with r <= observed, r >= observed,
 * |r| >= |observed|, and with r defined. dist may be NULL; otherwise it
 * receives every draw's r (reps values, or n! when exhaustive). */
int32_t net_qap_cor(const double *y, const double *x, int64_t n, int32_t mode,
                    int32_t diag, int32_t exhaustive, int64_t reps,
                    uint64_t seed, int32_t threads, double *summary,
                    double *dist, char *err, size_t err_len);

/* Multiple regression QAP of y on k predictors stored one after another in
 * x. The design has p = k + (intercept != 0) columns, intercept first.
 * method: 0 = double semi-partialling, 1 = Y-permutation (always used when
 * p == 1). statistic: 0 = t, 1 = coefficient. tol: rank tolerance (1e-7).
 * coef receives a p-by-7 column-major table: coefficient, standard error,
 * t, then counts of draws with stat <= observed, stat >= observed,
 * |stat| >= |observed|, and with a defined statistic. fit receives 10
 * values: dyads used, residual df, residual sum of squares, model sum of
 * squares, R-squared, adjusted R-squared, residual standard error, F,
 * method used, and draws. dist may be NULL; otherwise it receives a
 * draws-by-p column-major table of every draw's statistics. */
int32_t net_mrqap(const double *y, const double *x, int64_t n, int64_t k,
                  int32_t intercept, int32_t mode, int32_t diag,
                  int32_t method, int32_t statistic, double tol,
                  int32_t exhaustive, int64_t reps, uint64_t seed,
                  int32_t threads, double *coef, double *fit, double *dist,
                  char *err, size_t err_len);

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

/* ---- Ego networks and two-mode networks -------------------------------
 * Declarations to be merged into netcapi.h, after the topography section.
 * Codes: mode 1 = out, 2 = in, 3 = all; types holds one value per node,
 * 0 for the first mode and 1 for the second; two-mode measure
 * 0 = degree, 1 = closeness, 2 = betweenness; unreachable 0 = ignore,
 * 1 = Wasserman-Faust (NetworkX); projection weight 0 = shared partners,
 * 1 = Newman's collaboration weights, 2 = sum of weight products.
 * threads == 0 uses all cores. Weights are ignored except by the product
 * projection. */

/* Number of nodes within order steps of each ego and at least mindist
 * steps away (igraph's ego_size(); mindist 0 counts the ego). With
 * egos == NULL every node is an ego and out holds n doubles; otherwise
 * egos holds k node indices and out k doubles. */
int32_t net_ego_sizes(const GraphHandle *h, const int64_t *egos, int64_t k,
                      int64_t order, int32_t mode, int64_t mindist,
                      int32_t threads, double *out, char *err,
                      size_t err_len);

/* Neighborhood members of k egos in ascending node order: those of ego i
 * are members[offsets[i] .. offsets[i + 1]). offsets holds k + 1 values,
 * members capacity values (the sum of net_ego_sizes with the same
 * arguments); a smaller capacity is an error. */
int32_t net_ego_members(const GraphHandle *h, const int64_t *egos, int64_t k,
                        int64_t order, int32_t mode, int64_t mindist,
                        int32_t threads, int64_t *offsets, int64_t *members,
                        int64_t capacity, char *err, size_t err_len);

/* Egocentric measures of every node (alters within order steps): number
 * of alters, distinct ties among them, pairs of alters (ordered if
 * directed), density (NaN without pairs), weak components among the
 * alters without ego, and ego betweenness (skipped when ego_betweenness
 * is NULL). Each output holds n doubles. */
int32_t net_ego_summary(const GraphHandle *h, int64_t order, int32_t mode,
                        int32_t threads, double *size, double *ties,
                        double *pairs, double *density, double *components,
                        double *ego_betweenness, char *err, size_t err_len);

/* Distinct ties divided by n1 * n2; writes one double (NaN if a mode is
 * empty). Fails if an edge joins two nodes of the same mode. */
int32_t net_twomode_density(const GraphHandle *h, const int32_t *types,
                            double *out, char *err, size_t err_len);

/* Degree, closeness, or betweenness; normalized != 0 applies Borgatti and
 * Everett's (1997) two-mode normalization. unreachable is used by
 * closeness only (Wasserman-Faust requires normalized != 0). out holds n
 * doubles. */
int32_t net_twomode_centrality(const GraphHandle *h, const int32_t *types,
                               int32_t measure, int32_t normalized,
                               int32_t unreachable, int32_t threads,
                               double *out, char *err, size_t err_len);

/* Leading singular vectors of the binary incidence matrix, each mode
 * scaled to Euclidean length 1 (NaN for every node without ties). tol <= 0
 * and max_iter <= 0 select the defaults (1e-12, 10000). *eigenvalue
 * receives the largest singular value, *connected 1 if the network is
 * connected. */
int32_t net_twomode_eigenvector(const GraphHandle *h, const int32_t *types,
                                double tol, int64_t max_iter, int32_t threads,
                                double *out, double *eigenvalue,
                                int32_t *connected, char *err,
                                size_t err_len);

/* A one-mode projection, read in two steps like GroupsHandle. */
typedef struct ProjectionHandle ProjectionHandle;

/* Projection onto the first (second == 0) or second mode; loops != 0
 * keeps the diagonal as self-loops (not with Newman weights). Release
 * *out with net_projection_free. */
int32_t net_project(const GraphHandle *h, const int32_t *types, int32_t second,
                    int32_t weight, int32_t loops, int32_t threads,
                    ProjectionHandle **out, char *err, size_t err_len);

/* Number of edges, or -1 for a NULL handle. */
int64_t net_projection_count(const ProjectionHandle *p);

/* Copies the edges (0-based indices of the original network, from <= to,
 * ordered by from and then to) and their weights; each array holds
 * net_projection_count(p) values. */
int32_t net_projection_copy(const ProjectionHandle *p, int64_t *from,
                            int64_t *to, double *weights, char *err,
                            size_t err_len);

void net_projection_free(ProjectionHandle *p);

/* ---- Spectral centralities and mixing patterns -------------------------
 * Codes: mode and direction 1 = out, 2 = in, 3 = all; loops 0 = an
 * undirected self-loop counts twice (igraph), 1 = once (NetworkX),
 * 2 = not at all; dangling 0 = jump by the personalization vector,
 * 1 = jump uniformly; hits scale 0 = largest score 1, 1 = sum 1, 2 = unit
 * Euclidean norm; E-I scope 0 = network, 1 = groups, 2 = nodes; E-I
 * counting 0 = Krackhardt and Stern, 1 = the R package homophily.
 * Category codes are 0-based (-1 = missing, for mixing matrices only).
 * tol <= 0 and max_iter <= 0 select the defaults (1e-12, 10000).
 * threads == 0 uses all cores. */

/* PageRank; out holds n doubles summing to 1. personalization is NULL
 * (uniform) or n nonnegative values. directed == 0 ignores directions. */
int32_t net_pagerank(const GraphHandle *h, double damping, int32_t directed,
                     int32_t weighted, const double *personalization,
                     int32_t dangling, int32_t loops, double tol,
                     int64_t max_iter, int32_t threads, double *out,
                     char *err, size_t err_len);

/* Katz centrality x = alpha A' x + beta; beta holds 1 or n values. An
 * alpha not below 1 / lambda_max is an error. normalized != 0 scales to
 * unit Euclidean norm. */
int32_t net_katz(const GraphHandle *h, double alpha, const double *beta,
                 int64_t beta_len, int32_t mode, int32_t weighted,
                 int32_t loops, int32_t normalized, double tol,
                 int64_t max_iter, int32_t threads, double *out, char *err,
                 size_t err_len);

/* Bonacich power centrality alpha (I - exponent A)^-1 A 1, scaled so the
 * squares sum to n, or (rescale != 0) so the scores sum to 1. A singular
 * system is an error. */
int32_t net_bonacich_power(const GraphHandle *h, double exponent,
                           int32_t mode, int32_t weighted, int32_t loops,
                           int32_t rescale, double tol, int64_t max_iter,
                           int32_t threads, double *out, char *err,
                           size_t err_len);

/* Hub and authority scores (n doubles each). *value receives the largest
 * eigenvalue of A^T A, *unique 1 if the leading singular value is simple. */
int32_t net_hits(const GraphHandle *h, int32_t weighted, int32_t loops,
                 int32_t scale, double tol, int64_t max_iter, int32_t threads,
                 double *hubs, double *authorities, double *value,
                 int32_t *unique, char *err, size_t err_len);

/* Nominal assortativity of n category codes in 0..k; one double (NaN if
 * undefined). normalized == 0 gives the numerator (modularity). */
int32_t net_assortativity_nominal(const GraphHandle *h, const int64_t *types,
                                  int64_t k, int32_t directed, int32_t loops,
                                  int32_t normalized, double *out, char *err,
                                  size_t err_len);

/* Assortativity of numeric values at arc sources (values) and targets
 * (values_in, or values if NULL); normalized == 0 gives the covariance. */
int32_t net_assortativity_values(const GraphHandle *h, const double *values,
                                 const double *values_in, int32_t directed,
                                 int32_t loops, int32_t normalized,
                                 double *out, char *err, size_t err_len);

/* Degree assortativity (strengths if weighted != 0); one double. */
int32_t net_assortativity_degree(const GraphHandle *h, int32_t weighted,
                                 int32_t directed, int32_t loops,
                                 int32_t normalized, double *out, char *err,
                                 size_t err_len);

/* Mixing matrix of dim1 (codes 0..k1) by dim2 (0..k2); same != 0 when they
 * are one attribute. Writes k1 * k2 counts row-major; *transposed = 1 when
 * the matrix is k2 by k1 (rows dim2, directed networks, direction 2). */
int32_t net_mixing_matrix(const GraphHandle *h, const int64_t *dim1,
                          int64_t k1, const int64_t *dim2, int64_t k2,
                          int32_t same, int32_t direction, double *out,
                          int32_t *transposed, char *err, size_t err_len);

/* E-I index for the network (1 value), groups (k), or nodes (n); external
 * and internal receive the tie counts unless NULL. NaN without ties. */
int32_t net_ei_index(const GraphHandle *h, const int64_t *groups, int64_t k,
                     int32_t scope, int32_t direction, int32_t counting,
                     double *out, double *external, double *internal,
                     char *err, size_t err_len);

/* Within-group tie shares: k rows of 5 doubles (row-major): size, internal
 * ties, external ties, observed internal share, expected share
 * (size - 1) / (n - 1). */
int32_t net_group_shares(const GraphHandle *h, const int64_t *groups,
                         int64_t k, int32_t direction, double *out,
                         char *err, size_t err_len);

#ifdef __cplusplus
}
#endif

#endif
