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

/*
 * Network layouts. To be merged into netcapi.h.
 *
 * Every function writes node coordinates into x and y (n doubles each).
 * Direction and self-loops are ignored and parallel edges count once.
 * weights_as: 0 = distance, 1 = strength (lengths 1/w); only used when
 * weighted != 0, and weights must then be strictly positive. Zero selects
 * the default for pivots, max_iter, iterations, and sweeps, as does a
 * value <= 0 for tol and theta. threads == 0 uses all available cores;
 * results do not depend on the thread count.
 */

/* algorithm: 0 = all pairs up to 2000 nodes per component and sparse
 * stress beyond, 1 = all pairs (stress majorization), 2 = sparse stress
 * with pivots. Defaults: 50 pivots for the starting layout of all-pairs
 * stress and 100 for sparse stress; 500 (all pairs) or 100 (sparse)
 * iterations; tol 1e-4. seed drives the perturbation of the starting
 * layout and the sample of sparse-stress pivots. */
int32_t net_layout_stress(const GraphHandle *h, int32_t algorithm,
                          int32_t weighted, int32_t weights_as, int64_t pivots,
                          int64_t max_iter, double tol, uint64_t seed,
                          int32_t threads, double *x, double *y, char *err,
                          size_t err_len);

/* Fruchterman-Reingold. approximation: 0 = exact up to 1000 nodes per
 * component and Barnes-Hut beyond, 1 = exact, 2 = Barnes-Hut with opening
 * angle theta (default 1). Default 500 iterations. seed drives the
 * starting positions. */
int32_t net_layout_fr(const GraphHandle *h, int32_t weighted,
                      int32_t weights_as, int64_t iterations,
                      int32_t approximation, double theta, uint64_t seed,
                      int32_t threads, double *x, double *y, char *err,
                      size_t err_len);

/* kind: 0 = circle, 1 = star, 2 = grid. param: the star's center (a node
 * index, or -1 for the node with the most distinct neighbors) or the
 * grid's width (0 for the smallest square); ignored for the circle. */
int32_t net_layout_simple(const GraphHandle *h, int32_t kind, int64_t param,
                          double *x, double *y, char *err, size_t err_len);

/* Layered tree of each component. roots (may be NULL if n_roots is 0)
 * gives at most one root per component. */
int32_t net_layout_tree(const GraphHandle *h, const int64_t *roots,
                        int64_t n_roots, double *x, double *y, char *err,
                        size_t err_len);

/* Two rows for a two-mode network: types[v] == 0 on top. sweeps: number
 * of barycenter sweeps (0 for 24). */
int32_t net_layout_twomode(const GraphHandle *h, const int32_t *types,
                           int64_t sweeps, double *x, double *y, char *err,
                           size_t err_len);

/* Writes the stress of the layout (x, y) to out: the sum over pairs in
 * the same component of (||x_i - x_j|| - d_ij)^2 / d_ij^2. */
int32_t net_layout_stress_value(const GraphHandle *h, const double *x,
                                const double *y, int32_t weighted,
                                int32_t weights_as, int32_t threads,
                                double *out, char *err, size_t err_len);

/* ---- Key players, group centrality, flow betweenness ----------------
 * Declarations to be merged into netcapi.h.
 *
 * Node sets are 0-based indices without repeats.
 * measure (key players): 0 = fragmentation, 1 = distance-weighted
 *   fragmentation, 2 = m-reach, 3 = distance-weighted reach.
 * measure (group centrality): 0 = degree, 1 = closeness, 2 = betweenness.
 * mode: 1 = out, 2 = in, 3 = all. weights_as: 0 = distance, 1 = strength.
 * method: 0 = greedy swaps with random restarts, 1 = exact enumeration.
 * scale (flow betweenness): 0 = raw, 1 = normalized, 2 = fractional. */

/* Criterion value for the k nodes in set; m is the m-reach distance limit. */
int32_t net_kpp_score(const GraphHandle *h, const int64_t *set, int64_t k,
                      int32_t measure, double m, int32_t mode, int32_t weighted,
                      int32_t weights_as, double *out, char *err,
                      size_t err_len);

/* Finds k key players. Greedy restart r starts from stream r of seed. set
 * receives k indices in increasing order; score one double; evaluations
 * the number of sets evaluated. */
int32_t net_key_players(const GraphHandle *h, int64_t k, int32_t measure,
                        double m, int32_t mode, int32_t weighted,
                        int32_t weights_as, int32_t method, int64_t restarts,
                        uint64_t seed, int32_t threads, int64_t *set,
                        double *score, double *evaluations, char *err,
                        size_t err_len);

/* Group degree, closeness, or betweenness of the k nodes in group. */
int32_t net_group_centrality(const GraphHandle *h, const int64_t *group,
                             int64_t k, int32_t measure, int32_t mode,
                             int32_t weighted, int32_t weights_as,
                             int32_t normalized, int32_t threads, double *out,
                             char *err, size_t err_len);

/* Flow betweenness (Freeman, Borgatti, and White 1991) of every node into
 * out (n doubles). capacity != 0 uses edge weights as capacities. */
int32_t net_flow_betweenness(const GraphHandle *h, int32_t capacity,
                             int32_t scale, int32_t threads, double *out,
                             char *err, size_t err_len);

/*
 * Temporal networks: contacts with start and end times, time-respecting
 * paths, reachability, temporal closeness, slicing into windows,
 * aggregation, and inter-event statistics. To be merged into netcapi.h
 * (after the GraphHandle declarations; needs <stddef.h> and <stdint.h>).
 *
 * Codes:
 *   criterion: 0 = earliest arrival, 1 = latest departure, 2 = fastest,
 *              3 = shortest (fewest contacts).
 *   direction: 1 = walks from the roots, 2 = walks into the roots.
 *   traversal: 0 = the constant traversal_time for every contact,
 *              1 = each contact's duration (end - start).
 *   scope:     0 = nodes, 1 = node pairs; mode (nodes of directed
 *              networks): 1 = out, 2 = in, 3 = all.
 * Window bounds start and end may be NaN for the earliest start time and
 * the latest arrival of any contact; max_wait may be infinite. strict != 0
 * requires waiting times strictly greater than min_wait.
 */

typedef struct TemporalHandle TemporalHandle;

/* Numbered columns of doubles, for results whose length is not known in
 * advance. */
typedef struct TemporalResult TemporalResult;

/* end may be NULL for instantaneous contacts; otherwise end[i] >= start[i].
 * Times must be finite. Release the handle with net_temporal_free. */
int32_t net_temporal_new(int64_t n, const int64_t *src, const int64_t *tgt,
                         const double *start, const double *end, int64_t m,
                         int32_t directed, TemporalHandle **out, char *err,
                         size_t err_len);

void net_temporal_free(TemporalHandle *h);

int64_t net_temporal_node_count(const TemporalHandle *h);

/* out[r * n + j]: walks from roots[r] to j (direction 1) or from j to
 * roots[r] (direction 2). Unreachable: +Inf (-Inf for latest departure).
 * out must hold nroots * n doubles. */
int32_t net_temporal_paths(const TemporalHandle *h, int32_t criterion,
                           int32_t direction, const int64_t *roots,
                           int64_t nroots, double start, double end,
                           double min_wait, double max_wait, int32_t strict,
                           int32_t traversal, double traversal_time,
                           int32_t threads, double *out, char *err,
                           size_t err_len);

/* Number of other nodes each node reaches (direction 1) or is reached by
 * (direction 2). out must hold n doubles. */
int32_t net_temporal_reach(const TemporalHandle *h, int32_t direction,
                           double start, double end, double min_wait,
                           double max_wait, int32_t strict, int32_t traversal,
                           double traversal_time, int32_t threads, double *out,
                           char *err, size_t err_len);

/* Pan-Saramaki temporal closeness (n values) and, unless distances is NULL,
 * average temporal distances (n * n values, row = origin). */
int32_t net_temporal_closeness(const TemporalHandle *h, double start,
                               double end, double min_wait, double max_wait,
                               int32_t strict, int32_t traversal,
                               double traversal_time, int32_t threads,
                               double *closeness, double *distances, char *err,
                               size_t err_len);

/* Contacts active in [from[k], to[k]). Columns: 0 = offsets (k + 1),
 * 1 = 0-based contact indices. */
int32_t net_temporal_windows(const TemporalHandle *h, const double *from,
                             const double *to, int64_t k, TemporalResult **out,
                             char *err, size_t err_len);

/* Contacts active at time t. Column 0: 0-based contact indices. */
int32_t net_temporal_active_at(const TemporalHandle *h, double t,
                               TemporalResult **out, char *err,
                               size_t err_len);

/* Contacts grouped by node pair (within [from, to) when windowed != 0,
 * durations clipped). Columns: 0 = from, 1 = to, 2 = count, 3 = duration,
 * 4 = first start, 5 = last end, 6 = pair of each contact (-1 outside). */
int32_t net_temporal_aggregate(const TemporalHandle *h, int32_t windowed,
                               double from, double to, TemporalResult **out,
                               char *err, size_t err_len);

/* Columns: 0 = group of each inter-event time, 1 = inter-event times,
 * 2 = events per group. */
int32_t net_temporal_inter_event(const TemporalHandle *h, int32_t scope,
                                 int32_t mode, TemporalResult **out, char *err,
                                 size_t err_len);

/* Columns: 0 = events, 1 = burstiness B, 2 = memory M (Goh-Barabasi). */
int32_t net_temporal_burstiness(const TemporalHandle *h, int32_t scope,
                                int32_t mode, TemporalResult **out, char *err,
                                size_t err_len);

/* Length of column col, or -1 if there is no such column. */
int64_t net_temporal_result_length(const TemporalResult *h, int32_t col);

int32_t net_temporal_result_copy(const TemporalResult *h, int32_t col,
                                 double *out, char *err, size_t err_len);

void net_temporal_result_free(TemporalResult *h);

/*
 * Reach, edge betweenness, information centrality, bridging, regular and
 * automorphic equivalence, blockmodels, classical scaling, and geospatial
 * networks. To be merged into netcapi.h.
 *
 * mode: 1 = out, 2 = in, 3 = all. weights_as: 0 = distance, 1 = strength
 * (lengths 1/w). crs: 0 = longitude and latitude in degrees on a sphere of
 * radius `radius` (distances in its units), 1 = planar (Euclidean).
 * Matrices are column-major.
 */

typedef struct EdgeListHandle EdgeListHandle;

/* Number of other nodes within distance k (steps, or a length when
 * weighted != 0; k may be infinite), divided by n - 1 when
 * normalized != 0. */
int32_t net_reach(const GraphHandle *h, double k, int32_t mode,
                  int32_t weighted, int32_t weights_as, int32_t normalized,
                  int32_t threads, double *out, char *err, size_t err_len);

/* Shortest-path betweenness of every edge, in input order; out holds one
 * double per edge. normalized != 0 divides by n(n - 1) (directed) or
 * n(n - 1)/2 (undirected). */
int32_t net_edge_betweenness(const GraphHandle *h, int32_t weighted,
                             int32_t weights_as, int32_t normalized,
                             int32_t threads, double *out, char *err,
                             size_t err_len);

/* Stephenson-Zelen information centrality on the symmetrized network.
 * With weighted != 0, weights are tie strengths (strictly positive).
 * Isolates get 0; if the nodes with ties are not connected, they get NaN. */
int32_t net_information_centrality(const GraphHandle *h, int32_t weighted,
                                   double *out, char *err, size_t err_len);

/* Hwang et al.'s bridging coefficient; NaN for isolated nodes. */
int32_t net_bridging_coefficient(const GraphHandle *h, double *out, char *err,
                                 size_t err_len);

/* Maximal regular equivalence (CATREGE). seed may be NULL, or hold n_seed
 * (= n) labels of a partition to refine. out receives n 0-based class
 * numbers (in order of first node); rounds the number of refinement rounds
 * that split a class. With weighted != 0, total edge weights are tie
 * categories. */
int32_t net_regular_equivalence(const GraphHandle *h, const int64_t *seed,
                                int64_t n_seed, int32_t weighted,
                                int32_t *out, int32_t *rounds, char *err,
                                size_t err_len);

/* Automorphic equivalence classes (orbits of the automorphism group of the
 * binary ties), 0-based in order of first node. Fails rather than run past
 * max_steps refined colourings; max_steps <= 0 uses the default. */
int32_t net_automorphic_equivalence(const GraphHandle *h, int64_t max_steps,
                                    int32_t *out, char *err, size_t err_len);

/* Block densities of a partition with k blocks (labels numbered in order
 * of first appearance). density, sums, and pairs receive k-by-k matrices;
 * summary receives the overall density and the R-squared fit. */
int32_t net_block_model(const GraphHandle *h, const int64_t *membership,
                        int64_t n_membership, int64_t k, int32_t weighted,
                        int32_t loops, double *density, double *sums,
                        double *pairs, double *summary, char *err,
                        size_t err_len);

/* Classical multidimensional scaling of the n-by-n dissimilarities d into
 * at most k dimensions. points receives an n-by-k matrix whose first *dims
 * columns hold coordinates (the rest NaN); eigenvalues all n eigenvalues
 * in decreasing order; gof cmdscale's two goodness-of-fit values. */
int32_t net_classical_mds(const double *d, int64_t n, int64_t k,
                          double *points, int32_t *dims, double *eigenvalues,
                          double *gof, char *err, size_t err_len);

/* Kruskal's stress formula 1 of the n-by-p configuration points against
 * the n-by-n dissimilarities d; writes one double. */
int32_t net_mds_stress(const double *d, int64_t n, const double *points,
                       int64_t p, double *out, char *err, size_t err_len);

/* Length of every edge between its endpoints' coordinates (x, y hold
 * n_points = n values); out holds one double per edge. */
int32_t net_edge_lengths(const GraphHandle *h, const double *x,
                         const double *y, int64_t n_points, int32_t crs,
                         double radius, double *out, char *err,
                         size_t err_len);

/* Proximity network of n located points: with k > 0 the k nearest
 * neighbours of every point (directed edges, nearest first, ties by
 * index); otherwise every pair within distance r (undirected, smaller
 * index first). Release *out with net_edgelist_free. */
int32_t net_proximity(const double *x, const double *y, int64_t n,
                      int32_t crs, double radius, double r, int64_t k,
                      int32_t threads, EdgeListHandle **out, char *err,
                      size_t err_len);

/* Number of edges in an edge list, or -1 for NULL. */
int64_t net_edgelist_size(const EdgeListHandle *h);

/* Copies 0-based endpoints and distances (net_edgelist_size values each). */
int32_t net_edgelist_copy(const EdgeListHandle *h, int64_t *from, int64_t *to,
                          double *distance, char *err, size_t err_len);

void net_edgelist_free(EdgeListHandle *h);

/* Centroid and dispersion of the k groups given by groups (one label per
 * point, k distinct, numbered in order of first appearance). Each output
 * holds k doubles: sizes, centroid x and y, mean distance from the
 * centroid, and standard distance. */
int32_t net_group_geometry(const double *x, const double *y,
                           const int64_t *groups, int64_t n, int32_t crs,
                           double radius, int64_t k, double *sizes,
                           double *cx, double *cy, double *mean_distance,
                           double *standard_distance, char *err,
                           size_t err_len);

/*
 * Cohesive subsets, community partitions, topography (core-periphery,
 * compactness, small worlds), and random graphs. To be merged into
 * netcapi.h.
 */

/* A collection of node sets (bicomponents, n-cliques, k-plexes), read in
 * two steps like GroupsHandle: query the number of sets and members,
 * allocate, then copy. Members of set i are members[offsets[i] ..
 * offsets[i + 1]), 0-based and increasing; sets are in lexicographic
 * order. Release with net_node_sets_free. */
typedef struct NodeSetsHandle NodeSetsHandle;

/* Number of sets, or -1 for a NULL handle. */
int64_t net_node_sets_count(const NodeSetsHandle *h);
/* Total number of members over all sets, or -1 for a NULL handle. */
int64_t net_node_sets_total(const NodeSetsHandle *h);
/* offsets must hold count + 1 values and members total values. */
int32_t net_node_sets_copy(const NodeSetsHandle *h, int64_t *offsets,
                           int64_t *members, char *err, size_t err_len);
void net_node_sets_free(NodeSetsHandle *h);

/* Bicomponents (blocks: maximal subgraphs without a cut point) of the
 * undirected simple graph underlying the network (directions, loops,
 * and edge multiplicity ignored). A bridge is a block of two nodes;
 * isolated nodes belong to no block. */
int32_t net_bicomponents(const GraphHandle *h, NodeSetsHandle **out,
                         char *err, size_t err_len);

/* Maximal n-cliques (clans == 0) or n-clans (clans != 0) of the
 * undirected simple graph with at least min_size members: sets in which
 * every pair is within distance `order` in the whole network (n-clique)
 * or within the set itself (n-clan). threads == 0 uses all cores. */
int32_t net_n_cliques(const GraphHandle *h, int64_t order, int64_t min_size,
                      int32_t clans, int32_t threads, NodeSetsHandle **out,
                      char *err, size_t err_len);

/* Maximal k-plexes (every member adjacent to at least size - k others)
 * with at least min_size members. The search stops with an error after
 * max_steps steps (at least 1), a guard against exponential blow-up. */
int32_t net_k_plexes(const GraphHandle *h, int64_t k, int64_t min_size,
                     double max_steps, int32_t threads, NodeSetsHandle **out,
                     char *err, size_t err_len);

/* Edge connectivity (number of edge-disjoint paths) between every pair of
 * nodes of the undirected simple graph: an n-by-n symmetric matrix with
 * zero diagonal, written to out in column-major order. Lambda sets are
 * the classes of nodes joined by connectivity >= k. */
int32_t net_lambda_connectivity(const GraphHandle *h, double *out, char *err,
                                size_t err_len);

/* Factions: a partition into k non-empty groups (1 <= k <= n) by tabu
 * search from `starts` random partitions. fit: 0 = fewest errors (absent
 * ties within groups plus ties between groups), 1 = largest correlation
 * with the ideal pattern (needs 1 < k < n). A start ends after max_stall
 * steps without improvement; tenure == 0 chooses the default tabu
 * tenure. Writes 0-based labels (numbered by first appearance), the
 * number of errors, and the correlation (NaN where undefined) of the
 * partition. The same seed gives the same result on any thread count. */
int32_t net_factions(const GraphHandle *h, int64_t k, int32_t fit,
                     int64_t starts, int64_t max_stall, int64_t tenure,
                     uint64_t seed, int32_t threads, int64_t *membership,
                     double *errors, double *correlation, char *err,
                     size_t err_len);

/* Girvan-Newman divisive clustering of an undirected network (edges
 * removed in order of highest betweenness). Writes 0-based labels of the
 * partition of highest modularity to membership[0..n) and its modularity
 * to *q; the removed edges (smaller node first) to removed_from /
 * removed_to, the modularity before any removal and after each removal
 * to path, and the number of removals to *removed. removed_from and
 * removed_to must hold m values and path m + 1, where m is the number of
 * edges. */
int32_t net_girvan_newman(const GraphHandle *h, int32_t threads,
                          int64_t *membership, double *q, int64_t *removed_from,
                          int64_t *removed_to, double *path, int64_t *removed,
                          char *err, size_t err_len);

/* Newman's leading eigenvector communities of an undirected network (at
 * most 3000 nodes): 0-based labels to membership[0..n), modularity to *q,
 * number of divisions to *splits. Weights (weighted != 0) must be
 * non-negative. */
int32_t net_leading_eigen(const GraphHandle *h, int32_t weighted,
                          int64_t *membership, double *q, int64_t *splits,
                          char *err, size_t err_len);

/* Clauset-Newman-Moore greedy modularity communities of an undirected
 * network: 0-based labels to membership[0..n), modularity to *q.
 * resolution is the Reichardt-Bornholdt parameter (1 for plain
 * modularity). */
int32_t net_fast_greedy(const GraphHandle *h, double resolution,
                        int32_t weighted, int64_t *membership, double *q,
                        char *err, size_t err_len);

/* Discrete core-periphery fit (Borgatti and Everett 1999) by tabu search:
 * core[i] is 1 for core nodes and 0 otherwise, and the correlation
 * between the data and the pattern is written to *correlation. links:
 * 0 = core-periphery pairs left out of the fit, 1 = expected tied,
 * 2 = expected untied. weighted != 0 uses total weights between nodes as
 * values. Other arguments as in net_factions. */
int32_t net_core_periphery(const GraphHandle *h, int32_t links,
                           int32_t weighted, int64_t starts, int64_t max_stall,
                           int64_t tenure, uint64_t seed, int32_t threads,
                           int32_t *core, double *correlation, char *err,
                           size_t err_len);

/* Continuous core-periphery fit: coreness c (non-negative, unit Euclidean
 * norm) such that the pattern c_i c_j fits the data. fit: 0 = MINRES
 * (least squares, sum over i != j of (a_ij - c_i c_j)^2, UCINET's method;
 * deterministic coordinate descent with convergence tolerance tol on the
 * relative change of coreness, at most max_iter sweeps), 1 = maximum
 * correlation between the data and the pattern (Borgatti and Everett
 * 1999), by BFGS from `starts` points (the degree start plus seeded random
 * ones), at most max_iter iterations each. Writes coreness to
 * coreness[0..n), the correlation of the final pattern to *correlation,
 * the residual sum of squares to *residual (NaN for the correlation fit),
 * and the iterations (sweeps) made to *iterations. */
int32_t net_core_periphery_continuous(const GraphHandle *h, int32_t fit,
                                      int32_t weighted, int64_t max_iter,
                                      double tol, int64_t starts, uint64_t seed,
                                      int32_t threads, double *coreness,
                                      double *correlation, double *residual,
                                      int64_t *iterations, char *err,
                                      size_t err_len);

/* Compactness (Borgatti 2006): the mean of 1 / d over ordered pairs of
 * distinct nodes, unreachable pairs contributing 0; 1 minus the
 * distance-weighted fragmentation. mode: 1 = out, 2 = in, 3 = all.
 * weights_as: 0 = distance, 1 = strength. NaN below two nodes. */
int32_t net_compactness(const GraphHandle *h, int32_t mode, int32_t weighted,
                        int32_t weights_as, int32_t threads, double *out,
                        char *err, size_t err_len);

/* Small-world quotient (Watts and Strogatz 1998): clustering and mean
 * distance of the undirected simple graph against `reps` G(n, m) graphs
 * of the same size. Writes five values to out: clustering, mean distance,
 * mean random clustering, mean random distance, and the quotient
 * (C / C_random) / (L / L_random). local != 0 uses mean local clustering
 * instead of transitivity. */
int32_t net_small_world(const GraphHandle *h, int32_t local, int64_t reps,
                        uint64_t seed, int32_t threads, double *out, char *err,
                        size_t err_len);

/* ---- Random graphs -------------------------------------------------
 * Every generator returns `reps` edge lists on nodes 0..n in one handle;
 * replicate r uses random stream (seed, r), so a seed reproduces every
 * replicate on any platform and thread count. Edges are in increasing
 * (from, to) order with from < to in undirected graphs. */
typedef struct EdgeListsHandle EdgeListsHandle;

/* Number of lists, or -1 for a NULL handle. */
int64_t net_edge_lists_count(const EdgeListsHandle *h);
/* Total number of edges over all lists, or -1 for a NULL handle. */
int64_t net_edge_lists_total(const EdgeListsHandle *h);
/* The edges of list i are positions offsets[i] .. offsets[i + 1]) of from
 * and to. offsets must hold count + 1 values, from and to total values. */
int32_t net_edge_lists_copy(const EdgeListsHandle *h, int64_t *offsets,
                            int64_t *from, int64_t *to, char *err,
                            size_t err_len);
void net_edge_lists_free(EdgeListsHandle *h);

/* G(n, p): each possible pair tied independently with probability p.
 * loops != 0 makes self-loops possible. */
int32_t net_random_gnp(int64_t n, double p, int32_t directed, int32_t loops,
                       int64_t reps, uint64_t seed, int32_t threads,
                       EdgeListsHandle **out, char *err, size_t err_len);

/* G(n, m): m distinct pairs drawn uniformly. */
int32_t net_random_gnm(int64_t n, int64_t m, int32_t directed, int32_t loops,
                       int64_t reps, uint64_t seed, int32_t threads,
                       EdgeListsHandle **out, char *err, size_t err_len);

/* Configuration model: degrees holds n degrees of an undirected graph
 * (in_degrees NULL) or n out-degrees with in_degrees holding n
 * in-degrees (directed). simple == 0 returns the multigraph (loops and
 * parallel edges kept, a loop counting 2 toward an undirected degree);
 * simple != 0 repeats each matching up to max_tries times until it is
 * simple, and fails if none is. */
int32_t net_random_configuration(int64_t n, const int64_t *degrees,
                                 const int64_t *in_degrees, int32_t simple,
                                 int64_t max_tries, int64_t reps, uint64_t seed,
                                 int32_t threads, EdgeListsHandle **out,
                                 char *err, size_t err_len);

/* U|MAN: a directed graph drawn uniformly from those with `mutual`
 * mutual and `asym` asymmetric dyads (the rest null). */
int32_t net_random_uman(int64_t n, int64_t mutual, int64_t asym, int64_t reps,
                        uint64_t seed, int32_t threads, EdgeListsHandle **out,
                        char *err, size_t err_len);

/* Watts-Strogatz: a ring lattice joining each node to its k nearest
 * neighbors (k even), each edge rewired with probability p. */
int32_t net_random_small_world(int64_t n, int64_t k, double p, int64_t reps,
                               uint64_t seed, int32_t threads,
                               EdgeListsHandle **out, char *err,
                               size_t err_len);

/* Barabasi-Albert: from a star on m + 1 nodes, each new node joins m
 * distinct nodes drawn with probability proportional to degree. */
int32_t net_random_preferential(int64_t n, int64_t m, int64_t reps,
                                uint64_t seed, int32_t threads,
                                EdgeListsHandle **out, char *err,
                                size_t err_len);

/* Stochastic block model: k blocks of sizes[r] nodes (numbered block by
 * block) and tie probabilities probs[r * k + s] (row-major; symmetric for
 * undirected graphs). */
int32_t net_random_sbm(int64_t k, const int64_t *sizes, const double *probs,
                       int32_t directed, int32_t loops, int64_t reps,
                       uint64_t seed, int32_t threads, EdgeListsHandle **out,
                       char *err, size_t err_len);

#ifdef __cplusplus
}
#endif

#endif
