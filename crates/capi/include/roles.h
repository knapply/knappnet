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
