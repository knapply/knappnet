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
