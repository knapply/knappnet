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
