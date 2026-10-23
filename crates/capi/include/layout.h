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
