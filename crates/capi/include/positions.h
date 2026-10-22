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
