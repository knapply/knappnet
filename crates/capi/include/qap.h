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
