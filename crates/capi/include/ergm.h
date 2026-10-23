/* ---- Exponential random graph models ----------------------------------
 * A model is built once from a network, its covariates, and a list of
 * terms, then queried through an ErgmHandle. Simulated networks live in an
 * ErgmSample. Matrices are column-major, as R stores them. All functions
 * return 0 on success and write a message to err otherwise.
 *
 * Term codes (term_code), with the meaning of term_param and term_attr:
 *   0 edges; 1 mutual;
 *   2 nodecov, 3 nodeicov, 4 nodeocov        attr = numeric attribute;
 *   5 nodefactor, 6 nodeifactor, 7 nodeofactor
 *                                             attr = categorical attribute,
 *                                             levels = kept level codes;
 *   8 nodematch                               attr = categorical attribute;
 *   9 nodematch (diff = TRUE)                 attr = categorical, levels;
 *  10 absdiff                                 param = exponent, attr = numeric;
 *  11 edgecov                                 attr = dyadic covariate;
 *  12 kstar, 13 istar, 14 ostar               param = k;
 *  15 triangle, 16 ttriple, 17 ctriple;
 *  18 gwesp, 19 gwdsp, 20 gwnsp               param = decay (fixed);
 *  21 gwdegree, 22 gwidegree, 23 gwodegree    param = decay (fixed);
 *  24 degree, 25 idegree, 26 odegree          param = k;
 *  27 isolates.
 * The levels of term t are term_levels[term_levels_off[t] ..
 * term_levels_off[t] + term_levels_len[t]]. numeric holds k_numeric blocks
 * of n doubles, categorical k_categorical blocks of n level codes (0-based,
 * no missing values), dyadic k_dyadic column-major n-by-n matrices. fixed
 * may be NULL, or n*n row-major flags marking dyads that keep their
 * observed state (structural zeros and ones). Self-loops are an error and
 * parallel edges collapse to one tie. */

typedef struct ErgmHandle ErgmHandle;
typedef struct ErgmSample ErgmSample;

int32_t net_ergm_new(int64_t n, const int64_t *src, const int64_t *tgt,
                     int64_t m, int32_t directed, const int32_t *term_code,
                     const double *term_param, const int32_t *term_attr,
                     const int32_t *term_levels_off,
                     const int32_t *term_levels_len,
                     const int32_t *term_levels, int32_t nterms,
                     const double *numeric, int32_t k_numeric,
                     const int32_t *categorical, int32_t k_categorical,
                     const double *dyadic, int32_t k_dyadic,
                     const int32_t *fixed, ErgmHandle **out, char *err,
                     size_t err_len);

void net_ergm_free(ErgmHandle *h);

/* Number of statistics, number of free dyads, and whether every term is
 * dyad independent (1/0); -1 for a null handle. */
int64_t net_ergm_width(const ErgmHandle *h);
int64_t net_ergm_dyads(const ErgmHandle *h);
int32_t net_ergm_dyad_independent(const ErgmHandle *h);

/* The observed statistics (width doubles). */
int32_t net_ergm_stats(const ErgmHandle *h, double *out, char *err,
                       size_t err_len);

/* Change statistics of every free dyad of the observed network, in
 * canonical order (tail < head undirected; tail != head directed): x is
 * dyads-by-width, y the observed dyad states, tail and head 0-based. */
int32_t net_ergm_changestats(const ErgmHandle *h, double *x, double *y,
                             int64_t *tail, int64_t *head, char *err,
                             size_t err_len);

/* Logistic regression of y (0/1) on the m-by-p design x with optional case
 * weights w (NULL for unit weights); no intercept unless the design has
 * one. max_iter <= 0 and tol <= 0 use the defaults (50, 1e-12). fit
 * receives 5 values: deviance, null deviance (all coefficients 0), total
 * weight, iterations, converged. */
int32_t net_logit(const double *x, const double *y, const double *w,
                  int64_t m, int64_t p, int32_t max_iter, double tol,
                  double *coef, double *cov, double *fit, char *err,
                  size_t err_len);

/* Maximum pseudolikelihood estimate: coef (width), cov (width-by-width,
 * the logistic covariance, valid only for dyad-independent models), and
 * fit (6): deviance, null deviance, dyads, iterations, converged, dyad
 * independent. */
int32_t net_ergm_mple(const ErgmHandle *h, int32_t max_iter, double tol,
                      double *coef, double *cov, double *fit, char *err,
                      size_t err_len);

/* Metropolis-Hastings sample at theta: per chain, burnin proposals then
 * nsim samples interval proposals apart; chain c uses stream c of seed.
 * proposal: 0 = tie/no-tie, 1 = random dyad. start_empty != 0 starts from
 * the empty network instead of the observed one. stats receives a
 * (chains * nsim)-by-width matrix, rate the proposals accepted and made.
 * If out is not NULL and keep_networks != 0, the networks are kept in a
 * sample handle to release with net_ergm_sample_free. */
int32_t net_ergm_simulate(const ErgmHandle *h, const double *theta,
                          int64_t burnin, int64_t interval, int64_t nsim,
                          int32_t chains, int32_t proposal, uint64_t seed,
                          int32_t threads, int32_t keep_networks,
                          int32_t start_empty, double *stats, double *rate,
                          ErgmSample **out, char *err, size_t err_len);

void net_ergm_sample_free(ErgmSample *s);

/* Networks kept in a sample, ties of network k, and its tie list as
 * 0-based tail/head arrays. */
int64_t net_ergm_sample_count(const ErgmSample *s);
int64_t net_ergm_sample_ties(const ErgmSample *s, int64_t k);
int32_t net_ergm_sample_edges(const ErgmSample *s, int64_t k, int64_t *tail,
                              int64_t *head, char *err, size_t err_len);

/* Goodness-of-fit statistics. which: 0 = degree, 1 = in-degree,
 * 2 = out-degree (n values: nodes of degree 0..n-1), 3 = edgewise and
 * 4 = dyadwise shared partners (n-1 values: 0..n-2 partners), 5 = geodesic
 * distance (n values: pairs at distance 1..n-1, then unreachable pairs),
 * 6 = triad census (16 values, MAN order; directed only). The sample form
 * writes a count-by-len matrix, the observed form len values. */
int32_t net_ergm_gof_sample(const ErgmSample *s, int32_t which, double *out,
                            char *err, size_t err_len);
int32_t net_ergm_gof_observed(const ErgmHandle *h, int32_t which,
                              double *out, char *err, size_t err_len);

/* Monte Carlo maximum likelihood from theta0 with the sampling settings of
 * net_ergm_simulate per iteration, at most max_iter iterations, convex
 * hull margin (0.05), and chi-squared critical value for the convergence
 * test of the estimating equation. theta: estimate; cov: inverse Fisher
 * information; mcmc_cov: Monte Carlo term of the covariance; info (5):
 * iterations, converged, last step length, convergence statistic, rows of
 * stats; stats: the final sample, (chains * nsim)-by-width. */
int32_t net_ergm_mcmle(const ErgmHandle *h, const double *theta0,
                       int64_t burnin, int64_t interval, int64_t nsim,
                       int32_t chains, int32_t proposal, uint64_t seed,
                       int32_t threads, int32_t max_iter, double margin,
                       double critical, double *theta, double *cov,
                       double *mcmc_cov, double *info, double *stats,
                       char *err, size_t err_len);

/* Exact log-likelihood at theta by enumerating every network on the free
 * dyads (at most 22), with the exact mean (width) and covariance
 * (width-by-width) of the statistics. */
int32_t net_ergm_exact(const ErgmHandle *h, const double *theta,
                       double *loglik, double *mean, double *cov, char *err,
                       size_t err_len);

/* Exact maximum likelihood estimate by Newton's method from start; cov is
 * the inverse Fisher information and info (2) the log-likelihood at the
 * estimate and the iterations used. Fails when the estimate does not
 * exist. */
int32_t net_ergm_exact_mle(const ErgmHandle *h, const double *start,
                           int32_t max_iter, double *theta, double *cov,
                           double *info, char *err, size_t err_len);

/* Path-sampling estimate of l(theta) - l(0) with steps midpoint samples,
 * each drawn with the settings of net_ergm_simulate; out (2): estimate and
 * Monte Carlo standard error. l(0) = -D log 2 for D free dyads. */
int32_t net_ergm_loglik(const ErgmHandle *h, const double *theta,
                        int32_t steps, int64_t burnin, int64_t interval,
                        int64_t nsim, int32_t chains, int32_t proposal,
                        uint64_t seed, int32_t threads, double *out,
                        char *err, size_t err_len);

/* Pooled pseudolikelihood of a temporal ERGM: logistic regression of y on
 * the stacked m-by-p design x, period[r] the 0-based transition of row r
 * among periods, and reps bootstrap replicates resampling transitions with
 * replacement (replicate r uses stream r of seed). coef, cov, fit as in
 * net_logit; boot: reps-by-p replicate estimates (NaN where a replicate
 * did not converge). */
int32_t net_tergm_mple(const double *x, const double *y,
                       const int32_t *period, int64_t m, int64_t p,
                       int32_t periods, int64_t reps, uint64_t seed,
                       int32_t threads, int32_t max_iter, double tol,
                       double *coef, double *cov, double *fit, double *boot,
                       char *err, size_t err_len);
