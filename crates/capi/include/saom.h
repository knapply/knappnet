/*
 * Stochastic actor-oriented models for network panel data: panels of
 * directed binary networks with covariates, effect statistics and change
 * statistics, ministep simulation, and method-of-moments estimation by
 * the Robbins-Monro algorithm. To be merged into netcapi.h (after the
 * GraphHandle declarations; needs <stddef.h> and <stdint.h>).
 *
 * Effects are passed as three arrays of length k: codes (0 density,
 * 1 recip, 2 transTrip, 3 transTies, 4 transRecTrip, 5 cycle3, 6 gwespFF,
 * 7 gwespFB, 8 gwespBF, 9 gwespBB, 10 inPop, 11 inPopSqrt, 12 outPop,
 * 13 outPopSqrt, 14 inAct, 15 inActSqrt, 16 outAct, 17 outActSqrt,
 * 18 egoX, 19 altX, 20 simX, 21 sameX, 22 egoXaltX, 23 X), covariates
 * (0-based index of the actor covariate for 18-22 or of the dyadic
 * covariate for 23; -1 otherwise), and parameters (the internal effect
 * parameter; the GWESP weight is parameter / 100, RSiena's default is 69).
 *
 * Parameter vectors hold the M - 1 rate parameters first when estimation
 * is unconditional (conditional == 0), then the k effect parameters;
 * conditional estimation has only the k effect parameters.
 *
 * Matrices passed in are column-major as R stores them; matrices read out
 * of a fit are row-major.
 */

typedef struct SaomHandle SaomHandle;

/* An estimation result, read item by item (see the SAOM_FIT_* codes). */
typedef struct SaomFit SaomFit;

/* Number of doubles in the options vector of net_saom_estimate. */
#define SAOM_OPTIONS_LEN 10

/* Items of a fit for net_saom_fit_length and net_saom_fit_copy. p is the
 * number of parameters, M the number of waves, n3 the phase 3 iterations,
 * and A = 3n + 16 the number of auxiliary goodness-of-fit statistics
 * (in-degree counts 0..n-1, out-degree counts 0..n-1, triad census in MAN
 * order, pairs at geodesic distance 1..n-1 and unreachable pairs). */
#define SAOM_FIT_THETA 0                /* p estimates */
#define SAOM_FIT_SE 1                   /* p standard errors */
#define SAOM_FIT_COVARIANCE 2           /* p x p covariance of the estimates */
#define SAOM_FIT_DERIVATIVE 3           /* p x p derivative matrix */
#define SAOM_FIT_STATISTIC_COVARIANCE 4 /* p x p covariance of the statistics */
#define SAOM_FIT_TCONV 5                /* p convergence t-ratios */
#define SAOM_FIT_TCONV_MAX 6            /* 1: overall maximum convergence ratio */
#define SAOM_FIT_TARGETS 7              /* p observed statistics */
#define SAOM_FIT_MEANS 8                /* p mean simulated statistics */
#define SAOM_FIT_RATES 9                /* M - 1 rate parameters */
#define SAOM_FIT_RATE_SE 10             /* M - 1 standard errors of the rates */
#define SAOM_FIT_PHASE2_THETA 11        /* iterations x (1 + p): subphase, theta */
#define SAOM_FIT_PHASE2_DEVIATIONS 12   /* iterations x p deviations S - s */
#define SAOM_FIT_PHASE3_DEVIATIONS 13   /* n3 x p deviations S - s */
#define SAOM_FIT_PHASE3_TIMES 14        /* n3 x (M - 1) simulated period lengths */
#define SAOM_FIT_AUX_SIMULATED 15       /* n3 x A (empty unless keep_aux) */
#define SAOM_FIT_AUX_OBSERVED 16        /* A, summed over waves 2..M */
#define SAOM_FIT_META 17                /* 8: conditional, M - 1, k, n1, n2, n3, p, A */

/* Builds a panel from `waves` (>= 2) adjacency matrices, each n x n
 * column-major, one after another; any non-zero entry is a tie and the
 * diagonal is ignored. Missing (NaN) entries are refused. Release with
 * net_saom_free. */
int32_t net_saom_new(int64_t n, int64_t waves, const double *adjacency,
                     SaomHandle **out, char *err, size_t err_len);

void net_saom_free(SaomHandle *h);

/* Adds an actor covariate: len is n (constant) or n * waves (changing;
 * wave-major, i.e. an n x waves column-major matrix). Covariates are
 * numbered 0, 1, ... in the order added, separately for actor and dyadic
 * covariates. */
int32_t net_saom_add_actor_covariate(SaomHandle *h, const double *values,
                                     int64_t len, char *err, size_t err_len);

/* Adds a constant dyadic covariate (n x n column-major; diagonal ignored). */
int32_t net_saom_add_dyadic_covariate(SaomHandle *h, const double *values,
                                      char *err, size_t err_len);

/* Number of periods M - 1, or -1 for a null handle. */
int64_t net_saom_period_count(const SaomHandle *h);

/* Tie variables that changed in each period (M - 1 doubles). */
int32_t net_saom_distances(const SaomHandle *h, double *out, char *err,
                           size_t err_len);

/* Target statistics: for each period m, the k effect statistics of wave
 * m + 1 with the covariates of wave m; out has k * (M - 1) doubles,
 * period-major (out[m * k + kk]). */
int32_t net_saom_targets(const SaomHandle *h, int64_t k, const int32_t *codes,
                         const int32_t *covariates, const double *parameters,
                         double *out, char *err, size_t err_len);

/* Change statistics of every effect for actor `actor` on wave `wave`
 * (both 0-based) with the covariates of period `period`: n * k doubles,
 * alter-major (out[j * k + kk]); 0 for j == actor. */
int32_t net_saom_contributions(const SaomHandle *h, int64_t k,
                               const int32_t *codes, const int32_t *covariates,
                               const double *parameters, int64_t wave,
                               int64_t period, int64_t actor, double *out,
                               char *err, size_t err_len);

/* Standard initial parameter values: p doubles (k, or k + M - 1 when
 * conditional == 0). */
int32_t net_saom_initial(const SaomHandle *h, int64_t k, const int32_t *codes,
                         const int32_t *covariates, const double *parameters,
                         int32_t conditional, double *out, char *err,
                         size_t err_len);

/* Simulates every period once from the observed network at its start, with
 * effect parameters theta (k) and rates (M - 1; ignored and may be NULL
 * when conditional != 0), using stream `index` of `seed`. networks receives
 * the M - 1 end networks as n x n column-major 0/1 matrices, statistics the
 * k effect statistics per period (period-major), and times the simulated
 * period lengths (stopping times under rate 1 when conditional, else 1).
 * Each output pointer may be NULL. */
int32_t net_saom_simulate(const SaomHandle *h, int64_t k, const int32_t *codes,
                          const int32_t *covariates, const double *parameters,
                          const double *theta, const double *rates,
                          int32_t conditional, uint64_t seed, uint64_t index,
                          double *networks, double *statistics, double *times,
                          char *err, size_t err_len);

/* Method-of-moments estimation. options holds SAOM_OPTIONS_LEN doubles:
 * conditional (0/1), n1 (0 for 7 + 3p), nsub, n3, firstg, reduceg,
 * diagonalize, truncation, theta_bound, keep_aux (0/1). initial may be NULL
 * for the standard initial values, or hold p values with NaN for "standard".
 * threads == 0 uses all cores. Release the fit with net_saom_fit_free. */
int32_t net_saom_estimate(const SaomHandle *h, int64_t k, const int32_t *codes,
                          const int32_t *covariates, const double *parameters,
                          const double *options, const double *initial,
                          uint64_t seed, int32_t threads, SaomFit **out,
                          char *err, size_t err_len);

/* Length of fit item `which`, or -1 if it does not exist. */
int64_t net_saom_fit_length(const SaomFit *h, int32_t which);

/* Copies fit item `which` into out (net_saom_fit_length doubles). */
int32_t net_saom_fit_copy(const SaomFit *h, int32_t which, double *out,
                          char *err, size_t err_len);

void net_saom_fit_free(SaomFit *h);
