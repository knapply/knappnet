/*
 * R bindings for exponential random graph models. As in glue.c, R errors
 * are raised only after the engine call has returned. Engine-owned objects
 * (the model and samples of networks) are attached to external pointers
 * with finalizers before the engine creates them, so an R allocation
 * failure can never leak them.
 *
 * Terms arrive as parallel vectors (code, parameter, attribute index, level
 * offsets and lengths into a flat level vector); covariates as flat double
 * and integer vectors with one block per attribute. See ergm.h for the
 * codes. Matrices are returned in R's column-major layout.
 */
#include "glue.h"

static void ergm_finalizer(SEXP ptr) {
  ErgmHandle *h = (ErgmHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_ergm_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static void sample_finalizer(SEXP ptr) {
  ErgmSample *s = (ErgmSample *)R_ExternalPtrAddr(ptr);
  if (s != NULL) {
    net_ergm_sample_free(s);
    R_ClearExternalPtr(ptr);
  }
}

static ErgmHandle *ergm_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid model handle (internal error)");
  }
  return (ErgmHandle *)R_ExternalPtrAddr(ptr);
}

static ErgmSample *sample_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid sample handle (internal error)");
  }
  return (ErgmSample *)R_ExternalPtrAddr(ptr);
}

static const int32_t *as_int32(SEXP x) {
  if (TYPEOF(x) != INTSXP) Rf_error("expected an integer vector (internal error)");
  return (const int32_t *)INTEGER(x);
}

static const double *as_double(SEXP x) {
  if (TYPEOF(x) != REALSXP) Rf_error("expected a double vector (internal error)");
  return REAL(x);
}

static int ergm_width(ErgmHandle *h) {
  int64_t p = net_ergm_width(h);
  if (p < 0) Rf_error("invalid model handle (internal error)");
  return (int)p;
}

SEXP C_ergm_new(SEXP n, SEXP src, SEXP tgt, SEXP directed, SEXP term_code, SEXP term_param,
                SEXP term_attr, SEXP term_levels_off, SEXP term_levels_len, SEXP term_levels,
                SEXP numeric, SEXP k_numeric, SEXP categorical, SEXP k_categorical, SEXP dyadic,
                SEXP k_dyadic, SEXP fixed) {
  if (TYPEOF(src) != INTSXP || TYPEOF(tgt) != INTSXP || XLENGTH(src) != XLENGTH(tgt)) {
    Rf_error("edge indices must be integer vectors of equal length (internal error)");
  }
  R_xlen_t m = XLENGTH(src);
  const int64_t *s = as_int64(src);
  const int64_t *t = as_int64(tgt);
  int nterms = Rf_length(term_code);
  SEXP ptr = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, ergm_finalizer, TRUE);
  ErgmHandle *h = NULL;
  char err[ERR_LEN] = {0};
  int status =
      net_ergm_new((int64_t)Rf_asInteger(n), s, t, (int64_t)m, Rf_asLogical(directed),
                   as_int32(term_code), as_double(term_param), as_int32(term_attr),
                   as_int32(term_levels_off), as_int32(term_levels_len), as_int32(term_levels),
                   nterms, as_double(numeric), Rf_asInteger(k_numeric), as_int32(categorical),
                   Rf_asInteger(k_categorical), as_double(dyadic), Rf_asInteger(k_dyadic),
                   Rf_isNull(fixed) ? NULL : as_int32(fixed), &h, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(ptr, h);
  UNPROTECT(1);
  return ptr;
}

SEXP C_ergm_alive(SEXP ptr) {
  return Rf_ScalarLogical(TYPEOF(ptr) == EXTPTRSXP && R_ExternalPtrAddr(ptr) != NULL);
}

/* c(width, free dyads, dyad independent). */
SEXP C_ergm_info(SEXP ptr) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 3));
  REAL(out)[0] = (double)net_ergm_width(h);
  REAL(out)[1] = (double)net_ergm_dyads(h);
  REAL(out)[2] = (double)net_ergm_dyad_independent(h);
  UNPROTECT(1);
  return out;
}

SEXP C_ergm_stats(SEXP ptr) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, p));
  char err[ERR_LEN] = {0};
  int status = net_ergm_stats(h, REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* list(x = dyads-by-width matrix, y, tail, head) with 1-based endpoints. */
SEXP C_ergm_changestats(SEXP ptr) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  int64_t d = net_ergm_dyads(h);
  if (d < 0) Rf_error("invalid model handle (internal error)");
  SEXP x = PROTECT(Rf_allocMatrix(REALSXP, (int)d, p));
  SEXP y = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)d));
  SEXP tail = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)d));
  SEXP head = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)d));
  int64_t *t = (int64_t *)R_alloc(d > 0 ? (size_t)d : 1, sizeof(int64_t));
  int64_t *hd = (int64_t *)R_alloc(d > 0 ? (size_t)d : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  int status = net_ergm_changestats(h, REAL(x), REAL(y), t, hd, err, ERR_LEN);
  check_status(status, err, 4);
  for (int64_t i = 0; i < d; i++) {
    INTEGER(tail)[i] = (int)t[i] + 1;
    INTEGER(head)[i] = (int)hd[i] + 1;
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, x);
  SET_VECTOR_ELT(out, 1, y);
  SET_VECTOR_ELT(out, 2, tail);
  SET_VECTOR_ELT(out, 3, head);
  UNPROTECT(5);
  return out;
}

/* list(coef, cov, fit) where fit has 6 values (see ergm.h). */
SEXP C_ergm_mple(SEXP ptr, SEXP max_iter, SEXP tol) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  SEXP coef = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP cov = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP fit = PROTECT(Rf_allocVector(REALSXP, 6));
  char err[ERR_LEN] = {0};
  int status = net_ergm_mple(h, Rf_asInteger(max_iter), Rf_asReal(tol), REAL(coef), REAL(cov),
                             REAL(fit), err, ERR_LEN);
  check_status(status, err, 3);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, coef);
  SET_VECTOR_ELT(out, 1, cov);
  SET_VECTOR_ELT(out, 2, fit);
  UNPROTECT(4);
  return out;
}

/* list(stats = (chains * nsim)-by-width matrix, rate, sample handle or NULL). */
SEXP C_ergm_simulate(SEXP ptr, SEXP theta, SEXP burnin, SEXP interval, SEXP nsim, SEXP chains,
                     SEXP proposal, SEXP seed, SEXP threads, SEXP keep, SEXP start_empty) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  if (TYPEOF(theta) != REALSXP || XLENGTH(theta) != p) {
    Rf_error("coefficients must be a double vector with one value per statistic (internal error)");
  }
  double ns = Rf_asReal(nsim);
  int nchains = Rf_asInteger(chains);
  int keep_networks = Rf_asLogical(keep) == 1;
  R_xlen_t rows = (R_xlen_t)(ns * nchains);
  SEXP stats = PROTECT(Rf_allocMatrix(REALSXP, (int)rows, p));
  SEXP rate = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP sptr = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(sptr, sample_finalizer, TRUE);
  ErgmSample *s = NULL;
  char err[ERR_LEN] = {0};
  int status = net_ergm_simulate(
      h, REAL(theta), (int64_t)Rf_asReal(burnin), (int64_t)Rf_asReal(interval), (int64_t)ns,
      nchains, Rf_asInteger(proposal), as_seed(seed), Rf_asInteger(threads), keep_networks,
      Rf_asLogical(start_empty), REAL(stats), REAL(rate), keep_networks ? &s : NULL, err, ERR_LEN);
  check_status(status, err, 3);
  if (keep_networks) R_SetExternalPtrAddr(sptr, s);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, stats);
  SET_VECTOR_ELT(out, 1, rate);
  SET_VECTOR_ELT(out, 2, keep_networks ? sptr : R_NilValue);
  UNPROTECT(4);
  return out;
}

/* A list with one two-column integer matrix (1-based tail, head) per
 * sampled network. */
SEXP C_ergm_sample_networks(SEXP ptr) {
  ErgmSample *s = sample_from_ptr(ptr);
  int64_t count = net_ergm_sample_count(s);
  if (count < 0) Rf_error("invalid sample handle (internal error)");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)count));
  for (int64_t k = 0; k < count; k++) {
    int64_t m = net_ergm_sample_ties(s, k);
    if (m < 0) {
      UNPROTECT(1);
      Rf_error("sampled network %d is unavailable (internal error)", (int)k + 1);
    }
    SEXP mat = PROTECT(Rf_allocMatrix(INTSXP, (int)m, 2));
    int64_t *t = (int64_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int64_t));
    int64_t *hd = (int64_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int64_t));
    char err[ERR_LEN] = {0};
    int status = net_ergm_sample_edges(s, k, t, hd, err, ERR_LEN);
    check_status(status, err, 2);
    for (int64_t i = 0; i < m; i++) {
      INTEGER(mat)[i] = (int)t[i] + 1;
      INTEGER(mat)[m + i] = (int)hd[i] + 1;
    }
    SET_VECTOR_ELT(out, (R_xlen_t)k, mat);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

/* count-by-len matrix of one goodness-of-fit statistic over a sample; len
 * is supplied by R (n, n - 1, or 16). */
SEXP C_ergm_gof_sample(SEXP ptr, SEXP which, SEXP len) {
  ErgmSample *s = sample_from_ptr(ptr);
  int64_t count = net_ergm_sample_count(s);
  if (count < 0) Rf_error("invalid sample handle (internal error)");
  int l = Rf_asInteger(len);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)count, l));
  char err[ERR_LEN] = {0};
  int status = net_ergm_gof_sample(s, Rf_asInteger(which), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_ergm_gof_observed(SEXP ptr, SEXP which, SEXP len) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, Rf_asInteger(len)));
  char err[ERR_LEN] = {0};
  int status = net_ergm_gof_observed(h, Rf_asInteger(which), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* list(theta, cov, mcmc_cov, info, stats). */
SEXP C_ergm_mcmle(SEXP ptr, SEXP theta0, SEXP burnin, SEXP interval, SEXP nsim, SEXP chains,
                  SEXP proposal, SEXP seed, SEXP threads, SEXP max_iter, SEXP margin,
                  SEXP critical) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  if (TYPEOF(theta0) != REALSXP || XLENGTH(theta0) != p) {
    Rf_error(
        "starting values must be a double vector with one value per statistic (internal error)");
  }
  double ns = Rf_asReal(nsim);
  int nchains = Rf_asInteger(chains);
  R_xlen_t rows = (R_xlen_t)(ns * nchains);
  SEXP theta = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP cov = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP mc = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP info = PROTECT(Rf_allocVector(REALSXP, 5));
  SEXP stats = PROTECT(Rf_allocMatrix(REALSXP, (int)rows, p));
  char err[ERR_LEN] = {0};
  int status = net_ergm_mcmle(
      h, REAL(theta0), (int64_t)Rf_asReal(burnin), (int64_t)Rf_asReal(interval), (int64_t)ns,
      nchains, Rf_asInteger(proposal), as_seed(seed), Rf_asInteger(threads), Rf_asInteger(max_iter),
      Rf_asReal(margin), Rf_asReal(critical), REAL(theta), REAL(cov), REAL(mc), REAL(info),
      REAL(stats), err, ERR_LEN);
  check_status(status, err, 5);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, theta);
  SET_VECTOR_ELT(out, 1, cov);
  SET_VECTOR_ELT(out, 2, mc);
  SET_VECTOR_ELT(out, 3, info);
  SET_VECTOR_ELT(out, 4, stats);
  UNPROTECT(6);
  return out;
}

/* list(loglik, mean, cov) by exact enumeration. */
SEXP C_ergm_exact(SEXP ptr, SEXP theta) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  if (TYPEOF(theta) != REALSXP || XLENGTH(theta) != p) {
    Rf_error("coefficients must be a double vector with one value per statistic (internal error)");
  }
  SEXP ll = PROTECT(Rf_allocVector(REALSXP, 1));
  SEXP mean = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP cov = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  char err[ERR_LEN] = {0};
  int status = net_ergm_exact(h, REAL(theta), REAL(ll), REAL(mean), REAL(cov), err, ERR_LEN);
  check_status(status, err, 3);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, ll);
  SET_VECTOR_ELT(out, 1, mean);
  SET_VECTOR_ELT(out, 2, cov);
  UNPROTECT(4);
  return out;
}

/* list(theta, cov, info) with info = c(loglik, iterations). */
SEXP C_ergm_exact_mle(SEXP ptr, SEXP start, SEXP max_iter) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  if (TYPEOF(start) != REALSXP || XLENGTH(start) != p) {
    Rf_error(
        "starting values must be a double vector with one value per statistic (internal error)");
  }
  SEXP theta = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP cov = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP info = PROTECT(Rf_allocVector(REALSXP, 2));
  char err[ERR_LEN] = {0};
  int status = net_ergm_exact_mle(h, REAL(start), Rf_asInteger(max_iter), REAL(theta), REAL(cov),
                                  REAL(info), err, ERR_LEN);
  check_status(status, err, 3);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, theta);
  SET_VECTOR_ELT(out, 1, cov);
  SET_VECTOR_ELT(out, 2, info);
  UNPROTECT(4);
  return out;
}

/* c(log-likelihood ratio against theta = 0, Monte Carlo standard error). */
SEXP C_ergm_loglik(SEXP ptr, SEXP theta, SEXP steps, SEXP burnin, SEXP interval, SEXP nsim,
                   SEXP chains, SEXP proposal, SEXP seed, SEXP threads) {
  ErgmHandle *h = ergm_from_ptr(ptr);
  int p = ergm_width(h);
  if (TYPEOF(theta) != REALSXP || XLENGTH(theta) != p) {
    Rf_error("coefficients must be a double vector with one value per statistic (internal error)");
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  char err[ERR_LEN] = {0};
  int status = net_ergm_loglik(h, REAL(theta), Rf_asInteger(steps), (int64_t)Rf_asReal(burnin),
                               (int64_t)Rf_asReal(interval), (int64_t)Rf_asReal(nsim),
                               Rf_asInteger(chains), Rf_asInteger(proposal), as_seed(seed),
                               Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* Pooled pseudolikelihood with bootstrap: list(coef, cov, fit, boot). x is
 * an m-by-p double matrix, y a double vector, period an integer vector of
 * 0-based transition labels. */
SEXP C_tergm_mple(SEXP x, SEXP y, SEXP period, SEXP periods, SEXP reps, SEXP seed, SEXP threads,
                  SEXP max_iter, SEXP tol) {
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP) {
    Rf_error("the design and response must be double vectors (internal error)");
  }
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (Rf_isNull(dim) || Rf_length(dim) != 2)
    Rf_error("the design must be a matrix (internal error)");
  int m = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (XLENGTH(y) != m || XLENGTH(period) != m) {
    Rf_error("the response and period labels must have one value per row (internal error)");
  }
  double nreps = Rf_asReal(reps);
  SEXP coef = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP cov = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP fit = PROTECT(Rf_allocVector(REALSXP, 5));
  SEXP boot = PROTECT(Rf_allocMatrix(REALSXP, (int)nreps, p));
  char err[ERR_LEN] = {0};
  int status = net_tergm_mple(REAL(x), REAL(y), as_int32(period), (int64_t)m, (int64_t)p,
                              Rf_asInteger(periods), (int64_t)nreps, as_seed(seed),
                              Rf_asInteger(threads), Rf_asInteger(max_iter), Rf_asReal(tol),
                              REAL(coef), REAL(cov), REAL(fit), REAL(boot), err, ERR_LEN);
  check_status(status, err, 4);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, coef);
  SET_VECTOR_ELT(out, 1, cov);
  SET_VECTOR_ELT(out, 2, fit);
  SET_VECTOR_ELT(out, 3, boot);
  UNPROTECT(5);
  return out;
}
