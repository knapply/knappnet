/*
 * R bindings for stochastic actor-oriented models. As in glue.c, R errors
 * are raised only after the engine call has returned. The panel handle and
 * the estimation result are attached to external pointers with finalizers
 * before the engine creates them, so an R allocation failure can never
 * leak them; the estimation result is copied into R vectors and released
 * at once.
 *
 * Effects arrive as three parallel vectors: integer codes, integer
 * covariate indices (0-based, -1 for none), and double internal
 * parameters. Seeds arrive as doubles holding whole numbers in [0, 2^53),
 * checked in R.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static void saom_finalizer(SEXP ptr) {
  SaomHandle *h = (SaomHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_saom_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static void fit_finalizer(SEXP ptr) {
  SaomFit *f = (SaomFit *)R_ExternalPtrAddr(ptr);
  if (f != NULL) {
    net_saom_fit_free(f);
    R_ClearExternalPtr(ptr);
  }
}

static SaomHandle *saom_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid panel handle (internal error)");
  }
  return (SaomHandle *)R_ExternalPtrAddr(ptr);
}

static uint64_t as_u64(SEXP x) { return (uint64_t)Rf_asReal(x); }

/* The effect specification: codes and covariates integer, parameters
 * double, all of the same length k. */
typedef struct {
  int64_t k;
  const int32_t *codes;
  const int32_t *covariates;
  const double *parameters;
} Spec;

static Spec spec_of(SEXP codes, SEXP covariates, SEXP parameters) {
  if (TYPEOF(codes) != INTSXP || TYPEOF(covariates) != INTSXP || TYPEOF(parameters) != REALSXP ||
      XLENGTH(covariates) != XLENGTH(codes) || XLENGTH(parameters) != XLENGTH(codes)) {
    Rf_error("effects must be integer codes, integer covariates, and double parameters of one "
             "length (internal error)");
  }
  Spec s;
  s.k = (int64_t)XLENGTH(codes);
  s.codes = (const int32_t *)INTEGER(codes);
  s.covariates = (const int32_t *)INTEGER(covariates);
  s.parameters = REAL(parameters);
  return s;
}

/* adjacency: double vector of waves * n * n (column-major matrices one after
 * another); actor: list of double vectors (n or n * waves); dyadic: list of
 * double vectors (n * n). */
SEXP C_saom_new(SEXP n, SEXP waves, SEXP adjacency, SEXP actor, SEXP dyadic) {
  if (TYPEOF(adjacency) != REALSXP) Rf_error("adjacency must be a double vector (internal error)");
  int nn = Rf_asInteger(n);
  int w = Rf_asInteger(waves);
  if (XLENGTH(adjacency) != (R_xlen_t)nn * nn * w) {
    Rf_error("adjacency has the wrong length (internal error)");
  }
  SEXP ptr = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, saom_finalizer, TRUE);
  SaomHandle *h = NULL;
  char err[ERR_LEN] = {0};
  int status = net_saom_new((int64_t)nn, (int64_t)w, REAL(adjacency), &h, err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  R_SetExternalPtrAddr(ptr, h);
  for (R_xlen_t i = 0; i < Rf_xlength(actor); i++) {
    SEXP v = VECTOR_ELT(actor, i);
    if (TYPEOF(v) != REALSXP) Rf_error("actor covariates must be double vectors (internal error)");
    status = net_saom_add_actor_covariate(h, REAL(v), (int64_t)XLENGTH(v), err, ERR_LEN);
    if (status != 0) {
      UNPROTECT(1);
      Rf_error("%s", err);
    }
  }
  for (R_xlen_t i = 0; i < Rf_xlength(dyadic); i++) {
    SEXP v = VECTOR_ELT(dyadic, i);
    if (TYPEOF(v) != REALSXP || XLENGTH(v) != (R_xlen_t)nn * nn) {
      Rf_error("dyadic covariates must be double vectors of length n * n (internal error)");
    }
    status = net_saom_add_dyadic_covariate(h, REAL(v), err, ERR_LEN);
    if (status != 0) {
      UNPROTECT(1);
      Rf_error("%s", err);
    }
  }
  UNPROTECT(1);
  return ptr;
}

SEXP C_saom_alive(SEXP ptr) {
  return Rf_ScalarLogical(TYPEOF(ptr) == EXTPTRSXP && R_ExternalPtrAddr(ptr) != NULL);
}

SEXP C_saom_distances(SEXP ptr) {
  SaomHandle *h = saom_from_ptr(ptr);
  int64_t periods = net_saom_period_count(h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)periods));
  char err[ERR_LEN] = {0};
  if (net_saom_distances(h, REAL(out), err, ERR_LEN) != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns k * (M - 1) values, period-major. */
SEXP C_saom_targets(SEXP ptr, SEXP codes, SEXP covariates, SEXP parameters) {
  SaomHandle *h = saom_from_ptr(ptr);
  Spec s = spec_of(codes, covariates, parameters);
  int64_t periods = net_saom_period_count(h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)(s.k * periods)));
  char err[ERR_LEN] = {0};
  if (net_saom_targets(h, s.k, s.codes, s.covariates, s.parameters, REAL(out), err, ERR_LEN) !=
      0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns n * k values, alter-major; wave, period, and actor are 0-based. */
SEXP C_saom_contributions(SEXP ptr, SEXP codes, SEXP covariates, SEXP parameters, SEXP n,
                          SEXP wave, SEXP period, SEXP actor) {
  SaomHandle *h = saom_from_ptr(ptr);
  Spec s = spec_of(codes, covariates, parameters);
  R_xlen_t nn = (R_xlen_t)Rf_asInteger(n);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, nn * (R_xlen_t)s.k));
  char err[ERR_LEN] = {0};
  if (net_saom_contributions(h, s.k, s.codes, s.covariates, s.parameters,
                             (int64_t)Rf_asInteger(wave), (int64_t)Rf_asInteger(period),
                             (int64_t)Rf_asInteger(actor), REAL(out), err, ERR_LEN) != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_saom_initial(SEXP ptr, SEXP codes, SEXP covariates, SEXP parameters, SEXP conditional) {
  SaomHandle *h = saom_from_ptr(ptr);
  Spec s = spec_of(codes, covariates, parameters);
  int cond = Rf_asLogical(conditional);
  int64_t p = s.k + (cond ? 0 : net_saom_period_count(h));
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)p));
  char err[ERR_LEN] = {0};
  if (net_saom_initial(h, s.k, s.codes, s.covariates, s.parameters, cond, REAL(out), err,
                       ERR_LEN) != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns list(networks, statistics, times): networks holds the M - 1
 * simulated end networks as column-major n x n matrices one after another. */
SEXP C_saom_simulate(SEXP ptr, SEXP codes, SEXP covariates, SEXP parameters, SEXP n, SEXP theta,
                     SEXP rates, SEXP conditional, SEXP seed, SEXP index) {
  SaomHandle *h = saom_from_ptr(ptr);
  Spec s = spec_of(codes, covariates, parameters);
  if (TYPEOF(theta) != REALSXP || XLENGTH(theta) != (R_xlen_t)s.k) {
    Rf_error("theta must be a double vector with one value per effect (internal error)");
  }
  int cond = Rf_asLogical(conditional);
  int64_t periods = net_saom_period_count(h);
  if (!cond && (TYPEOF(rates) != REALSXP || XLENGTH(rates) != (R_xlen_t)periods)) {
    Rf_error("rates must be a double vector with one value per period (internal error)");
  }
  R_xlen_t nn = (R_xlen_t)Rf_asInteger(n);
  SEXP networks = PROTECT(Rf_allocVector(REALSXP, nn * nn * (R_xlen_t)periods));
  SEXP statistics = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)(s.k * periods)));
  SEXP times = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)periods));
  char err[ERR_LEN] = {0};
  if (net_saom_simulate(h, s.k, s.codes, s.covariates, s.parameters, REAL(theta),
                        cond ? NULL : REAL(rates), cond, as_u64(seed), as_u64(index),
                        REAL(networks), REAL(statistics), REAL(times), err, ERR_LEN) != 0) {
    UNPROTECT(3);
    Rf_error("%s", err);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, networks);
  SET_VECTOR_ELT(out, 1, statistics);
  SET_VECTOR_ELT(out, 2, times);
  UNPROTECT(4);
  return out;
}

#define FIT_ITEMS 18

/* Returns a list of the FIT_ITEMS fit items (double vectors), in the order
 * of the SAOM_FIT_* codes. options: SAOM_OPTIONS_LEN doubles; initial: a
 * double vector of p values or NULL. */
SEXP C_saom_estimate(SEXP ptr, SEXP codes, SEXP covariates, SEXP parameters, SEXP options,
                     SEXP initial, SEXP seed, SEXP threads) {
  SaomHandle *h = saom_from_ptr(ptr);
  Spec s = spec_of(codes, covariates, parameters);
  if (TYPEOF(options) != REALSXP || XLENGTH(options) != SAOM_OPTIONS_LEN) {
    Rf_error("options must be a double vector of length %d (internal error)", SAOM_OPTIONS_LEN);
  }
  if (!Rf_isNull(initial) && TYPEOF(initial) != REALSXP) {
    Rf_error("initial values must be a double vector (internal error)");
  }
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, fit_finalizer, TRUE);
  SaomFit *fit = NULL;
  char err[ERR_LEN] = {0};
  if (net_saom_estimate(h, s.k, s.codes, s.covariates, s.parameters, REAL(options),
                        Rf_isNull(initial) ? NULL : REAL(initial), as_u64(seed),
                        Rf_asInteger(threads), &fit, err, ERR_LEN) != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  R_SetExternalPtrAddr(holder, fit);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, FIT_ITEMS));
  for (int which = 0; which < FIT_ITEMS; which++) {
    int64_t len = net_saom_fit_length(fit, which);
    if (len < 0) {
      UNPROTECT(2);
      Rf_error("fit item %d is missing (internal error)", which);
    }
    SEXP item = Rf_allocVector(REALSXP, (R_xlen_t)len);
    SET_VECTOR_ELT(out, which, item);
    if (len > 0 && net_saom_fit_copy(fit, which, REAL(item), err, ERR_LEN) != 0) {
      UNPROTECT(2);
      Rf_error("%s", err);
    }
  }
  net_saom_fit_free(fit);
  R_ClearExternalPtr(holder);
  UNPROTECT(2);
  return out;
}
