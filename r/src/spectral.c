/*
 * R bindings for PageRank, Katz centrality, Bonacich power centrality, hubs
 * and authorities, assortativity, mixing matrices, the E-I index, and
 * within-group tie shares. As in glue.c, R errors are raised only after the
 * engine call has returned.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static GraphHandle *spectral_graph(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("Invalid network handle.");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

/* Category codes (an R integer vector) as 64-bit integers; NA becomes -1. */
static int64_t *spectral_codes(SEXP codes, R_xlen_t want) {
  if (TYPEOF(codes) != INTSXP) Rf_error("Category codes must be an integer vector.");
  if (XLENGTH(codes) != want) Rf_error("Category codes must have one value per node.");
  int64_t *out = (int64_t *)R_alloc(want > 0 ? (size_t)want : 1, sizeof(int64_t));
  const int *src = INTEGER(codes);
  for (R_xlen_t i = 0; i < want; i++) {
    out[i] = src[i] == NA_INTEGER ? -1 : (int64_t)src[i];
  }
  return out;
}

/* A double vector of length n, or NULL for R's NULL. */
static const double *spectral_optional(SEXP x, R_xlen_t n) {
  if (Rf_isNull(x)) return NULL;
  if (TYPEOF(x) != REALSXP || XLENGTH(x) != n) {
    Rf_error("Expected a numeric vector with one value per node.");
  }
  return REAL(x);
}

static SEXP spectral_vector_result(int status, SEXP out, const char *err) {
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_pagerank(SEXP ptr, SEXP damping, SEXP directed, SEXP weighted, SEXP personalized,
                SEXP dangling, SEXP loops, SEXP tol, SEXP max_iter, SEXP threads) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  const double *p = spectral_optional(personalized, n);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  char err[ERR_LEN] = {0};
  int status = net_pagerank(h, Rf_asReal(damping), Rf_asLogical(directed),
                            Rf_asLogical(weighted), p, Rf_asInteger(dangling),
                            Rf_asInteger(loops), Rf_asReal(tol), (int64_t)Rf_asReal(max_iter),
                            Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return spectral_vector_result(status, out, err);
}

SEXP C_katz(SEXP ptr, SEXP alpha, SEXP beta, SEXP mode, SEXP weighted, SEXP loops,
            SEXP normalized, SEXP tol, SEXP max_iter, SEXP threads) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  if (TYPEOF(beta) != REALSXP) Rf_error("`beta` must be numeric.");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  char err[ERR_LEN] = {0};
  int status = net_katz(h, Rf_asReal(alpha), REAL(beta), (int64_t)XLENGTH(beta),
                        Rf_asInteger(mode), Rf_asLogical(weighted), Rf_asInteger(loops),
                        Rf_asLogical(normalized), Rf_asReal(tol),
                        (int64_t)Rf_asReal(max_iter), Rf_asInteger(threads), REAL(out), err,
                        ERR_LEN);
  return spectral_vector_result(status, out, err);
}

SEXP C_power_centrality(SEXP ptr, SEXP exponent, SEXP mode, SEXP weighted, SEXP loops,
                        SEXP rescale, SEXP tol, SEXP max_iter, SEXP threads) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  char err[ERR_LEN] = {0};
  int status = net_bonacich_power(h, Rf_asReal(exponent), Rf_asInteger(mode),
                                  Rf_asLogical(weighted), Rf_asInteger(loops),
                                  Rf_asLogical(rescale), Rf_asReal(tol),
                                  (int64_t)Rf_asReal(max_iter), Rf_asInteger(threads),
                                  REAL(out), err, ERR_LEN);
  return spectral_vector_result(status, out, err);
}

/* Returns list(hubs, authorities, value, unique). */
SEXP C_hits(SEXP ptr, SEXP weighted, SEXP loops, SEXP scale, SEXP tol, SEXP max_iter,
            SEXP threads) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  SEXP hubs = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP auth = PROTECT(Rf_allocVector(REALSXP, n));
  double value = 0.0;
  int32_t unique = 0;
  char err[ERR_LEN] = {0};
  int status = net_hits(h, Rf_asLogical(weighted), Rf_asInteger(loops), Rf_asInteger(scale),
                        Rf_asReal(tol), (int64_t)Rf_asReal(max_iter), Rf_asInteger(threads),
                        REAL(hubs), REAL(auth), &value, &unique, err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(2);
    Rf_error("%s", err);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, hubs);
  SET_VECTOR_ELT(out, 1, auth);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(value));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(unique != 0));
  UNPROTECT(3);
  return out;
}

static SEXP spectral_scalar_result(int status, double value, const char *err) {
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarReal(value);
}

SEXP C_assortativity_nominal(SEXP ptr, SEXP types, SEXP k, SEXP directed, SEXP loops,
                             SEXP normalized) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *t = spectral_codes(types, n);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_assortativity_nominal(h, t, (int64_t)Rf_asInteger(k),
                                         Rf_asLogical(directed), Rf_asInteger(loops),
                                         Rf_asLogical(normalized), &value, err, ERR_LEN);
  return spectral_scalar_result(status, value, err);
}

SEXP C_assortativity_values(SEXP ptr, SEXP values, SEXP values_in, SEXP directed, SEXP loops,
                            SEXP normalized) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  const double *v = spectral_optional(values, n);
  const double *vin = spectral_optional(values_in, n);
  if (v == NULL) Rf_error("`values` must be numeric.");
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_assortativity_values(h, v, vin, Rf_asLogical(directed),
                                        Rf_asInteger(loops), Rf_asLogical(normalized), &value,
                                        err, ERR_LEN);
  return spectral_scalar_result(status, value, err);
}

SEXP C_assortativity_degree(SEXP ptr, SEXP weighted, SEXP directed, SEXP loops,
                            SEXP normalized) {
  GraphHandle *h = spectral_graph(ptr);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_assortativity_degree(h, Rf_asLogical(weighted), Rf_asLogical(directed),
                                        Rf_asInteger(loops), Rf_asLogical(normalized), &value,
                                        err, ERR_LEN);
  return spectral_scalar_result(status, value, err);
}

/* Returns list(counts, transposed): counts row-major. */
SEXP C_mixing_matrix(SEXP ptr, SEXP dim1, SEXP k1, SEXP dim2, SEXP k2, SEXP same,
                     SEXP direction) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *d1 = spectral_codes(dim1, n);
  int64_t *d2 = spectral_codes(dim2, n);
  int a = Rf_asInteger(k1);
  int b = Rf_asInteger(k2);
  if (a == NA_INTEGER || b == NA_INTEGER || a < 0 || b < 0) {
    Rf_error("Invalid number of categories.");
  }
  SEXP counts = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)a * (R_xlen_t)b));
  int32_t transposed = 0;
  char err[ERR_LEN] = {0};
  int status = net_mixing_matrix(h, d1, (int64_t)a, d2, (int64_t)b, Rf_asLogical(same),
                                 Rf_asInteger(direction), REAL(counts), &transposed, err,
                                 ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, counts);
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(transposed != 0));
  UNPROTECT(2);
  return out;
}

/* Returns list(index, external, internal). */
SEXP C_ei_index(SEXP ptr, SEXP groups, SEXP k, SEXP scope, SEXP direction, SEXP counting) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *g = spectral_codes(groups, n);
  int kk = Rf_asInteger(k);
  int sc = Rf_asInteger(scope);
  if (kk == NA_INTEGER || kk < 0) Rf_error("Invalid number of groups.");
  R_xlen_t size = sc == 0 ? 1 : (sc == 1 ? (R_xlen_t)kk : n);
  SEXP index = PROTECT(Rf_allocVector(REALSXP, size));
  SEXP external = PROTECT(Rf_allocVector(REALSXP, size));
  SEXP internal = PROTECT(Rf_allocVector(REALSXP, size));
  char err[ERR_LEN] = {0};
  int status = net_ei_index(h, g, (int64_t)kk, sc, Rf_asInteger(direction),
                            Rf_asInteger(counting), REAL(index), REAL(external),
                            REAL(internal), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(3);
    Rf_error("%s", err);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, index);
  SET_VECTOR_ELT(out, 1, external);
  SET_VECTOR_ELT(out, 2, internal);
  UNPROTECT(4);
  return out;
}

SEXP C_group_shares(SEXP ptr, SEXP groups, SEXP k, SEXP direction) {
  GraphHandle *h = spectral_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *g = spectral_codes(groups, n);
  int kk = Rf_asInteger(k);
  if (kk == NA_INTEGER || kk < 0) Rf_error("Invalid number of groups.");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5 * (R_xlen_t)kk));
  char err[ERR_LEN] = {0};
  int status = net_group_shares(h, g, (int64_t)kk, Rf_asInteger(direction), REAL(out), err,
                                ERR_LEN);
  return spectral_vector_result(status, out, err);
}
