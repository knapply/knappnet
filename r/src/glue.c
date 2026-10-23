/*
 * R bindings for the network engine. Each function converts R vectors to
 * plain C arrays, calls the engine through netcapi.h, and converts results
 * back. Errors are raised only after the engine call has returned, so R's
 * non-local error handling never crosses into compiled Rust code.
 *
 * Engine messages are plain lower-case clauses; the R wrapper
 * engine_call() prefixes them with the R function's name and turns edge
 * numbers into rows of the edge table. The engine runs on its own threads,
 * which may not call R_CheckUserInterrupt(), so a computation cannot be
 * interrupted until it returns.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static void graph_finalizer(SEXP ptr) {
  GraphHandle *h = (GraphHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_graph_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static GraphHandle *graph_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid network handle (internal error)");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

static int64_t *as_int64(SEXP x) {
  R_xlen_t n = XLENGTH(x);
  int64_t *out = (int64_t *)R_alloc(n > 0 ? n : 1, sizeof(int64_t));
  const int *src = INTEGER(x);
  for (R_xlen_t i = 0; i < n; i++) out[i] = (int64_t)src[i];
  return out;
}

/* The external pointer is allocated, with its finalizer, before the engine
 * builds the graph, so no R allocation can fail while the graph is owned
 * by nothing. */
SEXP C_graph_new(SEXP n, SEXP src, SEXP tgt, SEXP weights, SEXP directed) {
  if (TYPEOF(src) != INTSXP || TYPEOF(tgt) != INTSXP || XLENGTH(src) != XLENGTH(tgt)) {
    Rf_error("edge indices must be integer vectors of equal length (internal error)");
  }
  if (!Rf_isNull(weights) && (TYPEOF(weights) != REALSXP || XLENGTH(weights) != XLENGTH(src))) {
    Rf_error("weights must be a double vector with one value per edge (internal error)");
  }
  R_xlen_t m = XLENGTH(src);
  const int64_t *s = as_int64(src);
  const int64_t *t = as_int64(tgt);
  const double *w = Rf_isNull(weights) ? NULL : REAL(weights);
  SEXP ptr = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, graph_finalizer, TRUE);

  GraphHandle *h = NULL;
  char err[ERR_LEN] = {0};
  int status = net_graph_new((int64_t)Rf_asInteger(n), s, t, w, (int64_t)m,
                             Rf_asLogical(directed), &h, err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  R_SetExternalPtrAddr(ptr, h);
  UNPROTECT(1);
  return ptr;
}

SEXP C_graph_alive(SEXP ptr) {
  return Rf_ScalarLogical(TYPEOF(ptr) == EXTPTRSXP && R_ExternalPtrAddr(ptr) != NULL);
}

SEXP C_degree(SEXP ptr, SEXP mode, SEXP loops, SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_degree(h, Rf_asInteger(mode), Rf_asLogical(loops),
                          Rf_asLogical(weighted), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_betweenness(SEXP ptr, SEXP weighted, SEXP weights_as, SEXP normalized, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_betweenness(h, Rf_asLogical(weighted), Rf_asInteger(weights_as),
                               Rf_asLogical(normalized), Rf_asInteger(threads), REAL(out),
                               err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_triad_census(SEXP ptr, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 16));
  char err[ERR_LEN] = {0};
  int status = net_triad_census(h, Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}
