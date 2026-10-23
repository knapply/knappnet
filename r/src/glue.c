/*
 * The network handle: building a graph from an edge list, the degree,
 * betweenness, and triad census bindings, and the liveness check.
 */
#include "glue.h"

static void graph_finalizer(SEXP ptr) {
  GraphHandle *h = (GraphHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_graph_free(h);
    R_ClearExternalPtr(ptr);
  }
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
  int status = net_graph_new((int64_t)Rf_asInteger(n), s, t, w, (int64_t)m, Rf_asLogical(directed),
                             &h, err, ERR_LEN);
  check_status(status, err, 1);
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
  int status = net_degree(h, Rf_asInteger(mode), Rf_asLogical(loops), Rf_asLogical(weighted),
                          REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_betweenness(SEXP ptr, SEXP weighted, SEXP weights_as, SEXP normalized, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status =
      net_betweenness(h, Rf_asLogical(weighted), Rf_asInteger(weights_as), Rf_asLogical(normalized),
                      Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_triad_census(SEXP ptr, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 16));
  char err[ERR_LEN] = {0};
  int status = net_triad_census(h, Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}
