/*
 * R bindings for density, components, fragmentation, centralization,
 * closeness, harmonic centrality, and eigenvector centrality, and for the
 * adjacency-matrix reading of a network that the sna and UCINET
 * conventions use. As in glue.c, R errors are raised only after the engine
 * call has returned.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static GraphHandle *topo_graph(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("Invalid network handle.");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

static void topo_finalizer(SEXP ptr) {
  GraphHandle *h = (GraphHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_graph_free(h);
    R_ClearExternalPtr(ptr);
  }
}

/* keep: 0 = largest weight of parallel edges, 1 = smallest. */
SEXP C_graph_simplify(SEXP ptr, SEXP undirected, SEXP keep) {
  GraphHandle *h = topo_graph(ptr);
  GraphHandle *out = NULL;
  char err[ERR_LEN] = {0};
  int status = net_graph_simplify(h, Rf_asLogical(undirected), Rf_asInteger(keep), &out,
                                  err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);
  SEXP res = PROTECT(R_MakeExternalPtr(out, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(res, topo_finalizer, TRUE);
  UNPROTECT(1);
  return res;
}

SEXP C_density(SEXP ptr, SEXP loops) {
  GraphHandle *h = topo_graph(ptr);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_density(h, Rf_asLogical(loops), &value, err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarReal(value);
}

/* Returns 0-based component labels. */
SEXP C_components(SEXP ptr, SEXP kind) {
  GraphHandle *h = topo_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)net_graph_node_count(h)));
  int64_t count = 0;
  char err[ERR_LEN] = {0};
  int status = net_components(h, Rf_asInteger(kind), INTEGER(out), &count, err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_fragmentation(SEXP ptr, SEXP kind) {
  GraphHandle *h = topo_graph(ptr);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_fragmentation(h, Rf_asInteger(kind), &value, err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarReal(value);
}

SEXP C_centralization(SEXP ptr, SEXP measure, SEXP mode, SEXP loops, SEXP normalized,
                      SEXP convention, SEXP threads) {
  GraphHandle *h = topo_graph(ptr);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_centralization(h, Rf_asInteger(measure), Rf_asInteger(mode),
                                  Rf_asLogical(loops), Rf_asLogical(normalized),
                                  Rf_asInteger(convention), Rf_asInteger(threads), &value,
                                  err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarReal(value);
}

SEXP C_closeness(SEXP ptr, SEXP harmonic, SEXP mode, SEXP weighted, SEXP weights_as,
                 SEXP normalized, SEXP unreachable, SEXP threads) {
  GraphHandle *h = topo_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_closeness(h, Rf_asLogical(harmonic), Rf_asInteger(mode),
                             Rf_asLogical(weighted), Rf_asInteger(weights_as),
                             Rf_asLogical(normalized), Rf_asInteger(unreachable),
                             Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns list(values, eigenvalue, connected). */
SEXP C_eigenvector(SEXP ptr, SEXP mode, SEXP weighted, SEXP scale, SEXP loops, SEXP tol,
                   SEXP max_iter, SEXP threads) {
  GraphHandle *h = topo_graph(ptr);
  SEXP values = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  double eigenvalue = 0.0;
  int32_t connected = 0;
  char err[ERR_LEN] = {0};
  int status = net_eigenvector(h, Rf_asInteger(mode), Rf_asLogical(weighted),
                               Rf_asInteger(scale), Rf_asInteger(loops), Rf_asReal(tol),
                               (int64_t)Rf_asReal(max_iter), Rf_asInteger(threads),
                               REAL(values), &eigenvalue, &connected, err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, values);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(eigenvalue));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(connected != 0));
  UNPROTECT(2);
  return out;
}
