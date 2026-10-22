/*
 * R bindings for brokerage and positions: Burt's structural holes,
 * Gould-Fernandez brokerage, structural equivalence, and CONCOR. Errors are
 * raised only after the engine call has returned, so R's non-local error
 * handling never crosses into compiled Rust code.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static GraphHandle *positions_graph(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("Invalid network handle.");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

SEXP C_structural_holes(SEXP ptr, SEXP measure, SEXP weighted, SEXP threads) {
  GraphHandle *h = positions_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_structural_holes(h, Rf_asInteger(measure), Rf_asLogical(weighted),
                                    Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_brokerage(SEXP ptr, SEXP groups, SEXP threads) {
  GraphHandle *h = positions_graph(ptr);
  if (TYPEOF(groups) != INTSXP) Rf_error("Group codes must be an integer vector.");
  int64_t n = net_graph_node_count(h);
  R_xlen_t ng = XLENGTH(groups);
  int64_t *g = (int64_t *)R_alloc(ng > 0 ? (size_t)ng : 1, sizeof(int64_t));
  const int *src = INTEGER(groups);
  for (R_xlen_t i = 0; i < ng; i++) g[i] = (int64_t)src[i];
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)n, 6));
  char err[ERR_LEN] = {0};
  int status = net_brokerage(h, g, (int64_t)ng, Rf_asInteger(threads), REAL(out), err,
                             ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_structural_equivalence(SEXP ptr, SEXP method, SEXP weighted, SEXP threads) {
  GraphHandle *h = positions_graph(ptr);
  int64_t n = net_graph_node_count(h);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)n, (int)n));
  char err[ERR_LEN] = {0};
  int status = net_structural_equivalence(h, Rf_asInteger(method), Rf_asLogical(weighted),
                                          Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns list(membership, diagnostics) where diagnostics holds the counts
 * of unconverged and numerically fragile splits. */
SEXP C_concor(SEXP ptr, SEXP depth, SEXP weighted, SEXP tolerance, SEXP max_iter,
              SEXP threads) {
  GraphHandle *h = positions_graph(ptr);
  int64_t n = net_graph_node_count(h);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP membership = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)n));
  SEXP diagnostics = PROTECT(Rf_allocVector(INTSXP, 2));
  char err[ERR_LEN] = {0};
  int status = net_concor(h, Rf_asInteger(depth), Rf_asLogical(weighted),
                          Rf_asReal(tolerance), Rf_asInteger(max_iter),
                          Rf_asInteger(threads), INTEGER(membership),
                          INTEGER(diagnostics), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(3);
    Rf_error("%s", err);
  }
  SET_VECTOR_ELT(out, 0, membership);
  SET_VECTOR_ELT(out, 1, diagnostics);
  UNPROTECT(3);
  return out;
}
