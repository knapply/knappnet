/*
 * R bindings for network layouts. Each returns a double vector of length
 * 2n holding the x coordinates followed by the y coordinates; the R side
 * shapes it into a matrix. As in glue.c, R errors are raised only after
 * the engine call has returned.
 */
#include "glue.h"

static SEXP layout_alloc(GraphHandle *h, R_xlen_t *n) {
  *n = (R_xlen_t)net_graph_node_count(h);
  return Rf_allocVector(REALSXP, 2 * *n);
}

SEXP C_layout_stress(SEXP ptr, SEXP algorithm, SEXP weighted, SEXP weights_as, SEXP pivots,
                     SEXP max_iter, SEXP tol, SEXP seed, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  int64_t k = as_count(pivots, "pivots");
  int64_t iter = as_count(max_iter, "max_iter");
  R_xlen_t n = 0;
  SEXP out = PROTECT(layout_alloc(h, &n));
  char err[ERR_LEN] = {0};
  int status = net_layout_stress(h, Rf_asInteger(algorithm), Rf_asLogical(weighted),
                                 Rf_asInteger(weights_as), k, iter, Rf_asReal(tol), s,
                                 Rf_asInteger(threads), REAL(out), REAL(out) + n, err, ERR_LEN);
  UNPROTECT(1);
  check_status(status, err, 0);
  return out;
}

SEXP C_layout_fr(SEXP ptr, SEXP weighted, SEXP weights_as, SEXP iterations, SEXP approximation,
                 SEXP theta, SEXP seed, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  int64_t iter = as_count(iterations, "max_iter");
  R_xlen_t n = 0;
  SEXP out = PROTECT(layout_alloc(h, &n));
  char err[ERR_LEN] = {0};
  int status = net_layout_fr(h, Rf_asLogical(weighted), Rf_asInteger(weights_as), iter,
                             Rf_asInteger(approximation), Rf_asReal(theta), s,
                             Rf_asInteger(threads), REAL(out), REAL(out) + n, err, ERR_LEN);
  UNPROTECT(1);
  check_status(status, err, 0);
  return out;
}

/* kind: 0 = circle, 1 = star (param: 0-based center or -1), 2 = grid
 * (param: width or 0). */
SEXP C_layout_simple(SEXP ptr, SEXP kind, SEXP param) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = 0;
  SEXP out = PROTECT(layout_alloc(h, &n));
  char err[ERR_LEN] = {0};
  int status = net_layout_simple(h, Rf_asInteger(kind), (int64_t)Rf_asReal(param), REAL(out),
                                 REAL(out) + n, err, ERR_LEN);
  UNPROTECT(1);
  check_status(status, err, 0);
  return out;
}

/* roots: 0-based node indices as doubles (possibly empty). */
SEXP C_layout_tree(SEXP ptr, SEXP roots) {
  GraphHandle *h = graph_from_ptr(ptr);
  if (TYPEOF(roots) != REALSXP) Rf_error("roots must be a double vector (internal error)");
  R_xlen_t k = XLENGTH(roots);
  int64_t *r = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  for (R_xlen_t i = 0; i < k; i++) r[i] = (int64_t)REAL(roots)[i];
  R_xlen_t n = 0;
  SEXP out = PROTECT(layout_alloc(h, &n));
  char err[ERR_LEN] = {0};
  int status = net_layout_tree(h, r, (int64_t)k, REAL(out), REAL(out) + n, err, ERR_LEN);
  UNPROTECT(1);
  check_status(status, err, 0);
  return out;
}

/* types: integer vector of 0 (top row) and 1 (bottom row), one per node. */
SEXP C_layout_twomode(SEXP ptr, SEXP types, SEXP sweeps) {
  GraphHandle *h = graph_from_ptr(ptr);
  if (TYPEOF(types) != INTSXP) Rf_error("node types must be an integer vector (internal error)");
  int64_t sw = as_count(sweeps, "sweeps");
  R_xlen_t n = 0;
  SEXP out = PROTECT(layout_alloc(h, &n));
  if (XLENGTH(types) != n) {
    UNPROTECT(1);
    Rf_error("node types must have one value per node (internal error)");
  }
  char err[ERR_LEN] = {0};
  int status = net_layout_twomode(h, INTEGER(types), sw, REAL(out), REAL(out) + n, err, ERR_LEN);
  UNPROTECT(1);
  check_status(status, err, 0);
  return out;
}

/* xy: double vector of length 2n (x coordinates, then y). */
SEXP C_layout_stress_value(SEXP ptr, SEXP xy, SEXP weighted, SEXP weights_as, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  if (TYPEOF(xy) != REALSXP || XLENGTH(xy) != 2 * n) {
    Rf_error("the layout must have two numeric coordinates per node (internal error)");
  }
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_layout_stress_value(h, REAL(xy), REAL(xy) + n, Rf_asLogical(weighted),
                                       Rf_asInteger(weights_as), Rf_asInteger(threads), &value, err,
                                       ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}
