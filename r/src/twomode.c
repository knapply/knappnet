/*
 * R bindings for ego networks and two-mode networks. As in glue.c, R
 * errors are raised only after the engine call has returned and every
 * engine-owned object has been released.
 */
#include <stdlib.h>

#include "glue.h"

/* 0-based node indices from an R integer vector, in R-managed memory. */
static int64_t *ego_indices(SEXP egos, R_xlen_t *k) {
  if (TYPEOF(egos) != INTSXP) Rf_error("Node indices must be an integer vector.");
  *k = XLENGTH(egos);
  int64_t *out = (int64_t *)R_alloc(*k > 0 ? (size_t)*k : 1, sizeof(int64_t));
  const int *src = INTEGER(egos);
  for (R_xlen_t i = 0; i < *k; i++) out[i] = (int64_t)src[i];
  return out;
}

/* Node types (0 or 1 per node) from an R integer vector. */
static const int32_t *node_types(SEXP types, GraphHandle *h) {
  if (TYPEOF(types) != INTSXP) Rf_error("Node types must be an integer vector.");
  if ((int64_t)XLENGTH(types) != net_graph_node_count(h)) {
    Rf_error("There must be one node type per node.");
  }
  return (const int32_t *)INTEGER(types);
}

/* egos: NULL for every node, or 0-based indices. */
SEXP C_ego_sizes(SEXP ptr, SEXP egos, SEXP order, SEXP mode, SEXP mindist, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t k = 0;
  const int64_t *list = NULL;
  if (Rf_isNull(egos)) {
    k = (R_xlen_t)net_graph_node_count(h);
  } else {
    list = ego_indices(egos, &k);
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, k));
  char err[ERR_LEN] = {0};
  int status =
      net_ego_sizes(h, list, (int64_t)k, (int64_t)Rf_asInteger(order), Rf_asInteger(mode),
                    (int64_t)Rf_asInteger(mindist), Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* Returns a list with one integer vector of 1-based node positions per
 * ego, in ascending order. */
SEXP C_ego_members(SEXP ptr, SEXP egos, SEXP order, SEXP mode, SEXP mindist, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t k = 0;
  const int64_t *list = ego_indices(egos, &k);
  int64_t ord = (int64_t)Rf_asInteger(order);
  int32_t md = Rf_asInteger(mode);
  int64_t mind = (int64_t)Rf_asInteger(mindist);
  int32_t thr = Rf_asInteger(threads);
  double *sizes = (double *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(double));
  char err[ERR_LEN] = {0};
  int status = net_ego_sizes(h, list, (int64_t)k, ord, md, mind, thr, sizes, err, ERR_LEN);
  check_status(status, err, 0);
  double total = 0.0;
  for (R_xlen_t i = 0; i < k; i++) total += sizes[i];
  int64_t *offsets = (int64_t *)R_alloc((size_t)k + 1, sizeof(int64_t));
  int64_t *members = (int64_t *)R_alloc(total > 0 ? (size_t)total : 1, sizeof(int64_t));
  status = net_ego_members(h, list, (int64_t)k, ord, md, mind, thr, offsets, members,
                           (int64_t)total, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, k));
  for (R_xlen_t i = 0; i < k; i++) {
    R_xlen_t len = (R_xlen_t)(offsets[i + 1] - offsets[i]);
    SEXP group = PROTECT(Rf_allocVector(INTSXP, len));
    int *dst = INTEGER(group);
    for (R_xlen_t j = 0; j < len; j++) dst[j] = (int)members[offsets[i] + j] + 1;
    SET_VECTOR_ELT(out, i, group);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

/* Returns list(size, ties, pairs, density, components, ego_betweenness);
 * the last is NULL unless betweenness is TRUE. */
SEXP C_ego_summary(SEXP ptr, SEXP order, SEXP mode, SEXP betweenness, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int with_betweenness = Rf_asLogical(betweenness) == TRUE;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  for (int i = 0; i < 5; i++) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, n));
  if (with_betweenness) SET_VECTOR_ELT(out, 5, Rf_allocVector(REALSXP, n));
  char err[ERR_LEN] = {0};
  int status =
      net_ego_summary(h, (int64_t)Rf_asInteger(order), Rf_asInteger(mode), Rf_asInteger(threads),
                      REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                      REAL(VECTOR_ELT(out, 3)), REAL(VECTOR_ELT(out, 4)),
                      with_betweenness ? REAL(VECTOR_ELT(out, 5)) : NULL, err, ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_twomode_density(SEXP ptr, SEXP types) {
  GraphHandle *h = graph_from_ptr(ptr);
  const int32_t *t = node_types(types, h);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_twomode_density(h, t, &value, err, ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}

/* measure: 0 = degree, 1 = closeness, 2 = betweenness. unreachable:
 * 0 = ignore, 1 = Wasserman-Faust. */
SEXP C_twomode_centrality(SEXP ptr, SEXP types, SEXP measure, SEXP normalized, SEXP unreachable,
                          SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  const int32_t *t = node_types(types, h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_twomode_centrality(h, t, Rf_asInteger(measure), Rf_asLogical(normalized),
                                      Rf_asInteger(unreachable), Rf_asInteger(threads), REAL(out),
                                      err, ERR_LEN);
  return vector_result(status, out, err);
}

/* Returns list(values, eigenvalue, connected). */
SEXP C_twomode_eigenvector(SEXP ptr, SEXP types, SEXP tol, SEXP max_iter, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  const int32_t *t = node_types(types, h);
  SEXP values = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  double eigenvalue = 0.0;
  int32_t connected = 0;
  char err[ERR_LEN] = {0};
  int status = net_twomode_eigenvector(h, t, Rf_asReal(tol), (int64_t)Rf_asReal(max_iter),
                                       Rf_asInteger(threads), REAL(values), &eigenvalue, &connected,
                                       err, ERR_LEN);
  check_status(status, err, 1);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, values);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(eigenvalue));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(connected != 0));
  UNPROTECT(2);
  return out;
}

/* Returns list(from, to, weight): 1-based node positions and weights. The
 * edges are copied into malloc'd buffers and the engine handle is freed
 * before any R allocation, so an R error cannot leak the handle. */
SEXP C_project(SEXP ptr, SEXP types, SEXP second, SEXP weight, SEXP loops, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  const int32_t *t = node_types(types, h);
  char err[ERR_LEN] = {0};
  ProjectionHandle *p = NULL;
  int status = net_project(h, t, Rf_asLogical(second), Rf_asInteger(weight), Rf_asLogical(loops),
                           Rf_asInteger(threads), &p, err, ERR_LEN);
  check_status(status, err, 0);

  int64_t m = net_projection_count(p);
  size_t cells = (size_t)(m > 0 ? m : 1);
  int64_t *from = (int64_t *)malloc(cells * sizeof(int64_t));
  int64_t *to = (int64_t *)malloc(cells * sizeof(int64_t));
  double *w = (double *)malloc(cells * sizeof(double));
  if (from == NULL || to == NULL || w == NULL) {
    free(from);
    free(to);
    free(w);
    net_projection_free(p);
    Rf_error("Not enough memory for a projection with %.0f edges.", (double)m);
  }
  status = net_projection_copy(p, from, to, w, err, ERR_LEN);
  net_projection_free(p);
  if (status != 0) {
    free(from);
    free(to);
    free(w);
    Rf_error("%s", err);
  }

  /* From here on an R allocation error would leak the buffers, so copy
   * them into R-managed memory first. */
  int64_t *rfrom = (int64_t *)R_alloc(cells, sizeof(int64_t));
  int64_t *rto = (int64_t *)R_alloc(cells, sizeof(int64_t));
  double *rw = (double *)R_alloc(cells, sizeof(double));
  for (int64_t i = 0; i < m; i++) {
    rfrom[i] = from[i];
    rto[i] = to[i];
    rw[i] = w[i];
  }
  free(from);
  free(to);
  free(w);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP a = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)m));
  SEXP b = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)m));
  SEXP x = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)m));
  for (int64_t i = 0; i < m; i++) {
    INTEGER(a)[i] = (int)rfrom[i] + 1;
    INTEGER(b)[i] = (int)rto[i] + 1;
    REAL(x)[i] = rw[i];
  }
  SET_VECTOR_ELT(out, 0, a);
  SET_VECTOR_ELT(out, 1, b);
  SET_VECTOR_ELT(out, 2, x);
  UNPROTECT(4);
  return out;
}
