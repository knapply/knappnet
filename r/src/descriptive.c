/*
 * R bindings for reciprocity, transitivity and clustering, distances and
 * the measures built on them, reachability, cut points and bridges, and
 * the counts summary() reports. As in glue.c, R errors are raised only
 * after the engine call has returned.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static GraphHandle *descriptive_graph(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("Invalid network handle.");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

/* 0-based node positions from an integer vector, or NULL (every node) when
 * `sources` is NULL. The count goes to *count (-1 for every node). */
static const int64_t *descriptive_sources(SEXP sources, int64_t *count) {
  if (Rf_isNull(sources)) {
    *count = -1;
    return NULL;
  }
  if (TYPEOF(sources) != INTSXP) Rf_error("Source nodes must be an integer vector.");
  R_xlen_t k = XLENGTH(sources);
  int64_t *out = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  const int *src = INTEGER(sources);
  for (R_xlen_t i = 0; i < k; i++) out[i] = (int64_t)src[i];
  *count = (int64_t)k;
  return out;
}

static SEXP descriptive_scalar_result(int status, double value, const char *err) {
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarReal(value);
}

SEXP C_reciprocity(SEXP ptr, SEXP mode) {
  GraphHandle *h = descriptive_graph(ptr);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_reciprocity(h, Rf_asInteger(mode), &value, err, ERR_LEN);
  return descriptive_scalar_result(status, value, err);
}

SEXP C_transitivity(SEXP ptr, SEXP undefined, SEXP threads) {
  GraphHandle *h = descriptive_graph(ptr);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_transitivity(h, Rf_asInteger(undefined), Rf_asInteger(threads), &value,
                                err, ERR_LEN);
  return descriptive_scalar_result(status, value, err);
}

/* Local coefficients (average FALSE) or their mean (average TRUE). */
SEXP C_clustering(SEXP ptr, SEXP average, SEXP weighting, SEXP undefined, SEXP threads) {
  GraphHandle *h = descriptive_graph(ptr);
  int avg = Rf_asLogical(average) == TRUE;
  R_xlen_t len = avg ? 1 : (R_xlen_t)net_graph_node_count(h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, len));
  char err[ERR_LEN] = {0};
  int status = net_clustering(h, avg, Rf_asInteger(weighting), Rf_asInteger(undefined),
                              Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* A (number of sources) x n matrix of distances. */
SEXP C_distances(SEXP ptr, SEXP sources, SEXP mode, SEXP weighted, SEXP weights_as,
                 SEXP threads) {
  GraphHandle *h = descriptive_graph(ptr);
  int64_t n = net_graph_node_count(h);
  int64_t count = 0;
  const int64_t *src = descriptive_sources(sources, &count);
  int64_t rows = count < 0 ? n : count;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)rows, (int)n));
  char err[ERR_LEN] = {0};
  int status = net_distances(h, src, count, Rf_asInteger(mode), Rf_asLogical(weighted),
                             Rf_asInteger(weights_as), Rf_asInteger(threads), REAL(out),
                             err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Seven values: sources, pairs joined by a path, sum of their distances,
 * longest distance, sources not reaching every other node, mean distance,
 * and diameter. */
SEXP C_path_summary(SEXP ptr, SEXP sources, SEXP mode, SEXP weighted, SEXP weights_as,
                    SEXP unconnected, SEXP threads) {
  GraphHandle *h = descriptive_graph(ptr);
  int64_t count = 0;
  const int64_t *src = descriptive_sources(sources, &count);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 7));
  char err[ERR_LEN] = {0};
  int status = net_path_summary(h, src, count, Rf_asInteger(mode), Rf_asLogical(weighted),
                                Rf_asInteger(weights_as), Rf_asLogical(unconnected),
                                Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_eccentricity(SEXP ptr, SEXP mode, SEXP weighted, SEXP weights_as, SEXP threads) {
  GraphHandle *h = descriptive_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_eccentricity(h, Rf_asInteger(mode), Rf_asLogical(weighted),
                                Rf_asInteger(weights_as), Rf_asInteger(threads), REAL(out),
                                err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_reachability(SEXP ptr, SEXP mode, SEXP threads) {
  GraphHandle *h = descriptive_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_reachability(h, Rf_asInteger(mode), Rf_asInteger(threads), REAL(out), err,
                                ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns list(cut_point, bridge): logical vectors over nodes and edges. */
SEXP C_cut_structure(SEXP ptr) {
  GraphHandle *h = descriptive_graph(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  R_xlen_t m = (R_xlen_t)net_graph_edge_count(h);
  SEXP cut = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP bridge = PROTECT(Rf_allocVector(LGLSXP, m));
  int32_t *cut_flags = (int32_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int32_t));
  int32_t *bridge_flags = (int32_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int32_t));
  char err[ERR_LEN] = {0};
  int status = net_cut_structure(h, cut_flags, bridge_flags, err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(2);
    Rf_error("%s", err);
  }
  for (R_xlen_t i = 0; i < n; i++) LOGICAL(cut)[i] = cut_flags[i] != 0;
  for (R_xlen_t i = 0; i < m; i++) LOGICAL(bridge)[i] = bridge_flags[i] != 0;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cut);
  SET_VECTOR_ELT(out, 1, bridge);
  UNPROTECT(3);
  return out;
}

/* Three values: self-loops, multiple edges, and isolates. */
SEXP C_edge_counts(SEXP ptr) {
  GraphHandle *h = descriptive_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 3));
  char err[ERR_LEN] = {0};
  int status = net_edge_counts(h, REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}
