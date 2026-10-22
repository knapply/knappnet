/*
 * R bindings for cohesive subgroups: cores, cliques, and communities.
 * As in glue.c, errors are raised only after the engine call has returned
 * and every engine-owned object has been released.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "netcapi.h"

#define ERR_LEN 1024

static GraphHandle *subgroups_graph(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("Invalid network handle.");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

SEXP C_coreness(SEXP ptr, SEXP mode) {
  GraphHandle *h = subgroups_graph(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_coreness(h, Rf_asInteger(mode), REAL(out), err, ERR_LEN);
  if (status != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  UNPROTECT(1);
  return out;
}

/* Returns a list of integer vectors of 1-based node positions. The groups
 * are copied into malloc'd buffers and the engine handle is freed before
 * any R allocation, so an R error cannot leak the handle. */
SEXP C_max_cliques(SEXP ptr, SEXP min_size, SEXP threads) {
  GraphHandle *h = subgroups_graph(ptr);
  char err[ERR_LEN] = {0};
  GroupsHandle *groups = NULL;
  int status = net_max_cliques(h, (int64_t)Rf_asInteger(min_size), Rf_asInteger(threads),
                               &groups, err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);

  int64_t count = net_groups_count(groups);
  int64_t total = net_groups_total(groups);
  int64_t *offsets = (int64_t *)malloc((size_t)(count + 1) * sizeof(int64_t));
  int64_t *members = (int64_t *)malloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
  if (offsets == NULL || members == NULL) {
    free(offsets);
    free(members);
    net_groups_free(groups);
    Rf_error("Not enough memory for %.0f cliques.", (double)count);
  }
  status = net_groups_copy(groups, offsets, members, err, ERR_LEN);
  net_groups_free(groups);
  if (status != 0) {
    free(offsets);
    free(members);
    Rf_error("%s", err);
  }

  /* From here on an R allocation error would leak the two buffers, so
   * copy them into R-managed memory first. */
  int64_t *roff = (int64_t *)R_alloc((size_t)(count + 1), sizeof(int64_t));
  int64_t *rmem = (int64_t *)R_alloc((size_t)(total > 0 ? total : 1), sizeof(int64_t));
  for (int64_t i = 0; i <= count; i++) roff[i] = offsets[i];
  for (int64_t i = 0; i < total; i++) rmem[i] = members[i];
  free(offsets);
  free(members);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)count));
  for (int64_t i = 0; i < count; i++) {
    R_xlen_t len = (R_xlen_t)(roff[i + 1] - roff[i]);
    SEXP clique = PROTECT(Rf_allocVector(INTSXP, len));
    int *dst = INTEGER(clique);
    for (R_xlen_t j = 0; j < len; j++) dst[j] = (int)rmem[roff[i] + j] + 1;
    SET_VECTOR_ELT(out, (R_xlen_t)i, clique);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_clique_number(SEXP ptr, SEXP threads) {
  GraphHandle *h = subgroups_graph(ptr);
  char err[ERR_LEN] = {0};
  int64_t value = 0;
  int status = net_clique_number(h, Rf_asInteger(threads), &value, err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarInteger((int)value);
}

/* membership: 0-based integer labels, one per node. */
SEXP C_modularity(SEXP ptr, SEXP membership, SEXP resolution, SEXP weighted) {
  GraphHandle *h = subgroups_graph(ptr);
  if (TYPEOF(membership) != INTSXP) Rf_error("Membership must be an integer vector.");
  R_xlen_t n = XLENGTH(membership);
  int64_t *memb = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  const int *src = INTEGER(membership);
  for (R_xlen_t i = 0; i < n; i++) memb[i] = (int64_t)src[i];
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int status = net_modularity(h, memb, (int64_t)n, Rf_asReal(resolution),
                              Rf_asLogical(weighted), &q, err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);
  return Rf_ScalarReal(q);
}

/* Returns list(membership = 1-based integer labels, modularity = double).
 * seed is a whole number in [0, 2^53], passed as a double. */
SEXP C_louvain(SEXP ptr, SEXP resolution, SEXP weighted, SEXP seed) {
  GraphHandle *h = subgroups_graph(ptr);
  double s = Rf_asReal(seed);
  if (!R_FINITE(s) || s < 0 || s > 9007199254740992.0 || floor(s) != s) {
    Rf_error("`seed` must be a whole number between 0 and 2^53.");
  }
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int status = net_louvain(h, Rf_asReal(resolution), Rf_asLogical(weighted), (uint64_t)s,
                           memb, &q, err, ERR_LEN);
  if (status != 0) Rf_error("%s", err);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP m = PROTECT(Rf_allocVector(INTSXP, n));
  for (R_xlen_t i = 0; i < n; i++) INTEGER(m)[i] = (int)memb[i] + 1;
  SET_VECTOR_ELT(out, 0, m);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  UNPROTECT(2);
  return out;
}
