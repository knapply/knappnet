/*
 * R bindings for key players, group centrality, and flow betweenness. As
 * in glue.c, R errors are raised only after the engine call has returned.
 */
#include "glue.h"

SEXP C_kpp_score(SEXP ptr, SEXP set, SEXP measure, SEXP m, SEXP mode, SEXP weighted,
                 SEXP weights_as) {
  GraphHandle *h = graph_from_ptr(ptr);
  const int64_t *nodes = as_int64(set);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status = net_kpp_score(h, nodes, (int64_t)XLENGTH(set), Rf_asInteger(measure), Rf_asReal(m),
                             Rf_asInteger(mode), Rf_asLogical(weighted), Rf_asInteger(weights_as),
                             &value, err, ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}

/* Returns list(set (0-based integers), score, evaluations). */
SEXP C_key_players(SEXP ptr, SEXP k, SEXP measure, SEXP m, SEXP mode, SEXP weighted,
                   SEXP weights_as, SEXP method, SEXP restarts, SEXP seed, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  int64_t size = (int64_t)Rf_asReal(k);
  if (size < 1) Rf_error("the number of key players must be at least 1 (internal error)");
  int64_t *members = (int64_t *)R_alloc((size_t)size, sizeof(int64_t));
  double score = 0.0, evaluations = 0.0;
  char err[ERR_LEN] = {0};
  int status =
      net_key_players(h, size, Rf_asInteger(measure), Rf_asReal(m), Rf_asInteger(mode),
                      Rf_asLogical(weighted), Rf_asInteger(weights_as), Rf_asInteger(method),
                      (int64_t)Rf_asReal(restarts), (uint64_t)Rf_asReal(seed),
                      Rf_asInteger(threads), members, &score, &evaluations, err, ERR_LEN);
  check_status(status, err, 0);
  SEXP set = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)size));
  for (int64_t i = 0; i < size; i++) INTEGER(set)[i] = (int)members[i];
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, set);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(score));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(evaluations));
  UNPROTECT(2);
  return out;
}

SEXP C_group_centrality(SEXP ptr, SEXP group, SEXP measure, SEXP mode, SEXP weighted,
                        SEXP weights_as, SEXP normalized, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  const int64_t *nodes = as_int64(group);
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status =
      net_group_centrality(h, nodes, (int64_t)XLENGTH(group), Rf_asInteger(measure),
                           Rf_asInteger(mode), Rf_asLogical(weighted), Rf_asInteger(weights_as),
                           Rf_asLogical(normalized), Rf_asInteger(threads), &value, err, ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}

SEXP C_flow_betweenness(SEXP ptr, SEXP capacity, SEXP scale, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_flow_betweenness(h, Rf_asLogical(capacity), Rf_asInteger(scale),
                                    Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}
