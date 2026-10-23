/*
 * R bindings for reach, edge betweenness, information centrality, the
 * bridging coefficient, regular and automorphic equivalence, blockmodels,
 * classical scaling, and geospatial networks. As in glue.c, R errors are
 * raised only after the engine call has returned, so R's non-local error
 * handling never crosses into compiled Rust code.
 */
#include "glue.h"

static void check_reals(SEXP v, const char *what) {
  if (TYPEOF(v) != REALSXP) Rf_error("%s must be a double vector (internal error)", what);
}

SEXP C_reach(SEXP ptr, SEXP k, SEXP mode, SEXP weighted, SEXP weights_as, SEXP normalized,
             SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_reach(h, Rf_asReal(k), Rf_asInteger(mode), Rf_asLogical(weighted),
                         Rf_asInteger(weights_as), Rf_asLogical(normalized), Rf_asInteger(threads),
                         REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* One value per edge, in the order of x$edges. */
SEXP C_edge_betweenness(SEXP ptr, SEXP weighted, SEXP weights_as, SEXP normalized, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_edge_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_edge_betweenness(h, Rf_asLogical(weighted), Rf_asInteger(weights_as),
                                    Rf_asLogical(normalized), Rf_asInteger(threads), REAL(out), err,
                                    ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_information_centrality(SEXP ptr, SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_information_centrality(h, Rf_asLogical(weighted), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_bridging_coefficient(SEXP ptr) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_bridging_coefficient(h, REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* seed is NULL or an integer vector of class codes. Returns
 * list(membership, rounds) with 0-based classes. */
SEXP C_regular_equivalence(SEXP ptr, SEXP seed, SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  int64_t *s = NULL;
  R_xlen_t ns = 0;
  if (seed != R_NilValue) {
    s = as_int64(seed);
    ns = XLENGTH(seed);
  }
  SEXP membership = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)net_graph_node_count(h)));
  int32_t rounds = 0;
  char err[ERR_LEN] = {0};
  int status = net_regular_equivalence(h, s, (int64_t)ns, Rf_asLogical(weighted),
                                       INTEGER(membership), &rounds, err, ERR_LEN);
  check_status(status, err, 1);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, membership);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger((int)rounds));
  UNPROTECT(2);
  return out;
}

/* Returns 0-based orbit numbers. max_steps is a double (it may exceed the
 * range of an R integer). */
SEXP C_automorphic_equivalence(SEXP ptr, SEXP max_steps) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {0};
  int status =
      net_automorphic_equivalence(h, (int64_t)Rf_asReal(max_steps), INTEGER(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* membership holds 0-based block codes forming k blocks in order of first
 * appearance. Returns list(density, sums, pairs, summary): three k-by-k
 * matrices and c(overall density, R-squared). */
SEXP C_block_model(SEXP ptr, SEXP membership, SEXP k, SEXP weighted, SEXP loops) {
  GraphHandle *h = graph_from_ptr(ptr);
  int64_t *m = as_int64(membership);
  R_xlen_t n = XLENGTH(membership);
  int kk = Rf_asInteger(k);
  if (kk < 0) Rf_error("negative block count (internal error)");
  SEXP density = PROTECT(Rf_allocMatrix(REALSXP, kk, kk));
  SEXP sums = PROTECT(Rf_allocMatrix(REALSXP, kk, kk));
  SEXP pairs = PROTECT(Rf_allocMatrix(REALSXP, kk, kk));
  SEXP summary = PROTECT(Rf_allocVector(REALSXP, 2));
  char err[ERR_LEN] = {0};
  int status =
      net_block_model(h, m, (int64_t)n, (int64_t)kk, Rf_asLogical(weighted), Rf_asLogical(loops),
                      REAL(density), REAL(sums), REAL(pairs), REAL(summary), err, ERR_LEN);
  check_status(status, err, 4);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, density);
  SET_VECTOR_ELT(out, 1, sums);
  SET_VECTOR_ELT(out, 2, pairs);
  SET_VECTOR_ELT(out, 3, summary);
  UNPROTECT(5);
  return out;
}

/* d is an n-by-n double matrix. Returns list(points (n-by-k, unused
 * columns NaN), dims, eigenvalues, gof). */
SEXP C_classical_mds(SEXP d, SEXP n, SEXP k) {
  check_reals(d, "the dissimilarity matrix");
  int nn = Rf_asInteger(n);
  int kk = Rf_asInteger(k);
  if (nn < 0 || kk < 0) Rf_error("negative size (internal error)");
  if (XLENGTH(d) != (R_xlen_t)nn * nn) Rf_error("matrix size mismatch (internal error)");
  SEXP points = PROTECT(Rf_allocMatrix(REALSXP, nn, kk));
  SEXP eigenvalues = PROTECT(Rf_allocVector(REALSXP, nn));
  SEXP gof = PROTECT(Rf_allocVector(REALSXP, 2));
  int32_t dims = 0;
  char err[ERR_LEN] = {0};
  int status = net_classical_mds(REAL(d), (int64_t)nn, (int64_t)kk, REAL(points), &dims,
                                 REAL(eigenvalues), REAL(gof), err, ERR_LEN);
  check_status(status, err, 3);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, points);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger((int)dims));
  SET_VECTOR_ELT(out, 2, eigenvalues);
  SET_VECTOR_ELT(out, 3, gof);
  UNPROTECT(4);
  return out;
}

/* d is n-by-n and points n-by-p, both double matrices. */
SEXP C_mds_stress(SEXP d, SEXP n, SEXP points, SEXP p) {
  check_reals(d, "the dissimilarity matrix");
  check_reals(points, "the configuration");
  int nn = Rf_asInteger(n);
  int pp = Rf_asInteger(p);
  if (nn < 0 || pp < 0) Rf_error("negative size (internal error)");
  if (XLENGTH(d) != (R_xlen_t)nn * nn || XLENGTH(points) != (R_xlen_t)nn * pp) {
    Rf_error("matrix size mismatch (internal error)");
  }
  double value = 0.0;
  char err[ERR_LEN] = {0};
  int status =
      net_mds_stress(REAL(d), (int64_t)nn, REAL(points), (int64_t)pp, &value, err, ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}

/* x and y hold one coordinate per node. One length per edge, in the order
 * of x$edges. */
SEXP C_edge_lengths(SEXP ptr, SEXP x, SEXP y, SEXP crs, SEXP radius) {
  GraphHandle *h = graph_from_ptr(ptr);
  check_reals(x, "x coordinates");
  check_reals(y, "y coordinates");
  if (XLENGTH(x) != XLENGTH(y)) Rf_error("coordinate length mismatch (internal error)");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_edge_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_edge_lengths(h, REAL(x), REAL(y), (int64_t)XLENGTH(x), Rf_asInteger(crs),
                                Rf_asReal(radius), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* Returns list(from, to, distance) with 0-based endpoints as doubles (node
 * counts may exceed the range of an R integer only in theory, but doubles
 * cost nothing here). With k > 0 the k nearest neighbours of every point;
 * otherwise every pair within distance r. */
SEXP C_proximity(SEXP x, SEXP y, SEXP crs, SEXP radius, SEXP r, SEXP k, SEXP threads) {
  check_reals(x, "x coordinates");
  check_reals(y, "y coordinates");
  if (XLENGTH(x) != XLENGTH(y)) Rf_error("coordinate length mismatch (internal error)");
  EdgeListHandle *edges = NULL;
  char err[ERR_LEN] = {0};
  int status = net_proximity(REAL(x), REAL(y), (int64_t)XLENGTH(x), Rf_asInteger(crs),
                             Rf_asReal(radius), Rf_asReal(r), (int64_t)Rf_asReal(k),
                             Rf_asInteger(threads), &edges, err, ERR_LEN);
  check_status(status, err, 0);
  int64_t m = net_edgelist_size(edges);
  if (m < 0) {
    net_edgelist_free(edges);
    Rf_error("proximity network produced no edge list (internal error)");
  }
  int64_t *from = (int64_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int64_t));
  int64_t *to = (int64_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int64_t));
  SEXP distance = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)m));
  status = net_edgelist_copy(edges, from, to, REAL(distance), err, ERR_LEN);
  net_edgelist_free(edges);
  check_status(status, err, 1);
  SEXP from_r = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)m));
  SEXP to_r = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)m));
  for (int64_t i = 0; i < m; i++) {
    REAL(from_r)[i] = (double)from[i];
    REAL(to_r)[i] = (double)to[i];
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, from_r);
  SET_VECTOR_ELT(out, 1, to_r);
  SET_VECTOR_ELT(out, 2, distance);
  UNPROTECT(4);
  return out;
}

/* groups holds 0-based group codes forming k groups in order of first
 * appearance. Returns a k-by-5 matrix: size, centroid x, centroid y, mean
 * distance, standard distance. */
SEXP C_group_geometry(SEXP x, SEXP y, SEXP groups, SEXP crs, SEXP radius, SEXP k) {
  check_reals(x, "x coordinates");
  check_reals(y, "y coordinates");
  int64_t *g = as_int64(groups);
  R_xlen_t n = XLENGTH(groups);
  if (XLENGTH(x) != n || XLENGTH(y) != n) {
    Rf_error("coordinate length mismatch (internal error)");
  }
  int kk = Rf_asInteger(k);
  if (kk < 0) Rf_error("negative group count (internal error)");
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, kk, 5));
  double *o = REAL(out);
  char err[ERR_LEN] = {0};
  int status =
      net_group_geometry(REAL(x), REAL(y), g, (int64_t)n, Rf_asInteger(crs), Rf_asReal(radius),
                         (int64_t)kk, o, o + kk, o + 2 * kk, o + 3 * kk, o + 4 * kk, err, ERR_LEN);
  return vector_result(status, out, err);
}
