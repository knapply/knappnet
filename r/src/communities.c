/*
 * R bindings for community detection: Walktrap, label propagation,
 * Infomap and the map equation, Leiden, spinglass, exact maximum
 * modularity, fluid communities, Voronoi communities, partition
 * comparison, and dendrogram cutting. As in glue.c, R errors are raised
 * only after the engine call has returned.
 *
 * Memberships come back as 1-based integer vectors numbered by first
 * appearance. Dendrograms use igraph's `merges` numbering (leaves 1..n,
 * merge i creates cluster n + i) on the R side and 0-based indices in the
 * engine, so every index is shifted by one at this boundary.
 */
#include "glue.h"

/* Copies the 0-based labels in memb[0..n) into a 1-based integer vector. */
static SEXP labels_1based(const int64_t *memb, R_xlen_t n) {
  SEXP m = PROTECT(Rf_allocVector(INTSXP, n));
  for (R_xlen_t i = 0; i < n; i++) INTEGER(m)[i] = (int)memb[i] + 1;
  UNPROTECT(1);
  return m;
}

/* A buffer of `n` 64-bit values (at least one) that R reclaims. */
static int64_t *int64_buffer(R_xlen_t n) {
  return (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
}

/* An optional vector of n numbers: NULL passes through as NULL. */
static const double *optional_doubles(SEXP x, R_xlen_t n, const char *what) {
  if (Rf_isNull(x)) return NULL;
  if (TYPEOF(x) != REALSXP || XLENGTH(x) != n) {
    Rf_error("%s must be a numeric vector of the expected length (internal error)", what);
  }
  return REAL(x);
}

/* An optional vector of n integer codes as 64-bit values (NA becomes -1,
 * the engine's "unlabeled"): NULL passes through as NULL. */
static const int64_t *optional_codes(SEXP x, R_xlen_t n, const char *what) {
  if (Rf_isNull(x)) return NULL;
  if (TYPEOF(x) != INTSXP || XLENGTH(x) != n) {
    Rf_error("%s must be an integer vector of the expected length (internal error)", what);
  }
  int64_t *out = int64_buffer(n);
  const int *src = INTEGER(x);
  for (R_xlen_t i = 0; i < n; i++) out[i] = src[i] == NA_INTEGER ? -1 : (int64_t)src[i];
  return out;
}

/* Required integer codes of length n, as 64-bit values. */
static int64_t *required_codes(SEXP x, R_xlen_t n, const char *what) {
  if (TYPEOF(x) != INTSXP || XLENGTH(x) != n) {
    Rf_error("%s must be an integer vector of the expected length (internal error)", what);
  }
  int64_t *out = int64_buffer(n);
  const int *src = INTEGER(x);
  for (R_xlen_t i = 0; i < n; i++) {
    if (src[i] == NA_INTEGER) Rf_error("%s must not contain NA (internal error)", what);
    out[i] = (int64_t)src[i];
  }
  return out;
}

/* Returns list(membership, modularity, merges, path): the merges as an
 * integer matrix with two columns in igraph's numbering (1-based), the
 * modularity before any merge and after each as a numeric vector. */
SEXP C_walktrap(SEXP ptr, SEXP steps, SEXP weighted, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  int64_t t = as_count(steps, "steps");
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = int64_buffer(n);
  int64_t *ma = int64_buffer(n - 1);
  int64_t *mb = int64_buffer(n - 1);
  double *path = (double *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(double));
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int64_t k = 0;
  int status = net_walktrap(h, t, Rf_asLogical(weighted), Rf_asInteger(threads), memb, &q, ma, mb,
                            path, &k, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  SEXP merges = PROTECT(Rf_allocMatrix(INTSXP, (int)k, 2));
  for (R_xlen_t i = 0; i < (R_xlen_t)k; i++) {
    INTEGER(merges)[i] = (int)ma[i] + 1;
    INTEGER(merges)[i + k] = (int)mb[i] + 1;
  }
  SET_VECTOR_ELT(out, 2, merges);
  SEXP p = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)k + 1));
  for (R_xlen_t i = 0; i <= (R_xlen_t)k; i++) REAL(p)[i] = path[i];
  SET_VECTOR_ELT(out, 3, p);
  UNPROTECT(3);
  return out;
}

/* Returns the 1-based membership. `initial` holds integer codes with NA
 * for unlabeled nodes (or NULL), `fixed` logical flags (or NULL). */
SEXP C_label_propagation(SEXP ptr, SEXP variant, SEXP mode, SEXP weighted, SEXP initial, SEXP fixed,
                         SEXP seed) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  const int64_t *init = optional_codes(initial, n, "initial");
  int32_t *flags = NULL;
  if (!Rf_isNull(fixed)) {
    if (TYPEOF(fixed) != LGLSXP || XLENGTH(fixed) != n) {
      Rf_error("fixed must be a logical vector of the expected length (internal error)");
    }
    flags = (int32_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int32_t));
    for (R_xlen_t i = 0; i < n; i++) flags[i] = LOGICAL(fixed)[i] == TRUE;
  }
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  int status = net_label_propagation(h, Rf_asInteger(variant), Rf_asInteger(mode),
                                     Rf_asLogical(weighted), init, flags, s, memb, err, ERR_LEN);
  check_status(status, err, 0);
  return labels_1based(memb, n);
}

/* Returns list(membership, codelength). */
SEXP C_infomap(SEXP ptr, SEXP weighted, SEXP node_weights, SEXP trials, SEXP seed) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  int64_t tr = as_count(trials, "trials");
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  const double *nw = optional_doubles(node_weights, n, "node_weights");
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  double codelength = 0.0;
  int status = net_infomap(h, Rf_asLogical(weighted), nw, tr, s, memb, &codelength, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(codelength));
  UNPROTECT(1);
  return out;
}

/* The map equation of the partition `membership` (0-based integer codes). */
SEXP C_map_equation(SEXP ptr, SEXP membership, SEXP weighted, SEXP node_weights) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  const int64_t *labels = required_codes(membership, n, "membership");
  const double *nw = optional_doubles(node_weights, n, "node_weights");
  char err[ERR_LEN] = {0};
  double out = 0.0;
  int status = net_map_equation(h, labels, (int64_t)n, Rf_asLogical(weighted), nw, &out, err,
                                ERR_LEN);
  return scalar_result(status, out, err);
}

/* Returns list(membership, quality, iterations). `n_iterations` is a
 * whole number (negative runs until stable), passed as a double. */
SEXP C_leiden(SEXP ptr, SEXP objective, SEXP resolution, SEXP beta, SEXP n_iterations,
              SEXP weighted, SEXP node_weights, SEXP initial, SEXP seed) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  double iters = Rf_asReal(n_iterations);
  if (!R_FINITE(iters) || floor(iters) != iters || fabs(iters) > 9007199254740992.0) {
    Rf_error("n_iterations must be a whole number");
  }
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  const double *nw = optional_doubles(node_weights, n, "node_weights");
  const int64_t *init = NULL;
  if (!Rf_isNull(initial)) init = required_codes(initial, n, "initial");
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  double quality = 0.0;
  int64_t ran = 0;
  int status = net_leiden(h, Rf_asInteger(objective), Rf_asReal(resolution), Rf_asReal(beta),
                          (int64_t)iters, Rf_asLogical(weighted), nw, init, s, memb, &quality,
                          &ran, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(quality));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger((int)ran));
  UNPROTECT(1);
  return out;
}

/* Returns list(membership, modularity, energy). */
SEXP C_spinglass(SEXP ptr, SEXP spins, SEXP parallel, SEXP start_temp, SEXP stop_temp,
                 SEXP cool_fact, SEXP null_model, SEXP gamma, SEXP is_signed, SEXP gamma_minus,
                 SEXP weighted, SEXP seed) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  int64_t sp = as_count(spins, "spins");
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  double q = 0.0, energy = 0.0;
  int status = net_spinglass(h, sp, Rf_asLogical(parallel), Rf_asReal(start_temp),
                             Rf_asReal(stop_temp), Rf_asReal(cool_fact), Rf_asInteger(null_model),
                             Rf_asReal(gamma), Rf_asLogical(is_signed), Rf_asReal(gamma_minus),
                             Rf_asLogical(weighted), s, memb, &q, &energy, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(energy));
  UNPROTECT(1);
  return out;
}

/* Returns list(membership, modularity). */
SEXP C_optimal_modularity(SEXP ptr, SEXP resolution, SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int status = net_optimal_modularity(h, Rf_asReal(resolution), Rf_asLogical(weighted), memb, &q,
                                      err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  UNPROTECT(1);
  return out;
}

/* Returns the 1-based membership of k fluid communities. */
SEXP C_fluid_communities(SEXP ptr, SEXP k, SEXP seed) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  int64_t groups = as_count(k, "k");
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  int status = net_fluid_communities(h, groups, s, memb, err, ERR_LEN);
  check_status(status, err, 0);
  return labels_1based(memb, n);
}

/* Returns list(membership, generators, modularity, radius). `lengths` is
 * NULL or one number per edge; a negative radius chooses the radius of
 * highest modularity. */
SEXP C_voronoi(SEXP ptr, SEXP lengths, SEXP weighted, SEXP mode, SEXP radius, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  R_xlen_t m = (R_xlen_t)net_graph_edge_count(h);
  const double *len = optional_doubles(lengths, m, "lengths");
  int64_t *memb = int64_buffer(n);
  int64_t *gens = int64_buffer(n);
  char err[ERR_LEN] = {0};
  int64_t ngen = 0;
  double q = 0.0, used = 0.0;
  int status = net_voronoi(h, len, Rf_asLogical(weighted), Rf_asInteger(mode), Rf_asReal(radius),
                           Rf_asInteger(threads), memb, gens, &ngen, &q, &used, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, labels_1based(gens, (R_xlen_t)ngen));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(q));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(used));
  UNPROTECT(1);
  return out;
}

/* A distance or similarity between two partitions given as 0-based
 * integer codes of equal length. */
SEXP C_compare_partitions(SEXP a, SEXP b, SEXP method) {
  if (TYPEOF(a) != INTSXP || TYPEOF(b) != INTSXP || XLENGTH(a) != XLENGTH(b)) {
    Rf_error("partitions must be integer vectors of equal length (internal error)");
  }
  R_xlen_t n = XLENGTH(a);
  int64_t *x = required_codes(a, n, "a");
  int64_t *y = required_codes(b, n, "b");
  char err[ERR_LEN] = {0};
  double out = 0.0;
  int status = net_compare_partitions(x, y, (int64_t)n, Rf_asInteger(method), &out, err, ERR_LEN);
  return scalar_result(status, out, err);
}

/* Reads a dendrogram from two vectors of 1-based cluster indices and
 * returns it as 0-based 64-bit arrays of the same length. */
static void read_merges(SEXP a, SEXP b, int64_t **ma, int64_t **mb, R_xlen_t *len) {
  if (TYPEOF(a) != INTSXP || TYPEOF(b) != INTSXP || XLENGTH(a) != XLENGTH(b)) {
    Rf_error("merges must be two integer vectors of equal length (internal error)");
  }
  R_xlen_t k = XLENGTH(a);
  *ma = int64_buffer(k);
  *mb = int64_buffer(k);
  for (R_xlen_t i = 0; i < k; i++) {
    if (INTEGER(a)[i] == NA_INTEGER || INTEGER(b)[i] == NA_INTEGER) {
      Rf_error("merges must not contain NA (internal error)");
    }
    (*ma)[i] = (int64_t)INTEGER(a)[i] - 1;
    (*mb)[i] = (int64_t)INTEGER(b)[i] - 1;
  }
  *len = k;
}

/* Membership of n leaves after the first `steps` merges or, when
 * groups > 0, cut at exactly `groups` groups. */
SEXP C_cut_dendrogram(SEXP n, SEXP merge_a, SEXP merge_b, SEXP steps, SEXP groups) {
  int64_t leaves = as_count(n, "the number of leaves");
  int64_t st = as_count(steps, "steps");
  int64_t g = as_count(groups, "k");
  int64_t *ma, *mb;
  R_xlen_t k;
  read_merges(merge_a, merge_b, &ma, &mb, &k);
  int64_t *memb = int64_buffer((R_xlen_t)leaves);
  char err[ERR_LEN] = {0};
  int status = net_cut_dendrogram(leaves, ma, mb, (int64_t)k, st, g, memb, err, ERR_LEN);
  check_status(status, err, 0);
  return labels_1based(memb, (R_xlen_t)leaves);
}

/* Returns list(path, membership, modularity): the modularity of the
 * network after 0, 1, ... merges, and the first level of highest
 * modularity. */
SEXP C_dendrogram_modularity(SEXP ptr, SEXP merge_a, SEXP merge_b, SEXP resolution,
                             SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *ma, *mb;
  R_xlen_t k;
  read_merges(merge_a, merge_b, &ma, &mb, &k);
  double *path = (double *)R_alloc((size_t)k + 1, sizeof(double));
  int64_t *memb = int64_buffer(n);
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int status = net_dendrogram_modularity(h, ma, mb, (int64_t)k, Rf_asReal(resolution),
                                         Rf_asLogical(weighted), path, memb, &q, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP p = PROTECT(Rf_allocVector(REALSXP, k + 1));
  for (R_xlen_t i = 0; i <= k; i++) REAL(p)[i] = path[i];
  SET_VECTOR_ELT(out, 0, p);
  SET_VECTOR_ELT(out, 1, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(q));
  UNPROTECT(2);
  return out;
}
