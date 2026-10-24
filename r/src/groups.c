/*
 * R bindings for cohesive subsets (bicomponents, n-cliques, n-clans,
 * k-plexes, lambda sets, factions), community partitions (Girvan-Newman,
 * leading eigenvector, greedy modularity), whole-network structure
 * (core-periphery fits, compactness, the small-world quotient), and random
 * graphs. As in
 * glue.c, R errors are raised only after the engine call has returned and
 * every engine-owned object has been released.
 */
#include "glue.h"

/* Copies the 0-based labels in memb[0..n) into a 1-based integer vector. */
static SEXP labels_1based(const int64_t *memb, R_xlen_t n) {
  SEXP m = PROTECT(Rf_allocVector(INTSXP, n));
  for (R_xlen_t i = 0; i < n; i++) INTEGER(m)[i] = (int)memb[i] + 1;
  UNPROTECT(1);
  return m;
}

/* ---- Node sets ---------------------------------------------------------- */

static void node_sets_finalizer(SEXP ptr) {
  NodeSetsHandle *s = (NodeSetsHandle *)R_ExternalPtrAddr(ptr);
  if (s != NULL) {
    net_node_sets_free(s);
    R_ClearExternalPtr(ptr);
  }
}

/* Copies the node sets held by `holder` into a list of 1-based integer
 * vectors and releases them. The holder's finalizer guards against an R
 * allocation error leaking the engine's sets. */
static SEXP node_sets_to_list(SEXP holder) {
  NodeSetsHandle *sets = (NodeSetsHandle *)R_ExternalPtrAddr(holder);
  int64_t count = net_node_sets_count(sets);
  int64_t total = net_node_sets_total(sets);
  int64_t *offsets = (int64_t *)R_alloc((size_t)(count + 1), sizeof(int64_t));
  int64_t *members = (int64_t *)R_alloc((size_t)(total > 0 ? total : 1), sizeof(int64_t));
  char err[ERR_LEN] = {0};
  int status = net_node_sets_copy(sets, offsets, members, err, ERR_LEN);
  node_sets_finalizer(holder);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)count));
  for (int64_t i = 0; i < count; i++) {
    R_xlen_t len = (R_xlen_t)(offsets[i + 1] - offsets[i]);
    SEXP set = PROTECT(Rf_allocVector(INTSXP, len));
    int *dst = INTEGER(set);
    for (R_xlen_t j = 0; j < len; j++) dst[j] = (int)members[offsets[i] + j] + 1;
    SET_VECTOR_ELT(out, (R_xlen_t)i, set);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_bicomponents(SEXP ptr) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, node_sets_finalizer, TRUE);
  char err[ERR_LEN] = {0};
  NodeSetsHandle *sets = NULL;
  int status = net_bicomponents(h, &sets, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, sets);
  SEXP out = node_sets_to_list(holder);
  UNPROTECT(1);
  return out;
}

SEXP C_n_cliques(SEXP ptr, SEXP order, SEXP min_size, SEXP clans, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, node_sets_finalizer, TRUE);
  char err[ERR_LEN] = {0};
  NodeSetsHandle *sets = NULL;
  int status = net_n_cliques(h, as_count(order, "n"), as_count(min_size, "min_size"),
                             Rf_asLogical(clans), Rf_asInteger(threads), &sets, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, sets);
  SEXP out = node_sets_to_list(holder);
  UNPROTECT(1);
  return out;
}

SEXP C_k_plexes(SEXP ptr, SEXP k, SEXP min_size, SEXP max_steps, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, node_sets_finalizer, TRUE);
  char err[ERR_LEN] = {0};
  NodeSetsHandle *sets = NULL;
  int status = net_k_plexes(h, as_count(k, "k"), as_count(min_size, "min_size"),
                            Rf_asReal(max_steps), Rf_asInteger(threads), &sets, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, sets);
  SEXP out = node_sets_to_list(holder);
  UNPROTECT(1);
  return out;
}

/* An n-by-n matrix of pairwise edge connectivities. */
SEXP C_lambda_connectivity(SEXP ptr) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)n, (int)n));
  char err[ERR_LEN] = {0};
  int status = net_lambda_connectivity(h, REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* ---- Factions and partitions ------------------------------------------- */

/* Returns list(membership (1-based), errors, correlation). */
SEXP C_factions(SEXP ptr, SEXP k, SEXP fit, SEXP starts, SEXP max_stall, SEXP tenure, SEXP seed,
                SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  double errors = 0.0, correlation = 0.0;
  int status = net_factions(h, as_count(k, "k"), Rf_asInteger(fit), as_count(starts, "starts"),
                            as_count(max_stall, "max_stall"), as_count(tenure, "tenure"), s,
                            Rf_asInteger(threads), memb, &errors, &correlation, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(errors));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(correlation));
  UNPROTECT(1);
  return out;
}

/* Returns list(membership, modularity, removed_from, removed_to, path),
 * the removed edges as 1-based node positions. */
SEXP C_girvan_newman(SEXP ptr, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  R_xlen_t m = (R_xlen_t)net_graph_edge_count(h);
  int64_t *memb = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  int64_t *from = (int64_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int64_t));
  int64_t *to = (int64_t *)R_alloc(m > 0 ? (size_t)m : 1, sizeof(int64_t));
  double *path = (double *)R_alloc((size_t)m + 1, sizeof(double));
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int64_t removed = 0;
  int status =
      net_girvan_newman(h, Rf_asInteger(threads), memb, &q, from, to, path, &removed, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  SET_VECTOR_ELT(out, 2, labels_1based(from, (R_xlen_t)removed));
  SET_VECTOR_ELT(out, 3, labels_1based(to, (R_xlen_t)removed));
  SEXP p = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)removed + 1));
  for (R_xlen_t i = 0; i <= (R_xlen_t)removed; i++) REAL(p)[i] = path[i];
  SET_VECTOR_ELT(out, 4, p);
  UNPROTECT(2);
  return out;
}

/* Returns list(membership, modularity, splits). */
SEXP C_leading_eigen(SEXP ptr, SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int64_t splits = 0;
  int status = net_leading_eigen(h, Rf_asLogical(weighted), memb, &q, &splits, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger((int)splits));
  UNPROTECT(1);
  return out;
}

/* Returns list(membership, modularity). */
SEXP C_fast_greedy(SEXP ptr, SEXP resolution, SEXP weighted) {
  GraphHandle *h = graph_from_ptr(ptr);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int64_t *memb = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  double q = 0.0;
  int status =
      net_fast_greedy(h, Rf_asReal(resolution), Rf_asLogical(weighted), memb, &q, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, labels_1based(memb, n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  UNPROTECT(1);
  return out;
}

/* ---- Whole-network structure ------------------------------------------- */

/* Returns list(core (logical), correlation). */
SEXP C_core_periphery(SEXP ptr, SEXP links, SEXP weighted, SEXP starts, SEXP max_stall, SEXP tenure,
                      SEXP seed, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  int32_t *core = (int32_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int32_t));
  char err[ERR_LEN] = {0};
  double correlation = 0.0;
  int status =
      net_core_periphery(h, Rf_asInteger(links), Rf_asLogical(weighted), as_count(starts, "starts"),
                         as_count(max_stall, "max_stall"), as_count(tenure, "tenure"), s,
                         Rf_asInteger(threads), core, &correlation, err, ERR_LEN);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP c = PROTECT(Rf_allocVector(LGLSXP, n));
  for (R_xlen_t i = 0; i < n; i++) LOGICAL(c)[i] = core[i] != 0;
  SET_VECTOR_ELT(out, 0, c);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(correlation));
  UNPROTECT(2);
  return out;
}

/* Returns list(coreness, correlation, residual, iterations). fit: 0 = MINRES,
 * 1 = correlation. */
SEXP C_core_periphery_continuous(SEXP ptr, SEXP fit, SEXP weighted, SEXP max_iter, SEXP tol,
                                 SEXP starts, SEXP seed, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  R_xlen_t n = (R_xlen_t)net_graph_node_count(h);
  SEXP coreness = PROTECT(Rf_allocVector(REALSXP, n));
  char err[ERR_LEN] = {0};
  double correlation = 0.0, residual = 0.0;
  int64_t iterations = 0;
  int status = net_core_periphery_continuous(
      h, Rf_asInteger(fit), Rf_asLogical(weighted), as_count(max_iter, "max_iter"), Rf_asReal(tol),
      as_count(starts, "starts"), s, Rf_asInteger(threads), REAL(coreness), &correlation, &residual,
      &iterations, err, ERR_LEN);
  check_status(status, err, 1);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, coreness);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(correlation));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(residual));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal((double)iterations));
  UNPROTECT(2);
  return out;
}

SEXP C_compactness(SEXP ptr, SEXP mode, SEXP weighted, SEXP weights_as, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  char err[ERR_LEN] = {0};
  double value = 0.0;
  int status =
      net_compactness(h, Rf_asInteger(mode), Rf_asLogical(weighted), Rf_asInteger(weights_as),
                      Rf_asInteger(threads), &value, err, ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}

/* Five values: clustering, mean distance, random clustering, random mean
 * distance, quotient. */
SEXP C_small_world(SEXP ptr, SEXP local, SEXP reps, SEXP seed, SEXP threads) {
  GraphHandle *h = graph_from_ptr(ptr);
  uint64_t s = as_seed(seed);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  char err[ERR_LEN] = {0};
  int status = net_small_world(h, Rf_asLogical(local), as_count(reps, "reps"), s,
                               Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* ---- Random graphs ------------------------------------------------------ */

static void edge_lists_finalizer(SEXP ptr) {
  EdgeListsHandle *e = (EdgeListsHandle *)R_ExternalPtrAddr(ptr);
  if (e != NULL) {
    net_edge_lists_free(e);
    R_ClearExternalPtr(ptr);
  }
}

/* Copies the edge lists held by `holder` into a list of list(from, to)
 * (1-based node positions) and releases them. */
static SEXP edge_lists_to_list(SEXP holder) {
  EdgeListsHandle *lists = (EdgeListsHandle *)R_ExternalPtrAddr(holder);
  int64_t count = net_edge_lists_count(lists);
  int64_t total = net_edge_lists_total(lists);
  int64_t *offsets = (int64_t *)R_alloc((size_t)(count + 1), sizeof(int64_t));
  int64_t *from = (int64_t *)R_alloc((size_t)(total > 0 ? total : 1), sizeof(int64_t));
  int64_t *to = (int64_t *)R_alloc((size_t)(total > 0 ? total : 1), sizeof(int64_t));
  char err[ERR_LEN] = {0};
  int status = net_edge_lists_copy(lists, offsets, from, to, err, ERR_LEN);
  edge_lists_finalizer(holder);
  check_status(status, err, 0);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)count));
  for (int64_t i = 0; i < count; i++) {
    R_xlen_t len = (R_xlen_t)(offsets[i + 1] - offsets[i]);
    SEXP pair = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(pair, 0, labels_1based(from + offsets[i], len));
    SET_VECTOR_ELT(pair, 1, labels_1based(to + offsets[i], len));
    SET_VECTOR_ELT(out, (R_xlen_t)i, pair);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

/* Shared tail of every generator: the engine call has written `lists`
 * (or failed with `status`). */
static SEXP finish_edge_lists(SEXP holder, EdgeListsHandle *lists, int status, const char *err) {
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, lists);
  SEXP out = edge_lists_to_list(holder);
  UNPROTECT(1);
  return out;
}

static SEXP new_holder(void) {
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, edge_lists_finalizer, TRUE);
  return holder; /* still protected; finish_edge_lists unprotects it */
}

SEXP C_random_gnp(SEXP n, SEXP p, SEXP directed, SEXP loops, SEXP reps, SEXP seed, SEXP threads) {
  uint64_t s = as_seed(seed);
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status =
      net_random_gnp(as_count(n, "n"), Rf_asReal(p), Rf_asLogical(directed), Rf_asLogical(loops),
                     as_count(reps, "reps"), s, Rf_asInteger(threads), &lists, err, ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}

SEXP C_random_gnm(SEXP n, SEXP m, SEXP directed, SEXP loops, SEXP reps, SEXP seed, SEXP threads) {
  uint64_t s = as_seed(seed);
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status = net_random_gnm(as_count(n, "n"), as_count(m, "m"), Rf_asLogical(directed),
                              Rf_asLogical(loops), as_count(reps, "reps"), s, Rf_asInteger(threads),
                              &lists, err, ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}

/* degrees: integer vector of n degrees (or out-degrees); in_degrees: NULL
 * or n in-degrees. */
SEXP C_random_configuration(SEXP degrees, SEXP in_degrees, SEXP simple, SEXP max_tries, SEXP reps,
                            SEXP seed, SEXP threads) {
  uint64_t s = as_seed(seed);
  if (TYPEOF(degrees) != INTSXP) Rf_error("degrees must be integer (internal error)");
  R_xlen_t n = XLENGTH(degrees);
  int64_t *d = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  int64_t *din = NULL;
  for (R_xlen_t i = 0; i < n; i++) {
    if (INTEGER(degrees)[i] == NA_INTEGER || INTEGER(degrees)[i] < 0) {
      Rf_error("degrees must be non-negative whole numbers");
    }
    d[i] = INTEGER(degrees)[i];
  }
  if (in_degrees != R_NilValue) {
    if (TYPEOF(in_degrees) != INTSXP || XLENGTH(in_degrees) != n) {
      Rf_error("in_degrees must be an integer vector of the same length as degrees");
    }
    din = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
    for (R_xlen_t i = 0; i < n; i++) {
      if (INTEGER(in_degrees)[i] == NA_INTEGER || INTEGER(in_degrees)[i] < 0) {
        Rf_error("in_degrees must be non-negative whole numbers");
      }
      din[i] = INTEGER(in_degrees)[i];
    }
  }
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status = net_random_configuration((int64_t)n, d, din, Rf_asLogical(simple),
                                        as_count(max_tries, "max_tries"), as_count(reps, "reps"), s,
                                        Rf_asInteger(threads), &lists, err, ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}

SEXP C_random_uman(SEXP n, SEXP mutual, SEXP asym, SEXP reps, SEXP seed, SEXP threads) {
  uint64_t s = as_seed(seed);
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status =
      net_random_uman(as_count(n, "n"), as_count(mutual, "mutual"), as_count(asym, "asym"),
                      as_count(reps, "reps"), s, Rf_asInteger(threads), &lists, err, ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}

SEXP C_random_small_world(SEXP n, SEXP k, SEXP p, SEXP reps, SEXP seed, SEXP threads) {
  uint64_t s = as_seed(seed);
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status = net_random_small_world(as_count(n, "n"), as_count(k, "k"), Rf_asReal(p),
                                      as_count(reps, "reps"), s, Rf_asInteger(threads), &lists, err,
                                      ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}

SEXP C_random_preferential(SEXP n, SEXP m, SEXP reps, SEXP seed, SEXP threads) {
  uint64_t s = as_seed(seed);
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status = net_random_preferential(as_count(n, "n"), as_count(m, "m"), as_count(reps, "reps"),
                                       s, Rf_asInteger(threads), &lists, err, ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}

/* sizes: integer vector of k block sizes; probs: k-by-k numeric matrix
 * (column-major, as R stores it; probs[r, s] is read as row r, column s). */
SEXP C_random_sbm(SEXP sizes, SEXP probs, SEXP directed, SEXP loops, SEXP reps, SEXP seed,
                  SEXP threads) {
  uint64_t s = as_seed(seed);
  if (TYPEOF(sizes) != INTSXP) Rf_error("sizes must be integer (internal error)");
  if (TYPEOF(probs) != REALSXP) Rf_error("probs must be double (internal error)");
  R_xlen_t k = XLENGTH(sizes);
  if (XLENGTH(probs) != k * k) Rf_error("probs must be a k-by-k matrix (internal error)");
  int64_t *sz = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  double *pr = (double *)R_alloc(k * k > 0 ? (size_t)(k * k) : 1, sizeof(double));
  for (R_xlen_t i = 0; i < k; i++) {
    if (INTEGER(sizes)[i] == NA_INTEGER || INTEGER(sizes)[i] < 0) {
      Rf_error("sizes must be non-negative whole numbers");
    }
    sz[i] = INTEGER(sizes)[i];
  }
  /* The engine reads probs[r * k + s]: row-major. */
  for (R_xlen_t r = 0; r < k; r++) {
    for (R_xlen_t c = 0; c < k; c++) pr[r * k + c] = REAL(probs)[c * k + r];
  }
  SEXP holder = new_holder();
  char err[ERR_LEN] = {0};
  EdgeListsHandle *lists = NULL;
  int status =
      net_random_sbm((int64_t)k, sz, pr, Rf_asLogical(directed), Rf_asLogical(loops),
                     as_count(reps, "reps"), s, Rf_asInteger(threads), &lists, err, ERR_LEN);
  return finish_edge_lists(holder, lists, status, err);
}
