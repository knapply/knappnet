/*
 * R bindings for temporal networks. As in glue.c, R errors are raised only
 * after the engine call has returned. Engine-owned objects (the network
 * handle and variable-length results) are attached to external pointers
 * with finalizers before the engine creates them, so an R allocation
 * failure can never leak them.
 *
 * Path options arrive from R as one double vector:
 *   c(start, end, min_wait, max_wait, strict, traversal, traversal_time,
 *     threads)
 * with NA (NaN) start or end for the default window and Inf max_wait for
 * no limit; traversal is 0 (constant traversal_time) or 1 (durations).
 */
#include "glue.h"
#define OPTS_LEN 8

static void temporal_finalizer(SEXP ptr) {
  TemporalHandle *h = (TemporalHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_temporal_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static void result_finalizer(SEXP ptr) {
  TemporalResult *r = (TemporalResult *)R_ExternalPtrAddr(ptr);
  if (r != NULL) {
    net_temporal_result_free(r);
    R_ClearExternalPtr(ptr);
  }
}

static TemporalHandle *temporal_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid temporal network handle (internal error)");
  }
  return (TemporalHandle *)R_ExternalPtrAddr(ptr);
}

static const double *path_options(SEXP opts) {
  if (TYPEOF(opts) != REALSXP || XLENGTH(opts) != OPTS_LEN) {
    Rf_error("path options must be a double vector of length 8 (internal error)");
  }
  return REAL(opts);
}

/* Copies every column of the result held by `holder` into an R list of
 * double vectors, then releases the result. */
static SEXP result_columns(SEXP holder) {
  TemporalResult *r = (TemporalResult *)R_ExternalPtrAddr(holder);
  int ncol = 0;
  while (net_temporal_result_length(r, ncol) >= 0) ncol++;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, ncol));
  char err[ERR_LEN] = {0};
  for (int c = 0; c < ncol; c++) {
    int64_t k = net_temporal_result_length(r, c);
    SEXP col = Rf_allocVector(REALSXP, (R_xlen_t)k);
    SET_VECTOR_ELT(out, c, col);
    if (k > 0 && net_temporal_result_copy(r, c, REAL(col), err, ERR_LEN) != 0) {
      UNPROTECT(1);
      Rf_error("%s", err);
    }
  }
  net_temporal_result_free(r);
  R_ClearExternalPtr(holder);
  UNPROTECT(1);
  return out;
}

static SEXP new_holder(void) {
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, result_finalizer, TRUE);
  UNPROTECT(1);
  return holder;
}

/* src, tgt: 0-based integer endpoints; start: double; end: double or NULL. */
SEXP C_temporal_new(SEXP n, SEXP src, SEXP tgt, SEXP start, SEXP end, SEXP directed) {
  if (TYPEOF(src) != INTSXP || TYPEOF(tgt) != INTSXP || XLENGTH(src) != XLENGTH(tgt)) {
    Rf_error("edge indices must be integer vectors of equal length (internal error)");
  }
  R_xlen_t m = XLENGTH(src);
  if (TYPEOF(start) != REALSXP || XLENGTH(start) != m) {
    Rf_error("start times must be a double vector with one value per edge (internal error)");
  }
  if (!Rf_isNull(end) && (TYPEOF(end) != REALSXP || XLENGTH(end) != m)) {
    Rf_error("end times must be a double vector with one value per edge (internal error)");
  }
  const int64_t *s = as_int64(src);
  const int64_t *t = as_int64(tgt);
  SEXP ptr = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, temporal_finalizer, TRUE);
  TemporalHandle *h = NULL;
  char err[ERR_LEN] = {0};
  int status = net_temporal_new((int64_t)Rf_asInteger(n), s, t, REAL(start),
                                Rf_isNull(end) ? NULL : REAL(end), (int64_t)m,
                                Rf_asLogical(directed), &h, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(ptr, h);
  UNPROTECT(1);
  return ptr;
}

/* Returns a vector of length(roots) * n, row-major (row = root). */
SEXP C_temporal_paths(SEXP ptr, SEXP criterion, SEXP direction, SEXP roots, SEXP opts) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  const double *o = path_options(opts);
  const int64_t *r = as_int64(roots);
  R_xlen_t k = XLENGTH(roots);
  int64_t n = net_temporal_node_count(h);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, k * (R_xlen_t)n));
  char err[ERR_LEN] = {0};
  int status = net_temporal_paths(h, Rf_asInteger(criterion), Rf_asInteger(direction), r,
                                  (int64_t)k, o[0], o[1], o[2], o[3], (int32_t)o[4], (int32_t)o[5],
                                  o[6], (int32_t)o[7], REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

SEXP C_temporal_reach(SEXP ptr, SEXP direction, SEXP opts) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  const double *o = path_options(opts);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_temporal_node_count(h)));
  char err[ERR_LEN] = {0};
  int status = net_temporal_reach(h, Rf_asInteger(direction), o[0], o[1], o[2], o[3], (int32_t)o[4],
                                  (int32_t)o[5], o[6], (int32_t)o[7], REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}

/* Returns list(closeness, distances), distances NULL unless requested
 * (row-major, row = origin). */
SEXP C_temporal_closeness(SEXP ptr, SEXP opts, SEXP distances) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  const double *o = path_options(opts);
  R_xlen_t n = (R_xlen_t)net_temporal_node_count(h);
  int want = Rf_asLogical(distances) == 1;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP c = Rf_allocVector(REALSXP, n);
  SET_VECTOR_ELT(out, 0, c);
  SEXP d = R_NilValue;
  if (want) {
    d = Rf_allocVector(REALSXP, n * n);
    SET_VECTOR_ELT(out, 1, d);
  }
  char err[ERR_LEN] = {0};
  int status = net_temporal_closeness(h, o[0], o[1], o[2], o[3], (int32_t)o[4], (int32_t)o[5], o[6],
                                      (int32_t)o[7], REAL(c), want ? REAL(d) : NULL, err, ERR_LEN);
  return vector_result(status, out, err);
}

/* Returns list(offsets, indices) as doubles; indices are 0-based. */
SEXP C_temporal_windows(SEXP ptr, SEXP from, SEXP to) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  if (TYPEOF(from) != REALSXP || TYPEOF(to) != REALSXP || XLENGTH(from) != XLENGTH(to)) {
    Rf_error("window bounds must be double vectors of equal length (internal error)");
  }
  SEXP holder = PROTECT(new_holder());
  TemporalResult *r = NULL;
  char err[ERR_LEN] = {0};
  int status =
      net_temporal_windows(h, REAL(from), REAL(to), (int64_t)XLENGTH(from), &r, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, r);
  SEXP out = result_columns(holder);
  UNPROTECT(1);
  return out;
}

/* Returns 0-based indices (doubles) of the contacts active at time t. */
SEXP C_temporal_active_at(SEXP ptr, SEXP t) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  SEXP holder = PROTECT(new_holder());
  TemporalResult *r = NULL;
  char err[ERR_LEN] = {0};
  int status = net_temporal_active_at(h, Rf_asReal(t), &r, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, r);
  SEXP out = result_columns(holder);
  UNPROTECT(1);
  return VECTOR_ELT(out, 0);
}

/* window: NULL, or c(from, to). Returns list(from, to, count, duration,
 * first, last, pair_of) as doubles; endpoints and pairs are 0-based. */
SEXP C_temporal_aggregate(SEXP ptr, SEXP window) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  int windowed = !Rf_isNull(window);
  if (windowed && (TYPEOF(window) != REALSXP || XLENGTH(window) != 2)) {
    Rf_error("window must be NULL or a double vector of length 2 (internal error)");
  }
  double a = windowed ? REAL(window)[0] : 0.0;
  double b = windowed ? REAL(window)[1] : 0.0;
  SEXP holder = PROTECT(new_holder());
  TemporalResult *r = NULL;
  char err[ERR_LEN] = {0};
  int status = net_temporal_aggregate(h, windowed, a, b, &r, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, r);
  SEXP out = result_columns(holder);
  UNPROTECT(1);
  return out;
}

/* scope: 0 = nodes, 1 = pairs; mode: 1 = out, 2 = in, 3 = all. Returns
 * list(group, gap, events); groups are 0-based. */
SEXP C_temporal_inter_event(SEXP ptr, SEXP scope, SEXP mode) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  SEXP holder = PROTECT(new_holder());
  TemporalResult *r = NULL;
  char err[ERR_LEN] = {0};
  int status =
      net_temporal_inter_event(h, Rf_asInteger(scope), Rf_asInteger(mode), &r, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, r);
  SEXP out = result_columns(holder);
  UNPROTECT(1);
  return out;
}

/* Returns list(events, burstiness, memory). */
SEXP C_temporal_burstiness(SEXP ptr, SEXP scope, SEXP mode) {
  TemporalHandle *h = temporal_from_ptr(ptr);
  SEXP holder = PROTECT(new_holder());
  TemporalResult *r = NULL;
  char err[ERR_LEN] = {0};
  int status =
      net_temporal_burstiness(h, Rf_asInteger(scope), Rf_asInteger(mode), &r, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(holder, r);
  SEXP out = result_columns(holder);
  UNPROTECT(1);
  return out;
}
