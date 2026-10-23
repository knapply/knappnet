/*
 * Shared declarations for the C glue between R and the network engine.
 *
 * Every C file under r/src includes this header. Each exported C_* function
 * converts R vectors to plain C arrays, calls the engine through
 * netcapi.h, and converts the result back. Errors are raised only after
 * the engine call has returned, so R's non-local error handling never
 * crosses into compiled Rust code; the helpers below unprotect before
 * raising. The engine runs on its own threads, which may not call
 * R_CheckUserInterrupt(), so a computation cannot be interrupted until it
 * returns.
 *
 * Engine messages are plain lower-case clauses; the R wrapper
 * engine_call() prefixes them with the R function's name and turns edge
 * numbers into rows of the edge table. Messages marked "(internal error)"
 * report a mismatch between the R code and the glue, which users should
 * never see.
 *
 * Conventions:
 * - Graphs and other engine objects travel as external pointers whose
 *   finalizer releases the engine handle; the pointer is allocated, with
 *   its finalizer, before the engine builds the object, so no R
 *   allocation can fail while the object is owned by nothing.
 * - Node indices are 0-based in the engine; the R code converts.
 * - Seeds and counts arrive as doubles holding whole numbers.
 * - PROTECT every allocation; UNPROTECT (with the same count) before
 *   raising or returning.
 */
#ifndef KNAPPNET_GLUE_H
#define KNAPPNET_GLUE_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>

#include "netcapi.h"

/* Size of the buffer that receives an engine error message. */
#define ERR_LEN 1024

/* The engine graph behind a network's external pointer. */
static inline GraphHandle *graph_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid network handle (internal error)");
  }
  return (GraphHandle *)R_ExternalPtrAddr(ptr);
}

/* Raises the engine's message when `status` is non-zero, after
 * unprotecting `nprotect` objects. */
static inline void check_status(int status, const char *err, int nprotect) {
  if (status != 0) {
    if (nprotect > 0) UNPROTECT(nprotect);
    Rf_error("%s", err);
  }
}

/* A protected result vector: raises the engine's message when `status` is
 * non-zero, and unprotects it either way. */
static inline SEXP vector_result(int status, SEXP out, const char *err) {
  check_status(status, err, 1);
  UNPROTECT(1);
  return out;
}

/* A scalar result, or the engine's message when `status` is non-zero. */
static inline SEXP scalar_result(int status, double value, const char *err) {
  check_status(status, err, 0);
  return Rf_ScalarReal(value);
}

/* Copies an integer vector into a 64-bit array that R reclaims. */
static inline int64_t *as_int64(SEXP x) {
  if (TYPEOF(x) != INTSXP) Rf_error("indices must be an integer vector (internal error)");
  R_xlen_t n = XLENGTH(x);
  int64_t *out = (int64_t *)R_alloc(n > 0 ? (size_t)n : 1, sizeof(int64_t));
  const int *src = INTEGER(x);
  for (R_xlen_t i = 0; i < n; i++) out[i] = (int64_t)src[i];
  return out;
}

/* A seed: a whole number in [0, 2^53], passed as a double. */
static inline uint64_t as_seed(SEXP seed) {
  double s = Rf_asReal(seed);
  if (!R_FINITE(s) || s < 0 || s > 9007199254740992.0 || floor(s) != s) {
    Rf_error("seed must be a whole number between 0 and 2^53");
  }
  return (uint64_t)s;
}

/* A count (iterations, replications, sizes): a non-negative whole number
 * passed as a double, named `what` in the error. */
static inline int64_t as_count(SEXP value, const char *what) {
  double v = Rf_asReal(value);
  if (!R_FINITE(v) || v < 0 || v > 9007199254740992.0 || floor(v) != v) {
    Rf_error("%s must be a non-negative whole number", what);
  }
  return (int64_t)v;
}

#endif
