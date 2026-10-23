/*
 * R bindings for the quadratic assignment procedure. Matrices arrive as
 * double vectors in R's column-major layout. Errors are raised only after
 * the engine call has returned.
 */
#include "glue.h"

/* Seeds arrive as doubles holding whole numbers in [0, 2^53), checked in R. */

/* `draws` is the number of permutations whose statistics are returned:
 * reps, n! when exact, or 0 to return none. */
SEXP C_qap_cor(SEXP y, SEXP x, SEXP n, SEXP mode, SEXP diag, SEXP exact, SEXP reps, SEXP seed,
               SEXP threads, SEXP draws) {
  if (TYPEOF(y) != REALSXP || TYPEOF(x) != REALSXP) {
    Rf_error("matrices must be double vectors (internal error)");
  }
  double nd = Rf_asReal(draws);
  SEXP summary = PROTECT(Rf_allocVector(REALSXP, 7));
  SEXP dist = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nd));
  char err[ERR_LEN] = {0};
  int status =
      net_qap_cor(REAL(y), REAL(x), (int64_t)Rf_asInteger(n), Rf_asInteger(mode),
                  Rf_asLogical(diag), Rf_asLogical(exact), (int64_t)Rf_asReal(reps), as_seed(seed),
                  Rf_asInteger(threads), REAL(summary), nd > 0 ? REAL(dist) : NULL, err, ERR_LEN);
  check_status(status, err, 2);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, summary);
  SET_VECTOR_ELT(out, 1, dist);
  UNPROTECT(3);
  return out;
}

SEXP C_mrqap(SEXP y, SEXP x, SEXP n, SEXP k, SEXP intercept, SEXP mode, SEXP diag, SEXP method,
             SEXP statistic, SEXP tol, SEXP exact, SEXP reps, SEXP seed, SEXP threads, SEXP draws) {
  if (TYPEOF(y) != REALSXP || TYPEOF(x) != REALSXP) {
    Rf_error("matrices must be double vectors (internal error)");
  }
  /* `draws` is as in C_qap_cor. */
  int p = Rf_asInteger(k) + (Rf_asLogical(intercept) == 1);
  double nd = Rf_asReal(draws);
  SEXP coef = PROTECT(Rf_allocMatrix(REALSXP, p, 7));
  SEXP fit = PROTECT(Rf_allocVector(REALSXP, 10));
  SEXP dist = PROTECT(Rf_allocMatrix(REALSXP, (int)nd, p));
  char err[ERR_LEN] = {0};
  int status = net_mrqap(
      REAL(y), REAL(x), (int64_t)Rf_asInteger(n), (int64_t)Rf_asInteger(k), Rf_asLogical(intercept),
      Rf_asInteger(mode), Rf_asLogical(diag), Rf_asInteger(method), Rf_asInteger(statistic),
      Rf_asReal(tol), Rf_asLogical(exact), (int64_t)Rf_asReal(reps), as_seed(seed),
      Rf_asInteger(threads), REAL(coef), REAL(fit), nd > 0 ? REAL(dist) : NULL, err, ERR_LEN);
  check_status(status, err, 3);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, coef);
  SET_VECTOR_ELT(out, 1, fit);
  SET_VECTOR_ELT(out, 2, dist);
  UNPROTECT(4);
  return out;
}
