/*
 * R bindings for reading and writing network files. Tables cross the
 * boundary as lists of plain vectors; the R side turns them into data
 * frames and handles R-specific types (factors, dates) before writing.
 *
 * Text crosses the boundary as UTF-8 in both directions: strings from R
 * are translated with Rf_translateCharUTF8() (R on Windows before 4.2
 * stores them in the native code page), and strings from the engine are
 * marked CE_UTF8, so non-ASCII identifiers round-trip on every platform.
 * File paths are passed as UTF-8 too; the engine opens them through
 * Rust's standard library, which uses the wide-character API on Windows.
 *
 * Network handles are owned by an external pointer from the moment they
 * exist, so an R error (an allocation failure, or a string that cannot be
 * translated) never leaks one.
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "netcapi.h"

#define ERR_LEN 1024

static void network_finalizer(SEXP ptr) {
  NetworkHandle *h = (NetworkHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_network_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static SEXP network_holder(void) {
  SEXP holder = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(holder, network_finalizer, TRUE);
  UNPROTECT(1);
  return holder;
}

static SEXP table_to_list(NetworkHandle *h, int which) {
  int64_t nrow = net_table_nrow(h, which);
  int ncol = net_table_ncol(h, which);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, ncol));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, ncol));
  char err[ERR_LEN] = {0};

  for (int c = 0; c < ncol; c++) {
    const char *name = net_table_name(h, which, c);
    SET_STRING_ELT(names, c, Rf_mkCharCE(name == NULL ? "" : name, CE_UTF8));
    int kind = net_table_kind(h, which, c);
    SEXP col;
    if (kind == 3) {
      col = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)nrow));
      for (int64_t r = 0; r < nrow; r++) {
        const char *s = net_table_string(h, which, c, r);
        SET_STRING_ELT(col, (R_xlen_t)r, s == NULL ? NA_STRING : Rf_mkCharCE(s, CE_UTF8));
      }
    } else {
      double *vals = (double *)R_alloc(nrow > 0 ? nrow : 1, sizeof(double));
      int32_t *na = (int32_t *)R_alloc(nrow > 0 ? nrow : 1, sizeof(int32_t));
      if (net_table_numeric(h, which, c, vals, na, err, ERR_LEN) != 0) {
        UNPROTECT(2);
        Rf_error("%s", err);
      }
      if (kind == 0) {
        col = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)nrow));
        for (int64_t r = 0; r < nrow; r++)
          LOGICAL(col)[r] = na[r] ? NA_LOGICAL : (vals[r] != 0);
      } else {
        /* Integers that fit R's 32-bit integers stay integer; others become
         * doubles so no value is lost. */
        int fits = (kind == 1);
        for (int64_t r = 0; fits && r < nrow; r++)
          if (!na[r] && (vals[r] > 2147483647.0 || vals[r] < -2147483647.0)) fits = 0;
        if (fits) {
          col = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)nrow));
          for (int64_t r = 0; r < nrow; r++)
            INTEGER(col)[r] = na[r] ? NA_INTEGER : (int)vals[r];
        } else {
          col = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nrow));
          for (int64_t r = 0; r < nrow; r++) REAL(col)[r] = na[r] ? NA_REAL : vals[r];
        }
      }
    }
    SET_VECTOR_ELT(out, c, col);
    UNPROTECT(1);
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP C_read_network(SEXP path, SEXP format, SEXP directed) {
  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 ||
      (!Rf_isNull(format) && (TYPEOF(format) != STRSXP || XLENGTH(format) != 1))) {
    Rf_error("path and format must be single strings (internal error)");
  }
  const char *fmt = Rf_isNull(format) ? NULL : Rf_translateCharUTF8(STRING_ELT(format, 0));
  const char *file = Rf_translateCharUTF8(STRING_ELT(path, 0));
  int dir = Rf_isNull(directed) ? -1 : Rf_asLogical(directed);
  if (dir == NA_LOGICAL) dir = -1;
  SEXP holder = PROTECT(network_holder());
  NetworkHandle *h = NULL;
  char err[ERR_LEN] = {0};
  if (net_read(file, fmt, dir, &h, err, ERR_LEN) != 0) {
    UNPROTECT(1);
    Rf_error("%s", err);
  }
  R_SetExternalPtrAddr(holder, h);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(net_network_directed(h)));
  SET_VECTOR_ELT(out, 1, table_to_list(h, 0));
  SET_VECTOR_ELT(out, 2, table_to_list(h, 1));
  SET_VECTOR_ELT(out, 3, table_to_list(h, 2));
  network_finalizer(holder);

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, Rf_mkChar("directed"));
  SET_STRING_ELT(names, 1, Rf_mkChar("graph"));
  SET_STRING_ELT(names, 2, Rf_mkChar("nodes"));
  SET_STRING_ELT(names, 3, Rf_mkChar("edges"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(3);
  return out;
}

/* Adds every column of a named list to table `which`. Returns 0 on
 * success; on failure writes the message to err. */
static int add_columns(NetworkHandle *h, int which, SEXP cols, char *err) {
  if (TYPEOF(cols) != VECSXP) {
    snprintf(err, ERR_LEN, "table %d is not a list (internal error)", which);
    return 1;
  }
  SEXP names = Rf_getAttrib(cols, R_NamesSymbol);
  for (R_xlen_t c = 0; c < XLENGTH(cols); c++) {
    SEXP col = VECTOR_ELT(cols, c);
    if (TYPEOF(names) != STRSXP || STRING_ELT(names, c) == NA_STRING ||
        CHAR(STRING_ELT(names, c))[0] == '\0') {
      snprintf(err, ERR_LEN, "column %d of table %d has no name", (int)c + 1, which);
      return 1;
    }
    const char *name = Rf_translateCharUTF8(STRING_ELT(names, c));
    R_xlen_t n = XLENGTH(col);
    int status;
    if (TYPEOF(col) == STRSXP) {
      const char **vals = (const char **)R_alloc(n > 0 ? n : 1, sizeof(char *));
      for (R_xlen_t i = 0; i < n; i++) {
        SEXP s = STRING_ELT(col, i);
        vals[i] = (s == NA_STRING) ? NULL : Rf_translateCharUTF8(s);
      }
      status = net_table_add_string(h, which, name, vals, (int64_t)n, err, ERR_LEN);
    } else {
      double *vals = (double *)R_alloc(n > 0 ? n : 1, sizeof(double));
      int32_t *na = (int32_t *)R_alloc(n > 0 ? n : 1, sizeof(int32_t));
      int kind;
      if (TYPEOF(col) == LGLSXP) {
        kind = 0;
        for (R_xlen_t i = 0; i < n; i++) {
          na[i] = LOGICAL(col)[i] == NA_LOGICAL;
          vals[i] = na[i] ? 0 : LOGICAL(col)[i];
        }
      } else if (TYPEOF(col) == INTSXP) {
        kind = 1;
        for (R_xlen_t i = 0; i < n; i++) {
          na[i] = INTEGER(col)[i] == NA_INTEGER;
          vals[i] = na[i] ? 0 : INTEGER(col)[i];
        }
      } else if (TYPEOF(col) == REALSXP) {
        kind = 2;
        for (R_xlen_t i = 0; i < n; i++) {
          na[i] = R_IsNA(REAL(col)[i]);
          vals[i] = na[i] ? 0 : REAL(col)[i];
        }
      } else {
        snprintf(err, ERR_LEN, "column `%s` has an unsupported type", name);
        return 1;
      }
      status = net_table_add_numeric(h, which, name, kind, vals, na, (int64_t)n, err, ERR_LEN);
    }
    if (status != 0) return status;
  }
  return 0;
}

SEXP C_write_network(SEXP directed, SEXP graph, SEXP nodes, SEXP edges, SEXP path,
                     SEXP format, SEXP weight) {
  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 ||
      (!Rf_isNull(format) && (TYPEOF(format) != STRSXP || XLENGTH(format) != 1)) ||
      (!Rf_isNull(weight) && (TYPEOF(weight) != STRSXP || XLENGTH(weight) != 1))) {
    Rf_error("path, format, and weight must be single strings (internal error)");
  }
  SEXP holder = PROTECT(network_holder());
  NetworkHandle *h = net_network_new(Rf_asLogical(directed) == 1);
  if (h == NULL) {
    UNPROTECT(1);
    Rf_error("could not create a network");
  }
  R_SetExternalPtrAddr(holder, h);
  char err[ERR_LEN] = {0};
  int status = add_columns(h, 0, graph, err);
  if (status == 0) status = add_columns(h, 1, nodes, err);
  if (status == 0) status = add_columns(h, 2, edges, err);
  if (status == 0) {
    const char *fmt = Rf_isNull(format) ? NULL : Rf_translateCharUTF8(STRING_ELT(format, 0));
    const char *w = Rf_isNull(weight) ? NULL : Rf_translateCharUTF8(STRING_ELT(weight, 0));
    const char *file = Rf_translateCharUTF8(STRING_ELT(path, 0));
    status = net_write(h, file, fmt, w, err, ERR_LEN);
  }
  network_finalizer(holder);
  UNPROTECT(1);
  if (status != 0) Rf_error("%s", err);
  return R_NilValue;
}
