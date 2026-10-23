/*
 * The data model: an engine frame built from the role columns of the node
 * and edge tables, validation with structured diagnostics, the graph
 * handle, the edge index, the derived summary, and the parallel-edge
 * operations (multiplicity, simplicity, simplify).
 *
 * The R side converts classed columns before the call (factors and 64-bit
 * integers to character, dates and date-times to seconds, logical
 * identifiers to integers), so this file sees only integer, double,
 * logical, and character vectors. Integers and logicals cross with a
 * missing mask; doubles carry NaN for missing; text crosses as UTF-8 bytes
 * with offsets (R on Windows before 4.2 stores strings in the native code
 * page, so every element goes through Rf_translateCharUTF8()).
 *
 * A frame is owned by an external pointer from the moment it exists, so an
 * R error raised while filling it never leaks the handle. Rows and node
 * indices come back 0-based; the R side converts.
 */
#include "glue.h"

#include <string.h>

static void frame_finalizer(SEXP ptr) {
  FrameHandle *h = (FrameHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_frame_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static void graph_release(SEXP ptr) {
  GraphHandle *h = (GraphHandle *)R_ExternalPtrAddr(ptr);
  if (h != NULL) {
    net_graph_free(h);
    R_ClearExternalPtr(ptr);
  }
}

static FrameHandle *frame_from_ptr(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == NULL) {
    Rf_error("invalid frame handle (internal error)");
  }
  return (FrameHandle *)R_ExternalPtrAddr(ptr);
}

static const char *utf8(SEXP s) { return Rf_translateCharUTF8(s); }

/* Adds the R vector `col` as column `name` of table `table`. Returns the
 * engine's status; `err` holds the message on failure. */
static int add_column(FrameHandle *h, int table, const char *name, SEXP col, char *err) {
  R_xlen_t n = XLENGTH(col);
  size_t alloc = n > 0 ? (size_t)n : 1;
  switch (TYPEOF(col)) {
  case INTSXP: {
    int64_t *values = (int64_t *)R_alloc(alloc, sizeof(int64_t));
    uint8_t *missing = (uint8_t *)R_alloc(alloc, sizeof(uint8_t));
    const int *v = INTEGER(col);
    int any = 0;
    for (R_xlen_t i = 0; i < n; i++) {
      missing[i] = (v[i] == NA_INTEGER);
      values[i] = missing[i] ? 0 : (int64_t)v[i];
      any |= missing[i];
    }
    return net_frame_add_int(h, table, name, values, any ? missing : NULL, (int64_t)n, err,
                             ERR_LEN);
  }
  case REALSXP:
    return net_frame_add_double(h, table, name, REAL(col), (int64_t)n, err, ERR_LEN);
  case LGLSXP: {
    int32_t *values = (int32_t *)R_alloc(alloc, sizeof(int32_t));
    uint8_t *missing = (uint8_t *)R_alloc(alloc, sizeof(uint8_t));
    const int *v = LOGICAL(col);
    int any = 0;
    for (R_xlen_t i = 0; i < n; i++) {
      missing[i] = (v[i] == NA_LOGICAL);
      values[i] = missing[i] ? 0 : (v[i] != 0);
      any |= missing[i];
    }
    return net_frame_add_logical(h, table, name, values, any ? missing : NULL, (int64_t)n,
                                 err, ERR_LEN);
  }
  case STRSXP: {
    int64_t *offsets = (int64_t *)R_alloc(alloc + 1, sizeof(int64_t));
    uint8_t *missing = (uint8_t *)R_alloc(alloc, sizeof(uint8_t));
    int any = 0;
    int64_t total = 0;
    for (R_xlen_t i = 0; i < n; i++) {
      SEXP s = STRING_ELT(col, i);
      missing[i] = (s == NA_STRING);
      any |= missing[i];
      offsets[i] = total;
      if (!missing[i]) total += (int64_t)strlen(utf8(s));
    }
    offsets[n] = total;
    char *bytes = (char *)R_alloc(total > 0 ? (size_t)total : 1, sizeof(char));
    for (R_xlen_t i = 0; i < n; i++) {
      if (missing[i]) continue;
      const char *s = utf8(STRING_ELT(col, i));
      memcpy(bytes + offsets[i], s, (size_t)(offsets[i + 1] - offsets[i]));
    }
    return net_frame_add_string(h, table, name, offsets, bytes, any ? missing : NULL,
                                (int64_t)n, err, ERR_LEN);
  }
  default:
    snprintf(err, ERR_LEN, "column `%s` has an unsupported type (internal error)", name);
    return 1;
  }
}

/* Adds every column of a named list to a table. */
static int add_columns(FrameHandle *h, int table, SEXP cols, char *err) {
  if (TYPEOF(cols) != VECSXP) {
    snprintf(err, ERR_LEN, "table columns must be a list (internal error)");
    return 1;
  }
  SEXP names = Rf_getAttrib(cols, R_NamesSymbol);
  R_xlen_t k = XLENGTH(cols);
  if (k > 0 && (TYPEOF(names) != STRSXP || XLENGTH(names) != k)) {
    snprintf(err, ERR_LEN, "table columns must be named (internal error)");
    return 1;
  }
  for (R_xlen_t c = 0; c < k; c++) {
    int status = add_column(h, table, utf8(STRING_ELT(names, c)), VECTOR_ELT(cols, c), err);
    if (status != 0) return status;
  }
  return 0;
}

/* A frame from the role columns. `node_cols` is NULL when there is no
 * node table (nodes then come from the edges); `node_rows` gives the row
 * count of a node table that may have no columns. `extra_cols` is a list
 * of three tables (layers, node spells, node changes), each NULL or a
 * named list of columns. `roles` is a list of three parallel vectors:
 * table codes, role codes, and column names. `declarations` holds the
 * values of declarations 1 to 7 (loops, multiple, between_modes_only,
 * add_nodes, egos_are_nodes, couplings_only, overlapping_spells; -1 = not
 * declared). `geometry` is a list of two flags and a CRS string or NULL;
 * `coordinates` the coordinate kind or NULL. */
SEXP C_frame_new(SEXP directed, SEXP node_rows, SEXP node_cols, SEXP edge_cols, SEXP extra_cols,
                 SEXP roles, SEXP weight_type, SEXP declarations, SEXP geometry,
                 SEXP coordinates) {
  SEXP ptr = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, frame_finalizer, TRUE);
  FrameHandle *h = net_frame_new(Rf_asLogical(directed));
  if (h == NULL) {
    UNPROTECT(1);
    Rf_error("could not allocate a frame (internal error)");
  }
  R_SetExternalPtrAddr(ptr, h);
  char err[ERR_LEN] = {0};
  int status = 0;
  if (!Rf_isNull(node_cols)) {
    status = net_frame_node_rows(h, (int64_t)Rf_asReal(node_rows), err, ERR_LEN);
    check_status(status, err, 1);
    status = add_columns(h, 0, node_cols, err);
    check_status(status, err, 1);
  }
  status = add_columns(h, 1, edge_cols, err);
  check_status(status, err, 1);
  if (TYPEOF(extra_cols) != VECSXP || XLENGTH(extra_cols) != 3) {
    UNPROTECT(1);
    Rf_error("extra tables must be a list of three (internal error)");
  }
  for (int t = 0; t < 3; t++) {
    SEXP cols = VECTOR_ELT(extra_cols, t);
    if (Rf_isNull(cols)) continue;
    status = add_columns(h, t + 2, cols, err);
    check_status(status, err, 1);
  }

  if (TYPEOF(roles) != VECSXP || XLENGTH(roles) != 3) {
    UNPROTECT(1);
    Rf_error("roles must be a list of three vectors (internal error)");
  }
  SEXP tables = VECTOR_ELT(roles, 0);
  SEXP codes = VECTOR_ELT(roles, 1);
  SEXP names = VECTOR_ELT(roles, 2);
  R_xlen_t k = XLENGTH(names);
  if (TYPEOF(tables) != INTSXP || TYPEOF(codes) != INTSXP || TYPEOF(names) != STRSXP ||
      XLENGTH(tables) != k || XLENGTH(codes) != k) {
    UNPROTECT(1);
    Rf_error("roles must hold integer tables, integer codes, and names (internal error)");
  }
  for (R_xlen_t i = 0; i < k; i++) {
    status = net_frame_set_role(h, INTEGER(tables)[i], INTEGER(codes)[i],
                                utf8(STRING_ELT(names, i)), err, ERR_LEN);
    check_status(status, err, 1);
  }
  status = net_frame_set_weight_type(h, Rf_asInteger(weight_type), err, ERR_LEN);
  check_status(status, err, 1);

  if (TYPEOF(declarations) != INTSXP || XLENGTH(declarations) != 7) {
    UNPROTECT(1);
    Rf_error("declarations must be seven integers (internal error)");
  }
  for (int d = 0; d < 7; d++) {
    int value = INTEGER(declarations)[d];
    /* Declaration 4 (add_nodes) is never undeclared. */
    if (d == 3 && value < 0) value = 0;
    status = net_frame_declare(h, d + 1, value, err, ERR_LEN);
    check_status(status, err, 1);
  }

  if (!Rf_isNull(geometry)) {
    if (TYPEOF(geometry) != VECSXP || XLENGTH(geometry) != 3) {
      UNPROTECT(1);
      Rf_error("geometry must be a list of three elements (internal error)");
    }
    SEXP crs = VECTOR_ELT(geometry, 2);
    const char *crs_text = Rf_isNull(crs) ? NULL : utf8(STRING_ELT(crs, 0));
    for (int t = 0; t < 2; t++) {
      int present = Rf_asLogical(VECTOR_ELT(geometry, t));
      if (present == NA_LOGICAL) present = 0;
      if (present) {
        status = net_frame_set_geometry(h, t, 1, crs_text, err, ERR_LEN);
        check_status(status, err, 1);
      }
    }
  }
  if (!Rf_isNull(coordinates)) {
    status = net_frame_set_coordinates(h, Rf_asInteger(coordinates), err, ERR_LEN);
    check_status(status, err, 1);
  }
  UNPROTECT(1);
  return ptr;
}

static SEXP named_list(int k, const char **names) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, k));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, k));
  for (int i = 0; i < k; i++) SET_STRING_ELT(nm, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

/* Runs every check. Returns a list: errors, notes, and the diagnostics as
 * parallel vectors (code, message, severity, table, count) plus a list of
 * the listed 0-based rows of each. */
SEXP C_frame_validate(SEXP ptr, SEXP threads) {
  FrameHandle *h = frame_from_ptr(ptr);
  char err[ERR_LEN] = {0};
  int64_t errors = 0, notes = 0;
  int status = net_frame_validate(h, Rf_asInteger(threads), &errors, &notes, err, ERR_LEN);
  check_status(status, err, 0);
  int64_t k = net_frame_diagnostic_count(h);
  if (k < 0) k = 0;
  static const char *names[] = {"errors", "notes", "code", "message", "severity",
                                "table",  "count", "rows"};
  SEXP out = PROTECT(named_list(8, names));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal((double)errors));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal((double)notes));
  SEXP code = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)k));
  SEXP message = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)k));
  SEXP severity = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)k));
  SEXP table = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)k));
  SEXP count = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)k));
  SEXP rows = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)k));
  for (int64_t i = 0; i < k; i++) {
    const char *c = NULL, *m = NULL;
    int32_t sev = 0, tab = -1;
    int64_t cnt = 0, nr = 0;
    status = net_frame_diagnostic(h, i, &c, &m, &sev, &tab, &cnt, &nr, err, ERR_LEN);
    check_status(status, err, 7);
    SET_STRING_ELT(code, (R_xlen_t)i, Rf_mkCharCE(c == NULL ? "" : c, CE_UTF8));
    SET_STRING_ELT(message, (R_xlen_t)i, Rf_mkCharCE(m == NULL ? "" : m, CE_UTF8));
    INTEGER(severity)[i] = sev;
    INTEGER(table)[i] = tab;
    REAL(count)[i] = (double)cnt;
    SEXP r = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nr));
    if (nr > 0) {
      int64_t *buf = (int64_t *)R_alloc((size_t)nr, sizeof(int64_t));
      status = net_frame_diagnostic_rows(h, i, buf, err, ERR_LEN);
      check_status(status, err, 8);
      for (int64_t j = 0; j < nr; j++) REAL(r)[j] = (double)buf[j];
    }
    SET_VECTOR_ELT(rows, (R_xlen_t)i, r);
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(out, 2, code);
  SET_VECTOR_ELT(out, 3, message);
  SET_VECTOR_ELT(out, 4, severity);
  SET_VECTOR_ELT(out, 5, table);
  SET_VECTOR_ELT(out, 6, count);
  SET_VECTOR_ELT(out, 7, rows);
  UNPROTECT(7);
  return out;
}

/* Validates if needed and builds the graph handle. The external pointer
 * exists, with its finalizer, before the engine builds the graph. */
SEXP C_frame_build(SEXP ptr, SEXP threads) {
  FrameHandle *h = frame_from_ptr(ptr);
  SEXP out = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(out, graph_release, TRUE);
  GraphHandle *g = NULL;
  char err[ERR_LEN] = {0};
  int status = net_frame_build(h, Rf_asInteger(threads), &g, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(out, g);
  UNPROTECT(1);
  return out;
}

static SEXP int_vector_from(const int64_t *v, int64_t n, int offset) {
  SEXP out = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)n));
  for (int64_t i = 0; i < n; i++) INTEGER(out)[i] = (int)v[i] + offset;
  UNPROTECT(1);
  return out;
}

/* The built graph's node count, the 0-based node index of every edge
 * row's endpoints (m rows), for each added node the edge row it first
 * appeared in (0-based) and whether as the target, and the 0-based edge
 * row of every graph edge (the observed ties; NULL when every row is an
 * edge of the graph). */
SEXP C_frame_index(SEXP ptr, SEXP m) {
  FrameHandle *h = frame_from_ptr(ptr);
  int64_t n = net_frame_node_count(h);
  int64_t added = net_frame_added_count(h);
  int64_t g = net_frame_graph_edge_count(h);
  if (n < 0 || added < 0 || g < 0) Rf_error("the frame has not been built (internal error)");
  int64_t k = (int64_t)Rf_asReal(m);
  char err[ERR_LEN] = {0};
  int64_t *s = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  int64_t *t = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  int status = net_frame_edge_index(h, s, t, err, ERR_LEN);
  check_status(status, err, 0);
  int64_t *rows = (int64_t *)R_alloc(added > 0 ? (size_t)added : 1, sizeof(int64_t));
  int32_t *is_target = (int32_t *)R_alloc(added > 0 ? (size_t)added : 1, sizeof(int32_t));
  if (added > 0) {
    status = net_frame_added(h, rows, is_target, err, ERR_LEN);
    check_status(status, err, 0);
  }
  int64_t *graph_rows = (int64_t *)R_alloc(g > 0 ? (size_t)g : 1, sizeof(int64_t));
  if (g != k && g > 0) {
    status = net_frame_edge_rows(h, graph_rows, err, ERR_LEN);
    check_status(status, err, 0);
  }
  static const char *names[] = {"n", "src", "tgt", "added_row", "added_target", "graph_rows"};
  SEXP out = PROTECT(named_list(6, names));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal((double)n));
  SET_VECTOR_ELT(out, 1, int_vector_from(s, k, 0));
  SET_VECTOR_ELT(out, 2, int_vector_from(t, k, 0));
  SET_VECTOR_ELT(out, 3, int_vector_from(rows, added, 0));
  SEXP target = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)added));
  for (int64_t i = 0; i < added; i++) LOGICAL(target)[i] = is_target[i] != 0;
  SET_VECTOR_ELT(out, 4, target);
  SET_VECTOR_ELT(out, 5, g == k ? R_NilValue : int_vector_from(graph_rows, g, 0));
  UNPROTECT(2);
  return out;
}

static SEXP label_vector(const FrameHandle *h, int which, int64_t k) {
  SEXP out = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)k));
  for (int64_t i = 0; i < k; i++) {
    const char *s = net_frame_summary_label(h, which, i);
    SET_STRING_ELT(out, (R_xlen_t)i, Rf_mkCharCE(s == NULL ? "" : s, CE_UTF8));
  }
  UNPROTECT(1);
  return out;
}

static SEXP count_vector(const FrameHandle *h, int which, int64_t k, char *err, int nprotect) {
  int64_t *buf = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  if (k > 0) {
    int status = net_frame_summary_counts(h, which, buf, err, ERR_LEN);
    check_status(status, err, nprotect);
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)k));
  for (int64_t i = 0; i < k; i++) REAL(out)[i] = (double)buf[i];
  UNPROTECT(1);
  return out;
}

/* The derived properties: the NET_FRAME_SUMMARY_LEN numbers, the mode,
 * relation, layer, and ego labels with their counts (nodes per mode,
 * edges per relation, intralayer edges per layer, nodes and edge rows per
 * ego), the weight column, and the CRS. */
SEXP C_frame_summary(SEXP ptr, SEXP threads) {
  FrameHandle *h = frame_from_ptr(ptr);
  char err[ERR_LEN] = {0};
  SEXP values = PROTECT(Rf_allocVector(REALSXP, NET_FRAME_SUMMARY_LEN));
  int status = net_frame_summary(h, Rf_asInteger(threads), REAL(values), err, ERR_LEN);
  check_status(status, err, 1);
  int64_t modes = (int64_t)REAL(values)[3];
  int64_t relations = (int64_t)REAL(values)[9];
  int64_t layers = R_FINITE(REAL(values)[26]) ? (int64_t)REAL(values)[26] : 0;
  int64_t egos = R_FINITE(REAL(values)[32]) ? (int64_t)REAL(values)[32] : 0;
  static const char *names[] = {"values",       "modes",       "mode_counts", "relations",
                                "relation_counts", "weight",   "crs",         "layers",
                                "layer_counts", "egos",        "ego_nodes",   "ego_edges"};
  SEXP out = PROTECT(named_list(12, names));
  SET_VECTOR_ELT(out, 0, values);
  SET_VECTOR_ELT(out, 1, label_vector(h, 0, modes));
  SET_VECTOR_ELT(out, 2, count_vector(h, 0, modes, err, 2));
  SET_VECTOR_ELT(out, 3, label_vector(h, 1, relations));
  SET_VECTOR_ELT(out, 4, count_vector(h, 1, relations, err, 2));
  const char *w = net_frame_summary_label(h, 2, 0);
  SET_VECTOR_ELT(out, 5, w == NULL ? R_NilValue : Rf_mkString(w));
  const char *crs = net_frame_summary_label(h, 3, 0);
  SET_VECTOR_ELT(out, 6, crs == NULL ? R_NilValue : Rf_mkString(crs));
  SET_VECTOR_ELT(out, 7, label_vector(h, 4, layers));
  SET_VECTOR_ELT(out, 8, count_vector(h, 2, layers, err, 2));
  SET_VECTOR_ELT(out, 9, label_vector(h, 5, egos));
  SET_VECTOR_ELT(out, 10, count_vector(h, 3, egos, err, 2));
  SET_VECTOR_ELT(out, 11, count_vector(h, 4, egos, err, 2));
  UNPROTECT(2);
  return out;
}

/* Integer codes of the modes (which 0, one per node row), relations
 * (which 1, one per edge row), node egos (which 2, one per node row), or
 * edge egos (which 3, one per edge row), 0-based by first appearance. */
SEXP C_frame_codes(SEXP ptr, SEXP which, SEXP n) {
  FrameHandle *h = frame_from_ptr(ptr);
  int64_t k = (int64_t)Rf_asReal(n);
  int64_t *buf = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  int status = net_frame_codes(h, Rf_asInteger(which), buf, err, ERR_LEN);
  check_status(status, err, 0);
  return int_vector_from(buf, k, 0);
}

/* Number of edges sharing each edge's key; m values. */
SEXP C_frame_multiplicity(SEXP ptr, SEXP m) {
  FrameHandle *h = frame_from_ptr(ptr);
  int64_t k = (int64_t)Rf_asReal(m);
  int64_t *buf = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  char err[ERR_LEN] = {0};
  int status = net_frame_multiplicity(h, buf, err, ERR_LEN);
  check_status(status, err, 0);
  return int_vector_from(buf, k, 0);
}

SEXP C_frame_is_simple(SEXP ptr) {
  FrameHandle *h = frame_from_ptr(ptr);
  int32_t simple = 0;
  char err[ERR_LEN] = {0};
  int status = net_frame_is_simple(h, &simple, err, ERR_LEN);
  check_status(status, err, 0);
  return Rf_ScalarLogical(simple != 0);
}

/* Groups parallel edges: the kept edge each row became (0-based, -1 for a
 * dropped loop), the combined weight at every row (NULL without one), and
 * the number of kept edges. */
SEXP C_frame_simplify(SEXP ptr, SEXP m, SEXP combine, SEXP drop_loops, SEXP weighted) {
  FrameHandle *h = frame_from_ptr(ptr);
  int64_t k = (int64_t)Rf_asReal(m);
  int64_t *group = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  int want_weight = Rf_asLogical(weighted);
  SEXP weight = PROTECT(want_weight ? Rf_allocVector(REALSXP, (R_xlen_t)k) : R_NilValue);
  double *w = want_weight ? REAL(weight) : NULL;
  for (int64_t i = 0; w != NULL && i < k; i++) w[i] = NA_REAL;
  int64_t kept = 0;
  char err[ERR_LEN] = {0};
  int status = net_frame_simplify(h, Rf_asInteger(combine), Rf_asLogical(drop_loops), group, w,
                                  &kept, err, ERR_LEN);
  check_status(status, err, 1);
  static const char *names[] = {"group", "weight", "kept"};
  SEXP out = PROTECT(named_list(3, names));
  SET_VECTOR_ELT(out, 0, int_vector_from(group, k, 0));
  SET_VECTOR_ELT(out, 1, weight);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal((double)kept));
  UNPROTECT(2);
  return out;
}

/* Layer codes (0-based indices into the layer labels) of every edge row's
 * endpoints and each row's class (0 intralayer, 1 coupling, 2 interlayer);
 * m values each. */
SEXP C_frame_edge_layers(SEXP ptr, SEXP m) {
  FrameHandle *h = frame_from_ptr(ptr);
  int64_t k = (int64_t)Rf_asReal(m);
  size_t alloc = k > 0 ? (size_t)k : 1;
  int64_t *sl = (int64_t *)R_alloc(alloc, sizeof(int64_t));
  int64_t *tl = (int64_t *)R_alloc(alloc, sizeof(int64_t));
  int32_t *cl = (int32_t *)R_alloc(alloc, sizeof(int32_t));
  char err[ERR_LEN] = {0};
  int status = net_frame_edge_layers(h, sl, tl, cl, err, ERR_LEN);
  check_status(status, err, 0);
  static const char *names[] = {"source_layer", "target_layer", "class"};
  SEXP out = PROTECT(named_list(3, names));
  SET_VECTOR_ELT(out, 0, int_vector_from(sl, k, 0));
  SET_VECTOR_ELT(out, 1, int_vector_from(tl, k, 0));
  SEXP cls = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)k));
  for (int64_t i = 0; i < k; i++) INTEGER(cls)[i] = (int)cl[i];
  SET_VECTOR_ELT(out, 2, cls);
  UNPROTECT(2);
  return out;
}

/* The subgraph of layer i (which 0) or ego i (which 1), 0-based: a graph
 * handle, the built graph's 0-based node index of each subgraph node, and
 * the 0-based edge row of each subgraph edge. The external pointer exists,
 * with its finalizer, before the engine builds the graph. */
SEXP C_frame_subgraph(SEXP ptr, SEXP which, SEXP i) {
  FrameHandle *h = frame_from_ptr(ptr);
  int w = Rf_asInteger(which);
  int64_t idx = (int64_t)Rf_asReal(i);
  int64_t nn = 0, ne = 0;
  char err[ERR_LEN] = {0};
  int status = net_frame_subgraph_size(h, w, idx, &nn, &ne, err, ERR_LEN);
  check_status(status, err, 0);
  int64_t *node_index = (int64_t *)R_alloc(nn > 0 ? (size_t)nn : 1, sizeof(int64_t));
  int64_t *edge_rows = (int64_t *)R_alloc(ne > 0 ? (size_t)ne : 1, sizeof(int64_t));
  SEXP g = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(g, graph_release, TRUE);
  GraphHandle *out = NULL;
  status = net_frame_subgraph(h, w, idx, &out, node_index, edge_rows, err, ERR_LEN);
  check_status(status, err, 1);
  R_SetExternalPtrAddr(g, out);
  static const char *names[] = {"ptr", "nodes", "edge_rows"};
  SEXP res = PROTECT(named_list(3, names));
  SET_VECTOR_ELT(res, 0, g);
  SET_VECTOR_ELT(res, 1, int_vector_from(node_index, nn, 0));
  SET_VECTOR_ELT(res, 2, int_vector_from(edge_rows, ne, 0));
  UNPROTECT(2);
  return res;
}

/* The resolved presence spells of the n nodes, sorted by node and start:
 * offsets (n + 1 values, 0-based into the spells), start and end
 * (infinite when unbounded), and the 0-based row of each spell. */
SEXP C_frame_node_spells(SEXP ptr, SEXP n) {
  FrameHandle *h = frame_from_ptr(ptr);
  int64_t nodes = (int64_t)Rf_asReal(n);
  int64_t k = net_frame_node_spell_count(h);
  if (k < 0) Rf_error("the frame has not been built (internal error)");
  int64_t *offsets = (int64_t *)R_alloc((size_t)nodes + 1, sizeof(int64_t));
  int64_t *rows = (int64_t *)R_alloc(k > 0 ? (size_t)k : 1, sizeof(int64_t));
  SEXP start = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)k));
  SEXP end = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)k));
  char err[ERR_LEN] = {0};
  int status = net_frame_node_spells(h, offsets, k > 0 ? REAL(start) : NULL,
                                     k > 0 ? REAL(end) : NULL, rows, err, ERR_LEN);
  check_status(status, err, 2);
  static const char *names[] = {"offsets", "start", "end", "rows"};
  SEXP out = PROTECT(named_list(4, names));
  SET_VECTOR_ELT(out, 0, int_vector_from(offsets, nodes + 1, 0));
  SET_VECTOR_ELT(out, 1, start);
  SET_VECTOR_ELT(out, 2, end);
  SET_VECTOR_ELT(out, 3, int_vector_from(rows, k, 0));
  UNPROTECT(3);
  return out;
}
