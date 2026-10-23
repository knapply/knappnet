/* Compiles a C program against netcapi.h and the static library, checking
 * that the header matches the exported ABI. Run by
 * tools/check_frame_abi.sh. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "netcapi.h"

static int fail(const char *what, const char *err) {
  fprintf(stderr, "%s: %s\n", what, err);
  return 1;
}

int main(void) {
  char err[512];
  /* Nodes a, b, c with modes; edges a-b, b-c, b-c (parallel), c-d (d is
   * not in the node table). */
  FrameHandle *f = net_frame_new(0);
  const char *ids = "abc";
  int64_t id_off[] = {0, 1, 2, 3};
  if (net_frame_add_string(f, 0, "id", id_off, ids, NULL, 3, err, sizeof err) != 0)
    return fail("add id", err);
  int64_t modes[] = {1, 2, 1};
  if (net_frame_add_int(f, 0, "mode", modes, NULL, 3, err, sizeof err) != 0)
    return fail("add mode", err);
  const char *from = "abbc";
  const char *to = "bccd";
  int64_t e_off[] = {0, 1, 2, 3, 4};
  if (net_frame_add_string(f, 1, "from", e_off, from, NULL, 4, err, sizeof err) != 0)
    return fail("add from", err);
  if (net_frame_add_string(f, 1, "to", e_off, to, NULL, 4, err, sizeof err) != 0)
    return fail("add to", err);
  double w[] = {1.0, 2.0, 3.0, 4.0};
  if (net_frame_add_double(f, 1, "w", w, 4, err, sizeof err) != 0) return fail("add w", err);
  if (net_frame_set_role(f, 0, 0, "id", err, sizeof err) != 0) return fail("role id", err);
  if (net_frame_set_role(f, 0, 1, "mode", err, sizeof err) != 0) return fail("role mode", err);
  if (net_frame_set_role(f, 1, 0, "from", err, sizeof err) != 0) return fail("role from", err);
  if (net_frame_set_role(f, 1, 1, "to", err, sizeof err) != 0) return fail("role to", err);
  if (net_frame_set_role(f, 1, 2, "w", err, sizeof err) != 0) return fail("role w", err);
  if (net_frame_set_weight_type(f, 1, err, sizeof err) != 0) return fail("weight type", err);

  /* Without add_nodes, d is unknown: one error, and no note yet because
   * the structural checks need resolved endpoints. */
  int64_t errors = 0, notes = 0;
  if (net_frame_validate(f, 1, &errors, &notes, err, sizeof err) != 0)
    return fail("validate", err);
  if (errors != 1 || notes != 0) {
    fprintf(stderr, "expected 1 error and no note, got %lld and %lld\n", (long long)errors,
            (long long)notes);
    return 1;
  }
  const char *code = NULL, *message = NULL;
  int32_t severity = -9, table = -9;
  int64_t count = 0, n_rows = 0;
  if (net_frame_diagnostic(f, 0, &code, &message, &severity, &table, &count, &n_rows, err,
                           sizeof err) != 0)
    return fail("diagnostic", err);
  if (strcmp(code, "E_UNKNOWN_NODE") != 0 || severity != 0 || table != 1 || count != 1 ||
      n_rows != 1 || message == NULL)
    return fail("diagnostic fields", code);
  int64_t rows[1];
  if (net_frame_diagnostic_rows(f, 0, rows, err, sizeof err) != 0 || rows[0] != 3)
    return fail("diagnostic rows", err);
  GraphHandle *g = NULL;
  if (net_frame_build(f, 1, &g, err, sizeof err) == 0) return fail("build should fail", "");
  if (strstr(err, "E_") != NULL) return fail("error text should be prose", err);

  /* Declaring add_nodes makes it valid once the mode role is cleared (an
   * added node has no mode); the build then has 4 nodes. */
  if (net_frame_declare(f, 4, 1, err, sizeof err) != 0) return fail("declare", err);
  if (net_frame_diagnostic_count(f) != -1) return fail("stale diagnostics", "");
  if (net_frame_build(f, 1, &g, err, sizeof err) == 0) return fail("build without mode", "");
  if (net_frame_set_role(f, 0, 1, NULL, err, sizeof err) != 0) return fail("clear mode", err);
  if (net_frame_build(f, 1, &g, err, sizeof err) != 0) return fail("build", err);
  if (net_frame_validate(f, 1, &errors, &notes, err, sizeof err) != 0 || errors != 0 || notes != 2)
    return fail("notes after build", err);
  if (net_graph_node_count(g) != 4 || net_frame_node_count(f) != 4) return fail("node count", "");
  if (net_frame_added_count(f) != 1) return fail("added count", "");
  int64_t added_row[1];
  int32_t added_target[1];
  if (net_frame_added(f, added_row, added_target, err, sizeof err) != 0 || added_row[0] != 3 ||
      added_target[0] != 1)
    return fail("added", err);
  int64_t src[4], tgt[4];
  if (net_frame_edge_index(f, src, tgt, err, sizeof err) != 0 || src[3] != 2 || tgt[3] != 3)
    return fail("edge index", err);
  int64_t mult[4];
  if (net_frame_multiplicity(f, mult, err, sizeof err) != 0 || mult[1] != 2 || mult[0] != 1)
    return fail("multiplicity", err);
  int32_t simple = 1;
  if (net_frame_is_simple(f, &simple, err, sizeof err) != 0 || simple != 0)
    return fail("is simple", err);
  int64_t group[4], kept = 0;
  double cw[4] = {0, 0, 0, 0};
  if (net_frame_simplify(f, 0, 0, group, cw, &kept, err, sizeof err) != 0)
    return fail("simplify", err);
  if (kept != 3 || group[2] != 1 || cw[1] != 5.0 || cw[2] != 5.0) return fail("simplify values", "");
  net_graph_free(g);
  double s[NET_FRAME_SUMMARY_LEN];
  if (net_frame_summary(f, 1, s, err, sizeof err) != 0) return fail("summary", err);
  if (s[0] != 4 || s[1] != 4 || s[2] != 0 || s[5] != 1 || s[6] != 1 || s[7] != 1.0 ||
      s[8] != 4.0 || s[20] != 1 || s[22] != 1 || fabs(s[23] - 4.0 / 6.0) > 1e-12)
    return fail("summary values", "");
  const char *wname = net_frame_summary_label(f, 2, 0);
  if (wname == NULL || strcmp(wname, "w") != 0) return fail("weight label", "");
  if (net_frame_summary_label(f, 0, 0) != NULL) return fail("mode label without modes", "");
  net_frame_free(f);
  net_frame_free(NULL);

  /* Errors come back as messages, not crashes. */
  if (net_frame_add_int(NULL, 0, "x", NULL, NULL, 0, err, sizeof err) == 0) return 1;
  f = net_frame_new(1);
  if (net_frame_set_role(f, 1, 0, "nope", err, sizeof err) == 0) return 1;
  if (net_frame_declare(f, 9, 1, err, sizeof err) == 0) return 1;
  net_frame_free(f);
  puts("frame ABI smoke test passed");
  return 0;
}
