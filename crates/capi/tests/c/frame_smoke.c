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
  if (net_frame_set_role(f, 5, 0, "x", err, sizeof err) == 0) return 1;
  net_frame_free(f);

  /* A multilayer collection of ego networks with a node spells table, a
   * node changes table, and an unobserved pair: actors a, b, c (egos a
   * and c; b is a's alter); edges a-b in layer x, a-b in layer y, the
   * coupling a-a across x and y, and an unobserved a-b in layer x at a
   * later time. */
  f = net_frame_new(0);
  const char *ids3 = "abc";
  int64_t off3[] = {0, 1, 2, 3};
  if (net_frame_add_string(f, 0, "id", off3, ids3, NULL, 3, err, sizeof err) != 0)
    return fail("add id", err);
  const char *egos = "aac";
  if (net_frame_add_string(f, 0, "ego", off3, egos, NULL, 3, err, sizeof err) != 0)
    return fail("add ego", err);
  const char *from4 = "aaaa";
  const char *to4 = "bbab";
  int64_t off4[] = {0, 1, 2, 3, 4};
  if (net_frame_add_string(f, 1, "from", off4, from4, NULL, 4, err, sizeof err) != 0)
    return fail("add from", err);
  if (net_frame_add_string(f, 1, "to", off4, to4, NULL, 4, err, sizeof err) != 0)
    return fail("add to", err);
  const char *ls = "xyxx";
  const char *lt = "xyyx";
  if (net_frame_add_string(f, 1, "ls", off4, ls, NULL, 4, err, sizeof err) != 0)
    return fail("add ls", err);
  if (net_frame_add_string(f, 1, "lt", off4, lt, NULL, 4, err, sizeof err) != 0)
    return fail("add lt", err);
  double t4[] = {0.0, 0.0, 0.0, 5.0};
  if (net_frame_add_double(f, 1, "t", t4, 4, err, sizeof err) != 0) return fail("add t", err);
  int32_t na4[] = {0, 0, 0, 1};
  if (net_frame_add_logical(f, 1, "na", na4, NULL, 4, err, sizeof err) != 0)
    return fail("add na", err);
  const char *layer_ids = "yx";
  int64_t off2[] = {0, 1, 2};
  if (net_frame_add_string(f, 2, "layer", off2, layer_ids, NULL, 2, err, sizeof err) != 0)
    return fail("add layer", err);
  const char *spell_ids = "aa";
  double spell_s[] = {4.0, 0.0};
  double spell_e[] = {9.0, 2.0};
  if (net_frame_add_string(f, 3, "id", off2, spell_ids, NULL, 2, err, sizeof err) != 0)
    return fail("add spell id", err);
  if (net_frame_add_double(f, 3, "s", spell_s, 2, err, sizeof err) != 0)
    return fail("add spell s", err);
  if (net_frame_add_double(f, 3, "e", spell_e, 2, err, sizeof err) != 0)
    return fail("add spell e", err);
  const char *change_ids = "b";
  int64_t off1[] = {0, 1};
  double change_t[] = {1.0};
  if (net_frame_add_string(f, 4, "id", off1, change_ids, NULL, 1, err, sizeof err) != 0)
    return fail("add change id", err);
  if (net_frame_add_double(f, 4, "t", change_t, 1, err, sizeof err) != 0)
    return fail("add change t", err);
  if (net_frame_set_role(f, 0, 0, "id", err, sizeof err) != 0) return fail("role id", err);
  if (net_frame_set_role(f, 0, 7, "ego", err, sizeof err) != 0) return fail("role ego", err);
  if (net_frame_set_role(f, 1, 0, "from", err, sizeof err) != 0) return fail("role from", err);
  if (net_frame_set_role(f, 1, 1, "to", err, sizeof err) != 0) return fail("role to", err);
  if (net_frame_set_role(f, 1, 9, "ls", err, sizeof err) != 0) return fail("role ls", err);
  if (net_frame_set_role(f, 1, 10, "lt", err, sizeof err) != 0) return fail("role lt", err);
  if (net_frame_set_role(f, 1, 5, "t", err, sizeof err) != 0) return fail("role t", err);
  if (net_frame_set_role(f, 1, 8, "na", err, sizeof err) != 0) return fail("role na", err);
  if (net_frame_set_role(f, 2, 0, "layer", err, sizeof err) != 0) return fail("role layer", err);
  if (net_frame_set_role(f, 3, 0, "id", err, sizeof err) != 0) return fail("role spell id", err);
  if (net_frame_set_role(f, 3, 1, "s", err, sizeof err) != 0) return fail("role spell s", err);
  if (net_frame_set_role(f, 3, 2, "e", err, sizeof err) != 0) return fail("role spell e", err);
  if (net_frame_set_role(f, 4, 0, "id", err, sizeof err) != 0) return fail("role change id", err);
  if (net_frame_set_role(f, 4, 1, "t", err, sizeof err) != 0) return fail("role change t", err);
  if (net_frame_declare(f, 5, 1, err, sizeof err) != 0) return fail("declare egos", err);
  if (net_frame_declare(f, 6, 1, err, sizeof err) != 0) return fail("declare couplings", err);
  if (net_frame_validate(f, 1, &errors, &notes, err, sizeof err) != 0)
    return fail("validate", err);
  /* Only the unobserved-pair note (and c as an isolate). */
  if (errors != 0 || notes != 2) {
    fprintf(stderr, "expected no error and 2 notes, got %lld and %lld\n", (long long)errors,
            (long long)notes);
    return 1;
  }
  if (net_frame_build(f, 1, &g, err, sizeof err) != 0) return fail("build rich", err);
  if (net_graph_edge_count(g) != 3 || net_frame_graph_edge_count(f) != 3)
    return fail("observed edges", "");
  net_graph_free(g);
  int64_t erows[3];
  if (net_frame_edge_rows(f, erows, err, sizeof err) != 0 || erows[2] != 2)
    return fail("edge rows", err);
  int64_t msrc[1], mtgt[1], mrows[1];
  if (net_frame_missing_pairs(f, msrc, mtgt, mrows, err, sizeof err) != 0 || msrc[0] != 0 ||
      mtgt[0] != 1 || mrows[0] != 3)
    return fail("missing pairs", err);
  int64_t lsrc[4], ltgt[4];
  int32_t cls[4];
  if (net_frame_edge_layers(f, lsrc, ltgt, cls, err, sizeof err) != 0) return fail("layers", err);
  /* Layer labels follow the table: y is 0, x is 1. */
  if (lsrc[0] != 1 || ltgt[1] != 0 || cls[0] != 0 || cls[2] != 1) return fail("layer codes", "");
  int64_t egocodes[3];
  if (net_frame_codes(f, 2, egocodes, err, sizeof err) != 0 || egocodes[2] != 1)
    return fail("ego codes", err);
  if (net_frame_summary(f, 1, s, err, sizeof err) != 0) return fail("summary rich", err);
  /* The unobserved row still counts as an intralayer row of layer x. */
  if (s[1] != 4 || s[15] != 2 || s[25] != 1 || s[26] != 2 || s[27] != 1 || s[28] != 3 ||
      s[29] != 1 || s[30] != 0 || s[31] != 1 || s[32] != 2 || s[33] != 2 || s[34] != 1 ||
      s[35] != 2 || s[36] != 1)
    return fail("summary rich values", "");
  const char *lname = net_frame_summary_label(f, 4, 1);
  if (lname == NULL || strcmp(lname, "x") != 0) return fail("layer label", "");
  const char *ename = net_frame_summary_label(f, 5, 1);
  if (ename == NULL || strcmp(ename, "c") != 0) return fail("ego label", "");
  int64_t per_layer[2], per_ego[2];
  if (net_frame_summary_counts(f, 2, per_layer, err, sizeof err) != 0 || per_layer[0] != 1 ||
      per_layer[1] != 2)
    return fail("layer counts", err);
  if (net_frame_summary_counts(f, 3, per_ego, err, sizeof err) != 0 || per_ego[0] != 2 ||
      per_ego[1] != 1)
    return fail("ego node counts", err);
  int64_t sn = 0, se = 0;
  if (net_frame_subgraph_size(f, 0, 1, &sn, &se, err, sizeof err) != 0 || sn != 3 || se != 1)
    return fail("layer subgraph size", err);
  int64_t sub_nodes[3], sub_edges[1];
  if (net_frame_subgraph(f, 0, 1, &g, sub_nodes, sub_edges, err, sizeof err) != 0)
    return fail("layer subgraph", err);
  if (net_graph_edge_count(g) != 1 || sub_edges[0] != 0 || sub_nodes[2] != 2)
    return fail("layer subgraph values", "");
  net_graph_free(g);
  if (net_frame_subgraph_size(f, 1, 0, &sn, &se, err, sizeof err) != 0 || sn != 2 || se != 3)
    return fail("ego subgraph size", err);
  if (net_frame_subgraph(f, 1, 0, &g, NULL, NULL, err, sizeof err) != 0)
    return fail("ego subgraph", err);
  if (net_graph_node_count(g) != 2) return fail("ego subgraph nodes", "");
  net_graph_free(g);
  if (net_frame_subgraph_size(f, 1, 5, &sn, &se, err, sizeof err) == 0)
    return fail("ego out of range", "");
  if (net_frame_node_spell_count(f) != 2) return fail("spell count", "");
  int64_t soff[4], srows[2];
  double sstart[2], send[2];
  if (net_frame_node_spells(f, soff, sstart, send, srows, err, sizeof err) != 0)
    return fail("node spells", err);
  if (soff[1] != 2 || soff[3] != 2 || sstart[0] != 0.0 || send[1] != 9.0 || srows[0] != 1)
    return fail("node spell values", "");
  net_frame_free(f);
  puts("frame ABI smoke test passed");
  return 0;
}
