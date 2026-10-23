/*
 * Data model and validation. To be merged into netcapi.h as its own
 * section; until then include it after netcapi.h.
 */
#ifndef NETCAPI_FRAME_H
#define NETCAPI_FRAME_H

#include "netcapi.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Data model and validation ------------------------------------------
 * A frame is a node table and an edge table of typed columns (integer,
 * double, text, logical; each with missing values), the roles that give
 * columns structural meaning, and the declared properties of the
 * network. It is filled in steps (columns, roles, weight type,
 * declarations, geometry), validated, and built into a graph. Every
 * change after validation discards the diagnostics and the build.
 *
 * Codes: table 0 = nodes, 1 = edges. Node roles 0 = id, 1 = mode,
 * 2 = node_start, 3 = node_end, 4 = x, 5 = y, 6 = node_missing; edge
 * roles 0 = source, 1 = target, 2 = weight, 3 = relation, 4 = sign,
 * 5 = start, 6 = end, 7 = wave, 8 = edge_missing. weight_type 0 =
 * distance, 1 = strength, -1 = none. Declarations 0 = directed,
 * 1 = loops, 2 = multiple, 3 = between_modes_only, 4 = add_nodes, with
 * value 0 = false, 1 = true, -1 = not declared (1 to 3 only). Coordinate
 * kinds 0 = unknown, 1 = longitude and latitude, 2 = planar. Combine
 * 0 = sum, 1 = mean, 2 = min, 3 = max, 4 = first, 5 = count.
 *
 * Buffers: a column of n values is read from the caller's buffers and
 * copied; `missing` is n bytes (non-zero = missing) or NULL. Logical
 * values are int32_t (non-zero = true). Text is UTF-8 in one byte buffer
 * with n + 1 offsets (value i is bytes offsets[i]..offsets[i + 1]).
 * Strings returned by accessors are valid while the handle lives and
 * the frame is unchanged. */

/* A frame; release with net_frame_free. */
typedef struct FrameHandle FrameHandle;

/* Number of doubles net_frame_summary writes. */
#define NET_FRAME_SUMMARY_LEN 26

/* An empty frame; directed != 0 for a directed network. */
FrameHandle *net_frame_new(int32_t directed);

/* Releases a frame (NULL is allowed). */
void net_frame_free(FrameHandle *h);

/* Gives the frame a node table of n rows (only needed when no node
 * column is added; node columns must then have n rows). */
int32_t net_frame_node_rows(FrameHandle *h, int64_t n, char *err,
                            size_t err_len);

/* Adds an integer column of n values; missing may be NULL. */
int32_t net_frame_add_int(FrameHandle *h, int32_t table, const char *name,
                          const int64_t *values, const uint8_t *missing,
                          int64_t n, char *err, size_t err_len);

/* Adds a double column of n values; NaN is missing. */
int32_t net_frame_add_double(FrameHandle *h, int32_t table, const char *name,
                             const double *values, int64_t n, char *err,
                             size_t err_len);

/* Adds a logical column of n values (non-zero = true); missing may be
 * NULL. */
int32_t net_frame_add_logical(FrameHandle *h, int32_t table, const char *name,
                              const int32_t *values, const uint8_t *missing,
                              int64_t n, char *err, size_t err_len);

/* Adds a text column of n values from n + 1 offsets into bytes (UTF-8,
 * no terminators); missing may be NULL. */
int32_t net_frame_add_string(FrameHandle *h, int32_t table, const char *name,
                             const int64_t *offsets, const char *bytes,
                             const uint8_t *missing, int64_t n, char *err,
                             size_t err_len);

/* Assigns role of table to the column called name; NULL clears the role. */
int32_t net_frame_set_role(FrameHandle *h, int32_t table, int32_t role,
                           const char *name, char *err, size_t err_len);

/* Declares the meaning of the weights (0 distance, 1 strength, -1 none). */
int32_t net_frame_set_weight_type(FrameHandle *h, int32_t weight_type,
                                  char *err, size_t err_len);

/* Sets a declaration (see the codes above). */
int32_t net_frame_declare(FrameHandle *h, int32_t which, int32_t value,
                          char *err, size_t err_len);

/* Records whether table carries geometry and its coordinate reference
 * system (crs may be NULL). The two tables' systems must agree. */
int32_t net_frame_set_geometry(FrameHandle *h, int32_t table, int32_t present,
                               const char *crs, char *err, size_t err_len);

/* What the representative coordinates (roles x and y) are. */
int32_t net_frame_set_coordinates(FrameHandle *h, int32_t kind, char *err,
                                  size_t err_len);

/* Runs every check and writes the number of errors and notes (either
 * pointer may be NULL). Returns 0 even for an invalid frame; non-zero
 * only when the call itself fails. */
int32_t net_frame_validate(FrameHandle *h, int32_t threads, int64_t *errors,
                           int64_t *notes, char *err, size_t err_len);

/* Number of diagnostics, or -1 for NULL or before validation. */
int64_t net_frame_diagnostic_count(const FrameHandle *h);

/* Diagnostic i: code (such as "E_UNKNOWN_NODE"), message, severity
 * (0 error, 1 note), table (0 nodes, 1 edges, -1 network), the number of
 * affected rows in all, and how many net_frame_diagnostic_rows lists (at
 * most 50). Any output may be NULL. */
int32_t net_frame_diagnostic(const FrameHandle *h, int64_t i, const char **code,
                             const char **message, int32_t *severity,
                             int32_t *table, int64_t *count, int64_t *n_rows,
                             char *err, size_t err_len);

/* The listed 0-based rows of diagnostic i; rows holds n_rows values. */
int32_t net_frame_diagnostic_rows(const FrameHandle *h, int64_t i, int64_t *rows,
                                  char *err, size_t err_len);

/* Validates if needed and builds the graph (nodes: the node table's rows,
 * then nodes added from the edges in order of first appearance; edge e
 * is row e). Fails with the errors as the message when the frame is
 * invalid. Release *out with net_graph_free. */
int32_t net_frame_build(FrameHandle *h, int32_t threads, GraphHandle **out,
                        char *err, size_t err_len);

/* Nodes of the built graph, or -1 before a build. */
int64_t net_frame_node_count(const FrameHandle *h);

/* Nodes added from the edges, or -1 before a build. */
int64_t net_frame_added_count(const FrameHandle *h);

/* For each added node, in node order: the edge row its identifier first
 * appeared in and whether that was the target (1) or source (0). Both
 * buffers hold net_frame_added_count values. */
int32_t net_frame_added(const FrameHandle *h, int64_t *rows, int32_t *is_target,
                        char *err, size_t err_len);

/* The 0-based node index of every edge's endpoints (m values each). */
int32_t net_frame_edge_index(const FrameHandle *h, int64_t *src, int64_t *tgt,
                             char *err, size_t err_len);

/* Integer codes of the modes (which 0; one per node table row) or the
 * relations (which 1; one per edge), numbered by first appearance of
 * each label. Fails when the role is not set. */
int32_t net_frame_codes(const FrameHandle *h, int32_t which, int64_t *out,
                        char *err, size_t err_len);

/* Number of edges sharing each edge's key (endpoints, relation, spell);
 * m values. */
int32_t net_frame_multiplicity(const FrameHandle *h, int64_t *out, char *err,
                               size_t err_len);

/* 1 when the network has neither self-loops nor parallel edges. */
int32_t net_frame_is_simple(const FrameHandle *h, int32_t *out, char *err,
                            size_t err_len);

/* Groups parallel edges, dropping self-loops when drop_loops != 0.
 * group[e] (m values) is the kept edge that row e became (numbered by
 * first appearance; the kept edge is its group's first row) or -1 for a
 * removed loop; weight (m values or NULL) receives at every row the
 * combined weight of its group (only count writes it without a weight
 * role); *kept is the number of kept edges. */
int32_t net_frame_simplify(const FrameHandle *h, int32_t combine,
                           int32_t drop_loops, int64_t *group, double *weight,
                           int64_t *kept, char *err, size_t err_len);

/* Derived properties as NET_FRAME_SUMMARY_LEN doubles: 0 nodes, 1 edges,
 * 2 directed, 3 modes, 4 between_modes_only (-1 undeclared), 5 weighted,
 * 6 weight type (0 distance, 1 strength, -1 none), 7 smallest weight,
 * 8 largest weight, 9 relations, 10 signed, 11 negative edges, 12 time
 * kind (0 static, 1 contacts, 2 edge spells, 3 waves), 13 earliest
 * start, 14 latest end, 15 node spells, 16 node geometry, 17 edge
 * geometry, 18 coordinate kind (-1 none), 19 loops, 20 parallel edges,
 * 21 isolates, 22 weak components, 23 density, 24 unobserved nodes,
 * 25 unobserved edges. Undefined values are NaN. Builds if needed. */
int32_t net_frame_summary(FrameHandle *h, int32_t threads, double *out,
                          char *err, size_t err_len);

/* Label i of the modes (which 0) or relations (which 1), the weight
 * column (which 2, i 0), or the coordinate reference system (which 3,
 * i 0); NULL when out of range or before net_frame_summary. */
const char *net_frame_summary_label(const FrameHandle *h, int32_t which,
                                    int64_t i);

/* Nodes per mode (which 0) or edges per relation (which 1); out holds as
 * many values as the summary reports. */
int32_t net_frame_summary_counts(const FrameHandle *h, int32_t which,
                                 int64_t *out, char *err, size_t err_len);

#ifdef __cplusplus
}
#endif

#endif /* NETCAPI_FRAME_H */
