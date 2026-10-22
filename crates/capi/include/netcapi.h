/*
 * C interface to the network analysis engine.
 *
 * All functions return 0 on success. On failure they return non-zero and
 * write a NUL-terminated message into `err` (at most `err_len` bytes).
 * Node indices are 0-based. Output buffers are allocated by the caller.
 */
#ifndef NETCAPI_H
#define NETCAPI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GraphHandle GraphHandle;

int32_t net_graph_new(int64_t n, const int64_t *src, const int64_t *tgt,
                      const double *weights, int64_t m, int32_t directed,
                      GraphHandle **out, char *err, size_t err_len);

void net_graph_free(GraphHandle *h);

int64_t net_graph_node_count(const GraphHandle *h);

/* mode: 1 = out, 2 = in, 3 = all. weighted != 0 returns strength. */
int32_t net_degree(const GraphHandle *h, int32_t mode, int32_t loops,
                   int32_t weighted, double *out, char *err, size_t err_len);

/* weights_as: 0 = distance, 1 = strength (lengths 1/w); only used when
 * weighted != 0. threads == 0 uses all available cores. */
int32_t net_betweenness(const GraphHandle *h, int32_t weighted,
                        int32_t weights_as, int32_t normalized,
                        int32_t threads, double *out, char *err,
                        size_t err_len);

/* Writes 16 counts in MAN order: 003 012 102 021D 021U 021C 111D 111U
 * 030T 030C 201 120D 120U 120C 210 300. */
int32_t net_triad_census(const GraphHandle *h, int32_t threads, double *out,
                         char *err, size_t err_len);

/* ---- Files ----------------------------------------------------------
 * A network handle holds three tables: 0 = graph attributes (one row),
 * 1 = nodes (first column "id"), 2 = edges (first columns "from", "to").
 * Column kinds: 0 = logical, 1 = integer, 2 = double, 3 = text.
 * Formats: graphml, gexf, gml, pajek, dl, json, dynetml (ORA), dot (write
 * only). Reading a .xml file without a format chooses GraphML or DyNetML
 * from the root element; writing one produces GraphML. */

typedef struct NetworkHandle NetworkHandle;

/* format may be NULL to infer from the extension. directed is -1 (not
 * specified), 0, or 1; it only matters for formats without direction. */
int32_t net_read(const char *path, const char *format, int32_t directed,
                 NetworkHandle **out, char *err, size_t err_len);

NetworkHandle *net_network_new(int32_t directed);
void net_network_free(NetworkHandle *h);
int32_t net_network_directed(const NetworkHandle *h);

int64_t net_table_nrow(const NetworkHandle *h, int32_t which);
int32_t net_table_ncol(const NetworkHandle *h, int32_t which);
const char *net_table_name(const NetworkHandle *h, int32_t which, int32_t col);
int32_t net_table_kind(const NetworkHandle *h, int32_t which, int32_t col);

/* Logical, integer, and double columns as doubles; na[i] = 1 if missing. */
int32_t net_table_numeric(const NetworkHandle *h, int32_t which, int32_t col,
                          double *out, int32_t *na, char *err, size_t err_len);

/* NULL if missing or not a text column. Valid while the handle lives. */
const char *net_table_string(NetworkHandle *h, int32_t which, int32_t col,
                             int64_t row);

/* kind: 0 logical, 1 integer, 2 double. na may be NULL. */
int32_t net_table_add_numeric(NetworkHandle *h, int32_t which, const char *name,
                              int32_t kind, const double *vals,
                              const int32_t *na, int64_t n, char *err,
                              size_t err_len);

/* NULL entries are missing values. */
int32_t net_table_add_string(NetworkHandle *h, int32_t which, const char *name,
                             const char *const *vals, int64_t n, char *err,
                             size_t err_len);

/* format and weight may be NULL. */
int32_t net_write(const NetworkHandle *h, const char *path, const char *format,
                  const char *weight, char *err, size_t err_len);

#ifdef __cplusplus
}
#endif

#endif
