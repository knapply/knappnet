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

/* threads == 0 uses all available cores. */
int32_t net_betweenness(const GraphHandle *h, int32_t weighted,
                        int32_t normalized, int32_t threads, double *out,
                        char *err, size_t err_len);

/* Writes 16 counts in MAN order: 003 012 102 021D 021U 021C 111D 111U
 * 030T 030C 201 120D 120U 120C 210 300. */
int32_t net_triad_census(const GraphHandle *h, int32_t threads, double *out,
                         char *err, size_t err_len);

#ifdef __cplusplus
}
#endif

#endif
