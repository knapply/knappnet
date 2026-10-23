/*
 * Temporal networks: contacts with start and end times, time-respecting
 * paths, reachability, temporal closeness, slicing into windows,
 * aggregation, and inter-event statistics. To be merged into netcapi.h
 * (after the GraphHandle declarations; needs <stddef.h> and <stdint.h>).
 *
 * Codes:
 *   criterion: 0 = earliest arrival, 1 = latest departure, 2 = fastest,
 *              3 = shortest (fewest contacts).
 *   direction: 1 = walks from the roots, 2 = walks into the roots.
 *   traversal: 0 = the constant traversal_time for every contact,
 *              1 = each contact's duration (end - start).
 *   scope:     0 = nodes, 1 = node pairs; mode (nodes of directed
 *              networks): 1 = out, 2 = in, 3 = all.
 * Window bounds start and end may be NaN for the earliest start time and
 * the latest arrival of any contact; max_wait may be infinite. strict != 0
 * requires waiting times strictly greater than min_wait.
 */

typedef struct TemporalHandle TemporalHandle;

/* Numbered columns of doubles, for results whose length is not known in
 * advance. */
typedef struct TemporalResult TemporalResult;

/* end may be NULL for instantaneous contacts; otherwise end[i] >= start[i].
 * Times must be finite. Release the handle with net_temporal_free. */
int32_t net_temporal_new(int64_t n, const int64_t *src, const int64_t *tgt,
                         const double *start, const double *end, int64_t m,
                         int32_t directed, TemporalHandle **out, char *err,
                         size_t err_len);

void net_temporal_free(TemporalHandle *h);

int64_t net_temporal_node_count(const TemporalHandle *h);

/* out[r * n + j]: walks from roots[r] to j (direction 1) or from j to
 * roots[r] (direction 2). Unreachable: +Inf (-Inf for latest departure).
 * out must hold nroots * n doubles. */
int32_t net_temporal_paths(const TemporalHandle *h, int32_t criterion,
                           int32_t direction, const int64_t *roots,
                           int64_t nroots, double start, double end,
                           double min_wait, double max_wait, int32_t strict,
                           int32_t traversal, double traversal_time,
                           int32_t threads, double *out, char *err,
                           size_t err_len);

/* Number of other nodes each node reaches (direction 1) or is reached by
 * (direction 2). out must hold n doubles. */
int32_t net_temporal_reach(const TemporalHandle *h, int32_t direction,
                           double start, double end, double min_wait,
                           double max_wait, int32_t strict, int32_t traversal,
                           double traversal_time, int32_t threads, double *out,
                           char *err, size_t err_len);

/* Pan-Saramaki temporal closeness (n values) and, unless distances is NULL,
 * average temporal distances (n * n values, row = origin). */
int32_t net_temporal_closeness(const TemporalHandle *h, double start,
                               double end, double min_wait, double max_wait,
                               int32_t strict, int32_t traversal,
                               double traversal_time, int32_t threads,
                               double *closeness, double *distances, char *err,
                               size_t err_len);

/* Contacts active in [from[k], to[k]). Columns: 0 = offsets (k + 1),
 * 1 = 0-based contact indices. */
int32_t net_temporal_windows(const TemporalHandle *h, const double *from,
                             const double *to, int64_t k, TemporalResult **out,
                             char *err, size_t err_len);

/* Contacts active at time t. Column 0: 0-based contact indices. */
int32_t net_temporal_active_at(const TemporalHandle *h, double t,
                               TemporalResult **out, char *err,
                               size_t err_len);

/* Contacts grouped by node pair (within [from, to) when windowed != 0,
 * durations clipped). Columns: 0 = from, 1 = to, 2 = count, 3 = duration,
 * 4 = first start, 5 = last end, 6 = pair of each contact (-1 outside). */
int32_t net_temporal_aggregate(const TemporalHandle *h, int32_t windowed,
                               double from, double to, TemporalResult **out,
                               char *err, size_t err_len);

/* Columns: 0 = group of each inter-event time, 1 = inter-event times,
 * 2 = events per group. */
int32_t net_temporal_inter_event(const TemporalHandle *h, int32_t scope,
                                 int32_t mode, TemporalResult **out, char *err,
                                 size_t err_len);

/* Columns: 0 = events, 1 = burstiness B, 2 = memory M (Goh-Barabasi). */
int32_t net_temporal_burstiness(const TemporalHandle *h, int32_t scope,
                                int32_t mode, TemporalResult **out, char *err,
                                size_t err_len);

/* Length of column col, or -1 if there is no such column. */
int64_t net_temporal_result_length(const TemporalResult *h, int32_t col);

int32_t net_temporal_result_copy(const TemporalResult *h, int32_t col,
                                 double *out, char *err, size_t err_len);

void net_temporal_result_free(TemporalResult *h);
