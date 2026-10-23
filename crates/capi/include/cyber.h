/* ---- Key players, group centrality, flow betweenness ----------------
 * Declarations to be merged into netcapi.h.
 *
 * Node sets are 0-based indices without repeats.
 * measure (key players): 0 = fragmentation, 1 = distance-weighted
 *   fragmentation, 2 = m-reach, 3 = distance-weighted reach.
 * measure (group centrality): 0 = degree, 1 = closeness, 2 = betweenness.
 * mode: 1 = out, 2 = in, 3 = all. weights_as: 0 = distance, 1 = strength.
 * method: 0 = greedy swaps with random restarts, 1 = exact enumeration.
 * scale (flow betweenness): 0 = raw, 1 = normalized, 2 = fractional. */

/* Criterion value for the k nodes in set; m is the m-reach distance limit. */
int32_t net_kpp_score(const GraphHandle *h, const int64_t *set, int64_t k,
                      int32_t measure, double m, int32_t mode, int32_t weighted,
                      int32_t weights_as, double *out, char *err,
                      size_t err_len);

/* Finds k key players. Greedy restart r starts from stream r of seed. set
 * receives k indices in increasing order; score one double; evaluations
 * the number of sets evaluated. */
int32_t net_key_players(const GraphHandle *h, int64_t k, int32_t measure,
                        double m, int32_t mode, int32_t weighted,
                        int32_t weights_as, int32_t method, int64_t restarts,
                        uint64_t seed, int32_t threads, int64_t *set,
                        double *score, double *evaluations, char *err,
                        size_t err_len);

/* Group degree, closeness, or betweenness of the k nodes in group. */
int32_t net_group_centrality(const GraphHandle *h, const int64_t *group,
                             int64_t k, int32_t measure, int32_t mode,
                             int32_t weighted, int32_t weights_as,
                             int32_t normalized, int32_t threads, double *out,
                             char *err, size_t err_len);

/* Flow betweenness (Freeman, Borgatti, and White 1991) of every node into
 * out (n doubles). capacity != 0 uses edge weights as capacities. */
int32_t net_flow_betweenness(const GraphHandle *h, int32_t capacity,
                             int32_t scale, int32_t threads, double *out,
                             char *err, size_t err_len);
