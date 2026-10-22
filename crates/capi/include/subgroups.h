/* ---- Cohesive subgroups ----------------------------------------------
 * Cores, cliques, and communities. Declarations to be merged into
 * netcapi.h. */

/* Coreness (k-core number) of every node, written to out[0..n). mode:
 * 1 = out, 2 = in, 3 = all; ignored for undirected graphs. */
int32_t net_coreness(const GraphHandle *h, int32_t mode, double *out,
                     char *err, size_t err_len);

/* A collection of node groups, read in two steps: query the number of
 * groups and members, allocate, then copy. Members of group i are
 * members[offsets[i] .. offsets[i + 1]), 0-based and increasing. */
typedef struct GroupsHandle GroupsHandle;

int64_t net_groups_count(const GroupsHandle *h);
int64_t net_groups_total(const GroupsHandle *h);
/* offsets must hold count + 1 values and members total values. */
int32_t net_groups_copy(const GroupsHandle *h, int64_t *offsets,
                        int64_t *members, char *err, size_t err_len);
void net_groups_free(GroupsHandle *h);

/* Maximal cliques of the undirected simple graph underlying the network
 * (directions, loops, and edge multiplicity ignored) with at least
 * min_size members, in lexicographic order. Release the result with
 * net_groups_free. threads == 0 uses all available cores. */
int32_t net_max_cliques(const GraphHandle *h, int64_t min_size,
                        int32_t threads, GroupsHandle **out, char *err,
                        size_t err_len);

/* Size of the largest clique (0 for a graph with no nodes). */
int32_t net_clique_number(const GraphHandle *h, int32_t threads,
                          int64_t *out, char *err, size_t err_len);

/* Modularity of a partition: membership holds n non-negative labels.
 * Weights (weighted != 0) must be non-negative. Writes NaN when the
 * network has no edges or total weight zero. */
int32_t net_modularity(const GraphHandle *h, const int64_t *membership,
                       int64_t n, double resolution, int32_t weighted,
                       double *out, char *err, size_t err_len);

/* Louvain community detection on an undirected graph. Writes 0-based
 * group labels, numbered by first appearance, to membership[0..n) and the
 * partition's modularity to *q. The same seed gives the same result. */
int32_t net_louvain(const GraphHandle *h, double resolution, int32_t weighted,
                    uint64_t seed, int64_t *membership, double *q, char *err,
                    size_t err_len);
