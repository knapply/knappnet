/*
 * Cohesive subsets, community partitions, topography (core-periphery,
 * compactness, small worlds), and random graphs. To be merged into
 * netcapi.h.
 */

/* A collection of node sets (bicomponents, n-cliques, k-plexes), read in
 * two steps like GroupsHandle: query the number of sets and members,
 * allocate, then copy. Members of set i are members[offsets[i] ..
 * offsets[i + 1]), 0-based and increasing; sets are in lexicographic
 * order. Release with net_node_sets_free. */
typedef struct NodeSetsHandle NodeSetsHandle;

/* Number of sets, or -1 for a NULL handle. */
int64_t net_node_sets_count(const NodeSetsHandle *h);
/* Total number of members over all sets, or -1 for a NULL handle. */
int64_t net_node_sets_total(const NodeSetsHandle *h);
/* offsets must hold count + 1 values and members total values. */
int32_t net_node_sets_copy(const NodeSetsHandle *h, int64_t *offsets,
                           int64_t *members, char *err, size_t err_len);
void net_node_sets_free(NodeSetsHandle *h);

/* Bicomponents (blocks: maximal subgraphs without a cut point) of the
 * undirected simple graph underlying the network (directions, loops,
 * and edge multiplicity ignored). A bridge is a block of two nodes;
 * isolated nodes belong to no block. */
int32_t net_bicomponents(const GraphHandle *h, NodeSetsHandle **out,
                         char *err, size_t err_len);

/* Maximal n-cliques (clans == 0) or n-clans (clans != 0) of the
 * undirected simple graph with at least min_size members: sets in which
 * every pair is within distance `order` in the whole network (n-clique)
 * or within the set itself (n-clan). threads == 0 uses all cores. */
int32_t net_n_cliques(const GraphHandle *h, int64_t order, int64_t min_size,
                      int32_t clans, int32_t threads, NodeSetsHandle **out,
                      char *err, size_t err_len);

/* Maximal k-plexes (every member adjacent to at least size - k others)
 * with at least min_size members. The search stops with an error after
 * max_steps steps (at least 1), a guard against exponential blow-up. */
int32_t net_k_plexes(const GraphHandle *h, int64_t k, int64_t min_size,
                     double max_steps, int32_t threads, NodeSetsHandle **out,
                     char *err, size_t err_len);

/* Edge connectivity (number of edge-disjoint paths) between every pair of
 * nodes of the undirected simple graph: an n-by-n symmetric matrix with
 * zero diagonal, written to out in column-major order. Lambda sets are
 * the classes of nodes joined by connectivity >= k. */
int32_t net_lambda_connectivity(const GraphHandle *h, double *out, char *err,
                                size_t err_len);

/* Factions: a partition into k non-empty groups (1 <= k <= n) by tabu
 * search from `starts` random partitions. fit: 0 = fewest errors (absent
 * ties within groups plus ties between groups), 1 = largest correlation
 * with the ideal pattern (needs 1 < k < n). A start ends after max_stall
 * steps without improvement; tenure == 0 chooses the default tabu
 * tenure. Writes 0-based labels (numbered by first appearance), the
 * number of errors, and the correlation (NaN where undefined) of the
 * partition. The same seed gives the same result on any thread count. */
int32_t net_factions(const GraphHandle *h, int64_t k, int32_t fit,
                     int64_t starts, int64_t max_stall, int64_t tenure,
                     uint64_t seed, int32_t threads, int64_t *membership,
                     double *errors, double *correlation, char *err,
                     size_t err_len);

/* Girvan-Newman divisive clustering of an undirected network (edges
 * removed in order of highest betweenness). Writes 0-based labels of the
 * partition of highest modularity to membership[0..n) and its modularity
 * to *q; the removed edges (smaller node first) to removed_from /
 * removed_to, the modularity before any removal and after each removal
 * to path, and the number of removals to *removed. removed_from and
 * removed_to must hold m values and path m + 1, where m is the number of
 * edges. */
int32_t net_girvan_newman(const GraphHandle *h, int32_t threads,
                          int64_t *membership, double *q, int64_t *removed_from,
                          int64_t *removed_to, double *path, int64_t *removed,
                          char *err, size_t err_len);

/* Newman's leading eigenvector communities of an undirected network (at
 * most 3000 nodes): 0-based labels to membership[0..n), modularity to *q,
 * number of divisions to *splits. Weights (weighted != 0) must be
 * non-negative. */
int32_t net_leading_eigen(const GraphHandle *h, int32_t weighted,
                          int64_t *membership, double *q, int64_t *splits,
                          char *err, size_t err_len);

/* Clauset-Newman-Moore greedy modularity communities of an undirected
 * network: 0-based labels to membership[0..n), modularity to *q.
 * resolution is the Reichardt-Bornholdt parameter (1 for plain
 * modularity). */
int32_t net_fast_greedy(const GraphHandle *h, double resolution,
                        int32_t weighted, int64_t *membership, double *q,
                        char *err, size_t err_len);

/* Discrete core-periphery fit (Borgatti and Everett 1999) by tabu search:
 * core[i] is 1 for core nodes and 0 otherwise, and the correlation
 * between the data and the pattern is written to *correlation. links:
 * 0 = core-periphery pairs left out of the fit, 1 = expected tied,
 * 2 = expected untied. weighted != 0 uses total weights between nodes as
 * values. Other arguments as in net_factions. */
int32_t net_core_periphery(const GraphHandle *h, int32_t links,
                           int32_t weighted, int64_t starts, int64_t max_stall,
                           int64_t tenure, uint64_t seed, int32_t threads,
                           int32_t *core, double *correlation, char *err,
                           size_t err_len);

/* Continuous core-periphery fit: coreness c (non-negative, unit
 * Euclidean norm) maximizing the correlation between the data and
 * c_i c_j, by BFGS from `starts` points (the degree start plus seeded
 * random ones), at most max_iter iterations each. Writes coreness to
 * coreness[0..n), the correlation to *correlation, and the iterations of
 * the best start to *iterations. */
int32_t net_core_periphery_continuous(const GraphHandle *h, int32_t weighted,
                                      int64_t max_iter, int64_t starts,
                                      uint64_t seed, int32_t threads,
                                      double *coreness, double *correlation,
                                      int64_t *iterations, char *err,
                                      size_t err_len);

/* Compactness (Borgatti 2006): the mean of 1 / d over ordered pairs of
 * distinct nodes, unreachable pairs contributing 0; 1 minus the
 * distance-weighted fragmentation. mode: 1 = out, 2 = in, 3 = all.
 * weights_as: 0 = distance, 1 = strength. NaN below two nodes. */
int32_t net_compactness(const GraphHandle *h, int32_t mode, int32_t weighted,
                        int32_t weights_as, int32_t threads, double *out,
                        char *err, size_t err_len);

/* Small-world quotient (Watts and Strogatz 1998): clustering and mean
 * distance of the undirected simple graph against `reps` G(n, m) graphs
 * of the same size. Writes five values to out: clustering, mean distance,
 * mean random clustering, mean random distance, and the quotient
 * (C / C_random) / (L / L_random). local != 0 uses mean local clustering
 * instead of transitivity. */
int32_t net_small_world(const GraphHandle *h, int32_t local, int64_t reps,
                        uint64_t seed, int32_t threads, double *out, char *err,
                        size_t err_len);

/* ---- Random graphs -------------------------------------------------
 * Every generator returns `reps` edge lists on nodes 0..n in one handle;
 * replicate r uses random stream (seed, r), so a seed reproduces every
 * replicate on any platform and thread count. Edges are in increasing
 * (from, to) order with from < to in undirected graphs. */
typedef struct EdgeListsHandle EdgeListsHandle;

/* Number of lists, or -1 for a NULL handle. */
int64_t net_edge_lists_count(const EdgeListsHandle *h);
/* Total number of edges over all lists, or -1 for a NULL handle. */
int64_t net_edge_lists_total(const EdgeListsHandle *h);
/* The edges of list i are positions offsets[i] .. offsets[i + 1]) of from
 * and to. offsets must hold count + 1 values, from and to total values. */
int32_t net_edge_lists_copy(const EdgeListsHandle *h, int64_t *offsets,
                            int64_t *from, int64_t *to, char *err,
                            size_t err_len);
void net_edge_lists_free(EdgeListsHandle *h);

/* G(n, p): each possible pair tied independently with probability p.
 * loops != 0 makes self-loops possible. */
int32_t net_random_gnp(int64_t n, double p, int32_t directed, int32_t loops,
                       int64_t reps, uint64_t seed, int32_t threads,
                       EdgeListsHandle **out, char *err, size_t err_len);

/* G(n, m): m distinct pairs drawn uniformly. */
int32_t net_random_gnm(int64_t n, int64_t m, int32_t directed, int32_t loops,
                       int64_t reps, uint64_t seed, int32_t threads,
                       EdgeListsHandle **out, char *err, size_t err_len);

/* Configuration model: degrees holds n degrees of an undirected graph
 * (in_degrees NULL) or n out-degrees with in_degrees holding n
 * in-degrees (directed). simple == 0 returns the multigraph (loops and
 * parallel edges kept, a loop counting 2 toward an undirected degree);
 * simple != 0 repeats each matching up to max_tries times until it is
 * simple, and fails if none is. */
int32_t net_random_configuration(int64_t n, const int64_t *degrees,
                                 const int64_t *in_degrees, int32_t simple,
                                 int64_t max_tries, int64_t reps, uint64_t seed,
                                 int32_t threads, EdgeListsHandle **out,
                                 char *err, size_t err_len);

/* U|MAN: a directed graph drawn uniformly from those with `mutual`
 * mutual and `asym` asymmetric dyads (the rest null). */
int32_t net_random_uman(int64_t n, int64_t mutual, int64_t asym, int64_t reps,
                        uint64_t seed, int32_t threads, EdgeListsHandle **out,
                        char *err, size_t err_len);

/* Watts-Strogatz: a ring lattice joining each node to its k nearest
 * neighbors (k even), each edge rewired with probability p. */
int32_t net_random_small_world(int64_t n, int64_t k, double p, int64_t reps,
                               uint64_t seed, int32_t threads,
                               EdgeListsHandle **out, char *err,
                               size_t err_len);

/* Barabasi-Albert: from a star on m + 1 nodes, each new node joins m
 * distinct nodes drawn with probability proportional to degree. */
int32_t net_random_preferential(int64_t n, int64_t m, int64_t reps,
                                uint64_t seed, int32_t threads,
                                EdgeListsHandle **out, char *err,
                                size_t err_len);

/* Stochastic block model: k blocks of sizes[r] nodes (numbered block by
 * block) and tie probabilities probs[r * k + s] (row-major; symmetric for
 * undirected graphs). */
int32_t net_random_sbm(int64_t k, const int64_t *sizes, const double *probs,
                       int32_t directed, int32_t loops, int64_t reps,
                       uint64_t seed, int32_t threads, EdgeListsHandle **out,
                       char *err, size_t err_len);
