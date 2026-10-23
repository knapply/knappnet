# knappnet (development version)

First release. Networks are a pair of data frames, nodes and edges, and
every measure runs in a compiled, multithreaded engine written in Rust with
no dependencies. Each measure is validated against igraph, sna, NetworkX,
and a direct evaluation of its definition; where those tools disagree, the
`convention` argument selects their behavior and `?conventions` lists every
difference.

## Networks

* `net_frame()` builds a network from an edge table and an optional node
  table, with column roles for identifiers, endpoints, and weights.
  Directed and undirected networks, self-loops, and parallel edges are
  supported, and `print()` reports self-loops, parallel edges, and missing
  weights.
* Weighted and binary networks are explicit: a network is weighted only
  when it has a weight column, and `weight_type` says once whether weights
  are distances or strengths (tie intensities, read as lengths 1/w by
  path-based measures). The fastest exact algorithm is chosen
  automatically: breadth-first search for binary networks, Dial's bucket
  queue for small whole-number lengths, and Dijkstra's algorithm otherwise.
* Conversions to and from igraph (`net_to_igraph()`, `net_from()`),
  tidygraph (`net_to_tbl_graph()`), statnet's network (`net_to_network()`), and
  dense or sparse adjacency matrices (`as.matrix()`), keeping node, edge,
  and network attributes.
* `net_read()` and `net_write()` read and write GraphML, GEXF, GML,
  Pajek, UCINET DL, JSON node-link, ORA's DyNetML (multimode
  meta-networks), and CSV/TSV edge lists, and write DOT. Text is UTF-8 on
  every platform.

## Measures

* Centrality: `net_degree()`, `net_strength()`, `net_betweenness()` (Brandes),
  `net_closeness()`, `net_harmonic()`, and `net_eigen()`,
  with the conventions of igraph, NetworkX, sna, and UCINET for
  multigraphs, unreachable nodes, and normalization.
* Whole networks: `net_density()`, `net_components()` (weak and strong),
  `net_fragmentation()`, `net_centralization()` (Freeman), and `net_triad_census()`
  (Batagelj and Mrvar).
* Cohesive subgroups: `net_coreness()` (k-cores), `net_max_cliques()` and
  `net_clique_number()`, `net_modularity()` (with resolution and the directed
  form), and `net_louvain()` with a reproducible seed.
* Brokerage: Burt's `net_constraint()`, `net_effective_size()`, and `net_efficiency()`,
  and Gould-Fernandez `net_brokerage()` roles.
* Positions: `net_struct_equiv()` (Euclidean, correlation, Hamming)
  and `net_concor()`, which warns when rounding error may decide a split.
* Hypothesis tests: `net_qap_cor()` and `net_mrqap()` (double semi-partialling or
  Y-permutation), with exact tests over all permutations for small
  networks and seeded, thread-independent random permutations otherwise.

## Plotting

* `plot()` draws a netframe with a stress layout and appearance mapped from
  node and edge columns; spatial networks are drawn at their locations.
  Every plot uses one colour scheme: Okabe and Ito's colour-blind-safe
  palette for categories (HCL "Dark 3" beyond eight), viridis for scores,
  and a single-hue blue scale for matrices.
* Overlays on the network: `net_plot_centrality()` (nodes sized and colored by
  a score), `net_plot_partition()` (colored by a grouping, with modularity),
  `net_plot_coreness()` (concentric k-core shells), and `net_plot_ego()`.
* Matrices: `net_plot_matrix()` (adjacency or incidence matrix, sorted by a
  partition with block boundaries), `net_plot_mixing()`, and a `plot()`
  method for `net_blockmodel()` (permuted matrix and block densities with the
  image outlined).
* Distributions: `net_plot_degree()` (counts or a log-log CCDF),
  `net_plot_components()`, `net_plot_triads()` against the counts
  expected under Holland and Leinhardt's U|MAN model, and `net_plot_scores()`
  (ranked dot plot).
* Diagnostics: `plot()` methods for `net_qap_cor()` (permutation distribution
  with the observed statistic), `net_mrqap()` (coefficient plot with
  permutation intervals, or per-coefficient distributions), and
  `net_mds()` (configuration, Shepard diagram, scree plot; the result
  now has class `"classical_mds"` with a `print()` method).
* `net_plot_robustness()` tracks the largest component or fragmentation as
  nodes are removed by degree, betweenness, at random, or in a given order
  such as the key players.
* Temporal networks: `net_plot_timeline()` of contacts and `net_plot_rolling()` for
  measures over windows.

## Reliability

* Results are reproducible for a given number of threads, and integer
  results never depend on it. `threads = 0` uses all cores, but only 2
  when `_R_CHECK_LIMIT_CORES_` is set (as by `R CMD check --as-cran`) and
  at most `OMP_THREAD_LIMIT` when that is set.
* Inputs are checked before any computation by one set of shared checks,
  so the same kind of mistake gives the same message everywhere; errors
  name the argument and say how to fix the problem. A `seed` means the
  same thing in every function: a whole number between 0 and 2^53, or
  `NULL` to draw one from R's generator. Missing or empty identifiers, duplicated
  node identifiers, and edges to unknown nodes are errors naming the rows.
  Missing or non-finite weights are kept in the data but stop every
  measure, and zero or negative weights are checked by each measure that
  uses weights; errors name the row of the edge table. See "Missing values
  and invalid input" in `?netframe`.
