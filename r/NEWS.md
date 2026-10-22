# knappnet (development version)

First release. Networks are a pair of data frames, nodes and edges, and
every measure runs in a compiled, multithreaded engine written in Rust with
no dependencies. Each measure is validated against igraph, sna, NetworkX,
and a direct evaluation of its definition; where those tools disagree, the
`convention` argument selects their behavior and `?conventions` lists every
difference.

## Networks

* `netframe()` builds a network from an edge table and an optional node
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
* Conversions to and from igraph (`to_igraph()`, `as_netframe()`),
  tidygraph (`to_tbl_graph()`), statnet's network (`to_network()`), and
  dense or sparse adjacency matrices (`as.matrix()`), keeping node, edge,
  and network attributes.
* `read_network()` and `write_network()` read and write GraphML, GEXF, GML,
  Pajek, UCINET DL, JSON node-link, ORA's DyNetML (multimode
  meta-networks), and CSV/TSV edge lists, and write DOT. Text is UTF-8 on
  every platform.

## Measures

* Centrality: `degree()`, `strength()`, `betweenness()` (Brandes),
  `closeness()`, `harmonic_centrality()`, and `eigenvector_centrality()`,
  with the conventions of igraph, NetworkX, sna, and UCINET for
  multigraphs, unreachable nodes, and normalization.
* Whole networks: `graph_density()`, `components()` (weak and strong),
  `fragmentation()`, `centralization()` (Freeman), and `triad_census()`
  (Batagelj and Mrvar).
* Cohesive subgroups: `coreness()` (k-cores), `max_cliques()` and
  `clique_number()`, `modularity()` (with resolution and the directed
  form), and `cluster_louvain()` with a reproducible seed.
* Brokerage: Burt's `constraint()`, `effective_size()`, and `efficiency()`,
  and Gould-Fernandez `brokerage()` roles.
* Positions: `structural_equivalence()` (Euclidean, correlation, Hamming)
  and `concor()`, which warns when rounding error may decide a split.
* Hypothesis tests: `qap_cor()` and `mrqap()` (double semi-partialling or
  Y-permutation), with exact tests over all permutations for small
  networks and seeded, thread-independent random permutations otherwise.

## Reliability

* Results are reproducible for a given number of threads, and integer
  results never depend on it. `threads = 0` uses all cores, but only 2
  when `_R_CHECK_LIMIT_CORES_` is set (as by `R CMD check --as-cran`) and
  at most `OMP_THREAD_LIMIT` when that is set.
* Inputs are checked before any computation, and errors name the argument
  and say how to fix the problem. Missing or empty identifiers, duplicated
  node identifiers, and edges to unknown nodes are errors naming the rows.
  Missing or non-finite weights are kept in the data but stop every
  measure, and zero or negative weights are checked by each measure that
  uses weights; errors name the row of the edge table. See "Missing values
  and invalid input" in `?netframe`.
