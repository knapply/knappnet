# What is supported

The complete list, generated from NetworkX's registry of dispatchable
functions, is [docs/drafts/networkx_coverage.md](../../docs/drafts/networkx_coverage.md).
For each of NetworkX's functions it gives the status in knappnet:
**implemented** (every call runs on the engine), **partial** (the listed
arguments or inputs run in NetworkX), or **missing** (runs in NetworkX).

In summary, knappnet runs these on its engine:

- Centrality: degree, in- and out-degree, betweenness and edge
  betweenness, closeness, harmonic, eigenvector (power iteration and
  NumPy versions), Katz (both versions), PageRank, HITS, current-flow
  closeness (information centrality), group degree, closeness, and
  betweenness, and the bipartite degree, closeness, and betweenness
  centralities.
- Paths: single-source, single-target, and all-pairs shortest path
  lengths (unweighted and Dijkstra), `shortest_path_length`, `has_path`,
  descendants and ancestors, average shortest path length, Wiener index,
  global efficiency, eccentricity, diameter, radius, center, periphery.
- Structure: connected, strongly connected, and weakly connected
  components (and their counts and tests), `node_connected_component`,
  articulation points, bridges, biconnected components, core numbers,
  maximal cliques, triadic census, transitivity, clustering, average
  clustering, reciprocity, Burt's constraint and effective size.
- Mixing: degree, attribute, and numeric assortativity, degree Pearson
  correlation.
- Communities: modularity and Louvain communities.

Anything else, and any call knappnet declines, runs in NetworkX with the
same result it would always have had. [Differences](differences.md)
describes the few places where knappnet's results are not identical.
