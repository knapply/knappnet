# What is supported

The complete list, generated from NetworkX's registry of dispatchable
functions, is [docs/drafts/networkx_coverage.md](../../docs/drafts/networkx_coverage.md).
For each of NetworkX's functions it gives the status in knappnet:
**implemented** (every call runs on the engine), **partial** (the listed
arguments or inputs run in NetworkX), or **missing** (runs in NetworkX).

In summary, knappnet runs these on its engine (some only for the inputs
the coverage list names; the rest run in NetworkX):

- Centrality: degree, in- and out-degree, betweenness and edge
  betweenness (and their subset versions), load and edge load, closeness,
  harmonic, eigenvector (power iteration and NumPy versions), Katz (both
  versions), PageRank, HITS, percolation, current-flow betweenness (node,
  edge, and subset versions) and closeness, communicability betweenness,
  subgraph centrality and the Estrada index, second-order, Laplacian, and
  trophic measures, resistance distance, Kemeny constant, dispersion,
  VoteRank, reaching centralities (unweighted), group degree, closeness,
  and betweenness, and the bipartite degree, closeness, and betweenness
  centralities.
- Paths: single-source, single-target, and all-pairs shortest path
  lengths (unweighted, Dijkstra, and Bellman-Ford), `shortest_path_length`,
  Floyd-Warshall, negative-cycle detection, `has_path`, descendants and
  ancestors, average shortest path length, Wiener, Schultz, Gutman, and
  hyper-Wiener indices, global and local efficiency, harmonic diameter,
  closeness vitality, eccentricity, diameter, radius, center, periphery.
  Functions that return the paths themselves run on the engine when every
  shortest path is unique.
- Structure: connected, strongly connected, and weakly connected
  components (and their counts and tests), `node_connected_component`,
  articulation points, bridges, local bridges, biconnected components,
  condensation, attracting components, semiconnectedness, k-edge
  components, bridge components, k-components, core numbers and the k-core,
  k-shell, k-crust, k-corona, k-truss, and onion-layer subgraphs, maximal
  cliques, triadic census, transitivity, clustering, square clustering,
  generalized degree, average clustering, reciprocity, Burt's constraint
  and effective size.
- Flow and connectivity: maximum flow value, minimum cut and cut value,
  minimum cost flow cost, Gomory-Hu trees, node and edge connectivity
  (global and local), average and all-pairs node connectivity, minimum
  s-t node and edge cuts.
- Trees, DAGs, and cycles: minimum and maximum spanning trees and edges
  (Kruskal and Prim), topological sorts and generations, DAG longest
  paths, transitive closure and reduction, dominators, lowest common
  ancestors, cycle bases, girth, `find_cycle`, and the unnormalized rich
  club coefficient.
- Mixing and cuts: degree, attribute, and numeric assortativity, mixing
  dicts and matrices, average neighbor degree, average degree
  connectivity, cut size, volume, conductance, expansion measures.
- Similarity and isomorphism: the link prediction indices,
  Weisfeiler-Lehman hashes, SimRank, graph edit distance, the cheap
  isomorphism filters, and the VF2++ yes/no tests.
- Matching and colouring: maximal matching and the matching checks,
  greedy colouring with the largest-first and saturation strategies.
- Communities: modularity and Louvain communities.

Anything else, and any call knappnet declines, runs in NetworkX with the
same result it would always have had. [Differences](differences.md)
describes the few places where knappnet's results are not identical.
