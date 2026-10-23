# Scale notes: community detection

For each method as implemented in `crates/engine/src/metrics/partition/`:
memory per node and per edge, passes over the edges, whether it could
run on a memory-mapped or edge-streamed graph with O(n) RAM, and known
scalable variants. n nodes, m edges, k groups; "CSR" is the engine's
adjacency (about 12 bytes per arc, so 24 bytes per undirected edge with
weights).

- **Walktrap** (Pons and Latapy 2006): dense walk vectors, 8n bytes per
  node (8n^2 in all), plus the CSR and a heap of O(m) candidates. t
  sparse matrix-vector passes over the edges per node (n t passes in
  all), then one Ward update per merge costing O(n) per affected pair.
  Not O(n) RAM: the vectors are the method. Pons and Latapy (2006) bound
  memory by keeping only some vectors and recomputing the rest, still
  quadratic in the worst case; for large networks Louvain or Leiden is
  the practical substitute.
- **Label propagation** (Raghavan, Albert, and Kumara 2007): O(n) labels
  plus the CSR; one pass over the edges per sweep, with a per-node label
  tally. The asynchronous variants need random access to neighborhoods
  (CSR in memory or memory-mapped, O(n) RAM beyond the mapping); a
  synchronous sweep can run over an edge stream sorted by node with O(n)
  RAM for two label arrays. The fast variant (Traag and Subelj 2023)
  revisits only nodes whose neighborhood changed, usually far fewer than
  full sweeps. Scalable variants: parallel label propagation (Staudt and
  Meyerhenke 2016), balanced label propagation for partitioning (Ugander
  and Backstrom 2013).
- **Infomap** (Rosvall and Bergstrom 2008; Rosvall, Axelsson, and
  Bergstrom 2009): flows per node and per arc (about 16 bytes per arc
  for the flow lists in both directions), the PageRank iteration for
  directed networks (O(n) vectors, one pass over the edges per power
  step), then local moving (one pass per sweep) on each aggregation
  level; submodule movements rebuild subnetworks, briefly doubling the
  arc storage. Random access is needed (memory-mapped CSR possible, not
  streaming). Scalable variants: RelaxMap, a parallel relaxation of the
  core moves (Bae, Halperin, West, Rosvall, and Howe 2013).
- **Leiden** (Traag, Waltman, and van Eck 2019): CSR plus O(n) group
  totals, a node queue, and the refined partition; fast local moving
  touches each arc once per visit of its endpoint, and each aggregation
  shrinks the network. Needs random access (memory-mapped CSR with O(n)
  RAM beyond it is feasible for the first level; later levels are
  smaller and held in memory). Scalable variants: parallel Louvain-type
  local moving (Staudt and Meyerhenke 2016); Traag, Waltman, and van Eck
  (2019) report runs on networks with billions of edges.
- **Spinglass** (Reichardt and Bornholdt 2006; Traag and Bruggeman
  2009): O(n) spins, O(q) totals per spin state, and O(q) scratch per
  update; one pass over the edges per temperature (about 460 sweeps with
  the defaults) plus the final quench. The heat-bath update reads a
  node's neighborhood, so a memory-mapped CSR with O(n + q) RAM works; a
  sorted edge stream works for the parallel (synchronous) update. Cost is
  dominated by the schedule length times m + n q; no standard scalable
  variant, and the method is mainly of interest for moderate networks.
- **Exact maximum modularity** (Brandes et al. 2008; branch and bound
  with Russian doll search, Verfaillie, Lemaitre, and Schiex 1996):
  dense pair scores, 8 c^2 bytes for a component of c nodes, solved one
  component at a time; time exponential in c. Limited to components of
  at most 60 nodes by design. Scale is not a goal: exact solvers based on
  column generation (Aloise et al. 2010) reach a few thousand nodes on
  sparse networks.
- **Fluid communities** (Pares et al. 2018): O(n) community labels and
  O(k) densities, one pass over the edges per sweep; random access to
  neighborhoods (memory-mapped CSR with O(n + k) RAM). The paper presents
  the method as scalable (linear work per sweep) and it parallelizes in
  the same way as asynchronous label propagation.
- **Voronoi communities** (Deritei et al. 2014; Molnar et al. 2024):
  edge lengths (8 bytes per arc, both directions), O(n) densities, and
  for a fixed radius one truncated shortest-path search per generator
  plus one multi-source search, O(m log n) in all. Choosing the radius
  here runs one full shortest-path search per node (n m log n time) and
  keeps the distinct distance values (up to n^2 in the worst case with
  real-valued lengths), then evaluates up to 1,000 radii. A fixed radius
  could run on a memory-mapped CSR with O(n) RAM; the automatic radius
  cannot. Molnar et al. (2024) choose the radius with Brent's method,
  which needs only a few evaluations and is the scalable choice.
- **Partition comparison** (Meila 2007; Danon et al. 2005; van Dongen
  2000; Rand 1971; Hubert and Arabie 1985; Vinh, Epps, and Bailey 2010):
  one pass over the two label vectors builds a sparse contingency table
  (at most min(n, k_a k_b) cells), so any stream of label pairs works
  with O(k_a + k_b + cells) RAM. The adjusted mutual information's
  expected value costs O(n k_a k_b) in the worst case.
- **Dendrogram cutting**: O(n) union-find over the merges, one pass over
  the merge list (streamable); modularity along the dendrogram costs one
  pass over the edges plus O(1) per merge.

References not otherwise cited in the code (not yet in
`docs/references.bib`; verify before citing in a module): Staudt, C. L.
and Meyerhenke, H. (2016), Engineering parallel algorithms for community
detection in massive networks, IEEE Transactions on Parallel and
Distributed Systems 27(1); Ugander, J. and Backstrom, L. (2013), Balanced
label propagation for partitioning massive graphs, WSDM 2013; Bae, S.-H.,
Halperin, D., West, J., Rosvall, M., and Howe, B. (2013), Scalable
flow-based community detection for large-scale network analysis, IEEE
ICDM Workshops 2013; Aloise, D., Cafieri, S., Caporossi, G., Hansen, P.,
Perron, S., and Liberti, L. (2010), Column generation algorithms for
exact modularity maximization in networks, Physical Review E 82, 046112.
