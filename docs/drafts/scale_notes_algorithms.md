# Scale notes: community, position, layout, permutation, generator and dense kernels

Working notes for docs/SCALE.md (kernel ledger). Scope: Louvain, fast
greedy, leading eigenvector, Girvan-Newman, maximal cliques, factions and
core-periphery, structural equivalence and CONCOR, REGE (CATREGE), stress
layouts, QAP/MRQAP, random graphs, dense linear algebra. Path and
traversal kernels are in the paths workstream's notes; Leiden, Infomap,
label propagation, Walktrap and the other community methods in
scale_notes_communities.md.

Measurements: `cargo run --release --example bench_models` on a 2-vCPU
Xeon (2.1 GHz, AVX-512 capable), 7 GB RAM, one thread unless stated; heap
is the peak of live heap bytes from a counting allocator. Sparse = average
degree 8, medium = density 2%, dense = density 30%. Timings on this
virtual machine vary by up to 2x between runs; before/after figures are
the better of two alternating runs of each version (medians of 3), and
untouched workloads moved by -10% to +80% in the same runs, so ratios
under ~1.3x are noise. Thread scaling is not measurable here (2 vCPUs on
one core).

## Per kernel

### Louvain (Blondel et al. 2008; fast local move, Traag et al. 2019)

- Access: per level, local-move passes visit each node's adjacency once
  (sequential over CSR), with random reads and writes of O(n) community
  arrays (community id, community total degree, per-node scratch
  accumulator). Aggregation is one pass over edges, a sparse S' A S
  (bucket by community, merge in a dense accumulator; Gilbert, Moler and
  Schreiber 1992).
- Working set: CSR (8-12 bytes per arc) + ~5 arrays of n. Measured peak
  8.5 MB at n = 20k, m = 80k (106 B/edge incl. the input graph copy);
  14 MB at n = 1k, m = 150k (96 B/edge). Allocations per run fell from
  98k to 208 after the aggregation rewrite; time 0.059 s -> 0.040 s
  (sparse), within noise elsewhere.
- Memory-mapped CSR: viable for level 0 (sequential row reads; community
  arrays random but O(n)). Coarse levels are small and stay in RAM.
- Semi-external: yes. Node state O(n); each pass streams edges sorted by
  source. The fast local move queue makes later passes touch few rows:
  rows are read on demand, which on a mapped file is random I/O; a
  pass-ordered (plain sweep) variant keeps I/O sequential at the cost of
  more passes.
- Ordering: relabelling nodes by BFS or degree order (Wei et al. 2016,
  Gorder) improves locality of the community-array reads; the result then
  depends on the order (Louvain is order dependent), so any reordering must
  be applied to the visiting order explicitly and documented.
- Streaming alternative: one-pass edge-stream clustering with O(n) state
  (Hollocou et al. 2017, SCoDA) gives lower-quality partitions but needs no
  stored graph. Incremental updates: re-run local moves seeded with the
  previous partition (Aynaud and Guillaume 2010) on the affected nodes;
  the frontier queue of the fast local move already supports starting from
  an arbitrary node set.
- Billion edges: CSR of 10^9 arcs with 32-bit targets = 4 GB (+ weights 4
  GB in f32); community state for 10^8 nodes ~2-3 GB. Feasible on a 32-64
  GB machine in RAM, or semi-external from a mapped CSR on SSD.

### Fast greedy (Clauset, Newman and Moore 2004)

- Access: priority queue of candidate merges; per-group ordered maps
  (BTreeMap) merged smaller-into-larger. Random access throughout.
- Change: stale queue entries are re-evaluated lazily (a stored gain is an
  upper bound unless the tie weight changed, which always pushes a fresh
  entry), so only the absorbed group's neighbours are re-pushed per merge.
  Same merges, ties and output as before.
- Measured: sparse n = 10k, m = 40k: 2.50 s -> 0.077 s (33x), peak heap
  169 MB -> 4.4 MB; medium n = 2k, m = 40k: 1.53 s -> 0.039 s (39x); dense
  n = 600, m = 54k: 0.093 s -> 0.039 s.
- Working set: O(m) maps + O(m) queue. ~110 B/edge measured.
- Semi-external: no (global queue, arbitrary merges). Out-of-core variants
  are not practical; at scale use Louvain/Leiden.

### Leading eigenvector (Newman 2006)

- Access: dense n x n modularity matrix, dense symmetric eigensolver
  (EISPACK tred2/tql2, O(n^3), O(n^2) memory); limited to 3000 nodes.
- Change: the solver now works on the transpose of the EISPACK working
  matrix, so the Householder and rotation loops run over contiguous rows
  (same operations, same order, bit-identical). n = 1000 decomposition:
  3.24 s -> 1.44 s; leading eigenvector n = 1000: 4.1-4.5 s -> 1.9-2.0 s;
  classical MDS n = 1000: 2.92 s -> 1.29 s. Peak memory 3 -> 4 n^2
  doubles (the transposed copy).
  Compiling these loops for AVX2/AVX-512 made no difference (rotations are
  bound by memory traffic and dependent sums), so no wide copy is kept.
- Scalable variant: Lanczos on the implicit B x = A x - k (k'x)/2m costs
  O(m) per iteration (as igraph does with ARPACK); not adopted because it
  changes which vector is returned when the leading eigenvalue is
  degenerate, which changes partitions. Would remove the 3000-node limit
  and make the method semi-external (each iteration one edge pass).

### Girvan-Newman (Girvan and Newman 2002)

- Access: after each removal, Brandes accumulation from every node of the
  affected component; O(m^2 n) total. Each source pass is BFS over the
  remaining adjacency.
- Change: removed edges are deleted from adjacency rows (order kept), so
  BFS no longer tests a liveness flag per arc; distances in 32 bits.
  sparse n = 300: 6.58 s -> 4.14 s (allocations 1.47M -> 43k with the
  queue reuse).
- Not a large-graph method: quadratic in m by definition. Sampled
  betweenness (Riondato and Kornaropoulos 2016) would make each round
  cheaper but changes results.

### Maximal cliques (Bron-Kerbosch with Tomita pivoting, Eppstein et al. 2010 ordering)

- Access: degeneracy order, then per-vertex subproblems on the neighbourhood
  (local adjacency, bitsets for dense subproblems). Per-vertex subproblems
  are independent: natural for out-of-core processing, one vertex's
  later-neighbourhood at a time (2-hop working set, O(d^2) for degeneracy
  d).
- Measured: n = 20k sparse 0.035 s (1.2 MB); medium n = 3k 0.18 s (3.6 MB).
- Semi-external: yes with an edge stream sorted by degeneracy order, but
  the output can be exponential; streaming the output is required.

### Factions, core-periphery (tabu search)

- Access: dense or bit adjacency over n <= a few thousand; O(n) work per
  move evaluation. Measured: n = 400 sparse 0.08 s, 0.1 MB. Not large-graph
  methods (the objective is defined on all n^2 cells). Sampling-based
  objective estimates would be the only route; not planned.

### Structural equivalence and CONCOR (Breiger, Boorman and Arabie 1975)

- Access: all pairs of tie profiles. Sparse path: merge of sorted rows,
  O(n m) total. Dense binary path (new): bit rows, popcount of AND/XOR,
  O(n^2 * n/32 words); counts are exact so results are identical. Only
  i < j is evaluated (every measure is symmetric bit for bit).
- CONCOR iterations: b x b correlation matrices, O(b^3) per iteration;
  now a tiled product (8 x 4 register tiles over packed rows), compiled
  for AVX2 at run time with a baseline fallback; each sum is still in order
  so results are identical.
- Measured (depth 3, directed): dense n = 1k, m = 300k: 26.4 s -> 1.26 s
  (21x; the initial profile correlations were 75% of the old time).
  Sparse n = 2k, m = 16k: 124.7 s -> 8.4 s (15x). The AVX2 copy of the
  tile alone is worth 2x (20.4 s -> 10.9 s back to back); AVX2 copies of
  the eigensolver and matmul gained nothing and were dropped.
- Working set: 4-5 dense b x b matrices (145 MB at b = 2000). Inherently
  O(n^2) memory; at scale the iterated-correlation fixed point can be
  approximated by spectral bisection of the correlation matrix (not the
  same method; documented divergence would be needed).
- Billion edges: not applicable (O(n^2) output).

### REGE / CATREGE (Borgatti and Everett 1993)

- Access: rounds of class refinement; each round one pass over arcs,
  signatures hashed per node. This is colour refinement on (tie type,
  class) pairs: the same skeleton as 1-WL refinement and Weisfeiler-Lehman
  graph kernels (Shervashidze et al. 2011), so a WL hash/kernel can share
  the round structure. Measured: n = 2k sparse 0.011 s; dense n = 1k,
  m = 300k 0.06 s (26 MB).
- Semi-external: yes. O(n) class labels; each round streams arcs grouped by
  source (signature built per node) — sequential with a sorted CSR.
- Rounds are at most n, usually few; I/O = rounds x edge passes.

### Stress layouts (Gansner, Koren and North 2005; sparse stress, Ortmann, Klimenta and Brandes 2017)

- Full stress: all-pairs distances (n BFS/Dijkstra), dense n x n
  reciprocal distances (8n^2 bytes: 8 MB at n = 1k, 800 MB at 10k), then
  majorization passes over all pairs (87-94% of the time at n = 1k:
  0.47 s sparse, 0.64 s dense). Tried: four rows per pass in vector lanes
  (interleaved rows, branch-free selects, AVX2 copy); no measurable gain,
  reverted. A pivot-based start with fewer majorization passes, or a
  stochastic-gradient solver (Zheng, Pawar and Goodman 2019), are the
  remaining levers; both change results.
- Sparse stress: k pivots (k x n distances in f32), edge terms, colour
  classes for Gauss-Seidel updates. Working set O(k n + m); n = 10k:
  0.78 s, 19 MB. Semi-external: pivot BFS passes are edge streams; the
  update loop reads each node's arcs and the k pivot terms (row-major), so
  a mapped CSR with node-major term storage streams sequentially per
  colour class only if classes are stored contiguously (relabel by class).
- Billion edges: sparse stress with k = 200 needs 200 BFS passes and
  800 n bytes of terms: n = 10^8 -> 80 GB of terms, so terms must be
  recomputed or quantized (8-bit) at that size; multilevel (coarsening)
  layouts are the realistic route.

### QAP and MRQAP (Krackhardt 1988; Dekker, Krackhardt and Snijders 2007)

- Access: n x n matrices vectorized to dyads; per permutation a gather and
  a statistic (QAP) or a least-squares refit (MRQAP).
- Change: Householder reflectors of the fixed columns are computed once.
  Y-permutation refits only transform the permuted response (O(p N)
  instead of O(p^2 N) per draw); DSP refits only the permuted residual
  column. Operations and order as a fresh fit, so results are identical;
  residuals are no longer formed inside the loop. n = 100, 3 predictors,
  1000 permutations: DSP 1.43 s -> 0.73 s; Y-permutation 0.33 s -> 0.14 s.
- Working set: (p + 1) N doubles per thread for N = n^2 dyads. n = 10k
  -> N = 10^8, 0.8 GB per column: dyad streaming in blocks (the reflector
  application is a sequence of dot products and axpys over N, i.e. one
  streaming pass per column) makes MRQAP O(p) passes per permutation with
  O(n) RAM for the permutation itself.

### Random graphs (Batagelj and Brandes 2005 skipping; Floyd sampling)

- G(n, p) by geometric skipping: O(n + m), sequential output; trivially
  streaming (edges can be written as generated, O(1) state per row).
  Dense n = 2k, p = 0.3: 20 draws 1.1 s, 107 MB (output-bound).
- Generation at billion-edge scale is I/O bound; parallel per-row streams
  (seeded per row block) keep reproducibility independent of threads.

### Dense linear algebra

- Matmul now blocked (64 x 512 panels), same accumulation order: n = 1000
  0.38 s -> 0.25 s. Eigendecomposition as above. Cholesky unchanged
  (0.14 s at n = 1000). All O(n^2) memory; not used on large graphs except
  through p x p systems.

## Primitives worth sharing (for streaming and incremental work)

- Bucket-by-row + dense accumulator (Louvain aggregation) is a general
  sparse "S' A S" / SpGEMM-with-dense-accumulator kernel; usable for
  coarsening, quotient graphs and WL-style label aggregation.
- Bit rows with popcount (structural equivalence) generalize to common
  neighbour counts, Jaccard/cosine similarity on dense graphs, and
  triangle counting by masked products (Azad, Buluc and Gilbert 2015).
- Lazy-priority merging (fast greedy) works for any agglomeration whose
  merge scores only fall when a group grows.
- QR prefix reuse applies to any repeated least squares with a fixed
  design block (permutation tests, jackknife over one column).

## What would make billion-edge graphs feasible on one machine

1. Mapped CSR with 64-bit offsets and 32-bit targets (4-8 bytes per arc)
   plus f32 weights; kernels audited for O(n) RAM state.
2. Louvain/Leiden, REGE/WL refinement, label propagation, components and
   k-core as edge-stream (semi-external) kernels: passes x 4-8 GB per pass
   at SSD speeds (2-7 GB/s) -> seconds per pass.
3. Sketches where exactness is not required: HyperLogLog neighbourhood
   functions (Boldi, Rosa and Vigna 2011) for distance-based measures;
   MinHash/bottom-k sketches of neighbour sets for structural-equivalence
   similarity (approximate, O(n k) memory instead of O(n^2) output).
4. Compression (WebGraph, Boldi and Vigna 2004) to bring 10^9 arcs to
   ~1-3 GB when sequential decoding suffices (Louvain sweeps, WL rounds).
5. Methods inherently O(n^2) (CONCOR, structural-equivalence matrices,
   full stress, QAP on all dyads) need blocking and output streaming, or
   sampling; they do not reach billion-edge scale as exact methods.

## Flows, connectivity, trees, DAGs, and cycles

Kernels added with the flow and connectivity work: `crates/engine/src/flow/`,
`metrics/connectivity/`, `metrics/trees.rs`, `metrics/dag.rs`,
`metrics/cycles.rs`. Measurements: `cargo run --release --example
bench_flows` on the btower bench runner (8 threads visible), median of 3,
random networks with average degree 8 (connectivity: 20; cycles: 6;
minimum cycle bases: 4) and weights 1..10; peak heap from the counting
allocator, including the input graph's own copy. NetworkX 3.7 figures were
taken on the 2-vCPU development VM on comparable random graphs, where the
engine runs about 2x slower than on btower, so the ratios quoted are
conservative.

### Maximum flow and minimum cut (Dinic 1970)

- One residual network (`flow::FlowNetwork`) serves every flow: arcs in
  reverse pairs, per-pair flow in one f64 array, CSR of arcs by tail built
  once; BFS levels plus DFS blocking flow with current-arc pointers.
  Working set ~50 B per arc (tail-indexed CSR, heads, capacities, flows)
  plus O(n) scratch. Measured: n = 10^5, m = 8 x 10^5: 0.12 s, 47 MB;
  n = 10^4: 8 ms (NetworkX `maximum_flow_value`, preflow-push: 0.57 s on
  the VM, ~40x).
- Unit-capacity networks (connectivity) need O(m min(n^(2/3), m^(1/2)))
  per flow; general capacities O(n^2 m) worst case, rarely approached on
  sparse graphs. A highest-label push-relabel with gap and global
  relabelling (Goldberg and Tarjan 1988; Cherkassky and Goldberg 1997) is
  the next step for large dense or adversarial instances; the API
  (`FlowNetwork::max_flow`, `reachable_from`, `reaching`) would not change.
- Semi-external: no; flows touch arcs repeatedly in data-dependent order.
  Memory-mapped CSR for the arc lists is viable (read-only), with the
  flow array in RAM (8 B per arc pair).
- Node capacities (node splitting) double the node count and add n arcs;
  several sources or sinks add one node and O(|S|) arcs.

### Minimum-cost flow (network simplex)

- Primal network simplex, strongly feasible trees, block pivot search
  (block ~ sqrt(m)), parent/child-list tree with subtree re-rooting.
  n = 10^4, m = 8 x 10^4 (max flow of least cost): 47 ms, 6.4 MB; n = 10^3:
  1.9 ms (NetworkX `max_flow_min_cost`: 0.11 s, ~30x).
- Allocations scale with pivots (child lists); a thread/depth array
  representation (Ahuja, Magnanti and Orlin 1993, ch. 11) would remove
  them. Floating-point data: exact for integers below 2^53.

### Gomory-Hu trees (Gusfield 1990)

- n - 1 maximum flows on the original network, reusing one residual
  network (reset is O(m)). n = 3000, m = 12000: 1.7 s, 3 MB; n = 300: 15 ms
  (NetworkX `gomory_hu_tree`: 0.31 s, ~10x before the btower factor).
- Sequential by construction: each flow's terminals depend on the tree so
  far. Parallel options: speculative flows for independent leaves, or the
  recent cut-equivalent tree algorithms with polylogarithmic max-flow
  calls (Abboud et al. 2021, unit capacities). For undirected unit
  capacities, k-edge-components reuse this tree.

### Node and edge connectivity, cuts, disjoint paths

- Even's reduction for kappa (pairs at a least-degree node and among its
  neighbours) and a dominating set for undirected lambda (Esfahanian and
  Hakimi 1984); flows stop at the best value so far. Pairs run in
  parallel, one network copy per thread (so memory grows with threads:
  23.7 MB at 8 threads vs 5.1 MB at 1 for n = 3000, m = 30000).
- kappa, n = 3000, m = 30000: 2.24 s on 1 thread, 0.33 s on 8 (6.8x);
  lambda: 0.96 s / 0.15 s before the dominating set, which cut the number
  of flows from n - 1 to |D| (n = 1000: 0.19 s -> 0.03 s on the VM).
  NetworkX on n = 300, m = 3000: kappa 0.48 s (engine 18 ms, 1 thread),
  lambda 33 ms (engine 2 ms).
- All-pairs node connectivity is n(n - 1)/2 unit flows on the split
  network: n = 300: 2.9 s on 1 thread, 0.41 s on 8 (7x); NetworkX n = 100:
  2.2 s vs 0.09 s. For large n the right tool is a cut-equivalent tree
  for vertex connectivity where one exists, or sampling pairs.
- Global cuts and k-components need cut sets, not just values, so they
  compute full flows for every candidate pair.

### k-components and k-edge-components

- k-components: recursion on minimum node cuts of induced subgraphs, from
  the blocks; each level rebuilds the induced subgraph and its split
  network. n = 300, m = 1200: 0.57 s on 8 threads (NetworkX `k_components`
  on n = 100: 3.8 s vs 0.05 s). Worst case exponential in the number of
  nested cuts, as for Moody and White's procedure.
- k-edge-components: Gomory-Hu tree (undirected) or representative-based
  pairwise flows with cutoff k inside strong components (directed).
  n = 5000, k = 3: 3.3 s. NetworkX's `k_edge_components` on n = 1000
  (average degree 8) fails with RecursionError after 26 s; the engine
  takes 0.12 s.

### Spanning trees and branchings

- Kruskal (sort by (weight, position), union-find): n = 10^6, m = 4 x 10^6:
  0.24 s, 82 MB. Boruvka 0.51 s, 57 MB (rounds parallelizable: the
  per-round lightest-edge scan is a natural parallel reduction with the
  total order guaranteeing the same forest). Prim with a binary heap of
  all arcs: 4.4 s, 447 MB, kept for completeness. NetworkX Kruskal on
  n = 10^4: 0.12 s vs 1.2 ms.
- Optimum branchings (Edmonds/Tarjan with leftist heaps, lazy costs,
  rollback union-find): n = 10^6, m = 8 x 10^6: 2.4 s; 32-bit heap links
  reduced the arena to ~76 B per arc. NetworkX `maximum_branching` on
  n = 10^4: 12 s vs 3.7 ms (~3000x). Heap merges recurse to O(log m)
  depth only.

### DAG kernels

- Topological sort (Kahn by generations), dominators (Cooper-Harvey-
  Kennedy over reverse postorder): n = 10^6, m = 8 x 10^6: 0.35 s and 1.5 s
  (296 MB and 400 MB, mostly deduplicated neighbour rows; the dominator
  rows could reuse the CSR when no input-order is requested). For very
  large flow graphs Lengauer-Tarjan (1979) has better worst-case bounds;
  CHK's iterations stayed small here. NetworkX `immediate_dominators`
  n = 10^4: 0.15 s vs 3.3 ms.
- Transitive closure: one bit row per strong component, unioned in reverse
  topological order (Purdom 1970); output-bound (n = 10^4: 0.18 s, 214 MB
  for ~25 x 10^6 pairs). Transitive reduction reuses the reachability bits:
  n = 10^4: 21 ms. NetworkX on n = 1000: closure 2.2 s vs 1.7 ms,
  reduction 0.14 s vs 0.7 ms.
- Lowest common ancestors in a DAG use ancestor bitsets (n^2/8 bytes);
  beyond ~10^5 nodes, pairs should be answered by targeted searches
  instead (Bender et al. 2005).

### Cycles

- Girth: BFS from every node with pruning at the best cycle so far;
  n = 10^4, m = 3 x 10^4: 0.4 ms (the bound prunes almost everything).
- Bounded simple cycles (DFS above the start node): n = 10^4, m = 3 x 10^4,
  length <= 5: 1.9 s on btower with a strong component search per start
  node (O(nm) in all), 0.30 s on the slower VM after dropping it (it is
  only needed for Johnson's unbounded search); NetworkX on n = 1000:
  11.8 s vs 20 ms.
- Chordless cycles (path extension with neighbour counters): length <= 6,
  n = 10^4: 0.34 s; NetworkX n = 1000: 5.4 s vs 27 ms.
- Minimum cycle basis (de Pina): one Dijkstra per witness edge endpoint
  per basis cycle, in parallel. Restricting starts to the witness's
  endpoints made n = 300 22x faster (4.4 s -> 0.2 s on the VM); n = 1000,
  m = 2000: 2.6 s on 1 thread, 0.61 s on 8. NetworkX `minimum_cycle_basis`
  n = 300: 85 s vs 0.09 s (~900x). Further: Mehlhorn and Michail's
  O(m^2 n / log n) variants and isometric-cycle candidate sets.
