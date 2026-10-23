# Scale notes: community, position, layout, permutation, generator and dense kernels

Working notes for docs/SCALE.md (kernel ledger). Scope: Louvain, fast
greedy, leading eigenvector, Girvan-Newman, maximal cliques, factions and
core-periphery, structural equivalence and CONCOR, REGE (CATREGE), stress
layouts, QAP/MRQAP, random graphs, dense linear algebra. Path and
traversal kernels are in the paths workstream's notes.

Measurements: `cargo run --release --example bench_models` on a 2-vCPU
Xeon (2.1 GHz, AVX-512 capable), 7 GB RAM, one thread unless stated; heap
is the peak of live heap bytes from a counting allocator. Sparse = average
degree 8, medium = density 2%, dense = density 30%. Timings on this
virtual machine vary by up to 2x between runs; ratios below come from runs
of both versions back to back.

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
  98k to 208 after the aggregation rewrite.
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
- Measured: sparse n = 10k, m = 40k: 4.46 s -> 0.086 s (52x), peak heap
  169 MB -> 4.4 MB; medium n = 2k, m = 40k: 2.61 s -> 0.042 s; dense n =
  600, m = 54k: 0.166 s -> 0.056 s.
- Working set: O(m) maps + O(m) queue. ~110 B/edge measured.
- Semi-external: no (global queue, arbitrary merges). Out-of-core variants
  are not practical; at scale use Louvain/Leiden.

### Leading eigenvector (Newman 2006)

- Access: dense n x n modularity matrix, dense symmetric eigensolver
  (EISPACK tred2/tql2, O(n^3), O(n^2) memory); limited to 3000 nodes.
- Change: the solver now works on the transpose of the EISPACK working
  matrix, so the Householder and rotation loops run over contiguous rows
  (same operations, same order, bit-identical). n = 1000 decomposition:
  5.0 s -> 1.5 s; leading eigenvector n = 1000: 7-8 s -> 2.0 s.
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
  sparse n = 300: 6.6 s -> 5.2 s in back-to-back runs (allocations 1.47M
  -> 43k with the earlier queue reuse).
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
- Measured (depth 3, directed): dense n = 1k, m = 300k: 29.5 s -> 1.5-2.5 s;
  initial correlations 75% of the old time, now a few percent. Sparse
  n = 2k, m = 16k: 149 s -> 11-14 s. AVX2 copy of the tile alone: 2x
  (20.4 s -> 10.9 s back to back).
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
  majorization passes over all pairs. Change: each pass processes four
  rows at once in independent lanes (AVX2 copy at run time); per-row sums
  in order, identical results.
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
  1000 permutations: DSP 2.19 s -> 0.78 s; Y-permutation 0.54 s -> 0.17 s.
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
  0.42 s -> 0.25-0.30 s. Eigendecomposition as above. Cholesky unchanged
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
