# Scale notes: graph storage and traversal kernels

Working notes from the traversal/path optimization pass, for docs/SCALE.md.
Measurements are from `crates/engine/examples/bench_paths.rs` (counting
allocator, `VmHWM`) on a 2-core x86-64 machine; n nodes, m edges, arcs = m
(directed) or 2m (undirected).

## Storage (`graph.rs`)

- Layout today: edge list (`sources`, `targets`: u32 each, 8 B/edge), plus
  out-CSR (`offsets`: usize, 8 B/node; `targets`: u32, 4 B/arc), plus an
  in-CSR for directed graphs, plus f64 weights per edge and per arc when
  weighted. Unweighted undirected: about 16 B/edge + 8 B/node resident;
  weighted adds 24 B/edge (one copy per edge, two per arc).
- Construction (`from_edges`): two stable counting sorts; measured peak
  heap 2.6 MB for m = 50k undirected (about 52 B/edge transient). Two
  sequential passes over the edges plus one scattered write per arc.
- Memory-mapped CSR: feasible as is. Every kernel reads `offsets` and
  `targets` only through `Csr::neighbors`/`arc_range`; a read-only mmap of
  the two arrays (with u64 offsets) would serve every read-only measure.
  Dropping the edge list once the CSR exists would save 8 B/edge but the
  edge list is needed for edge-indexed outputs (edge betweenness, bridges).
- u32 offsets when arcs < 2^32 would save 4 B/node; negligible next to
  arcs (rejected for now: n is small next to m in every realistic case).
- Compression: neighbor lists are sorted, so gap + varint coding
  (Boldi and Vigna 2004, WebGraph) would cut `targets` to about 1-2 B/arc
  on locality-ordered graphs; decoding is sequential per list, which fits
  BFS/MS-BFS (they scan whole lists) but not `BitMatrix`-style membership
  tests. Needs a node ordering with locality first (BFS order, Cuthill and
  McKee 1969; layered label propagation, Boldi et al. 2011; Gorder, Wei et
  al. 2016) to pay off.

## Breadth-first search, one source (`descriptive`, `ego`, `resilience`)

- Access: queue scanned sequentially; per arc one random read of the
  target's state (4-8 B). Working set: O(n) state + the CSR.
- Edge streaming: level-synchronous BFS is semi-external: O(n) state in RAM,
  one sequential pass over all edges per level (diameter D passes). Good
  for low-diameter graphs (social, web: D < 30), bad for road networks.
- Direction-optimizing BFS (Beamer, Asanovic and Patterson 2012) reduces
  arcs examined by 5-20x on low-diameter graphs; needs the in-CSR (already
  built for directed graphs) and a bitmap of visited nodes (n/8 bytes).

## Multi-source BFS (closeness, harmonic, eccentricity, mean distance, reach, neighborhood sizes)

- Now runs 64 sources per pass with one u64 per node (`bitset::MultiBfs`,
  Then et al. 2014). State: 3 x 8 B/node per worker; arcs scanned once per
  level per batch instead of once per source, so all-pairs cost drops from
  n searches to n/64 batches of word-wide ORs.
- Access: per level a sequential scan of the active list and random ORs
  into `next[w]`; then either a sorted touched list (small frontiers) or a
  sequential scan of all n words (large frontiers).
- Semi-external form is natural: each level is one pass over the edge
  stream (or only the edges of active nodes, with an index), with O(n)
  words in RAM per 64 sources. With W-word lanes (64W sources) per pass,
  all-pairs needs n / 64W x D passes over the edges: for n = 10^8, W = 8,
  D = 20, about 4 x 10^6 passes, still infeasible; for the source-sampled
  estimators (Eppstein and Wang 2004; Brandes and Pich 2007) with k = 512
  sources it is one batch: D passes over the edges.
- Direction-optimizing levels (Beamer, Asanovic and Patterson 2012, applied
  to source sets): a level whose frontier holds more than 1/8 of the arcs
  is expanded bottom-up, a sequential sweep over all nodes that reads each
  unfinished node's in-list and stops early once all 64 sources have
  reached it. Measured (1 thread, CPU s, identical results): closeness
  dense1k 0.0088 -> 0.0019, dense2k 0.043 -> 0.008, medium3k 0.060 ->
  0.022, sparse10k 0.37 -> 0.29; eccentricity similar. A bottom-up level
  is one sequential pass over the in-CSR, the streaming-friendly pattern;
  top-down levels are random-access. Cost: needs the in-CSR (directed
  graphs already keep one: +4 B/arc +8 B/node).
- Direction views: `In` and `All` searches on directed graphs read the out-
  and in-CSR in place (`bitset::View`) instead of building a reversed or
  symmetrized copy (about 50 B/edge transient: edge list + two counting
  sorts). Neighborhood sizes, order 2, mode all: densedir1k (300k arcs)
  0.0110 -> 0.0030 s, peak heap 18.3 MB -> 55 KB; meddir3k 11.1 MB -> 186 KB.
- When batching pays: a 64-source batch saves work only where the searches
  overlap. Depth-limited searches (neighborhood sizes, reach k) on sparse
  graphs have balls of about d^r nodes; when 64 d^r < n the batch only adds
  3 words/node of state traffic, so one plain search per source
  (`bitset::SingleBfs`, 4 B/node state) is used instead
  (`bitset::batches_pay`). Measured, neighborhood sizes order 2: sparse100k
  0.47 s batched -> 0.10 s per source (main: 0.11 s).
- Blocking: partitioning nodes into cache-sized ranges and bucketing the
  pushes by target range (X-Stream, Roy, Mihailovic and Zwaenepoel 2013;
  GridGraph, Zhu, Han and Chen 2015) turns the random ORs into sequential
  writes per block.

## Dijkstra and Dial (weighted closeness, eccentricity, betweenness, reach)

- Dijkstra now uses an indexed 4-ary heap (Johnson 1975): one entry per
  node instead of one per relaxation; heap memory O(n) per worker (12 B per
  entry plus 4 B/node position), versus O(m) worst case for lazy deletion.
- Access: heap operations are random within O(frontier) entries; arcs are
  read sequentially per settled node; per relaxed arc one random read of
  the target's distance.
- Dial: bucket array of width max-length + 1; O(n) nodes in buckets.
- Not streamable in label-setting form; semi-external variants exist
  (delta-stepping, Meyer and Sanders 2003, processes buckets of width delta
  with one edge pass each, the same pattern as level-synchronous BFS).
- Edge weights add 8 B/arc to the working set; f32 lengths (or u8/u16
  integer lengths for Dial, which only needs lengths up to `DIAL_MAX`)
  would halve or quarter that.

## Brandes betweenness (node and edge)

- Per source: forward search (as above) plus a reverse pass over the
  settled nodes in reverse order, re-reading each node's arcs (two passes
  over the reachable arcs per source). State per worker: sigma, coefficient
  (f64), hops (u32) or distance, rank: about 24-32 B/node, plus one f64
  accumulator per node (node) or per edge (edge betweenness).
- Measured: 586 KB peak heap for n = 10k (sparse10k), about 59 B/node; the
  per-thread accumulators make it O(threads x n).
- Reduction: sources split into 64 fixed groups, each summed into its own
  accumulator and merged in group order (`parallel::sum_in_order`), so
  results are thread-count independent. Live accumulators: at most about
  2T + 1 (a thread more than 2T groups ahead waits), each n f64 (node) or
  m f64 (edge betweenness: 8 B/edge each, so 3T+ copies of an m-vector is
  the dominant memory term for edge betweenness on large m). Workers no
  longer hold their own result vector.
- Semi-external: the reverse pass needs the settle order (4 B/node) and the
  forward level structure; per source two passes over the edge stream per
  level. Exact all-pairs is O(nm) regardless of storage, so the feasible
  route at billions of edges is source sampling (Brandes and Pich 2007;
  Riondato and Kornaropoulos 2016) with MS-BFS-style batching of the
  forward phase.

## Triangles and clustering (`descriptive`)

- Access: per node, sorted-list intersections (or `BitMatrix` probes for
  dense graphs). Working set: O(m) for the oriented adjacency.
- Degree ordering / forward algorithm (Schank and Wagner 2005; Latapy
  2008) bounds work by O(m^1.5) and halves the lists; with lists in
  degree order this is the masked SpGEMM formulation (Azad, Buluc and
  Gilbert 2015; Wolf et al. 2017).
- External memory: triangle listing needs O(m^1.5 / (M^0.5 B)) I/Os
  (Pagh and Silvestri 2014); simple chunked forward algorithm with node
  ranges that fit in RAM is the practical route.

## Linear-time kernels (degree, coreness, PageRank, eigenvector, Katz, HITS)

- One sequential pass over the CSR per iteration plus a random read of the
  source's value per arc (pull form). Fully streamable: O(n) values in RAM,
  edges on disk, one pass per iteration (GraphChi, Kyrola, Blelloch and
  Guestrin 2012). Hilbert or degree ordering of the edge stream improves
  the cache behaviour of the random reads (McSherry, Isard and Murray 2015,
  "Scalability! But at what COST?").
