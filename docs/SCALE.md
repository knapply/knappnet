# Analysis beyond memory

Running notes on what knappnet needs to analyse networks larger than RAM,
straight from columnar files and databases, and on algorithms that make
otherwise infeasible analyses practical on one machine. Kept up to date as
the optimization work measures each kernel. Citations are verified before
implementation.

## Goal

Billions of edges on one workstation or one large cloud machine, with the
same validated results as the in-memory engine. A single well-engineered
machine outperforms many distributed graph systems on analytical work
(McSherry, Isard and Murray 2015), so the priority is efficient
single-node execution, not a cluster.

## Storage layers

**Memory-mapped CSR.** A versioned on-disk format (header, offsets,
targets, optional weights and edge ids, little-endian, 64-byte aligned)
that the engine maps read-only. Kernels that already read CSR run
unchanged; the OS pages data in and out. Needs: 64-bit offsets throughout,
no per-edge allocation, and kernels audited for random-access hot spots.

**Compressed graphs.** Gap and variable-length coding of sorted adjacency
lists with reference compression (WebGraph; Boldi and Vigna 2004) brings
web and social graphs to a few bits per edge. Decoding is sequential per
node, which suits traversal and streaming kernels. Node reordering
(degree, BFS, Layered Label Propagation; Boldi et al. 2011; Gorder, Wei et
al. 2016) improves both compression and cache behaviour.

**Arrow and Parquet.** Read edge and node tables zero-copy through the
Arrow C data interface (no Arrow dependency in the engine: the C ABI takes
column pointers). Build CSR from Parquet in a streaming external sort when
the edge list exceeds RAM. Write results (node and edge measures,
partitions) back as Arrow columns.

**Databases.** A DuckDB extension that builds a graph from a query and
exposes measures as table functions; SQL/PGQ (SQL:2023 property graph
queries) as the target query surface. Postgres and SQLite through the same
Arrow path. Graph databases (Neo4j, Kùzu) as import sources via Arrow or
CSV export, not live drivers.

## Execution models

**Semi-external.** Node state in RAM (O(n)), edges streamed from disk in
sequential passes. Fits BFS-family, label propagation, PageRank and other
iterative linear algebra, connected components, k-core peeling, and many
community methods. Abello, Buchsbaum and Westbrook (2002) for the model.

**Edge-centric streaming and grid partitioning.** X-Stream (Roy et al.
2013), GraphChi parallel sliding windows (Kyrola et al. 2012), GridGraph
2-D blocking (Zhu et al. 2015): sequential disk reads, selective block
scheduling to skip inactive regions.

**Sketches and sampling.** HyperBall / HyperANF (Boldi and Vigna 2013;
Boldi, Rosa and Vigna 2011) for distance distributions, closeness and
harmonic centrality on billion-edge graphs; sampled betweenness with error
guarantees (Riondato and Kornaropoulos 2016; KADABRA, Borassi and Natale
2019); triangle and motif estimators by edge sampling; count-min and
HyperLogLog neighbourhood sketches.

**Streaming and incremental updates.** Dynamic connectivity, incremental
k-core, and incremental PageRank for graphs that change while being
analysed.

## Kernel ledger

One line per kernel: passes over edges, working set, semi-external or
not, best scalable variant. Filled in from the optimization and community
workstreams (docs/drafts/scale_notes_*.md until merged).

| Kernel | Passes | Working set | Semi-external | Scalable variant |
|---|---|---|---|---|
| | | | | |

## Statistical modelling at scale

Dyad design matrices streamed in blocks with case-control subsampling of
non-edges; ERGM estimation by equilibrium expectation on large networks
(Byshkin et al. 2018); relational event models with sampled risk sets
(Lerner and Lomi 2020); latent space models with case-control likelihoods
(Raftery et al. 2012). These are the routes to fitting models people
currently do not attempt on networks of realistic size.

## Open questions

- On-disk format: own format, or adopt an existing one (WebGraph, Kùzu
  storage, GraphAr on Parquet) for interchange.
- Whether DuckDB integration ships as a separate extension crate.
- Minimum disk throughput assumptions for the streaming kernels.
