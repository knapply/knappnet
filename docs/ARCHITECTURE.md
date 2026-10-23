# Architecture

This document records how the project is built and the rules that keep it
correct, fast, and maintainable as it grows.

## Layers

```
crates/engine     Rust. Graph storage and algorithms. No dependencies.
                  Published to crates.io under the project name.
crates/capi       Rust. A small C ABI over the engine (crates/capi/include/netcapi.h).
r/                The R package. C glue (r/src/*.c) calls the C ABI.
validation/       Fixtures, golden results from reference tools, comparison
                  (validation/README.md).
bench/            Benchmarks against other tools (bench.yml; the quick tier
                  runs in ci.yml).
data-raw/         Scripts that rebuild the bundled example data from its
                  published sources (data-raw.yml).
tools/            Maintenance scripts, listed in CONTRIBUTING.md.
docs/             Architecture notes, the release checklist, and the
                  bibliography with its generated index.
```

Inside the engine, shared machinery lives in its own modules: `graph`
(CSR storage and derived views such as `Graph::simple_undirected`),
`matrix` (`Dense`), `linalg` (all dense linear algebra, one implementation
per operation), `union_find`, `rng`, `parallel`, `weights` (path lengths
and weight checks), and `conventions`. Measures live in `metrics`,
statistical tests in `qap`, and dynamics, space,
drawing, and generation in `temporal`, `spatial`, `layout`, and
`generate`. The crate documentation (`lib.rs`) and `metrics/mod.rs` carry
the full map. The C ABI header `netcapi.h` is divided into commented
sections by area after a preamble on return codes, ownership, and buffer
sizes, and the R glue shares its helpers through `r/src/glue.h`.

Each area of functionality (descriptive statistics, positions, spectral
measures, ...) is one module per layer with the same name: an engine
module, a C ABI module, a C glue file, an R file, a test file, a manifest
file, a runner, and oracle scripts. `tools/new_measure.py` creates them.

The R package (`r/`) keeps user-facing functions in `r/R/<area>.R`,
shared helpers in `r/R/utils-*.R` (argument checks, the engine handle and
error rewriting, node lookups), and plotting in `r/R/plot-*.R`. Help pages
are hand-written Rd files in `r/man/`, the only source of documentation
(there are no roxygen comments); `tools/check_r_docs.py` keeps them, the
function signatures, `NAMESPACE`, and the pkgdown index consistent.

Language bindings never call the engine directly; they go through the C ABI.
Python bindings will use the same ABI, so every language shares one tested
surface.

### Why a C ABI instead of a binding generator

Linking a Rust static library into an R package through a thin C layer
(`ooms_hellorust`) needs nothing on the R side beyond R's own headers. With
no external crates, CRAN builds are fully offline and there is nothing to
vendor. The R package carries a copy of the engine under `r/src/rust/`,
maintained by `tools/sync_r_sources.py`; CI fails if the copy drifts.

## Data model

A network is a pair of tables:

- **nodes**: one row per node, with an identifier column and any attributes.
- **edges**: one row per edge, with two endpoint columns and any attributes.

Column roles record which columns carry structural meaning. The engine's
`frame` module (`crates/engine/src/frame/`) defines them and their
meaning, with sources, in its module documentation; in short:

| Role | Table | Meaning |
|---|---|---|
| `id`, `source`, `target` | nodes, edges | Identifiers: integer, double, or text; the endpoint columns share one type and the identifier column matches it. |
| `weight` + weight type | edges | A numeric column whose meaning is declared once: `strength` or `distance`. No default. |
| `mode` | nodes | Any discrete column; with `between_modes_only`, a two-mode or k-partite network. |
| `relation` | edges | The relation of a multiplex network; one long edge table. |
| `sign` | edges | `-1`/`1` or logical; the weight is then a non-negative magnitude. Provisional. |
| `start`, `end`, `wave` | edges | Half-open spells `[start, end)`; `start` alone is a contact; wave `k` is `[k, k + 1)`. |
| `node_start`, `node_end` | nodes | Presence spells; an edge is active only while both endpoints are present. |
| `x`, `y` + geometry flags | nodes | Representative coordinates; geometry and its CRS stay with the binding, which reports the CRS of each table so they can be checked for agreement. |
| `node_missing`, `edge_missing` | nodes, edges | Logical marks for non-respondents and unobserved dyads. Provisional. |

Declared properties (`directed`, and optionally `loops`, `multiple`,
`between_modes_only`, `add_nodes`) are enforced: violating rows are
errors, contradictory declarations are errors, and an undeclared property
is observed and reported as a note. Two edges are parallel when they have
the same endpoints and the same relation and spell, so a multiplex or
temporal network with one edge per relation and spell is simple.

Attributes stay in the data frames. The engine receives the role columns
as typed buffers (integer, double, text, logical, each with missing
values), validates them, and builds the CSR graph; results come back as
vectors aligned with the node table.

### Validation and construction

Every binding builds its network object through the frame, so the checks
and their messages are shared. `frame::validate` runs every check in
O(n + m) time (identifiers are hashed once; parallel edges are found by
sorting the edge key; column scans run in parallel chunks merged in row
order, so the result does not depend on the thread count) and returns a
report rather than stopping at the first problem. Each diagnostic has a
stable code (`E_UNKNOWN_NODE`, `E_WEIGHT_TYPE_MISSING`, `N_LOOPS`, ...;
the full list is on `frame::Code`), a severity (errors stop construction,
notes are advisory), the table it concerns, a message, and the first 50
offending 0-based rows with the total count, so a binding phrases row
lists its own way. The checks fall into four groups, modelled on the
validation of `as.network.data.frame()` in the network package
(`butts2008network`): structure (endpoint and identifier types, missing,
duplicate, and unknown identifiers), declared properties, attribute
sanity (weights, modes, relations, signs, times, coordinates, geometry),
and advisory notes (nodes added from edges, isolates, 0-based-looking
identifiers, coerced identifiers, loops and parallel edges when nothing
was declared).

`frame::build` validates and builds the graph: node `i` is row `i` of the
node table, followed by nodes added from the edges in order of first
appearance; edge `e` is row `e` of the edge table; the mapping is
returned. `frame::multiplicity` implements the parallel-edge rule (edge
multiplicity, the simplicity test, and `simplify`, which returns an index
mapping plus combined weights so bindings aggregate other attributes
themselves), and `frame::summary` computes everything a printed summary
shows in one call. The C ABI section in `crates/capi/include/frame.h`
exposes all of it over typed column buffers.

## Files and interoperability

File readers and writers live in the engine (`crates/engine/src/io/`), so
every language binding shares one implementation. Each format reads into,
and writes from, the same neutral model: a node table starting with `id`, an
edge table starting with `from` and `to`, typed attribute columns
(logical, integer, double, text, each with missing values), and
network-level attributes. The XML and JSON parsers are small, dependency-free,
and never resolve external entities.

Conversions to and from other R packages (igraph, tidygraph, network) live in
the R package, since those objects only exist in R. Formats are validated in
both directions against NetworkX (`validation/io_check.py`) and igraph
(`r/tests/testthat/test-io.R`).

## Storage

Graphs are stored in compressed sparse row (CSR) form. Outgoing arcs of node
`v` are `targets[offsets[v]..offsets[v+1]]`, sorted by target, with weights in
a parallel array. Directed graphs also store incoming arcs in the same
layout. Undirected graphs store each edge once per endpoint. Rows are built
with two stable counting sorts in O(n + m).

Node indices are 32-bit, which supports roughly four billion nodes. Offsets
are machine words. A 64-bit index option is planned for very large graphs.

## Weighted and binary networks

Whether a network is weighted is a property of the network (`weight =` in
`netframe()`), not something each function guesses. Measures follow it by
default and can be overridden per call.

- **Binary networks** (an edge exists or it does not) never touch weight
  arrays. Path-based measures use breadth-first search with exact integer
  hop counts. Dense binary networks use a bit matrix (`engine::bitset`) for
  constant-time arc tests when it is no larger than adjacency lists.
- **Weight meaning** is declared once: `weight_type = "distance"` (cost,
  length, latency; the igraph and NetworkX convention) or `"strength"`
  (contact frequency, tie intensity), which path-based measures convert to
  lengths 1/w (`newman2001scientific`, `opsahl2010node`). Additive measures
  such as strength (`barrat2004architecture`) use weights as given.
- **Algorithm choice** follows the weights (`engine::weights`): uniform
  weights reduce exactly to the binary case; positive integer lengths up to
  1024 use Dial's bucket queue with exact integer distances
  (`dial1969algorithm`); anything else uses Dijkstra's algorithm
  (`dijkstra1959note`) with a relative tolerance for ties.
- Measures defined only on binary ties (the triad census) ignore weights and
  say so in their documentation.

Property tests check that every algorithm returns the same results on the
same lengths, that uniform weights reproduce unweighted results, and that
strengths w give the same results as distances 1/w. The validation corpus
includes uniform-weight, integer-weight, and real-weight networks so each
path is exercised against the oracles.

## Conventions

Every algorithm documents its conventions in its source file and in the R
help page: treatment of self-loops and parallel edges, whether weights are
distances or strengths, normalization, and handling of disconnected graphs.
Defaults follow igraph and NetworkX where they agree. Where reference tools
disagree, the manifest records the difference and the reason.

## Parallelism

`engine::parallel` assigns work items to threads in a fixed pattern and
merges results in thread order, so results are reproducible for a given
thread count. Integer results are identical for any thread count;
floating-point results can differ in the last bits across thread counts.
Threads come from the standard library, which avoids OpenMP (unavailable on
some CRAN build machines).

## Validation

`validation/manifest.d/<area>.yaml` lists every measure, its engine call,
its tolerance, and its equivalent in each reference tool, with every known
convention difference and its reason. Oracles:

| Oracle | Where it runs | Notes |
|---|---|---|
| definition | locally and CI | Direct evaluation of the definition with exact arithmetic on small fixtures (`validation/oracles/<area>_definition_oracle.py`). Shares no code with anything else. |
| networkx | locally and CI | `validation/oracles/<area>_networkx_oracle.py` |
| igraph, sna, and other R packages | CI | `validation/oracles/r_<area>.R` |
| UCINET | manual | Golden files produced by hand on Windows (`validation/UCINET.md`) and committed. |

`validation/compare.py` classifies each comparison as a match, a documented
convention difference, or a failure, and writes `validation/REPORT.md`.
Golden files are committed along with the oracle versions that produced
them, and are regenerated when an oracle version changes. Checks that do
not fit the manifest (file formats, QAP, random graphs, layouts) have their
own scripts, and `validation/scale/` repeats the comparison on larger
networks. `validation/README.md` has the details.

Engine invariants (relabeling nodes permutes results, unit weights match
unweighted results, thread count does not change results, census totals)
are checked with randomized tests in `crates/engine/tests/*_properties.rs`.

## Rules

1. No measure is exported without a manifest entry, at least one passing
   oracle, and a citation in `docs/references.bib` (`tools/check_bib.py`
   checks the bibliography and generates `docs/references.md` from it).
2. Every algorithm cites its source in a comment at the top of its file and
   is implemented from the published source, never from the code of a
   copyleft reference tool (see CONTRIBUTING.md, "Licensing").
3. The engine stays dependency-free unless a dependency is clearly worth its
   cost for CRAN builds.
4. The project name never appears in identifiers, C symbols, or R function
   names; it lives in the package metadata (including the published crate
   name, which workspace code imports as `engine`) and the docs.
5. Approximate algorithms are opt-in, named as such, and report their error
   guarantees.

## Designing for accelerators

GPU backends are out of scope for now, but these choices keep the door open:

- Heavy operations sit behind small, batch-friendly functions over flat
  arrays, so a device implementation can replace them without restructuring.
- Storage is flat and contiguous (CSR, parallel weight arrays), which copies
  directly to device memory.
- Inherently sequential algorithms (for example ERGM edge toggles) are kept
  separate from the bulk kernels so they can stay on the CPU.
- Randomized algorithms will use counter-based generators, which give
  reproducible parallel streams on any backend.
- The validation harness runs per backend, so a future accelerator path is
  checked against the same golden files.

A future backend layer will dispatch between implementations based on graph
size and density. Sparse linear-algebra formulations of graph algorithms
(`kepner2016mathematical`, `davis2019algorithm`) and direction-optimizing
traversal (`beamer2012direction`) are the main influences.
