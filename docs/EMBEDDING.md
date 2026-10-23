# The engine as a component

The engine is meant to sit underneath anything that analyses networks:
language packages, databases, desktop and browser tools, services, and
embedded devices. These are the properties it must have, with status.

## Interface stability

- [ ] `net_abi_version()` and a header version macro; semantic versioning
      of the C ABI with a written compatibility policy.
- [ ] Every function returns a status code and writes a message into a
      caller buffer (done for most); a single error-code enum documented
      in the header.
- [ ] A frozen, documented subset ("stable") distinct from experimental
      entry points.
- [ ] A Rust API for direct use (`engine` crate) with the same stability
      tiers.

## Runtime behaviour

- [x] No global mutable state; handles are independent and thread-safe to
      read concurrently.
- [x] Caller-allocated output buffers; deterministic results for a given
      seed, independent of thread count.
- [ ] Cancellation and progress: an optional callback polled by long
      kernels (betweenness, community detection, layouts, simulation).
- [ ] Caller-supplied thread pool or thread budget per call (done) and a
      way to cap memory.
- [ ] Optional allocator hooks for hosts that manage memory (databases,
      game engines, R/Python arenas).

## Data interchange

- [ ] Arrow C data interface in and out (zero-copy columns), so any
      Arrow-speaking host (Python, R, DuckDB, Polars, Spark, JavaScript)
      shares buffers without conversion.
- [ ] A storage abstraction so kernels run on in-memory CSR, memory-mapped
      or compressed graphs, and database-backed adjacency without code
      changes (see SCALE.md).
- [ ] Versioned on-disk graph format.

## Extensibility

- [ ] User-supplied edge weights and costs via callbacks or columns, and
      generic semirings for path and propagation kernels.
- [ ] Custom statistics for models and null-model tests registered
      through the ABI.

## Targets

- [x] Linux, macOS, Windows; MSRV 1.70; std only.
- [ ] WebAssembly build (browser and server runtimes).
- [ ] Bindings: R (done), Python, JavaScript/WASM, then Julia, Java/JVM,
      .NET, Go as demand appears; a DuckDB extension.
- [ ] Small static library and shared library builds with symbol
      visibility limited to `net_*`.
