# Contributing

Accuracy is the first priority of this project: every measure must be
correct, cited, documented, and validated against reference tools before it
ships. This page is the recipe for adding one, the rules that keep the
library consistent, and the checks that CI runs. Read
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) first for how the layers fit
together, and [validation/README.md](validation/README.md) for the
validation framework.

## Setting up

- Rust (stable; the engine supports Rust 1.70 and newer and has no
  dependencies, so `cargo` works offline).
- Python 3.10 or newer with `networkx`, `numpy`, `scipy`, `pandas`, and
  `pyyaml` for the validation and the Python benchmarks.
- R 4.1 or newer with the packages in `r/DESCRIPTION` (`Suggests` holds the
  reference tools) for the R package, its tests, and the R oracles. A C
  compiler is needed to build the R package; the Rust engine is compiled by
  R's build through `r/src/Makevars`.

## Adding a measure

`python3 tools/new_measure.py <area> <name>` writes a skeleton for every
step below (or prints the snippet to paste when the area's files exist),
marks every place that needs the real measure with `TODO`, and lists the
checks to run. `--dry-run` only prints. The steps, in order:

1. **Engine function** in `crates/engine/src/metrics/<area>.rs` (or the
   matching module for layouts, files, models). The module comment cites
   the sources of the definition and the algorithm by BibTeX key and lists
   every convention (see "Conventions" below). Unit tests with hand-checked
   values go in the module. Use `crate::parallel` for per-node work,
   `crate::rng::Stream` for randomness, `crate::matrix::Dense` for dense
   matrices, `crate::linalg` for every dense linear-algebra operation
   (Cholesky and SPD solves, inverses, pseudo-inverses, symmetric
   eigendecomposition; add new operations there, never in a measure's
   module), `crate::union_find::UnionFind` for disjoint sets,
   `Graph::simple_undirected` for the simple undirected neighbour lists,
   `crate::weights` for path lengths and `weights::additive_weights` for
   weights that are added as given. `crates/engine/src/lib.rs` and
   `metrics/mod.rs` carry a map of where each kind of code lives.
2. **Property tests** in `crates/engine/tests/<area>_properties.rs`:
   randomized invariants (relabeling nodes permutes results, the thread
   count does not change them, unit weights reproduce the unweighted
   result, totals and symmetries the measure must satisfy), plus an
   independent brute-force evaluation on small random networks where one
   is easy to write.
3. **C ABI** function `net_<name>` in `crates/capi/src/<area>.rs`, declared
   in `crates/capi/include/netcapi.h` with a comment. It is
   `#[no_mangle] pub unsafe extern "C"`, has a `# Safety` section, returns
   0 on success, and writes its error message through `guard`. Options are
   small integer codes documented at the top of the module.
4. **C glue** `C_<name>` in `r/src/<area>.c`. Include `"glue.h"`, which
   brings in `netcapi.h` and the shared helpers (`graph_from_ptr()` for the
   graph handle, `as_count()`, `as_seed()`, `as_int64()`, `check_status()`
   to turn an engine error into an R error after unprotecting,
   `vector_result()` and `scalar_result()`); raise R errors only after the
   engine call has returned. Register it in `r/src/init.c`: an `extern`
   line and a `R_CallMethodDef` row in the group for its file, sorted by
   name. `netcapi.h` is organized in commented sections; add the
   declaration to its area's section.
5. **R function** in `r/R/<area>.R`, its help page `r/man/<name>.Rd`,
   tests in `r/tests/testthat/test-<area>.R`, `export(<name>)` in
   `r/NAMESPACE`, and the help topic in the reference index of
   `r/_pkgdown.yml`. The R package has no roxygen comments: the Rd files
   are the single, hand-written source of documentation, each function
   carries a one-line comment naming its Rd file, and
   `tools/check_r_docs.py` checks that every `\usage` matches the
   function's formals, every export is defined and documented, every
   topic is in the [pkgdown](https://pkgdown.r-lib.org/) index, and Rd sections come in the package's
   order. Shared helpers live in `r/R/utils-*.R`: `check_netframe()`,
   `engine_graph()`, `engine_call()`, `resolve_threads()`,
   `resolve_weighting()`, and `resolve_additive_weights()` in
   `utils-engine.R`; `check_flag()`, `check_choice()`, `check_count()`, and
   the other argument checks in `utils-checks.R`; node lookups and
   attribute columns in `utils-nodes.R`. Plotting lives in `plot-*.R`.
   Tests use `skip_if_not_installed()` for every suggested package and at
   most 2 threads.
6. **Validation**: a manifest entry in `validation/manifest.d/<area>.yaml`,
   a runner in `validation/bindings/<area>.py`, a definition oracle in
   `validation/oracles/<area>_definition_oracle.py`, and an oracle for every
   reference tool that implements the measure ([NetworkX](https://networkx.org/) in Python; [igraph](https://r.igraph.org/),
   [sna](https://cran.r-project.org/package=sna), and others in `validation/oracles/r_<area>.R`). Regenerate the
   Python golden files and commit them.
7. **Bibliography**: the BibTeX entries in `docs/references.bib`, verified
   against Crossref (`https://api.crossref.org/works/<DOI>`) or the
   publisher, then `python3 tools/check_bib.py --update`, which sorts the
   file and regenerates `docs/references.md`.
8. **R sources**: `python3 tools/sync_r_sources.py` copies the engine into
   `r/src/rust/`; CI fails when the copy is stale.
9. **Checks**: everything under "Checks" below. If the README shows the new
   measure, run `python3 tools/link_tools.py README.md` and
   `Rscript tools/check_readme.R`.

A measure is *validated* when its manifest entry has at least one passing
oracle on every fixture it applies to and no failure; a definition oracle
(a slow, obviously correct implementation of the textbook formula sharing
no code with the engine) is expected wherever one is practical.

## Conventions

- **Document every convention** at the top of the engine module and in the
  R help page, in this order: self-loops, parallel edges, direction,
  weights, normalization, disconnected networks and undefined values, and
  which reference tools the result matches.
- **Weighted versus binary is explicit.** Path-based measures use
  `crate::weights::{path_lengths, Lengths, WeightKind, PathAlgorithm}`, so
  binary networks use breadth-first search, small integer lengths Dial's
  bucket queue, and the rest Dijkstra's algorithm; `weights_as` says
  whether weights are distances or strengths. Measures defined only on
  binary ties ignore weights and say so. Additive measures use weights as
  given. In R, `weighted = NULL` follows the network's weight role.
- **`convention` argument** where the reference tools disagree (`"igraph"`,
  `"networkx"`, `"sna"`, `"ucinet"`): the default follows [igraph](https://r.igraph.org/), the
  differences are tabulated in `crates/engine/src/conventions.rs` and in
  `?conventions`, and the help page's examples show one call per
  convention on a small network where the results visibly differ, with a
  comment naming the tool each reproduces. The same goes for any other
  argument that exists to match a particular tool (`mode` variants,
  `scale`, `unreachable` rules, E-I `scope`).
- **Threads**: every per-node computation takes `threads` (0 means all
  cores; `resolve_threads()` applies CRAN's limit). Results are
  reproducible for a given thread count; integer results for any count.
  Examples and tests use at most 2 threads.
- **Seeds**: randomized methods take a `seed` and draw from
  `crate::rng::Stream`, so results are identical on every machine and
  thread count. In R, a missing seed is drawn from R's generator so
  `set.seed()` works.
- **Undefined values** are `NaN`, never an error, and the help page says
  when they occur. Each validation area states how its runner and oracles
  encode them.
- **Errors**: engine errors are `Error::InvalidArgument` with a message a
  user can act on; `engine_call()` prefixes the R function name and
  rewrites edge positions as rows of `x$edges`. Invalid R arguments are
  checked in R (`check_flag()`, `check_choice()`) before the engine runs.
- **Results** come back aligned with `x$nodes` (node measures), as a named
  vector (census-like results), or as a data frame.
- **Style**: `cargo fmt` and `cargo clippy -D warnings`; base R only in
  package code; comments explain why and cite sources rather than narrate
  the code. The project name appears only in package metadata and docs,
  never in identifiers, C symbols, or R function names.

## Licensing

The project is proprietary (all rights reserved; see LICENSE), and
outside contributions are not accepted at present. Implement every method from its published source
(papers, books, official documentation). Reference tools under the GPL or
another copyleft license ([igraph](https://r.igraph.org/), [sna](https://cran.r-project.org/package=sna), network, [statnet](https://statnet.org/), [ergm](https://cran.r-project.org/package=ergm), [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/),
tsna, xucinet, and others) are used only as black-box oracles: run them,
never read, copy, or translate their source. Permissively licensed
projects ([NetworkX](https://networkx.org/), BSD; [NetworKit](https://networkit.github.io/), MIT) may be read but not copied; cite
them when they informed a design.

## Checks

CI runs all of these on every push (see `.github/workflows/ci.yml`). Run
them before opening a pull request:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc -p knappnet --no-deps
./tools/check_c_abi.sh                      # netcapi.h matches the ABI
python3 tools/check_bib.py                  # bibliography and docs/references.md
python3 tools/check_r_functions.py          # no R function defined twice
python3 tools/check_r_docs.py               # Rd usage, NAMESPACE, pkgdown index
python3 tools/sync_r_sources.py --check     # r/src/rust/ is current
python3 tools/new_measure.py probe probe_measure --dry-run
python3 validation/check_manifest.py        # manifest entries well-formed
python3 validation/fixtures.py              # fixtures are reproducible
python3 validation/run_oracles.py           # --python-only without R
python3 validation/compare.py               # 0 FAIL
python3 validation/io_check.py              # file formats against NetworkX
python3 validation/qap_check.py             # QAP against exact enumeration
python3 validation/cug_check.py             # random graphs and CUG tests
python3 validation/layout_check.py          # layout properties
Rscript tools/check_readme.R                # README output matches R
R CMD build r && R CMD check --as-cran knappnet_*.tar.gz
```

## Tools

| Tool | Purpose |
|---|---|
| `tools/map.py` | Cross-layer index: for each exported R function, its R file, help page, C glue, C ABI functions, engine calls, and validation variants (`--find NAME` searches every layer; writes the ignored `docs/map.json`). |
| `tools/new_measure.py <area> <name>` | Scaffold a measure through every layer; `--dry-run` prints only. |
| `tools/check_bib.py` | Check `docs/references.bib` (unique sorted keys, fields, DOIs, cited keys exist) and that `docs/references.md` is current; `--update` sorts and regenerates. |
| `tools/check_c_abi.sh` | Compile and run a C program against `netcapi.h`. |
| `tools/check_r_functions.py` | Fail if two files in `r/R` define the same function. |
| `tools/check_r_docs.py` | Static checks of the R package without R: balanced brackets, `\usage` against function formals, NAMESPACE against definitions and Rd aliases, [pkgdown](https://pkgdown.r-lib.org/) index completeness, Rd section order. |
| `tools/check_readme.R` | Run the README's R chunks and compare their printed output; `--update` rewrites it. |
| `tools/sync_r_sources.py` | Copy the engine into `r/src/rust/`; `--check` for CI. |
| `tools/link_tools.py <files>` | Link the first mention of each reference tool in a Markdown file to its website. |
| `tools/ci_annotate.sh <title> <logs>` | Show the important part of log files as GitHub annotations (used by the workflows). |

Every tool prints its usage with `--help`.

## Releasing

See [docs/RELEASING.md](docs/RELEASING.md).
