# Validation

Every measure the engine exposes is compared, on every fixture, with every
reference implementation that computes it and with a direct evaluation of
its definition. This directory holds the corpus, the oracles that produce
golden values, and the comparison. CI runs all of it on Linux, macOS, and
Windows on every push; `REPORT.md` is the latest result.

## Layout

```
fixtures.py              generates the corpus (deterministic; CI checks it)
fixtures/<name>/         edges.csv (0-based endpoints, weight), meta.json
fixtures-extra/          datasets for the R package's examples, not oracles
manifest.d/<area>.yaml   what is compared and how (see "Manifest" below)
bindings/<area>.py       ctypes runners: one engine call per manifest `fn`
oracles/<area>_definition_oracle.py   definition oracles (Python, exact)
oracles/<area>_networkx_oracle.py     NetworkX oracles
oracles/r_<area>.R                    igraph, sna, and other R oracles
derived/                 data derived from the fixtures that several
                         oracles share (contact times, SAOM panels)
golden/<oracle>/<fixture>/<variant>.csv   oracle answers; VERSION per oracle
engine.py                loads the C ABI and dispatches `fn` to a runner
common.py                fixtures, manifests, golden files: shared helpers
run_oracles.py           runs every oracle script
compare.py               engine vs golden values; writes REPORT.md
check_manifest.py        shape check of every manifest entry
io_check.py, qap_check.py, cug_check.py, layout_check.py
                         checks that do not fit the manifest (below)
scale/                   the same comparison on larger generated networks
UCINET.md                how the UCINET golden files are produced by hand
```

The validation goes through the same C ABI the R package uses, so a
passing comparison also exercises the exported interface.

## Running it

```sh
cargo build --release --offline -p capi        # the library compare.py loads
python3 validation/check_manifest.py           # manifest entries well-formed
python3 validation/run_oracles.py --python-only   # definition and NetworkX
python3 validation/compare.py                  # must print FAIL 0
```

`run_oracles.py` without `--python-only` also runs the R scripts, which
need R with [igraph](https://r.igraph.org/), [sna](https://cran.r-project.org/package=sna), and yaml (and, for some areas, network, [ergm](https://cran.r-project.org/package=ergm),
btergm, [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/), tsna, [keyplayer](https://cran.r-project.org/package=keyplayer), and homophily; a script whose package is
missing says so and is skipped). `ENGINE_LIB=path/to/libcapi.so` points
`compare.py` at a library built elsewhere. `compare.py --strict` (CI) also
fails when an oracle never ran.

The golden files of the Python oracles are committed with the oracle
version that produced them (`golden/<oracle>/VERSION`) and regenerated when
that version changes. The R oracles run in CI and their golden files are
uploaded as a build artifact, so locally their comparisons show as "not
run".

## Outcomes

`compare.py` classifies each (measure, oracle, fixture) as:

| Outcome | Meaning |
|---|---|
| match | values agree within the manifest tolerance |
| documented | values differ on a fixture the manifest declares as a known convention difference for that oracle, with the reason |
| FAIL | any other disagreement (CI fails) |
| not run | no golden file: the oracle was not available |
| undefined | the manifest declares that the oracle has no answer on that fixture, with the reason |
| manual | no golden file yet for an oracle run by hand ([UCINET](https://sites.google.com/site/ucinetsoftware/home)) |

A measure is *validated* when it has a manifest entry, at least one
passing oracle, and no FAIL, and every documented difference and undefined
case carries a reason that a reader can check against the oracle's
documentation. A definition oracle is strongly preferred: it is the
strongest evidence, and when a reference tool disagrees with both the
engine and the definition oracle, the difference is the tool's.

## Manifest

One file per area in `manifest.d/`; `core.yaml` is the reference example.
Keys are unique across all files and every entry has this shape:

```yaml
betweenness_normalized:
  engine: {fn: betweenness, weighted: false, normalized: true}
  applies_to: any            # any | directed | undirected
  requires_weights: false    # optional; true skips unweighted fixtures
  compare: values            # optional; values | sorted | partition | sets
  tolerance: {rel: 1.0e-9, abs: 1.0e-12}
  oracles:
    definition:
      call: "direct evaluation of the definition (core_definition_oracle.py)"
    networkx:
      call: "nx.betweenness_centrality(G, normalized = True)"
      differs_when: [multi_edges]
      reason: NetworkX counts parallel edges as a single shortest path.
    igraph:
      call: "igraph::betweenness(g, weights = NA, normalized = TRUE)"
      differs_when: [edge_case]
      reason: >-
        With fewer than three nodes the normalizing constant is zero.
        igraph returns NaN; this engine, like NetworkX, returns 0.
```

- `engine.fn` names a runner registered with `@register("fn")` in
  `bindings/<area>.py`; the other `engine` keys are its options and are
  documented in the comment at the top of the area's manifest.
- `compare`: `values` (node-aligned vectors, scalars as one-element
  lists), `sorted` (order does not matter), `partition` (group labels are
  arbitrary), `sets` (a collection of groups, such as cliques).
- `differs_when` (fixture tags) and `differs_on` (fixture names) declare
  where an oracle follows a different convention; `reason` says why, in
  terms a user could verify. On those fixtures a mismatch is "documented";
  anywhere else it fails.
- `undefined_when` and `undefined_on` declare fixtures where the oracle
  gives no answer (it errors, or the value is undefined), with
  `undefined_reason`. Without them a missing golden file is "not run" and
  fails `--strict`.
- Fixture tags: `directed`, `undirected`, `weighted`, `integer_weights`,
  `loops`, `multi_edges`, `disconnected`, `edge_case`, `twomode`
  (`fixtures.py` computes them; `check_manifest.py` rejects others).
- Undefined values (NaN) in results: each area states in its manifest
  comment how its runner and oracles encode them (an empty line, `-1`
  for measures that are never negative, or `-999`/`-1000` markers), so
  that undefined only matches undefined.

## Oracles

A Python oracle is a script `oracles/<area>_<oracle>_oracle.py` that reads
its area's entries with `common.load_area("<area>")`, computes each
variant it implements for every fixture it applies to, writes the result
with `common.write_values(oracle, fixture, variant, values)` (a value is a
number or a tuple of numbers, one group per line), and records its version
with `common.write_version(oracle, ...)`. `run_oracles.py` runs every file
matching `*_oracle.py`, then every `*.R`.

A **definition oracle** is a slow, obviously correct implementation of the
textbook definition that shares no code or algorithmic idea with the
engine: exact rational arithmetic (`fractions.Fraction`), brute-force
enumeration, dense linear algebra. It runs on fixtures small enough for
that (`MAX_NODES`). A **reference-tool oracle** calls the tool's own
function with the arguments written in the manifest's `call`, so a reader
can reproduce any golden value by hand.

R oracle scripts (`oracles/r_<area>.R`) read the same manifest and write
the same files. They run from the repository root, parse `meta.json`
without a JSON package, and each carries its own small helpers because R
is not available where the Python side is developed; keep them
self-contained. Licensing: reference tools under the GPL ([igraph](https://r.igraph.org/), [sna](https://cran.r-project.org/package=sna),
[statnet](https://statnet.org/), [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/), ...) are run as black boxes; their source is never read
or translated.

## Checks outside the manifest

Some things are not node-level values on a fixture and have their own
scripts, each a definition check independent of the engine:

- `io_check.py`: every file format, written by the engine and read by
  [NetworkX](https://networkx.org/) and the reverse, on every fixture and an attribute-rich network.
- `qap_check.py`: QAP and MRQAP against numpy least squares and an exact
  enumeration of all permutations for small n.
- `cug_check.py`: random graph generators and the null models of
  `cug_test()` against their defining constraints and distributions.
- `layout_check.py`: layouts by their properties (finite, reproducible
  across seeds and threads, components separated, stress no worse than
  NetworkX's Kamada-Kawai).
- `scale/`: the manifest comparison repeated on generated networks of up
  to 50,000 nodes against NetworkX, [igraph](https://r.igraph.org/), and [sna](https://cran.r-project.org/package=sna) (weekly, `bench.yml`).

## Adding a measure

`python3 tools/new_measure.py <area> <name>` writes the manifest entry,
runner, and oracle stubs; CONTRIBUTING.md has the full recipe. In short:
a manifest entry in `manifest.d/<area>.yaml`, a runner in
`bindings/<area>.py`, a definition oracle, every reference tool that
implements the measure as a further oracle (with each known difference
documented), `run_oracles.py --python-only`, `compare.py` with `FAIL 0`,
and the regenerated golden files committed.
