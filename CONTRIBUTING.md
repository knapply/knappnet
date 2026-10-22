# Contributing

Thanks for your interest. Accuracy is the first priority of this project, so
changes follow a few firm rules. Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)
first.

## Adding a measure

1. Implement it in `crates/engine/src/metrics/`. Start the file with a comment
   that cites the source of the definition and the algorithm (add the entries
   to `docs/references.bib`) and lists every convention: loops, parallel
   edges, weights, normalization, disconnected graphs.
2. Add unit tests with hand-checked values, and property tests where an
   invariant applies.
3. Expose it in `crates/capi` and `crates/capi/include/netcapi.h`.
4. Add a manifest entry in `validation/manifest.yaml` and oracle code for each
   reference tool that implements it.
5. Expose it in R (`r/R/`, `r/src/glue.c`, `r/src/init.c`, `r/man/`) with the
   same conventions and references in the help page.
6. Run the checks below.

## Checks

```sh
cargo test --workspace                     # engine unit and property tests
./tools/check_c_abi.sh                     # C header matches the ABI
python3 tools/sync_r_sources.py            # copy engine sources into the R package
python3 validation/fixtures.py             # regenerate fixtures (deterministic)
python3 validation/oracles/definition_oracle.py
python3 validation/oracles/networkx_oracle.py
Rscript validation/oracles/r_oracles.R     # needs igraph, sna, yaml
python3 validation/compare.py              # must report 0 failures
python3 tools/rename.py --self-test        # renaming still works
R CMD build r && R CMD check --as-cran knappnet_*.tar.gz
```

## Style

- Rust: `cargo fmt` and `cargo clippy`.
- R: base R only in package code; no new hard dependencies without discussion.
- Comments explain *why* and cite sources; they do not narrate the code.

## Renaming

The project name lives in `project.toml`. To change it, or after renaming the
GitHub repository, run `python3 tools/rename.py` (see `--help`). To change
the name or email on past commits, add an entry to `.mailmap`.
