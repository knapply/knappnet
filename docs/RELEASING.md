# Releasing

The checklist for a release of the R package to CRAN and of the engine
crate to crates.io. Work through it in order; each step names the command
or the place to look. A release is cut from `main` with CI green.

## Blockers for the first CRAN release

- [ ] **Maintainer email.** `Authors@R` in `r/DESCRIPTION` gives
      `knapply@users.noreply.github.com`, which does not receive mail. CRAN
      writes to the maintainer (and sends a confirmation link on every
      submission), so it must be an address the maintainer reads. Change it
      in `r/DESCRIPTION` and `CITATION.cff` before submitting.
- [ ] **First version number.** `r/DESCRIPTION`, `Cargo.toml`
      (`workspace.package.version`), and `CITATION.cff` say `0.0.1`. Decide
      the first released version (for example `0.1.0`) and set it in all
      three, then run `python3 tools/sync_r_sources.py`.

## 1. Versions and news

- [ ] Choose the version. The R package and the crate share
      `workspace.package.version`; bump `r/DESCRIPTION` to the same number.
- [ ] Update `CITATION.cff` (`version`, and `date-released` if present).
- [ ] In `r/NEWS.md`, rename the top heading from
      `(development version)` to the version, and check that every
      user-visible change since the last release is listed: new functions,
      changed defaults, new errors or warnings, and bug fixes that change
      results. After the release, add a new `(development version)` heading.
- [ ] `python3 tools/sync_r_sources.py`, then commit. CI fails if the copy
      of the engine in `r/src/rust/` is stale.

## 2. Local checks

Run from the repository root (see `CONTRIBUTING.md` for what each does):

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc -p knappnet --no-deps
./tools/check_c_abi.sh
python3 tools/sync_r_sources.py --check
python3 validation/run_oracles.py          # needs R with igraph, sna, yaml
python3 validation/compare.py --strict     # 0 failures
Rscript tools/check_readme.R               # README output matches R
R CMD build r
R CMD check --as-cran knappnet_*.tar.gz
```

- [ ] `R CMD check --as-cran` gives no errors, no warnings other than the
      known Linux warning about `abort`/`exit` (explained in
      `r/cran-comments.md`), and no notes other than "New submission" (and
      possibly the installed size note, also explained there).
- [ ] The check took less than 10 minutes in total, and no example takes
      more than 5 seconds (`R CMD check` reports slow examples in
      `knappnet.Rcheck/knappnet-Ex.timings`).
- [ ] The tarball is under 5 MB (`ls -l knappnet_*.tar.gz`) and contains no
      build products (`tar tzf knappnet_*.tar.gz | grep -E 'target|\.o$|\.so$'`
      prints nothing).
- [ ] Nothing was written outside the temporary directory: after the check,
      `git status` is clean and no files appeared in the home directory.

## 3. Checks on every platform

- [ ] CI on `main` is green: engine tests, validation on Linux, macOS, and
      Windows, R CMD check on R release (three platforms) and R-devel,
      oldrel-1, and 4.1 (Linux), the AddressSanitizer/UBSan and valgrind
      jobs, and the README check.
- [ ] Run the "R CMD check (extra platforms)" workflow (`r-check.yml`,
      from the Actions tab, or by pushing the release tag; see step 6):
      R-devel and oldrel-1 on macOS and Windows, and Linux on ARM.
- [ ] win-builder, R-devel and R release. From R, in `r/`:
      `devtools::check_win_devel()` and `devtools::check_win_release()`, or
      upload the tarball at <https://win-builder.r-project.org/upload.aspx>.
      Results arrive by email at the maintainer address, which is one more
      reason it must be real.
- [ ] macOS builder: <https://mac.r-project.org/macbuilder/submit.html>.
- [ ] R-hub v2: `rhub::rhub_check()` from `r/` (runs on GitHub Actions in
      a fork set up with `rhub::rhub_setup()`). At least the `linux`,
      `macos-arm64`, `windows`, `clang-asan`, `valgrind`, `nosuggests`, and
      `rchk` platforms. `rchk` checks PROTECT balance in the C glue; fix
      anything it reports.
- [ ] Check without Rust installed: on a machine or container without
      `cargo`, `R CMD INSTALL knappnet_*.tar.gz` must stop with the message
      from `r/src/Makevars` telling the user how to install Rust.
- [ ] Reverse dependencies: none yet. From the second release on, run
      `revdepcheck::revdep_check()` and report the results in
      `cran-comments.md`.

## 4. cran-comments.md

- [ ] Update `r/cran-comments.md`: the test environments actually used and
      the check results from the final runs (every section marked
      "confirm"). Keep the explanations of the compiled-code warning, the
      Rust toolchain, threads, and size up to date.
- [ ] For a resubmission, add a "Resubmission" section at the top that
      answers each point the CRAN reviewers raised.

## 5. Submit to CRAN

- [ ] `R CMD build r` on a clean checkout, then submit the tarball at
      <https://cran.r-project.org/submit.html> (or `devtools::submit_cran()`
      from `r/`). Confirm the submission from the email CRAN sends.
- [ ] Reply to reviewers from the maintainer address; keep replies short
      and fix what they ask before resubmitting with an incremented
      version only if they request it.

## 6. After acceptance

- [ ] Tag the release on the commit that was submitted:
      `git tag -a vX.Y.Z -m "knappnet X.Y.Z"` and `git push origin vX.Y.Z`.
      The tag runs `r-check.yml`.
- [ ] Publish the engine crate: the `publish-crate.yml` workflow, or
      `cargo publish -p knappnet` (run `cargo package -p knappnet` first).
- [ ] Create a GitHub release from the tag with the `r/NEWS.md` section.
- [ ] Add a new `# knappnet (development version)` heading to `r/NEWS.md`
      and bump the version to `X.Y.Z.9000` in `r/DESCRIPTION`.

## Known limitations to revisit

- A long computation cannot be interrupted from R: the engine runs on its
  own threads, which may not call `R_CheckUserInterrupt()`. Supporting it
  needs a cancellation flag in the C ABI that the engine's workers poll
  between work items, set by the C glue from a check on the R thread
  (for example by running the engine call on a helper thread and polling
  `R_CheckUserInterrupt()` through `R_ToplevelExec()` while waiting).
- A panic inside the engine is caught and returned as an R error, but
  Rust's default panic hook also prints the panic message to standard
  error first. CRAN asks compiled code not to write to the console; a
  silent panic hook installed by the C ABI would avoid it.
- Writing CSV and TSV edge lists goes through `utils::write.table()`, which
  cannot write characters outside the native encoding in R before 4.2 on
  Windows and in other non-UTF-8 locales. Every other format is UTF-8 on
  every platform.
