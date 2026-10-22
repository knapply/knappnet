## Submission

This is a new submission.

<!-- Sections marked "confirm" are filled in from the release checklist in
     docs/RELEASING.md before each submission. -->

## Test environments (confirm)

* GitHub Actions: R release on Ubuntu, macOS, and Windows; R-devel,
  oldrel-1, and 4.1 (the minimum version) on Ubuntu; R-devel on macOS and
  Windows, oldrel-1 on macOS and Windows, and R release on Linux ARM.
* R-devel built with AddressSanitizer and UndefinedBehaviorSanitizer (gcc),
  and R-devel with valgrind instrumentation (wch1/r-debug): no reports.
* win-builder: R-devel and R release.
* R-hub: <!-- list the platforms run -->.

## R CMD check results (confirm)

0 errors | 1 warning | 1 note

* The warning, on Linux only, is about `abort`, `exit`, and `_exit` in the
  compiled code; see below.
* The note is "New submission".

## Compiled code: `abort`, `exit`, `_exit`

On Linux, R CMD check reports that the shared library contains references
to `abort()`, `exit()`, and `_exit()`. They come from the Rust standard
library, which is linked statically into the package, not from this
package's code. Every C entry point runs the engine inside
`std::panic::catch_unwind` and returns any failure to R as an R error, and
no code calls `std::process::exit` or `std::process::abort`. As in all Rust
code, the one remaining path to `abort()` is a failed memory allocation,
where R itself would also fail. Other CRAN packages that compile Rust
report the same symbols for the same reason.

## Rust

The package follows CRAN's notes on using Rust:

* `SystemRequirements` names Cargo and rustc (>= 1.70). The build looks for
  `cargo` on the `PATH` and in `~/.cargo/bin`, prints the versions of
  `rustc` and `cargo`, and stops with installation instructions if either
  is missing or rustc is older than 1.70.
* The Rust code has no dependencies outside the standard library, so
  nothing is vendored and the build never uses the network (`cargo build
  --offline --locked`). All Rust sources are in `src/rust/`, written for
  this package by its author.
* Cargo runs with `-j 2`, and `CARGO_HOME` is set inside the package's
  build directory and removed afterwards, so nothing is written to the
  user's home directory.
* Intermediate build files are removed after linking.

## Threads

The engine uses the Rust standard library's threads (no OpenMP). Functions
default to all available cores, but use 2 when `_R_CHECK_LIMIT_CORES_` is
set and at most `OMP_THREAD_LIMIT` when that is set, so examples and tests
use at most 2 cores during `R CMD check --as-cran`. Tests that compare
thread counts use 1 and 2.

## Size

The source package is small (the Rust sources are about 0.7 MB). The
installed shared library is about 6 MB on Linux when not stripped, most of
it debugging information from the Rust standard library, and about 2 MB
stripped, so an "installed package size" note may appear on platforms that
do not strip shared libraries.
