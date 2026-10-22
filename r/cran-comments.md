## R CMD check notes

### Compiled code: `abort`, `exit`, `_exit`

On Linux, R CMD check reports that the shared library contains references
to `abort()`, `exit()`, and `_exit()`. They come from the Rust standard
library, which is linked statically into `src/rust/target/release/libcapi.a`,
not from this package's code. The package never terminates R: every C entry
point runs the engine inside `std::panic::catch_unwind` and returns any
failure to R as an error condition, and no code calls `std::process::exit`
or `std::process::abort`. Other CRAN packages that compile Rust report the
same symbols for the same reason.
