# Installing from source

knappnet is built from the repository; it is not on PyPI.

## Requirements

- CPython 3.10 or later (one build serves every later version).
- NetworkX 3.5 or later. The backend is tested with NetworkX 3.7.
- A Rust toolchain, 1.70 or later (<https://rustup.rs>). The engine has no
  third-party crates, so the build needs no network access for Rust.
- A C compiler: GCC or Clang on Linux, the Xcode command line tools on
  macOS, Microsoft C++ Build Tools on Windows.

## Install

From a clone of the repository:

```sh
python -m pip install ./python
```

The build compiles the engine's C interface (`crates/capi`) with Cargo and
links it into the extension module `knappnet._core`. Cargo's target
directory (`target/`, or `CARGO_TARGET_DIR`) is shared with the rest of the
repository, so a second build is quick.

To run the backend's tests as well:

```sh
python -m pip install "./python[test]"
python -m pytest python/tests
```

## Check the installation

```python
import networkx as nx
print("knappnet" in nx.utils.backends.backends)   # True
print(nx.config.backends.knappnet)                # Config(threads=0, min_nodes=64)
```
