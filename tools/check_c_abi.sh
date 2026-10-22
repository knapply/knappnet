#!/usr/bin/env sh
# Builds the static library and links a C program against netcapi.h.
set -eu
cd "$(dirname "$0")/.."
cargo build --release -p capi --lib
out="${TMPDIR:-/tmp}/netcapi-smoke"
cc -std=c99 -Wall -Wextra -Werror -Icrates/capi/include crates/capi/tests/c/smoke.c \
  target/release/libcapi.a -lpthread -ldl -lm -o "$out"
"$out"
