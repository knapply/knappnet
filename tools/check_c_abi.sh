#!/usr/bin/env sh
# Builds the static library, then compiles, links, and runs a C program
# against crates/capi/include/netcapi.h with the system C compiler, so the
# header is checked against the exported ABI. CI runs it on Linux and
# macOS; on Windows the validation job and R CMD check exercise the ABI.
#
#   ./tools/check_c_abi.sh
set -eu
case "${1:-}" in -h|--help) sed -n '2,7p' "$0" | sed 's/^# \{0,1\}//'; exit 0;; esac
cd "$(dirname "$0")/.."
cargo build --release -p capi --lib
out="${TMPDIR:-/tmp}/netcapi-smoke"
cc -std=c99 -Wall -Wextra -Werror -Icrates/capi/include crates/capi/tests/c/smoke.c \
  target/release/libcapi.a -lpthread -ldl -lm -o "$out"
"$out"
