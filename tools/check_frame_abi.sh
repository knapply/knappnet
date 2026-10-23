#!/usr/bin/env sh
# Builds the static library, then compiles, links, and runs a C program
# against the data model section of crates/capi/include/netcapi.h, so that
# section is checked against the exported ABI.
#
#   ./tools/check_frame_abi.sh
set -eu
case "${1:-}" in -h|--help) sed -n '2,7p' "$0" | sed 's/^# \{0,1\}//'; exit 0;; esac
cd "$(dirname "$0")/.."
cargo build --release -p capi --lib
out="${TMPDIR:-/tmp}/netcapi-frame-smoke"
cc -std=c99 -Wall -Wextra -Werror -Icrates/capi/include crates/capi/tests/c/frame_smoke.c \
  target/release/libcapi.a -lpthread -ldl -lm -o "$out"
"$out"
