#!/usr/bin/env bash
# Runs a command, showing its output and keeping a copy in ci-logs/<name>.log
# (or $CI_LOG_DIR/<name>.log) for tools/ci_annotate.sh to summarize when a
# job fails. The exit status is the command's own.
#
#   tools/ci_log.sh compare python validation/compare.py --strict
set -u
case "${1:-}" in -h|--help|"") sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 0;; esac
name="$1"
shift
dir="${CI_LOG_DIR:-ci-logs}"
mkdir -p "$dir"
set -o pipefail
"$@" 2>&1 | tee "$dir/${name}.log"
