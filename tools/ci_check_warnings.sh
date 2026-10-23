#!/usr/bin/env bash
# Fails if an R CMD check log contains any WARNING other than the known one
# about compiled code: on Linux the statically linked Rust standard library
# contains references to abort() and exit(), which R CMD check reports even
# though the engine never calls them (see r/cran-comments.md). The
# workflows let check-r-package stop only on errors there and run this
# afterwards, so every other warning still fails the job.
#
#   tools/ci_check_warnings.sh r/check/*.Rcheck/00check.log
set -u
case "${1:-}" in -h|--help|"") sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0;; esac
status=0
for log in "$@"; do
  other=$(grep -E '\.\.\. (\[[^]]*\] )?WARNING$' "$log" | grep -v 'checking compiled code' || true)
  if [ -n "$other" ]; then
    echo "Unexpected warnings in $log:"
    echo "$other"
    status=1
  fi
done
exit "$status"
