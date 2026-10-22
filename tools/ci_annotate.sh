#!/usr/bin/env bash
# Shows the end of each given log file as a CI annotation (an error unless
# LEVEL says otherwise), so it is visible on the run's summary page.
#
#   tools/ci_annotate.sh "R CMD check" check/*.Rcheck/00check.log
#   LEVEL=notice tools/ci_annotate.sh Benchmarks bench.log
set -u
title="$1"
shift
for f in "$@"; do
  [ -f "$f" ] || continue
  body=$(tail -n 150 "$f" | sed -e 's/%/%25/g' -e 's/\r/%0D/g' | awk '{printf "%s%%0A", $0}')
  echo "::${LEVEL:-error} title=${title}: $(basename "$f")::${body}"
done
