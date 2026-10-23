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
  case "$f" in
    *.Rout.fail)
      # Test logs can contain very long diffs; show each failure's header
      # and first lines, then the summary at the end.
      text=$( (grep -n -A8 -E '^── (Failure|Error|Warning)' "$f" | head -n 200; echo '...'; tail -n 15 "$f") )
      ;;
    *)
      text=$(tail -n 150 "$f")
      ;;
  esac
  body=$(printf '%s\n' "$text" | sed -e 's/%/%25/g' -e 's/\r/%0D/g' | awk '{printf "%s%%0A", $0}')
  echo "::${LEVEL:-error} title=${title}: $(basename "$f")::${body}"
done
