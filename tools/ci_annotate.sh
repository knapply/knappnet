#!/usr/bin/env bash
# Shows the important part of each given log file as CI annotations (errors
# unless LEVEL says otherwise), so it is visible on the run's summary page.
# GitHub truncates long annotation messages, so the text is split into
# numbered parts.
#
#   tools/ci_annotate.sh "R CMD check" check/*.Rcheck/00check.log
#   LEVEL=notice tools/ci_annotate.sh Benchmarks bench.log
set -u
case "${1:-}" in -h|--help|"") sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'; exit 0;; esac
title="$1"
shift
for f in "$@"; do
  [ -f "$f" ] || continue
  case "$f" in
    *.Rout.fail)
      # Failures and errors with their first lines, warnings as one line
      # each, then the summary at the end.
      # Every failure's heading first (GitHub shows at most ten
      # annotations per step, so the details that follow are cut short).
      text=$( (echo '--- failures ---'
               grep -E '^── (Failure|Error)' "$f" | sed 's/ ─*$//' | head -n 170
               echo '--- details ---'
               grep -n -A6 -E '^── (Failure|Error)' "$f" | head -n 140
               echo '--- warnings ---'
               grep -n -A1 -E '^── Warning' "$f" | grep -v '^--$' | head -n 60
               echo '--- end of log ---'
               tail -n 12 "$f") )
      ;;
    *-Ex.Rout)
      text=$(tail -n 60 "$f")
      ;;
    *)
      text=$(tail -n 120 "$f")
      ;;
  esac
  part=1
  # Split into chunks of at most 60 lines.
  printf '%s\n' "$text" | split -l 60 - "${TMPDIR:-/tmp}/annot_" 
  for chunk in "${TMPDIR:-/tmp}"/annot_*; do
    body=$(sed -e 's/%/%25/g' -e 's/\r/%0D/g' "$chunk" | awk '{printf "%s%%0A", $0}')
    echo "::${LEVEL:-error} title=${title}: $(basename "$f") (part ${part})::${body}"
    part=$((part + 1))
    rm -f "$chunk"
  done
done
