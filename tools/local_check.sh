#!/usr/bin/env bash
# Runs every check CI runs, on one machine, and prints a compact summary:
# one line per step with its status and time, then only the lines that
# explain each failure. Full logs go to check-logs/.
#
#   tools/local_check.sh                 # everything
#   tools/local_check.sh engine valid    # only steps whose name contains a word
#
# Steps: lint, engine, abi, docs, oracles, compare, checks, python,
# rcheck. Needs Rust, Python with validation/requirements.txt, and R with
# the oracle packages (the btower runner image has all of them).
set -u
cd "$(dirname "$0")/.."
LOG=check-logs
mkdir -p "$LOG"
want=("$@")
declare -a SUMMARY
fails=0

selected() {
  [ ${#want[@]} -eq 0 ] && return 0
  for w in "${want[@]}"; do [[ "$1" == *"$w"* ]] && return 0; done
  return 1
}

step() {
  local name="$1"; shift
  selected "$name" || return 0
  local t0=$SECONDS
  if bash -c "$*" >"$LOG/$name.log" 2>&1; then
    SUMMARY+=("ok    $(printf '%-11s' "$name") $((SECONDS - t0))s")
  else
    SUMMARY+=("FAIL  $(printf '%-11s' "$name") $((SECONDS - t0))s  -> $LOG/$name.log")
    fails=$((fails + 1))
  fi
}

step lint 'cargo fmt --all --check && cargo clippy -q --workspace --all-targets -- -D warnings \
  && python3 tools/sync_r_sources.py --check && python3 tools/check_bib.py \
  && python3 validation/check_manifest.py && python3 tools/check_r_docs.py'
step engine 'cargo test -q --workspace --release'
step abi 'bash tools/check_c_abi.sh && bash tools/check_frame_abi.sh'
step docs 'RUSTDOCFLAGS="-D warnings" cargo doc -q --no-deps --workspace'
step oracles 'python3 validation/run_oracles.py'
step compare 'python3 validation/compare.py --strict'
step checks 'python3 -W ignore validation/io_check.py && python3 validation/qap_check.py \
  && python3 validation/cug_check.py && python3 validation/layout_check.py'
# The NetworkX backend: built into .venv-python, its tests, then NetworkX's
# whole test suite through it (tools/python_check.sh).
step python 'bash tools/python_check.sh'
step rcheck 'cd r && rm -f knappnet_*.tar.gz && R CMD build --no-build-vignettes . \
  && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --ignore-vignettes knappnet_*.tar.gz \
  && ! grep -q "Status:.*WARNING" knappnet.Rcheck/00check.log'

echo "== summary"
printf '%s\n' "${SUMMARY[@]}"
for f in "$LOG"/*.log; do
  name=$(basename "$f" .log)
  printf '%s\n' "${SUMMARY[@]}" | grep -q "^FAIL  $name " || continue
  echo "== $name"
  case "$name" in
    compare) grep -E '^(FAIL|match )' "$f" | head -60 ;;
    engine) grep -E 'panicked|FAILED|^failures:|error(\[|:)' "$f" | head -40 ;;
    python) grep -E '^(FAILED|ERROR)|error:|Error:|[0-9]+ (passed|failed)' "$f" | head -40 ;;
    rcheck) sed -n '/Status:/p;/^\* checking .*\.\.\. \(ERROR\|WARNING\)/,/^\* /p' r/knappnet.Rcheck/00check.log 2>/dev/null | head -60
            grep -E 'Failure|Error' r/knappnet.Rcheck/tests/testthat.Rout.fail 2>/dev/null | head -30 ;;
    *) grep -E -i 'error|warning|fail|diff' "$f" | head -40 ;;
  esac
done
exit $fails
