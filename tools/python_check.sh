#!/usr/bin/env bash
# Builds the Python package (the NetworkX backend) into a virtual environment
# in the workspace, runs its tests, then runs NetworkX's test suite through
# the backend.
#
#   tools/python_check.sh                          # latest NetworkX, full suite
#   tools/python_check.sh --networkx 'networkx==3.5' --python 3.10
#
# --python picks a CPython version (installed with uv when it is not the
# system's python3); without it, python3 is used, or a uv-managed CPython
# when python3 has no development headers. The environment is
# .venv-python (removed and rebuilt on every run).
set -eu
cd "$(dirname "$0")/.."
NX=networkx
WANT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --networkx) NX="$2"; shift 2 ;;
    --python) WANT="$2"; shift 2 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done

rm -rf .venv-python
has_headers() {
  "$1" -c 'import os, sys, sysconfig; sys.exit(not os.path.exists(os.path.join(sysconfig.get_paths()["include"], "Python.h")))'
}
if [ -z "$WANT" ] && has_headers python3; then
  python3 -m venv .venv-python
else
  # A standalone CPython (with headers) from uv.
  [ -x .venv-uv/bin/uv ] || { python3 -m venv .venv-uv && .venv-uv/bin/python -m pip install -q uv; }
  .venv-uv/bin/uv python install "${WANT:-3.12}"
  .venv-uv/bin/uv venv --seed --python "${WANT:-3.12}" .venv-python
fi
PY=.venv-python/bin/python
[ -x "$PY" ] || PY=.venv-python/Scripts/python.exe
"$PY" --version
"$PY" -m pip install -q "$NX" pytest numpy scipy
"$PY" -m pip install -q ./python
"$PY" -c 'import networkx as nx, knappnet; print("networkx", nx.__version__, "knappnet", knappnet.__version__)'
"$PY" -m pytest -q python/tests
NETWORKX_TEST_BACKEND=knappnet NETWORKX_FALLBACK_TO_NX=True \
  "$PY" -m pytest -q -p no:cacheprovider --pyargs networkx
