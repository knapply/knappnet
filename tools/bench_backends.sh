#!/usr/bin/env bash
# Times NetworkX and its backends with bench/networkx_backend.py in a fresh
# virtual environment.
#
#   tools/bench_backends.sh cpu    # knappnet, nx-parallel, graphblas-algorithms
#   tools/bench_backends.sh gpu    # knappnet and nx-cugraph (needs an NVIDIA GPU)
#   tools/bench_backends.sh cpu --quick --only pagerank
#
# On the btower runner: dispatch with runner=bench for cpu, runner=gpu for gpu.
set -eu
cd "$(dirname "$0")/.."
KIND="${1:-cpu}"
[ $# -gt 0 ] && shift
tools/python_check.sh --build-only
PY=.venv-python/bin/python
case "$KIND" in
  cpu)
    "$PY" -m pip install -q nx-parallel graphblas-algorithms
    BACKENDS="knappnet parallel graphblas" ;;
  gpu)
    "$PY" -m pip install -q nx-cugraph-cu12 --extra-index-url https://pypi.nvidia.com
    BACKENDS="knappnet cugraph" ;;
  *) echo "usage: $0 cpu|gpu [bench options]" >&2; exit 2 ;;
esac
"$PY" -m pip list 2>/dev/null | grep -i -E "^(networkx|nx-|graphblas|python-graphblas|cugraph|knappnet)" || true
NX_CUGRAPH_AUTOCONFIG=False "$PY" bench/networkx_backend.py --backends $BACKENDS "$@"
