"""Engine runners for key players, group centrality, and flow betweenness.

Undefined results (NaN) are written as -1, since none of these measures can
be negative; every oracle for them uses the same marker.

Manifest engine options:
  fn: kpp_single      the key player criterion of every single node, as a
                      node-aligned vector (the set {v} for each node v)
  fn: kpp_optimum     the best value over all sets of size `k` (exact
                      enumeration), a scalar
  fn: group           a group centrality of the fixture's standard group
                      (group_members() below)
  fn: flow_betweenness
  measure, m, mode, weighted, weights_as, normalized, scale, capacity as
  in the C interface (crates/capi/include/cyber.h).
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

UNDEFINED = -1.0

KPP = {"fragmentation": 0, "distance_fragmentation": 1, "m_reach": 2, "distance_reach": 3}
GROUP = {"degree": 0, "closeness": 1, "betweenness": 2}
MODES = {"out": 1, "in": 2, "all": 3}
WEIGHTS_AS = {"distance": 0, "strength": 1}
SCALE = {"raw": 0, "normalized": 1, "fractional": 2}

_declared = False


def group_members(n: int) -> list[int]:
    """The group used for group centrality on a fixture with n nodes: the
    first three nodes, leaving at least two non-members when possible."""
    return list(range(min(3, max(1, n - 2))))


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    f64p, i64p = ctypes.POINTER(f64), ctypes.POINTER(i64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_kpp_score.argtypes = [h, i64p, i64, i32, f64, i32, i32, i32, f64p, cp, sz]
    lib.net_key_players.argtypes = [h, i64, i32, f64, i32, i32, i32, i32, i64, u64, i32,
                                    i64p, f64p, f64p, cp, sz]
    lib.net_group_centrality.argtypes = [h, i64p, i64, i32, i32, i32, i32, i32, i32, f64p,
                                         cp, sz]
    lib.net_flow_betweenness.argtypes = [h, i32, i32, i32, f64p, cp, sz]
    for fn in ("net_kpp_score", "net_key_players", "net_group_centrality",
               "net_flow_betweenness"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _mark(x: float) -> float:
    return UNDEFINED if math.isnan(x) else x


def _kpp_args(spec):
    return (KPP[spec["measure"]], float(spec.get("m", 1)), MODES[spec.get("mode", "out")],
            int(spec.get("weighted", False)), WEIGHTS_AS[spec.get("weights_as", "distance")])


def kpp_score(g, members, spec) -> float:
    lib = _declare()
    k = len(members)
    arr = (ctypes.c_int64 * max(k, 1))(*members)
    out = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_kpp_score(g.handle, arr, k, *_kpp_args(spec), ctypes.byref(out), err,
                             ERR_LEN), err)
    return out.value


def key_players(g, k, spec, method="exact", restarts=1, seed=1, threads=0):
    lib = _declare()
    members = (ctypes.c_int64 * max(k, 1))()
    score, evals = ctypes.c_double(), ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_key_players(g.handle, k, *_kpp_args(spec), {"greedy": 0, "exact": 1}[method],
                               restarts, seed, threads, members, ctypes.byref(score),
                               ctypes.byref(evals), err, ERR_LEN), err)
    return list(members)[:k], score.value


@register("kpp_single")
def _kpp_single(g, spec, fx):
    return [_mark(kpp_score(g, [v], spec)) for v in range(g.n)]


@register("kpp_optimum")
def _kpp_optimum(g, spec, fx):
    k = int(spec["k"])
    if g.n < k:
        return [UNDEFINED]
    return [_mark(key_players(g, k, spec)[1])]


@register("group")
def _group(g, spec, fx):
    lib = _declare()
    members = group_members(g.n)
    arr = (ctypes.c_int64 * len(members))(*members)
    out = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_group_centrality(
        g.handle, arr, len(members), GROUP[spec["measure"]], MODES[spec.get("mode", "out")],
        int(spec.get("weighted", False)), WEIGHTS_AS[spec.get("weights_as", "distance")],
        int(spec.get("normalized", True)), 0, ctypes.byref(out), err, ERR_LEN), err)
    return [_mark(out.value)]


@register("flow_betweenness")
def _flow_betweenness(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_flow_betweenness(g.handle, int(spec.get("capacity", False)),
                                    SCALE[spec.get("scale", "raw")], 0, out, err, ERR_LEN), err)
    return [_mark(x) for x in list(out)[:g.n]]
