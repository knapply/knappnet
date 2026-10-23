"""Engine runners for temporal networks (validation/manifest.d/temporal.yaml).

Each fixture's edges are read as contacts with the times of
validation/derived/temporal.py. Infinities (no walk) are reported as
UNREACHABLE = -1000 and NaN (undefined) as UNDEFINED = -999, as the oracle
does.
"""

from __future__ import annotations

import ctypes
import math
import sys
from pathlib import Path

from engine import ERR_LEN, EngineError, load, register

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from derived.temporal import contact_times  # noqa: E402

UNREACHABLE = -1000.0
UNDEFINED = -999.0

CRITERIA = {"earliest_arrival": 0, "latest_departure": 1, "fastest": 2, "shortest": 3}
DIRECTIONS = {"out": 1, "in": 2}
MODES = {"out": 1, "in": 2, "all": 3}
SCOPES = {"node": 0, "pair": 1}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    i64p, f64p = ctypes.POINTER(i64), ctypes.POINTER(f64)
    h = ctypes.c_void_p
    hp = ctypes.POINTER(h)
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    opts = [f64, f64, f64, f64, i32, i32, f64, i32]
    lib.net_temporal_new.argtypes = [i64, i64p, i64p, f64p, f64p, i64, i32, hp, cp, sz]
    lib.net_temporal_free.argtypes = [h]
    lib.net_temporal_free.restype = None
    lib.net_temporal_paths.argtypes = [h, i32, i32, i64p, i64] + opts + [f64p, cp, sz]
    lib.net_temporal_reach.argtypes = [h, i32] + opts + [f64p, cp, sz]
    lib.net_temporal_closeness.argtypes = [h] + opts + [f64p, f64p, cp, sz]
    lib.net_temporal_windows.argtypes = [h, f64p, f64p, i64, hp, cp, sz]
    lib.net_temporal_active_at.argtypes = [h, f64, hp, cp, sz]
    lib.net_temporal_aggregate.argtypes = [h, i32, f64, f64, hp, cp, sz]
    lib.net_temporal_inter_event.argtypes = [h, i32, i32, hp, cp, sz]
    lib.net_temporal_burstiness.argtypes = [h, i32, i32, hp, cp, sz]
    lib.net_temporal_result_length.argtypes = [h, i32]
    lib.net_temporal_result_length.restype = i64
    lib.net_temporal_result_copy.argtypes = [h, i32, f64p, cp, sz]
    lib.net_temporal_result_free.argtypes = [h]
    lib.net_temporal_result_free.restype = None
    for fn in ("net_temporal_new", "net_temporal_paths", "net_temporal_reach",
               "net_temporal_closeness", "net_temporal_windows", "net_temporal_active_at",
               "net_temporal_aggregate", "net_temporal_inter_event",
               "net_temporal_burstiness", "net_temporal_result_copy"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _mark(v: float) -> float:
    if math.isnan(v):
        return UNDEFINED
    if math.isinf(v):
        return UNREACHABLE
    return v


class Temporal:
    """A temporal network handle for one fixture."""

    def __init__(self, fx, instantaneous: bool = False):
        self.lib = _declare()
        self.n = fx.n
        self.start, self.end = contact_times(fx, instantaneous)
        m = len(fx.src)
        s = (ctypes.c_int64 * max(m, 1))(*fx.src)
        t = (ctypes.c_int64 * max(m, 1))(*fx.tgt)
        a = (ctypes.c_double * max(m, 1))(*self.start)
        b = (ctypes.c_double * max(m, 1))(*self.end)
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_temporal_new(fx.n, s, t, a, b, m, int(fx.directed),
                                         ctypes.byref(self.handle), err, ERR_LEN), err)

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_temporal_free(self.handle)
            self.handle = None

    def result(self, fn, *args) -> list[list[float]]:
        """Calls a function returning a TemporalResult; returns its columns."""
        res = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(fn(self.handle, *args, ctypes.byref(res), err, ERR_LEN), err)
        try:
            cols = []
            col = 0
            while (k := self.lib.net_temporal_result_length(res, col)) >= 0:
                out = (ctypes.c_double * max(k, 1))()
                _check(self.lib.net_temporal_result_copy(res, col, out, err, ERR_LEN), err)
                cols.append(list(out)[:k])
                col += 1
            return cols
        finally:
            self.lib.net_temporal_result_free(res)


def _options(spec) -> list:
    traversal = spec.get("traversal", "zero")
    if traversal == "zero":
        code, value = 0, 0.0
    elif traversal == "duration":
        code, value = 1, 0.0
    else:
        code, value = 0, float(traversal)
    return [float(spec.get("start", math.nan)), float(spec.get("end", math.nan)),
            float(spec.get("min_wait", 0.0)), float(spec.get("max_wait", math.inf)),
            int(spec.get("strict", False)), code, value, 0]


def _paths(t: Temporal, spec) -> list[float]:
    """All-pairs matrix, row = origin."""
    n = t.n
    crit = CRITERIA[spec["criterion"]]
    default = "in" if spec["criterion"] == "latest_departure" else "out"
    direction = spec.get("direction", default)
    roots = (ctypes.c_int64 * max(n, 1))(*range(n))
    out = (ctypes.c_double * max(n * n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(t.lib.net_temporal_paths(t.handle, crit, DIRECTIONS[direction], roots, n,
                                    *_options(spec), out, err, ERR_LEN), err)
    m = list(out)[:n * n]
    if direction == "in":
        m = [m[j * n + i] for i in range(n) for j in range(n)]
    return [_mark(v) for v in m]


def _reach(t: Temporal, spec, direction: str) -> list[float]:
    out = (ctypes.c_double * max(t.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(t.lib.net_temporal_reach(t.handle, DIRECTIONS[direction], *_options(spec), out,
                                    err, ERR_LEN), err)
    return list(out)[:t.n]


def _closeness(t: Temporal, spec, distances: bool) -> list[float]:
    n = t.n
    c = (ctypes.c_double * max(n, 1))()
    d = (ctypes.c_double * max(n * n, 1))() if distances else None
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(t.lib.net_temporal_closeness(t.handle, *_options(spec), c, d, err, ERR_LEN), err)
    values = list(d)[:n * n] if distances else list(c)[:n]
    return [_mark(v) for v in values]


def series(start: list[float], end: list[float], width: float, step: float):
    """Window starts and ends: a = first start, + step, ... while a <= last end."""
    if not start:
        return [], []
    lo, hi = min(start), max(end)
    a, k = [], 0
    while lo + k * step <= hi:
        a.append(lo + k * step)
        k += 1
    return a, [x + width for x in a]


def grid(start: list[float], end: list[float]) -> list[float]:
    """Times from one unit before the first start to one after the last end, by halves."""
    if not start:
        return [0.0]
    lo, hi = min(start) - 1, max(end) + 1
    return [lo + k / 2 for k in range(int(2 * (hi - lo)) + 1)]


@register("temporal_paths")
def _run_paths(g, spec, fx):
    return _paths(Temporal(fx, spec.get("instantaneous", False)), spec)


@register("temporal_reach")
def _run_reach(g, spec, fx):
    return _reach(Temporal(fx), spec, spec["direction"])


@register("temporal_reachability_ratio")
def _run_ratio(g, spec, fx):
    if fx.n < 2:
        return [UNDEFINED]
    return [sum(_reach(Temporal(fx), spec, "out")) / (fx.n * (fx.n - 1))]


@register("temporal_closeness")
def _run_closeness(g, spec, fx):
    return _closeness(Temporal(fx), spec, False)


@register("temporal_average_distance")
def _run_average_distance(g, spec, fx):
    return _closeness(Temporal(fx), spec, True)


@register("temporal_active")
def _run_active(g, spec, fx):
    t = Temporal(fx)
    return [tuple(t.result(t.lib.net_temporal_active_at, x)[0]) for x in grid(t.start, t.end)]


@register("temporal_windows")
def _run_windows(g, spec, fx):
    t = Temporal(fx)
    a, b = series(t.start, t.end, float(spec["width"]), float(spec["step"]))
    k = len(a)
    fa = (ctypes.c_double * max(k, 1))(*a)
    fb = (ctypes.c_double * max(k, 1))(*b)
    offsets, idx = t.result(t.lib.net_temporal_windows, fa, fb, k)
    return [tuple(idx[int(offsets[w]):int(offsets[w + 1])]) for w in range(k)]


@register("temporal_aggregate")
def _run_aggregate(g, spec, fx):
    t = Temporal(fx)
    w = spec.get("window")
    cols = t.result(t.lib.net_temporal_aggregate, int(w is not None),
                    float(w[0]) if w else 0.0, float(w[1]) if w else 0.0)
    return [tuple(c[i] for c in cols[:6]) for i in range(len(cols[0]))]


@register("temporal_burstiness")
def _run_burstiness(g, spec, fx):
    t = Temporal(fx)
    scope = spec["scope"]
    events, b, m = t.result(t.lib.net_temporal_burstiness, SCOPES[scope],
                            MODES[spec.get("mode", "all")])
    rows = [(e, _mark(x), _mark(y)) for e, x, y in zip(events, b, m)]
    if scope == "pair":
        cols = t.result(t.lib.net_temporal_aggregate, 0, 0.0, 0.0)
        rows = [(cols[0][i], cols[1][i]) + r for i, r in enumerate(rows)]
    return rows


@register("temporal_inter_event")
def _run_inter_event(g, spec, fx):
    t = Temporal(fx)
    group, gap, _ = t.result(t.lib.net_temporal_inter_event, SCOPES[spec["scope"]],
                             MODES[spec.get("mode", "all")])
    return list(zip(group, gap))
