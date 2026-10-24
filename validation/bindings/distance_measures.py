"""Engine runners for the distance measures (crates/capi/src/distance.rs):
center, periphery, barycenter, the Wiener family, harmonic diameter,
efficiency, and closeness vitality.

Markers: NaN is written as UNDEFINED = -999, infinity as INFINITE = -1000,
and minus infinity as NEG_INFINITE = -1001. Node sets (center, periphery,
barycenter) are returned as their sorted node indices (empty when the
measure needs a connected network and the network is not).

Manifest engine options:
  fn: dm_center, dm_periphery   reachable: require | ignore
  fn: dm_barycenter             tie: relative | exact
  fn: dm_index                  index: wiener | schultz | gutman |
                                hyper_wiener | hyper_wiener_networkx |
                                harmonic_diameter | global_efficiency
  fn: dm_local_efficiency       loop_is_neighbor
  fn: dm_closeness_vitality
  mode: out | in | all; weighted; weights_as: distance | strength.
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

UNDEFINED = -999.0
INFINITE = -1000.0
NEG_INFINITE = -1001.0
MODES = {"out": 1, "in": 2, "all": 3}
WEIGHTS_AS = {"distance": 0, "strength": 1}
INDEX = {"wiener": 0, "schultz": 1, "gutman": 2, "hyper_wiener": 3,
         "hyper_wiener_networkx": 4, "harmonic_diameter": 5, "global_efficiency": 6}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p = ctypes.POINTER(f64), ctypes.POINTER(i64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    sig = {
        "net_center_periphery": [h, i32, i32, i32, i32, i32, i32, i64p, i64p, cp, sz],
        "net_barycenter": [h, i32, i32, i32, i32, i32, i64p, i64p, f64p, cp, sz],
        "net_distance_index": [h, i32, i32, i32, i32, i32, f64p, cp, sz],
        "net_local_efficiency": [h, i32, i32, f64p, cp, sz],
        "net_closeness_vitality": [h, i32, i32, i32, i32, f64p, cp, sz],
    }
    for name, args in sig.items():
        fn = getattr(lib, name)
        fn.argtypes = args
        fn.restype = i32
    _declared = True
    return lib


def _mark(x: float) -> float:
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE if x > 0 else NEG_INFINITE
    return x


def _path(spec):
    return (MODES[spec.get("mode", "out")], int(spec.get("weighted", False)),
            WEIGHTS_AS[spec.get("weights_as", "distance")])


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def _node_set(g, fn, *args):
    out = (ctypes.c_int64 * max(g.n, 1))()
    count = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(g.handle, *args, out, ctypes.byref(count), err, ERR_LEN), err)
    return sorted(float(x) for x in list(out)[:count.value])


@register("dm_center")
def _center(g, spec, fx):
    lib = _declare()
    reach = {"require": 0, "ignore": 1}[spec.get("reachable", "require")]
    return _node_set(g, lib.net_center_periphery, 0, reach, *_path(spec), 0)


@register("dm_periphery")
def _periphery(g, spec, fx):
    lib = _declare()
    reach = {"require": 0, "ignore": 1}[spec.get("reachable", "require")]
    return _node_set(g, lib.net_center_periphery, 1, reach, *_path(spec), 0)


@register("dm_barycenter")
def _barycenter(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_int64 * max(g.n, 1))()
    count = ctypes.c_int64()
    totals = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    tie = {"relative": 0, "exact": 1}[spec.get("tie", "exact")]
    _check(lib.net_barycenter(g.handle, *_path(spec), tie, 0, out, ctypes.byref(count), totals,
                              err, ERR_LEN), err)
    return sorted(float(x) for x in list(out)[:count.value])


def _scalar(fn, g, *args):
    out = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(g.handle, *args, ctypes.byref(out), err, ERR_LEN), err)
    return [_mark(out.value)]


@register("dm_index")
def _index(g, spec, fx):
    lib = _declare()
    return _scalar(lib.net_distance_index, g, INDEX[spec["index"]], *_path(spec), 0)


@register("dm_local_efficiency")
def _local_efficiency(g, spec, fx):
    lib = _declare()
    return _scalar(lib.net_local_efficiency, g, int(spec.get("loop_is_neighbor", False)), 0)


@register("dm_closeness_vitality")
def _vitality(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_closeness_vitality(g.handle, *_path(spec), 0, out, err, ERR_LEN), err)
    return [_mark(x) for x in list(out)[:g.n]]
