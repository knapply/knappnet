"""Engine runners for cut measures (crates/capi/src/cut_measures.rs).

Derived inputs, which the oracles build the same way:
  S   the nodes v with v % 2 == 0
  T   the complement of S (`t: complement`), or the nodes v with
      v % 3 == 1 (`t: mod3`, which overlaps S)

Encodings: nine values in this order: cut size, volume of S, volume of T,
normalized cut size, conductance, edge expansion, mixing expansion, node
expansion, boundary expansion; NaN is -999 and infinity -1000 (where
NetworkX divides by zero).

Manifest engine options:
  fn: cut_measures   t (complement | mod3), weighted
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    lib.net_cut_measures.argtypes = [ctypes.c_void_p, ctypes.POINTER(i64), i64,
                                     ctypes.POINTER(i64), i64, i32, ctypes.POINTER(f64),
                                     ctypes.c_char_p, ctypes.c_size_t]
    lib.net_cut_measures.restype = i32
    _declared = True
    return lib


def encode(x: float) -> float:
    if math.isnan(x):
        return -999.0
    if math.isinf(x):
        return -1000.0
    return x


def sets(fx, spec):
    s = [v for v in range(fx.n) if v % 2 == 0]
    t = None if spec.get("t", "complement") == "complement" else [v for v in range(fx.n)
                                                                  if v % 3 == 1]
    return s, t


@register("cut_measures")
def _cut_measures(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    s, t = sets(fx, spec)
    sa = (ctypes.c_int64 * max(len(s), 1))(*s)
    ta = (ctypes.c_int64 * max(len(t), 1))(*t) if t is not None else None
    out = (ctypes.c_double * 9)()
    if lib.net_cut_measures(g.handle, sa, len(s), ta, len(t) if t is not None else 0,
                            int(bool(spec.get("weighted"))), out, err, ERR_LEN):
        raise EngineError(err.value.decode())
    return [encode(out[i]) for i in range(9)]
