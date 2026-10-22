"""Engine runners for brokerage and positions: Burt's structural holes,
Gould-Fernandez brokerage, structural equivalence, and CONCOR.

Golden-file conventions shared with the oracles in
validation/oracles/positions_*_oracle.py and validation/oracles/r_positions.R:

- An undefined value (NaN, such as the constraint of an isolate) is written
  as an empty line and read back as an empty group, so undefined matches
  undefined and nothing else.
- Brokerage writes one line per node with six counts: coordinator,
  itinerant, representative, gatekeeper, liaison, total.
- Structural equivalence writes one line per row of the n-by-n matrix.
- Fixtures have no group attribute, so brokerage uses a deterministic
  grouping named by the manifest (`groups: mod3` is node index modulo 3).
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

_declared = False


def _lib():
    global _declared
    lib = load()
    if not _declared:
        f64p = ctypes.POINTER(ctypes.c_double)
        i64p = ctypes.POINTER(ctypes.c_int64)
        i32p = ctypes.POINTER(ctypes.c_int32)
        h = ctypes.c_void_p
        i32 = ctypes.c_int32
        lib.net_structural_holes.argtypes = [h, i32, i32, i32, f64p, ctypes.c_char_p,
                                             ctypes.c_size_t]
        lib.net_structural_holes.restype = i32
        lib.net_brokerage.argtypes = [h, i64p, ctypes.c_int64, i32, f64p, ctypes.c_char_p,
                                      ctypes.c_size_t]
        lib.net_brokerage.restype = i32
        lib.net_structural_equivalence.argtypes = [h, i32, i32, i32, f64p, ctypes.c_char_p,
                                                   ctypes.c_size_t]
        lib.net_structural_equivalence.restype = i32
        lib.net_concor.argtypes = [h, i32, i32, ctypes.c_double, i32, i32, i32p, i32p,
                                   ctypes.c_char_p, ctypes.c_size_t]
        lib.net_concor.restype = i32
        _declared = True
    return lib


def groups_for(rule: str, n: int) -> list[int]:
    """Deterministic node groups for fixtures without a group attribute."""
    if rule == "mod3":
        return [v % 3 for v in range(n)]
    raise ValueError(f"unknown grouping rule {rule!r}")


def undefined_as_empty(values):
    return [() if isinstance(v, float) and math.isnan(v) else v for v in values]


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _holes(g, spec, code: int):
    lib = _lib()
    out = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_structural_holes(g.handle, code, int(spec["weighted"]), 0, out, err,
                                    ERR_LEN), err)
    return undefined_as_empty(list(out)[:g.n])


@register("constraint")
def _constraint(g, spec, fx):
    return _holes(g, spec, 0)


@register("effective_size")
def _effective_size(g, spec, fx):
    return _holes(g, spec, 1)


@register("efficiency")
def _efficiency(g, spec, fx):
    return _holes(g, spec, 2)


@register("brokerage")
def _brokerage(g, spec, fx):
    lib = _lib()
    n = g.n
    groups = (ctypes.c_int64 * max(n, 1))(*groups_for(spec["groups"], n))
    out = (ctypes.c_double * max(6 * n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_brokerage(g.handle, groups, n, 0, out, err, ERR_LEN), err)
    return [tuple(out[r + n * c] for c in range(6)) for r in range(n)]


METHODS = {"euclidean": 0, "correlation": 1, "hamming": 2}


@register("structural_equivalence")
def _structural_equivalence(g, spec, fx):
    lib = _lib()
    n = g.n
    out = (ctypes.c_double * max(n * n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_structural_equivalence(g.handle, METHODS[spec["method"]],
                                          int(spec.get("weighted", False)), 0, out, err,
                                          ERR_LEN), err)
    return [tuple(out[r + n * c] for c in range(n)) for r in range(n)]


@register("concor")
def _concor(g, spec, fx):
    lib = _lib()
    n = g.n
    out = (ctypes.c_int32 * max(n, 1))()
    diagnostics = (ctypes.c_int32 * 2)()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_concor(g.handle, spec["depth"], int(spec.get("weighted", False)),
                          spec.get("tolerance", 1e-7), spec.get("max_iter", 100), 0, out,
                          diagnostics, err, ERR_LEN), err)
    return [float(x) for x in list(out)[:n]]
