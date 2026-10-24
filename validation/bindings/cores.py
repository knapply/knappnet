"""Engine runners for k-shells, k-crusts, k-coronas, k-trusses, and onion
layers (crates/capi/src/cores.rs).

Encodings: a node or edge set is its increasing indices, one per line
(no lines for an empty set); truss numbers are edge-aligned (0 for
self-loops); onion layers are node-aligned.

Manifest engine options:
  fn: core_subset   which (core | shell | crust | corona), k (omit for
                    NetworkX's default)
  fn: truss         the truss number of every edge
  fn: k_truss       k, part (edges | nodes)
  fn: onion         every node's layer
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

WHICH = {"core": 0, "shell": 1, "crust": 2, "corona": 3}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64 = ctypes.c_int32, ctypes.c_int64
    i64p = ctypes.POINTER(i64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_core_subset.argtypes = [h, i32, i64, i32, i64p, i64p, cp, sz]
    lib.net_truss_numbers.argtypes = [h, i64p, cp, sz]
    lib.net_k_truss.argtypes = [h, i64, i64p, i64p, i64p, i64p, cp, sz]
    lib.net_onion_layers.argtypes = [h, i64p, cp, sz]
    for fn in ("net_core_subset", "net_truss_numbers", "net_k_truss", "net_onion_layers"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


@register("core_subset")
def _core_subset(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nodes = (ctypes.c_int64 * max(fx.n, 1))()
    count = ctypes.c_int64()
    k = spec.get("k")
    _check(lib.net_core_subset(g.handle, WHICH[spec["which"]], int(k or 0), int(k is not None),
                               nodes, ctypes.byref(count), err, ERR_LEN), err)
    return [nodes[i] for i in range(count.value)]


@register("truss")
def _truss(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    m = len(fx.src)
    out = (ctypes.c_int64 * max(m, 1))()
    _check(lib.net_truss_numbers(g.handle, out, err, ERR_LEN), err)
    return [out[i] for i in range(m)]


@register("k_truss")
def _k_truss(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    m = len(fx.src)
    edges = (ctypes.c_int64 * max(m, 1))()
    nodes = (ctypes.c_int64 * max(fx.n, 1))()
    ne, nn = ctypes.c_int64(), ctypes.c_int64()
    _check(lib.net_k_truss(g.handle, int(spec["k"]), edges, ctypes.byref(ne), nodes,
                           ctypes.byref(nn), err, ERR_LEN), err)
    if spec["part"] == "edges":
        return [edges[i] for i in range(ne.value)]
    return [nodes[i] for i in range(nn.value)]


@register("onion")
def _onion(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    out = (ctypes.c_int64 * max(fx.n, 1))()
    _check(lib.net_onion_layers(g.handle, out, err, ERR_LEN), err)
    return [out[i] for i in range(fx.n)]
