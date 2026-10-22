"""Engine runners for cohesive subgroups: cores, cliques, and communities."""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

_declared = False


def _lib():
    global _declared
    lib = load()
    if _declared:
        return lib
    i64p = ctypes.POINTER(ctypes.c_int64)
    f64p = ctypes.POINTER(ctypes.c_double)
    handle = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_coreness.argtypes = [handle, ctypes.c_int32, f64p, cp, sz]
    lib.net_coreness.restype = ctypes.c_int32
    lib.net_max_cliques.argtypes = [handle, ctypes.c_int64, ctypes.c_int32,
                                    ctypes.POINTER(handle), cp, sz]
    lib.net_max_cliques.restype = ctypes.c_int32
    lib.net_groups_count.argtypes = [handle]
    lib.net_groups_count.restype = ctypes.c_int64
    lib.net_groups_total.argtypes = [handle]
    lib.net_groups_total.restype = ctypes.c_int64
    lib.net_groups_copy.argtypes = [handle, i64p, i64p, cp, sz]
    lib.net_groups_copy.restype = ctypes.c_int32
    lib.net_groups_free.argtypes = [handle]
    lib.net_groups_free.restype = None
    lib.net_clique_number.argtypes = [handle, ctypes.c_int32, i64p, cp, sz]
    lib.net_clique_number.restype = ctypes.c_int32
    lib.net_modularity.argtypes = [handle, i64p, ctypes.c_int64, ctypes.c_double,
                                   ctypes.c_int32, f64p, cp, sz]
    lib.net_modularity.restype = ctypes.c_int32
    lib.net_louvain.argtypes = [handle, ctypes.c_double, ctypes.c_int32, ctypes.c_uint64,
                                i64p, f64p, cp, sz]
    lib.net_louvain.restype = ctypes.c_int32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def block_partition(n: int, k: int = 3) -> list[int]:
    """The fixed partition used to validate modularity: nodes split into k
    contiguous blocks by index, label floor(k * v / n)."""
    return [k * v // n for v in range(n)]


def coreness(g, mode: str = "all") -> list[float]:
    lib = _lib()
    out = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_coreness(g.handle, {"out": 1, "in": 2, "all": 3}[mode], out, err, ERR_LEN),
           err)
    return list(out)[:g.n]


def max_cliques(g, min_size: int = 0, threads: int = 0) -> list[tuple[int, ...]]:
    lib = _lib()
    err = ctypes.create_string_buffer(ERR_LEN)
    h = ctypes.c_void_p()
    _check(lib.net_max_cliques(g.handle, min_size, threads, ctypes.byref(h), err, ERR_LEN), err)
    try:
        count = lib.net_groups_count(h)
        total = lib.net_groups_total(h)
        offsets = (ctypes.c_int64 * (count + 1))()
        members = (ctypes.c_int64 * max(total, 1))()
        _check(lib.net_groups_copy(h, offsets, members, err, ERR_LEN), err)
    finally:
        lib.net_groups_free(h)
    return [tuple(members[offsets[i]:offsets[i + 1]]) for i in range(count)]


def clique_number(g, threads: int = 0) -> int:
    lib = _lib()
    err = ctypes.create_string_buffer(ERR_LEN)
    out = ctypes.c_int64()
    _check(lib.net_clique_number(g.handle, threads, ctypes.byref(out), err, ERR_LEN), err)
    return out.value


def modularity(g, membership, resolution: float = 1.0, weighted: bool = False) -> float:
    lib = _lib()
    err = ctypes.create_string_buffer(ERR_LEN)
    memb = (ctypes.c_int64 * max(len(membership), 1))(*membership)
    out = ctypes.c_double()
    _check(lib.net_modularity(g.handle, memb, len(membership), resolution, int(weighted),
                              ctypes.byref(out), err, ERR_LEN), err)
    return out.value


def louvain(g, resolution: float = 1.0, weighted: bool = False, seed: int = 1):
    """Returns (membership, modularity)."""
    lib = _lib()
    err = ctypes.create_string_buffer(ERR_LEN)
    memb = (ctypes.c_int64 * max(g.n, 1))()
    q = ctypes.c_double()
    _check(lib.net_louvain(g.handle, resolution, int(weighted), seed, memb, ctypes.byref(q),
                           err, ERR_LEN), err)
    return list(memb)[:g.n], q.value


@register("coreness")
def _coreness(g, spec, fx):
    return coreness(g, spec["mode"])


@register("max_cliques")
def _max_cliques(g, spec, fx):
    return max_cliques(g, spec.get("min_size", 0))


@register("clique_number")
def _clique_number(g, spec, fx):
    return [clique_number(g)]


@register("modularity")
def _modularity(g, spec, fx):
    assert spec["partition"] == "blocks3"
    return [modularity(g, block_partition(g.n), spec.get("resolution", 1.0),
                       spec.get("weighted", False))]


@register("louvain_modularity")
def _louvain_modularity(g, spec, fx):
    _, q = louvain(g, spec.get("resolution", 1.0), spec.get("weighted", False),
                   spec.get("seed", 1))
    return [q]


@register("louvain_quality")
def _louvain_quality(g, spec, fx):
    seeds = range(1, spec["seeds"] + 1)
    qs = [louvain(g, spec.get("resolution", 1.0), spec.get("weighted", False), s)[1]
          for s in seeds]
    return [sum(qs) / len(qs)]
