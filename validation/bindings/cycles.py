"""Engine runners for cycle bases, minimum cycle bases, simple and
chordless cycles, girth, and find_cycle (crates/capi/src/cycles.rs).

Encodings: a cycle is one tuple per line. Sets of cycles are canonical:
rotated to start at their lowest node and, for undirected networks,
turned so the second node is the lower of its two neighbours on the
cycle. Cycle bases keep the engine's order and orientation (compared with
NetworkX exactly). An infinite girth is -2. find_cycle reports one
(edge position, reversed) pair per line, nothing without a cycle. Results
that cannot be compared exactly (bases, which are not unique) also come
as [count, valid], valid being checked here: every list is a cycle of the
network, the cycles are independent modulo 2, and there are m - n + c.

Manifest engine options:
  fn: cycle_basis        the basis; order: sorted | input
  fn: basis_check        [count, valid] of cycle_basis
  fn: min_basis_weights  each cycle's weight (compare sorted); weighted
  fn: min_basis_check    [count, valid] of minimum_cycle_basis; weighted
  fn: simple_cycles      canonical cycles; length_bound
  fn: chordless_cycles   canonical cycles; length_bound
  fn: girth              [girth]
  fn: find_cycle         (edge, reversed) pairs; orientation
  skip: [fixtures]         simple or chordless cycles are not enumerated
                           on these (too many cycles); the runner returns
                           nothing and every oracle is declared undefined
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

ORDER = {"sorted": 0, "input": 1}
ORIENTATION = {"original": 0, "reverse": 1, "ignore": 2}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    hp = ctypes.POINTER(h)
    lib.net_cycle_basis.argtypes = [h, i32, i64, hp, cp, sz]
    lib.net_minimum_cycle_basis.argtypes = [h, i32, i32, hp, cp, sz]
    lib.net_simple_cycles.argtypes = [h, i64, i64, hp, cp, sz]
    lib.net_chordless_cycles.argtypes = [h, i64, i64, hp, cp, sz]
    lib.net_girth.argtypes = [h, f64p, cp, sz]
    lib.net_find_cycle.argtypes = [h, i64, i32, i32, i64p, i32p, i64p, cp, sz]
    for fn in ("net_cycle_basis", "net_minimum_cycle_basis", "net_simple_cycles",
               "net_chordless_cycles", "net_girth", "net_find_cycle"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def _lists(lib, handle):
    from bindings_connectivity import read_lists
    return read_lists(lib, handle)


def canonical(cycle, directed):
    c = list(cycle)
    i = c.index(min(c))
    c = c[i:] + c[:i]
    if not directed and len(c) > 2 and c[1] > c[-1]:
        c = [c[0]] + c[1:][::-1]
    return tuple(float(v) for v in c)


def basis_valid(fx, cycles):
    """Whether `cycles` (node lists) is a cycle basis: each is a cycle of
    the multigraph (consecutive nodes joined, a loop for [v], two parallel
    edges for [u, v]), there are m - n + c of them, and their edge sets
    (with parallel edges told apart greedily) are independent modulo 2."""
    edges = {}
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        edges.setdefault((min(a, b), max(a, b)), []).append(e)
    parent = list(range(fx.n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(fx.src, fx.tgt):
        parent[find(a)] = find(b)
    comps = len({find(v) for v in range(fx.n)})
    if len(cycles) != len(fx.src) - fx.n + comps:
        return False
    vectors = []
    seen_short = {}
    for c in cycles:
        if len(c) <= 2:
            # The k-th loop [v] is v's k-th loop edge; the k-th 2-cycle
            # [u, v] is the first and the (k + 1)-th edge between them.
            key = (min(c), max(c))
            k = seen_short.get(key, 0)
            seen_short[key] = k + 1
            lst = edges.get(key, [])
            if len(c) == 1:
                if k >= len(lst):
                    return False
                vectors.append(1 << lst[k])
            else:
                if k + 1 >= len(lst):
                    return False
                vectors.append((1 << lst[0]) ^ (1 << lst[k + 1]))
            continue
        if len(set(c)) != len(c):
            return False
        bits = 0
        for a, b in zip(c, c[1:] + c[:1]):
            key = (min(a, b), max(a, b))
            if key not in edges:
                return False
            bits ^= 1 << edges[key][0]
        vectors.append(bits)
    # Gaussian elimination over GF(2).
    basis = {}
    for v in vectors:
        while v:
            top = v.bit_length() - 1
            if top not in basis:
                basis[top] = v
                break
            v ^= basis[top]
        else:
            return False
    return True


@register("cycle_basis")
def _run_cycle_basis(g, spec, fx):
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_cycle_basis(g.handle, ORDER[spec.get("order", "input")], -1, ctypes.byref(hd),
                               err, ERR_LEN), err)
    lists, _ = _lists(lib, hd)
    return [tuple(float(v) for v in c) for c in lists]


@register("basis_check")
def _run_basis_check(g, spec, fx):
    cycles = [[int(v) for v in c] for c in _run_cycle_basis(g, spec, fx)]
    return [float(len(cycles)), 1.0 if basis_valid(fx, cycles) else 0.0]


def _min_basis(g, spec, fx):
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    weighted = spec.get("weighted", False) and fx.weighted
    _check(lib.net_minimum_cycle_basis(g.handle, int(weighted), 0, ctypes.byref(hd), err,
                                       ERR_LEN), err)
    return _lists(lib, hd)


@register("min_basis_weights")
def _run_min_basis_weights(g, spec, fx):
    return _min_basis(g, spec, fx)[1]


@register("min_basis_check")
def _run_min_basis_check(g, spec, fx):
    cycles, _ = _min_basis(g, spec, fx)
    return [float(len(cycles)), 1.0 if basis_valid(fx, cycles) else 0.0]


def _cycles(fn, g, spec, fx):
    if fx.name in spec.get("skip", []):
        return []
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    bound = spec.get("length_bound", -1)
    _check(fn(g.handle, bound, 0, ctypes.byref(hd), err, ERR_LEN), err)
    lists, _ = _lists(lib, hd)
    return [canonical(c, fx.directed) for c in lists]


@register("simple_cycles")
def _run_simple(g, spec, fx):
    return _cycles(_declare().net_simple_cycles, g, spec, fx)


@register("chordless_cycles")
def _run_chordless(g, spec, fx):
    return _cycles(_declare().net_chordless_cycles, g, spec, fx)


@register("girth")
def _run_girth(g, spec, fx):
    lib = _declare()
    out = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_girth(g.handle, ctypes.byref(out), err, ERR_LEN), err)
    return [-2.0 if out.value == float("inf") else out.value]


@register("find_cycle")
def _run_find_cycle(g, spec, fx):
    lib = _declare()
    m = max(len(fx.src), 1)
    edges = (ctypes.c_int64 * m)()
    rev = (ctypes.c_int32 * m)()
    size = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_find_cycle(g.handle, -1, ORIENTATION[spec.get("orientation", "original")],
                              ORDER["input"], edges, rev, ctypes.byref(size), err, ERR_LEN), err)
    return [(float(edges[i]), float(rev[i])) for i in range(size.value)]
