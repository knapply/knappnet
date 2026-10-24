"""Engine runners for edge swaps, reference graphs, the small-world
coefficients, and the rich-club coefficient (crates/capi/src/rewiring.rs).

Swaps and reference graphs are random, so the engine and NetworkX draw
different graphs; what they must share is the contract, which the runner
checks on the engine's graph and the NetworkX oracle on NetworkX's:

  double_edge_swap      (degrees kept, simple, number of edges)
  connected_swap        (degrees kept, connected, number of edges)
  directed_swap         (in-degrees kept, out-degrees kept, simple, arcs)
  random_reference,     (degrees kept, connected if the input is, number
  lattice_reference      of edges)

Each tuple is [-1] where the function fails (too few nodes or edges, no
swap possible, a disconnected input). sigma, omega, and the normalized
rich-club coefficient are means over random graphs and are compared with
a tolerance; the unnormalized rich-club coefficient is exact. NaN is -999.

Manifest engine options:
  fn: swap_contract   procedure (double | connected | directed | random |
                      lattice)
  fn: small_world     coefficient (sigma | omega), niter, nrand; the mean
                      over seeds 0 .. seeds - 1
  fn: rich_club       normalized, q, n_samples
"""

from __future__ import annotations

import ctypes
import math
from collections import Counter

from engine import ERR_LEN, load, register

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    i64p, f64p = ctypes.POINTER(i64), ctypes.POINTER(f64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_double_edge_swap.argtypes = [h, i64, i64, u64, i64p, i64p, i64p, cp, sz]
    lib.net_connected_double_edge_swap.argtypes = [h, i64, u64, i64p, i64p, i64p, i64p, cp, sz]
    lib.net_directed_edge_swap.argtypes = [h, i64, i64, u64, i64p, i64p, i64p, cp, sz]
    lib.net_random_reference.argtypes = [h, i64, i32, u64, i64p, i64p, i64p, cp, sz]
    lib.net_lattice_reference.argtypes = [h, i64, f64p, i32, u64, i64p, i64p, i64p, cp, sz]
    lib.net_small_world_sigma.argtypes = [h, i64, i64, u64, i32, f64p, cp, sz]
    lib.net_small_world_omega.argtypes = [h, i64, i64, u64, i32, f64p, cp, sz]
    lib.net_rich_club.argtypes = [h, i32, i64, i64, u64, i32, f64p, i64p, cp, sz]
    for fn in ("net_double_edge_swap", "net_connected_double_edge_swap",
               "net_directed_edge_swap", "net_random_reference", "net_lattice_reference",
               "net_small_world_sigma", "net_small_world_omega", "net_rich_club"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def simple_pairs(fx):
    if fx.directed:
        return set(zip(fx.src, fx.tgt))
    return {(min(a, b), max(a, b)) for a, b in zip(fx.src, fx.tgt)}


def degrees(n, pairs, directed):
    if directed:
        return (Counter(a for a, _ in pairs), Counter(b for _, b in pairs))
    d = Counter()
    for a, b in pairs:
        d[a] += 1
        d[b] += 1
    return d


def connected(n, pairs):
    if n == 0:
        return False
    adj = {v: set() for v in range(n)}
    for a, b in pairs:
        adj[a].add(b)
        adj[b].add(a)
    seen, stack = {0}, [0]
    while stack:
        for w in adj[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def contract(fx, procedure, out_pairs):
    """The contract tuple of a rewired graph (see the module docstring)."""
    before = simple_pairs(fx)
    simple = len(set(out_pairs)) == len(out_pairs)
    if procedure == "directed":
        di, do = degrees(fx.n, before, True)
        ai, ao = degrees(fx.n, out_pairs, True)
        return (int(di == ai), int(do == ao), int(simple), len(out_pairs))
    same = degrees(fx.n, before, False) == degrees(fx.n, out_pairs, False)
    if procedure == "double":
        return (int(same), int(simple), len(out_pairs))
    if procedure == "connected":
        return (int(same), int(connected(fx.n, out_pairs)), len(out_pairs))
    keeps = (not connected(fx.n, before)) or connected(fx.n, out_pairs)
    return (int(same), int(keeps), len(out_pairs))


@register("swap_contract")
def _swap_contract(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    m = max(len(fx.src), 1)
    f = (ctypes.c_int64 * m)()
    t = (ctypes.c_int64 * m)()
    k = ctypes.c_int64()
    p = spec["procedure"]
    if p == "double":
        st = lib.net_double_edge_swap(g.handle, 10, 1000, 7, f, t, ctypes.byref(k), err, ERR_LEN)
    elif p == "connected":
        sw = ctypes.c_int64()
        st = lib.net_connected_double_edge_swap(g.handle, 20, 7, f, t, ctypes.byref(k),
                                                ctypes.byref(sw), err, ERR_LEN)
    elif p == "directed":
        st = lib.net_directed_edge_swap(g.handle, 5, 500, 7, f, t, ctypes.byref(k), err, ERR_LEN)
    elif p == "random":
        st = lib.net_random_reference(g.handle, 2, 1, 7, f, t, ctypes.byref(k), err, ERR_LEN)
    else:
        st = lib.net_lattice_reference(g.handle, 2, None, 1, 7, f, t, ctypes.byref(k), err,
                                       ERR_LEN)
    if st:
        return [-1]
    return [contract(fx, p, [(f[i], t[i]) for i in range(k.value)])]


@register("small_world")
def _small_world(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    out = ctypes.c_double()
    fn = lib.net_small_world_sigma if spec["coefficient"] == "sigma" else lib.net_small_world_omega
    values = []
    for seed in range(int(spec.get("seeds", 1))):
        if fn(g.handle, int(spec["niter"]), int(spec["nrand"]), seed, 0, ctypes.byref(out), err,
              ERR_LEN):
            return [-1]
        values.append(out.value)
    mean = sum(values) / len(values)
    return [-999.0 if math.isnan(mean) else mean]


@register("rich_club")
def _rich_club(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    out = (ctypes.c_double * max(fx.n, 1))()
    c = ctypes.c_int64()
    if lib.net_rich_club(g.handle, int(bool(spec.get("normalized"))), int(spec.get("q", 100)),
                         int(spec.get("n_samples", 1)), 1, 0, out, ctypes.byref(c), err,
                         ERR_LEN):
        return [-1]
    return [out[i] for i in range(c.value)]
