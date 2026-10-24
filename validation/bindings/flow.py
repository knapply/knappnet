"""Engine runners for maximum flows, minimum cuts, minimum-cost flows, and
Gomory-Hu trees (crates/capi/src/flow.rs).

Every flow problem on a fixture uses the same terminals and data, which
the oracles reproduce (see flow_problem() below):
  source 0 and sink n - 1 (sources 0, 1 and sinks n - 2, n - 1 for the
  `multi` variants); edge capacities the fixture's weights, or 1 on an
  unweighted fixture; edge costs ((7 e) mod 11) - 3 for edge e in input
  order; node capacities (v mod 3) + 1 where a variant weights nodes.
Results that do not exist (too few nodes for the terminals, or adjacent
terminals for a node cut) are written as -1 by the runner and every
oracle.

Manifest engine options:
  fn: max_flow             [value]
  fn: max_flow_check       [value, valid]: the value recomputed from the
                           edge flows, and 1 if those flows respect every
                           capacity and conserve flow at every other node
  fn: min_cut              node-aligned 1/0 source side; side: smallest |
                           largest
  fn: min_node_cut         [value]; weighted: node capacities as above
  fn: max_flow_min_cost    [value, cost]
  fn: gomory_hu            tree edges (u, v, weight) with u < v
  fn: gomory_hu_values     minimum cut value of every pair u < v, read
                           off the tree
  multi: true              several sources and sinks
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

MISSING = -1.0
SIDE = {"smallest": 0, "largest": 1}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_max_flow.argtypes = [h, i64p, i64, i64p, i64, f64p, f64p, f64p, cp, sz]
    lib.net_min_cut.argtypes = [h, i64p, i64, i64p, i64, f64p, i32, f64p, i32p, i32p, cp, sz]
    lib.net_min_node_cut_capacitated.argtypes = [h, i64p, i64, i64p, i64, f64p, f64p, i32,
                                                 f64p, i32p, i32p, cp, sz]
    lib.net_min_cost_flow.argtypes = [h, f64p, f64p, f64p, f64p, f64p, cp, sz]
    lib.net_max_flow_min_cost.argtypes = [h, i64, i64, f64p, f64p, f64p, f64p, f64p, cp, sz]
    lib.net_gomory_hu_tree.argtypes = [h, f64p, i64p, f64p, cp, sz]
    for fn in ("net_max_flow", "net_min_cut", "net_min_node_cut_capacitated",
               "net_min_cost_flow", "net_max_flow_min_cost", "net_gomory_hu_tree"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def flow_problem(fx, multi=False):
    """(sources, sinks, capacities, costs) of the fixture's flow problem, or
    None when it has too few nodes."""
    if multi:
        if fx.n < 4:
            return None
        s, t = [0, 1], [fx.n - 2, fx.n - 1]
    else:
        if fx.n < 2:
            return None
        s, t = [0], [fx.n - 1]
    m = len(fx.src)
    cap = list(fx.weights) if fx.weighted else [1.0] * m
    cost = [float((7 * e) % 11 - 3) for e in range(m)]
    return s, t, cap, cost


def node_capacities(n):
    return [float(v % 3 + 1) for v in range(n)]


def _arr(kind, values):
    return (kind * max(len(values), 1))(*values)


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def _max_flow(g, s, t, cap):
    lib = _declare()
    m = len(cap)
    value = ctypes.c_double()
    flows = (ctypes.c_double * max(m, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_max_flow(g.handle, _arr(ctypes.c_int64, s), len(s), _arr(ctypes.c_int64, t),
                            len(t), _arr(ctypes.c_double, cap), ctypes.byref(value), flows,
                            err, ERR_LEN), err)
    return value.value, list(flows)[:m]


@register("max_flow")
def _run_max_flow(g, spec, fx):
    p = flow_problem(fx, spec.get("multi", False))
    if p is None:
        return [MISSING]
    s, t, cap, _ = p
    return [_max_flow(g, s, t, cap)[0]]


@register("max_flow_check")
def _run_max_flow_check(g, spec, fx):
    p = flow_problem(fx, spec.get("multi", False))
    if p is None:
        return [MISSING, MISSING]
    s, t, cap, _ = p
    value, flows = _max_flow(g, s, t, cap)
    tol = 1e-9 * max([1.0] + [c for c in cap if math.isfinite(c)])
    net = [0.0] * fx.n
    ok = True
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        f = flows[e]
        if a == b:
            ok &= f == 0.0
            continue
        if fx.directed:
            ok &= -tol <= f <= cap[e] + tol
        else:
            ok &= abs(f) <= cap[e] + tol
        net[a] -= f
        net[b] += f
    terminals = set(s) | set(t)
    ok &= all(abs(net[v]) <= tol * len(fx.src) for v in range(fx.n) if v not in terminals)
    out = -sum(net[v] for v in s)
    return [out, 1.0 if ok else 0.0]


@register("min_cut")
def _run_min_cut(g, spec, fx):
    p = flow_problem(fx, spec.get("multi", False))
    if p is None:
        return [MISSING]
    s, t, cap, _ = p
    lib = _declare()
    value = ctypes.c_double()
    side = (ctypes.c_int32 * max(fx.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_min_cut(g.handle, _arr(ctypes.c_int64, s), len(s), _arr(ctypes.c_int64, t),
                           len(t), _arr(ctypes.c_double, cap), SIDE[spec["side"]],
                           ctypes.byref(value), side, None, err, ERR_LEN), err)
    return [float(x) for x in list(side)[:fx.n]]


@register("min_node_cut")
def _run_min_node_cut(g, spec, fx):
    p = flow_problem(fx)
    if p is None:
        return [MISSING]
    s, t, _, _ = p
    lib = _declare()
    value = ctypes.c_double()
    nodes = _arr(ctypes.c_double, node_capacities(fx.n)) if spec.get("weighted") else None
    err = ctypes.create_string_buffer(ERR_LEN)
    status = lib.net_min_node_cut_capacitated(
        g.handle, _arr(ctypes.c_int64, s), len(s), _arr(ctypes.c_int64, t), len(t), nodes, None,
        0, ctypes.byref(value), None, None, err, ERR_LEN)
    if status:
        if "no node cut" in err.value.decode():
            return [MISSING]
        raise EngineError(err.value.decode())
    return [value.value]


@register("max_flow_min_cost")
def _run_max_flow_min_cost(g, spec, fx):
    p = flow_problem(fx)
    if p is None:
        return [MISSING, MISSING]
    s, t, cap, cost = p
    lib = _declare()
    value, total = ctypes.c_double(), ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_max_flow_min_cost(g.handle, s[0], t[0], _arr(ctypes.c_double, cap),
                                     _arr(ctypes.c_double, cost), ctypes.byref(value),
                                     ctypes.byref(total), None, err, ERR_LEN), err)
    return [value.value, total.value]


def _gomory_hu(g, fx):
    lib = _declare()
    m = len(fx.src)
    cap = list(fx.weights) if fx.weighted else [1.0] * m
    parent = (ctypes.c_int64 * max(fx.n, 1))()
    weight = (ctypes.c_double * max(fx.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_gomory_hu_tree(g.handle, _arr(ctypes.c_double, cap), parent, weight, err,
                                  ERR_LEN), err)
    return list(parent)[:fx.n], list(weight)[:fx.n]


@register("gomory_hu")
def _run_gomory_hu(g, spec, fx):
    parent, weight = _gomory_hu(g, fx)
    return [(min(v, parent[v]), max(v, parent[v]), weight[v]) for v in range(1, fx.n)]


@register("gomory_hu_values")
def _run_gomory_hu_values(g, spec, fx):
    parent, weight = _gomory_hu(g, fx)
    adj = {v: [] for v in range(fx.n)}
    for v in range(1, fx.n):
        adj[v].append((parent[v], weight[v]))
        adj[parent[v]].append((v, weight[v]))
    out = []
    for u in range(fx.n):
        low = {u: math.inf}
        stack = [u]
        while stack:
            x = stack.pop()
            for y, w in adj[x]:
                if y not in low:
                    low[y] = min(low[x], w)
                    stack.append(y)
        out += [low[v] for v in range(u + 1, fx.n)]
    return out
