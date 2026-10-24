#!/usr/bin/env python3
"""Definition oracle for maximum flows, minimum cuts, minimum-cost flows,
and Gomory-Hu trees.

By the max-flow min-cut theorem (Ford and Fulkerson 1956) the maximum flow
value is the least capacity of a cut, so it is found here by evaluating
the capacity of every node subset S (the arcs leaving S, or the edges with
one end in S if undirected) with numpy, on fixtures of at most MAX_NODES
nodes. The minimum cuts of a network are closed under union and
intersection, so the smallest and largest source sides are the
intersection and the union of every minimum cut's source side. Node cuts
are found by trying node sets in increasing size. Minimum-cost flows are
solved as linear programs (scipy's HiGHS), independently of any flow
algorithm. Nothing here shares code with the engine.

The flow problem of a fixture (terminals, capacities, costs, node
capacities) is the one documented in validation/bindings/flow.py.
Results that do not exist are written as -1.

    python3 validation/oracles/flow_definition_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 20
GH_MAX_NODES = 20
MISSING = -1.0


def problem(fx, multi=False):
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


def cut_capacities(fx, cap):
    """Capacity of every node subset S (bit v of the index = v in S)."""
    n = fx.n
    idx = np.arange(1 << n, dtype=np.int64)
    bits = [((idx >> v) & 1).astype(bool) for v in range(n)]
    total = np.zeros(1 << n)
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b or cap[e] == 0:
            continue
        if fx.directed:
            crossing = bits[a] & ~bits[b]
        else:
            crossing = bits[a] ^ bits[b]
        total += np.where(crossing, cap[e], 0.0)
    return idx, total


def feasible_mask(idx, sources, sinks):
    mask = np.ones(len(idx), dtype=bool)
    for v in sources:
        mask &= ((idx >> v) & 1).astype(bool)
    for v in sinks:
        mask &= ~((idx >> v) & 1).astype(bool)
    return mask


def min_cuts(fx, sources, sinks, cap):
    idx, total = cut_capacities(fx, cap)
    mask = feasible_mask(idx, sources, sinks)
    best = total[mask].min()
    tol = 1e-9 * max(1.0, abs(best))
    sides = idx[mask & (total <= best + tol)]
    union = np.bitwise_or.reduce(sides)
    inter = np.bitwise_and.reduce(sides)
    return best, int(union), int(inter)


def flags(bits, n):
    return [float((bits >> v) & 1) for v in range(n)]


def separated(fx, removed, s, t):
    """Whether no path leads from s to t once `removed` is deleted."""
    adj = [[] for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        adj[a].append(b)
        if not fx.directed:
            adj[b].append(a)
    seen = {s} | set(removed)
    stack = [s]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w == t:
                return False
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return True


def min_node_cut(fx, weights):
    s, t = 0, fx.n - 1
    adjacent = any(a == s and b == t or (not fx.directed and a == t and b == s)
                   for a, b in zip(fx.src, fx.tgt))
    if adjacent:
        return MISSING
    others = [v for v in range(fx.n) if v not in (s, t)]
    best = None
    for k in range(len(others) + 1):
        if best is not None and k * min(weights) >= best:
            break
        for comb in itertools.combinations(others, k):
            w = sum(weights[v] for v in comb)
            if (best is None or w < best) and separated(fx, comb, s, t):
                best = w
    return float(best)


def lp_max_flow_min_cost(fx, cap, cost):
    from scipy.optimize import linprog
    n, m = fx.n, len(fx.src)
    if m == 0:
        return [0.0, 0.0]
    s, t = 0, n - 1
    inc = np.zeros((n, m))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a != b:
            inc[a, e] -= 1.0
            inc[b, e] += 1.0
    # Stage 1: maximise F subject to inflow - outflow = F at t, -F at s.
    col = np.zeros((n, 1))
    col[s, 0], col[t, 0] = 1.0, -1.0
    a_eq = np.hstack([inc, col])
    c = np.zeros(m + 1)
    c[-1] = -1.0
    bounds = [(0.0, cap[e]) for e in range(m)] + [(0.0, None)]
    r1 = linprog(c, A_eq=a_eq, b_eq=np.zeros(n), bounds=bounds, method="highs")
    flow = -r1.fun
    b = np.zeros(n)
    b[s], b[t] = -flow, flow
    r2 = linprog(np.array(cost), A_eq=inc, b_eq=b, bounds=bounds[:m], method="highs")
    return [flow, r2.fun]


def gh_values(fx, cap):
    idx, total = cut_capacities(fx, cap)
    out = []
    for u in range(fx.n):
        bu = ((idx >> u) & 1).astype(bool)
        for v in range(u + 1, fx.n):
            bv = ((idx >> v) & 1).astype(bool)
            out.append(float(total[bu & ~bv].min()))
    return out


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    if fn == "gomory_hu_values":
        if fx.n > GH_MAX_NODES:
            return None
        m = len(fx.src)
        return gh_values(fx, list(fx.weights) if fx.weighted else [1.0] * m)
    if fn == "max_flow_min_cost":
        p = problem(fx)
        if p is None:
            return [MISSING, MISSING]
        return lp_max_flow_min_cost(fx, p[2], p[3])
    if fx.n > MAX_NODES:
        return None
    if fn == "min_node_cut":
        if fx.n < 2:
            return [MISSING]
        w = [float(v % 3 + 1) for v in range(fx.n)] if e.get("weighted") else [1.0] * fx.n
        return [min_node_cut(fx, w)]
    p = problem(fx, e.get("multi", False))
    if p is None:
        return [MISSING]
    s, t, cap, _ = p
    best, union, inter = min_cuts(fx, s, t, cap)
    if fn == "max_flow":
        return [best]
    if fn == "min_cut":
        return flags(union if e["side"] == "largest" else inter, fx.n)
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("flow")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", []):
                continue
            values = compute(spec, fx)
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (flow): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
