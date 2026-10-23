#!/usr/bin/env python3
"""Golden values for ego networks and two-mode networks, computed directly
from their definitions.

This oracle shares no code or algorithmic strategy with the engine.
Neighborhoods come from repeated expansion of adjacency sets; ties among
alters are checked pair by pair; components come from the boolean
transitive closure; ego betweenness and two-mode betweenness count
shortest paths pair by pair from Floyd-Warshall distances with exact
integers and fractions (the engine uses Brandes' accumulation over
breadth-first searches). Two-mode scores are evaluated on the incidence
matrix: projections as exact matrix products, eigenvector scores from a
dense singular value decomposition with numpy (the engine uses power
iteration), and normalizations from the formulas of Borgatti and Everett
(1997).

Undefined values (NaN) are written as -1, as the engine runner does.

    python3 validation/oracles/twomode_definition_oracle.py
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, FIXTURES, applies, fixture_names, load_area,
                    load_fixture, write_values, write_version)  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80
UNDEFINED = -1.0
INF = None


def manifest() -> dict:
    return load_area("twomode")


def node_types(fx):
    return json.loads((FIXTURES / fx.name / "meta.json").read_text()).get("type")


# ---- Adjacency ----------------------------------------------------------------

def arc_set(fx):
    """Distinct arcs between distinct nodes; both directions if undirected."""
    arcs = set()
    for a, b in zip(fx.src, fx.tgt):
        if a == b:
            continue
        arcs.add((a, b))
        if not fx.directed:
            arcs.add((b, a))
    return arcs


def step_sets(fx, mode):
    """For each node, the nodes one step away in the direction of `mode`."""
    nxt = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if not fx.directed or mode in ("out", "all"):
            nxt[a].add(b)
        if not fx.directed or mode in ("in", "all"):
            nxt[b].add(a)
    return nxt


def within(fx, ego, order, mode):
    """{node: steps} for every node within `order` steps of ego."""
    nxt = step_sets(fx, mode)
    steps = {ego: 0}
    frontier = {ego}
    for k in range(1, order + 1):
        new = set()
        for v in frontier:
            new |= nxt[v]
        new -= set(steps)
        for v in new:
            steps[v] = k
        frontier = new
    return steps


# ---- Ego networks -------------------------------------------------------------

def components_count(nodes, arcs):
    """Weak components among `nodes` by the transitive closure of the
    symmetrized ties."""
    nodes = sorted(nodes)
    reach = {v: {v} for v in nodes}
    for (a, b) in arcs:
        if a in reach and b in reach:
            reach[a].add(b)
            reach[b].add(a)
    changed = True
    while changed:
        changed = False
        for v in nodes:
            grown = set().union(*(reach[u] for u in reach[v]))
            if grown != reach[v]:
                reach[v] = grown
                changed = True
    return len({frozenset(reach[v]) for v in nodes})


def shortest_path_counts(nodes, arcs):
    """Floyd-Warshall distances and exact numbers of shortest paths among
    `nodes` using only `arcs` between them."""
    idx = {v: i for i, v in enumerate(nodes)}
    k = len(nodes)
    d = [[INF] * k for _ in range(k)]
    for i in range(k):
        d[i][i] = 0
    for (a, b) in arcs:
        if a in idx and b in idx:
            d[idx[a]][idx[b]] = 1
    for m in range(k):
        for i in range(k):
            if d[i][m] is INF:
                continue
            for j in range(k):
                if d[m][j] is INF:
                    continue
                if d[i][j] is INF or d[i][m] + d[m][j] < d[i][j]:
                    d[i][j] = d[i][m] + d[m][j]
    succ = [[idx[b] for (a, b) in arcs if a == v and b in idx] for v in nodes]
    sigma = [[0] * k for _ in range(k)]
    for s in range(k):
        sigma[s][s] = 1
        order = sorted((t for t in range(k) if d[s][t] is not INF), key=lambda t: d[s][t])
        for u in order:
            for t in succ[u]:
                if d[s][t] is not INF and d[s][u] + 1 == d[s][t]:
                    sigma[s][t] += sigma[s][u]
    return idx, d, sigma


def node_betweenness(nodes, arcs, v, directed):
    """Betweenness of v: sum over pairs s, t (neither v) with a path of
    sigma_st(v) / sigma_st, exact."""
    idx, d, sigma = shortest_path_counts(nodes, arcs)
    x = idx[v]
    total = Fraction(0)
    k = len(nodes)
    for s in range(k):
        for t in range(k):
            if s == t or s == x or t == x or d[s][t] is INF:
                continue
            if d[s][x] is INF or d[x][t] is INF or d[s][x] + d[x][t] != d[s][t]:
                continue
            total += Fraction(sigma[s][x] * sigma[x][t], sigma[s][t])
    return total if directed else total / 2


def ego_measure(fx, e):
    arcs = arc_set(fx)
    order, mode = e.get("order", 1), e.get("mode", "all")
    measure = e["measure"]
    out = []
    for ego in range(fx.n):
        steps = within(fx, ego, order, mode)
        alters = [v for v in steps if v != ego]
        k = len(alters)
        among = {(a, b) for (a, b) in arcs if a in steps and b in steps and ego not in (a, b)}
        ties = len(among) if fx.directed else len(among) // 2
        pairs = k * (k - 1) if fx.directed else k * (k - 1) // 2
        if measure == "size":
            out.append(float(k))
        elif measure == "ties":
            out.append(float(ties))
        elif measure == "pairs":
            out.append(float(pairs))
        elif measure == "density":
            out.append(float(Fraction(ties, pairs)) if pairs else UNDEFINED)
        elif measure == "components":
            out.append(float(components_count(alters, among)))
        elif measure == "betweenness":
            members = sorted(steps)
            sub = {(a, b) for (a, b) in arcs if a in steps and b in steps}
            out.append(float(node_betweenness(members, sub, ego, fx.directed)))
        else:
            raise ValueError(measure)
    return out


# ---- Two-mode networks --------------------------------------------------------

def incidence(fx, types, weighted=False):
    """Rows: first-mode nodes in node order; columns: second-mode nodes.
    Binary (distinct ties) or valued (weights of parallel edges summed)."""
    rows = [v for v in range(fx.n) if types[v] == 0]
    cols = [v for v in range(fx.n) if types[v] == 1]
    r = {v: i for i, v in enumerate(rows)}
    c = {v: j for j, v in enumerate(cols)}
    b = [[Fraction(0)] * len(cols) for _ in rows]
    for e, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if types[s] == types[t]:
            raise ValueError(f"edge {e} joins one mode")
        u, v = (s, t) if types[s] == 0 else (t, s)
        w = Fraction(fx.weights[e]) if weighted else Fraction(1)
        if weighted:
            b[r[u]][c[v]] += w
        else:
            b[r[u]][c[v]] = Fraction(1)
    return rows, cols, b


def bipartite_arcs(fx):
    arcs = set()
    for a, b in zip(fx.src, fx.tgt):
        arcs.add((a, b))
        arcs.add((b, a))
    return arcs


def b_max(own, other):
    s, t = divmod(own - 1, other)
    return Fraction(other ** 2 * (s + 1) ** 2 + other * (s + 1) * (2 * t - s - 1)
                    - t * (2 * s - t + 3), 2)


def twomode_measure(fx, types, e):
    measure = e["measure"]
    rows, cols, b = incidence(fx, types)
    n1, n2 = len(rows), len(cols)
    size = {0: n1, 1: n2}
    if measure == "density":
        ties = sum(1 for row in b for x in row if x)
        return [float(Fraction(ties, n1 * n2))] if n1 * n2 else [UNDEFINED]
    deg = {}
    for i, v in enumerate(rows):
        deg[v] = sum(1 for x in b[i] if x)
    for j, v in enumerate(cols):
        deg[v] = sum(1 for i in range(n1) if b[i][j])
    if measure == "degree":
        if e.get("normalized", True):
            return [float(Fraction(deg[v], size[1 - types[v]])) for v in range(fx.n)]
        return [float(deg[v]) for v in range(fx.n)]
    arcs = bipartite_arcs(fx)
    nodes = list(range(fx.n))
    if measure == "closeness":
        idx, d, _ = shortest_path_counts(nodes, arcs)
        out = []
        for v in nodes:
            reach = [d[v][u] for u in nodes if u != v and d[v][u] is not INF]
            r, total = len(reach), sum(reach)
            if r == 0:
                out.append(0.0 if e.get("unreachable") == "wasserman_faust" else UNDEFINED)
                continue
            own, other = size[types[v]], size[1 - types[v]]
            c = (Fraction(other + 2 * (own - 1), total) if e.get("normalized", True)
                 else Fraction(1, total))
            if e.get("unreachable") == "wasserman_faust":
                c *= Fraction(r, fx.n - 1)
            out.append(float(c))
        return out
    if measure == "betweenness":
        out = []
        for v in nodes:
            raw = node_betweenness(nodes, arcs, v, directed=False)
            if e.get("normalized", True):
                m = b_max(size[types[v]], size[1 - types[v]])
                out.append(float(raw / m) if m > 0 else UNDEFINED)
            else:
                out.append(float(raw))
        return out
    if measure == "eigenvector":
        a = np.array([[float(x) for x in row] for row in b])
        u, _, vt = np.linalg.svd(a)
        left, right = np.abs(u[:, 0]), np.abs(vt[0, :])
        left, right = left / np.linalg.norm(left), right / np.linalg.norm(right)
        out = [0.0] * fx.n
        for i, v in enumerate(rows):
            out[v] = float(left[i])
        for j, v in enumerate(cols):
            out[v] = float(right[j])
        return out
    raise ValueError(measure)


def projection(fx, types, e):
    weight = e.get("weight", "count")
    rows, cols, b = incidence(fx, types)
    if e["second"]:
        nodes = cols
        b = [list(col) for col in zip(*b)]
    else:
        nodes = rows
    valued = None
    if weight == "product":
        r2, c2, w = incidence(fx, types, weighted=True)
        valued = [list(col) for col in zip(*w)] if e["second"] else w
    partners_of = [sum(1 for i in range(len(b)) if b[i][j]) for j in range(len(b[0]) if b else 0)]
    out = []
    for i in range(len(nodes)):
        for j in range(i if e.get("loops") else i + 1, len(nodes)):
            shared = [k for k in range(len(b[i])) if b[i][k] and b[j][k]]
            if not shared:
                continue
            if weight == "count":
                w_ij = Fraction(len(shared))
            elif weight == "newman":
                w_ij = sum((Fraction(1, partners_of[k] - 1) for k in shared
                            if partners_of[k] > 1), Fraction(0))
            else:
                w_ij = sum((valued[i][k] * valued[j][k] for k in shared), Fraction(0))
            out.append((nodes[i], nodes[j], float(w_ij)))
    return sorted(out)


# ---- Driver -------------------------------------------------------------------

def compute(fx, e):
    fn = e["fn"]
    if fn == "ego_size":
        return [float(sum(1 for v, k in within(fx, ego, e.get("order", 1), e.get("mode", "all"))
                          .items() if k >= e.get("mindist", 0)))
                for ego in range(fx.n)]
    if fn == "ego_members":
        out = []
        for ego in range(fx.n):
            steps = within(fx, ego, e.get("order", 1), e.get("mode", "all"))
            out.append(tuple(float(v) for v in sorted(v for v, k in steps.items()
                                                      if k >= e.get("mindist", 0))))
        return out
    if fn == "ego_summary":
        return ego_measure(fx, e)
    types = node_types(fx)
    if types is None:
        return None
    if fn == "twomode":
        return twomode_measure(fx, types, e)
    if fn == "project":
        return projection(fx, types, e)
    return None


def main() -> None:
    spec = manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (two-mode and ego): wrote {count} golden files")


if __name__ == "__main__":
    main()
