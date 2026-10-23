#!/usr/bin/env python3
"""Golden values for key players, group centrality, and flow betweenness,
computed directly from their definitions.

This oracle shares no code or algorithmic strategy with the engine:
- key player criteria are evaluated set by set from distance matrices
  (breadth-first search per node, or Floyd-Warshall in exact arithmetic for
  weighted networks), and optima by enumerating every set; the engine uses
  union-find, masked searches, and a parallel enumeration;
- group betweenness counts shortest paths pair by pair from the distance
  matrix, with exact integers or fractions;
- flow betweenness runs Edmonds-Karp (Edmonds and Karp 1972,
  `edmonds1972theoretical`) in exact arithmetic for every pair, and again for every pair with every other node deleted, with no
  shortcuts (the engine uses Dinic's algorithm and skips nodes that carry
  no flow).

Undefined values (NaN) are written as -1, as the engine runner does.

    python3 validation/oracles/resilience_definition_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (GOLDEN, applies, fixture_names, load_area, load_fixture,
                    write_values)  # noqa: E402

ORACLE = "definition"
UNDEFINED = -1.0
INF = None


def manifest() -> dict:
    return load_area("resilience")


def group_members(n):
    return list(range(min(3, max(1, n - 2))))


def exact(w):
    """Integers stay integers; other weights become exact fractions."""
    f = Fraction(w)
    return int(f) if f.denominator == 1 else f


def arcs(fx, mode="out", weighted=False):
    """(tail, head, length) for every arc along which paths run, loops
    dropped. Undirected edges and mode 'all' give both directions."""
    out = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        w = exact(fx.weights[i]) if weighted else 1
        if not fx.directed or mode == "all":
            out += [(a, b, w), (b, a, w)]
        elif mode == "out":
            out.append((a, b, w))
        else:
            out.append((b, a, w))
    return out


def distance_matrix(n, arc_list, removed=frozenset(), weighted=False):
    """d[i][j], or INF. Breadth-first search per node when unweighted,
    Floyd-Warshall otherwise."""
    alive = [v for v in range(n) if v not in removed]
    if not weighted:
        nbrs = [[] for _ in range(n)]
        for a, b, _ in arc_list:
            if a not in removed and b not in removed:
                nbrs[a].append(b)
        d = [[INF] * n for _ in range(n)]
        for s in alive:
            d[s][s] = 0
            q = deque([s])
            while q:
                v = q.popleft()
                for w in nbrs[v]:
                    if d[s][w] is INF:
                        d[s][w] = d[s][v] + 1
                        q.append(w)
        return d
    d = [[INF] * n for _ in range(n)]
    for v in alive:
        d[v][v] = 0
    for a, b, w in arc_list:
        if a in removed or b in removed:
            continue
        if d[a][b] is INF or w < d[a][b]:
            d[a][b] = w
    for k in alive:
        dk = d[k]
        for i in alive:
            dik = d[i][k]
            if dik is INF:
                continue
            di = d[i]
            for j in alive:
                if dk[j] is INF:
                    continue
                alt = dik + dk[j]
                if di[j] is INF or alt < di[j]:
                    di[j] = alt
    return d


# ---- Key players -----------------------------------------------------------

def kpp(fx, spec, members, cache):
    n = fx.n
    measure = spec["measure"]
    mode = spec.get("mode", "out")
    weighted = bool(spec.get("weighted", False))
    arc_list = arcs(fx, mode, weighted)
    if measure in ("fragmentation", "distance_fragmentation"):
        removed = frozenset(members)
        rest = [v for v in range(n) if v not in removed]
        if len(rest) < 2:
            return UNDEFINED
        d = distance_matrix(n, arc_list, removed, weighted)
        total = Fraction(0)
        for i in rest:
            for j in rest:
                if i == j or d[i][j] is INF:
                    continue
                total += 1 if measure == "fragmentation" else Fraction(1) / d[i][j]
        pairs = len(rest) * (len(rest) - 1)
        return float(1 - total / pairs)
    key = (mode, weighted)
    if key not in cache:
        cache[key] = distance_matrix(n, arc_list, frozenset(), weighted)
    d = cache[key]
    near = []
    for j in range(n):
        if j in members:
            continue
        ds = [d[i][j] for i in members if d[i][j] is not INF]
        near.append(min(ds) if ds else INF)
    if measure == "m_reach":
        m = spec.get("m", 1)
        return float(sum(1 for x in near if x is not INF and x <= m))
    total = Fraction(len(members)) + sum((Fraction(1) / x for x in near if x is not INF),
                                         Fraction(0))
    return float(total / n)


def kpp_single(fx, spec):
    cache = {}
    return [kpp(fx, spec, [v], cache) for v in range(fx.n)]


def kpp_optimum(fx, spec):
    k = int(spec["k"])
    if fx.n < k:
        return [UNDEFINED]
    cache = {}
    best = None
    for s in itertools.combinations(range(fx.n), k):
        x = kpp(fx, spec, list(s), cache)
        if x == UNDEFINED:
            continue
        best = x if best is None else max(best, x)
    return [UNDEFINED if best is None else best]


# ---- Group centrality --------------------------------------------------------

def group(fx, spec):
    n = fx.n
    members = set(group_members(n))
    k = len(members)
    measure = spec["measure"]
    mode = spec.get("mode", "out")
    normalized = spec.get("normalized", True)
    weighted = bool(spec.get("weighted", False))
    if measure == "degree":
        hit = set()
        for a, b, _ in arcs(fx, mode):
            if a in members:
                hit.add(b)
        count = len(hit - members)
        if not normalized:
            return [float(count)]
        return [UNDEFINED if n == k else count / (n - k)]
    if measure == "closeness":
        d = distance_matrix(n, arcs(fx, mode, weighted), frozenset(), weighted)
        total = 0
        for v in range(n):
            if v in members:
                continue
            ds = [d[c][v] for c in members if d[c][v] is not INF]
            if ds:
                total += min(ds)
        if total == 0:
            return [0.0]
        return [float(Fraction(n - k) / total if normalized else Fraction(1) / total)]
    # Betweenness: path counts by dynamic programming over the distance
    # matrix, all paths and paths avoiding the group. Parallel edges are
    # distinct paths.
    arc_list = arcs(fx, "out", weighted)
    d = distance_matrix(n, arc_list, frozenset(), weighted)
    incoming = [[] for _ in range(n)]
    for a, b, w in arc_list:
        incoming[b].append((a, w))
    total = Fraction(0)
    others = [v for v in range(n) if v not in members]
    for s in others:
        reach = sorted((v for v in range(n) if d[s][v] is not INF), key=lambda v: d[s][v])
        sigma = [0] * n
        avoid = [0] * n
        sigma[s] = avoid[s] = 1
        for t in reach:
            if t == s:
                continue
            for a, w in incoming[t]:
                if d[s][a] is not INF and d[s][a] + w == d[s][t]:
                    sigma[t] += sigma[a]
                    if a not in members:
                        avoid[t] += avoid[a]
        for t in others:
            if t != s and d[s][t] is not INF:
                total += 1 - Fraction(avoid[t], sigma[t])
    if not fx.directed:
        total /= 2
    if not normalized:
        return [float(total)]
    nk = n - k
    if nk < 2:
        return [UNDEFINED]
    pairs = nk * (nk - 1) if fx.directed else Fraction(nk * (nk - 1), 2)
    return [float(total / pairs)]


# ---- Flow betweenness ----------------------------------------------------------

def capacities(fx, use_weights):
    """Residual capacity matrix as a dict of dicts. Parallel edges add;
    undirected edges give capacity both ways; loops carry nothing."""
    cap = [dict() for _ in range(fx.n)]
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        c = exact(fx.weights[i]) if use_weights else 1
        cap[a][b] = cap[a].get(b, 0) + c
        cap[b].setdefault(a, 0)
        if not fx.directed:
            cap[b][a] = cap[b].get(a, 0) + c
    return cap


def max_flow(cap, s, t, removed=None):
    """Edmonds-Karp: shortest augmenting paths found by breadth-first
    search, in exact arithmetic."""
    res = [dict(row) for row in cap]
    flow = 0
    while True:
        parent = {s: None}
        q = deque([s])
        while q and t not in parent:
            v = q.popleft()
            for w, c in res[v].items():
                if c > 0 and w not in parent and w != removed:
                    parent[w] = v
                    q.append(w)
        if t not in parent:
            return flow
        path = []
        w = t
        while parent[w] is not None:
            path.append((parent[w], w))
            w = parent[w]
        b = min(res[u][v] for u, v in path)
        for u, v in path:
            res[u][v] -= b
            res[v][u] = res[v].get(u, 0) + b
        flow += b


def flow_table(fx, use_weights):
    n = fx.n
    cap = capacities(fx, use_weights)
    raw = [Fraction(0)] * n
    frac = [Fraction(0)] * n
    involved = [Fraction(0)] * n
    total = Fraction(0)
    pairs = (itertools.permutations(range(n), 2) if fx.directed
             else itertools.combinations(range(n), 2))
    for i, j in pairs:
        f = max_flow(cap, i, j)
        total += f
        involved[i] += f
        involved[j] += f
        for v in range(n):
            if v in (i, j):
                continue
            lost = f - max_flow(cap, i, j, removed=v)
            raw[v] += lost
            if f > 0:
                frac[v] += Fraction(lost) / f
    norm = []
    for v in range(n):
        denom = total - involved[v]
        norm.append(float(raw[v] / denom) if denom > 0 else UNDEFINED)
    return {"raw": [float(x) for x in raw], "fractional": [float(x) for x in frac],
            "normalized": norm}


def main() -> int:
    spec_all = manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        flows = {}
        for variant, spec in spec_all.items():
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            if name in spec["oracles"][ORACLE].get("undefined_on", []):
                continue
            e = spec["engine"]
            fn = e["fn"]
            if fn == "kpp_single":
                values = kpp_single(fx, e)
            elif fn == "kpp_optimum":
                values = kpp_optimum(fx, e)
            elif fn == "group":
                values = group(fx, e)
            elif fn == "flow_betweenness":
                key = bool(e.get("capacity", False))
                if key not in flows:
                    flows[key] = flow_table(fx, key)
                values = flows[key][e.get("scale", "raw")]
            else:
                raise ValueError(f"no definition for {fn}")
            write_values(ORACLE, name, variant, values)
            count += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    print(f"definition (resilience): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
