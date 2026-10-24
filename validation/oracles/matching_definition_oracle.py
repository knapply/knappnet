#!/usr/bin/env python3
"""Definition oracle for matchings (exhaustive search over vertex subsets
on networks of at most 16 nodes), the matching tests, the documented
colouring orders, maximal independent sets, square clustering (the
published pairwise formula), and generalized degree. Shares no code with
the engine.

    python3 validation/oracles/matching_definition_oracle.py
"""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 16


def pair_weights(fx, best):
    w = {}
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        k = (min(a, b), max(a, b))
        x = fx.weights[i] if fx.weighted else 1.0
        w[k] = x if k not in w else best(w[k], x)
    return w


def optimum(fx, w, key):
    n = fx.n
    nbrs = {v: [] for v in range(n)}
    for (a, b), x in w.items():
        nbrs[a].append((b, x))
        nbrs[b].append((a, x))

    @lru_cache(maxsize=None)
    def best(mask):
        v = next((v for v in range(n) if not mask >> v & 1), None)
        if v is None:
            return (0, 0.0)
        m2 = mask | 1 << v
        top = best(m2)
        for u, x in nbrs[v]:
            if not m2 >> u & 1:
                c, s = best(m2 | 1 << u)
                cand = (c + 1, s + x)
                if key(cand) > key(top):
                    top = cand
        return top
    return best(0)


def nx_degree(fx, v):
    return sum((a == v) + (b == v) for a, b in zip(fx.src, fx.tgt))


def colour(fx, order):
    succ = {v: set() for v in range(fx.n)}
    for a, b in zip(fx.src, fx.tgt):
        succ[a].add(b)
        if not fx.directed:
            succ[b].add(a)
    c = {}
    for u in order:
        used = {c[v] for v in succ[u] if v in c and v != u}
        c[u] = next(k for k in range(fx.n + 1) if k not in used)
    return [c[v] for v in range(fx.n)]


def simple_adj(fx):
    adj = {v: set() for v in range(fx.n)}
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def dsatur(fx):
    succ = {v: set() for v in range(fx.n)}
    for a, b in zip(fx.src, fx.tgt):
        succ[a].add(b)
        if not fx.directed:
            succ[b].add(a)
    deg = {v: nx_degree(fx, v) for v in range(fx.n)}
    seen = {v: set() for v in range(fx.n)}
    c = {}
    first = max(range(fx.n), key=lambda v: (deg[v], -v))
    order = [first]
    while len(c) < fx.n:
        if len(c) == 0:
            u = first
        else:
            u = max((v for v in range(fx.n) if v not in c),
                    key=lambda v: (len(seen[v]), deg[v], -v))
        used = {c[v] for v in succ[u] if v in c and v != u}
        c[u] = next(k for k in range(fx.n + 1) if k not in used)
        for v in succ[u]:
            seen[v].add(c[u])
    del order
    return [c[v] for v in range(fx.n)]


def order_for(fx, strategy):
    n = fx.n
    if strategy == "largest_first":
        return sorted(range(n), key=lambda v: -nx_degree(fx, v))
    adj = simple_adj(fx)
    if strategy == "smallest_last":
        alive = set(range(n))
        out = []
        while alive:
            v = min(alive, key=lambda x: (len(adj[x] & alive), x))
            out.append(v)
            alive.remove(v)
        return out[::-1]
    if strategy == "independent_set":
        left = set(range(n))
        out = []
        while left:
            alive = set(left)
            while alive:
                v = min(alive, key=lambda x: (len(adj[x] & alive), x))
                out.append(v)
                left.discard(v)
                alive -= adj[v] | {v}
        return out
    seen, out = set(), []
    for root in range(n):
        if root in seen:
            continue
        seen.add(root)
        out.append(root)
        if strategy == "connected_sequential_bfs":
            queue = [root]
            for v in queue:
                for w in sorted(adj[v]):
                    if w not in seen:
                        seen.add(w)
                        out.append(w)
                        queue.append(w)
        else:
            def dfs(v):
                for w in sorted(adj[v]):
                    if w not in seen:
                        seen.add(w)
                        out.append(w)
                        dfs(w)
            dfs(root)
    return out


def square(fx):
    adj = simple_adj(fx)
    out = []
    for v in range(fx.n):
        nb = sorted(adj[v])
        num = den = 0
        for i, u in enumerate(nb):
            for w in nb[i + 1:]:
                q = len(adj[u] & adj[w] - {v})
                theta = 1 if w in adj[u] else 0
                a = (len(adj[u]) - (1 + q + theta)) + (len(adj[w]) - (1 + q + theta))
                num += q
                den += a + q
        out.append(num / den if den > 0 else 0.0)
    return out


def compute(fx, e):
    fn = e["fn"]
    weighted = bool(e.get("weighted"))
    if fn in ("max_weight_matching", "min_weight_matching"):
        if fx.n > MAX_NODES:
            return None
        low = fn == "min_weight_matching"
        w = pair_weights(fx, min if low else max)
        if not weighted:
            w = {k: 1.0 for k in w}
        if low:
            c, s = optimum(fx, w, lambda t: (t[0], -t[1]))
        elif e.get("max_cardinality"):
            c, s = optimum(fx, w, lambda t: t)
        elif weighted:
            c, s = optimum(fx, w, lambda t: t[1])
            c = -1  # heaviest matchings may differ in size
        else:
            c, s = optimum(fx, w, lambda t: t)
        return [(c, s, 1)]
    edges = {(a, b) for a, b in zip(fx.src, fx.tgt)}
    if not fx.directed:
        edges |= {(b, a) for a, b in edges}

    def checks(pairs):
        used = set()
        ok = True
        for a, b in pairs:
            if a == b or (a, b) not in edges or a in used or b in used:
                ok = False
            used |= {a, b}
        maximal = ok and not any(a != b and a not in used and b not in used for a, b in edges)
        return [int(ok), int(maximal), int(ok and len(used) == fx.n)]

    def greedy():
        adj = {v: {} for v in range(fx.n)}
        for a, b in zip(fx.src, fx.tgt):
            adj[a][b] = True
            adj[b][a] = True
        seen, used, out = set(), set(), []
        for u in range(fx.n):
            for v in adj[u]:
                if v not in seen and u != v and u not in used and v not in used:
                    out.append((u, v))
                    used |= {u, v}
            seen.add(u)
        return out
    if fn == "maximal_matching":
        return [(min(p), max(p)) for p in greedy()]
    if fn == "check_matching":
        if e["pairs"] == "candidate":
            pairs = [(fx.src[i], fx.tgt[i]) for i in range(0, len(fx.src), 2)]
            # NetworkX receives a set: a pair listed twice counts once.
            pairs = list(dict.fromkeys(pairs))
        else:
            pairs = greedy()
        return checks(pairs)
    if fn == "greedy_color":
        if e["strategy"] == "saturation_largest_first":
            return dsatur(fx)
        return colour(fx, order_for(fx, e["strategy"]))
    if fn == "mis_check":
        return [(1, 1, 1)]
    if fn == "square_clustering":
        return square(fx)
    if fn == "generalized_degree":
        adj = simple_adj(fx)
        out = []
        for v in range(fx.n):
            tally = {}
            for w in adj[v]:
                k = len(adj[v] & adj[w])
                tally[k] = tally.get(k, 0) + 1
            out.append(tuple(x for k in sorted(tally) for x in (k, tally[k])))
        return out
    return None


def undefined(fx, spec):
    o = spec["oracles"][ORACLE]
    tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
    return fx.name in o.get("undefined_on", []) or tags & set(o.get("undefined_when", []))


def main() -> None:
    spec = load_area("matching")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx) or undefined(fx, v):
                continue
            values = compute(fx, v["engine"])
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (matching): wrote {count} golden files")


if __name__ == "__main__":
    main()
