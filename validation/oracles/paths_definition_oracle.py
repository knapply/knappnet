#!/usr/bin/env python3
"""Golden values for shortest paths computed directly from the definitions.

This oracle shares no code or algorithmic strategy with the engine:
- distances come from rounds of relaxing every arc (the Bellman-Ford
  recurrence) in exact rational arithmetic, from every source separately;
  a change in round n means a negative cycle reachable from the source
  (the engine uses breadth-first search, Dijkstra's algorithm with a heap,
  a FIFO label-correcting queue, Goldberg-Radzik scans, Johnson
  reweighting, or Floyd-Warshall, in floating point);
- predecessors are the arcs u -> v with d(u) + len(u, v) = d(v) exactly;
- numbers of shortest paths come from the recurrence
  sigma(v) = sum of sigma(u) over predecessors u, and lists of shortest
  paths from following predecessor arcs;
- k shortest simple paths come from enumerating every simple path
  (networks with at most 10 nodes).

The edge lengths repeat validation/bindings/paths.py. Markers:
unreachable -1e9, negative cycle -2e9.

    python3 validation/oracles/paths_definition_oracle.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 80
MAX_NODES_ENUMERATION = 10
INFINITE = -1e9
NEGATIVE_CYCLE = -2e9


def edge_lengths(fx, kind):
    m = len(fx.src)
    w = fx.weights if fx.weights is not None else [1.0] * m
    if kind == "unit":
        return [Fraction(1)] * m
    if kind == "weights":
        return [Fraction(x) for x in w]
    if kind == "strength":
        return [1 / Fraction(x) for x in w]
    if kind == "potential":
        return [Fraction(w[e]) + fx.src[e] % 5 - fx.tgt[e] % 5 for e in range(m)]
    if kind == "shift":
        return [Fraction(w[e]) - (Fraction(3, 2) if e % 3 == 0 else 0) for e in range(m)]
    raise ValueError(kind)


def arcs(fx, e):
    """(u, v, length) for every arc the paths may follow."""
    lens = edge_lengths(fx, e.get("lengths", "unit"))
    mode = e.get("mode", "out")
    out = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if not fx.directed or mode == "all":
            out.append((a, b, lens[i]))
            if a != b:
                out.append((b, a, lens[i]))
        elif mode == "out":
            out.append((a, b, lens[i]))
        else:
            out.append((b, a, lens[i]))
    return out


def relax(n, arc_list, sources):
    """Exact distances from the sources (all at 0), or None when a
    negative cycle is reachable."""
    d = [None] * n
    for s in sources:
        d[s] = Fraction(0)
    for _ in range(n + 1):
        changed = False
        for u, v, length in arc_list:
            if d[u] is None:
                continue
            alt = d[u] + length
            if d[v] is None or alt < d[v]:
                d[v] = alt
                changed = True
        if not changed:
            return d
    return None


def mark(x):
    return INFINITE if x is None else float(x)


def predecessors(n, arc_list, d):
    preds = [set() for _ in range(n)]
    for u, v, length in arc_list:
        if u != v and d[u] is not None and d[v] is not None and d[u] + length == d[v]:
            preds[v].add(u)
    return preds


def simple_paths(n, arc_list, s, t):
    best = {}
    for u, v, length in arc_list:
        if u != v and ((u, v) not in best or length < best[(u, v)]):
            best[(u, v)] = length
    nbrs = [[] for _ in range(n)]
    for (u, v) in sorted(best):
        nbrs[u].append(v)
    out = []
    path = [s]

    def go():
        v = path[-1]
        if v == t:
            out.append(sum((best[(a, b)] for a, b in zip(path, path[1:])), Fraction(0)))
            return
        for w in nbrs[v]:
            if w not in path:
                path.append(w)
                go()
                path.pop()

    go()
    return sorted(out)


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    a = arcs(fx, e)
    cutoff = e.get("cutoff")

    def cut(x):
        return None if x is None or (cutoff is not None and x > Fraction(cutoff)) else x

    if fn in ("paths_sssp", "paths_all_pairs", "paths_bidirectional", "paths_astar"):
        return [tuple(mark(cut(x)) for x in relax(n, a, [s])) for s in range(n)]
    if fn == "paths_multi_source":
        d0, d1 = relax(n, a, [0]), relax(n, a, [n - 1])
        return [mark(min((x for x in (p, q) if x is not None), default=None))
                for p, q in zip(d0, d1)]
    if fn == "paths_predecessors":
        d = relax(n, a, [0])
        preds = predecessors(n, a, d)
        return [tuple(v * (n + 1) + u for u in sorted(preds[v])) or (v * (n + 1) + n,)
                for v in range(n)]
    if fn == "paths_negative":
        d = relax(n, a, [0])
        return [NEGATIVE_CYCLE] if d is None else [mark(x) for x in d]
    if fn == "paths_negative_cycle":
        return [1.0 if relax(n, a, list(range(n))) is None else 0.0]
    if fn in ("paths_johnson", "paths_floyd_warshall"):
        if relax(n, a, list(range(n))) is None:
            return [NEGATIVE_CYCLE]
        return [tuple(mark(x) for x in relax(n, a, [s])) for s in range(n)]
    if fn == "paths_count":
        d = relax(n, a, [0])
        preds = predecessors(n, a, d)
        sigma = [0] * n
        sigma[0] = 1
        for v in sorted((v for v in range(n) if d[v] is not None and v != 0),
                        key=lambda v: d[v]):
            sigma[v] = sum(sigma[u] for u in preds[v])
        return [float(x) for x in sigma]
    if fn == "paths_all_shortest":
        d = relax(n, a, [0])
        if d[n - 1] is None:
            return []
        preds = predecessors(n, a, d)
        paths = []

        def back(path):
            v = path[-1]
            if v == 0:
                paths.append(list(reversed(path)))
                return
            for u in sorted(preds[v]):
                if u not in path:
                    back(path + [u])

        back([n - 1])
        return [tuple(i * n + v for i, v in enumerate(p)) for p in paths]
    if fn == "paths_k_shortest":
        k = int(e.get("k", 10))
        return [float(x) for x in simple_paths(n, a, 0, n - 1)[:k]]
    raise ValueError(f"no definition oracle for {fn}")


def main() -> int:
    manifest = load_area("paths")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, spec in manifest.items():
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            if spec["engine"]["fn"] == "paths_k_shortest" and fx.n > MAX_NODES_ENUMERATION:
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (paths): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
