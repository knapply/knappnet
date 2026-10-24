#!/usr/bin/env python3
"""Definition oracle for cycle bases, minimum cycle bases, simple and
chordless cycles, and girth.

Independent of the engine: the size of a cycle basis is m - n + c
(Kirchhoff); the weights of a minimum cycle basis come from Horton's
(1987) method (every cycle made of an edge and the two shortest paths from
some node to its ends, sorted by weight, kept greedily while independent
modulo 2), with every edge of the multigraph counted; simple cycles are
enumerated by depth-first search over paths whose first node is the
cycle's lowest (small fixtures only), and chordless cycles are those among
them whose nodes induce exactly as many edges as the cycle is long; the
girth is the least, over all edges u-v, of one plus the distance from v
back to u without that edge. Encodings follow
validation/bindings/cycles.py.

    python3 validation/oracles/cycles_definition_oracle.py
"""

from __future__ import annotations

import heapq
import sys
from fractions import Fraction
from collections import deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
ENUM_MAX_NODES = 15


def components(fx):
    parent = list(range(fx.n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(fx.src, fx.tgt):
        parent[find(a)] = find(b)
    return len({find(v) for v in range(fx.n)})


def horton_weights(fx, weighted):
    n, m = fx.n, len(fx.src)
    w = list(fx.weights) if weighted and fx.weighted else [1.0] * m
    # Exact perturbed weights w * K + 2^e make every path weight distinct
    # (so shortest paths are unique, as Horton's theorem assumes) while
    # keeping every strict comparison of the original totals: the
    # perturbations of a whole basis sum to less than K.
    K = (m + 1) * (1 << m) + 1
    wp = [Fraction(x) * K + (1 << e) for e, x in enumerate(w)]
    adj = [[] for _ in range(n)]
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        adj[a].append((b, e))
        if a != b:
            adj[b].append((a, e))
    candidates = []
    for v in range(n):
        # Shortest path tree from v (Dijkstra, ties by node index).
        dist = {v: Fraction(0)}
        via = {v: None}
        heap = [(Fraction(0), v)]
        done = set()
        while heap:
            d, x = heapq.heappop(heap)
            if x in done:
                continue
            done.add(x)
            for y, e in sorted(adj[x]):
                nd = d + wp[e]
                if y not in dist or nd < dist[y]:
                    dist[y] = nd
                    via[y] = (x, e)
                    heapq.heappush(heap, (nd, y))

        def path_edges(x):
            out = []
            while via[x] is not None:
                p, e = via[x]
                out.append(e)
                x = p
            return out
        for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            if a not in dist or b not in dist:
                continue
            pa, pb = path_edges(a), path_edges(b)
            if e in pa or e in pb:
                continue
            bits = 1 << e
            for x in pa + pb:
                bits ^= 1 << x
            key = sum(wp[x] for x in range(m) if bits >> x & 1)
            weight = sum(w[x] for x in range(m) if bits >> x & 1)
            candidates.append((key, weight, bits))
    candidates.sort()
    need = m - n + components(fx)
    basis, chosen = {}, []
    for _, weight, bits in candidates:
        if len(chosen) == need:
            break
        v = bits
        while v:
            top = v.bit_length() - 1
            if top not in basis:
                basis[top] = v
                chosen.append(weight)
                break
            v ^= basis[top]
    return chosen


def all_cycles(fx, bound=None):
    """Node-sequence cycles of the multigraph, as NetworkX lists them."""
    n = fx.n
    succ = [set() for _ in range(n)]
    loops = set()
    mult = {}
    for a, b in zip(fx.src, fx.tgt):
        if a == b:
            loops.add(a)
            continue
        succ[a].add(b)
        if not fx.directed:
            succ[b].add(a)
        key = (a, b) if fx.directed else (min(a, b), max(a, b))
        mult[key] = mult.get(key, 0) + 1
    out = [(v,) for v in sorted(loops)]
    if bound == 1:
        return out
    if not fx.directed:
        out += [k for k, c in mult.items() if c > 1]
    for s in range(n):
        stack = [(s, [s])]
        while stack:
            v, path = stack.pop()
            for w in succ[v]:
                if w == s:
                    if fx.directed or len(path) >= 3:
                        if fx.directed or path[1] < path[-1]:
                            out.append(tuple(path))
                elif w > s and w not in path and (bound is None or len(path) < bound):
                    stack.append((w, path + [w]))
    return out


def chordless(fx, cycle):
    nodes = set(cycle)
    induced = sum(1 for a, b in zip(fx.src, fx.tgt) if a in nodes and b in nodes)
    return induced == len(cycle)


def canonical(cycle, directed):
    c = list(cycle)
    i = c.index(min(c))
    c = c[i:] + c[:i]
    if not directed and len(c) > 2 and c[1] > c[-1]:
        c = [c[0]] + c[1:][::-1]
    return tuple(float(v) for v in c)


def girth(fx):
    best = float("inf")
    adj = [[] for _ in range(fx.n)]
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        adj[a].append((b, e))
        if not fx.directed and a != b:
            adj[b].append((a, e))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            return 1.0
        # Distance from b to a without edge e.
        dist = {b: 0}
        q = deque([b])
        while q:
            x = q.popleft()
            for y, f in adj[x]:
                if f != e and y not in dist:
                    dist[y] = dist[x] + 1
                    q.append(y)
        if a in dist:
            best = min(best, dist[a] + 1)
    return -2.0 if best == float("inf") else float(best)


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    if fn in ("basis_check", "min_basis_check"):
        return [float(len(fx.src) - fx.n + components(fx)), 1.0]
    if fn == "min_basis_weights":
        return horton_weights(fx, e.get("weighted", False))
    if fn == "girth":
        return [girth(fx)]
    if fx.n > ENUM_MAX_NODES:
        return None
    bound = e.get("length_bound")
    cycles = all_cycles(fx, bound)
    if fn == "chordless_cycles":
        cycles = [c for c in cycles if chordless(fx, c)]
    return [canonical(c, fx.directed) for c in cycles]


def main() -> int:
    manifest = load_area("cycles")
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
    print(f"definition (cycles): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
