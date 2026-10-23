#!/usr/bin/env python3
"""Golden values computed directly from textbook definitions.

This oracle shares no code or algorithmic strategy with the engine or the
other oracles. It is deliberately slow and simple, so it only runs on small
fixtures. For betweenness it evaluates Freeman's definition (Freeman 1977)
over all ordered pairs, counting shortest paths with exact rational
arithmetic. Parallel edges count as distinct paths, which arbitrates between
the igraph and NetworkX conventions on multigraphs.

    python3 validation/oracles/core_definition_oracle.py
"""

from __future__ import annotations

import heapq
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_fixture,
                    load_manifest, write_values, write_version)  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80


def arcs(fx, weighted: bool, strength: bool = False):
    """Adjacency list of (target, length) with parallel arcs kept."""
    adj = [[] for _ in range(fx.n)]
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        if weighted and strength:
            w = 1 / w
        adj[a].append((b, w))
        if not fx.directed:
            adj[b].append((a, w))
    return adj


def shortest(adj, s):
    """Exact distances and path counts from s (parallel arcs count)."""
    n = len(adj)
    dist = [None] * n
    dist[s] = Fraction(0)
    heap = [(Fraction(0), s)]
    done = [False] * n
    order = []
    while heap:
        d, v = heapq.heappop(heap)
        if done[v] or d != dist[v]:
            continue
        done[v] = True
        order.append(v)
        for w, length in adj[v]:
            alt = d + length
            if dist[w] is None or alt < dist[w]:
                dist[w] = alt
                heapq.heappush(heap, (alt, w))
    sigma = [0] * n
    sigma[s] = 1
    for v in order:
        for w, length in adj[v]:
            if dist[w] is not None and dist[v] + length == dist[w]:
                sigma[w] += sigma[v]
    return dist, sigma


def betweenness(fx, weighted: bool, normalized: bool, strength: bool = False) -> list[float]:
    adj = arcs(fx, weighted, strength)
    n = fx.n
    table = [shortest(adj, s) for s in range(n)]
    bc = [Fraction(0)] * n
    for s in range(n):
        ds, ss = table[s]
        for t in range(n):
            if t == s or ds[t] is None:
                continue
            for v in range(n):
                if v in (s, t) or ds[v] is None:
                    continue
                dv, sv = table[v]
                if dv[t] is not None and ds[v] + dv[t] == ds[t]:
                    bc[v] += Fraction(ss[v] * sv[t], ss[t])
    scale = Fraction(1) if fx.directed else Fraction(1, 2)
    if normalized and n > 2:
        pairs = (n - 1) * (n - 2)
        scale *= Fraction(1, pairs) if fx.directed else Fraction(2, pairs)
    return [float(b * scale) for b in bc]


def degree(fx, mode: str, loops: bool) -> list[float]:
    out = [0] * fx.n
    for a, b in zip(fx.src, fx.tgt):
        if a == b and not loops:
            continue
        if not fx.directed or mode != "in":
            out[a] += 1
        if not fx.directed or mode != "out":
            out[b] += 1
    return [float(x) for x in out]


def triad_census(fx) -> list[float]:
    """Classifies every triple by its (mutual, asymmetric, null) dyad counts
    and, where those are ambiguous, by orientation (Holland and Leinhardt 1976)."""
    arc = set()
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            arc.add((a, b))
            if not fx.directed:
                arc.add((b, a))
    names = ["003", "012", "102", "021D", "021U", "021C", "111D", "111U",
             "030T", "030C", "201", "120D", "120U", "120C", "210", "300"]
    counts = dict.fromkeys(names, 0)
    n = fx.n
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                counts[classify((i, j, k), arc)] += 1
    return [float(counts[x]) for x in names]


def classify(nodes, arc) -> str:
    i, j, k = nodes
    pairs = [(i, j), (i, k), (j, k)]
    m = sum((a, b) in arc and (b, a) in arc for a, b in pairs)
    a_ = sum(((a, b) in arc) != ((b, a) in arc) for a, b in pairs)
    nul = 3 - m - a_
    base = f"{m}{a_}{nul}"
    outdeg = {v: sum((v, w) in arc for w in nodes if w != v) for v in nodes}
    indeg = {v: sum((w, v) in arc for w in nodes if w != v) for v in nodes}
    if base == "021":
        if any(outdeg[v] == 2 for v in nodes):
            return "021D"
        if any(indeg[v] == 2 for v in nodes):
            return "021U"
        return "021C"
    if base == "111":
        # 111D (A<->B<-C): the asymmetric arc points into the mutual dyad.
        # 111U (A<->B->C): it points out of the mutual dyad.
        mutual_nodes = {v for a, b in pairs for v in (a, b)
                        if (a, b) in arc and (b, a) in arc}
        for a, b in pairs:
            if ((a, b) in arc) != ((b, a) in arc):
                src, dst = (a, b) if (a, b) in arc else (b, a)
                return "111D" if dst in mutual_nodes else "111U"
    if base == "030":
        return "030C" if all(outdeg[v] == 1 for v in nodes) else "030T"
    if base == "120":
        mutual_nodes = {v for a, b in pairs for v in (a, b)
                        if (a, b) in arc and (b, a) in arc}
        third = next(v for v in nodes if v not in mutual_nodes)
        if outdeg[third] == 2:
            return "120D"
        if indeg[third] == 2:
            return "120U"
        return "120C"
    return base


def main() -> None:
    manifest = load_manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, spec in manifest.items():
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            e = spec["engine"]
            if e["fn"] == "betweenness":
                values = betweenness(fx, e["weighted"], e["normalized"],
                                     e.get("weights_as") == "strength")
            elif e["fn"] == "degree":
                values = degree(fx, e["mode"], e["loops"])
            elif e["fn"] == "triad_census":
                values = triad_census(fx)
            else:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition: wrote {count} golden files")


if __name__ == "__main__":
    main()
