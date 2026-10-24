#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
shortest paths: distances from every source and of every pair, predecessor
lists, multi-source distances, bidirectional and A* search, Bellman-Ford,
Goldberg-Radzik, negative cycles, Johnson, Floyd-Warshall, counts and
lists of shortest paths, and k shortest simple paths.

Networks are built as MultiGraph or MultiDiGraph with each edge's length
in the attribute "len" (NetworkX takes the lightest of parallel edges).
The edge lengths repeat validation/bindings/paths.py: unit, the weights,
strengths 1/w, potential w + (u mod 5) - (v mod 5), and shift w - 1.5 on
every third edge. Markers: unreachable -1e9, negative cycle -2e9.

    python3 validation/oracles/paths_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"
INFINITE = -1e9
NEGATIVE_CYCLE = -2e9


def edge_lengths(fx, kind):
    m = len(fx.src)
    w = fx.weights if fx.weights is not None else [1.0] * m
    if kind == "unit":
        return [1.0] * m
    if kind == "weights":
        return list(w)
    if kind == "strength":
        return [1.0 / x for x in w]
    if kind == "potential":
        return [w[e] + fx.src[e] % 5 - fx.tgt[e] % 5 for e in range(m)]
    if kind == "shift":
        return [w[e] - (1.5 if e % 3 == 0 else 0.0) for e in range(m)]
    raise ValueError(kind)


def build(fx, e):
    kind = e.get("lengths", "unit")
    lens = edge_lengths(fx, kind)
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, len=lens[i])
    mode = e.get("mode", "out")
    if fx.directed and mode == "in":
        g = g.reverse(copy=True)
    elif fx.directed and mode == "all":
        g = g.to_undirected()
    weight = None if kind == "unit" else "len"
    return g, weight


def simple(g):
    """The simple graph keeping the lightest of parallel edges, no loops."""
    h = nx.DiGraph() if g.is_directed() else nx.Graph()
    h.add_nodes_from(g)
    for a, b, d in g.edges(data=True):
        if a == b:
            continue
        if not h.has_edge(a, b) or d["len"] < h[a][b]["len"]:
            h.add_edge(a, b, len=d["len"])
    return h


def row(n, dist):
    return tuple(float(dist[v]) if v in dist else INFINITE for v in range(n))


def pair_length(fn):
    try:
        return float(fn())
    except nx.NetworkXNoPath:
        return INFINITE


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    g, weight = build(fx, e)
    n = fx.n
    w = weight or None
    if fn == "paths_sssp":
        cutoff = e.get("cutoff")
        rows = []
        for s in range(n):
            if w is None:
                d = nx.single_source_shortest_path_length(g, s, cutoff=cutoff)
            else:
                d = nx.single_source_dijkstra_path_length(g, s, cutoff=cutoff, weight=w)
            rows.append(row(n, d))
        return rows
    if fn == "paths_multi_source":
        d = nx.multi_source_dijkstra_path_length(g, {0, n - 1}, weight=w)
        return list(row(n, d))
    if fn == "paths_predecessors":
        if w is None:
            pred = nx.predecessor(g, 0)
        else:
            pred = nx.dijkstra_predecessor_and_distance(g, 0, weight=w)[0]
        return [tuple(v * (n + 1) + u for u in pred.get(v, [])) or (v * (n + 1) + n,)
                for v in range(n)]
    if fn == "paths_all_pairs":
        if w is None:
            d = dict(nx.all_pairs_shortest_path_length(g))
        else:
            d = dict(nx.all_pairs_dijkstra_path_length(g, weight=w))
        return [row(n, d[s]) for s in range(n)]
    if fn == "paths_bidirectional":
        rows = []
        for s in range(n):
            r = []
            for t in range(n):
                if w is None:
                    r.append(pair_length(lambda: len(nx.bidirectional_shortest_path(g, s, t)) - 1))
                else:
                    r.append(pair_length(lambda: nx.bidirectional_dijkstra(g, s, t, weight=w)[0]))
            rows.append(tuple(r))
        return rows
    if fn == "paths_astar":
        factor = float(e.get("factor", 0.5))
        rows = [[INFINITE] * n for _ in range(n)]
        rev = g.reverse(copy=False) if g.is_directed() else g
        for t in range(n):
            if w is None:
                to_t = nx.single_source_shortest_path_length(rev, t)
            else:
                to_t = nx.single_source_dijkstra_path_length(rev, t, weight=w)

            def h(u, _t, to_t=to_t):
                return factor * to_t.get(u, 0.0)

            for s in range(n):
                rows[s][t] = pair_length(
                    lambda: nx.astar_path_length(g, s, t, heuristic=h, weight=w or (lambda *a: 1)))
        return [tuple(r) for r in rows]
    if fn == "paths_negative":
        try:
            if e.get("method", "bellman_ford") == "bellman_ford":
                d = nx.single_source_bellman_ford_path_length(g, 0, weight="len")
            else:
                d = nx.goldberg_radzik(g, 0, weight="len")[1]
        except nx.NetworkXUnbounded:
            return [NEGATIVE_CYCLE]
        return list(row(n, d))
    if fn == "paths_negative_cycle":
        return [1.0 if nx.negative_edge_cycle(g, weight="len") else 0.0]
    if fn == "paths_johnson":
        try:
            paths = nx.johnson(g, weight="len")
        except nx.NetworkXUnbounded:
            return [NEGATIVE_CYCLE]
        h = simple_lengths(g)
        rows = []
        for s in range(n):
            p = paths.get(s, {s: [s]})
            rows.append(tuple(path_sum(h, p[v]) if v in p else INFINITE for v in range(n)))
        return rows
    if fn == "paths_floyd_warshall":
        try:
            d = nx.floyd_warshall(g, weight=w or "len")
        except nx.NetworkXUnbounded:
            return [NEGATIVE_CYCLE]
        return [tuple(INFINITE if math.isinf(d[s][v]) else float(d[s][v]) for v in range(n))
                for s in range(n)]
    if fn == "paths_count":
        out = []
        for t in range(n):
            try:
                out.append(float(sum(1 for _ in nx.all_shortest_paths(g, 0, t, weight=w))))
            except nx.NetworkXNoPath:
                out.append(0.0)
        return out
    if fn == "paths_all_shortest":
        try:
            paths = list(nx.all_shortest_paths(g, 0, n - 1, weight=w))
        except nx.NetworkXNoPath:
            paths = []
        return [tuple(i * n + v for i, v in enumerate(p)) for p in paths]
    if fn == "paths_k_shortest":
        k = int(e.get("k", 10))
        h = simple(g)
        out = []
        try:
            for p in nx.shortest_simple_paths(h, 0, n - 1, weight=w):
                out.append(float(nx.path_weight(h, p, "len")) if w else float(len(p) - 1))
                if len(out) == k:
                    break
        except nx.NetworkXNoPath:
            pass
        return out
    raise ValueError(f"no NetworkX oracle for {fn}")


def simple_lengths(g):
    """Lightest length of each ordered pair, as NetworkX's weight function
    reads a multigraph."""
    best = {}
    for a, b, d in g.edges(data=True):
        for p in ((a, b), (b, a)) if not g.is_directed() else ((a, b),):
            if p not in best or d["len"] < best[p]:
                best[p] = d["len"]
    return best


def path_sum(h, path):
    total = 0.0
    for a, b in zip(path, path[1:]):
        total += h[(a, b)]
    return total


def main() -> int:
    manifest = load_area("paths")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", []):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (paths): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
