#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for group
degree, closeness, and betweenness (Everett and Borgatti 1999).

NetworkX measures group closeness by distances *to* the group in directed
graphs, so the out-mode variants run on the reversed graph. Group
betweenness runs on the simple graph (NetworkX counts parallel edges as one
path). Undefined values are written as -1, as the engine runner does.

    python3 validation/oracles/resilience_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values)  # noqa: E402

ORACLE = "networkx"


def group_members(n):
    return list(range(min(3, max(1, n - 2))))


def path_graph(fx):
    """A simple graph keeping the shortest of any parallel edges and no
    self-loops."""
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        w = fx.weights[i] if fx.weighted else 1.0
        if not g.has_edge(a, b) or w < g[a][b]["weight"]:
            g.add_edge(a, b, weight=w)
    return g


def compute(spec, fx):
    e = spec["engine"]
    g = path_graph(fx)
    c = group_members(fx.n)
    measure, mode = e["measure"], e.get("mode", "out")
    weight = "weight" if e.get("weighted") else None
    if measure == "degree":
        if fx.directed:
            fn = {"out": nx.group_out_degree_centrality,
                  "in": nx.group_in_degree_centrality}[mode]
            return [fn(g, c)]
        return [nx.group_degree_centrality(g, c)]
    if measure == "closeness":
        h = g.reverse(copy=True) if fx.directed and mode == "out" else g
        return [nx.group_closeness_centrality(h, c, weight=weight)]
    return [nx.group_betweenness_centrality(g, c, normalized=e.get("normalized", True),
                                            weight=weight)]


def main() -> int:
    manifest = load_area("resilience")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", []):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    print(f"networkx (resilience): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
