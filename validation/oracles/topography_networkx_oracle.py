#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
density, components, closeness, harmonic centrality, and eigenvector
centrality.

NetworkX measures closeness and harmonic centrality on paths *to* a node
in directed graphs, so the out-mode variants run on the reversed graph.
Undefined values are written as -1, as the engine runner does.

    python3 validation/oracles/topography_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,
                    write_version)  # noqa: E402

ORACLE = "networkx"


def multigraph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def simple_graph(fx):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def path_graph(fx):
    """A simple graph keeping the shortest of any parallel edges, which is
    what NetworkX's path searches use on a multigraph."""
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if fx.weighted else 1.0
        if not g.has_edge(a, b) or w < g[a][b]["weight"]:
            g.add_edge(a, b, weight=w)
    return g


def oriented(g, mode):
    if not g.is_directed() or mode == "in":
        return g
    if mode == "out":
        return g.reverse(copy=True)
    raise ValueError("NetworkX closeness has no mode 'all'")


def labels(groups, n):
    out = [0.0] * n
    for c in groups:
        first = float(min(c))
        for v in c:
            out[v] = first
    return out


def compute(variant: str, spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    nodes = range(fx.n)
    if fn == "density":
        if e["loops"]:
            return None
        g = simple_graph(fx) if e.get("simple") else multigraph(fx)
        return [nx.density(g)]
    if fn == "components":
        g = simple_graph(fx)
        if not fx.directed:
            groups = nx.connected_components(g)
        elif e["kind"] == "weak":
            groups = nx.weakly_connected_components(g)
        else:
            groups = nx.strongly_connected_components(g)
        return labels(groups, fx.n)
    if fn in ("closeness", "harmonic"):
        mode = e.get("mode", "out")
        if fx.directed and mode == "all":
            return None
        if e.get("weights_as") == "strength":
            return None
        g = oriented(path_graph(fx), mode)
        dist = "weight" if e.get("weighted") else None
        if fn == "harmonic":
            if e.get("normalized"):
                return None
            h = nx.harmonic_centrality(g, distance=dist)
            return [float(h[v]) for v in nodes]
        if e.get("unreachable") != "wasserman_faust":
            return None
        c = nx.closeness_centrality(g, distance=dist, wf_improved=True)
        return [c[v] for v in nodes]
    if fn == "betweenness_simple":
        b = nx.betweenness_centrality(multigraph(fx), normalized=e["normalized"],
                                      weight="weight" if e["weighted"] else None)
        return [b[v] for v in nodes]
    if fn == "eigenvector":
        if fx.directed and e.get("mode", "in") != "in":
            return None
        if e.get("loops") == "twice" and any(a == b for a, b in zip(fx.src, fx.tgt)) \
                and not fx.directed:
            return None
        g = multigraph(fx)
        try:
            x = nx.eigenvector_centrality_numpy(
                g, weight="weight" if e.get("weighted") else None, max_iter=100000, tol=0)
        except (nx.AmbiguousSolution, nx.NetworkXException, TypeError):
            # Not (strongly) connected, or too small for ARPACK.
            return None
        values = [x[v] for v in nodes]
        if e.get("scale", "max") == "max":
            top = max(values)
            values = [v / top for v in values]
        return values
    return None


def main() -> None:
    spec = load_area("topography")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(variant, v, fx)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (topography): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
