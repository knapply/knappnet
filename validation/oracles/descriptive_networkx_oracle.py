#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
reciprocity, transitivity and clustering, distances, mean distance,
diameter, eccentricity, reachability, cut points, and bridges.

NetworkX defines several of these only on (strongly) connected graphs.
Where it does, the oracle checks that its own function agrees with the
value computed from NetworkX's shortest-path lengths, which are what it
writes, so disconnected fixtures are covered too. NaN is written as -999
and infinity as -1000, as the engine runner does.

    python3 validation/oracles/descriptive_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import networkx as nx
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (GOLDEN, MANIFEST_DIR, applies, fixture_names,  # noqa: E402
                    load_fixture, write_values)

ORACLE = "networkx"
UNDEFINED = -999.0
INFINITE = -1000.0


def multigraph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
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


def simple_undirected(fx):
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def weighted_undirected(fx):
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if g.has_edge(a, b):
            raise ValueError("weighted clustering fixtures must not have parallel edges")
        g.add_edge(a, b, weight=fx.weights[i])
    return g


def lengths(g, weighted):
    if weighted:
        return dict(nx.all_pairs_dijkstra_path_length(g, weight="weight"))
    return dict(nx.all_pairs_shortest_path_length(g))


def connected(g):
    if g.number_of_nodes() == 0:
        return False
    return nx.is_strongly_connected(g) if g.is_directed() else nx.is_connected(g)


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    nodes = range(fx.n)
    if fn == "reciprocity":
        if e["mode"] != "default":
            return None
        g = nx.DiGraph()
        g.add_nodes_from(nodes)
        g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
        if g.number_of_edges() == 0:
            return None
        return [nx.overall_reciprocity(g)]
    if fn == "transitivity":
        return [nx.transitivity(simple_undirected(fx))]
    if fn in ("clustering", "average_clustering"):
        if e.get("undefined") != "zero":
            return None
        weighting = e.get("weighting", "none")
        if weighting == "barrat":
            return None
        if weighting == "onnela":
            g, weight = weighted_undirected(fx), "weight"
        else:
            g, weight = simple_undirected(fx), None
        if fn == "average_clustering":
            return [nx.average_clustering(g, weight=weight)]
        c = nx.clustering(g, weight=weight)
        return [float(c[v]) for v in nodes]
    if fn in ("distances", "mean_distance", "diameter", "eccentricity"):
        if e.get("weights_as") == "strength" or e.get("unconnected", True) is False:
            return None
        g = path_graph(fx)
        mode = e.get("mode", "out")
        if fx.directed and mode == "all":
            g = g.to_undirected()
        elif fx.directed and mode == "in":
            return None
        weighted = bool(e.get("weighted"))
        d = lengths(g, weighted)
        if fn == "distances":
            return [tuple(float(d[i][j]) if j in d[i] else INFINITE for j in nodes)
                    for i in nodes]
        finite = [d[i][j] for i in nodes for j in d[i] if j != i]
        if fn == "mean_distance":
            if not finite:
                return None
            mean = sum(finite) / len(finite)
            if connected(g) and fx.n > 1:
                check = nx.average_shortest_path_length(
                    g, weight="weight" if weighted else None)
                assert math.isclose(check, mean, rel_tol=1e-12), (check, mean)
            return [mean]
        if fn == "diameter":
            longest = float(max(finite, default=0))
            if connected(g) and not weighted:
                assert nx.diameter(g) == longest
            return [longest]
        ecc = [float(max(d[v].values())) for v in nodes]
        if connected(g):
            check = nx.eccentricity(g, weight="weight" if weighted else None)
            assert all(math.isclose(check[v], ecc[v], rel_tol=1e-12) for v in nodes)
        return ecc
    if fn == "reachability":
        g = multigraph(fx)
        mode = e["mode"]
        if not fx.directed:
            return [float(len(nx.node_connected_component(g, v)) - 1) for v in nodes]
        if mode == "out":
            return [float(len(nx.descendants(g, v))) for v in nodes]
        if mode == "in":
            return [float(len(nx.ancestors(g, v))) for v in nodes]
        u = g.to_undirected()
        return [float(len(nx.node_connected_component(u, v)) - 1) for v in nodes]
    if fn in ("cut_points", "bridges"):
        g = nx.MultiGraph()
        g.add_nodes_from(nodes)
        g.add_edges_from(zip(fx.src, fx.tgt))
        if fn == "cut_points":
            return sorted(float(v) for v in nx.articulation_points(g))
        return [(float(min(a, b)), float(max(a, b))) for a, b in nx.bridges(g)]
    return None


def main() -> None:
    spec = yaml.safe_load((MANIFEST_DIR / "descriptive.yaml").read_text())
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(v, fx)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (descriptive): wrote {count} golden files")
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(f"networkx {nx.__version__}\n", newline="")


if __name__ == "__main__":
    main()
