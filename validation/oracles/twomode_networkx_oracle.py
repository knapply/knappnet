#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for ego
networks (`ego_graph`) and two-mode networks (the `bipartite` module).

Ego measures use the simple graph (parallel edges collapsed, self-loops
dropped), as the engine does. Undefined values are written as -1, as the
engine runner does.

    python3 validation/oracles/twomode_networkx_oracle.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import networkx as nx
from networkx.algorithms import bipartite

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (FIXTURES, applies, fixture_names, load_area, load_fixture,
                    write_values, write_version)  # noqa: E402

ORACLE = "networkx"
UNDEFINED = -1.0


def manifest() -> dict:
    return load_area("twomode")


def node_types(fx):
    return json.loads((FIXTURES / fx.name / "meta.json").read_text()).get("type")


def simple_graph(fx):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def weighted_graph(fx):
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def ego_nodes(g, v, order, mode, center):
    """Nodes of NetworkX's ego graph: successors for `out`, predecessors
    (on the reversed graph) for `in`, both directions for `all`."""
    if g.is_directed() and mode == "in":
        h = nx.ego_graph(g.reverse(copy=True), v, radius=order, center=center)
    else:
        h = nx.ego_graph(g, v, radius=order, center=center,
                         undirected=g.is_directed() and mode == "all")
    return set(h.nodes())


def ego_values(fx, e):
    g = simple_graph(fx)
    fn, order, mode = e["fn"], e.get("order", 1), e.get("mode", "all")
    out = []
    for v in range(fx.n):
        if fn == "ego_size":
            nodes = ego_nodes(g, v, order, mode, center=e.get("mindist", 0) == 0)
            out.append(float(len(nodes)))
            continue
        if fn == "ego_members":
            nodes = ego_nodes(g, v, order, mode, center=e.get("mindist", 0) == 0)
            out.append(tuple(float(u) for u in sorted(nodes)))
            continue
        measure = e["measure"]
        if measure == "betweenness":
            h = g.subgraph(ego_nodes(g, v, order, mode, center=True))
            out.append(nx.betweenness_centrality(h, normalized=False)[v])
            continue
        h = g.subgraph(ego_nodes(g, v, order, mode, center=False))
        if measure == "ties":
            out.append(float(h.number_of_edges()))
        elif measure == "density":
            out.append(nx.density(h) if h.number_of_nodes() > 1 else UNDEFINED)
        elif measure == "components":
            count = (nx.number_weakly_connected_components(h) if h.is_directed()
                     else nx.number_connected_components(h))
            out.append(float(count))
        else:
            return None
    return out


def projection_edges(p):
    out = []
    for a, b, w in p.edges(data="weight"):
        a, b = min(a, b), max(a, b)
        out.append((float(a), float(b), float(w)))
    return sorted(out)


def product_weight(g, u, v):
    return sum(g[u][e]["weight"] * g[v][e]["weight"] for e in set(g[u]) & set(g[v]))


def twomode_values(fx, types, e):
    g = weighted_graph(fx)
    first = [v for v in range(fx.n) if types[v] == 0]
    second = [v for v in range(fx.n) if types[v] == 1]
    if e["fn"] == "project":
        if e.get("loops"):
            return None
        nodes = second if e["second"] else first
        weight = e.get("weight", "count")
        if weight == "count":
            return projection_edges(bipartite.weighted_projected_graph(g, nodes))
        if weight == "newman":
            return projection_edges(bipartite.collaboration_weighted_projected_graph(g, nodes))
        return projection_edges(bipartite.generic_weighted_projected_graph(
            g, nodes, weight_function=lambda b, u, v: product_weight(b, u, v)))
    measure = e["measure"]
    if measure == "density":
        return [bipartite.density(g, first)]
    if measure == "degree":
        d = bipartite.degree_centrality(g, first)
    elif measure == "closeness":
        if not e.get("normalized", True):
            return None
        d = bipartite.closeness_centrality(
            g, first, normalized=e.get("unreachable") == "wasserman_faust")
    elif measure == "betweenness":
        if e.get("normalized", True):
            d = bipartite.betweenness_centrality(g, first)
        else:
            d = nx.betweenness_centrality(g, normalized=False)
    else:
        return None
    return [float(d[v]) for v in range(fx.n)]


def compute(fx, e):
    if e["fn"] in ("ego_size", "ego_members", "ego_summary"):
        return ego_values(fx, e)
    types = node_types(fx)
    if types is None:
        return None
    return twomode_values(fx, types, e)


def main() -> None:
    spec = manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx {nx.__version__} (two-mode and ego): wrote {count} golden files")


if __name__ == "__main__":
    main()
