#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for node
and edge connectivity, minimum cuts, disjoint paths, k-components,
k-edge-components, local bridges, condensation, attracting components,
and semiconnectedness.

Every function runs on NetworkX's simple graph of the fixture (nodes
0..n - 1 added first, then the edges in file order, so parallel edges
collapse and neighbours are stored in input order); self-loops are dropped
where NetworkX would otherwise count them (the engine ignores them
throughout). Encodings follow validation/bindings/connectivity.py.

    python3 validation/oracles/connectivity_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx
from networkx.algorithms import connectivity as nxc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture,  # noqa: E402
                    write_values, write_version)

ORACLE = "networkx"
MISSING = -1.0


def simple(fx, loops=False, order="input"):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    edges = list(zip(fx.src, fx.tgt))
    if order == "sorted":
        if fx.directed:
            edges = sorted(edges)
        else:
            # Sorted rows: add each node's neighbours in increasing order.
            edges = sorted({(min(a, b), max(a, b)) for a, b in edges}
                           | {(max(a, b), min(a, b)) for a, b in edges})
    g.add_edges_from((a, b) for a, b in edges if loops or a != b)
    return g


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    if fn == "local_connectivity":
        if n < 2:
            return [MISSING]
        g = simple(fx)
        f = nxc.local_node_connectivity if e["kind"] == "node" else nxc.local_edge_connectivity
        return [f(g, 0, n - 1)]
    if fn == "connectivity":
        g = simple(fx)
        f = nx.node_connectivity if e["kind"] == "node" else nx.edge_connectivity
        return [f(g)]
    if fn == "all_pairs_node":
        d = nx.all_pairs_node_connectivity(simple(fx))
        return [0.0 if u == v else d[u][v] for u in range(n) for v in range(n)]
    if fn == "average_node":
        return [nx.average_node_connectivity(simple(fx))]
    if fn == "st_edge_cut":
        if n < 2:
            return []
        return sorted(nxc.minimum_st_edge_cut(simple(fx), 0, n - 1))
    if fn == "st_node_cut":
        if n < 2:
            return [MISSING]
        return sorted(float(v) for v in nxc.minimum_st_node_cut(simple(fx), 0, n - 1))
    if fn == "st_node_cut_check":
        g = simple(fx)
        if n < 2 or g.has_edge(0, n - 1):
            return [MISSING, MISSING]
        return [len(nxc.minimum_st_node_cut(g, 0, n - 1)), 1.0]
    if fn == "node_cut_check":
        return [len(nx.minimum_node_cut(simple(fx))), 1.0]
    if fn == "edge_cut_check":
        return [len(nx.minimum_edge_cut(simple(fx))), 1.0]
    if fn == "disjoint_paths_check":
        if n < 2:
            return [MISSING, MISSING]
        f = nx.node_disjoint_paths if e["kind"] == "node" else nx.edge_disjoint_paths
        try:
            count = len(list(f(simple(fx), 0, n - 1)))
        except nx.NetworkXNoPath:
            count = 0
        return [count, 1.0]
    if fn == "k_components":
        out = []
        for k, comps in nx.k_components(simple(fx)).items():
            for c in comps:
                out.append(tuple(float(k * 1000 + v) for v in sorted(c)))
        return out
    if fn == "k_edge_components":
        label = [0] * n
        for i, comp in enumerate(nx.k_edge_components(simple(fx), e["k"])):
            for v in comp:
                label[v] = i
        return label
    if fn == "local_bridges":
        g = simple(fx)
        span = {}
        for u, v, s in nx.local_bridges(g, with_span=True):
            span[(u, v)] = span[(v, u)] = -2.0 if s == float("inf") else float(s)
        return [span.get((a, b), -1.0) if a != b else -1.0 for a, b in zip(fx.src, fx.tgt)]
    if fn == "condensation":
        mapping = nx.condensation(simple(fx, loops=True, order=e["order"])).graph["mapping"]
        return [mapping[v] for v in range(n)]
    if fn == "attracting":
        g = simple(fx, loops=True, order=e["order"])
        mapping = nx.condensation(g).graph["mapping"]
        label = [MISSING] * n
        for comp in nx.attracting_components(g):
            for v in comp:
                label[v] = mapping[v]
        return label
    if fn == "semiconnected":
        return [1.0 if nx.is_semiconnected(simple(fx)) else 0.0]
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("connectivity")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if (ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", [])
                    or tags & set(ospec.get("undefined_when", []))):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (connectivity): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
