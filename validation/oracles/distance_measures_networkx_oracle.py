#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for center,
periphery, barycenter, the Wiener, Schultz, Gutman, and hyper-Wiener
indices, harmonic diameter, global and local efficiency, and closeness
vitality.

Networks are MultiGraph or MultiDiGraph objects with the fixture weights in
"weight". A NetworkX error on a node set means the network is not
connected and is written as an empty set; on a number, as -999 (NaN).
Infinity is -1000 and minus infinity -1001.

    python3 validation/oracles/distance_measures_networkx_oracle.py
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
UNDEFINED = -999.0
INFINITE = -1000.0
NEG_INFINITE = -1001.0


def mark(x):
    x = float(x)
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE if x > 0 else NEG_INFINITE
    return x


def graph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if fx.weighted:
            g.add_edge(a, b, weight=fx.weights[i])
        else:
            g.add_edge(a, b)
    return g


def simple(fx):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if fx.weighted:
            g.add_edge(a, b, weight=fx.weights[i])
        else:
            g.add_edge(a, b)
    return g


def node_set(fn):
    try:
        return sorted(float(v) for v in fn())
    except (nx.NetworkXError, nx.NetworkXNoPath, nx.NetworkXException):
        return []


def number(fn):
    try:
        return mark(fn())
    except (nx.NetworkXError, nx.NetworkXException, ZeroDivisionError):
        return UNDEFINED


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    g = graph(fx)
    w = "weight" if e.get("weighted") else None
    if fn == "dm_center":
        return node_set(lambda: nx.center(g, weight=w))
    if fn == "dm_periphery":
        return node_set(lambda: nx.periphery(g, weight=w))
    if fn == "dm_barycenter":
        return node_set(lambda: nx.barycenter(g, weight=w))
    if fn == "dm_index":
        index = e["index"]
        if index == "wiener":
            return [number(lambda: nx.wiener_index(g, weight=w))]
        if index == "schultz":
            return [number(lambda: nx.schultz_index(simple(fx), weight=w))]
        if index == "gutman":
            return [number(lambda: nx.gutman_index(simple(fx), weight=w))]
        if index == "hyper_wiener_networkx":
            return [number(lambda: nx.hyper_wiener_index(simple(fx), weight=w))]
        if index == "harmonic_diameter":
            return [number(lambda: nx.harmonic_diameter(g, weight=w))]
        if index == "global_efficiency":
            return [number(lambda: nx.global_efficiency(g))]
    if fn == "dm_local_efficiency":
        return [number(lambda: nx.local_efficiency(g))]
    if fn == "dm_closeness_vitality":
        try:
            r = nx.closeness_vitality(g, weight=w)
        except (nx.NetworkXError, nx.NetworkXException):
            return [UNDEFINED] * fx.n
        return [mark(r[v]) for v in range(fx.n)]
    raise ValueError(f"no NetworkX oracle for {fn} {e}")


def main() -> int:
    manifest = load_area("distance_measures")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx):
                continue
            if name in ospec.get("undefined_on", []) or tags & set(ospec.get("undefined_when",
                                                                              [])):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (distance_measures): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
