#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008).

    python3 validation/oracles/networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (GOLDEN, applies, fixture_names, load_fixture,  # noqa: E402
                    load_manifest, write_values)

ORACLE = "networkx"
TRIADS = ["003", "012", "102", "021D", "021U", "021C", "111D", "111U",
          "030T", "030C", "201", "120D", "120U", "120C", "210", "300"]


def multigraph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def compute(variant: str, fx) -> list[float] | None:
    g = multigraph(fx)
    nodes = range(fx.n)
    if variant == "degree_all":
        d = dict(g.degree())
        return [d[v] for v in nodes]
    if variant == "degree_out":
        d = dict(g.out_degree())
        return [d[v] for v in nodes]
    if variant == "degree_in":
        d = dict(g.in_degree())
        return [d[v] for v in nodes]
    if variant == "strength_all":
        d = dict(g.degree(weight="weight"))
        return [d[v] for v in nodes]
    if variant in ("betweenness", "betweenness_normalized", "betweenness_weighted"):
        b = nx.betweenness_centrality(
            g,
            normalized=variant == "betweenness_normalized",
            weight="weight" if variant == "betweenness_weighted" else None,
        )
        return [b[v] for v in nodes]
    if variant == "triad_census":
        if not fx.directed:
            return None
        c = nx.triadic_census(nx.DiGraph(g))
        return [c[t] for t in TRIADS]
    return None


def main() -> None:
    manifest = load_manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            values = compute(variant, fx)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__}: wrote {count} golden files to {GOLDEN / ORACLE}")
    (GOLDEN / ORACLE / "VERSION").write_text(f"networkx {nx.__version__}\n")


if __name__ == "__main__":
    main()
