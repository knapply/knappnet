#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for cycle
bases, minimum cycle bases, simple and chordless cycles, girth, and
find_cycle.

Simple cycles, chordless cycles, and find_cycle run on NetworkX's
multigraph of the fixture (nodes 0..n - 1, then the edges in file order,
keyed by their position); cycle bases and girth, which NetworkX refuses
for multigraphs, run on its simple graph built the same way. Encodings
follow validation/bindings/cycles.py.

    python3 validation/oracles/cycles_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture,  # noqa: E402
                    write_values, write_version)

ORACLE = "networkx"


def canonical(cycle, directed):
    c = list(cycle)
    i = c.index(min(c))
    c = c[i:] + c[:i]
    if not directed and len(c) > 2 and c[1] > c[-1]:
        c = [c[0]] + c[1:][::-1]
    return tuple(float(v) for v in c)


def multigraph(fx):
    g = (nx.MultiDiGraph if fx.directed else nx.MultiGraph)()
    g.add_nodes_from(range(fx.n))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, key=e)
    return g


def simple(fx, weighted=False):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if weighted:
            w = fx.weights[e]
            if g.has_edge(a, b):
                w = min(w, g[a][b]["weight"])
            g.add_edge(a, b, weight=w)
        else:
            g.add_edge(a, b)
    return g


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    if fn == "cycle_basis":
        return [tuple(float(v) for v in c) for c in nx.cycle_basis(simple(fx))]
    if fn in ("basis_check", "min_basis_check"):
        g = simple(fx)
        return [float(g.number_of_edges() - g.number_of_nodes()
                      + nx.number_connected_components(g)), 1.0]
    if fn == "min_basis_weights":
        weighted = e.get("weighted", False) and fx.weighted
        g = simple(fx, weighted)
        out = []
        for c in nx.minimum_cycle_basis(g, weight="weight" if weighted else None):
            if len(c) == 1:
                pairs = [(c[0], c[0])]
            else:
                pairs = list(zip(c, c[1:] + c[:1]))
            out.append(sum(g[a][b]["weight"] if weighted else 1.0 for a, b in pairs))
        return out
    if fn == "simple_cycles":
        return [canonical(c, fx.directed)
                for c in nx.simple_cycles(multigraph(fx), length_bound=e.get("length_bound"))]
    if fn == "chordless_cycles":
        return [canonical(c, fx.directed)
                for c in nx.chordless_cycles(multigraph(fx), length_bound=e.get("length_bound"))]
    if fn == "girth":
        gi = nx.girth(simple(fx))
        return [-2.0 if gi == float("inf") else float(gi)]
    if fn == "find_cycle":
        orient = {"original": None, "reverse": "reverse", "ignore": "ignore"}[
            e.get("orientation", "original")]
        try:
            c = nx.find_cycle(multigraph(fx), orientation=orient)
        except nx.NetworkXNoCycle:
            return []
        return [(float(x[2]), 1.0 if len(x) > 3 and x[3] == "reverse" else 0.0) for x in c]
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("cycles")
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
    print(f"networkx (cycles): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
