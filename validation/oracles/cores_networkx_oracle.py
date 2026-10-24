#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
k-cores, k-shells, k-crusts, k-coronas, k-trusses, and onion layers (see
../bindings/cores.py for the encodings).

    python3 validation/oracles/cores_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"


def simple(fx, directed):
    g = nx.DiGraph() if directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def truss_numbers(fx):
    g = simple(fx, False)
    number = {}
    k = 2
    while True:
        h = nx.k_truss(g, k)
        if h.number_of_edges() == 0:
            break
        for u, v in h.edges:
            number[frozenset((u, v))] = k
        k += 1
    return [0 if a == b else number.get(frozenset((a, b)), 2) for a, b in zip(fx.src, fx.tgt)]


def compute(fx, e):
    fn = e["fn"]
    if fn == "core_subset":
        g = simple(fx, fx.directed)
        k = e.get("k")
        f = {"core": nx.k_core, "shell": nx.k_shell, "crust": nx.k_crust,
             "corona": nx.k_corona}[e["which"]]
        try:
            h = f(g, k) if k is not None else f(g)
        except ValueError:
            return []  # max() of an empty graph's core numbers
        return sorted(h.nodes)
    if fn == "truss":
        return truss_numbers(fx)
    if fn == "k_truss":
        h = nx.k_truss(simple(fx, False), int(e["k"]))
        if e["part"] == "nodes":
            return sorted(h.nodes)
        return [i for i, (a, b) in enumerate(zip(fx.src, fx.tgt))
                if a != b and h.has_edge(a, b)]
    if fn == "onion":
        layers = nx.onion_layers(simple(fx, False))
        return [layers[v] for v in range(fx.n)]
    return None


def main() -> None:
    spec = load_area("cores")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (cores): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
