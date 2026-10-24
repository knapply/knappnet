#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for edge
swaps, reference graphs, the small-world coefficients, and the rich-club
coefficient (see ../bindings/rewiring.py for the contracts and encodings).

    python3 validation/oracles/rewiring_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
import warnings
from collections import Counter
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"


def simple(fx):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from(zip(fx.src, fx.tgt))
    return g


def contract(fx, procedure, before, after):
    simple_ok = True  # a NetworkX Graph or DiGraph has no parallel edges
    if procedure == "directed":
        return (int(dict(before.in_degree) == dict(after.in_degree)),
                int(dict(before.out_degree) == dict(after.out_degree)), int(simple_ok),
                after.number_of_edges())
    same = Counter(dict(before.degree)) == Counter(dict(after.degree)) and \
        dict(before.degree) == dict(after.degree)
    m = after.number_of_edges()
    if procedure == "double":
        return (int(same), int(simple_ok), m)
    if procedure == "connected":
        return (int(same), int(nx.is_connected(after)), m)
    keeps = (not nx.is_connected(before)) or nx.is_connected(after)
    return (int(same), int(keeps), m)


def compute(fx, e):
    g = simple(fx)
    fn = e["fn"]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if fn == "swap_contract":
                p = e["procedure"]
                h = g.copy()
                if p == "double":
                    nx.double_edge_swap(h, nswap=10, max_tries=1000, seed=7)
                elif p == "connected":
                    nx.connected_double_edge_swap(h, nswap=20, seed=7)
                elif p == "directed":
                    nx.directed_edge_swap(h, nswap=5, max_tries=500, seed=7)
                elif p == "random":
                    h = nx.random_reference(g, niter=2, seed=7)
                else:
                    h = nx.lattice_reference(g, niter=2, seed=7)
                return [contract(fx, p, g, h)]
            if fn == "small_world":
                f = nx.sigma if e["coefficient"] == "sigma" else nx.omega
                vals = [f(g, niter=e["niter"], nrand=e["nrand"], seed=s)
                        for s in range(e.get("seeds", 1))]
                r = sum(vals) / len(vals)
                return [-999.0 if math.isnan(r) else r]
            if fn == "rich_club":
                rc = nx.rich_club_coefficient(g, normalized=bool(e.get("normalized")),
                                              Q=e.get("q", 100), n_samples=e.get("n_samples", 1),
                                              seed=7)
                return [rc[k] for k in range(len(rc))]
    except (nx.NetworkXError, nx.NetworkXAlgorithmError, nx.NetworkXUnfeasible, Exception):
        return [-1]
    return None


def undefined(fx, spec):
    o = spec["oracles"][ORACLE]
    tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
    return fx.name in o.get("undefined_on", []) or tags & set(o.get("undefined_when", []))


def main() -> None:
    spec = load_area("rewiring")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx) or undefined(fx, v):
                continue
            values = compute(fx, v["engine"])
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (rewiring): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
