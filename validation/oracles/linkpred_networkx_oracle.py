#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for link
prediction scores: NetworkX's own functions on the undirected simple graph
of each fixture, community labels node % 3.

    python3 validation/oracles/linkpred_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,
                    write_version)  # noqa: E402

ORACLE = "networkx"


def simple_graph(fx):
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from(zip(fx.src, fx.tgt))
    for v in g:
        g.nodes[v]["community"] = v % 3
    return g


def scores(g, e):
    if e["pairs"] == "non_edges":
        pairs = sorted(tuple(sorted(p)) for p in nx.non_edges(g))
    else:
        pairs = [(u, v) for u in range(len(g)) for v in range(u + 1, len(g))]
    m = e["method"]
    if m == "common_neighbors":
        return [(u, v, len(nx.common_neighbors(g, u, v))) for u, v in pairs]
    funcs = {
        "jaccard": nx.jaccard_coefficient,
        "adamic_adar": nx.adamic_adar_index,
        "resource_allocation": nx.resource_allocation_index,
        "preferential_attachment": nx.preferential_attachment,
        "cnc": lambda G, eb: nx.common_neighbor_centrality(G, eb, alpha=e["alpha"]),
        "cn_sh": lambda G, eb: nx.cn_soundarajan_hopcroft(G, eb, community="community"),
        "ra_sh": lambda G, eb: nx.ra_index_soundarajan_hopcroft(G, eb, community="community"),
        "wic": lambda G, eb: nx.within_inter_cluster(G, eb, delta=e["delta"],
                                                     community="community"),
    }
    out = list(funcs[m](g, pairs))
    if [(u, v) for u, v, _ in out] != pairs:
        raise RuntimeError("NetworkX returned the pairs in another order")
    return out


def main() -> None:
    spec = load_area("linkpred")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        g = simple_graph(fx)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            write_values(ORACLE, name, variant, scores(g, v["engine"]))
            count += 1
    print(f"networkx {nx.__version__} (linkpred): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
