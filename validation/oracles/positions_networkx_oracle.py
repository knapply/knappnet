#!/usr/bin/env python3
"""Golden values for Burt's structural-hole measures from NetworkX
(Hagberg, Schult, and Swart 2008).

Constraint uses NetworkX's matrix path (`nodes=None`), which sums parallel
edges and evaluates Burt's formula directly.

Effective size cannot be taken from NetworkX 3.7 as-is on these fixtures:
- the matrix path divides each mutual weight z_jq by the largest tie of q
  rather than of j (the vector of row maxima broadcasts across columns),
  which changes results whenever nodes' strongest ties differ;
- the per-node path reads multigraph weights as 1, and in a directed graph
  gives NaN to a node with in-arcs but no out-arcs.
Burt's effective size depends only on the mutual ties z = A + A^T, so the
per-node path is called on an undirected simple graph whose edge weight is
z (self-loops included, as NetworkX would see them). That evaluates
NetworkX's own redundancy formula on exactly the ties it is defined on.

NetworkX has no efficiency function; its documentation defines efficiency
as effective size divided by degree, which is computed here with the
number of distinct contacts (in- or out-neighbors other than the node).
Undefined values (NaN) are written as empty lines.

    python3 validation/oracles/positions_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_fixture, load_manifest,
                    write_values, write_version)  # noqa: E402

ORACLE = "networkx"
FUNCTIONS = {"constraint", "effective_size", "efficiency"}


def multigraph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def mutual_graph(fx, weighted: bool):
    """Undirected simple graph with weight z_ij = a_ij + a_ji."""
    z = {}
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if weighted else 1.0
        key = (min(a, b), max(a, b))
        z[key] = z.get(key, 0.0) + w * (1 if fx.directed or a == b else 2)
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    for (a, b), w in z.items():
        g.add_edge(a, b, z=w)
    return g


def compute(fn: str, fx, weighted: bool) -> list:
    if fn == "constraint":
        d = nx.constraint(multigraph(fx), weight="weight" if weighted else None)
    else:
        g = mutual_graph(fx, weighted)
        d = nx.effective_size(g, nodes=list(range(fx.n)), weight="z")
        if fn == "efficiency":
            k = {v: len(set(g[v]) - {v}) for v in g}
            d = {v: d[v] / k[v] if k[v] else math.nan for v in d}
    return [() if math.isnan(d[v]) else d[v] for v in range(fx.n)]


def main() -> None:
    manifest = load_manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            e = spec["engine"]
            if e["fn"] not in FUNCTIONS:
                continue
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            write_values(ORACLE, name, variant, compute(e["fn"], fx, e["weighted"]))
            count += 1
    print(f"networkx {nx.__version__} (positions): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
