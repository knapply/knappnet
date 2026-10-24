#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
mixing dicts, assortativity coefficients over node subsets, average
neighbour degree, and average degree connectivity (see
../bindings/assortativity.py for the derived inputs and encodings).

    python3 validation/oracles/assortativity_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"
UNDEFINED = -999.0


def build(fx):
    multi = len(set(zip(fx.src, fx.tgt))) < len(fx.src) or (
        not fx.directed and len({frozenset(p) for p in zip(fx.src, fx.tgt)}) < len(fx.src))
    cls = {(True, True): nx.MultiDiGraph, (True, False): nx.DiGraph,
           (False, True): nx.MultiGraph, (False, False): nx.Graph}[(fx.directed, multi)]
    g = cls()
    for v in range(fx.n):
        g.add_node(v, label=v % 3, value=float((v * 7) % 5))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def flat(d):
    return sorted((x, y, c) for x, row in d.items() for y, c in row.items() if c)


def compute(fx, e):
    g = build(fx)
    nodes = [v for v in range(fx.n) if v % 2 == 0] if e.get("subset") else None
    w = "weight" if e.get("weighted") else None
    fn = e["fn"]
    if fn == "mixing":
        if e["kind"] == "degree":
            return flat(nx.degree_mixing_dict(g, x=e.get("x", "out"), y=e.get("y", "in"),
                                              weight=w, nodes=nodes))
        return flat(nx.attribute_mixing_dict(g, "label", nodes=nodes))
    if fn == "assort":
        c = e["coefficient"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                if c == "degree":
                    if "pearson" in e.get("variant", ""):
                        r = nx.degree_pearson_correlation_coefficient(g)
                    else:
                        r = nx.degree_assortativity_coefficient(
                            g, x=e.get("x", "out"), y=e.get("y", "in"), weight=w, nodes=nodes)
                elif c == "attribute":
                    r = nx.attribute_assortativity_coefficient(g, "label", nodes=nodes)
                else:
                    r = nx.numeric_assortativity_coefficient(g, "value", nodes=nodes)
            except (ZeroDivisionError, ValueError, IndexError):
                r = math.nan
        return [UNDEFINED if math.isnan(r) else r]
    if fn == "avg_nbr_deg":
        kw = {}
        if fx.directed:
            kw = dict(source=e.get("source", "out"), target=e.get("target", "out"))
        d = nx.average_neighbor_degree(g, nodes=nodes, weight=w, **kw)
        return [d[v] for v in (nodes if nodes is not None else range(fx.n))]
    if fn == "avg_deg_conn":
        kw = {}
        if fx.directed:
            kw = dict(source=e.get("source", "in+out"), target=e.get("target", "in+out"))
        d = nx.average_degree_connectivity(g, nodes=nodes, weight=w, **kw)
        return sorted((float(k), v) for k, v in d.items())
    return None


def main() -> None:
    spec = load_area("assortativity")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            e = dict(v["engine"])
            if variant.endswith("_pearson"):
                e["variant"] = "pearson"
            values = compute(fx, e)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (assortativity): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
