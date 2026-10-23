#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
PageRank, Katz centrality, hubs and authorities, and assortativity
(validation/manifest.d/spectral.yaml).

Networks are built as MultiGraph or MultiDiGraph, so parallel edges add
up as they do in NetworkX's adjacency matrices. Node attributes come from
the rules described in validation/bindings/spectral.py (restated here).
Undefined values are written as empty lines.

    python3 validation/oracles/spectral_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (GOLDEN, MANIFEST_DIR, applies, fixture_names,  # noqa: E402
                    load_fixture, write_values)

ORACLE = "networkx"


def rule(name, n):
    return {
        "mod3": [v % 3 for v in range(n)],
        "mod3plus1": [v % 3 + 1 for v in range(n)],
        "num7": [7 * v % 11 for v in range(n)],
        "num3": [3 * v % 5 for v in range(n)],
    }[name]


def multigraph(fx, weighted=True):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if (fx.weighted and weighted) else 1.0)
    return g


def summed(fx, weighted):
    """A simple graph whose edge weights add up parallel edges (their
    weights, or their number), for functions that refuse multigraphs; its
    adjacency matrix is the multigraph's."""
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if (fx.weighted and weighted) else 1.0
        if g.has_edge(a, b):
            g[a][b]["weight"] += w
        else:
            g.add_edge(a, b, weight=w)
    return g


def defined(values):
    return [() if math.isnan(v) else v for v in values]


def loops_matter(fx):
    """Whether the count of an undirected self-loop changes the result."""
    return not fx.directed and any(a == b for a, b in zip(fx.src, fx.tgt))


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    nodes = range(n)
    g = multigraph(fx)
    weight = "weight" if e.get("weighted") else None
    if fn == "pagerank":
        if (e.get("loops") != "once" and loops_matter(fx)) or not e.get("directed", True):
            return None
        kwargs = {}
        if e.get("personalization"):
            kwargs["personalization"] = dict(zip(nodes, rule(e["personalization"], n)))
        if e.get("dangling") == "uniform":
            kwargs["dangling"] = {v: 1.0 for v in nodes}
        x = nx.pagerank(g, alpha=e["damping"], weight=weight, tol=1e-14, max_iter=20000,
                        **kwargs)
        return [x[v] for v in nodes]
    if fn == "katz":
        if (e.get("loops") != "once" and loops_matter(fx)) or n < 2:
            # katz_centrality_numpy fails on a single node (its result
            # squeezes to a scalar).
            return None
        s = summed(fx, bool(weight))
        a = nx.to_numpy_array(s, nodelist=list(nodes), weight="weight")
        lam = max(abs(np.linalg.eigvals(a))) if n else 0.0
        alpha = e["alpha_fraction"] / lam if lam > 0 else e["alpha_fraction"]
        x = nx.katz_centrality_numpy(s, alpha=alpha, beta=1.0, weight="weight",
                                     normalized=e.get("normalized", True))
        return [x[v] for v in nodes]
    if fn == "hits":
        if (e.get("loops") != "once" and loops_matter(fx)) or e.get("scale") != "sum":
            return None
        try:
            h, a = nx.hits(summed(fx, bool(weight)), max_iter=200000, tol=1e-13,
                           normalized=True)
        except nx.PowerIterationFailedConvergence:
            return None
        x = h if e["which"] == "hubs" else a
        return [x[v] for v in nodes]
    if fn == "assortativity_degree":
        if (e.get("loops") != "once" and loops_matter(fx)) or not e.get("directed", True):
            return None
        return [nx.degree_assortativity_coefficient(g, weight=weight)]
    if fn == "assortativity_nominal":
        if (e.get("loops") != "once" and loops_matter(fx)):
            return None
        nx.set_node_attributes(g, dict(zip(nodes, rule(e["types"], n))), "group")
        return [nx.attribute_assortativity_coefficient(g, "group")]
    if fn == "assortativity_values":
        if (e.get("loops") != "once" and loops_matter(fx)) or e.get("values_in"):
            return None
        nx.set_node_attributes(g, dict(zip(nodes, rule(e["values"], n))), "value")
        return [nx.numeric_assortativity_coefficient(g, "value")]
    return None


def main() -> None:
    spec = yaml.safe_load((MANIFEST_DIR / "spectral.yaml").read_text())
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            if name in v["oracles"][ORACLE].get("undefined_on", []):
                continue
            with np.errstate(all="ignore"):
                values = compute(v, fx)
            if values is not None:
                write_values(ORACLE, name, variant, defined([float(x) for x in values]))
                count += 1
    print(f"networkx {nx.__version__} (spectral): wrote {count} golden files")
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(f"networkx {nx.__version__}\n", newline="")


if __name__ == "__main__":
    main()
