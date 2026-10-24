#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
Weisfeiler-Lehman graph and subgraph hashes, SimRank, Panther similarity,
and graph edit distance, on the simple Graph or DiGraph of each fixture
(see ../bindings/similarity.py for the derived inputs and encodings).

    python3 validation/oracles/similarity_networkx_oracle.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,
                    write_version)  # noqa: E402

ORACLE = "networkx"
GED_MAX_NODES = 10


def hex_ints(text):
    return tuple(int(text[i:i + 8], 16) for i in range(0, len(text), 8))


def simple(fx, labels=False):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    for v in range(fx.n):
        g.add_node(v, label=str(v % 3))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, label=str(e % 2))
    return g


def multigraph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def hub(fx):
    nbrs = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    return max(range(fx.n), key=lambda v: (len(nbrs[v]), -v))


def ged_graphs(fx, labels):
    n = fx.n
    edges = list(zip(fx.src, fx.tgt))[1:]
    present = set(edges) | (set() if fx.directed else {(b, a) for a, b in edges})
    if (0, n - 1) not in present:
        edges.append((0, n - 1))
    edges.append((n, min(1, n - 1)))
    out = []
    for nn, ee in ((n, list(zip(fx.src, fx.tgt))), (n + 1, edges)):
        g = nx.DiGraph() if fx.directed else nx.Graph()
        for v in range(nn):
            g.add_node(v, label=v % 2)
        for e, (a, b) in enumerate(ee):
            g.add_edge(a, b, label=e % 2)
        out.append(g)
    return out


def compute(fx, e):
    fn = e["fn"]
    if fn in ("wl_hash", "wl_subgraph"):
        g = simple(fx)
        kw = dict(iterations=int(e.get("iterations", 3)), digest_size=int(e.get("digest_size", 16)),
                  node_attr="label" if e.get("node_labels") else None,
                  edge_attr="label" if e.get("edge_labels") else None)
        if fn == "wl_hash":
            return [hex_ints(nx.weisfeiler_lehman_graph_hash(g, **kw))]
        hashes = nx.weisfeiler_lehman_subgraph_hashes(
            g, include_initial_labels=bool(e.get("include_initial")), **kw)
        return [hex_ints("".join(hashes[v])) for v in range(fx.n)]
    if fn == "simrank":
        tol = float(e.get("tolerance", 1e-4))
        if e.get("stop", "allclose") == "relative":
            s = nx.similarity._simrank_similarity_python(simple(fx), tolerance=tol)
        elif e.get("weighted"):
            s = nx.simrank_similarity(multigraph(fx), tolerance=tol)
        else:
            g = nx.DiGraph() if fx.directed else nx.Graph()
            g.add_nodes_from(range(fx.n))
            g.add_edges_from(zip(fx.src, fx.tgt))
            s = nx.simrank_similarity(g, tolerance=tol)
        return [float(s[u][v]) for u in range(fx.n) for v in range(fx.n)]
    if fn == "panther":
        g = nx.Graph()
        g.add_nodes_from(range(fx.n))
        g.add_edges_from(zip(fx.src, fx.tgt))
        src = hub(fx)
        kept = [v for v in g if g.degree(v) > 0]
        r = nx.panther_similarity(g, src, k=len(kept) - 1, eps=float(e["eps"]),
                                  seed=int(e.get("seed", 1)))
        out = [-1.0] * fx.n
        for v, x in r.items():
            out[v] = x
        return out
    if fn == "ged":
        if fx.n > GED_MAX_NODES:
            return None
        g1, g2 = ged_graphs(fx, e["labels"])
        lab = e.get("labels", "none")
        same = lambda a, b: a["label"] == b["label"]  # noqa: E731
        d = nx.graph_edit_distance(g1, g2, node_match=same if lab in ("nodes", "both") else None,
                                   edge_match=same if lab == "both" else None)
        return [d]
    return None


def undefined(fx, spec):
    o = spec["oracles"][ORACLE]
    tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
    return fx.name in o.get("undefined_on", []) or tags & set(o.get("undefined_when", []))


def main() -> None:
    warnings.simplefilter("ignore")
    spec = load_area("similarity")
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
    print(f"networkx {nx.__version__} (similarity): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
