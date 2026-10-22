#!/usr/bin/env python3
"""Validate file formats against NetworkX in both directions.

For every fixture, plus an attribute-rich network with awkward values
(Unicode, quotes, missing values, isolates, self-loops, parallel edges):

  ours -> NetworkX   write with this engine, read with NetworkX
  NetworkX -> ours   write with NetworkX, read with this engine

and compare node identifiers, edges (as a multiset, with direction),
weights, and node attributes. Formats that cannot carry a feature (for
example, attributes in Pajek) are only checked for what they can carry.

    python3 validation/io_check.py
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
from collections import Counter
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import fixture_names, load_fixture  # noqa: E402
from engine import read_network, write_network  # noqa: E402

# Formats NetworkX can both read and write, and what each carries.
# Pajek has no attribute columns, treats a missing link value as 1, and has
# no way to write a double quote inside a label.
FORMATS = {
    "graphml": {"attrs": True, "ext": "graphml"},
    "gexf": {"attrs": True, "ext": "gexf"},
    "gml": {"attrs": True, "ext": "gml"},
    "pajek": {"attrs": False, "ext": "net", "default_weight": 1.0, "no_quotes": True},
    "json": {"attrs": True, "ext": "json"},
}
IGNORED_NODE_ATTRS = {"label", "id", "x", "y", "z", "shape"}


def rich_network() -> dict:
    return {
        "directed": True,
        "graph": {"name": "rich"},
        "nodes": {
            "id": ["ana", "bo \"b\"", "céline", "d&d", "isolate"],
            "age": [31, 45, None, 22, 60],
            "score": [0.5, None, 1.25, -3.0, 1e-7],
            "member": [True, False, None, True, False],
            "team": ["red", "blue", "<blue>", None, "green & gold"],
        },
        "edges": {
            "from": ["ana", "bo \"b\"", "céline", "ana", "d&d", "ana"],
            "to": ["bo \"b\"", "céline", "céline", "d&d", "ana", "bo \"b\""],
            "weight": [1.0, 2.5, 0.75, 3.0, 4.0, 2.0],
        },
    }


def without_quotes(net: dict) -> dict:
    fix = lambda v: v.replace('"', "'") if isinstance(v, str) else v  # noqa: E731
    out = {"directed": net["directed"], "graph": net["graph"]}
    for t in ("nodes", "edges"):
        out[t] = {k: ([fix(x) for x in v] if k in ("id", "from", "to") else v)
                  for k, v in net[t].items()}
    return out


def fixture_network(name: str) -> dict:
    fx = load_fixture(name)
    edges = {"from": [str(s) for s in fx.src], "to": [str(t) for t in fx.tgt]}
    if fx.weighted:
        edges["weight"] = list(fx.weights)
    return {"directed": fx.directed, "graph": {},
            "nodes": {"id": [str(i) for i in range(fx.n)]}, "edges": edges}


def to_networkx(net: dict) -> nx.Graph:
    pairs = list(zip(net["edges"]["from"], net["edges"]["to"]))
    multi = len(set(pairs if net["directed"] else [tuple(sorted(p)) for p in pairs])) < len(pairs)
    cls = {(True, True): nx.MultiDiGraph, (True, False): nx.DiGraph,
           (False, True): nx.MultiGraph, (False, False): nx.Graph}[(net["directed"], multi)]
    g = cls(**net.get("graph", {}))
    cols = net["nodes"]
    for i, node in enumerate(cols["id"]):
        g.add_node(node, **{k: v[i] for k, v in cols.items() if k != "id" and v[i] is not None})
    e = net["edges"]
    for i, (a, b) in enumerate(pairs):
        g.add_edge(a, b, **{k: v[i] for k, v in e.items() if k not in ("from", "to") and v[i] is not None})
    return g


def edge_multiset(pairs, weights, directed, default_weight=None) -> Counter:
    out = Counter()
    for (a, b), w in zip(pairs, weights):
        w = default_weight if w is None else w
        key = (a, b) if directed else tuple(sorted((a, b)))
        out[key + (None if w is None else round(float(w), 9),)] += 1
    return out


def summarize_ours(net: dict, default_weight=None):
    nodes = net["nodes"]["id"]
    pairs = list(zip(net["edges"]["from"], net["edges"]["to"]))
    weights = net["edges"].get("weight", [None] * len(pairs))
    attrs = {}
    for i, node in enumerate(nodes):
        attrs[node] = {k: v[i] for k, v in net["nodes"].items()
                       if k not in IGNORED_NODE_ATTRS and v[i] is not None}
    return set(nodes), edge_multiset(pairs, weights, net["directed"], default_weight), attrs


def summarize_nx(g: nx.Graph, default_weight=None):
    nodes = {str(n) for n in g.nodes}
    pairs, weights = [], []
    for a, b, d in g.edges(data=True):
        pairs.append((str(a), str(b)))
        weights.append(d.get("weight"))
    attrs = {str(n): {k: v for k, v in d.items() if k not in IGNORED_NODE_ATTRS}
             for n, d in g.nodes(data=True)}
    return nodes, edge_multiset(pairs, weights, g.is_directed(), default_weight), attrs


def same_value(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool):
        return int(a) == int(b)
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)
    return str(a) == str(b)


def compare(expected, actual, check_attrs: bool) -> list[str]:
    (en, ee, ea), (an, ae, aa) = expected, actual
    problems = []
    if en != an:
        problems.append(f"nodes differ: missing {sorted(en - an)[:3]}, extra {sorted(an - en)[:3]}")
    if ee != ae:
        problems.append(f"edges differ: missing {list((ee - ae).items())[:2]}, "
                        f"extra {list((ae - ee).items())[:2]}")
    if check_attrs:
        for node, want in ea.items():
            got = aa.get(node, {})
            for k, v in want.items():
                if k not in got or not same_value(v, got[k]):
                    problems.append(f"node {node!r} attribute {k!r}: want {v!r}, got {got.get(k)!r}")
                    break
    return problems


def nx_read(fmt: str, path: Path) -> nx.Graph:
    if fmt == "graphml":
        return nx.read_graphml(path)
    if fmt == "gexf":
        return nx.read_gexf(path)
    if fmt == "gml":
        return nx.read_gml(path, label="label")
    if fmt == "pajek":
        return nx.read_pajek(path)
    return nx.node_link_graph(json.loads(path.read_text(encoding="utf-8")))


def nx_write(fmt: str, g: nx.Graph, path: Path) -> None:
    if fmt == "graphml":
        nx.write_graphml(g, path)
    elif fmt == "gexf":
        nx.write_gexf(g, path)
    elif fmt == "gml":
        nx.write_gml(g, path)
    elif fmt == "pajek":
        nx.write_pajek(g, path)
    else:
        path.write_text(json.dumps(nx.node_link_data(g)), encoding="utf-8")


def main() -> int:
    cases = {"rich": rich_network()}
    cases.update({name: fixture_network(name) for name in fixture_names()})
    failures, checks = [], 0
    with tempfile.TemporaryDirectory() as tmp:
        for case, net in cases.items():
            for fmt, spec in FORMATS.items():
                if spec.get("no_quotes"):
                    net_f = without_quotes(net)
                else:
                    net_f = net
                dw = spec.get("default_weight")
                expected = summarize_ours(net_f, dw)
                path = Path(tmp) / f"{case}.{spec['ext']}"
                write_network(net_f, path, fmt)
                got = summarize_nx(nx_read(fmt, path), dw)
                for p in compare(expected, got, spec["attrs"]):
                    failures.append(f"ours -> networkx [{fmt}] {case}: {p}")
                checks += 1

                path2 = Path(tmp) / f"{case}.nx.{spec['ext']}"
                g = to_networkx(net_f)
                nx_write(fmt, g, path2)
                ours = summarize_ours(read_network(path2, fmt), dw)
                for p in compare(summarize_nx(g, dw), ours, spec["attrs"]):
                    failures.append(f"networkx -> ours [{fmt}] {case}: {p}")
                checks += 1

    for f in failures:
        print("FAIL", f)
    print(f"{checks} format checks, {len(failures)} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
