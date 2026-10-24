#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for cut
measures (see ../bindings/cut_measures.py for the sets and encodings).

    python3 validation/oracles/cut_measures_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"


def build(fx):
    pairs = list(zip(fx.src, fx.tgt))
    keys = {p if fx.directed else frozenset(p) for p in pairs}
    multi = len(keys) < len(pairs)
    cls = {(True, True): nx.MultiDiGraph, (True, False): nx.DiGraph,
           (False, True): nx.MultiGraph, (False, False): nx.Graph}[(fx.directed, multi)]
    g = cls()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(pairs):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def ieee(f, num_den):
    """f() or, where NetworkX divides by zero, the IEEE value of num/den."""
    try:
        return float(f())
    except ZeroDivisionError:
        num, den = num_den()
        if den != 0:
            return num / den
        return math.nan if num == 0 or math.isnan(num) else math.copysign(math.inf, num)


def encode(x):
    return -999.0 if math.isnan(x) else (-1000.0 if math.isinf(x) else x)


def compute(fx, e):
    g = build(fx)
    s = {v for v in range(fx.n) if v % 2 == 0}
    t = (set(range(fx.n)) - s) if e.get("t", "complement") == "complement" else \
        {v for v in range(fx.n) if v % 3 == 1}
    w = "weight" if e.get("weighted") else None
    cut = nx.cut_size(g, s, t, weight=w)
    vs, vt = nx.volume(g, s, weight=w), nx.volume(g, t, weight=w)
    m = g.number_of_edges()
    nb = set()
    for v in s:
        nb |= set(g.neighbors(v))
    vals = [
        cut, vs, vt,
        # Replaced below when a volume is zero.
        ieee(lambda: nx.normalized_cut_size(g, s, t, weight=w), lambda: (math.nan, 1)),
        ieee(lambda: nx.conductance(g, s, t, weight=w), lambda: (cut, min(vs, vt))),
        ieee(lambda: nx.edge_expansion(g, s, t, weight=w), lambda: (cut, min(len(s), len(t)))),
        ieee(lambda: nx.mixing_expansion(g, s, t, weight=w), lambda: (cut, 2 * m)),
        ieee(lambda: nx.node_expansion(g, s), lambda: (len(nb), len(s))),
        ieee(lambda: nx.boundary_expansion(g, s), lambda: (len(nb - s), len(s))),
    ]
    if vs == 0 or vt == 0:
        # cut * (1/vs + 1/vt) with IEEE division.
        a = math.inf if vs == 0 else 1 / vs
        b = math.inf if vt == 0 else 1 / vt
        vals[3] = cut * (a + b) if cut != 0 else math.nan
    return [encode(x) for x in vals]


def main() -> None:
    spec = load_area("cut_measures")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            write_values(ORACLE, name, variant, compute(fx, v["engine"]))
            count += 1
    print(f"networkx {nx.__version__} (cut_measures): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
