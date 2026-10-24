#!/usr/bin/env python3
"""Definition oracle for link prediction scores: each formula evaluated
directly on the symmetric 0/1 adjacency matrix of the fixture's
undirected simple graph (loops on the diagonal), in exact rational
arithmetic where the score is rational. Shares no code with the engine.

    python3 validation/oracles/linkpred_definition_oracle.py
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 80
INF = None


def adjacency(fx):
    a = [[0] * fx.n for _ in range(fx.n)]
    for s, t in zip(fx.src, fx.tgt):
        a[s][t] = a[t][s] = 1
    return a


def hops(a):
    n = len(a)
    d = [[0 if i == j else (1 if a[i][j] else INF) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if d[i][k] is INF:
                continue
            for j in range(n):
                if d[k][j] is not INF and (d[i][j] is INF or d[i][k] + d[k][j] < d[i][j]):
                    d[i][j] = d[i][k] + d[k][j]
    return d


def compute(fx, e, a, dist):
    n = fx.n
    deg = [sum(a[v]) + a[v][v] for v in range(n)]
    comm = [v % 3 for v in range(n)]
    if e["pairs"] == "non_edges":
        pairs = [(u, v) for u in range(n) for v in range(u + 1, n) if not a[u][v]]
    else:
        pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    m = e["method"]
    out = []
    for u, v in pairs:
        cn = [w for w in range(n) if w not in (u, v) and a[u][w] and a[v][w]]
        if m == "common_neighbors":
            s = len(cn)
        elif m == "jaccard":
            union = sum(1 for w in range(n) if a[u][w] or a[v][w])
            s = Fraction(len(cn), union) if union else 0
        elif m == "adamic_adar":
            s = math.fsum(1 / math.log(deg[w]) for w in cn)
        elif m == "resource_allocation":
            s = sum((Fraction(1, deg[w]) for w in cn), Fraction(0))
        elif m == "preferential_attachment":
            s = deg[u] * deg[v]
        elif m == "cnc":
            alpha = Fraction(e["alpha"]).limit_denominator(1000)
            d = dist[u][v]
            s = alpha * len(cn) + (0 if d is INF else (1 - alpha) * Fraction(n, d))
        elif m == "cn_sh":
            s = len(cn) + (sum(1 for w in cn if comm[w] == comm[u]) if comm[u] == comm[v] else 0)
        elif m == "ra_sh":
            s = (sum((Fraction(1, deg[w]) for w in cn if comm[w] == comm[u]), Fraction(0))
                 if comm[u] == comm[v] else 0)
        elif m == "wic":
            if comm[u] != comm[v]:
                s = 0
            else:
                within = sum(1 for w in cn if comm[w] == comm[u])
                s = Fraction(within) / (len(cn) - within + Fraction(1, 1000))
        else:
            return None
        out.append((u, v, float(s)))
    return out


def main() -> None:
    spec = load_area("linkpred")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        a = adjacency(fx)
        dist = hops(a)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"], a, dist)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (linkpred): wrote {count} golden files")


if __name__ == "__main__":
    main()
