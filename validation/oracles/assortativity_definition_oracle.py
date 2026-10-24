#!/usr/bin/env python3
"""Definition oracle for mixing counts, assortativity over node subsets,
average neighbour degree, and average degree connectivity, written from
the definitions (Newman 2002, 2003; Barrat et al. 2004; Pastor-Satorras et
al. 2001) over the fixture's edge list: every edge is an arc (both
directions when undirected, a loop once), degrees count edge ends (a loop
twice when undirected), and pair sums are exact fractions before the
final square root. Shares no code with the engine.

    python3 validation/oracles/assortativity_definition_oracle.py
"""

from __future__ import annotations

import math
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
UNDEFINED = -999.0


def arcs(fx):
    out = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if fx.weighted else Fraction(1)
        out.append((a, b, w))
        if not fx.directed and a != b:
            out.append((b, a, w))
    return out


def degrees(fx, kind, weighted):
    d = [Fraction(0)] * fx.n
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if (weighted and fx.weighted) else Fraction(1)
        if not fx.directed or kind == "in+out":
            d[a] += w
            d[b] += w
        elif kind == "out":
            d[a] += w
        else:
            d[b] += w
    return d


def correlation(pairs):
    if not pairs:
        return math.nan
    m = len(pairs)
    mx = sum(p[0] for p in pairs) / m
    my = sum(p[1] for p in pairs) / m
    cxy = sum((x - mx) * (y - my) for x, y in pairs)
    cxx = sum((x - mx) ** 2 for x, _ in pairs)
    cyy = sum((y - my) ** 2 for _, y in pairs)
    if cxx == 0 or cyy == 0:
        return math.nan
    return float(cxy) / math.sqrt(float(cxx) * float(cyy))


def nominal(pairs):
    if not pairs:
        return 0.0  # NetworkX's formula on an empty mixing matrix
    m = len(pairs)
    c = Counter(pairs)
    cats = {x for p in pairs for x in p}
    trace = sum(Fraction(c[(k, k)], m) for k in cats)
    ab = sum(Fraction(sum(c[(k, j)] for j in cats), m) * Fraction(sum(c[(j, k)] for j in cats), m)
             for k in cats)
    if ab == 1:
        return math.nan
    return float((trace - ab) / (1 - ab))


def compute(fx, e):
    n = fx.n
    inside = [v % 2 == 0 for v in range(n)] if e.get("subset") else [True] * n
    weighted = bool(e.get("weighted"))
    label = [v % 3 for v in range(n)]
    value = [(v * 7) % 5 for v in range(n)]
    fn = e["fn"]
    A = arcs(fx)
    if fn in ("mixing", "assort") and (e.get("kind") == "degree" or e.get("coefficient") == "degree"):
        dx = degrees(fx, e.get("x", "out"), weighted)
        dy = degrees(fx, e.get("y", "in"), weighted)
        pairs = [(dx[a], dy[b]) for a, b, _ in A if inside[a] and inside[b]]
        if fn == "mixing":
            return sorted((float(x), float(y), c) for (x, y), c in Counter(pairs).items())
        r = correlation(pairs)
        return [UNDEFINED if math.isnan(r) else r]
    if fn == "mixing":
        pairs = [(label[a], label[b]) for a, b, _ in A if inside[a]]
        return sorted((x, y, c) for (x, y), c in Counter(pairs).items())
    if fn == "assort":
        if e["coefficient"] == "attribute":
            r = nominal([(label[a], label[b]) for a, b, _ in A if inside[a]])
        else:
            known = {value[v] for v in range(n) if inside[v]}
            pairs = [(value[a], value[b]) for a, b, _ in A
                     if inside[a] and value[a] in known and value[b] in known]
            r = correlation(pairs)
        return [UNDEFINED if math.isnan(r) else r]
    if fn == "avg_nbr_deg":
        src = e.get("source", "out") if fx.directed else "in+out"
        tgt = e.get("target", "out") if fx.directed else "in+out"
        sdeg = degrees(fx, src, weighted)
        tdeg = degrees(fx, tgt, False)
        out = []
        for v in range(n):
            if not inside[v]:
                continue
            if sdeg[v] == 0:
                out.append(0.0)
                continue
            total = Fraction(0)
            for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
                w = Fraction(fx.weights[i]) if weighted else Fraction(1)
                if not fx.directed:
                    if a == v:
                        total += w * tdeg[b]
                    if b == v and a != b:
                        total += w * tdeg[a]
                else:
                    if a == v and "out" in src:
                        total += w * tdeg[b]
                    if b == v and "in" in src:
                        total += w * tdeg[a]
            out.append(float(total / sdeg[v]))
        return out
    if fn == "avg_deg_conn":
        src = e.get("source", "in+out") if fx.directed else "in+out"
        tgt = e.get("target", "in+out") if fx.directed else "in+out"
        k = degrees(fx, src, False)
        kw = degrees(fx, src, weighted)
        tdeg = degrees(fx, tgt, False)
        dsum, dnorm = {}, {}
        for v in range(n):
            if not inside[v]:
                continue
            s = Fraction(0)
            for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
                w = Fraction(fx.weights[i]) if weighted else Fraction(1)
                if not fx.directed:
                    if a == v:
                        s += w * tdeg[b]
                    if b == v and a != b:
                        s += w * tdeg[a]
                elif src == "in":
                    if b == v:
                        s += w * tdeg[a]
                elif a == v:
                    # "out", and "in+out" (successors, as NetworkX).
                    s += w * tdeg[b]
            dsum[k[v]] = dsum.get(k[v], 0) + s
            dnorm[k[v]] = dnorm.get(k[v], 0) + kw[v]
        return sorted((float(x), float(dsum[x] if dnorm[x] == 0 else dsum[x] / dnorm[x]))
                      for x in dsum)
    return None


def main() -> None:
    spec = load_area("assortativity")
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
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (assortativity): wrote {count} golden files")


if __name__ == "__main__":
    main()
