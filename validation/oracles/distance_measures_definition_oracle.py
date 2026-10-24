#!/usr/bin/env python3
"""Golden values for the distance measures computed directly from their
definitions: exact distances by the Floyd-Warshall recurrence in rational
arithmetic (the engine uses breadth-first search, Dial's buckets, or
Dijkstra's algorithm), then each measure as a sum over pairs; closeness
vitality recomputes the distances of the network without each node; local
efficiency recomputes the distances inside each neighborhood.

Markers: NaN -999, infinity -1000, minus infinity -1001; node sets are
sorted node indices (empty when the network is not strongly connected).

    python3 validation/oracles/distance_measures_definition_oracle.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 80
UNDEFINED = -999.0
INFINITE = -1000.0
NEG_INFINITE = -1001.0


def lengths(fx, e):
    out = []
    for i in range(len(fx.src)):
        if not e.get("weighted"):
            out.append(Fraction(1))
        elif e.get("weights_as", "distance") == "strength":
            out.append(1 / Fraction(fx.weights[i]))
        else:
            out.append(Fraction(fx.weights[i]))
    return out


def distances(n, arcs, removed=None):
    d = [[None] * n for _ in range(n)]
    for v in range(n):
        if v != removed:
            d[v][v] = Fraction(0)
    for a, b, L in arcs:
        if a == b or removed in (a, b):
            continue
        if d[a][b] is None or L < d[a][b]:
            d[a][b] = L
    for k in range(n):
        for i in range(n):
            if d[i][k] is None:
                continue
            for j in range(n):
                if d[k][j] is not None and (d[i][j] is None or d[i][k] + d[k][j] < d[i][j]):
                    d[i][j] = d[i][k] + d[k][j]
    return d


def arcs_of(fx, e):
    lens = lengths(fx, e)
    mode = e.get("mode", "out")
    out = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if not fx.directed or mode == "all":
            out += [(a, b, lens[i]), (b, a, lens[i])]
        elif mode == "out":
            out.append((a, b, lens[i]))
        else:
            out.append((b, a, lens[i]))
    return out


def wiener(n, d, nodes, halve):
    total = Fraction(0)
    for i in nodes:
        for j in nodes:
            if i == j:
                continue
            if d[i][j] is None:
                return None
            total += d[i][j]
    return total / 2 if halve else total


_CACHE: dict = {}


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    arcs = arcs_of(fx, e)
    key = (fx.name, e.get("mode", "out"), bool(e.get("weighted")),
           e.get("weights_as", "distance"))
    if key not in _CACHE:
        _CACHE[key] = distances(n, arcs)
    d = _CACHE[key]
    connected = all(d[i][j] is not None for i in range(n) for j in range(n))
    undirected_reading = not fx.directed or e.get("mode", "out") == "all"
    if fn in ("dm_center", "dm_periphery"):
        if e.get("reachable", "require") == "require":
            if not connected:
                return []
            ecc = [max(d[i]) for i in range(n)]
        else:
            ecc = [max(x for x in d[i] if x is not None) for i in range(n)]
        best = min(ecc) if fn == "dm_center" else max(ecc)
        return [float(v) for v in range(n) if ecc[v] == best]
    if fn == "dm_barycenter":
        if not connected:
            return []
        tot = [sum(d[i]) for i in range(n)]
        return [float(v) for v in range(n) if tot[v] == min(tot)]
    if fn == "dm_index":
        index = e["index"]
        deg = [0] * n
        for a, b in zip(fx.src, fx.tgt):
            deg[a] += 1
            deg[b] += 1
        if index == "harmonic_diameter":
            inv = sum((1 / d[i][j] for i in range(n) for j in range(n)
                       if i != j and d[i][j] is not None and d[i][j] != 0), Fraction(0))
            if inv != 0:
                return [float(n * (n - 1) / inv)]
            return [INFINITE if n > 1 else UNDEFINED]
        if index == "global_efficiency":
            if n < 2:
                return [0.0]
            inv = sum((1 / d[i][j] for i in range(n) for j in range(n)
                       if i != j and d[i][j] is not None), Fraction(0))
            return [float(inv / (n * (n - 1)))]
        if n == 0:
            return [UNDEFINED]
        if not connected:
            return [INFINITE]
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if index == "wiener":
            return [float(wiener(n, d, range(n), undirected_reading))]
        if index == "schultz":
            return [float(sum((deg[i] + deg[j]) * d[i][j] for i, j in pairs))]
        if index == "gutman":
            return [float(sum(deg[i] * deg[j] * d[i][j] for i, j in pairs))]
        if index == "hyper_wiener":
            return [float(sum(d[i][j] + d[i][j] ** 2 for i, j in pairs) / 2)]
        if index == "hyper_wiener_networkx":
            return [float(sum(d[i][j] + d[i][j] ** 2 for i, j in pairs))]
    if fn == "dm_local_efficiency":
        if n == 0:
            return [UNDEFINED]
        nbrs = [set() for _ in range(n)]
        for a, b in zip(fx.src, fx.tgt):
            if a != b:
                nbrs[a].add(b)
                nbrs[b].add(a)
            elif e.get("loop_is_neighbor"):
                nbrs[a].add(a)
        total = Fraction(0)
        for v in range(n):
            members = sorted(nbrs[v])
            k = len(members)
            if k < 2:
                continue
            sub = [(a, b, Fraction(1)) for a in members for b in members
                   if a != b and (b in nbrs[a])]
            idx = {x: i for i, x in enumerate(members)}
            dd = distances(k, [(idx[a], idx[b], L) for a, b, L in sub])
            inv = sum((1 / dd[i][j] for i in range(k) for j in range(k)
                       if i != j and dd[i][j] is not None), Fraction(0))
            total += inv / (k * (k - 1))
        return [float(total / n)]
    if fn == "dm_closeness_vitality":
        if n <= 1:
            return [UNDEFINED] * n
        whole = wiener(n, d, range(n), undirected_reading)
        out = []
        for v in range(n):
            dv = distances(n, arcs, removed=v)
            after = wiener(n, dv, [u for u in range(n) if u != v], undirected_reading)
            if whole is None and after is None:
                out.append(UNDEFINED)
            elif whole is None:
                out.append(INFINITE)
            elif after is None:
                out.append(NEG_INFINITE)
            else:
                out.append(float(whole - after))
        return out
    raise ValueError(f"no definition oracle for {fn} {e}")


def main() -> int:
    manifest = load_area("distance_measures")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx):
                continue
            if name in ospec.get("undefined_on", []) or tags & set(ospec.get("undefined_when",
                                                                              [])):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (distance_measures): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
