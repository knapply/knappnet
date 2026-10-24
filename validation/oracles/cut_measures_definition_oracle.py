#!/usr/bin/env python3
"""Definition oracle for cut measures: every formula evaluated over the
fixture's edge list in exact fractions, with IEEE results where a
denominator is zero. Shares no code with the engine.

    python3 validation/oracles/cut_measures_definition_oracle.py
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


def div(a, b):
    if b == 0:
        return math.nan if a == 0 else math.inf
    return float(Fraction(a) / Fraction(b))


def encode(x):
    return -999.0 if math.isnan(x) else (-1000.0 if math.isinf(x) else x)


def compute(fx, e):
    n = fx.n
    S = {v for v in range(n) if v % 2 == 0}
    T = (set(range(n)) - S) if e.get("t", "complement") == "complement" else \
        {v for v in range(n) if v % 3 == 1}
    weighted = bool(e.get("weighted"))
    cut = Fraction(0)
    vol = {"S": Fraction(0), "T": Fraction(0)}
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        if fx.directed:
            # Arcs S -> T and T -> S, each read as NetworkX's edge_boundary
            # reads the out-edges of its first set.
            if a in S and (b in T or (b in S and a in T)):
                cut += w
            if a in T and (b in S or (b in T and a in S)):
                cut += w
        elif (a in S and b in T) or (b in S and a in T):
            cut += w
        for key, X in (("S", S), ("T", T)):
            if a in X:
                vol[key] += w
            if not fx.directed and b in X:
                vol[key] += w
    nb = set()
    for a, b in zip(fx.src, fx.tgt):
        if a in S:
            nb.add(b)
        if not fx.directed and b in S:
            nb.add(a)
    vs, vt = vol["S"], vol["T"]
    inv = lambda x: math.inf if x == 0 else 1 / Fraction(x)  # noqa: E731
    if vs == 0 or vt == 0:
        ncut = math.nan if cut == 0 else math.inf
    else:
        ncut = float(cut * (inv(vs) + inv(vt)))
    vals = [float(cut), float(vs), float(vt), ncut, div(cut, min(vs, vt)),
            div(cut, min(len(S), len(T))), div(cut, 2 * len(fx.src)),
            div(len(nb), len(S)), div(len(nb - S), len(S))]
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
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (cut_measures): wrote {count} golden files")


if __name__ == "__main__":
    main()
