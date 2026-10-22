#!/usr/bin/env python3
"""Golden values for brokerage and positions, computed from the definitions.

Shares no code with the engine. Each measure is evaluated the slow, obvious
way on a dense matrix:

- Constraint and effective size (Burt 1992) follow Burt's formulas term by
  term in exact rational arithmetic, over each node's contacts.
- Efficiency is effective size divided by the number of contacts.
- Gould-Fernandez brokerage (Gould and Fernandez 1989) enumerates every
  ordered triple (i, j, k) of distinct nodes with i -> j -> k and no
  i -> k, and classifies it by the groups of i, j, and k.
- Structural equivalence compares the stacked row and column of each pair
  of nodes over k != i, j (Burt 1976; Wasserman and Faust 1994), with exact
  sums; correlation is evaluated from exact central moments.
- CONCOR (Breiger, Boorman, and Arabie 1975) iterates numpy correlation
  matrices; see the conventions in crates/engine/src/metrics/equivalence.rs.
  Where a correlation with the reference node comes within 1e-6 of zero
  after the first iteration, rounding decides the partition, so no golden
  file is written for that fixture (it is reported as skipped).

Undefined values (NaN) are written as empty lines.

    python3 validation/oracles/positions_definition_oracle.py
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, applies, fixture_names, load_fixture, write_values  # noqa: E402
from common import load_manifest  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80
NEAR_ZERO = 1e-6
STALL = 1e-12
FUNCTIONS = {"constraint", "effective_size", "efficiency", "brokerage",
             "structural_equivalence", "concor"}


def groups_for(rule: str, n: int) -> list[int]:
    if rule == "mod3":
        return [v % 3 for v in range(n)]
    raise ValueError(rule)


def tie_matrix(fx, weighted: bool):
    """a[i][k]: total weight (or count, or 0/1) of edges i -> k, no loops.
    Undirected edges are entered in both directions."""
    n = fx.n
    a = [[Fraction(0)] * n for _ in range(n)]
    for e, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        w = Fraction(fx.weights[e]) if weighted else Fraction(1)
        a[s][t] += w
        if not fx.directed:
            a[t][s] += w
    return a


def binary(a):
    return [[Fraction(1) if x != 0 else Fraction(0) for x in row] for row in a]


# ---- Structural holes --------------------------------------------------------

def holes(fx, weighted: bool, which: str) -> list:
    """Burt (1992): p_ij = z_ij / sum_k z_ik with z_ij = a_ij + a_ji, where
    parallel edges add up (weights as given, or 1 per edge)."""
    n = fx.n
    a = [[Fraction(0)] * n for _ in range(n)]
    for e, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        w = Fraction(fx.weights[e]) if weighted else Fraction(1)
        a[s][t] += w
        if not fx.directed:
            a[t][s] += w
    z = [[a[i][j] + a[j][i] for j in range(n)] for i in range(n)]
    contacts = [[j for j in range(n) if j != i and z[i][j] != 0] for i in range(n)]
    total = [sum(z[i]) for i in range(n)]
    p = [[z[i][j] / total[i] if total[i] else Fraction(0) for j in range(n)] for i in range(n)]
    out = []
    for i in range(n):
        c = contacts[i]
        if not c:
            out.append(())
            continue
        if which == "constraint":
            val = Fraction(0)
            for j in c:
                term = p[i][j] + sum((p[i][q] * p[q][j] for q in c if q != j), Fraction(0))
                val += term * term
            out.append(float(val))
            continue
        es = Fraction(0)
        for j in c:
            mx = max(z[j][k] for k in range(n) if k != j)
            red = sum((p[i][q] * z[j][q] / mx for q in c if q != j), Fraction(0))
            es += 1 - red
        out.append(float(es) if which == "effective_size" else float(es / len(c)))
    return out


# ---- Brokerage ---------------------------------------------------------------

def brokerage(fx, rule: str) -> list:
    n = fx.n
    arc = set()
    for s, t in zip(fx.src, fx.tgt):
        if s != t:
            arc.add((s, t))
            if not fx.directed:
                arc.add((t, s))
    grp = groups_for(rule, n)
    counts = [[0] * 5 for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if len({i, j, k}) < 3:
                    continue
                if (i, j) not in arc or (j, k) not in arc or (i, k) in arc:
                    continue
                gi, gj, gk = grp[i], grp[j], grp[k]
                if gi == gj == gk:
                    role = 0  # coordinator: A A A
                elif gi == gk:
                    role = 1  # itinerant: A B A
                elif gi == gj:
                    role = 2  # representative: A A B
                elif gj == gk:
                    role = 3  # gatekeeper: A B B
                else:
                    role = 4  # liaison: A B C
                counts[j][role] += 1
    return [tuple(c + [sum(c)]) for c in counts]


# ---- Structural equivalence --------------------------------------------------

def profile_pairs(a, i, j):
    n = len(a)
    ks = [k for k in range(n) if k != i and k != j]
    x = [v for k in ks for v in (a[i][k], a[k][i])]
    y = [v for k in ks for v in (a[j][k], a[k][j])]
    return x, y


def exact_correlation(x, y) -> float:
    m = len(x)
    if m < 2:
        return 0.0
    mx = sum(x, Fraction(0)) / m
    my = sum(y, Fraction(0)) / m
    sxy = sum(((u - mx) * (v - my) for u, v in zip(x, y)), Fraction(0))
    sxx = sum(((u - mx) ** 2 for u in x), Fraction(0))
    syy = sum(((v - my) ** 2 for v in y), Fraction(0))
    if sxx == 0 or syy == 0:
        return 0.0
    # sxy / sqrt(sxx * syy), with a single rounding of the square root.
    return float(sxy) / math.sqrt(float(sxx * syy)) if sxy else 0.0


def structural_equivalence(fx, method: str, weighted: bool) -> list:
    a = tie_matrix(fx, weighted)
    if not weighted:
        a = binary(a)
    n = fx.n
    rows = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                if method != "correlation":
                    row.append(0.0)
                    continue
                # sna's diagonal: the profile against itself over k != i.
                ks = [k for k in range(n) if k != i]
                x = [v for k in ks for v in (a[i][k], a[k][i])]
                row.append(1.0 if len(x) >= 2 and len(set(x)) > 1 else 0.0)
                continue
            x, y = profile_pairs(a, i, j)
            if method == "euclidean":
                row.append(math.sqrt(sum((u - v) ** 2 for u, v in zip(x, y))))
            elif method == "hamming":
                row.append(float(sum(abs(u - v) for u, v in zip(x, y))))
            else:
                row.append(exact_correlation(x, y))
        rows.append(tuple(row))
    return rows


# ---- CONCOR ------------------------------------------------------------------

def concor(fx, depth: int, weighted: bool, tol: float, max_iter: int) -> list:
    a = tie_matrix(fx, weighted)
    if not weighted:
        a = binary(a)
    n = fx.n
    A = np.array([[float(x) for x in row] for row in a]) if n else np.zeros((0, 0))
    # Profiles with the diagonal missing: row then column of each node.
    np.fill_diagonal(A, np.nan) if n else None
    P = np.vstack([A, A.T])  # column v is node v's profile
    constant = {}
    active = []
    for v in range(n):
        vals = P[:, v][~np.isnan(P[:, v])]
        if vals.size == 0 or np.all(vals == vals[0]):
            constant.setdefault(float(vals[0]) if vals.size else 0.0, []).append(v)
        else:
            active.append(v)
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mask = ~np.isnan(P[:, i]) & ~np.isnan(P[:, j])
            x, y = P[mask, i], P[mask, j]
            if x.size >= 2 and x.std() > 0 and y.std() > 0:
                C[i, j] = np.corrcoef(x, y)[0, 1]

    fragile = False

    def converge(M):
        nonlocal fragile
        for _ in range(max_iter):
            if np.all(np.abs(M) >= 1 - tol):
                break
            with np.errstate(invalid="ignore", divide="ignore"):
                N = np.nan_to_num(np.corrcoef(M, rowvar=False), nan=0.0)
            step = np.abs(N - M).max()
            M = N
            fragile |= bool(np.any(np.abs(M[0]) < NEAR_ZERO))
            if step <= STALL:
                break  # a fixed point that is not +1 or -1
        return M

    blocks = [active]
    for _ in range(depth):
        nxt = []
        for b in blocks:
            if len(b) < 2:
                if b:
                    nxt.append(b)
                continue
            M = converge(C[np.ix_(b, b)])
            first = [v for c, v in enumerate(b) if M[0, c] > 0]
            second = [v for c, v in enumerate(b) if not M[0, c] > 0]
            nxt += [part for part in (first, second) if part]
        blocks = nxt
    if fragile:
        return None
    blocks += [constant[k] for k in constant]
    member = [0] * n
    for label, b in enumerate(blocks, start=1):
        for v in b:
            member[v] = label
    return [float(m) for m in member]


def main() -> None:
    manifest = load_manifest()
    count = 0
    skipped = []
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, spec in manifest.items():
            e = spec["engine"]
            if e["fn"] not in FUNCTIONS:
                continue
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            fn = e["fn"]
            if fn in ("constraint", "effective_size", "efficiency"):
                values = holes(fx, e["weighted"], fn)
            elif fn == "brokerage":
                values = brokerage(fx, e["groups"])
            elif fn == "structural_equivalence":
                values = structural_equivalence(fx, e["method"], e.get("weighted", False))
            else:
                values = concor(fx, e["depth"], e.get("weighted", False),
                                e.get("tolerance", 1e-7), e.get("max_iter", 100))
            if values is None:
                skipped.append(f"{variant} on {name}")
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    version = GOLDEN / ORACLE / "VERSION"
    if not version.exists():
        version.parent.mkdir(parents=True, exist_ok=True)
        version.write_text("definition (direct evaluation, exact arithmetic)\n", newline="")
    print(f"definition (positions): wrote {count} golden files")
    if skipped:
        print("   numerically fragile CONCOR splits, no golden written: " + ", ".join(skipped))


if __name__ == "__main__":
    main()
