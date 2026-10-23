#!/usr/bin/env python3
"""Golden values for temporal networks, computed directly from the
definitions.

Every fixture is read as a contact sequence with the times of
validation/temporal_fixtures.py. A time-respecting walk is a path in the
*arc graph*: its nodes are the temporal arcs (one per directed contact,
two per undirected contact, none for self-loops), and arc b may follow arc
a when b leaves the node where a arrives and the waiting time
dep(b) - arr(a) satisfies the waiting rule. Optimal walks are found by
breadth-first search in that graph from every possible first arc, in
exact rational arithmetic, so walks that revisit nodes are included and
nothing depends on the order in which contacts are processed. The engine
instead makes one pass over contacts sorted by time (Wu et al. 2014).

Temporal closeness (Pan and Saramaki 2011) integrates the latency over
starting times piece by piece between consecutive departures, rather than
through the engine's fastest-path profile and Eq. 1. Burstiness and memory
(Goh and Barabasi 2008) use exact rational moments.

Infinities are written as -1000 and NaN as -999, as the engine runner does.

    python3 validation/oracles/temporal_definition_oracle.py
"""

from __future__ import annotations

import math
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, MANIFEST_DIR, applies, fixture_names, load_fixture  # noqa: E402
from common import write_values  # noqa: E402
from temporal_fixtures import contact_times  # noqa: E402

ORACLE = "definition"
UNREACHABLE = -1000.0
UNDEFINED = -999.0


def manifest() -> dict:
    return yaml.safe_load((MANIFEST_DIR / "temporal.yaml").read_text())


def mark(v) -> float:
    if v is None:
        return UNREACHABLE
    if isinstance(v, float) and math.isnan(v):
        return UNDEFINED
    return float(v)


# ---- Contacts, arcs, and the waiting rule ------------------------------------

class Model:
    def __init__(self, fx, e):
        start, end = contact_times(fx, e.get("instantaneous", False))
        self.fx = fx
        self.start = [Fraction(int(s)) for s in start]
        self.end = [Fraction(int(x)) for x in end]
        traversal = e.get("traversal", "zero")
        self.min_wait = Fraction(e.get("min_wait", 0))
        self.max_wait = None if "max_wait" not in e else Fraction(e["max_wait"])
        self.strict = e.get("strict", False)

        def lam(c):
            if traversal == "zero":
                return Fraction(0)
            if traversal == "duration":
                return self.end[c] - self.start[c]
            return Fraction(traversal)

        m = len(fx.src)
        self.arcs = []  # (from, to, dep, arr)
        for c in range(m):
            u, v = fx.src[c], fx.tgt[c]
            if u == v:
                continue
            arr = self.start[c] + lam(c)
            self.arcs.append((u, v, self.start[c], arr))
            if not fx.directed:
                self.arcs.append((v, u, self.start[c], arr))
        if "start" in e:
            self.lo = Fraction(e["start"])
        else:
            self.lo = min(self.start) if m else Fraction(0)
        if "end" in e:
            self.hi = Fraction(e["end"])
        else:
            self.hi = max(self.start[c] + lam(c) for c in range(m)) if m else self.lo
        by_tail = {}
        for j, a in enumerate(self.arcs):
            by_tail.setdefault(a[0], []).append(j)
        self.next = []
        for a in self.arcs:
            nxt = []
            for j in by_tail.get(a[1], []):
                w = self.arcs[j][2] - a[3]
                ok = w > self.min_wait if self.strict else w >= self.min_wait
                if ok and (self.max_wait is None or w <= self.max_wait):
                    nxt.append(j)
            self.next.append(nxt)
        self._reach = {}

    def bfs(self, starts):
        """Hop count (1 for the starts) of every arc reachable from `starts`."""
        depth = {}
        q = deque()
        for s in starts:
            if s not in depth:
                depth[s] = 1
                q.append(s)
        while q:
            a = q.popleft()
            for b in self.next[a]:
                if b not in depth:
                    depth[b] = depth[a] + 1
                    q.append(b)
        return depth

    def starts(self, i):
        return [k for k, a in enumerate(self.arcs) if a[0] == i and a[2] >= self.lo]

    def reach(self, a):
        """For walks whose first arc is a: earliest arrival at each node within the window."""
        if a not in self._reach:
            best = {}
            for b in self.bfs([a]):
                _, v, _, arr = self.arcs[b]
                if arr <= self.hi and (v not in best or arr < best[v]):
                    best[v] = arr
            self._reach[a] = best
        return self._reach[a]


# ---- Paths -------------------------------------------------------------------

def matrix(md: Model, criterion: str):
    n = md.fx.n
    out = []
    for i in range(n):
        starts = md.starts(i)
        row = [None] * n
        if criterion == "shortest":
            for b, d in md.bfs(starts).items():
                _, v, _, arr = md.arcs[b]
                if arr <= md.hi and (row[v] is None or d < row[v]):
                    row[v] = d
        else:
            for a in starts:
                dep = md.arcs[a][2]
                for v, arr in md.reach(a).items():
                    val = {"earliest_arrival": arr, "fastest": arr - dep,
                           "latest_departure": dep}[criterion]
                    better = (row[v] is None
                              or (val > row[v] if criterion == "latest_departure"
                                  else val < row[v]))
                    if better:
                        row[v] = val
        row[i] = {"earliest_arrival": md.lo, "latest_departure": md.hi}.get(criterion, 0)
        out.extend(row)
    return out


def reach_counts(md: Model, direction: str):
    n = md.fx.n
    m = matrix(md, "shortest")
    if direction == "out":
        return [sum(1 for j in range(n) if j != i and m[i * n + j] is not None)
                for i in range(n)]
    return [sum(1 for i in range(n) if i != j and m[i * n + j] is not None)
            for j in range(n)]


def average_distance(md: Model, i: int, j: int):
    """Pan and Saramaki's average temporal distance from i to j."""
    pairs = []
    for a in md.starts(i):
        r = md.reach(a)
        if j in r:
            pairs.append((md.arcs[a][2], r[j]))
    if not pairs:
        return None
    span = md.hi - md.lo
    first = min(p[1] for p in pairs)
    cuts = sorted({md.lo, md.hi} | {p[0] for p in pairs})
    total = Fraction(0)
    for x, y in zip(cuts, cuts[1:]):
        # Starting in (x, y], the earliest arrival comes from walks leaving
        # at or after y, or else from the first walk of the next period.
        arrive = min([p[1] for p in pairs if p[0] >= y] + [first + span])
        total += (y - x) * (arrive - (x + y) / 2)
    return total / span


def closeness(md: Model, distances: bool):
    n = md.fx.n
    if n < 2 or md.hi <= md.lo:
        if distances:
            return [0.0 if i == j else math.nan for i in range(n) for j in range(n)]
        return [math.nan] * n
    tau = [[Fraction(0) if i == j else average_distance(md, i, j) for j in range(n)]
           for i in range(n)]
    if distances:
        return [tau[i][j] for i in range(n) for j in range(n)]
    return [sum((1 / tau[i][j] for j in range(n) if j != i and tau[i][j] is not None),
                Fraction(0)) / (n - 1) for i in range(n)]


# ---- Slicing and aggregation -------------------------------------------------

def active_at(s, e, t):
    return (s == e == t) or (s <= t < e)


def active_in(s, e, a, b):
    if s == e:
        return a <= s < b
    return s < b and e > a


def grid(md: Model):
    if not md.start:
        return [Fraction(0)]
    lo, hi = min(md.start) - 1, max(md.end) + 1
    return [lo + Fraction(k, 2) for k in range(int(2 * (hi - lo)) + 1)]


def series(md: Model, width, step):
    if not md.start:
        return []
    lo, hi = min(md.start), max(md.end)
    out, k = [], 0
    while lo + k * step <= hi:
        out.append((lo + k * step, lo + k * step + width))
        k += 1
    return out


def aggregate(md: Model, window):
    fx = md.fx
    groups = {}
    for c in range(len(fx.src)):
        s, e = md.start[c], md.end[c]
        if window is not None and not active_in(s, e, *window):
            continue
        u, v = fx.src[c], fx.tgt[c]
        key = (u, v) if fx.directed else (min(u, v), max(u, v))
        if key not in groups:
            groups[key] = {"from": u, "to": v, "count": 0, "duration": Fraction(0),
                           "first": s, "last": e, "order": c}
        g = groups[key]
        g["count"] += 1
        if window is None:
            g["duration"] += e - s
        else:
            g["duration"] += max(Fraction(0), min(e, window[1]) - max(s, window[0]))
        g["first"] = min(g["first"], s)
        g["last"] = max(g["last"], e)
    return sorted(groups.values(), key=lambda g: g["order"])


# ---- Events ------------------------------------------------------------------

def bm(times):
    """Goh and Barabasi's B and M of a sorted event sequence."""
    gaps = [b - a for a, b in zip(times, times[1:])]
    k = len(gaps)

    def mean(x):
        return sum(x, Fraction(0)) / len(x)

    def var(x, mu):
        return sum(((v - mu) ** 2 for v in x), Fraction(0)) / len(x)

    b = math.nan
    if k >= 2:
        mu = mean(gaps)
        sigma = math.sqrt(var(gaps, mu))
        if sigma + mu > 0:
            b = (sigma - float(mu)) / (sigma + float(mu))
    m = math.nan
    if k >= 3:
        x, y = gaps[:-1], gaps[1:]
        m1, m2 = mean(x), mean(y)
        v1, v2 = var(x, m1), var(y, m2)
        if v1 > 0 and v2 > 0:
            cov = sum(((p - m1) * (q - m2) for p, q in zip(x, y)), Fraction(0)) / (k - 1)
            m = float(cov) / math.sqrt(v1 * v2)
    return b, m


def node_events(md: Model, mode: str):
    fx = md.fx
    ev = [[] for _ in range(fx.n)]
    for c in range(len(fx.src)):
        u, v, t = fx.src[c], fx.tgt[c], md.start[c]
        if not fx.directed or mode == "all":
            ev[u].append(t)
            if v != u:
                ev[v].append(t)
        elif mode == "out":
            ev[u].append(t)
        else:
            ev[v].append(t)
    return [sorted(x) for x in ev]


def burstiness(md: Model, scope: str, mode: str):
    if scope == "node":
        rows = []
        for times in node_events(md, mode):
            b, m = bm(times)
            rows.append((len(times), mark(b), mark(m)))
        return rows
    fx = md.fx
    rows = []
    for g in aggregate(md, None):
        key = (g["from"], g["to"])
        times = sorted(md.start[c] for c in range(len(fx.src))
                       if (fx.src[c], fx.tgt[c]) == key
                       or (not fx.directed and (fx.tgt[c], fx.src[c]) == key))
        b, m = bm(times)
        rows.append((g["from"], g["to"], len(times), mark(b), mark(m)))
    return rows


# ---- Dispatch ----------------------------------------------------------------

def compute(fx, e):
    fn = e["fn"]
    md = Model(fx, e)
    if fn == "temporal_paths":
        return [mark(v) for v in matrix(md, e["criterion"])]
    if fn == "temporal_reach":
        return reach_counts(md, e["direction"])
    if fn == "temporal_reachability_ratio":
        if fx.n < 2:
            return [UNDEFINED]
        return [Fraction(sum(reach_counts(md, "out")), fx.n * (fx.n - 1))]
    if fn == "temporal_closeness":
        return [mark(v) for v in closeness(md, False)]
    if fn == "temporal_average_distance":
        return [mark(v) for v in closeness(md, True)]
    if fn == "temporal_active":
        return [tuple(c for c in range(len(fx.src)) if active_at(md.start[c], md.end[c], t))
                for t in grid(md)]
    if fn == "temporal_windows":
        return [tuple(c for c in range(len(fx.src)) if active_in(md.start[c], md.end[c], a, b))
                for a, b in series(md, Fraction(e["width"]), Fraction(e["step"]))]
    if fn == "temporal_aggregate":
        w = e.get("window")
        w = None if w is None else (Fraction(w[0]), Fraction(w[1]))
        return [(g["from"], g["to"], g["count"], g["duration"], g["first"], g["last"])
                for g in aggregate(md, w)]
    if fn == "temporal_burstiness":
        return burstiness(md, e["scope"], e.get("mode", "all"))
    if fn == "temporal_inter_event":
        rows = []
        for v, times in enumerate(node_events(md, e.get("mode", "all"))):
            rows.extend((v, b - a) for a, b in zip(times, times[1:]))
        return rows
    raise ValueError(f"unknown temporal function {fn!r}")


def main() -> None:
    spec = manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            write_values(ORACLE, name, variant, compute(fx, v["engine"]))
            count += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(
        "definition (direct evaluation, exact arithmetic)\n", newline="")
    print(f"definition (temporal): wrote {count} golden files")


if __name__ == "__main__":
    main()
