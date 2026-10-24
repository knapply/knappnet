#!/usr/bin/env python3
"""Definition oracle for spanning trees, branchings and arborescences, and
directed acyclic graphs.

Independent evaluations of the definitions, sharing no code with the
engine: spanning-forest weights by Jarnik's (Prim's) algorithm on a dense
matrix of the lightest (or heaviest) edge between each pair; optimum
branchings and arborescences by trying every choice of at most one
incoming arc per node on small networks; generations as longest-path
depths; acyclicity, closures, and reductions from breadth-first
reachability; longest paths by enumerating every path; dominators from
the definition (d dominates v when v is unreachable from the start once d
is removed) and dominance frontiers from theirs (Cytron et al. 1991);
lowest common ancestors from ancestor sets. Derived networks and
encodings follow validation/bindings/trees.py.

    python3 validation/oracles/trees_definition_oracle.py
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MISSING = -1.0
BRANCHING_MAX_CHOICES = 300_000
PATHS_MAX_NODES = 20


def weights_of(fx, weighted):
    return list(fx.weights) if weighted and fx.weighted else [1.0] * len(fx.src)


def forest_weight(fx, w, maximum):
    n = fx.n
    sign = -1.0 if maximum else 1.0
    best = {}
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key not in best or sign * w[e] < sign * best[key]:
            best[key] = w[e]
    inside = [False] * n
    total = 0.0
    for root in range(n):
        if inside[root]:
            continue
        inside[root] = True
        dist = {}
        for (a, b), x in best.items():
            if root in (a, b):
                other = b if a == root else a
                dist[other] = min(dist.get(other, math.inf), sign * x)
        while True:
            cand = [(d, v) for v, d in dist.items() if not inside[v]]
            if not cand:
                break
            d, v = min(cand)
            inside[v] = True
            total += sign * d
            for (a, b), x in best.items():
                if v in (a, b):
                    other = b if a == v else a
                    if not inside[other]:
                        dist[other] = min(dist.get(other, math.inf), sign * x)
    return total


def branching_weight(fx, w, kind):
    n = fx.n
    arcs_in = [[None] for _ in range(n)]
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a != b:
            arcs_in[b].append(e)
    choices = 1
    for c in arcs_in:
        choices *= len(c)
    if choices > BRANCHING_MAX_CHOICES:
        return None
    if kind == "minimal" and w:
        hi, lo = max(w), min(w)
        score = [2 * hi - lo + 1 - x for x in w]
    else:
        score = w
    best = None
    for pick in itertools.product(*arcs_in):
        roots = sum(1 for p in pick if p is None)
        if kind in ("max_arborescence", "min_arborescence") and roots != 1:
            continue
        parent = [None if p is None else fx.src[p] for p in pick]
        ok = True
        for v in range(n):
            seen, x = set(), v
            while x is not None:
                if x in seen:
                    ok = False
                    break
                seen.add(x)
                x = parent[x]
            if not ok:
                break
        if not ok:
            continue
        s = sum(score[p] for p in pick if p is not None)
        real = sum(w[p] for p in pick if p is not None)
        if kind in ("min_branching", "min_arborescence"):
            better = best is None or s < best[0]
        else:
            better = best is None or s > best[0]
        if better:
            best = (s, real)
    if best is None:
        return MISSING
    return best[1]


def dag_succ(fx):
    succ = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a < b:
            succ[a].add(b)
    return succ


def reach_from(succ, s, removed=None):
    seen = set()
    q = deque([s])
    while q:
        v = q.popleft()
        for w in succ[v]:
            if w != removed and w not in seen:
                seen.add(w)
                q.append(w)
    return seen


def simple_succ(fx):
    succ = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        succ[a].add(b)
        if not fx.directed:
            succ[b].add(a)
    return succ


def dominators(fx):
    n = fx.n
    succ = [set() for _ in range(n)]
    pred = [set() for _ in range(n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            succ[a].add(b)
            pred[b].add(a)
    reachable = {0} | reach_from(succ, 0)
    dom = {}
    for v in reachable:
        dom[v] = {0, v}
        for d in reachable:
            if d not in (0, v) and v not in reach_from_start(succ, d):
                dom[v].add(d)
    idom = {}
    for v in reachable:
        if v == 0:
            idom[v] = 0
            continue
        strict = dom[v] - {v}
        idom[v] = next(d for d in strict if all(o in dom[d] for o in strict))
    return reachable, dom, idom, pred


def reach_from_start(succ, removed):
    seen = {0}
    q = deque([0])
    while q:
        v = q.popleft()
        for w in succ[v]:
            if w != removed and w not in seen:
                seen.add(w)
                q.append(w)
    return seen


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    if fn == "spanning_weight":
        return [forest_weight(fx, weights_of(fx, e.get("weighted")), e.get("maximum", False))]
    if fn == "branching":
        v = branching_weight(fx, weights_of(fx, e.get("weighted")), e["kind"])
        return None if v is None else [v]
    if fn == "generations":
        succ = dag_succ(fx)
        depth = [0] * n
        for v in range(n):  # arcs go from lower to higher index
            for w in succ[v]:
                depth[w] = max(depth[w], depth[v] + 1)
        return depth
    if fn == "is_dag":
        succ = simple_succ(fx)
        cyclic = any(v in reach_from(succ, v) for v in range(n))
        return [0.0 if cyclic or not fx.directed else 1.0]
    if fn == "closure":
        succ = simple_succ(fx)
        loops = {a for a, b in zip(fx.src, fx.tgt) if a == b}
        pairs = []
        for u in range(n):
            r = reach_from(succ, u)
            for v in range(n):
                if u != v and v in r:
                    pairs.append((u, v))
                elif u == v:
                    mode = e["reflexive"]
                    if (mode == "always" or (mode == "cycles" and u in r)
                            or (mode == "never" and u in loops)):
                        pairs.append((u, v))
        return pairs
    if fn == "reduction":
        succ = dag_succ(fx)
        pairs = []
        for u in range(n):
            longer = set()
            for w in succ[u]:
                longer |= reach_from(succ, w)
            pairs += [(u, v) for v in sorted(succ[u]) if v not in longer]
        return pairs
    if fn == "longest_length":
        if n > PATHS_MAX_NODES:
            return None
        weighted = e.get("weighted", False) and fx.weighted
        best_arc = {}
        for k, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            if a < b:
                x = fx.weights[k] if weighted else 1.0
                best_arc[(a, b)] = max(best_arc.get((a, b), -math.inf), x)
        succ = dag_succ(fx)
        best = 0.0

        def walk(v, total):
            nonlocal best
            best = max(best, total)
            for w in succ[v]:
                walk(w, total + best_arc[(v, w)])

        for v in range(n):
            walk(v, 0.0)
        return [best]
    if fn == "idom":
        if n == 0:
            return []
        reachable, _, idom, _ = dominators(fx)
        return [float(idom[v]) if v in reachable else MISSING for v in range(n)]
    if fn == "frontiers":
        if n == 0:
            return []
        reachable, dom, _, pred = dominators(fx)
        loops = {a for a, b in zip(fx.src, fx.tgt) if a == b}
        out = []
        for x in range(n):
            if x not in reachable:
                out.append((MISSING,))
                continue
            f = set()
            for y in reachable:
                preds = {p for p in pred[y] if p in reachable} | ({y} if y in loops else set())
                if any(x in dom[p] for p in preds) and not (x in dom[y] and x != y):
                    f.add(y)
            out.append(tuple(sorted(f)))
        return out
    if fn == "dag_lca_check":
        succ = dag_succ(fx)
        anc = [{v} for v in range(n)]
        for v in range(n):
            for w in succ[v]:
                anc[w] |= anc[v]
        pairs = [(u, v) for u in range(n) for v in range(u, n)]
        return [float(sum(1 for u, v in pairs if anc[u] & anc[v])), 1.0]
    if fn == "tree_lca":
        adj = [set() for _ in range(n)]
        for a, b in zip(fx.src, fx.tgt):
            if a != b:
                adj[a].add(b)
                adj[b].add(a)
        parent = {0: None} if n else {}
        q = deque([0] if n else [])
        while q:
            v = q.popleft()
            for w in sorted(adj[v]):
                if w not in parent:
                    parent[w] = v
                    q.append(w)
        nodes = sorted(parent)

        def chain(v):
            out = []
            while v is not None:
                out.append(v)
                v = parent[v]
            return out[::-1]

        res = []
        for i, u in enumerate(nodes):
            for v in nodes[i:]:
                cu, cv = chain(u), chain(v)
                k = 0
                while k < min(len(cu), len(cv)) and cu[k] == cv[k]:
                    k += 1
                res.append(float(cu[k - 1]))
        return res
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("trees")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", []):
                continue
            values = compute(spec, fx)
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (trees): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
