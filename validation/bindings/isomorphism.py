"""Engine runners for VF2++ isomorphism, induced subgraph isomorphism, and
monomorphism, the quick isomorphism checks, and tree isomorphism
(crates/capi/src/isomorphism.rs).

Derived inputs, which the oracles build the same way:
  rev(fx)         the fixture with node v renamed n - 1 - v and its edges
                  in reverse order (isomorphic to the fixture)
  perturbed(fx)   the fixture without its first edge, plus edge (0, n - 1)
  PATTERNS        small graphs searched for in each fixture (below), on
                  nodes 0..k - 1, directed for directed fixtures
  labels          v % 2 on every graph (rev(fx) keeps each node's label)
  bfs_tree(fx)    the breadth-first tree from node 0 over the simple
                  undirected graph (neighbours in increasing order), on the
                  nodes it reaches renamed 0, 1, ... in order of discovery
  dfs_tree(fx)    the same with depth-first search (a different tree)

Encodings: counts are capped at CAP; a search result is (number of
mappings, 1 if the first mapping returned is valid, 0 if there is none);
booleans are 0 or 1; an engine error the manifest declares undefined is
[-1].

Manifest engine options:
  fn: vf2pp           problem (iso | induced | mono), pattern (a key of
                      PATTERNS; iso compares the fixture with rev(fx)),
                      labels (true / false)
  fn: could_be_iso    properties (d | dt | dtc), against (rev | perturbed)
  fn: tree_iso        (free: rev tree, free: dfs tree, rooted with the
                      roots mapped, rooted with a wrong root)
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, Graph, load, register

CAP = 5000
PROBLEM = {"iso": 0, "induced": 1, "mono": 2}
PATTERNS = {
    "path3": (3, [(0, 1), (1, 2)]),
    "triangle": (3, [(0, 1), (1, 2), (2, 0)]),
    "star3": (4, [(0, 1), (0, 2), (0, 3)]),
    "cycle4": (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "feedforward": (3, [(0, 1), (1, 2), (0, 2)]),
}
PROPERTIES = {"d": 1, "dt": 3, "dtc": 7}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64 = ctypes.c_int32, ctypes.c_int64
    i64p, i32p = ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    hp = ctypes.POINTER(h)
    lib.net_vf2pp_new.argtypes = [h, h, i64p, i64p, i32, i32, hp, cp, sz]
    lib.net_vf2pp_next.argtypes = [h, i64, i64p, i64p, cp, sz]
    lib.net_vf2pp_free.argtypes = [h]
    lib.net_vf2pp_free.restype = None
    lib.net_vf2pp_count.argtypes = [h, h, i64p, i64p, i32, i32, i64, i64p, cp, sz]
    lib.net_could_be_isomorphic.argtypes = [h, h, i32, i32, i32p, cp, sz]
    lib.net_tree_isomorphism.argtypes = [h, h, i64, i64, i64p, i32p, cp, sz]
    for fn in ("net_vf2pp_new", "net_vf2pp_next", "net_vf2pp_count",
               "net_could_be_isomorphic", "net_tree_isomorphism"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def rev(fx):
    n = fx.n
    pairs = [(n - 1 - a, n - 1 - b) for a, b in zip(fx.src, fx.tgt)][::-1]
    return n, pairs


def perturbed(fx):
    pairs = list(zip(fx.src, fx.tgt))[1:] + [(0, fx.n - 1)]
    return fx.n, pairs


def _graph(n, pairs, directed):
    return Graph(n, [a for a, _ in pairs], [b for _, b in pairs], None, directed)


def _arr(values):
    return (ctypes.c_int64 * max(len(values), 1))(*values)


def edge_counts(n, pairs, directed):
    c = {}
    for a, b in pairs:
        k = (a, b) if directed else (min(a, b), max(a, b))
        c[k] = c.get(k, 0) + 1
    return c


def valid(mapping, sg, fg, problem, directed, labels):
    """Whether mapping (sg node -> fg node) solves the problem."""
    (ns, ps), (nf, pf) = sg, fg
    if len(set(mapping)) != ns:
        return False
    if labels and any(u % 2 != mapping[u] % 2 for u in range(ns)):
        return False
    cs, cf = edge_counts(ns, ps, directed), edge_counts(nf, pf, directed)
    key = (lambda a, b: (a, b)) if directed else (lambda a, b: (min(a, b), max(a, b)))
    for u in range(ns):
        for w in range(ns):
            if not directed and w < u:
                continue
            s = cs.get(key(u, w), 0)
            f = cf.get(key(mapping[u], mapping[w]), 0)
            if (problem == "mono" and s > f) or (problem != "mono" and s != f):
                return False
    return True


@register("vf2pp")
def _vf2pp(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    problem = spec["problem"]
    fg = (fx.n, list(zip(fx.src, fx.tgt)))
    if problem == "iso":
        sg = rev(fx)
    else:
        sg = PATTERNS[spec["pattern"]]
    hs = _graph(sg[0], sg[1], fx.directed)
    lf = ls = None
    if spec.get("labels"):
        lf = _arr([v % 2 for v in range(fx.n)])
        ls = _arr([v % 2 for v in range(sg[0])])
    m = ctypes.c_void_p()
    _check(lib.net_vf2pp_new(g.handle, hs.handle, lf, ls, PROBLEM[problem], 0,
                             ctypes.byref(m), err, ERR_LEN), err)
    try:
        ns = sg[0]
        batch = 256
        buf = (ctypes.c_int64 * max(batch * ns, 1))()
        found = ctypes.c_int64()
        total, first_ok = 0, 0
        while total < CAP:
            _check(lib.net_vf2pp_next(m, batch, buf, ctypes.byref(found), err, ERR_LEN), err)
            k = found.value
            if total == 0 and k > 0:
                first = [buf[i] for i in range(ns)]
                first_ok = int(valid(first, sg, fg, problem, fx.directed, spec.get("labels")))
            total += k
            if k < batch:
                break
        count = ctypes.c_int64()
        _check(lib.net_vf2pp_count(g.handle, hs.handle, lf, ls, PROBLEM[problem], 0, CAP,
                                   ctypes.byref(count), err, ERR_LEN), err)
        if min(total, CAP) != count.value:
            raise EngineError("iterator and count disagree")
        return [(min(total, CAP), first_ok)]
    finally:
        lib.net_vf2pp_free(m)


@register("could_be_iso")
def _could_be_iso(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    other = rev(fx) if spec["against"] == "rev" else perturbed(fx)
    h2 = _graph(other[0], other[1], fx.directed)
    out = ctypes.c_int32()
    status = lib.net_could_be_isomorphic(g.handle, h2.handle, PROPERTIES[spec["properties"]], 0,
                                         ctypes.byref(out), err, ERR_LEN)
    if status and fx.directed:
        return [-1]
    _check(status, err)
    return [out.value]


def _tree(fx, depth_first):
    nbrs = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    name = {0: 0}
    edges = []
    if depth_first:
        stack = [0]
        seen = {0}
        while stack:
            v = stack[-1]
            nxt = next((w for w in sorted(nbrs[v]) if w not in seen), None)
            if nxt is None:
                stack.pop()
                continue
            seen.add(nxt)
            name[nxt] = len(name)
            edges.append((name[v], name[nxt]))
            stack.append(nxt)
    else:
        queue = [0]
        for v in queue:
            for w in sorted(nbrs[v]):
                if w not in name:
                    name[w] = len(name)
                    edges.append((name[v], name[w]))
                    queue.append(w)
    return len(name), edges


def trees(fx):
    t1 = _tree(fx, False)
    k = t1[0]
    t2 = (k, [(k - 1 - a, k - 1 - b) for a, b in t1[1]][::-1])
    t3 = _tree(fx, True)
    return t1, t2, t3


@register("tree_iso")
def _tree_iso(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    t1, t2, t3 = trees(fx)
    k = t1[0]
    h1, h2, h3 = (_graph(t[0], t[1], False) for t in (t1, t2, t3))
    mapping = (ctypes.c_int64 * k)()
    found = ctypes.c_int32()

    def run(ha, hb, r1, r2, target):
        _check(lib.net_tree_isomorphism(ha.handle, hb.handle, r1, r2, mapping,
                                        ctypes.byref(found), err, ERR_LEN), err)
        if not found.value:
            return 0, 0
        m = [mapping[i] for i in range(k)]
        ok = valid(m, t1, target, "iso", False, False)
        return 1, int(ok)

    free_rev = run(h1, h2, -1, -1, t2)
    free_dfs = run(h1, h3, -1, -1, t3) if t3[0] == k else (0, 0)
    rooted = run(h1, h2, 0, k - 1, t2)
    # Node k - 2 of t2 is the image of node 1 of t1, a child of the root.
    wrong = run(h1, h2, 0, k - 2, t2) if k > 2 else (0, 0)
    return [free_rev, free_dfs, rooted, (wrong[0], wrong[1])]
