"""ctypes binding to the engine's C ABI.

Validation goes through the same C interface the R package uses, so a
passing comparison here also exercises the exported ABI.
"""

from __future__ import annotations

import ctypes
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ERR_LEN = 1024

_lib = None


def _library_path() -> Path:
    if env := os.environ.get("ENGINE_LIB"):
        return Path(env)
    name = {"darwin": "libcapi.dylib", "win32": "capi.dll"}.get(sys.platform, "libcapi.so")
    return ROOT / "target" / "release" / name


def load() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    path = _library_path()
    if "ENGINE_LIB" not in os.environ:
        subprocess.run(["cargo", "build", "--release", "--offline", "-p", "capi"],
                       cwd=ROOT, check=True, capture_output=True)
    lib = ctypes.CDLL(str(path))
    i64p = ctypes.POINTER(ctypes.c_int64)
    f64p = ctypes.POINTER(ctypes.c_double)
    handle = ctypes.c_void_p
    lib.net_graph_new.argtypes = [ctypes.c_int64, i64p, i64p, f64p, ctypes.c_int64,
                                  ctypes.c_int32, ctypes.POINTER(handle),
                                  ctypes.c_char_p, ctypes.c_size_t]
    lib.net_graph_new.restype = ctypes.c_int32
    lib.net_graph_free.argtypes = [handle]
    lib.net_degree.argtypes = [handle, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
                               f64p, ctypes.c_char_p, ctypes.c_size_t]
    lib.net_degree.restype = ctypes.c_int32
    lib.net_betweenness.argtypes = [handle, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
                                    ctypes.c_int32, f64p, ctypes.c_char_p, ctypes.c_size_t]
    lib.net_betweenness.restype = ctypes.c_int32
    lib.net_triad_census.argtypes = [handle, ctypes.c_int32, f64p, ctypes.c_char_p,
                                     ctypes.c_size_t]
    lib.net_triad_census.restype = ctypes.c_int32
    cp = ctypes.c_char_p
    lib.net_read.argtypes = [cp, cp, ctypes.c_int32, ctypes.POINTER(handle), cp, ctypes.c_size_t]
    lib.net_read.restype = ctypes.c_int32
    lib.net_network_new.argtypes = [ctypes.c_int32]
    lib.net_network_new.restype = handle
    lib.net_network_free.argtypes = [handle]
    lib.net_network_directed.argtypes = [handle]
    lib.net_network_directed.restype = ctypes.c_int32
    lib.net_table_nrow.argtypes = [handle, ctypes.c_int32]
    lib.net_table_nrow.restype = ctypes.c_int64
    lib.net_table_ncol.argtypes = [handle, ctypes.c_int32]
    lib.net_table_ncol.restype = ctypes.c_int32
    lib.net_table_name.argtypes = [handle, ctypes.c_int32, ctypes.c_int32]
    lib.net_table_name.restype = cp
    lib.net_table_kind.argtypes = [handle, ctypes.c_int32, ctypes.c_int32]
    lib.net_table_kind.restype = ctypes.c_int32
    lib.net_table_numeric.argtypes = [handle, ctypes.c_int32, ctypes.c_int32, f64p,
                                      ctypes.POINTER(ctypes.c_int32), cp, ctypes.c_size_t]
    lib.net_table_numeric.restype = ctypes.c_int32
    lib.net_table_string.argtypes = [handle, ctypes.c_int32, ctypes.c_int32, ctypes.c_int64]
    lib.net_table_string.restype = cp
    lib.net_table_add_numeric.argtypes = [handle, ctypes.c_int32, cp, ctypes.c_int32, f64p,
                                          ctypes.POINTER(ctypes.c_int32), ctypes.c_int64,
                                          cp, ctypes.c_size_t]
    lib.net_table_add_numeric.restype = ctypes.c_int32
    lib.net_table_add_string.argtypes = [handle, ctypes.c_int32, cp, ctypes.POINTER(cp),
                                         ctypes.c_int64, cp, ctypes.c_size_t]
    lib.net_table_add_string.restype = ctypes.c_int32
    lib.net_write.argtypes = [handle, cp, cp, cp, cp, ctypes.c_size_t]
    lib.net_write.restype = ctypes.c_int32
    _lib = lib
    return lib


class EngineError(RuntimeError):
    pass


class Graph:
    def __init__(self, n: int, src: list[int], tgt: list[int],
                 weights: list[float] | None, directed: bool):
        self.lib = load()
        self.n = n
        m = len(src)
        s = (ctypes.c_int64 * max(m, 1))(*src)
        t = (ctypes.c_int64 * max(m, 1))(*tgt)
        w = (ctypes.c_double * max(m, 1))(*weights) if weights is not None else None
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        status = self.lib.net_graph_new(n, s, t, w, m, int(directed),
                                        ctypes.byref(self.handle), err, ERR_LEN)
        if status:
            raise EngineError(err.value.decode())

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_graph_free(self.handle)
            self.handle = None

    def _call(self, fn, size: int, *args) -> list[float]:
        out = (ctypes.c_double * max(size, 1))()
        err = ctypes.create_string_buffer(ERR_LEN)
        if fn(self.handle, *args, out, err, ERR_LEN):
            raise EngineError(err.value.decode())
        return list(out)[:size]

    def degree(self, mode: str = "all", loops: bool = True, weighted: bool = False):
        code = {"out": 1, "in": 2, "all": 3}[mode]
        return self._call(self.lib.net_degree, self.n, code, int(loops), int(weighted))

    def betweenness(self, weighted: bool = False, normalized: bool = False, threads: int = 0,
                    weights_as: str = "distance"):
        kind = {"distance": 0, "strength": 1}[weights_as]
        return self._call(self.lib.net_betweenness, self.n, int(weighted), kind,
                          int(normalized), threads)

    def triad_census(self, threads: int = 0):
        return self._call(self.lib.net_triad_census, 16, threads)


def run_variant(g: Graph, spec: dict) -> list[float]:
    fn = spec["fn"]
    if fn == "degree":
        return g.degree(spec["mode"], spec["loops"], weighted=False)
    if fn == "strength":
        return g.degree(spec["mode"], spec["loops"], weighted=True)
    if fn == "betweenness":
        return g.betweenness(spec["weighted"], spec["normalized"],
                             weights_as=spec.get("weights_as", "distance"))
    if fn == "triad_census":
        return g.triad_census()
    raise ValueError(f"no engine binding for {fn!r}")


# ---- Files ----------------------------------------------------------------
# A network is {"directed": bool, "graph": {name: value},
#               "nodes": {column: [values]}, "edges": {column: [values]}},
# with None for missing values. Column order is preserved.

KINDS = {0: bool, 1: int, 2: float, 3: str}


def read_network(path, fmt=None, directed=None) -> dict:
    lib = load()
    h = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    d = -1 if directed is None else int(directed)
    status = lib.net_read(str(path).encode(), fmt.encode() if fmt else None, d,
                          ctypes.byref(h), err, ERR_LEN)
    if status:
        raise EngineError(err.value.decode())
    try:
        tables = []
        for which in range(3):
            nrow = lib.net_table_nrow(h, which)
            cols = {}
            for c in range(lib.net_table_ncol(h, which)):
                name = lib.net_table_name(h, which, c).decode()
                kind = lib.net_table_kind(h, which, c)
                if kind == 3:
                    vals = []
                    for r in range(nrow):
                        p = lib.net_table_string(h, which, c, r)
                        vals.append(None if p is None else p.decode())
                else:
                    out = (ctypes.c_double * max(nrow, 1))()
                    na = (ctypes.c_int32 * max(nrow, 1))()
                    if lib.net_table_numeric(h, which, c, out, na, err, ERR_LEN):
                        raise EngineError(err.value.decode())
                    cast = KINDS[kind]
                    vals = [None if na[r] else cast(out[r]) for r in range(nrow)]
                cols[name] = vals
            tables.append(cols)
        graph = {k: v[0] for k, v in tables[0].items()}
        return {"directed": bool(lib.net_network_directed(h)), "graph": graph,
                "nodes": tables[1], "edges": tables[2]}
    finally:
        lib.net_network_free(h)


def _kind_of(values) -> int:
    present = [v for v in values if v is not None]
    if present and all(isinstance(v, bool) for v in present):
        return 0
    if present and all(isinstance(v, int) and not isinstance(v, bool) for v in present):
        return 1
    if present and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in present):
        return 2
    return 3


def write_network(net: dict, path, fmt=None, weight=None) -> None:
    lib = load()
    h = lib.net_network_new(int(net["directed"]))
    err = ctypes.create_string_buffer(ERR_LEN)
    try:
        graph = {k: [v] for k, v in net.get("graph", {}).items()}
        for which, table in enumerate([graph, net["nodes"], net["edges"]]):
            for name, values in table.items():
                n = len(values)
                kind = 3 if name in ("id", "from", "to") else _kind_of(values)
                if kind == 3:
                    arr = (ctypes.c_char_p * max(n, 1))(
                        *[None if v is None else str(v).encode() for v in values])
                    status = lib.net_table_add_string(h, which, name.encode(), arr, n, err, ERR_LEN)
                else:
                    vals = (ctypes.c_double * max(n, 1))(*[0.0 if v is None else float(v) for v in values])
                    na = (ctypes.c_int32 * max(n, 1))(*[int(v is None) for v in values])
                    status = lib.net_table_add_numeric(h, which, name.encode(), kind, vals, na, n,
                                                       err, ERR_LEN)
                if status:
                    raise EngineError(err.value.decode())
        status = lib.net_write(h, str(path).encode(), fmt.encode() if fmt else None,
                               weight.encode() if weight else None, err, ERR_LEN)
        if status:
            raise EngineError(err.value.decode())
    finally:
        lib.net_network_free(h)
