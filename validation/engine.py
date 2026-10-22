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
                                    f64p, ctypes.c_char_p, ctypes.c_size_t]
    lib.net_betweenness.restype = ctypes.c_int32
    lib.net_triad_census.argtypes = [handle, ctypes.c_int32, f64p, ctypes.c_char_p,
                                     ctypes.c_size_t]
    lib.net_triad_census.restype = ctypes.c_int32
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

    def betweenness(self, weighted: bool = False, normalized: bool = False, threads: int = 0):
        return self._call(self.lib.net_betweenness, self.n, int(weighted), int(normalized),
                          threads)

    def triad_census(self, threads: int = 0):
        return self._call(self.lib.net_triad_census, 16, threads)


def run_variant(g: Graph, spec: dict) -> list[float]:
    fn = spec["fn"]
    if fn == "degree":
        return g.degree(spec["mode"], spec["loops"], weighted=False)
    if fn == "strength":
        return g.degree(spec["mode"], spec["loops"], weighted=True)
    if fn == "betweenness":
        return g.betweenness(spec["weighted"], spec["normalized"])
    if fn == "triad_census":
        return g.triad_census()
    raise ValueError(f"no engine binding for {fn!r}")
