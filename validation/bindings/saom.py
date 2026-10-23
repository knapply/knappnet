"""Engine runners for stochastic actor-oriented models
(validation/manifest.d/saom.yaml).

Every directed fixture becomes the three-wave panel of
validation/saom_fixtures.py. The runners return the deterministic parts
of the model: target statistics per period, change statistics for every
actor and alter on wave 1, and the number of changed ties per period.

Engine options:
  fn          saom_targets | saom_contributions | saom_distances
  effect      an RSiena short name (density, recip, ..., X)
  covariate   constant | changing | binary | dyadic (covariate effects)
  parameter   internal effect parameter (GWESP: 100 alpha; default 69)
"""

from __future__ import annotations

import ctypes
import sys
from pathlib import Path

from engine import ERR_LEN, EngineError, load, register

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import saom_fixtures as sf  # noqa: E402

EFFECT_CODES = {
    "density": 0, "recip": 1, "transTrip": 2, "transTies": 3, "transRecTrip": 4, "cycle3": 5,
    "gwespFF": 6, "gwespFB": 7, "gwespBF": 8, "gwespBB": 9, "inPop": 10, "inPopSqrt": 11,
    "outPop": 12, "outPopSqrt": 13, "inAct": 14, "inActSqrt": 15, "outAct": 16,
    "outActSqrt": 17, "egoX": 18, "altX": 19, "simX": 20, "sameX": 21, "egoXaltX": 22, "X": 23,
}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    i32p, f64p = ctypes.POINTER(i32), ctypes.POINTER(f64)
    h = ctypes.c_void_p
    hp = ctypes.POINTER(h)
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    spec = [i64, i32p, i32p, f64p]
    lib.net_saom_new.argtypes = [i64, i64, f64p, hp, cp, sz]
    lib.net_saom_free.argtypes = [h]
    lib.net_saom_free.restype = None
    lib.net_saom_add_actor_covariate.argtypes = [h, f64p, i64, cp, sz]
    lib.net_saom_add_dyadic_covariate.argtypes = [h, f64p, cp, sz]
    lib.net_saom_period_count.argtypes = [h]
    lib.net_saom_period_count.restype = i64
    lib.net_saom_distances.argtypes = [h, f64p, cp, sz]
    lib.net_saom_targets.argtypes = [h] + spec + [f64p, cp, sz]
    lib.net_saom_contributions.argtypes = [h] + spec + [i64, i64, i64, f64p, cp, sz]
    lib.net_saom_initial.argtypes = [h] + spec + [i32, f64p, cp, sz]
    lib.net_saom_simulate.argtypes = [h] + spec + [f64p, f64p, i32, u64, u64, f64p, f64p, f64p,
                                                  cp, sz]
    lib.net_saom_estimate.argtypes = [h] + spec + [f64p, f64p, u64, i32, hp, cp, sz]
    lib.net_saom_fit_length.argtypes = [h, i32]
    lib.net_saom_fit_length.restype = i64
    lib.net_saom_fit_copy.argtypes = [h, i32, f64p, cp, sz]
    lib.net_saom_fit_free.argtypes = [h]
    lib.net_saom_fit_free.restype = None
    for fn in ("net_saom_new", "net_saom_add_actor_covariate", "net_saom_add_dyadic_covariate",
               "net_saom_distances", "net_saom_targets", "net_saom_contributions",
               "net_saom_initial", "net_saom_simulate", "net_saom_estimate",
               "net_saom_fit_copy"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _doubles(values) -> ctypes.Array:
    return (ctypes.c_double * max(len(values), 1))(*values)


class SaomPanel:
    """The engine's panel for one fixture, with the four covariates."""

    def __init__(self, fx):
        self.lib = _declare()
        self.n = n = fx.n
        waves = sf.waves(fx)
        # Column-major per wave: entry (i, j) at j * n + i.
        flat = []
        for w in waves:
            flat.extend(w[i][j] for j in range(n) for i in range(n))
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_saom_new(n, len(waves), _doubles(flat), ctypes.byref(self.handle),
                                     err, ERR_LEN), err)
        self.actor_index = {}
        self.dyadic_index = {}
        for name, (kind, make) in sf.COVARIATES.items():
            values = make(n)
            if kind == "actor":
                if name == "changing":
                    flat = [v for wave in values for v in wave]
                else:
                    flat = list(values)
                _check(self.lib.net_saom_add_actor_covariate(self.handle, _doubles(flat),
                                                             len(flat), err, ERR_LEN), err)
                self.actor_index[name] = len(self.actor_index)
            else:
                flat = [values[i][j] for j in range(n) for i in range(n)]
                _check(self.lib.net_saom_add_dyadic_covariate(self.handle, _doubles(flat), err,
                                                              ERR_LEN), err)
                self.dyadic_index[name] = len(self.dyadic_index)
        self.periods = self.lib.net_saom_period_count(self.handle)

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_saom_free(self.handle)
            self.handle = None

    def spec(self, e: dict):
        name = e["effect"]
        code = EFFECT_CODES[name]
        cov = -1
        if name in ("egoX", "altX", "simX", "sameX", "egoXaltX"):
            cov = self.actor_index[e["covariate"]]
        elif name == "X":
            cov = self.dyadic_index[e["covariate"]]
        par = float(e.get("parameter", 69))
        return (1, (ctypes.c_int32 * 1)(code), (ctypes.c_int32 * 1)(cov),
                (ctypes.c_double * 1)(par))

    def targets(self, e: dict) -> list[float]:
        out = (ctypes.c_double * max(self.periods, 1))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_saom_targets(self.handle, *self.spec(e), out, err, ERR_LEN), err)
        return list(out)[: self.periods]

    def contributions(self, e: dict) -> list[float]:
        """Change statistics for the first CONTRIBUTION_ACTORS actors i and
        every alter j on wave 1, in i-major order (0 for j = i)."""
        n = self.n
        values = []
        err = ctypes.create_string_buffer(ERR_LEN)
        for i in range(min(n, sf.CONTRIBUTION_ACTORS)):
            out = (ctypes.c_double * max(n, 1))()
            _check(self.lib.net_saom_contributions(self.handle, *self.spec(e), 0, 0, i, out, err,
                                                   ERR_LEN), err)
            values.extend(list(out)[:n])
        return values

    def distances(self) -> list[float]:
        out = (ctypes.c_double * max(self.periods, 1))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_saom_distances(self.handle, out, err, ERR_LEN), err)
        return list(out)[: self.periods]


_cache: dict[str, SaomPanel] = {}


def panel_for(fx) -> SaomPanel:
    if fx.name not in _cache:
        _cache.clear()
        _cache[fx.name] = SaomPanel(fx)
    return _cache[fx.name]


@register("saom_targets")
def _targets(g, spec, fx):
    return panel_for(fx).targets(spec)


@register("saom_contributions")
def _contributions(g, spec, fx):
    return panel_for(fx).contributions(spec)


@register("saom_distances")
def _distances(g, spec, fx):
    return panel_for(fx).distances()
