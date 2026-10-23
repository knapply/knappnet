#!/usr/bin/env python3
"""Golden values for the stochastic actor-oriented model, computed
directly from the definitions of the effects.

Every directed fixture is read as the three-wave panel of
validation/saom_fixtures.py. Each effect's actor statistic s_i(x) is
written out as the formula in chapter 13 of the RSiena manual (Ripley et
al. 2026), with plain loops over actors and no incremental bookkeeping:
transitive triplets as sum_{j,h} x_ij x_ih x_hj, GWESP as the sum over
alters of e^a (1 - (1 - e^-a)^esp_ij) with esp_ij counted from the
definition of each two-path type, and so on. Covariates are centered as
the manual describes (section 5.3.2): actor covariates by their mean over
all actors and waves, dyadic covariates by the mean of the off-diagonal
cells, similarity 1 - |v_i - v_j| / range by its mean over ordered pairs of
distinct actors and waves; `sameX` is not centered. Means are exact
rationals, converted to floating point at the end.

Target statistics for period m are sum_i s_i of wave m + 1 with the
covariates of wave m. Change statistics are s_i(x with i -> j) - s_i(x
without i -> j), evaluated by recomputing s_i on both networks (for the
first CONTRIBUTION_ACTORS actors i and every alter j), except for the GWESP
effects, which RSiena defines as elementary effects whose change statistic
is the weight of the tie i -> j itself.

    python3 validation/oracles/saom_definition_oracle.py
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, MANIFEST_DIR, applies, fixture_names, load_fixture  # noqa: E402
from common import write_values  # noqa: E402
import saom_fixtures as sf  # noqa: E402

ORACLE = "definition"


def manifest() -> dict:
    return yaml.safe_load((MANIFEST_DIR / "saom.yaml").read_text())


# ---- Covariates ---------------------------------------------------------------

class Covariates:
    def __init__(self, n: int):
        self.n = n
        self.actor = {}
        for name, (kind, make) in sf.COVARIATES.items():
            values = make(n)
            if kind == "actor":
                # One column per period: the last wave's values of a
                # changing covariate are not used (as in RSiena).
                waves = values[:sf.WAVES - 1] if name == "changing" else [values] * (sf.WAVES - 1)
                waves = [[Fraction(v) for v in w] for w in waves]
                self.actor[name] = waves
            else:
                self.dyadic = [[Fraction(v) for v in row] for row in values]
        cells = [self.dyadic[i][j] for i in range(n) for j in range(n) if i != j]
        self.dyadic_mean = sum(cells) / len(cells) if cells else Fraction(0)
        # Means, ranges, and similarity means, once per covariate.
        self.means, self.ranges, self.sim_means = {}, {}, {}
        for name, waves in self.actor.items():
            flat = [v for w in waves for v in w]
            self.means[name] = sum(flat) / len(flat)
            r = max(flat) - min(flat)
            self.ranges[name] = r
            if r > 0 and n > 1:
                sims = [1 - abs(w[i] - w[j]) / r for w in waves
                        for i in range(n) for j in range(n) if i != j]
                self.sim_means[name] = sum(sims) / len(sims)

    def centered(self, name: str, period: int) -> list[Fraction]:
        m = self.means[name]
        return [v - m for v in self.actor[name][period]]

    def raw(self, name: str, period: int) -> list[Fraction]:
        return self.actor[name][period]

    def similarity(self, name: str, period: int, i: int, j: int) -> Fraction:
        v = self.actor[name][period]
        return 1 - abs(v[i] - v[j]) / self.ranges[name] - self.sim_means[name]


# ---- Effects, from their formulas --------------------------------------------

def gwesp_weight(alpha: float, esp: int) -> float:
    if esp == 0:
        return 0.0
    return math.exp(alpha) * (1 - (1 - math.exp(-alpha)) ** esp)


def esp(x, i: int, j: int, kind: str) -> int:
    """Shared partners h of the tie i -> j for the four GWESP types."""
    n = len(x)
    count = 0
    for h in range(n):
        if h in (i, j):
            continue
        if kind == "gwespFF":
            count += x[i][h] * x[h][j]
        elif kind == "gwespFB":
            count += x[i][h] * x[j][h]
        elif kind == "gwespBF":
            count += x[h][i] * x[h][j]
        else:  # gwespBB: two-path i <- h <- j
            count += x[h][i] * x[j][h]
    return count


def actor_statistic(e: dict, x, cov: Covariates, period: int, i: int):
    n = len(x)
    name = e["effect"]
    alters = [j for j in range(n) if x[i][j]]
    out_i = len(alters)
    if name == "density":
        return out_i
    if name == "recip":
        return sum(x[j][i] for j in alters)
    if name == "transTrip":
        return sum(x[i][j] * x[i][h] * x[h][j] for j in range(n) for h in range(n))
    if name == "transTies":
        return sum(max(x[i][h] * x[h][j] for h in range(n)) for j in alters)
    if name == "transRecTrip":
        return sum(x[i][j] * x[j][i] * x[i][h] * x[h][j] for j in range(n) for h in range(n))
    if name == "cycle3":
        return sum(x[i][j] * x[j][h] * x[h][i] for j in range(n) for h in range(n))
    if name.startswith("gwesp"):
        alpha = float(e.get("parameter", 69)) / 100
        return sum(gwesp_weight(alpha, esp(x, i, j, name)) for j in alters)
    indeg = lambda j: sum(x[h][j] for h in range(n))  # noqa: E731
    outdeg = lambda j: sum(x[j])  # noqa: E731
    if name == "inPop":
        return sum(indeg(j) for j in alters)
    if name == "inPopSqrt":
        return sum(math.sqrt(indeg(j)) for j in alters)
    if name == "outPop":
        return sum(outdeg(j) for j in alters)
    if name == "outPopSqrt":
        return sum(math.sqrt(outdeg(j)) for j in alters)
    if name == "inAct":
        return out_i * indeg(i)
    if name == "inActSqrt":
        return out_i * math.sqrt(indeg(i))
    if name == "outAct":
        return out_i ** 2
    if name == "outActSqrt":
        return out_i ** 1.5
    c = e.get("covariate")
    if name == "egoX":
        return cov.centered(c, period)[i] * out_i
    if name == "altX":
        v = cov.centered(c, period)
        return sum(v[j] for j in alters)
    if name == "simX":
        return sum(cov.similarity(c, period, i, j) for j in alters)
    if name == "sameX":
        v = cov.raw(c, period)
        return sum(1 for j in alters if v[i] == v[j])
    if name == "egoXaltX":
        v = cov.centered(c, period)
        return v[i] * sum(v[j] for j in alters)
    if name == "X":
        return sum(cov.dyadic[i][j] - cov.dyadic_mean for j in alters)
    raise ValueError(name)


def statistic(e: dict, x, cov: Covariates, period: int) -> float:
    total = sum(actor_statistic(e, x, cov, period, i) for i in range(len(x)))
    # RSiena reports the number of 3-cycles; the actor sum counts each three times.
    if e["effect"] == "cycle3":
        total = Fraction(total) / 3
    return float(total)


def change_statistics(e: dict, x, cov: Covariates, period: int) -> list[float]:
    n = len(x)
    out = []
    for i in range(min(n, sf.CONTRIBUTION_ACTORS)):
        for j in range(n):
            if i == j:
                out.append(0.0)
                continue
            if e["effect"].startswith("gwesp"):
                alpha = float(e.get("parameter", 69)) / 100
                out.append(gwesp_weight(alpha, esp(x, i, j, e["effect"])))
                continue
            with_tie = [row[:] for row in x]
            with_tie[i][j] = 1
            without = [row[:] for row in x]
            without[i][j] = 0
            out.append(float(actor_statistic(e, with_tie, cov, period, i)
                             - actor_statistic(e, without, cov, period, i)))
    return out


def distances(waves) -> list[float]:
    n = len(waves[0])
    return [float(sum(1 for i in range(n) for j in range(n) if waves[m][i][j] != waves[m + 1][i][j]))
            for m in range(len(waves) - 1)]


def main() -> None:
    spec = manifest()
    for name in fixture_names():
        fx = load_fixture(name)
        if not fx.directed:
            continue
        waves = sf.waves(fx)
        cov = Covariates(fx.n)
        for variant, v in spec.items():
            if not applies(v, fx) or "definition" not in v["oracles"]:
                continue
            e = v["engine"]
            if e["fn"] == "saom_distances":
                values = distances(waves)
            elif e["fn"] == "saom_targets":
                values = [statistic(e, waves[m + 1], cov, m) for m in range(len(waves) - 1)]
            elif e["fn"] == "saom_contributions":
                values = change_statistics(e, waves[0], cov, 0)
            else:
                raise ValueError(e["fn"])
            write_values(ORACLE, name, variant, values)
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    print("definition oracle (saom): done")


if __name__ == "__main__":
    main()
