"""Network panels for the SAOM validation.

Stochastic actor-oriented models need two or more observations of the
same actors. Every directed fixture of the ordinary corpus
(validation/fixtures/) is turned into a three-wave panel: wave 1 is the
fixture read as an adjacency matrix (self-loops dropped, parallel edges
collapsed), and each later wave toggles a fixed, pseudo-random set of tie
variables of the wave before it. Four covariates are derived from the
actor indices: a constant one with values 0..4, a changing one (three
waves, values 0..3), a binary one, and a dyadic one with values 0, 0.5,
..., 3. Everything is a small integer or half-integer, so every oracle sees
exactly the same data.

Both the engine runners (validation/bindings/saom.py) and the oracles
(validation/oracles/saom_definition_oracle.py, validation/oracles/r_saom.R)
build the panel from here or from the identical rules written in R.
"""

from __future__ import annotations

WAVES = 3

# Change statistics are validated for the first actors of each fixture (and
# every alter), which keeps the brute-force oracle fast on the larger ones.
CONTRIBUTION_ACTORS = 12


def adjacency(fx) -> list[list[int]]:
    """Wave 1: the fixture as a 0/1 matrix without loops."""
    n = fx.n
    a = [[0] * n for _ in range(n)]
    for i, j in zip(fx.src, fx.tgt):
        if i != j:
            a[i][j] = 1
    return a


def toggle(a: list[list[int]], mult: int, add: int, mod: int) -> list[list[int]]:
    """A copy of `a` in which tie (i, j) is flipped when
    (mult * (i n + j) + add) mod `mod` == 0."""
    n = len(a)
    b = [row[:] for row in a]
    for i in range(n):
        for j in range(n):
            if i != j and (mult * (i * n + j) + add) % mod == 0:
                b[i][j] ^= 1
    return b


def waves(fx) -> list[list[list[int]]]:
    w1 = adjacency(fx)
    w2 = toggle(w1, 37, 11, 7)
    w3 = toggle(w2, 53, 5, 9)
    return [w1, w2, w3]


def constant_covariate(n: int) -> list[float]:
    return [float((i * 7) % 5) for i in range(n)]


def changing_covariate(n: int) -> list[list[float]]:
    """One list per wave."""
    return [[float((i * 3 + 2 * w) % 4) for i in range(n)] for w in range(WAVES)]


def binary_covariate(n: int) -> list[float]:
    return [float(i % 2) for i in range(n)]


def dyadic_covariate(n: int) -> list[list[float]]:
    return [[((i * 5 + j * 3) % 7) / 2 for j in range(n)] for i in range(n)]


# Covariate names used in the manifest, with the actor/dyadic kind.
COVARIATES = {
    "constant": ("actor", constant_covariate),
    "changing": ("actor", changing_covariate),
    "binary": ("actor", binary_covariate),
    "dyadic": ("dyadic", dyadic_covariate),
}
