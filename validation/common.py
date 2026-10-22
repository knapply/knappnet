"""Shared helpers for loading fixtures, golden files, and the manifest."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"
GOLDEN = HERE / "golden"
MANIFEST = HERE / "manifest.yaml"
MANIFEST_DIR = HERE / "manifest.d"


@dataclass
class Fixture:
    name: str
    n: int
    directed: bool
    weighted: bool
    tags: list[str]
    src: list[int]
    tgt: list[int]
    weights: list[float] | None


def load_fixture(name: str) -> Fixture:
    d = FIXTURES / name
    meta = json.loads((d / "meta.json").read_text())
    src, tgt, w = [], [], []
    with open(d / "edges.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            src.append(int(row["from"]))
            tgt.append(int(row["to"]))
            if meta["weighted"]:
                w.append(float(row["weight"]))
    return Fixture(name, meta["n"], meta["directed"], meta["weighted"], meta["tags"],
                   src, tgt, w if meta["weighted"] else None)


def fixture_names() -> list[str]:
    return sorted(p.name for p in FIXTURES.iterdir() if (p / "meta.json").exists())


def load_manifest() -> dict:
    """The main manifest plus every file in manifest.d/, in name order."""
    manifest = yaml.safe_load(MANIFEST.read_text()) or {}
    for extra in sorted(MANIFEST_DIR.glob("*.yaml")):
        for key, spec in (yaml.safe_load(extra.read_text()) or {}).items():
            if key in manifest:
                raise ValueError(f"measure {key!r} is defined twice (in {extra.name})")
            manifest[key] = spec
    return manifest


def applies(variant: dict, fx: Fixture) -> bool:
    kind = variant.get("applies_to", "any")
    if kind == "directed" and not fx.directed:
        return False
    if kind == "undirected" and fx.directed:
        return False
    if variant.get("requires_weights") and not fx.weighted:
        return False
    return True


def write_values(oracle: str, fixture: str, variant: str, values) -> None:
    """Writes one value per line. A value may be a number or a sequence of
    numbers (written space-separated, e.g. one clique per line)."""
    d = GOLDEN / oracle / fixture
    d.mkdir(parents=True, exist_ok=True)
    with open(d / f"{variant}.csv", "w", newline="") as fh:
        fh.write("value\n")
        for v in values:
            if isinstance(v, (list, tuple)):
                fh.write(" ".join(f"{float(x):.17g}" for x in v) + "\n")
            else:
                fh.write(f"{float(v):.17g}\n")


def read_values(oracle: str, fixture: str, variant: str):
    """Numbers, or tuples of numbers for lines holding several values."""
    p = GOLDEN / oracle / fixture / f"{variant}.csv"
    if not p.exists():
        return None
    out = []
    with open(p) as fh:
        next(fh)
        for line in fh:
            # R writes missing values as NA; they are undefined, like NaN.
            parts = ["nan" if x == "NA" else x for x in line.split()]
            if not parts:
                # An empty line is an empty group (e.g. no members).
                out.append(())
            elif len(parts) == 1:
                out.append(float(parts[0]))
            else:
                out.append(tuple(float(x) for x in parts))
    return out
