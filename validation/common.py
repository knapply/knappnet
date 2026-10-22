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
    return yaml.safe_load(MANIFEST.read_text())


def applies(variant: dict, fx: Fixture) -> bool:
    kind = variant.get("applies_to", "any")
    if kind == "directed" and not fx.directed:
        return False
    if kind == "undirected" and fx.directed:
        return False
    if variant.get("requires_weights") and not fx.weighted:
        return False
    return True


def write_values(oracle: str, fixture: str, variant: str, values: list[float]) -> None:
    d = GOLDEN / oracle / fixture
    d.mkdir(parents=True, exist_ok=True)
    with open(d / f"{variant}.csv", "w", newline="") as fh:
        fh.write("value\n")
        for v in values:
            fh.write(f"{float(v):.17g}\n")


def read_values(oracle: str, fixture: str, variant: str) -> list[float] | None:
    p = GOLDEN / oracle / fixture / f"{variant}.csv"
    if not p.exists():
        return None
    with open(p) as fh:
        next(fh)
        return [float(line) for line in fh if line.strip()]
