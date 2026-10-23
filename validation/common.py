"""Shared helpers for fixtures, golden files, and the manifest.

Every Python oracle and engine binding imports from here, so the layout of
the corpus is defined once:

    fixtures/<name>/edges.csv, meta.json   the corpus (validation/fixtures.py)
    manifest.d/<area>.yaml                 what is compared, and how
    golden/<oracle>/<fixture>/<variant>.csv  one oracle's answers
    golden/<oracle>/VERSION                the oracle version that wrote them
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"
GOLDEN = HERE / "golden"
MANIFEST_DIR = HERE / "manifest.d"

# Fixture tags that manifest entries may name in `differs_when` and
# `undefined_when`. `directed` and `undirected` are added from meta.json;
# the rest are written by validation/fixtures.py.
TAGS = {"directed", "undirected", "weighted", "integer_weights", "loops",
        "multi_edges", "disconnected", "edge_case", "twomode"}

# Comparison modes a manifest entry may ask for (see compare.py).
COMPARE_MODES = {"values", "sorted", "partition", "sets"}


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


def manifest_files() -> list[Path]:
    """Every manifest file, in name order."""
    return sorted(MANIFEST_DIR.glob("*.yaml"))


def load_area(area: str) -> dict:
    """The manifest entries of one area (manifest.d/<area>.yaml)."""
    return yaml.safe_load((MANIFEST_DIR / f"{area}.yaml").read_text()) or {}


def load_manifest() -> dict:
    """Every manifest entry, in file order; a key may appear only once."""
    manifest: dict = {}
    for path in manifest_files():
        for key, spec in (yaml.safe_load(path.read_text()) or {}).items():
            if key in manifest:
                raise ValueError(f"measure {key!r} is defined twice (in {path.name})")
            manifest[key] = spec
    return manifest


def load_manifest_by_area() -> dict[str, dict]:
    """Area (manifest file stem) -> its entries, in file order."""
    areas: dict[str, dict] = {}
    seen: set[str] = set()
    for path in manifest_files():
        entries = yaml.safe_load(path.read_text()) or {}
        for key in entries:
            if key in seen:
                raise ValueError(f"measure {key!r} is defined twice (in {path.name})")
            seen.add(key)
        areas[path.stem] = entries
    return areas


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


def write_version(oracle: str, version: str) -> None:
    """Records which oracle version produced golden/<oracle>/ (shown in
    REPORT.md). Every oracle script calls this once, after writing."""
    (GOLDEN / oracle).mkdir(parents=True, exist_ok=True)
    (GOLDEN / oracle / "VERSION").write_text(version.rstrip("\n") + "\n", newline="")


DEFINITION_VERSION = "definition (direct evaluation, exact arithmetic)"


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
