#!/usr/bin/env python3
"""Copy the Rust engine into the R package source tree.

An R source package must be self-contained, so the R package carries its own
copy of the engine under r/src/rust/. The canonical sources live in crates/.
Run this script after changing anything in crates/ or the root Cargo files.

    python3 tools/sync_r_sources.py           # update the copy
    python3 tools/sync_r_sources.py --check   # exit 1 if the copy is stale (CI)
"""

from __future__ import annotations

import filecmp
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "r" / "src" / "rust"

# (source, destination) pairs relative to ROOT and DEST.
FILES = [
    ("Cargo.toml", "Cargo.toml"),
    ("Cargo.lock", "Cargo.lock"),
    ("crates/engine/Cargo.toml", "crates/engine/Cargo.toml"),
    ("crates/capi/Cargo.toml", "crates/capi/Cargo.toml"),
]
DIRS = [
    ("crates/engine/src", "crates/engine/src"),
    ("crates/capi/src", "crates/capi/src"),
]
HEADER = ("crates/capi/include/netcapi.h", ROOT / "r" / "src" / "netcapi.h")


def expected() -> dict[Path, Path]:
    pairs = {DEST / dst: ROOT / src for src, dst in FILES}
    for src, dst in DIRS:
        for f in sorted((ROOT / src).rglob("*")):
            if f.is_file():
                pairs[DEST / dst / f.relative_to(ROOT / src)] = f
    pairs[HEADER[1]] = ROOT / HEADER[0]
    return pairs


def stale(pairs: dict[Path, Path]) -> list[str]:
    problems = []
    for dst, src in pairs.items():
        if not dst.exists() or not filecmp.cmp(src, dst, shallow=False):
            problems.append(str(dst.relative_to(ROOT)))
    for sub in ("crates/engine/src", "crates/capi/src"):
        for f in (DEST / sub).rglob("*"):
            if f.is_file() and f not in pairs:
                problems.append(f"{f.relative_to(ROOT)} (not in canonical sources)")
    return problems


def main() -> int:
    pairs = expected()
    if "--check" in sys.argv[1:]:
        problems = stale(pairs)
        if problems:
            print("R package sources are out of date. Run tools/sync_r_sources.py.")
            for p in problems:
                print(f"  {p}")
            return 1
        print("R package sources are in sync.")
        return 0

    for sub in ("crates/engine/src", "crates/capi/src"):
        shutil.rmtree(DEST / sub, ignore_errors=True)
    for dst, src in pairs.items():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    print(f"Copied {len(pairs)} files into {DEST.relative_to(ROOT)}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
