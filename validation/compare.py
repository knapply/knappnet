#!/usr/bin/env python3
"""Compare the engine against every oracle's golden values.

    python3 validation/compare.py            # compare and write REPORT.md
    python3 validation/compare.py --strict   # also fail if an oracle never ran

Outcomes per (fixture, measure, oracle):
  match        values agree within the manifest tolerance
  documented   values differ on a fixture tagged with a convention difference
               the manifest declares for that oracle
  FAIL         any other disagreement (CI fails)
  not run      no golden file (oracle not available when goldens were made)
  undefined    no golden file because the oracle gives no answer on that
               fixture (manifest: undefined_on / undefined_when)
  manual       no golden file yet for an oracle run by hand (UCINET)

--strict fails only on "not run".
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (GOLDEN, HERE, applies, fixture_names, load_fixture,  # noqa: E402
                    load_manifest, read_values)
from engine import Graph, run_variant  # noqa: E402


# Oracles that cannot run unattended (UCINET is a Windows GUI program).
# Their golden files are produced by hand; until they exist, the missing
# comparisons are reported separately and do not fail --strict.
MANUAL_ORACLES = {"ucinet"}


def close(a: float, b: float, rel: float, abs_: float) -> bool:
    return abs(a - b) <= abs_ + rel * max(abs(a), abs(b))


def canonical(values, mode: str):
    """Puts results in a form where equal answers compare equal.

    values     as given (node-aligned vectors, scalars as one-element lists)
    sorted     order does not matter
    partition  group labels are arbitrary: relabel by first appearance
    sets       a collection of groups (such as cliques): order within and
               between groups does not matter
    """
    if mode == "values":
        return list(values)
    if mode == "sorted":
        return sorted(values)
    if mode == "partition":
        seen = {}
        return [float(seen.setdefault(v, len(seen))) for v in values]
    if mode == "sets":
        groups = [v if isinstance(v, tuple) else (v,) for v in values]
        return sorted(tuple(sorted(g)) for g in groups)
    raise ValueError(f"unknown comparison mode {mode!r}")


def agree(ours, theirs, rel: float, abs_: float) -> bool:
    if len(ours) != len(theirs):
        return False
    for a, b in zip(ours, theirs):
        if isinstance(a, tuple) or isinstance(b, tuple):
            a = a if isinstance(a, tuple) else (a,)
            b = b if isinstance(b, tuple) else (b,)
            if len(a) != len(b) or not all(close(x, y, rel, abs_) for x, y in zip(a, b)):
                return False
        elif not close(a, b, rel, abs_):
            return False
    return True


def worst_diff(xs: list[float], ys: list[float]) -> tuple[int, float, float]:
    i = max(range(len(xs)), key=lambda k: abs(xs[k] - ys[k]))
    return i, xs[i], ys[i]


def main() -> int:
    manifest = load_manifest()
    oracles = sorted({o for spec in manifest.values() for o in spec["oracles"]})
    tally: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    failures, documented, missing = [], [], []

    for name in fixture_names():
        fx = load_fixture(name)
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        g = Graph(fx.n, fx.src, fx.tgt, fx.weights, fx.directed)
        for variant, spec in manifest.items():
            if not applies(spec, fx):
                continue
            mode = spec.get("compare", "values")
            ours = canonical(run_variant(g, spec["engine"], fx), mode)
            tol = spec.get("tolerance", {})
            rel, abs_ = float(tol.get("rel", 0)), float(tol.get("abs", 0))
            for oracle, ospec in spec["oracles"].items():
                if (name in ospec.get("undefined_on", [])
                        or tags & set(ospec.get("undefined_when", []))):
                    # No well-defined answer to compare on this fixture, for
                    # the reason recorded in the manifest (whether or not the
                    # oracle wrote a value).
                    tally[(variant, oracle)]["undefined"] += 1
                    continue
                theirs = read_values(oracle, name, variant)
                if theirs is None:
                    key = "manual" if oracle in MANUAL_ORACLES else "not run"
                    tally[(variant, oracle)][key] += 1
                    if key == "not run":
                        missing.append((variant, oracle, name))
                    continue
                theirs = canonical(theirs, mode)
                if agree(ours, theirs, rel, abs_):
                    tally[(variant, oracle)]["match"] += 1
                    continue
                known = tags & set(ospec.get("differs_when", []))
                if name in ospec.get("differs_on", []):
                    known = known | {"listed"}
                numeric = all(not isinstance(v, tuple) for v in ours + theirs)
                if len(theirs) == len(ours) and numeric:
                    i, a, b = worst_diff(ours, theirs)
                    detail = f"index {i}: engine {a!r}, {oracle} {b!r}"
                elif len(theirs) == len(ours):
                    i = next(k for k, (a, b) in enumerate(zip(ours, theirs))
                             if not agree([a], [b], rel, abs_))
                    detail = f"item {i}: engine {ours[i]!r}, {oracle} {theirs[i]!r}"
                else:
                    detail = f"length {len(ours)} vs {len(theirs)}"
                row = (variant, oracle, name, detail)
                if known:
                    tally[(variant, oracle)]["documented"] += 1
                    documented.append(row + (", ".join(sorted(known)),))
                else:
                    tally[(variant, oracle)]["FAIL"] += 1
                    failures.append(row)

    lines = [
        "# Validation report",
        "",
        "Generated by `validation/compare.py`. Each cell counts fixtures as",
        "match / documented convention difference / FAIL / not run.",
        "",
    ]
    versions = []
    for o in oracles:
        vf = GOLDEN / o / "VERSION"
        versions.append(vf.read_text().strip() if vf.exists() else f"{o} (no goldens yet)")
    lines += ["Oracle versions: " + "; ".join(versions) + ".", ""]
    lines += ["| Measure | " + " | ".join(oracles) + " |",
              "|---|" + "---|" * len(oracles)]
    for variant in manifest:
        cells = []
        for o in oracles:
            if o not in manifest[variant]["oracles"]:
                cells.append("")
                continue
            t = tally[(variant, o)]
            cells.append(f"{t['match']} / {t['documented']} / {t['FAIL']} / {t['not run']}")
        lines.append(f"| `{variant}` | " + " | ".join(cells) + " |")

    if failures:
        lines += ["", "## Failures", ""]
        lines += [f"- `{v}` vs {o} on `{f}`: {d}" for v, o, f, d in failures]
    if documented:
        lines += ["", "## Documented convention differences", ""]
        lines += [f"- `{v}` vs {o} on `{f}` ({tags}): {d}"
                  for v, o, f, d, tags in documented]
    (HERE / "REPORT.md").write_text("\n".join(lines) + "\n", newline="")

    total = defaultdict(int)
    for t in tally.values():
        for k, v in t.items():
            total[k] += v
    print(f"match {total['match']}, documented {total['documented']}, "
          f"FAIL {total['FAIL']}, not run {total['not run']}, "
          f"awaiting a manual run {total['manual']}, "
          f"undefined for the oracle {total['undefined']}")
    for v, o, f, d in failures:
        print(f"FAIL {v} vs {o} on {f}: {d}")
    if "--strict" in sys.argv[1:] and total["not run"]:
        for v, o, f in missing:
            print(f"not run: {v} vs {o} on {f}")
    if failures:
        return 1
    if "--strict" in sys.argv[1:] and total["not run"]:
        print("Some oracles have not been run (--strict).")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
