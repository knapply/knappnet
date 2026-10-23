#!/usr/bin/env python3
"""Check every manifest entry (validation/manifest.d/*.yaml) for shape.

Each entry must have `engine` (with `fn`), `applies_to`, `tolerance` (with
`rel` and `abs`), and at least one oracle with a `call`. An oracle entry
that declares a known convention difference (`differs_when`, `differs_on`)
needs a `reason`; one that declares fixtures without an answer
(`undefined_on`, `undefined_when`) needs an `undefined_reason`. Fixture
tags and names must exist, `compare` must be a known mode, every
`engine.fn` must have a registered runner (validation/bindings/), and a
key may appear in only one file. See README.md for the format.

    python3 validation/check_manifest.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (COMPARE_MODES, TAGS, fixture_names,  # noqa: E402
                    load_manifest_by_area)

APPLIES = {"any", "directed", "undirected"}
ORACLE_KEYS = {"call", "differs_when", "differs_on", "reason", "undefined_on",
               "undefined_when", "undefined_reason", "note"}
ENTRY_KEYS = {"engine", "applies_to", "requires_weights", "tolerance", "compare",
              "oracles", "note"}


def check() -> list[str]:
    problems = []
    fixtures = set(fixture_names())
    try:
        areas = load_manifest_by_area()
    except ValueError as e:
        return [str(e)]
    try:
        from engine import RUNNERS, _load_bindings
        _load_bindings()
        runners = set(RUNNERS)
    except Exception as e:  # the bindings need the built library only to run
        problems.append(f"could not load the engine bindings: {e}")
        runners = None

    for area, entries in areas.items():
        where = f"manifest.d/{area}.yaml"
        if not isinstance(entries, dict):
            problems.append(f"{where}: not a mapping of measures")
            continue
        for key, spec in entries.items():
            p = f"{where}: {key}"
            if not isinstance(spec, dict):
                problems.append(f"{p}: not a mapping")
                continue
            for k in spec:
                if k not in ENTRY_KEYS:
                    problems.append(f"{p}: unknown key {k!r}")
            engine = spec.get("engine")
            if not isinstance(engine, dict) or "fn" not in engine:
                problems.append(f"{p}: needs engine: {{fn: ...}}")
            elif runners is not None and engine["fn"] not in runners:
                problems.append(f"{p}: no runner registered for engine fn {engine['fn']!r}")
            if spec.get("applies_to") not in APPLIES:
                problems.append(f"{p}: applies_to must be one of {sorted(APPLIES)}")
            tol = spec.get("tolerance")
            if not isinstance(tol, dict) or {"rel", "abs"} - set(tol):
                problems.append(f"{p}: needs tolerance: {{rel: ..., abs: ...}}")
            else:
                for k in ("rel", "abs"):
                    try:
                        if float(tol[k]) < 0:
                            problems.append(f"{p}: tolerance.{k} is negative")
                    except (TypeError, ValueError):
                        problems.append(f"{p}: tolerance.{k} is not a number")
            if spec.get("compare", "values") not in COMPARE_MODES:
                problems.append(f"{p}: compare must be one of {sorted(COMPARE_MODES)}")
            oracles = spec.get("oracles")
            if not isinstance(oracles, dict) or not oracles:
                problems.append(f"{p}: needs at least one oracle")
                continue
            for oracle, ospec in oracles.items():
                q = f"{p}: oracle {oracle}"
                if not isinstance(ospec, dict) or "call" not in ospec:
                    problems.append(f"{q}: needs call")
                    continue
                for k in ospec:
                    if k not in ORACLE_KEYS:
                        problems.append(f"{q}: unknown key {k!r}")
                if ("differs_when" in ospec or "differs_on" in ospec) and not ospec.get("reason"):
                    problems.append(f"{q}: differs_when/differs_on needs a reason")
                if "reason" in ospec and "differs_when" not in ospec and "differs_on" not in ospec:
                    problems.append(f"{q}: reason without differs_when or differs_on")
                if (("undefined_on" in ospec or "undefined_when" in ospec)
                        and not ospec.get("undefined_reason")):
                    problems.append(f"{q}: undefined_on/undefined_when needs an undefined_reason")
                if ("undefined_reason" in ospec and "undefined_on" not in ospec
                        and "undefined_when" not in ospec):
                    problems.append(f"{q}: undefined_reason without undefined_on or undefined_when")
                for k in ("differs_when", "undefined_when"):
                    for tag in ospec.get(k) or []:
                        if tag not in TAGS:
                            problems.append(f"{q}: unknown fixture tag {tag!r} in {k}")
                for k in ("differs_on", "undefined_on"):
                    for name in ospec.get(k) or []:
                        if name not in fixtures:
                            problems.append(f"{q}: unknown fixture {name!r} in {k}")
    return problems


def main() -> int:
    problems = check()
    for p in problems:
        print(p)
    n = sum(len(v) for v in load_manifest_by_area().values())
    print(f"{n} manifest entries, {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
