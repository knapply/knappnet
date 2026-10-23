#!/usr/bin/env python3
"""Validation at scale: the engine against NetworkX (and, through
check_scale.R, igraph and sna) on larger networks.

    python3 validation/scale/check_scale.py                 # NetworkX
    python3 validation/scale/check_scale.py --prepare-r DIR # inputs for check_scale.R
    python3 validation/scale/check_scale.py --compare-r DIR # igraph and sna results

The networks come from validation/scale/generate.py (seeded, a few
hundred to 50,000 nodes). Nothing about the comparison is specific to
this tier:

- the measures, engine calls, tolerances, and documented convention
  differences are read from validation/manifest.yaml and manifest.d/, the
  single source of truth also used by validation/compare.py;
- the oracle values come from the existing oracle scripts
  (validation/oracles/*networkx_oracle.py, and the R scripts run by
  check_scale.R), pointed at the scale fixtures instead of the main ones;
- results are compared with compare.py's own functions.

The definition oracles are not run: they evaluate definitions exactly and
are far too slow for networks of this size.

Each measure runs only on fixtures within its cost budget (BUDGETS below),
chosen so the whole check takes a few minutes on a 2-core machine:
shortest-path measures in NetworkX (pure Python, O(nm)) run on networks
of up to about 1,500 nodes, linear-time measures on every network.

Exit status: 1 if any comparison fails (a disagreement on a fixture
whose tags the manifest does not declare as a known convention
difference), 0 otherwise. Comparisons for which an oracle gives no value
(it declines the network, or does not implement that variant) are
counted as "not run" and listed, but do not fail.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import runpy
import shutil
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
VALIDATION = HERE.parent
sys.path.insert(0, str(VALIDATION))
sys.path.insert(0, str(HERE))

import common  # noqa: E402
from compare import agree, canonical, worst_diff  # noqa: E402
from engine import Graph, run_variant  # noqa: E402

import generate  # noqa: E402

WORK = HERE / "work"
ORACLES = VALIDATION / "oracles"

# Cost budgets. A measure runs on a fixture when the fixture is within its
# budget: ("nm", x) means n * m <= x, ("n+m", x) means n + m <= x, ("n", x)
# means n <= x. A budget may instead be a dict with "undirected" and
# "directed" entries (a missing entry means the measure never runs on
# that kind of network), for NetworkX routines whose cost differs between
# the two. BUDGETS is keyed by engine function; VARIANT_BUDGETS overrides
# it for the variants whose engine parameters match, and functions listed
# in neither get the shortest-path budget, the most conservative.
PATHS = ("nm", 1.5e7)
LINEAR = ("n+m", 1e7)
BUDGETS = {
    # Near-linear in NetworkX: every fixture.
    "degree": LINEAR,
    "strength": LINEAR,
    "degree_simple": LINEAR,
    "density": LINEAR,
    "components": LINEAR,
    "fragmentation": LINEAR,
    "coreness": LINEAR,
    "max_cliques": LINEAR,
    "clique_number": LINEAR,
    "modularity": LINEAR,
    "louvain_modularity": LINEAR,
    "eigenvector": LINEAR,
    "triad_census": LINEAR,
    "transitivity": LINEAR,
    "clustering": LINEAR,
    "average_clustering": LINEAR,
    "reciprocity": LINEAR,
    "bridges": LINEAR,
    "cut_points": LINEAR,
    "pagerank": LINEAR,
    "hits": LINEAR,
    "assortativity_degree": LINEAR,
    "assortativity_nominal": LINEAR,
    "assortativity_values": LINEAR,
    "project": LINEAR,
    # Dense n x n matrices in the oracle (Katz: an eigenvalue problem).
    "katz": ("n", 1500),
    "constraint": ("n", 1500),
    "effective_size": ("n", 1500),
    "efficiency": ("n", 1500),
    # Twenty runs of NetworkX's Louvain.
    "louvain_quality": ("n+m", 15000),
    # Ego networks: nx.ego_graph is linear per node on undirected
    # networks, but converts the whole network for every node when it
    # reads a directed one in both directions.
    "ego_size": {"undirected": ("n", 1500), "directed": ("n", 600)},
    "ego_members": {"undirected": ("n", 1500), "directed": ("n", 600)},
    "ego_summary": {"undirected": ("n", 1500), "directed": ("n", 600)},
    # O(nm) shortest paths; two-mode centralities are shortest-path based.
    "betweenness": PATHS,
    "betweenness_simple": PATHS,
    "closeness": PATHS,
    "harmonic": PATHS,
    "centralization": PATHS,
    "distances": PATHS,
    "diameter": PATHS,
    "eccentricity": PATHS,
    "mean_distance": PATHS,
    "reachability": PATHS,
    "twomode": PATHS,
}
VARIANT_BUDGETS = [
    # Betweenness within every order-2 ego network: on a directed network
    # that is minutes even at 600 nodes.
    ({"fn": "ego_summary", "measure": "betweenness", "order": 2},
     {"undirected": ("n", 1000)}),
]


def budget_of(spec: dict) -> dict | tuple:
    eng = spec["engine"]
    for match, budget in VARIANT_BUDGETS:
        if all(eng.get(k) == v for k, v in match.items()):
            return budget
    return BUDGETS.get(eng["fn"], PATHS)


# Convention differences that only networks of the scale tier exhibit,
# declared here as the manifest declares the others (the tags come from
# generate.py). Each entry: the engine functions, whether weighted, the
# oracles, the fixture tags, and the reason.
SCALE_DIFFERENCES = [
    ({"betweenness", "betweenness_simple"}, True, {"networkx"}, {"near_ties"},
     "the engine, like igraph, counts shortest paths whose lengths agree to a relative "
     "tolerance of 1e-10 as tied; NetworkX compares lengths exactly, so on networks with "
     "nearly equal path lengths (generate.py's near_ties tag) it counts fewer shortest paths"),
]


def scale_difference(spec: dict, oracle: str, tags: set) -> set:
    """Tags under which `oracle` is known to differ on `spec` in this tier."""
    eng = spec["engine"]
    out = set()
    for fns, weighted, oracles, when, _ in SCALE_DIFFERENCES:
        if eng["fn"] in fns and bool(eng.get("weighted")) == weighted and oracle in oracles:
            out |= tags & when
    return out


# The R oracle scripts build dense adjacency matrices for sna, so the R
# tier uses networks of up to R_MAX_NODES nodes. igraph runs every
# measure there; sna's slower matrix routines get their own budgets.
R_MAX_NODES = 5000
SNA_BUDGETS = {
    "triad_census": ("n", 600),
    "brokerage": ("n", 600),
    "structural_equivalence": ("n", 600),
    "concor": ("n", 600),
    "centralization": ("n", 1500),
    "closeness": ("n", 1500),
    "harmonic": ("n", 1500),
    "betweenness": ("n", 1500),
    "betweenness_simple": ("n", 1500),
    "eigenvector": ("n", 1500),
    "degree": ("n", 5000),
    "degree_simple": ("n", 5000),
    "density": ("n", 5000),
}
SNA_DEFAULT = ("n", 600)


def undefined_by_name(spec: dict, ospec: dict, tags: set) -> bool:
    """The manifests declare an oracle undefined on main fixtures by name
    (`undefined_on`). The one such rule that carries over to the scale
    fixtures is that the two-mode measures are undefined on one-mode
    networks, which have no node types; here the `twomode` tag tells the
    two apart. The others name degenerate main fixtures (empty, pair)."""
    return bool(ospec.get("undefined_on")) and spec["engine"]["fn"] in {"twomode", "project"} \
        and "twomode" not in tags


def within(budget, fx) -> bool:
    if isinstance(budget, dict):
        budget = budget.get("directed" if fx.directed else "undirected")
        if budget is None:
            return False
    kind, limit = budget
    n, m = fx.n, len(fx.src)
    size = {"nm": n * m, "n+m": n + m, "n": n}[kind]
    return size <= limit


def in_budget(spec: dict, fx) -> bool:
    return within(budget_of(spec), fx)


def generate_fixtures(out: Path) -> None:
    t = time.time()
    generate.write(generate.build(), out)
    print(f"generated {len(list(out.iterdir()))} scale fixtures in {time.time() - t:.1f} s")


def run_oracle(script: str, fixtures: str, golden: str, name: str) -> tuple[str, str, float]:
    """Runs one NetworkX oracle script on one scale fixture (in a worker
    process). The scripts take their paths, fixture list, and `applies`
    test from validation/common.py, which is pointed at that fixture of the
    scale tier (with the cost budgets added) before the script is loaded."""
    started = time.time()
    common.FIXTURES = Path(fixtures)
    common.GOLDEN = Path(golden)
    plain_applies, plain_names = common.applies, common.fixture_names
    common.applies = lambda spec, fx: plain_applies(spec, fx) and in_budget(spec, fx)
    common.fixture_names = lambda: [name]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            runpy.run_path(script, run_name="__main__")
    finally:
        common.applies, common.fixture_names = plain_applies, plain_names
    return Path(script).name, name, time.time() - started


def run_python_oracles(fixtures: Path, golden: Path) -> None:
    """Runs every NetworkX oracle script on every scale fixture, writing
    golden files to `golden`. Each (script, fixture) pair is a separate
    task, run on all cores, the largest fixtures first."""
    shutil.rmtree(golden, ignore_errors=True)
    common.FIXTURES = fixtures
    sizes = {}
    for name in common.fixture_names():
        meta = json.loads((fixtures / name / "meta.json").read_text())
        sizes[name] = meta["n"] + sum(1 for _ in open(fixtures / name / "edges.csv"))
    scripts = sorted(str(p) for p in ORACLES.glob("*networkx_oracle.py"))
    tasks = [(s, str(fixtures), str(golden), name)
             for name in sorted(sizes, key=sizes.get, reverse=True) for s in scripts]
    per_script = defaultdict(float)
    started = time.time()
    with ProcessPoolExecutor(max_workers=os.cpu_count() or 1) as pool:
        for script, name, seconds in pool.map(run_oracle, *zip(*tasks)):
            per_script[script] += seconds
    for script, seconds in sorted(per_script.items()):
        print(f"   {script}: {seconds:.0f} s")
    print(f"NetworkX oracles took {time.time() - started:.0f} s on "
          f"{os.cpu_count() or 1} cores")


def compare_tree(fixtures: Path, golden: Path, oracles: set[str], budget, tally, failures,
                 documented, missing) -> None:
    """Compares the engine with each oracle's golden files under `golden`
    for every fixture in `fixtures`, as validation/compare.py does."""
    common.FIXTURES = fixtures
    common.GOLDEN = golden
    manifest = common.load_manifest()
    for name in common.fixture_names():
        fx = common.load_fixture(name)
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        g = None
        for variant, spec in manifest.items():
            wanted = [o for o in spec["oracles"] if o in oracles and budget(o, spec, fx)]
            if not wanted or not common.applies(spec, fx):
                continue
            mode = spec.get("compare", "values")
            if g is None:
                g = Graph(fx.n, fx.src, fx.tgt, fx.weights, fx.directed)
            ours = canonical(run_variant(g, spec["engine"], fx), mode)
            tol = spec.get("tolerance", {})
            rel, abs_ = float(tol.get("rel", 0)), float(tol.get("abs", 0))
            for oracle in wanted:
                ospec = spec["oracles"][oracle]
                theirs = common.read_values(oracle, name, variant)
                if theirs is None:
                    if tags & set(ospec.get("undefined_when", [])) or \
                            undefined_by_name(spec, ospec, tags):
                        tally[(variant, oracle)]["undefined"] += 1
                    else:
                        tally[(variant, oracle)]["not run"] += 1
                        missing.append((variant, oracle, name))
                    continue
                theirs = canonical(theirs, mode)
                if agree(ours, theirs, rel, abs_):
                    tally[(variant, oracle)]["match"] += 1
                    continue
                numeric = all(not isinstance(v, tuple) for v in ours + theirs)
                if len(theirs) == len(ours) and numeric:
                    i, a, b = worst_diff(ours, theirs)
                    detail = f"index {i}: engine {a!r}, {oracle} {b!r}"
                else:
                    detail = f"length {len(ours)} vs {len(theirs)}" \
                        if len(theirs) != len(ours) else "groups differ"
                known = (tags & set(ospec.get("differs_when", []))) \
                    | scale_difference(spec, oracle, tags)
                if known:
                    tally[(variant, oracle)]["documented"] += 1
                    documented.append((variant, oracle, name, detail, ", ".join(sorted(known))))
                else:
                    tally[(variant, oracle)]["FAIL"] += 1
                    failures.append((variant, oracle, name, detail))


def report(tally, failures, documented, missing, started: float) -> int:
    total = defaultdict(int)
    print(f"\n{'measure':<40} {'oracle':<9} match documented FAIL notrun")
    for (variant, oracle), t in sorted(tally.items()):
        for k, v in t.items():
            total[k] += v
        print(f"{variant:<40} {oracle:<9} {t['match']:>5} {t['documented']:>10} "
              f"{t['FAIL']:>4} {t['not run']:>6}")
    for v, o, f, d, tags in documented:
        print(f"documented: {v} vs {o} on {f} ({tags}): {d}")
    for v, o, f in missing:
        print(f"not run: {v} vs {o} on {f}")
    for v, o, f, d in failures:
        print(f"FAIL {v} vs {o} on {f}: {d}")
    print(f"\nmatch {total['match']}, documented {total['documented']}, FAIL {total['FAIL']}, "
          f"not run {total['not run']}, undefined for the oracle {total['undefined']}; "
          f"{time.time() - started:.0f} s")
    return 1 if failures else 0


def networkx_check() -> int:
    started = time.time()
    fixtures, golden = WORK / "fixtures", WORK / "golden"
    generate_fixtures(fixtures)
    run_python_oracles(fixtures, golden)
    tally = defaultdict(lambda: defaultdict(int))
    failures, documented, missing = [], [], []
    compare_tree(fixtures, golden, {"networkx"}, lambda o, spec, fx: in_budget(spec, fx),
                 tally, failures, documented, missing)
    return report(tally, failures, documented, missing, started)


def sna_ok(spec: dict, fx) -> bool:
    return within(SNA_BUDGETS.get(spec["engine"]["fn"], SNA_DEFAULT), fx)


def r_budget(oracle: str, spec: dict, fx) -> bool:
    if fx.n > R_MAX_NODES:
        return False
    return oracle != "sna" or sna_ok(spec, fx)


def prepare_r(out: Path) -> int:
    """Writes one directory per fixture (up to R_MAX_NODES nodes) laid out
    like the repository's validation/ folder, holding that fixture and the
    manifests restricted to what igraph and sna run on it, so the R oracle
    scripts can run there unchanged."""
    shutil.rmtree(out, ignore_errors=True)
    fixtures = out / "_fixtures"
    generate_fixtures(fixtures)
    common.FIXTURES = fixtures
    sources = [VALIDATION / "manifest.yaml"] + sorted((VALIDATION / "manifest.d").glob("*.yaml"))
    for name in common.fixture_names():
        fx = common.load_fixture(name)
        if fx.n > R_MAX_NODES:
            continue
        root = out / name / "validation"
        shutil.copytree(fixtures / name, root / "fixtures" / name)
        for src in sources:
            specs = yaml.safe_load(src.read_text()) or {}
            kept = {}
            for variant, spec in specs.items():
                oracles = {o: v for o, v in spec["oracles"].items()
                           if (o == "igraph" or (o == "sna" and sna_ok(spec, fx)))}
                if oracles:
                    kept[variant] = dict(spec, oracles=oracles)
            dest = root / src.relative_to(VALIDATION)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(yaml.safe_dump(kept, sort_keys=False) if kept else "{}\n")
        (root / "golden").mkdir()
    print(f"prepared R inputs in {out}")
    return 0


def compare_r(out: Path) -> int:
    started = time.time()
    tally = defaultdict(lambda: defaultdict(int))
    failures, documented, missing = [], [], []
    trees = sorted(p for p in out.iterdir() if (p / "validation" / "fixtures").is_dir())
    for tree in trees:
        compare_tree(tree / "validation" / "fixtures", tree / "validation" / "golden",
                     {"igraph", "sna"}, r_budget, tally, failures, documented, missing)
    for oracle in ("igraph", "sna"):
        versions = {(t / "validation" / "golden" / oracle / "VERSION").read_text().strip()
                    for t in trees if (t / "validation" / "golden" / oracle / "VERSION").exists()}
        print(f"{oracle}: {', '.join(sorted(versions)) or 'no golden files'}")
    if not trees:
        print(f"no prepared fixtures in {out}")
        return 1
    status = report(tally, failures, documented, missing, started)
    if not any(t["match"] for (_, o), t in tally.items() if o == "igraph"):
        print("no igraph comparison matched: did the R oracle scripts run?")
        return 1
    return status


def main() -> int:
    args = sys.argv[1:]
    if args[:1] == ["--prepare-r"] and len(args) == 2:
        return prepare_r(Path(args[1]))
    if args[:1] == ["--compare-r"] and len(args) == 2:
        return compare_r(Path(args[1]))
    if args:
        print(__doc__)
        return 2
    return networkx_check()


if __name__ == "__main__":
    sys.exit(main())
