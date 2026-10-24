"""Writes docs/drafts/networkx_coverage.md: every NetworkX dispatchable
function with its status in the knappnet backend.

    python python/tools/networkx_coverage.py            # write the file
    python python/tools/networkx_coverage.py --check    # fail if it is stale

Status:
- implemented: every call NetworkX accepts runs on the engine (inputs on
  which NetworkX raises are left to NetworkX, which raises).
- partial: some arguments or graph types are left to NetworkX; they are
  listed.
- missing: always runs in NetworkX. Where the engine already has a
  related function, the last column says which, as a lead for future work.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

import networkx as nx

from knappnet import _algorithms  # noqa: F401
from knappnet._graph import IMPLEMENTATIONS, NOTES

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "drafts" / "networkx_coverage.md"

# Missing functions for which the C ABI already has a related routine, and
# why it is not wired yet.
ENGINE_LEADS = {
    "label_propagation_communities": "net_label_propagation (semi-synchronous); NetworkX "
    "visits nodes in set order, so results differ on ties",
    "asyn_lpa_communities": "net_label_propagation (retention); different random stream",
    "fast_label_propagation_communities": "net_label_propagation (fast); different random "
    "stream",
    "greedy_modularity_communities": "net_fast_greedy (Clauset-Newman-Moore); ties are "
    "broken differently, so partitions can differ",
    "leiden_communities": "net_leiden; different random stream",
    "asyn_fluidc": "net_fluid_communities; different random stream",
    "girvan_newman": "net_girvan_newman (removal order and best partition); NetworkX "
    "yields every level",
    "edge_betweenness_partition": "net_girvan_newman removal order",
    "k_core": "net_coreness (returns a graph)",
    "k_shell": "net_coreness (returns a graph)",
    "k_crust": "net_coreness (returns a graph)",
    "k_corona": "net_coreness (returns a graph)",
    "ego_graph": "net_ego_members (returns a graph)",
    "projected_graph": "net_project (returns a graph)",
    "weighted_projected_graph": "net_project, shared partners (returns a graph)",
    "collaboration_weighted_projected_graph": "net_project, Newman weights (returns a graph)",
    "generic_weighted_projected_graph": "net_project, weight products, for the default "
    "weight function only (returns a graph)",
    "density": "net_twomode_density (bipartite.density is O(1) in NetworkX)",
    "local_edge_connectivity": "net_lambda_connectivity (all pairs)",
    "edge_connectivity": "net_lambda_connectivity (all pairs)",
    "attribute_mixing_matrix": "net_mixing_matrix",
    "degree_mixing_matrix": "net_mixing_matrix",
    "sigma": "net_small_world (different random reference graphs)",
    "omega": "net_small_world (different random reference graphs)",
    "spring_layout": "net_layout_fr (different random start)",
    "gnp_random_graph": "net_random_gnp (different random stream)",
    "fast_gnp_random_graph": "net_random_gnp (different random stream)",
    "gnm_random_graph": "net_random_gnm (different random stream)",
    "dense_gnm_random_graph": "net_random_gnm (different random stream)",
    "configuration_model": "net_random_configuration (different random stream)",
    "watts_strogatz_graph": "net_random_small_world (different random stream)",
    "barabasi_albert_graph": "net_random_preferential (different random stream)",
    "stochastic_block_model": "net_random_sbm (different random stream)",
    "harmonic_diameter": "net_distances",
    "descendants_at_distance": "net_distances",
    "triangles": "triangle counts inside net_transitivity / net_clustering (not exported)",
    "max_weight_clique": "net_max_cliques (unweighted)",
    "node_clique_number": "net_max_cliques",
}


def group_of(module: str) -> str:
    parts = module.split(".")[1:]
    if parts[0] == "algorithms" and len(parts) > 2:
        return ".".join(parts[:2])
    if parts[0] in ("algorithms", "generators", "linalg", "readwrite", "drawing",
                    "classes", "convert", "convert_matrix", "relabel"):
        return ".".join(parts[:2]) if len(parts) > 1 else parts[0]
    return parts[0]


def render() -> str:
    algos = nx.utils.backends._registered_algorithms
    groups = collections.defaultdict(list)
    for name, func in algos.items():
        module = func.orig_func.__module__
        groups[group_of(module)].append((name, f"{module}.{func.orig_func.__name__}"))
    counts = collections.Counter()
    sections = []
    for group in sorted(groups):
        rows = []
        for name, path in sorted(groups[group]):
            if name in IMPLEMENTATIONS:
                declined = NOTES.get(name, "")
                status = "partial" if declined else "implemented"
                note = f"Left to NetworkX: {declined}" if declined else ""
            else:
                status = "missing"
                lead = ENGINE_LEADS.get(name)
                note = f"Engine: {lead}" if lead else ""
            counts[status] += 1
            rows.append(f"| `{name}` | `{path.removeprefix('networkx.')}` | {status} | "
                        f"{note} |")
        n_impl = sum(1 for r in rows if "| implemented |" in r or "| partial |" in r)
        sections.append(f"### {group} ({n_impl} of {len(rows)})\n\n"
                        "| Dispatch name | Function | Status | Notes |\n"
                        "|---|---|---|---|\n" + "\n".join(rows) + "\n")
    total = sum(counts.values())
    head = f"""# NetworkX coverage of the knappnet backend

Generated by `python/tools/networkx_coverage.py` from NetworkX
{nx.__version__}'s registry of dispatchable functions
(`networkx.utils.backends._registered_algorithms`) and the backend's table
of implementations. Do not edit by hand.

- **implemented**: every call NetworkX accepts runs on the engine. Inputs on
  which NetworkX raises an exception are left to NetworkX, which raises it.
- **partial**: the listed arguments, graph types, or inputs are left to
  NetworkX (the call still succeeds, in NetworkX's own implementation).
- **missing**: always runs in NetworkX. "Engine:" names a related routine
  of the C ABI where one exists.

"Weights <= 0" means graphs whose edge attribute has a zero or negative
value; "weight functions" means a callable passed as `weight`.

| Status | Functions |
|---|---|
| implemented | {counts['implemented']} |
| partial | {counts['partial']} |
| missing | {counts['missing']} |
| total | {total} |

"""
    return head + "\n".join(sections)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = render()
    if args.check:
        if not OUT.exists() or OUT.read_text() != text:
            print(f"{OUT} is out of date; run python/tools/networkx_coverage.py",
                  file=sys.stderr)
            return 1
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
