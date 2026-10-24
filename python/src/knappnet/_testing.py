"""Known differences from NetworkX's own test suite.

NetworkX's tests can be run with the backend (NETWORKX_TEST_BACKEND=knappnet,
NETWORKX_FALLBACK_TO_NX=True; see python/docs/testing.md). Tests listed here
check behaviour the backend deliberately does not share, such as the exact
partition NetworkX's random number stream produces for a seed. Each entry
maps a test node id suffix to the reason.
"""

_LOUVAIN = ("expects the partition NetworkX's random number stream gives for a seed; "
            "the backend uses the engine's stream")

KNOWN_DIFFERENCES = {
    "community/tests/test_louvain.py::test_karate_club_partition": _LOUVAIN,
    "community/tests/test_louvain.py::test_undirected_selfloops": _LOUVAIN,
    "community/tests/test_louvain.py::test_none_weight_param": _LOUVAIN,
    "community/tests/test_louvain.py::test_multigraph": _LOUVAIN,
    "community/tests/test_louvain.py::test_petersen": _LOUVAIN,
}


def mark_known_differences(items):
    import pytest

    for item in items:
        for suffix, reason in KNOWN_DIFFERENCES.items():
            if item.nodeid.endswith(suffix):
                item.add_marker(pytest.mark.xfail(reason=reason, strict=False))
