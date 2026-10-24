"""Backend metadata for NetworkX (the ``networkx.backend_info`` entry point).

NetworkX reads this when it is imported, to list the backend's functions
in its documentation and to set up ``nx.config.backends.knappnet``.
"""


def get_info():
    from . import _algorithms  # noqa: F401  (registers the implementations)
    from ._graph import NOTES

    return {
        "backend_name": "knappnet",
        "project": "knappnet",
        "package": "knappnet",
        "url": "https://github.com/knapply/knappnet",
        "short_summary": "Compiled, multithreaded engine; results match NetworkX.",
        "default_config": {
            # Threads for engine calls; 0 uses every core. Results do not
            # depend on the thread count.
            "threads": 0,
            # With automatic dispatch (backend_priority), graphs with fewer
            # nodes than this are left to NetworkX, where conversion would
            # cost more than it saves. Explicit backend="knappnet" ignores it.
            "min_nodes": 64,
        },
        "functions": {
            name: {"additional_docs": notes} if notes else {}
            for name, notes in sorted(NOTES.items())
        },
    }
