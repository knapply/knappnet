"""Argument checks for ``can_run``.

Each check takes the call's arguments bound to NetworkX's signature and
returns None (the backend can run the call) or a string saying why not;
NetworkX then runs its own implementation. A check declines whatever the
backend cannot reproduce exactly: options it does not implement, graph
types it does not handle, and inputs on which NetworkX raises an exception
(so that NetworkX raises it).
"""

from __future__ import annotations

import networkx as nx


def graph_of(a, name="G"):
    G = a.get(name)
    if getattr(G, "__networkx_backend__", None) == "knappnet":
        return G.graph
    return G


def all_of(*checks):
    def check(a):
        for c in checks:
            reason = c(a)
            if reason:
                return reason
        return None

    return check


def no_multigraph(a, name="G"):
    if graph_of(a, name).is_multigraph():
        return "multigraphs are not supported"
    return None


def undirected(a, name="G"):
    if graph_of(a, name).is_directed():
        return "NetworkX does not define this for directed graphs"
    return None


def directed(a, name="G"):
    if not graph_of(a, name).is_directed():
        return "NetworkX does not define this for undirected graphs"
    return None


def nonempty(a, name="G"):
    if len(graph_of(a, name)) == 0:
        return "the graph has no nodes"
    return None


def no_selfloops(a, name="G"):
    if nx.number_of_selfloops(graph_of(a, name)) > 0:
        return "self-loops are handled by NetworkX"
    return None


def weight_attr(param="weight"):
    """The weight argument must be None or an attribute name (not a function)."""

    def check(a):
        w = a.get(param)
        if w is None or isinstance(w, str):
            return None
        return f"{param} must be None or an attribute name"

    return check


def equals(param, *values):
    """The argument must be one of the given values (booleans only match booleans)."""

    def check(a):
        v = a.get(param)
        for x in values:
            if v is x:
                return None
            if isinstance(v, bool) or isinstance(x, bool) or v is None or x is None:
                continue
            try:
                if v == x:
                    return None
            except Exception:
                pass
        return f"{param}={v!r} is not supported"

    return check


def is_none(*params):
    def check(a):
        for p in params:
            if a.get(p) is not None:
                return f"{p} is not supported"
        return None

    return check


def at_least(param, minimum):
    def check(a):
        v = a.get(param)
        try:
            if v >= minimum:
                return None
        except TypeError:
            pass
        return f"{param} below {minimum} is left to NetworkX"

    return check


def node_in_graph(param, name="G"):
    def check(a):
        G = graph_of(a, name)
        try:
            if a.get(param) in G:
                return None
        except TypeError:
            pass
        return f"{param} is not a node of the graph"

    return check
