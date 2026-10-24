"""Argument checks for ``can_run``.

Each check takes the call's arguments bound to NetworkX's signature and
returns None (the backend can run the call) or a string saying why not;
NetworkX then runs its own implementation. A check declines whatever the
backend cannot reproduce exactly: options it does not implement, graph
types it does not handle, and inputs on which NetworkX raises an exception
(so that NetworkX raises it).
"""

from __future__ import annotations

from collections.abc import Collection

import networkx as nx


def graph_of(a, name="G"):
    G = a.get(name)
    if G is None and name == "G":
        # Functions whose graph parameter has another name.
        G = next((a[k] for k in ("flowG", "FG", "G1") if a.get(k) is not None), None)
    if getattr(G, "__networkx_backend__", None) == "knappnet":
        return G.graph
    return G


# Methods through which a subclass could present a graph other than the one
# stored in its adjacency dicts (NetworkX's approximation module has such a
# class for the complement graph). Subclasses overriding any are declined.
_STORAGE_METHODS = ("__getitem__", "__iter__", "__len__", "__contains__", "adj", "succ",
                    "pred", "neighbors", "successors", "predecessors", "edges", "degree",
                    "nodes", "has_edge", "has_node", "number_of_edges", "_adj", "_succ",
                    "_pred", "_node")
_BASES = None


def plain_storage(G):
    """True if G's graph is exactly what its adjacency dicts hold."""
    global _BASES
    if _BASES is None:
        _BASES = (nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph)
    cls = type(G)
    if cls in _BASES:
        return True
    base = next((b for b in cls.__mro__ if b in _BASES), None)
    if base is None:
        return False
    for klass in cls.__mro__:
        if klass is base:
            break
        if any(name in vars(klass) for name in _STORAGE_METHODS):
            return False
    return True


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


def nodes_in_graph(*params, name="G"):
    """Every named argument is a node of the graph."""
    checks = [node_in_graph(p, name) for p in params]
    return all_of(*checks)


def nodes_all_in_graph(param, name="G", allow_empty=False):
    """The argument is a non-empty container of nodes of the graph."""

    def check(a):
        G = graph_of(a, name)
        value = a.get(param)
        if not isinstance(value, Collection) or isinstance(value, (str, bytes)):
            return f"{param} must be a collection of nodes"
        items = list(value)
        if not items and not allow_empty:
            return f"{param} is empty"
        try:
            if all(v in G for v in items):
                return None
        except TypeError:
            pass
        return f"{param} has an item that is not a node of the graph"

    return check


def numeric_or_none(param):
    def check(a):
        v = a.get(param)
        if v is None or (isinstance(v, (int, float)) and not isinstance(v, bool)):
            return None
        return f"{param} must be a number or None"

    return check


def optional_nodes(param, name="G"):
    """None, or a non-empty collection of distinct nodes of the graph."""
    inner = nodes_all_in_graph(param, name)

    def check(a):
        value = a.get(param)
        if value is None:
            return None
        reason = inner(a)
        if reason:
            return reason
        if len(set(value)) != len(value):
            return f"{param} repeats a node"
        return None

    return check


def choice(param, *values):
    """The argument is one of the given strings."""

    def check(a):
        v = a.get(param)
        if isinstance(v, str) and v in values:
            return None
        return f"{param} is not one of {', '.join(values)}"

    return check


def integer(param, minimum=None):
    """The argument is an int (not a bool), at least ``minimum`` if given."""

    def check(a):
        v = a.get(param)
        if isinstance(v, int) and not isinstance(v, bool) and (minimum is None or v >= minimum):
            return None
        return f"{param} must be an integer" + (f" of at least {minimum}" if minimum is not None else "")

    return check


def optional_integer(param, minimum=None):
    inner = integer(param, minimum)

    def check(a):
        return None if a.get(param) is None else inner(a)

    return check
