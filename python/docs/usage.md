# Using the backend

knappnet is a NetworkX backend: NetworkX decides, call by call, whether a
function runs on knappnet's engine or in NetworkX itself. There is nothing
to import from knappnet and no graph class to learn. It works on NetworkX's
`Graph`, `DiGraph`, `MultiGraph`, and `MultiDiGraph` (and their views).

## One call

Pass `backend="knappnet"` to any supported function:

```python
import networkx as nx

G = nx.karate_club_graph()
bc = nx.betweenness_centrality(G, weight="weight", backend="knappnet")
```

If knappnet cannot run that call exactly as NetworkX would (an argument it
does not implement, such as `k=` sampling), NetworkX raises
`NotImplementedError` saying so, because the backend was requested
explicitly. To have such calls run in NetworkX instead, use the backend
priority described next.

## Every call

Make knappnet the first choice for all dispatchable functions with an
environment variable (read when NetworkX is imported):

```sh
NETWORKX_BACKEND_PRIORITY=knappnet python analysis.py
```

or in code:

```python
nx.config.backend_priority.algos = ["knappnet"]
```

NetworkX then tries knappnet first and falls back to its own implementation
whenever knappnet declines: for functions knappnet does not implement, for
arguments or inputs it does not handle, and for graphs too small to
benefit (see `min_nodes` below). Linear-time functions (components,
degree centrality, core numbers, single-source path lengths, and the
like), which NetworkX runs in about the time a conversion takes, run on
the engine only once the graph has been converted by an earlier call.
Results are the same either way.

`backend=` on a call still overrides the priority list, so
`backend="networkx"` forces NetworkX's own implementation.

## Converted graphs are cached

The first call converts the graph (one pass over its edges, reading at most
one edge attribute); NetworkX keeps the result in `G.__networkx_cache__`,
so later calls on the same graph reuse it. Methods such as `G.add_edge()`
clear the cache. Changing edge data in place (`G[u][v]["weight"] = 2`)
does not, and NetworkX warns about this when it reuses a cached
conversion; clear it with `G.__networkx_cache__.clear()`, or silence the
warning with:

```python
nx.config.warnings_to_ignore.add("cache")
```

Caching is controlled by NetworkX (`nx.config.cache_converted_graphs`,
or the environment variable `NETWORKX_CACHE_CONVERTED_GRAPHS`).

## Configuration

`nx.config.backends.knappnet` holds two settings:

| Setting | Default | Meaning |
|---|---|---|
| `threads` | `0` | Threads for engine calls; 0 uses every core. Results do not depend on the thread count. |
| `min_nodes` | `64` | With automatic dispatch (backend priority), graphs with fewer nodes run in NetworkX, where conversion would cost more than it saves. An explicit `backend="knappnet"` ignores it. |

```python
nx.config.backends.knappnet.threads = 4
```

## Seeing what ran where

NetworkX logs every dispatch decision, including the reason knappnet gave
for declining a call:

```python
import logging
logging.basicConfig()
logging.getLogger("networkx").setLevel(logging.DEBUG)
```
