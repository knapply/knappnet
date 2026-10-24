# knappnet for Python: a NetworkX backend

knappnet runs [NetworkX](https://networkx.org/) algorithms on knappnet's
compiled, multithreaded engine. It adds no functions of its own: NetworkX
code stays as it is, and the functions knappnet implements return what
NetworkX returns, keyed by your node labels.

```python
import networkx as nx

G = nx.les_miserables_graph()
nx.betweenness_centrality(G, backend="knappnet")   # one call

# or every supported call, with no code changes:
#   NETWORKX_BACKEND_PRIORITY=knappnet python analysis.py
```

- [Installing from source](docs/install.md)
- [Using the backend](docs/usage.md)
- [What is supported](docs/supported.md)
- [Differences from NetworkX](docs/differences.md)
- [Testing](docs/testing.md)

knappnet is proprietary software; see [LICENSE](../LICENSE). It is not
published on PyPI.
