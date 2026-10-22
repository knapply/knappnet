Fast, validated network analysis in Rust.

This crate is the engine behind the knappnet R package. It stores graphs in
compressed sparse row form and implements network measures with parallel,
dependency-free code. Every measure is validated against established
reference implementations (igraph, sna, NetworkX) and against a direct
evaluation of its definition; see the
[validation report](https://github.com/knapply/knappnet/blob/main/validation/REPORT.md).

```rust
use knappnet::{betweenness, degree, triad_census, BetweennessOptions, DegreeMode, Graph, WeightKind};

// A path 0 - 1 - 2 - 3. Nodes are 0-based indices; edges are parallel arrays.
let g = Graph::from_edges(4, &[0, 1, 2], &[1, 2, 3], None, false)?;

assert_eq!(degree(&g, DegreeMode::All, true), vec![1.0, 2.0, 2.0, 1.0]);
assert_eq!(betweenness(&g, BetweennessOptions::default())?, vec![0.0, 2.0, 2.0, 0.0]);

// Weighted, directed, normalized, on all cores. Weights are distances by
// default; use `weights_as: WeightKind::Strength` for tie strengths.
let w = Graph::from_edges(3, &[0, 1, 0], &[1, 2, 2], Some(&[1.0, 1.0, 5.0]), true)?;
let opts = BetweennessOptions { weighted: true, normalized: true, ..Default::default() };
assert_eq!(betweenness(&w, opts)?, vec![0.0, 0.5, 0.0]);

// The same weights read as strengths: the strong direct tie 0 -> 2 wins.
let opts = BetweennessOptions { weighted: true, weights_as: WeightKind::Strength, ..Default::default() };
assert_eq!(betweenness(&w, opts)?, vec![0.0, 0.0, 0.0]);

// The 16 Holland-Leinhardt triad types, in MAN order.
let census = triad_census(&g, 0);
assert_eq!(census.iter().sum::<f64>(), 4.0);
# Ok::<(), knappnet::Error>(())
```

## Conventions

Each function documents how it treats self-loops, parallel edges, weights,
and normalization. Defaults follow igraph and NetworkX where they agree;
where reference tools disagree, the difference is recorded in the project's
validation manifest.

## Features

- Weighted and unweighted networks each get the fastest exact algorithm:
  breadth-first search for unweighted (or uniformly weighted) networks,
  Dial's bucket queue for small integer weights, and Dijkstra's algorithm
  otherwise. Weights can be read as distances or as tie strengths.
- No dependencies. Builds anywhere Rust's standard library does, including
  WebAssembly (where work runs on a single thread).
- Deterministic parallelism: results are reproducible for a given thread
  count, and integer results never depend on it.
- A C interface for other languages lives in the same repository.

## License

MIT
