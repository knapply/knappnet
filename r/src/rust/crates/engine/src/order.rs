//! The order in which traversals visit neighbours.
//!
//! Most measures do not depend on the order in which a node's neighbours
//! are visited, but some results do: which spanning tree a depth-first
//! search finds, the order of strongly connected components, a
//! topological order among several, the cycles of a cycle basis, the
//! cycle `find_cycle` reports. For these the engine offers two orders:
//!
//! - [`NeighborOrder::Sorted`]: neighbours in increasing index, the order
//!   of the graph's own adjacency lists (the default).
//! - [`NeighborOrder::Input`]: neighbours in the order of the edges that
//!   join them, as given. NetworkX (`hagberg2008exploring`) stores a
//!   node's neighbours in the order they were first added, so a NetworkX
//!   graph built by adding nodes 0..n and then the same edges in the same
//!   order visits neighbours in exactly this order; with it, traversal-
//!   dependent results reproduce NetworkX's.
//!
//! Rows are built in O(n + m) with counting sorts; parallel edges either
//! repeat a neighbour or, when rows are deduplicated (the simple-graph
//! reading NetworkX's `Graph` and `DiGraph` give), keep only its first
//! occurrence, which is where NetworkX keeps it.

use crate::error::{Error, Result};
use crate::graph::Graph;

/// Neighbour order for traversal-dependent results. See the module
/// documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NeighborOrder {
    /// Increasing neighbour index.
    #[default]
    Sorted,
    /// Order of the joining edges in the input (NetworkX's insertion
    /// order).
    Input,
}

impl NeighborOrder {
    /// Decodes the C interface codes: 0 = sorted, 1 = input order.
    pub fn from_code(code: i32) -> Result<NeighborOrder> {
        match code {
            0 => Ok(NeighborOrder::Sorted),
            1 => Ok(NeighborOrder::Input),
            _ => Err(Error::InvalidArgument(format!(
                "unknown neighbor order code {code}"
            ))),
        }
    }
}

/// Which arcs a row holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Arcs leaving the node (all incident edges if undirected).
    Out,
    /// Arcs entering the node (all incident edges if undirected).
    In,
    /// Arcs in either direction.
    Both,
}

/// Adjacency rows with the edge behind each entry.
#[derive(Debug, Clone)]
pub struct Rows {
    offsets: Vec<usize>,
    nodes: Vec<u32>,
    edges: Vec<u32>,
}

impl Rows {
    /// Neighbours of `v` in the chosen order.
    #[inline]
    pub fn neighbors(&self, v: usize) -> &[u32] {
        &self.nodes[self.offsets[v]..self.offsets[v + 1]]
    }

    /// Edges (input positions) behind `neighbors(v)`, aligned with it.
    #[inline]
    pub fn edges(&self, v: usize) -> &[u32] {
        &self.edges[self.offsets[v]..self.offsets[v + 1]]
    }

    /// Number of rows.
    pub fn node_count(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Number of entries in `v`'s row.
    #[inline]
    pub fn degree(&self, v: usize) -> usize {
        self.offsets[v + 1] - self.offsets[v]
    }
}

/// Rows of `g` in `order`. Self-loops appear once in their node's row;
/// with `skip_loops` they are left out. With `dedupe`, a neighbour appears
/// only at its first position (in `order`).
pub fn rows(
    g: &Graph,
    direction: Direction,
    order: NeighborOrder,
    skip_loops: bool,
    dedupe: bool,
) -> Rows {
    let n = g.node_count();
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let m = src.len();
    // Entries (row, neighbour, edge) in input order.
    let mut entries: Vec<(u32, u32, u32)> = Vec::with_capacity(2 * m);
    let both_ways = !g.is_directed() || direction == Direction::Both;
    for e in 0..m {
        let (a, b) = (src[e], tgt[e]);
        if a == b {
            if !skip_loops {
                entries.push((a, a, e as u32));
            }
            continue;
        }
        if both_ways {
            entries.push((a, b, e as u32));
            entries.push((b, a, e as u32));
        } else if direction == Direction::Out {
            entries.push((a, b, e as u32));
        } else {
            entries.push((b, a, e as u32));
        }
    }
    if order == NeighborOrder::Sorted {
        // Stable counting sort by neighbour keeps edge order among
        // parallel edges.
        entries = counting_sort(entries, n, |x| x.1 as usize);
    }
    let entries = counting_sort(entries, n, |x| x.0 as usize);
    let mut offsets = vec![0usize; n + 1];
    let mut nodes = Vec::with_capacity(entries.len());
    let mut edges = Vec::with_capacity(entries.len());
    let mut mark = vec![u32::MAX; if dedupe { n } else { 0 }];
    let mut i = 0;
    for v in 0..n {
        while i < entries.len() && entries[i].0 as usize == v {
            let (_, w, e) = entries[i];
            i += 1;
            if dedupe {
                if mark[w as usize] == v as u32 {
                    continue;
                }
                mark[w as usize] = v as u32;
            }
            nodes.push(w);
            edges.push(e);
        }
        offsets[v + 1] = nodes.len();
    }
    Rows {
        offsets,
        nodes,
        edges,
    }
}

fn counting_sort<T: Copy>(items: Vec<T>, n: usize, key: impl Fn(&T) -> usize) -> Vec<T> {
    let mut count = vec![0usize; n + 1];
    for x in &items {
        count[key(x) + 1] += 1;
    }
    for i in 0..n {
        count[i + 1] += count[i];
    }
    let mut out = items.clone();
    for x in items {
        let k = key(&x);
        out[count[k]] = x;
        count[k] += 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orders_and_dedupe() {
        // Edges: 0-2, 0-1, 2-0 again, 1-1.
        let g = Graph::from_edges(3, &[0, 0, 2, 1], &[2, 1, 0, 1], None, false).unwrap();
        let s = rows(&g, Direction::Out, NeighborOrder::Sorted, false, false);
        assert_eq!(s.neighbors(0), &[1, 2, 2]);
        assert_eq!(s.edges(0), &[1, 0, 2]);
        let i = rows(&g, Direction::Out, NeighborOrder::Input, false, true);
        assert_eq!(i.neighbors(0), &[2, 1]);
        assert_eq!(i.neighbors(1), &[0, 1]);
        let i = rows(&g, Direction::Out, NeighborOrder::Input, true, true);
        assert_eq!(i.neighbors(1), &[0]);
        let d = Graph::from_edges(3, &[2, 0, 1], &[0, 2, 0], None, true).unwrap();
        let r = rows(&d, Direction::In, NeighborOrder::Input, true, false);
        assert_eq!(r.neighbors(0), &[2, 1]);
        let b = rows(&d, Direction::Both, NeighborOrder::Sorted, true, true);
        assert_eq!(b.neighbors(0), &[1, 2]);
    }
}
