//! Strongly connected components in traversal order, condensation,
//! attracting components, semiconnectedness, and local bridges.
//!
//! Sources:
//! - Strongly connected components by Tarjan's (1972, `tarjan1972depth`)
//!   depth-first algorithm with an explicit stack. Components are numbered
//!   in the order the algorithm completes them, which is a reverse
//!   topological order of the condensation; NetworkX's
//!   `strongly_connected_components` (`hagberg2008exploring`), a variant
//!   due to Pearce (2016, `pearce2016space`) that completes components in
//!   the same order, numbers `condensation` nodes this way. With the input
//!   neighbour order ([`crate::order`]) the numbering is NetworkX's.
//! - The condensation contracts each strongly connected component to one
//!   node (Harary, Norman, and Cartwright 1965, `harary1965structural`).
//!   An attracting component is one the condensation has no arc out of: a
//!   walk that enters it never leaves (Markov chains: a closed class;
//!   Kemeny and Snell 1960, `kemeny1960finite`).
//! - A directed network is semiconnected (unilaterally connected) when for
//!   every pair one node reaches the other, which holds exactly when its
//!   condensation has a path through all its nodes, that is when
//!   consecutive components in a topological order are joined by an arc
//!   (Harary, Norman, and Cartwright 1965).
//! - A local bridge is an edge whose endpoints have no common neighbour,
//!   so that it is not part of any triangle; its span is the distance
//!   between its endpoints once it is removed (Granovetter 1973,
//!   `granovetter1973strength`; Easley and Kleinberg 2010,
//!   `easley2010networks`, section 3.2).
//!
//! Conventions: condensation, attracting components, and semiconnectedness
//! are for directed networks (NetworkX refuses undirected graphs); parallel
//! edges and loops do not matter to them. Local bridges are for undirected
//! networks and use the simple graph (NetworkX refuses directed graphs and
//! multigraphs); a self-loop is never a local bridge, parallel edges share
//! their pair's answer, spans count edges (weights are not used), and a
//! span is infinite when no other path joins the endpoints.

use super::Simple;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::order::{rows, Direction, NeighborOrder};

/// Strongly connected component of every node, numbered in the order
/// Tarjan's algorithm completes them (roots tried in node order,
/// neighbours in `order`). For an undirected network, the connected
/// components in the same numbering.
pub fn strong_components(g: &Graph, order: NeighborOrder) -> Vec<u32> {
    let n = g.node_count();
    let adj = rows(g, Direction::Out, order, true, false);
    const UNSEEN: u32 = u32::MAX;
    let mut index = vec![UNSEEN; n];
    let mut low = vec![0u32; n];
    let mut on_stack = vec![false; n];
    let mut comp = vec![UNSEEN; n];
    let mut stack: Vec<u32> = Vec::new();
    let mut frames: Vec<(u32, usize)> = Vec::new();
    let mut counter = 0u32;
    let mut ncomp = 0u32;
    for root in 0..n {
        if index[root] != UNSEEN {
            continue;
        }
        index[root] = counter;
        low[root] = counter;
        counter += 1;
        stack.push(root as u32);
        on_stack[root] = true;
        frames.push((root as u32, 0));
        while let Some(frame) = frames.last_mut() {
            let v = frame.0 as usize;
            let nbrs = adj.neighbors(v);
            if frame.1 < nbrs.len() {
                let w = nbrs[frame.1] as usize;
                frame.1 += 1;
                if index[w] == UNSEEN {
                    index[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push(w as u32);
                    on_stack[w] = true;
                    frames.push((w as u32, 0));
                } else if on_stack[w] {
                    low[v] = low[v].min(index[w]);
                }
                continue;
            }
            frames.pop();
            if let Some(&(p, _)) = frames.last() {
                let p = p as usize;
                low[p] = low[p].min(low[v]);
            }
            if low[v] == index[v] {
                loop {
                    let w = stack.pop().expect("component on the stack") as usize;
                    on_stack[w] = false;
                    comp[w] = ncomp;
                    if w == v {
                        break;
                    }
                }
                ncomp += 1;
            }
        }
    }
    comp
}

fn directed_only(g: &Graph, what: &str) -> Result<()> {
    if g.is_directed() {
        Ok(())
    } else {
        Err(Error::InvalidArgument(format!(
            "{what} is defined for directed networks"
        )))
    }
}

/// Component of every node in the condensation, numbered as
/// [`strong_components`] does (NetworkX's `condensation` mapping when
/// `order` is the input order). The condensation's arcs are the distinct
/// pairs (label of u, label of v) over the arcs u -> v between different
/// components.
pub fn condensation(g: &Graph, order: NeighborOrder) -> Result<Vec<u32>> {
    directed_only(g, "the condensation")?;
    Ok(strong_components(g, order))
}

/// Whether each strongly connected component (numbered as in
/// [`condensation`]) has an arc leaving it, as a vector indexed by label.
fn has_exit(g: &Graph, comp: &[u32]) -> Vec<bool> {
    let k = comp.iter().map(|&c| c as usize + 1).max().unwrap_or(0);
    let mut exit = vec![false; k];
    for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
        let (ca, cb) = (comp[a as usize], comp[b as usize]);
        if ca != cb {
            exit[ca as usize] = true;
        }
    }
    exit
}

/// Label of each node's attracting component (in the numbering of
/// [`condensation`]), or `None` for a node outside every attracting
/// component.
pub fn attracting_components(g: &Graph, order: NeighborOrder) -> Result<Vec<Option<u32>>> {
    directed_only(g, "an attracting component")?;
    let comp = strong_components(g, order);
    let exit = has_exit(g, &comp);
    Ok(comp
        .iter()
        .map(|&c| (!exit[c as usize]).then_some(c))
        .collect())
}

/// Whether a directed network is semiconnected. A network with no nodes is
/// an error (NetworkX raises `NetworkXPointlessConcept`).
pub fn is_semiconnected(g: &Graph) -> Result<bool> {
    directed_only(g, "semiconnectedness")?;
    let n = g.node_count();
    if n == 0 {
        return Err(Error::InvalidArgument(
            "connectivity is not defined for a network with no nodes".into(),
        ));
    }
    if !Simple::new(g).is_connected() {
        return Ok(false);
    }
    let comp = strong_components(g, NeighborOrder::Sorted);
    let k = comp.iter().map(|&c| c as usize + 1).max().unwrap_or(0);
    // Labels are a reverse topological order: component i + 1 must have an
    // arc into component i.
    let mut joined = vec![false; k];
    for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
        let (ca, cb) = (comp[a as usize], comp[b as usize]);
        if ca == cb + 1 {
            joined[cb as usize] = true;
        }
    }
    Ok(joined[..k - 1].iter().all(|&j| j))
}

/// Local bridges of an undirected network, per edge in input order.
#[derive(Debug, Clone, PartialEq)]
pub struct LocalBridges {
    /// Whether each edge is a local bridge.
    pub is_local_bridge: Vec<bool>,
    /// Span of each local bridge (infinite when no other path joins its
    /// endpoints); NaN for other edges.
    pub span: Vec<f64>,
}

/// Local bridges and their spans. See the module documentation.
pub fn local_bridges(g: &Graph) -> Result<LocalBridges> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(
            "local bridges are defined for undirected networks".into(),
        ));
    }
    let simple = Simple::new(g);
    let n = simple.n;
    let m = g.edge_count();
    let mut is_local_bridge = vec![false; m];
    let mut span = vec![f64::NAN; m];
    let mut dist = vec![u32::MAX; n];
    let mut queue: Vec<usize> = Vec::new();
    let mut cache: std::collections::HashMap<(u32, u32), f64> = std::collections::HashMap::new();
    for e in 0..m {
        let (a, b) = (g.edge_sources()[e], g.edge_targets()[e]);
        if a == b {
            continue;
        }
        let (u, v) = (a.min(b), a.max(b));
        let (ru, rv) = (&simple.out[u as usize], &simple.out[v as usize]);
        let (mut i, mut j) = (0, 0);
        let mut common = false;
        while i < ru.len() && j < rv.len() {
            match ru[i].cmp(&rv[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    common = true;
                    break;
                }
            }
        }
        if common {
            continue;
        }
        is_local_bridge[e] = true;
        if let Some(&s) = cache.get(&(u, v)) {
            span[e] = s;
            continue;
        }
        // Breadth-first search from u that does not use the edge u-v.
        let (u, v) = (u as usize, v as usize);
        queue.clear();
        dist[u] = 0;
        queue.push(u);
        let mut found = f64::INFINITY;
        let mut head = 0;
        'bfs: while head < queue.len() {
            let x = queue[head];
            head += 1;
            for &y in &simple.out[x] {
                let y = y as usize;
                if x == u && y == v {
                    continue;
                }
                if dist[y] == u32::MAX {
                    dist[y] = dist[x] + 1;
                    if y == v {
                        found = dist[y] as f64;
                        break 'bfs;
                    }
                    queue.push(y);
                }
            }
        }
        for &x in &queue {
            dist[x] = u32::MAX;
        }
        dist[v] = u32::MAX;
        cache.insert((u as u32, v as u32), found);
        span[e] = found;
    }
    Ok(LocalBridges {
        is_local_bridge,
        span,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn directed_structure() {
        // 0 <-> 1 -> 2 <-> 3, 4 isolated.
        let g = graph(5, &[(0, 1), (1, 0), (1, 2), (2, 3), (3, 2)], true);
        let c = condensation(&g, NeighborOrder::Sorted).unwrap();
        assert_eq!(c, vec![1, 1, 0, 0, 2]);
        let a = attracting_components(&g, NeighborOrder::Sorted).unwrap();
        assert_eq!(a, vec![None, None, Some(0), Some(0), Some(2)]);
        assert!(!is_semiconnected(&g).unwrap());
        let h = graph(4, &[(0, 1), (1, 0), (1, 2), (2, 3), (3, 2)], true);
        assert!(is_semiconnected(&h).unwrap());
        assert!(condensation(&graph(2, &[(0, 1)], false), NeighborOrder::Sorted).is_err());
    }

    #[test]
    fn local_bridge_spans() {
        // Square 0-1-2-3-0 with a triangle 0-1-4 and a pendant 2-5.
        let g = graph(
            6,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                (0, 4),
                (1, 4),
                (2, 5),
                (5, 5),
            ],
            false,
        );
        let b = local_bridges(&g).unwrap();
        assert_eq!(
            b.is_local_bridge,
            vec![false, true, true, true, false, false, true, false]
        );
        assert_eq!(b.span[1], 3.0);
        assert_eq!(b.span[6], f64::INFINITY);
        assert!(b.span[0].is_nan());
    }
}
