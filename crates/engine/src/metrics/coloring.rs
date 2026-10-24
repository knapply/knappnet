//! Greedy node colouring and random maximal independent sets.
//!
//! [`greedy_color`] visits the nodes in an order chosen by a strategy and
//! gives each the smallest colour (0, 1, ...) not used by its already
//! coloured neighbours, as NetworkX 3.7's `greedy_color`
//! (`hagberg2008exploring`). The strategies (Kosowski and Manuszewski
//! 2004, `kosowski2004classical`):
//!
//! - largest first: decreasing degree (Welsh and Powell 1967,
//!   `welsh1967upper`), ties in node order (NetworkX's stable sort);
//! - random sequential: a random order from the counter-based stream
//!   (`seed`), not NetworkX's Python generator;
//! - smallest last: the reverse of an order that repeatedly removes a node
//!   of smallest degree in what remains (Matula and Beck 1983,
//!   `matula1983smallest`), ties by smallest index;
//! - independent set: repeatedly a maximal independent set of the
//!   uncoloured nodes, built by taking a node of smallest degree among
//!   those left (ties by smallest index) and removing it and its
//!   neighbours; each set takes the next colour classes;
//! - connected sequential (breadth- or depth-first): component by
//!   component (in order of their smallest node), from the smallest node,
//!   in breadth- or depth-first order with neighbours in increasing index
//!   order;
//! - saturation largest first (DSATUR; Brélaz 1979, `brelaz1979new`):
//!   first a node of largest degree, then always the node with the most
//!   distinct colours among its coloured neighbours, ties by degree and
//!   then node order.
//!
//! Largest first and DSATUR reproduce NetworkX's colourings exactly.
//! For the other strategies NetworkX breaks ties by Python set order and
//! traverses neighbours in insertion order; the engine's choices are
//! documented above and give a valid colouring of the same kind.
//! NetworkX's `interchange` option is not provided.
//!
//! Conventions: degrees count parallel edges and an undirected self-loop
//! twice (NetworkX's `G.degree`), except that smallest last and
//! independent set use the degree in the simple graph without loops
//! (NetworkX fails on loops there); a self-loop never forbids a colour.
//! Directed graphs: as in NetworkX, a node avoids the colours of its
//! successors, degrees are in- plus out-degrees, and DSATUR counts the
//! colours of predecessors; the strategies that need components or
//! removals (smallest last, independent set, connected sequential) need an
//! undirected graph (NetworkX raises or misbehaves there).
//!
//! [`maximal_independent_set`] (NetworkX's `maximal_independent_set`):
//! starting from a given independent set (or one random node), add nodes
//! drawn uniformly from those not yet excluded until none is left. Adding
//! nodes in a uniformly random order and skipping excluded ones draws from
//! the same distribution, which is how it is computed. As in NetworkX, a
//! node with a self-loop may be added (only the starting set is checked
//! for independence); the graph must be undirected (direction is ignored).

use std::collections::BinaryHeap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::rng::Stream;

/// Node order for [`greedy_color`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColoringStrategy {
    /// Decreasing degree.
    LargestFirst,
    /// A seeded random order.
    RandomSequential,
    /// Smallest-last (degeneracy) order.
    SmallestLast,
    /// Successive maximal independent sets.
    IndependentSet,
    /// Breadth-first within components.
    ConnectedSequentialBfs,
    /// Depth-first within components.
    ConnectedSequentialDfs,
    /// DSATUR.
    SaturationLargestFirst,
}

impl ColoringStrategy {
    /// Decodes the C interface codes: 0 = largest first, 1 = random
    /// sequential, 2 = smallest last, 3 = independent set, 4 = connected
    /// sequential breadth-first, 5 = depth-first, 6 = saturation largest
    /// first (DSATUR).
    pub fn from_code(code: i32) -> Result<ColoringStrategy> {
        use ColoringStrategy::*;
        Ok(match code {
            0 => LargestFirst,
            1 => RandomSequential,
            2 => SmallestLast,
            3 => IndependentSet,
            4 => ConnectedSequentialBfs,
            5 => ConnectedSequentialDfs,
            6 => SaturationLargestFirst,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown colouring strategy code {code}"
                )))
            }
        })
    }
}

fn nx_degree(g: &Graph, v: usize) -> usize {
    if g.is_directed() {
        g.out().degree(v) + g.inc().degree(v)
    } else {
        let row = g.out().neighbors(v);
        row.len() + row.iter().filter(|&&u| u as usize == v).count()
    }
}

/// Colours in the given order: the smallest colour unused by coloured
/// neighbours (successors in a directed graph).
fn colour_in_order(g: &Graph, order: &[u32]) -> Vec<u32> {
    let n = g.node_count();
    let mut colour = vec![u32::MAX; n];
    let mut mark = vec![usize::MAX; n + 1];
    for (step, &u) in order.iter().enumerate() {
        let u = u as usize;
        for &v in g.out().neighbors(u) {
            let c = colour[v as usize];
            if v as usize != u && c != u32::MAX && (c as usize) <= n {
                mark[c as usize] = step;
            }
        }
        let c = (0..).find(|&c| mark[c] != step).expect("a free colour");
        colour[u] = c as u32;
    }
    colour
}

/// A greedy colouring (colour of every node, from 0). `seed` drives the
/// random sequential strategy. See the module documentation.
pub fn greedy_color(g: &Graph, strategy: ColoringStrategy, seed: u64) -> Result<Vec<u32>> {
    use ColoringStrategy::*;
    let n = g.node_count();
    let needs_undirected = matches!(
        strategy,
        SmallestLast | IndependentSet | ConnectedSequentialBfs | ConnectedSequentialDfs
    );
    if needs_undirected && g.is_directed() {
        return Err(Error::InvalidArgument(
            "this colouring strategy needs an undirected graph".into(),
        ));
    }
    let order: Vec<u32> = match strategy {
        LargestFirst => {
            let mut o: Vec<u32> = (0..n as u32).collect();
            o.sort_by_key(|&v| std::cmp::Reverse(nx_degree(g, v as usize)));
            o
        }
        RandomSequential => {
            let mut p = vec![0usize; n];
            Stream::new(seed, 0).permutation(&mut p);
            p.into_iter().map(|v| v as u32).collect()
        }
        SmallestLast => smallest_last(g),
        IndependentSet => independent_sets(g),
        ConnectedSequentialBfs | ConnectedSequentialDfs => {
            connected_sequential(g, strategy == ConnectedSequentialDfs)
        }
        SaturationLargestFirst => return Ok(dsatur(g)),
    };
    Ok(colour_in_order(g, &order))
}

fn smallest_last(g: &Graph) -> Vec<u32> {
    let s = g.simple_undirected(false);
    let n = s.node_count();
    let mut deg: Vec<usize> = (0..n).map(|v| s.degree(v)).collect();
    let mut heap: BinaryHeap<std::cmp::Reverse<(usize, u32)>> = (0..n)
        .map(|v| std::cmp::Reverse((deg[v], v as u32)))
        .collect();
    let mut removed = vec![false; n];
    let mut order = Vec::with_capacity(n);
    while let Some(std::cmp::Reverse((d, v))) = heap.pop() {
        let v = v as usize;
        if removed[v] || d != deg[v] {
            continue;
        }
        removed[v] = true;
        order.push(v as u32);
        for &w in s.neighbors(v) {
            let w = w as usize;
            if !removed[w] {
                deg[w] -= 1;
                heap.push(std::cmp::Reverse((deg[w], w as u32)));
            }
        }
    }
    order.reverse();
    order
}

fn independent_sets(g: &Graph) -> Vec<u32> {
    let s = g.simple_undirected(false);
    let n = s.node_count();
    let mut coloured = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut left = n;
    while left > 0 {
        // One maximal independent set of the uncoloured nodes.
        let mut alive: Vec<bool> = coloured.iter().map(|&c| !c).collect();
        let mut deg: Vec<usize> = (0..n)
            .map(|v| {
                if alive[v] {
                    s.neighbors(v)
                        .iter()
                        .filter(|&&w| alive[w as usize])
                        .count()
                } else {
                    0
                }
            })
            .collect();
        let mut heap: BinaryHeap<std::cmp::Reverse<(usize, u32)>> = (0..n)
            .filter(|&v| alive[v])
            .map(|v| std::cmp::Reverse((deg[v], v as u32)))
            .collect();
        while let Some(std::cmp::Reverse((d, v))) = heap.pop() {
            let v = v as usize;
            if !alive[v] || d != deg[v] {
                continue;
            }
            order.push(v as u32);
            coloured[v] = true;
            left -= 1;
            // Remove v and its neighbours from this round's graph.
            let mut gone = vec![v];
            gone.extend(
                s.neighbors(v)
                    .iter()
                    .map(|&w| w as usize)
                    .filter(|&w| alive[w]),
            );
            for &x in &gone {
                alive[x] = false;
            }
            for &x in &gone {
                for &y in s.neighbors(x) {
                    let y = y as usize;
                    if alive[y] {
                        deg[y] -= 1;
                        heap.push(std::cmp::Reverse((deg[y], y as u32)));
                    }
                }
            }
        }
    }
    order
}

fn connected_sequential(g: &Graph, depth_first: bool) -> Vec<u32> {
    let s = g.simple_undirected(false);
    let n = s.node_count();
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    for root in 0..n {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        order.push(root as u32);
        if depth_first {
            // Preorder of NetworkX's dfs_edges: (node, next neighbour).
            let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
            while let Some(top) = stack.last_mut() {
                let (v, i) = *top;
                let row = s.neighbors(v);
                if i == row.len() {
                    stack.pop();
                    continue;
                }
                top.1 += 1;
                let w = row[i] as usize;
                if !seen[w] {
                    seen[w] = true;
                    order.push(w as u32);
                    stack.push((w, 0));
                }
            }
        } else {
            let mut head = order.len() - 1;
            while head < order.len() {
                let v = order[head] as usize;
                head += 1;
                for &w in s.neighbors(v) {
                    if !seen[w as usize] {
                        seen[w as usize] = true;
                        order.push(w);
                    }
                }
            }
        }
    }
    order
}

fn dsatur(g: &Graph) -> Vec<u32> {
    let n = g.node_count();
    if n == 0 {
        return Vec::new();
    }
    let deg: Vec<usize> = (0..n).map(|v| nx_degree(g, v)).collect();
    // Colours seen by each node: from predecessors in a directed graph
    // (NetworkX adds a node's colour to its successors' sets).
    let mut seen: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut colour = vec![u32::MAX; n];
    let mut heap: BinaryHeap<(usize, usize, std::cmp::Reverse<u32>)> = BinaryHeap::new();
    let first = (0..n)
        .max_by(|&a, &b| deg[a].cmp(&deg[b]).then(b.cmp(&a)))
        .expect("a node");
    heap.push((usize::MAX, deg[first], std::cmp::Reverse(first as u32)));
    for (v, &d) in deg.iter().enumerate() {
        if v != first {
            heap.push((0, d, std::cmp::Reverse(v as u32)));
        }
    }
    let mut mark = vec![usize::MAX; n + 1];
    let mut step = 0;
    while let Some((sat, _, std::cmp::Reverse(u))) = heap.pop() {
        let u = u as usize;
        if colour[u] != u32::MAX || (sat != usize::MAX && sat != seen[u].len()) {
            continue;
        }
        for &v in g.out().neighbors(u) {
            let c = colour[v as usize];
            if v as usize != u && c != u32::MAX {
                mark[c as usize] = step;
            }
        }
        let c = (0..).find(|&c| mark[c] != step).expect("a free colour") as u32;
        step += 1;
        colour[u] = c;
        for &v in g.out().neighbors(u) {
            let v = v as usize;
            if colour[v] == u32::MAX {
                if let Err(pos) = seen[v].binary_search(&c) {
                    seen[v].insert(pos, c);
                    heap.push((seen[v].len(), deg[v], std::cmp::Reverse(v as u32)));
                }
            }
        }
    }
    colour
}

/// A random maximal independent set containing `start` (which must be
/// independent), or containing one random node when `start` is empty;
/// nodes increasing. See the module documentation.
pub fn maximal_independent_set(g: &Graph, start: &[u32], seed: u64) -> Result<Vec<u32>> {
    let n = g.node_count();
    let mut in_set = vec![false; n];
    let mut excluded = vec![false; n];
    let neighbours = |v: usize| {
        let a = g.out().neighbors(v).iter();
        let b: &[u32] = if g.is_directed() {
            g.inc().neighbors(v)
        } else {
            &[]
        };
        a.chain(b.iter())
            .map(|&w| w as usize)
            .collect::<Vec<usize>>()
    };
    for &v in start {
        if v as usize >= n {
            return Err(Error::InvalidArgument(format!(
                "node {v} of the starting set is outside 0..{n}"
            )));
        }
        in_set[v as usize] = true;
    }
    for &v in start {
        for w in neighbours(v as usize) {
            if in_set[w] {
                return Err(Error::InvalidArgument(
                    "the starting nodes are not an independent set".into(),
                ));
            }
            excluded[w] = true;
        }
        excluded[v as usize] = true;
    }
    let mut order = vec![0usize; n];
    Stream::new(seed, 0).permutation(&mut order);
    for v in order {
        if excluded[v] {
            continue;
        }
        in_set[v] = true;
        excluded[v] = true;
        for w in neighbours(v) {
            excluded[w] = true;
        }
    }
    Ok((0..n).filter(|&v| in_set[v]).map(|v| v as u32).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colourings_of_small_graphs() {
        // A 5-cycle needs 3 colours; every strategy finds a proper one.
        let g = Graph::from_edges(5, &[0, 1, 2, 3, 4], &[1, 2, 3, 4, 0], None, false).unwrap();
        for code in 0..7 {
            let c = greedy_color(&g, ColoringStrategy::from_code(code).unwrap(), 3).unwrap();
            for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
                assert_ne!(c[a as usize], c[b as usize]);
            }
            assert_eq!(*c.iter().max().unwrap(), 2);
        }
        // Largest first on a star: the centre first.
        let s = Graph::from_edges(4, &[1, 1, 1], &[0, 2, 3], None, false).unwrap();
        assert_eq!(
            greedy_color(&s, ColoringStrategy::LargestFirst, 0).unwrap(),
            vec![1, 0, 1, 1]
        );
        let m = maximal_independent_set(&s, &[0], 9).unwrap();
        assert_eq!(m, vec![0, 2, 3]);
        assert!(maximal_independent_set(&s, &[0, 1], 9).is_err());
    }
}
