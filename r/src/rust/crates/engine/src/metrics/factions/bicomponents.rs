//! Bicomponents (biconnected components, or blocks).
//!
//! A bicomponent is a maximal set of nodes that stays connected when any
//! one of its nodes is removed; two bicomponents share at most one node,
//! a cut point (Wasserman and Faust 1994, `wasserman1994social`, section
//! 4.2). They are found with one
//! depth-first search computing low points (Hopcroft and Tarjan 1973,
//! `hopcroft1973algorithm`), written with an explicit stack so deep
//! networks cannot overflow the call stack: when the search returns from a
//! child w of v with low(w) >= disc(v), the nodes discovered since w,
//! together with v, form a bicomponent.
//!
//! Conventions (matching igraph's `biconnected_components()` and
//! NetworkX's `biconnected_components()`):
//! - The undirected simple reading of the network is used (see
//!   [`crate::metrics::factions`]).
//! - A bridge and its two endpoints form a bicomponent of two nodes, so
//!   every edge lies in exactly one bicomponent. Isolated nodes belong to
//!   none.
//! - Members are listed in increasing order and bicomponents in
//!   lexicographic order.

use super::{canonical_sets, simple_rows};
use crate::graph::Graph;

/// Node sets of the bicomponents of `g`. See the module documentation.
pub fn bicomponents(g: &Graph) -> Vec<Vec<u32>> {
    let rows = simple_rows(g);
    let n = rows.len();
    const UNSEEN: usize = usize::MAX;
    let mut disc = vec![UNSEEN; n];
    let mut low = vec![0usize; n];
    let mut stack: Vec<u32> = Vec::new();
    let mut out: Vec<Vec<u32>> = Vec::new();
    // Frames: node, its parent (or itself for a root), next row position.
    let mut frames: Vec<(usize, usize, usize)> = Vec::new();
    let mut timer = 0usize;
    for root in 0..n {
        if disc[root] != UNSEEN {
            continue;
        }
        disc[root] = timer;
        low[root] = timer;
        timer += 1;
        stack.push(root as u32);
        frames.push((root, root, 0));
        while let Some(frame) = frames.last_mut() {
            let (v, parent, pos) = *frame;
            if pos < rows[v].len() {
                frame.2 += 1;
                let w = rows[v][pos] as usize;
                if disc[w] == UNSEEN {
                    disc[w] = timer;
                    low[w] = timer;
                    timer += 1;
                    stack.push(w as u32);
                    frames.push((w, v, 0));
                } else if w != parent {
                    low[v] = low[v].min(disc[w]);
                }
                continue;
            }
            frames.pop();
            if v == parent {
                // The root: everything on the stack was assigned already.
                stack.pop();
                continue;
            }
            low[parent] = low[parent].min(low[v]);
            if low[v] >= disc[parent] {
                let mut block = vec![parent as u32];
                loop {
                    let x = stack.pop().expect("bicomponent stack underflow");
                    block.push(x);
                    if x as usize == v {
                        break;
                    }
                }
                out.push(block);
            }
        }
    }
    canonical_sets(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn blocks_of_a_bowtie_and_a_tail() {
        // Triangles {0,1,2} and {2,3,4} share cut point 2; 4-5 is a bridge;
        // 6 is isolated.
        let g = und(7, &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2), (4, 5)]);
        assert_eq!(
            bicomponents(&g),
            vec![vec![0, 1, 2], vec![2, 3, 4], vec![4, 5]]
        );
    }

    #[test]
    fn cycles_paths_and_loops() {
        let cycle = und(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
        assert_eq!(bicomponents(&cycle), vec![vec![0, 1, 2, 3, 4]]);
        let path = und(4, &[(0, 1), (1, 2), (2, 3), (1, 1), (1, 2)]);
        assert_eq!(
            bicomponents(&path),
            vec![vec![0, 1], vec![1, 2], vec![2, 3]]
        );
        assert!(bicomponents(&und(3, &[])).is_empty());
        // Direction is ignored.
        let d = Graph::from_edges(3, &[0, 1, 2], &[1, 2, 0], None, true).unwrap();
        assert_eq!(bicomponents(&d), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn deep_path_does_not_overflow() {
        let n = 100_000;
        let e: Vec<(i64, i64)> = (0..n as i64 - 1).map(|i| (i, i + 1)).collect();
        assert_eq!(bicomponents(&und(n, &e)).len(), n - 1);
    }
}
