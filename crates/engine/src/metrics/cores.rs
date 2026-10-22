//! k-cores and coreness.
//!
//! The k-core of a graph (Seidman 1983, `seidman1983network`) is its largest
//! subgraph in which every node has degree at least k. The coreness (core
//! number) of a node is the largest k for which it belongs to the k-core.
//! Coreness is computed with the O(n + m) bucket algorithm of Batagelj and
//! Zaversnik (2003, `batagelj2003algorithm`): nodes are removed in order of
//! non-decreasing current degree, and each removal lowers the degree of the
//! node's remaining neighbors.
//!
//! Conventions (matching igraph's `coreness()`):
//! - Directed graphs take a mode: `All` counts in-arcs plus out-arcs (a
//!   mutual pair counts twice), `Out` builds cores on out-degree, and `In`
//!   on in-degree. Undirected graphs ignore the mode.
//! - Parallel edges each count, as in degree.
//! - A self-loop adds 2 to a node's degree in an undirected graph, and 1 to
//!   out-degree and 1 to in-degree in a directed graph. A loop stays inside
//!   every subgraph that contains its node, so it never stops counting.
//!   NetworkX's `core_number` rejects self-loops and multigraphs; on simple
//!   graphs without loops the two agree.
//! - Weights are ignored: cores are defined on binary ties.
//! - Isolated nodes have coreness 0; disconnected graphs need no special
//!   treatment.

use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;

/// Coreness of every node, in node order. See the module documentation.
pub fn coreness(g: &Graph, mode: DegreeMode) -> Vec<usize> {
    let n = g.node_count();
    let mode = if g.is_directed() {
        mode
    } else {
        DegreeMode::All
    };
    let deg: Vec<usize> = (0..n)
        .map(|v| {
            if !g.is_directed() {
                // Undirected rows store a loop once; it contributes 2.
                let row = g.out().neighbors(v);
                row.len() + row.iter().filter(|&&u| u as usize == v).count()
            } else {
                match mode {
                    DegreeMode::Out => g.out().degree(v),
                    DegreeMode::In => g.inc().degree(v),
                    DegreeMode::All => g.out().degree(v) + g.inc().degree(v),
                }
            }
        })
        .collect();
    // Removing v lowers the degree of every node whose counted arcs reach v:
    // v's in-neighbors for out-cores, its out-neighbors for in-cores.
    peel(deg, |v, buf| {
        buf.clear();
        if !g.is_directed() {
            buf.extend_from_slice(g.out().neighbors(v));
            return;
        }
        if mode != DegreeMode::Out {
            buf.extend_from_slice(g.out().neighbors(v));
        }
        if mode != DegreeMode::In {
            buf.extend_from_slice(g.inc().neighbors(v));
        }
    })
    .0
}

/// The bucket algorithm of Batagelj and Zaversnik. `deg` holds initial
/// degrees; `affected(v, buf)` fills `buf` with the nodes whose degree drops
/// by one per entry when `v` is removed. Returns the core numbers and the
/// removal order, which is a smallest-last (degeneracy) ordering (Matula
/// and Beck 1983, `matula1983smallest`).
pub(crate) fn peel<F>(mut deg: Vec<usize>, affected: F) -> (Vec<usize>, Vec<u32>)
where
    F: Fn(usize, &mut Vec<u32>),
{
    let n = deg.len();
    let max_deg = deg.iter().copied().max().unwrap_or(0);
    // bin[d] = position in `vert` of the first node with current degree d.
    let mut bin = vec![0usize; max_deg + 1];
    for &d in &deg {
        bin[d] += 1;
    }
    let mut start = 0;
    for b in bin.iter_mut() {
        let count = *b;
        *b = start;
        start += count;
    }
    let mut pos = vec![0usize; n];
    let mut vert = vec![0u32; n];
    {
        let mut next = bin.clone();
        for v in 0..n {
            pos[v] = next[deg[v]];
            vert[pos[v]] = v as u32;
            next[deg[v]] += 1;
        }
    }
    let mut buf = Vec::new();
    for i in 0..n {
        let v = vert[i] as usize;
        affected(v, &mut buf);
        for &u in &buf {
            let u = u as usize;
            if deg[u] > deg[v] {
                // Move u to the front of its bin, then shrink the bin by one.
                let du = deg[u];
                let pu = pos[u];
                let pw = bin[du];
                let w = vert[pw] as usize;
                if u != w {
                    pos[u] = pw;
                    vert[pu] = w as u32;
                    pos[w] = pu;
                    vert[pw] = u as u32;
                }
                bin[du] += 1;
                deg[u] -= 1;
            }
        }
    }
    (deg, vert)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn hand_checked_undirected() {
        // A triangle 0-1-2 with a pendant 3 on node 2 and an isolate 4.
        let g = und(5, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 2, 2, 1, 0]);
        // The NetworkX documentation example (Havel-Hakimi graph of
        // degrees [0, 1, 2, 2, 2, 2, 3]) has core numbers
        // {0: 1, 1: 2, 2: 2, 3: 2, 4: 1, 5: 2, 6: 0}; a K4 is a 3-core.
        let k4 = und(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]);
        assert_eq!(coreness(&k4, DegreeMode::All), vec![3; 4]);
    }

    #[test]
    fn loops_and_parallel_edges_count() {
        // A loop alone gives degree 2 and keeps it.
        let g = und(2, &[(0, 0), (0, 1)]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 1]);
        // A doubled edge is a 2-core.
        let g = und(2, &[(0, 1), (0, 1)]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 2]);
    }

    #[test]
    fn directed_modes() {
        // 0 -> 1, 0 -> 2, 1 -> 2, 2 -> 0
        let g = Graph::from_edges(3, &[0, 0, 1, 2], &[1, 2, 2, 0], None, true).unwrap();
        // Out-degrees 2, 1, 1. Removing 1 drops nothing below 1; the whole
        // graph is a 1-out-core.
        assert_eq!(coreness(&g, DegreeMode::Out), vec![1, 1, 1]);
        // In-degrees 1, 1, 2.
        assert_eq!(coreness(&g, DegreeMode::In), vec![1, 1, 1]);
        // Total degrees 3, 2, 3: removing 1 leaves 0 and 2 joined by two
        // arcs, a 2-core.
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 2, 2]);
    }

    #[test]
    fn empty_graph() {
        let g = und(0, &[]);
        assert!(coreness(&g, DegreeMode::All).is_empty());
        let g = und(3, &[]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![0, 0, 0]);
    }
}
