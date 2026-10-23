//! Simple layouts: circle, star, grid, layered trees, and two rows for
//! two-mode networks.
//!
//! - Circle: node i at angle 2 pi i / n on the unit circle, counterclockwise
//!   from (1, 0), in node order (as igraph's `layout_in_circle()`).
//! - Star: the center at the origin and the other nodes on the unit circle
//!   in node order. The default center is the lowest-numbered node with the
//!   most distinct neighbors.
//! - Grid: rows of `width` nodes (by default the smallest width with
//!   width^2 >= n), filled left to right and top to bottom in node order,
//!   one unit apart.
//! - Tree: each component is drawn as a breadth-first spanning tree from a
//!   root, with depth d at y = -d. Leaves are one unit apart in depth-first
//!   order (children in node order), and each parent is centered over its
//!   first and last child. This satisfies the aesthetic rules of Reingold
//!   and Tilford (1981, `reingold1981tidier`): nodes at the same depth lie
//!   on a line, parents are centered over their children, subtrees do not
//!   overlap, and isomorphic subtrees are drawn identically; it is not as
//!   narrow as their algorithm, since every leaf gets its own column. The
//!   default root of a component is the middle node of the path between
//!   the two ends of a double sweep (a breadth-first search from the
//!   component's lowest-numbered node finds the farthest node u, and a
//!   search from u finds the farthest node v), which is a center of any
//!   tree. Components are packed as in [`crate::layout`].
//! - Two-mode: nodes of the first type in a top row and nodes of the second
//!   in a bottom row, each row spread over the same width and the rows half
//!   that width apart (at least one unit). The order within rows comes from
//!   the barycenter heuristic of Sugiyama, Tagawa, and Toda (1981,
//!   `sugiyama1981methods`): each sweep sorts one row by the mean relative
//!   position of each node's neighbors in the other row (nodes without such
//!   neighbors keep their own relative position; ties keep the previous
//!   order), alternating rows. Starting from node order, the ordering with
//!   the fewest edge crossings seen over `sweeps` sweeps is kept. Crossings
//!   are counted exactly as inversions of the edge list sorted by top
//!   position. Ties between nodes of the same type are drawn but do not
//!   affect the order.
//!
//! Direction and weights are ignored, parallel edges count once, and
//! self-loops are ignored. None of these layouts is random.

use super::{pack, parts, Layout, Part};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::weights::WeightKind;

/// Nodes on the unit circle in node order.
pub fn circle_layout(g: &Graph) -> Layout {
    let n = g.node_count();
    let mut l = Layout {
        x: vec![0.0; n],
        y: vec![0.0; n],
    };
    if n == 1 {
        return l;
    }
    for i in 0..n {
        let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        l.x[i] = a.cos();
        l.y[i] = a.sin();
    }
    l
}

/// One node at the center and the others on the unit circle in node order.
/// `center` defaults to the lowest-numbered node with the most distinct
/// neighbors.
pub fn star_layout(g: &Graph, center: Option<usize>) -> Result<Layout> {
    let n = g.node_count();
    let mut l = Layout {
        x: vec![0.0; n],
        y: vec![0.0; n],
    };
    if n == 0 {
        return Ok(l);
    }
    let center = match center {
        Some(c) if c >= n => {
            return Err(Error::InvalidArgument(format!(
                "center {c} is not a node index (0..{n})"
            )))
        }
        Some(c) => c,
        None => {
            let deg = distinct_degrees(g)?;
            (0..n).fold(0, |best, v| if deg[v] > deg[best] { v } else { best })
        }
    };
    let rim = n - 1;
    let mut k = 0;
    for v in 0..n {
        if v == center {
            continue;
        }
        let a = 2.0 * std::f64::consts::PI * k as f64 / rim as f64;
        l.x[v] = a.cos();
        l.y[v] = a.sin();
        k += 1;
    }
    Ok(l)
}

fn distinct_degrees(g: &Graph) -> Result<Vec<usize>> {
    let (ps, _) = parts(g, false, WeightKind::Distance)?;
    let mut deg = vec![0usize; g.node_count()];
    for p in &ps {
        for (i, &v) in p.nodes.iter().enumerate() {
            deg[v as usize] = p.degree(i);
        }
    }
    Ok(deg)
}

/// Nodes on a grid in node order, `width` per row (0 for the smallest
/// width whose square is at least n).
pub fn grid_layout(g: &Graph, width: usize) -> Layout {
    let n = g.node_count();
    let width = if width == 0 {
        let mut w = (n as f64).sqrt() as usize;
        while w * w < n {
            w += 1;
        }
        w.max(1)
    } else {
        width
    };
    let mut l = Layout {
        x: vec![0.0; n],
        y: vec![0.0; n],
    };
    for i in 0..n {
        l.x[i] = (i % width) as f64;
        l.y[i] = -((i / width) as f64);
    }
    l
}

/// Layered tree drawings of each component. `roots` may give at most one
/// root per component; the others use the default root.
pub fn tree_layout(g: &Graph, roots: &[usize]) -> Result<Layout> {
    let n = g.node_count();
    let (ps, _) = parts(g, false, WeightKind::Distance)?;
    let mut comp_of = vec![0usize; n];
    let mut local = vec![0usize; n];
    for (k, p) in ps.iter().enumerate() {
        for (i, &v) in p.nodes.iter().enumerate() {
            comp_of[v as usize] = k;
            local[v as usize] = i;
        }
    }
    let mut chosen: Vec<Option<usize>> = vec![None; ps.len()];
    for &r in roots {
        if r >= n {
            return Err(Error::InvalidArgument(format!(
                "root {r} is not a node index (0..{n})"
            )));
        }
        let k = comp_of[r];
        if chosen[k].is_some() {
            return Err(Error::InvalidArgument(format!(
                "more than one root given for the component of node {r}"
            )));
        }
        chosen[k] = Some(local[r]);
    }
    let mut pieces = Vec::with_capacity(ps.len());
    for (p, root) in ps.into_iter().zip(chosen) {
        let root = root.unwrap_or_else(|| center(&p));
        let xy = layered(&p, root);
        pieces.push((p.nodes, xy));
    }
    Ok(pack(n, &pieces, 1.0))
}

/// Breadth-first order and parents from `s`.
fn bfs(p: &Part, s: usize) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let c = p.len();
    let mut parent = vec![usize::MAX; c];
    let mut depth = vec![0usize; c];
    let mut order = Vec::with_capacity(c);
    parent[s] = s;
    order.push(s);
    let mut head = 0;
    while head < order.len() {
        let v = order[head];
        head += 1;
        for &t in &p.targets[p.arcs(v)] {
            let t = t as usize;
            if parent[t] == usize::MAX {
                parent[t] = v;
                depth[t] = depth[v] + 1;
                order.push(t);
            }
        }
    }
    (order, parent, depth)
}

/// The middle of the path between the ends of a double sweep.
fn center(p: &Part) -> usize {
    let farthest = |depth: &[usize]| -> usize {
        (0..depth.len()).fold(0, |b, v| if depth[v] > depth[b] { v } else { b })
    };
    let (_, _, d0) = bfs(p, 0);
    let u = farthest(&d0);
    let (_, parent, du) = bfs(p, u);
    let v = farthest(&du);
    let mut path = vec![v];
    let mut x = v;
    while x != u {
        x = parent[x];
        path.push(x);
    }
    path[(path.len() - 1) / 2]
}

fn layered(p: &Part, root: usize) -> Vec<[f64; 2]> {
    let c = p.len();
    let (order, parent, depth) = bfs(p, root);
    // Children in node order (BFS visits each row's targets in order).
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); c];
    for &v in order.iter().skip(1) {
        children[parent[v]].push(v);
    }
    // Leaves in depth-first preorder get consecutive columns.
    let mut x = vec![0.0; c];
    let mut next_leaf = 0.0;
    let mut stack = vec![root];
    while let Some(v) = stack.pop() {
        if children[v].is_empty() {
            x[v] = next_leaf;
            next_leaf += 1.0;
        }
        for &ch in children[v].iter().rev() {
            stack.push(ch);
        }
    }
    // Parents centered over their first and last child, deepest first.
    for &v in order.iter().rev() {
        if let (Some(&a), Some(&b)) = (children[v].first(), children[v].last()) {
            x[v] = (x[a] + x[b]) / 2.0;
        }
    }
    (0..c).map(|v| [x[v], -(depth[v] as f64)]).collect()
}

/// Two rows for a two-mode network: nodes with `types[v] == false` on top.
/// See the module documentation.
pub fn twomode_layout(g: &Graph, types: &[bool], sweeps: usize) -> Result<Layout> {
    let n = g.node_count();
    if types.len() != n {
        return Err(Error::LengthMismatch {
            what: "types",
            expected: n,
            found: types.len(),
        });
    }
    // Distinct ties between the two rows, as (top node, bottom node).
    let (ps, _) = parts(g, false, WeightKind::Distance)?;
    let mut ties: Vec<(usize, usize)> = Vec::new();
    for p in &ps {
        for i in 0..p.len() {
            let v = p.nodes[i] as usize;
            if types[v] {
                continue;
            }
            for &t in &p.targets[p.arcs(i)] {
                let w = p.nodes[t as usize] as usize;
                if types[w] {
                    ties.push((v, w));
                }
            }
        }
    }
    let top: Vec<usize> = (0..n).filter(|&v| !types[v]).collect();
    let bottom: Vec<usize> = (0..n).filter(|&v| types[v]).collect();
    let mut rank = vec![0usize; n];
    for row in [&top, &bottom] {
        for (r, &v) in row.iter().enumerate() {
            rank[v] = r;
        }
    }
    let mut nbr: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b) in &ties {
        nbr[a].push(b);
        nbr[b].push(a);
    }

    let mut rows = [top, bottom];
    let mut best = rank.clone();
    let mut best_cross = crossings(&ties, &rank);
    for sweep in 0..sweeps {
        if best_cross == 0 {
            break;
        }
        let which = if sweep % 2 == 0 { 1 } else { 0 };
        let other_len = rows[1 - which].len();
        let own_len = rows[which].len();
        let rel = |r: usize, len: usize| -> f64 {
            if len > 1 {
                r as f64 / (len - 1) as f64
            } else {
                0.5
            }
        };
        let key: Vec<f64> = rows[which]
            .iter()
            .map(|&v| {
                if nbr[v].is_empty() {
                    rel(rank[v], own_len)
                } else {
                    nbr[v].iter().map(|&u| rel(rank[u], other_len)).sum::<f64>()
                        / nbr[v].len() as f64
                }
            })
            .collect();
        let mut idx: Vec<usize> = (0..own_len).collect();
        idx.sort_by(|&a, &b| key[a].total_cmp(&key[b]));
        let reordered: Vec<usize> = idx.iter().map(|&i| rows[which][i]).collect();
        for (r, &v) in reordered.iter().enumerate() {
            rank[v] = r;
        }
        rows[which] = reordered;
        let c = crossings(&ties, &rank);
        if c < best_cross {
            best_cross = c;
            best.copy_from_slice(&rank);
        }
    }

    let (nt, nb) = (rows[0].len(), rows[1].len());
    let width = (nt.max(nb).max(2) - 1) as f64;
    let gap = (width / 2.0).max(1.0);
    let mut l = Layout {
        x: vec![0.0; n],
        y: vec![0.0; n],
    };
    for v in 0..n {
        let len = if types[v] { nb } else { nt };
        l.x[v] = if len > 1 {
            width * (best[v] as f64 / (len - 1) as f64 - 0.5)
        } else {
            0.0
        };
        l.y[v] = if types[v] { -gap / 2.0 } else { gap / 2.0 };
    }
    Ok(l)
}

/// Number of pairs of ties (a, b), (c, d) with a above c and b below d in
/// the orders given by `rank`: inversions of the bottom ranks after sorting
/// ties by (top rank, bottom rank), counted by merge sort.
fn crossings(ties: &[(usize, usize)], rank: &[usize]) -> u64 {
    let mut pairs: Vec<(usize, usize)> = ties.iter().map(|&(a, b)| (rank[a], rank[b])).collect();
    pairs.sort_unstable();
    let mut seq: Vec<usize> = pairs.into_iter().map(|p| p.1).collect();
    let mut buf = vec![0usize; seq.len()];
    inversions(&mut seq, &mut buf)
}

fn inversions(v: &mut [usize], buf: &mut [usize]) -> u64 {
    let n = v.len();
    if n < 2 {
        return 0;
    }
    let mid = n / 2;
    let mut count = {
        let (a, b) = v.split_at_mut(mid);
        let (ba, bb) = buf.split_at_mut(mid);
        inversions(a, ba) + inversions(b, bb)
    };
    let (mut i, mut j, mut k) = (0, mid, 0);
    while i < mid && j < n {
        if v[i] <= v[j] {
            buf[k] = v[i];
            i += 1;
        } else {
            buf[k] = v[j];
            j += 1;
            count += (mid - i) as u64;
        }
        k += 1;
    }
    buf[k..k + mid - i].copy_from_slice(&v[i..mid]);
    k += mid - i;
    buf[k..k + n - j].copy_from_slice(&v[j..n]);
    v.copy_from_slice(&buf[..n]);
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, edges: &[(i64, i64)]) -> Graph {
        let s: Vec<i64> = edges.iter().map(|e| e.0).collect();
        let t: Vec<i64> = edges.iter().map(|e| e.1).collect();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn circle_star_and_grid() {
        let g = graph(4, &[(0, 1), (1, 2), (1, 3)]);
        let c = circle_layout(&g);
        assert!((c.x[1]).abs() < 1e-15 && (c.y[1] - 1.0).abs() < 1e-15);
        let s = star_layout(&g, None).unwrap();
        assert_eq!((s.x[1], s.y[1]), (0.0, 0.0));
        assert!((s.x[0] - 1.0).abs() < 1e-15);
        assert!(star_layout(&g, Some(4)).is_err());
        let gr = grid_layout(&graph(5, &[]), 0);
        assert_eq!(gr.x, vec![0.0, 1.0, 2.0, 0.0, 1.0]);
        assert_eq!(gr.y, vec![0.0, 0.0, 0.0, -1.0, -1.0]);
    }

    #[test]
    fn trees_are_layered_and_centered() {
        // 0 - 1 - 2 - 3 - 4 with 2 also joined to 5: center is 2.
        let g = graph(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)]);
        let l = tree_layout(&g, &[]).unwrap();
        // Root on top, children one level down.
        assert!(l.y[2] > l.y[1] && (l.y[1] - l.y[3]).abs() < 1e-12);
        assert!((l.y[1] - l.y[5]).abs() < 1e-12 && (l.y[2] - l.y[1] - 1.0).abs() < 1e-12);
        // Parent centered over first and last child (1 and 5).
        assert!((l.x[2] - (l.x[1] + l.x[5]) / 2.0).abs() < 1e-12);
        let r = tree_layout(&g, &[0]).unwrap();
        assert!(r.y[0] > r.y[1]);
        assert!(tree_layout(&g, &[0, 4]).is_err());
    }

    #[test]
    fn twomode_rows_and_crossings() {
        // Top 0, 1; bottom 2, 3 with ties 0-3 and 1-2: one crossing in
        // node order, none after reordering.
        let g = graph(4, &[(0, 3), (1, 2)]);
        let types = [false, false, true, true];
        let ties = [(0, 3), (1, 2)];
        assert_eq!(crossings(&ties, &[0, 1, 0, 1]), 1);
        let l = twomode_layout(&g, &types, 8).unwrap();
        assert!(l.y[0] > 0.0 && l.y[2] < 0.0);
        let order_top = l.x[0] < l.x[1];
        let order_bottom = l.x[3] < l.x[2];
        assert_eq!(order_top, order_bottom);
        assert!(twomode_layout(&g, &types[..3], 8).is_err());
    }

    #[test]
    fn inversion_count_matches_brute_force() {
        let v = [3usize, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
        let mut brute = 0;
        for i in 0..v.len() {
            for j in i + 1..v.len() {
                if v[i] > v[j] {
                    brute += 1;
                }
            }
        }
        let mut w = v.to_vec();
        let mut buf = vec![0; v.len()];
        assert_eq!(inversions(&mut w, &mut buf), brute);
    }
}
