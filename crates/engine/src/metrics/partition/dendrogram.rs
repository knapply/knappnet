//! Dendrograms of hierarchical community methods: conversion to a common
//! representation, cutting at a number of groups, and cutting at the level
//! of highest modularity.
//!
//! The representation is that of igraph's `merges` matrix (and, shifted by
//! one and with leaves negated, of R's `hclust`): the n leaves are the
//! nodes 0..n, and merge i joins two existing clusters into the new
//! cluster n + i. A dendrogram of a disconnected network from a method that
//! only merges connected groups is incomplete: it has n - c merges for c
//! components, and cutting it into fewer than c groups is an error, as in
//! igraph's `cut_at()`.
//!
//! - [`merges_from_node_pairs`] converts a merge list whose groups are
//!   named by any of their nodes (as [`super::fast_greedy()`] reports them).
//! - [`merges_from_edge_removals`] builds the dendrogram of a divisive
//!   method from its edge removals (as [`super::girvan_newman()`] reports
//!   them): the edges are put back in reverse order of removal, and every
//!   edge that joins two components records a merge, which is how igraph
//!   derives the merges of `cluster_edge_betweenness()`.
//! - [`cut_dendrogram`] performs the first `steps` merges (igraph's
//!   `community_to_membership()`), and [`cut_at_groups`] cuts at k groups.
//! - [`dendrogram_modularity`] is the modularity of
//!   [`crate::metrics::community::modularity`] after 0, 1, ... merges, and
//!   [`cut_at_max_modularity`] cuts at the first level of highest
//!   modularity. Undirected networks only, as for the hierarchical methods.
//!
//! Groups are labelled 0, 1, ... in order of their first node.

use super::{additive_weights, relabel, require_undirected};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::union_find::UnionFind;

/// Checks that `merges` is a valid (possibly incomplete) dendrogram on n
/// leaves: merge i names two distinct clusters that exist before it
/// (below n + i) and that no earlier merge consumed.
pub fn check_merges(n: usize, merges: &[(usize, usize)]) -> Result<()> {
    if n > 0 && merges.len() > n - 1 {
        return Err(Error::InvalidArgument(format!(
            "a dendrogram on {n} leaves has at most {} merges, not {}",
            n - 1,
            merges.len()
        )));
    }
    if n == 0 && !merges.is_empty() {
        return Err(Error::InvalidArgument(
            "a dendrogram without leaves has no merges".into(),
        ));
    }
    let mut used = vec![false; n + merges.len()];
    for (i, &(a, b)) in merges.iter().enumerate() {
        for c in [a, b] {
            if c >= n + i {
                return Err(Error::InvalidArgument(format!(
                    "merge {} names cluster {c}, which does not exist yet (clusters \
                     0..{} exist before it)",
                    i + 1,
                    n + i
                )));
            }
            if used[c] {
                return Err(Error::InvalidArgument(format!(
                    "merge {} names cluster {c}, which an earlier merge already joined",
                    i + 1
                )));
            }
            used[c] = true;
        }
        if a == b {
            return Err(Error::InvalidArgument(format!(
                "merge {} joins cluster {a} with itself",
                i + 1
            )));
        }
    }
    Ok(())
}

/// Converts mergers of groups named by any of their member nodes into
/// igraph's numbering.
pub fn merges_from_node_pairs(n: usize, pairs: &[(u32, u32)]) -> Result<Vec<(usize, usize)>> {
    let mut sets = UnionFind::new(n);
    // Cluster id of each union-find root.
    let mut cluster: Vec<usize> = (0..n).collect();
    let mut out = Vec::with_capacity(pairs.len());
    for (i, &(a, b)) in pairs.iter().enumerate() {
        let (a, b) = (a as usize, b as usize);
        if a >= n || b >= n {
            return Err(Error::InvalidArgument(format!(
                "merge {} names node {}, but there are {n} nodes",
                i + 1,
                a.max(b)
            )));
        }
        let (ra, rb) = (sets.find(a), sets.find(b));
        if ra == rb {
            return Err(Error::InvalidArgument(format!(
                "merge {} joins nodes {a} and {b}, which are already in one group",
                i + 1
            )));
        }
        out.push((cluster[ra], cluster[rb]));
        sets.union(ra, rb);
        let r = sets.find(ra);
        cluster[r] = n + out.len() - 1;
    }
    Ok(out)
}

/// The dendrogram of a divisive method from its edge removals, earliest
/// removal first: edges are restored in reverse order and every edge that
/// joins two components is a merge.
pub fn merges_from_edge_removals(n: usize, removed: &[(u32, u32)]) -> Result<Vec<(usize, usize)>> {
    let mut sets = UnionFind::new(n);
    let mut pairs = Vec::new();
    for &(a, b) in removed.iter().rev() {
        if a as usize >= n || b as usize >= n {
            return Err(Error::InvalidArgument(format!(
                "an edge removal names node {}, but there are {n} nodes",
                a.max(b)
            )));
        }
        if sets.union(a as usize, b as usize).is_some() {
            pairs.push((a, b));
        }
    }
    merges_from_node_pairs(n, &pairs)
}

/// Membership after the first `steps` merges (labels 0, 1, ... in order of
/// first node).
pub fn cut_dendrogram(n: usize, merges: &[(usize, usize)], steps: usize) -> Result<Vec<usize>> {
    check_merges(n, merges)?;
    if steps > merges.len() {
        return Err(Error::InvalidArgument(format!(
            "cannot perform {steps} merges: the dendrogram has {}",
            merges.len()
        )));
    }
    let mut sets = UnionFind::new(n);
    // A node representing each cluster id.
    let mut rep: Vec<usize> = (0..n).collect();
    for &(a, b) in &merges[..steps] {
        let (ra, rb) = (rep[a], rep[b]);
        sets.union(ra, rb);
        rep.push(ra);
    }
    let roots: Vec<usize> = (0..n).map(|v| sets.find(v)).collect();
    Ok(relabel(&roots))
}

/// Membership with exactly `k` groups (1 <= k <= n), cutting the
/// dendrogram after n - k merges.
pub fn cut_at_groups(n: usize, merges: &[(usize, usize)], k: usize) -> Result<Vec<usize>> {
    check_merges(n, merges)?;
    if k == 0 || k > n.max(1) {
        return Err(Error::InvalidArgument(format!(
            "the number of groups must be between 1 and {}, not {k}",
            n.max(1)
        )));
    }
    let steps = n - k.min(n);
    if steps > merges.len() {
        return Err(Error::InvalidArgument(format!(
            "the dendrogram is incomplete (the network has disconnected parts): it \
             has {} merges, so it cannot be cut into fewer than {} groups",
            merges.len(),
            n - merges.len()
        )));
    }
    cut_dendrogram(n, merges, steps)
}

/// Modularity before any merge (entry 0) and after each merge, on the
/// undirected network `g` whose nodes are the leaves. NaN throughout when
/// the network has no edges (or total weight zero).
pub fn dendrogram_modularity(
    g: &Graph,
    merges: &[(usize, usize)],
    resolution: f64,
    weighted: bool,
) -> Result<Vec<f64>> {
    require_undirected(g, "Modularity of a dendrogram")?;
    if !(resolution.is_finite() && resolution >= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "resolution must be finite and non-negative, not {resolution}"
        )));
    }
    let n = g.node_count();
    check_merges(n, merges)?;
    let w = additive_weights(g, weighted)?;
    // Adjacency rows of (neighbor, weight), loops kept separately.
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut inside = vec![0.0f64; n + merges.len()];
    let mut deg = vec![0.0f64; n + merges.len()];
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    for e in 0..g.edge_count() {
        let x = w.map_or(1.0, |w| w[e]);
        let (u, v) = (src[e] as usize, tgt[e] as usize);
        deg[u] += x;
        deg[v] += x;
        if u == v {
            inside[u] += x;
        } else {
            rows[u].push((v, x));
            rows[v].push((u, x));
        }
    }
    let m: f64 = deg[..n].iter().sum::<f64>() / 2.0;
    if m == 0.0 {
        return Ok(vec![f64::NAN; merges.len() + 1]);
    }
    let term = |inside: f64, deg: f64| inside / m - resolution * (deg / (2.0 * m)).powi(2);
    let mut q: f64 = (0..n).map(|c| term(inside[c], deg[c])).sum();
    let mut path = Vec::with_capacity(merges.len() + 1);
    path.push(q);
    // Members of each live cluster, and the cluster of each node; the
    // smaller cluster's members move into the larger's list.
    let mut members: Vec<Vec<usize>> = (0..n).map(|v| vec![v]).collect();
    members.resize(n + merges.len(), Vec::new());
    let mut label: Vec<usize> = (0..n).collect();
    for (i, &(a, b)) in merges.iter().enumerate() {
        let c = n + i;
        let (small, large) = if members[a].len() <= members[b].len() {
            (a, b)
        } else {
            (b, a)
        };
        let mut cross = 0.0;
        for &v in &members[small] {
            for &(u, x) in &rows[v] {
                if label[u] == large {
                    cross += x;
                }
            }
        }
        q -= term(inside[a], deg[a]) + term(inside[b], deg[b]);
        inside[c] = inside[a] + inside[b] + cross;
        deg[c] = deg[a] + deg[b];
        q += term(inside[c], deg[c]);
        let mut list = std::mem::take(&mut members[large]);
        let moved = std::mem::take(&mut members[small]);
        list.extend_from_slice(&moved);
        for &v in &list {
            label[v] = c;
        }
        members[c] = list;
        path.push(q);
    }
    Ok(path)
}

/// Cuts at the first level of highest modularity. Returns the membership
/// and its modularity (NaN, with one group per node, when the network has
/// no edges).
pub fn cut_at_max_modularity(
    g: &Graph,
    merges: &[(usize, usize)],
    resolution: f64,
    weighted: bool,
) -> Result<(Vec<usize>, f64)> {
    let path = dendrogram_modularity(g, merges, resolution, weighted)?;
    let mut best = 0;
    for (i, &q) in path.iter().enumerate() {
        if q > path[best] {
            best = i;
        }
    }
    let membership = cut_dendrogram(g.node_count(), merges, best)?;
    Ok((membership, path[best]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::community::modularity;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn conversion_and_cuts() {
        // Merge {0,1}, {2,3}, then everything.
        let m = merges_from_node_pairs(4, &[(1, 0), (3, 2), (2, 0)]).unwrap();
        assert_eq!(m, vec![(1, 0), (3, 2), (5, 4)]);
        assert_eq!(cut_dendrogram(4, &m, 0).unwrap(), vec![0, 1, 2, 3]);
        assert_eq!(cut_at_groups(4, &m, 2).unwrap(), vec![0, 0, 1, 1]);
        assert_eq!(cut_at_groups(4, &m, 1).unwrap(), vec![0, 0, 0, 0]);
        assert!(cut_at_groups(4, &m[..2], 1).is_err());
        assert!(check_merges(3, &[(0, 3)]).is_err());
        assert!(check_merges(3, &[(0, 1), (1, 2)]).is_err());
        assert!(merges_from_node_pairs(3, &[(0, 1), (1, 0)]).is_err());
    }

    #[test]
    fn edge_removals() {
        // Path 0-1-2-3 cut in the middle first, then at the ends.
        let m = merges_from_edge_removals(4, &[(1, 2), (0, 1), (2, 3)]).unwrap();
        assert_eq!(m, vec![(2, 3), (0, 1), (5, 4)]);
        assert_eq!(cut_at_groups(4, &m, 2).unwrap(), vec![0, 0, 1, 1]);
    }

    #[test]
    fn modularity_path_matches_direct_evaluation() {
        let g = und(
            7,
            &[
                (0, 1),
                (0, 2),
                (1, 2),
                (3, 4),
                (3, 5),
                (4, 5),
                (2, 3),
                (6, 6),
                (5, 6),
                (5, 6),
            ],
        );
        let m =
            merges_from_node_pairs(7, &[(0, 1), (0, 2), (3, 4), (5, 6), (3, 5), (0, 3)]).unwrap();
        let path = dendrogram_modularity(&g, &m, 1.0, false).unwrap();
        for (s, &q) in path.iter().enumerate() {
            let memb: Vec<i64> = cut_dendrogram(7, &m, s)
                .unwrap()
                .iter()
                .map(|&c| c as i64)
                .collect();
            assert!((q - modularity(&g, &memb, 1.0, false).unwrap()).abs() < 1e-14);
        }
        let (memb, q) = cut_at_max_modularity(&g, &m, 1.0, false).unwrap();
        // The loop and the doubled edge at node 6 make {5, 6} a group.
        assert_eq!(memb, vec![0, 0, 0, 1, 1, 2, 2]);
        assert!((q - path[4]).abs() < 1e-15);
        assert!(path.iter().all(|&x| x <= q));
    }
}
