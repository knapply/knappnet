//! VoteRank, Laplacian centrality, and dispersion.
//!
//! Sources:
//! - **VoteRank** (Zhang, Chen, Dong, and Zhao 2016, `zhang2016identifying`)
//!   picks a set of spreaders one at a time. Every node starts with voting
//!   ability 1; in each round a node's score is the sum of the abilities of
//!   the nodes that vote for it (its neighbors; in directed networks the
//!   heads of its out-arcs, as in NetworkX's `voterank`), the node with the
//!   highest score is chosen, its own ability drops to 0, and each node it
//!   reaches by an edge loses 1/k of ability (k the average degree), down to
//!   0. The selection stops early when the best score is 0.
//! - **Laplacian centrality** (Qi, Fuller, Wu, Wu, and Zhang 2012,
//!   `qi2012laplacian`): the drop in Laplacian energy (the sum of squared
//!   entries of the weighted Laplacian L = D - A) when a node is removed,
//!   where removing node i deletes its row and column and lowers each
//!   neighbor's diagonal entry by the weight of its edge to i. In closed
//!   form the drop is L_ii^2 + 2 sum_j w_ij^2 + sum_j (L_jj^2 - (L_jj - |w_ij|)^2),
//!   O(m) for all nodes. Normalized, it is divided by the energy.
//! - **Dispersion** (Backstrom and Kleinberg 2014, `backstrom2014romantic`)
//!   of the tie between u and v: the number of pairs s, t of common
//!   neighbors of u and v that are not adjacent and have no common neighbor
//!   among u's neighbors other than u and v. Normalized,
//!   (dispersion + b)^alpha / (embeddedness + c), the embeddedness being the
//!   number of common neighbors; NetworkX's defaults are alpha = 1, b = 0,
//!   c = 0.
//!
//! Conventions:
//! - VoteRank: every edge votes, so parallel edges vote more than once and a
//!   self-loop lets a node vote for itself (twice in an undirected network,
//!   whose average degree counts a loop twice), exactly as NetworkX's loop
//!   over its edges; votes are summed in input order of the edges, the order
//!   NetworkX's `G.edges()` yields when the network was built from it, and
//!   among equal scores the node with the smallest index wins (NetworkX's
//!   `max` over its node order).
//! - Laplacian centrality uses the undirected reading: direction is ignored,
//!   parallel edges add their weights, and self-loops cancel out of the
//!   Laplacian, as in NetworkX's `laplacian_centrality` for undirected
//!   graphs (for directed graphs NetworkX uses the random-walk Laplacian of
//!   Chung, which is not implemented here). Weights are used as given. The
//!   normalized value is NaN when the energy is 0.
//! - Dispersion reads neighbor sets as NetworkX's adjacency does: distinct
//!   neighbors (out-neighbors in directed networks), including the node
//!   itself when it has a self-loop. In directed networks the pair test is
//!   not symmetric; pairs s, t are taken with s before t in node order.

use crate::error::{Error, Result};
use crate::graph::Graph;

/// VoteRank: up to `count` nodes (all when `None`), in the order chosen.
pub fn voterank(g: &Graph, count: Option<usize>) -> Vec<u32> {
    let n = g.node_count();
    let mut chosen = Vec::new();
    if n == 0 {
        return chosen;
    }
    let k = count.map_or(n, |c| c.min(n));
    let m = g.edge_count() as f64;
    let directed = g.is_directed();
    let avg = if directed {
        m / n as f64
    } else {
        2.0 * m / n as f64
    };
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut ability = vec![1.0f64; n];
    let mut score = vec![0.0f64; n];
    let mut is_chosen = vec![false; n];
    for _ in 0..k {
        score.iter_mut().for_each(|x| *x = 0.0);
        for (&s, &t) in src.iter().zip(tgt) {
            let (s, t) = (s as usize, t as usize);
            score[s] += ability[t];
            if !directed {
                score[t] += ability[s];
            }
        }
        for v in 0..n {
            if is_chosen[v] {
                score[v] = 0.0;
            }
        }
        let mut best = 0;
        for v in 1..n {
            if score[v] > score[best] {
                best = v;
            }
        }
        if score[best] == 0.0 {
            break;
        }
        chosen.push(best as u32);
        is_chosen[best] = true;
        ability[best] = 0.0;
        for &w in g.out().neighbors(best) {
            let w = w as usize;
            ability[w] = (ability[w] - 1.0 / avg).max(0.0);
        }
    }
    chosen
}

/// Laplacian centrality of every node. See the module documentation.
pub fn laplacian_centrality(g: &Graph, weighted: bool, normalized: bool) -> Result<Vec<f64>> {
    let n = g.node_count();
    // Parallel edges add up whether or not the network is weighted.
    let w: Vec<f64> = if weighted {
        g.edge_weights().ok_or(Error::MissingWeights)?.to_vec()
    } else {
        vec![1.0; g.edge_count()]
    };
    let src: Vec<i64> = g.edge_sources().iter().map(|&v| v as i64).collect();
    let tgt: Vec<i64> = g.edge_targets().iter().map(|&v| v as i64).collect();
    let s = Graph::from_edges(n, &src, &tgt, Some(&w), false)?.simple_undirected(true);
    // L_jj: the weighted degree without self-loops.
    let w_of = |v: usize, j: usize| s.weights(v).map_or(1.0, |w| w[j]);
    let diag: Vec<f64> = (0..n)
        .map(|v| (0..s.degree(v)).map(|j| w_of(v, j)).sum())
        .collect();
    let mut energy: f64 = diag.iter().map(|d| d * d).sum();
    for v in 0..n {
        for j in 0..s.degree(v) {
            energy += w_of(v, j) * w_of(v, j);
        }
    }
    Ok((0..n)
        .map(|i| {
            let mut drop = diag[i] * diag[i];
            for (j, &u) in s.neighbors(i).iter().enumerate() {
                let w = w_of(i, j);
                let lj = diag[u as usize];
                drop += 2.0 * w * w + lj * lj - (lj - w.abs()) * (lj - w.abs());
            }
            if normalized {
                if energy == 0.0 {
                    f64::NAN
                } else {
                    drop / energy
                }
            } else {
                drop
            }
        })
        .collect())
}

/// Options for [`dispersion`].
#[derive(Debug, Clone, Copy)]
pub struct DispersionOptions {
    /// Use the normalized form.
    pub normalized: bool,
    /// Exponent of the normalized form.
    pub alpha: f64,
    /// Added to the dispersion in the normalized form.
    pub b: f64,
    /// Added to the embeddedness in the normalized form.
    pub c: f64,
}

impl Default for DispersionOptions {
    fn default() -> Self {
        DispersionOptions {
            normalized: true,
            alpha: 1.0,
            b: 0.0,
            c: 0.0,
        }
    }
}

/// Dispersion of each pair (u, v) in `pairs`. See the module
/// documentation.
pub fn dispersion(g: &Graph, pairs: &[(u32, u32)], opts: DispersionOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    if let Some(&(u, v)) = pairs
        .iter()
        .find(|&&(u, v)| u as usize >= n || v as usize >= n)
    {
        return Err(Error::InvalidArgument(format!(
            "pair ({u}, {v}) refers to a node outside the network of {n} nodes"
        )));
    }
    // Distinct out-neighbors, the node itself included when looped.
    let nbrs: Vec<Vec<u32>> = (0..n)
        .map(|v| {
            let mut row = g.out().neighbors(v).to_vec();
            row.dedup();
            row
        })
        .collect();
    let has = |v: usize, x: u32| nbrs[v].binary_search(&x).is_ok();
    Ok(pairs
        .iter()
        .map(|&(u, v)| {
            let (u, v) = (u as usize, v as usize);
            let st: Vec<u32> = nbrs[v].iter().copied().filter(|&x| has(u, x)).collect();
            let mut total = 0u64;
            for (i, &s) in st.iter().enumerate() {
                // u's neighbors adjacent from s, other than u and v.
                let nbrs_s: Vec<u32> = nbrs[s as usize]
                    .iter()
                    .copied()
                    .filter(|&x| x as usize != u && x as usize != v && has(u, x))
                    .collect();
                for &t in &st[i + 1..] {
                    if nbrs_s.binary_search(&t).is_ok() {
                        continue;
                    }
                    if nbrs_s.iter().all(|&x| !has(t as usize, x)) {
                        total += 1;
                    }
                }
            }
            let emb = st.len() as f64;
            if opts.normalized {
                let x = (total as f64 + opts.b).powf(opts.alpha);
                if emb + opts.c != 0.0 {
                    x / (emb + opts.c)
                } else {
                    x
                }
            } else {
                total as f64
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn voterank_on_a_star_with_a_tail() {
        // Star centre 0 with leaves 1..4, and 4 - 5 - 6.
        let g =
            Graph::from_edges(7, &[0, 0, 0, 0, 4, 5], &[1, 2, 3, 4, 5, 6], None, false).unwrap();
        let v = voterank(&g, None);
        assert_eq!(v[0], 0);
        assert!(v.len() <= 7);
        assert_eq!(voterank(&g, Some(1)), vec![0]);
    }

    #[test]
    fn laplacian_centrality_of_a_path() {
        // Path 0-1-2 (unit weights): L energy = 1 + 4 + 1 + 4 * 1 = 10.
        // Removing 1 leaves two isolated nodes (energy 0): drop 10.
        // Removing 0 leaves the edge 1-2 (energy 4): drop 6.
        let g = Graph::from_edges(3, &[0, 1], &[1, 2], None, false).unwrap();
        let c = laplacian_centrality(&g, false, false).unwrap();
        assert_eq!(c, vec![6.0, 10.0, 6.0]);
        let c = laplacian_centrality(&g, false, true).unwrap();
        assert!((c[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dispersion_of_a_tie_in_two_groups() {
        // u = 0, v = 1 share neighbors 2, 3 (adjacent) and 4 (apart).
        let g = Graph::from_edges(
            5,
            &[0, 0, 0, 0, 1, 1, 1, 2],
            &[1, 2, 3, 4, 2, 3, 4, 3],
            None,
            false,
        )
        .unwrap();
        let d = dispersion(
            &g,
            &[(0, 1)],
            DispersionOptions {
                normalized: false,
                ..Default::default()
            },
        )
        .unwrap();
        // Pairs (2, 4) and (3, 4) are apart with no common neighbor among
        // 0's neighbors other than 0 and 1; (2, 3) are adjacent.
        assert_eq!(d, vec![2.0]);
        let d = dispersion(&g, &[(0, 1)], DispersionOptions::default()).unwrap();
        assert!((d[0] - 2.0 / 3.0).abs() < 1e-12);
    }
}
