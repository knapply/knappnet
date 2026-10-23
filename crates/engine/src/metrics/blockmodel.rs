//! Blockmodels: block densities, image matrices, and fit for a partition.
//!
//! Given a partition of the nodes into blocks (positions), the network's
//! adjacency matrix is divided into submatrices, one per ordered pair of
//! blocks. A blockmodel reduces each submatrix to one number, its density,
//! and the image matrix codes each block as a 1-block or a 0-block
//! (White, Boorman, and Breiger 1976, `white1976social`; Wasserman and
//! Faust 1994, `wasserman1994social`, ch. 10).
//!
//! - **Density** of block (r, s): the sum of a_ij over i in r and j in s
//!   (i != j), divided by the number of such pairs, n_r n_s, or
//!   n_r (n_r - 1) when r = s. With `loops = true` the diagonal counts:
//!   a_ii is included and the diagonal blocks have n_r^2 pairs.
//! - **Image** at a cutoff alpha: 1 where the block density is at least
//!   alpha, 0 elsewhere (the alpha density criterion of Arabie, Boorman,
//!   and Levitt 1978, `arabie1978constructing`; Wasserman and Faust 1994,
//!   ch. 10). The usual cutoff is the density of the whole network
//!   (with the same treatment of the diagonal), which [`BlockModel`]
//!   reports.
//! - **Fit**: the squared Pearson correlation between the observed a_ij
//!   and the density of the block each pair falls in, over the same pairs
//!   (Wasserman and Faust 1994, ch. 10): the share of the variance in
//!   the ties that the block densities account for.
//!
//! Conventions:
//! - Without weights a_ij is 1 when there is an edge from i to j and 0
//!   otherwise, so parallel edges count once; with weights a_ij is the
//!   total weight of the edges from i to j, used as given (density is then
//!   the mean tie value, as in sna's `blockmodel(block.content =
//!   "density")` on a valued matrix). An undirected edge is a tie in both
//!   directions, so the density matrix of an undirected network is
//!   symmetric.
//! - Self-loops are ignored unless `loops = true`.
//! - Blocks are numbered from 0 in order of their first node; only equality
//!   of the supplied labels matters.
//! - A block with no pairs (a one-node block on the diagonal without
//!   loops) has density NaN and image NaN. The fit is NaN when either the
//!   ties or the fitted densities have no variance, or there are no pairs.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::matrix::Dense;

/// Options for [`block_model`].
#[derive(Debug, Clone, Copy, Default)]
pub struct BlockModelOptions {
    /// Use total edge weights instead of 0/1 ties.
    pub weighted: bool,
    /// Count self-loops and the diagonal of the adjacency matrix.
    pub loops: bool,
}

/// Block sums, pair counts, and densities for a partition.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockModel {
    /// Block of every node, numbered from 0 in order of first node.
    pub membership: Vec<u32>,
    /// Number of nodes in each block.
    pub sizes: Vec<usize>,
    /// Sum of tie values in each block (k by k).
    pub sums: Dense,
    /// Number of pairs of nodes in each block (k by k).
    pub pairs: Dense,
    /// Block densities, sums / pairs (k by k).
    pub density: Dense,
    /// Density of the whole network, with the same treatment of the
    /// diagonal: total sum / total pairs.
    pub overall: f64,
    /// Squared correlation between ties and block densities.
    pub r_squared: f64,
}

impl BlockModel {
    /// Number of blocks.
    pub fn blocks(&self) -> usize {
        self.sizes.len()
    }

    /// The binary image: 1 where a block's density is at least `cutoff`,
    /// 0 where it is less, NaN where the density is undefined.
    pub fn image(&self, cutoff: f64) -> Dense {
        let k = self.blocks();
        let mut out = Dense::zeros(k, k);
        for r in 0..k {
            for s in 0..k {
                let d = self.density.get(r, s);
                out.set(
                    r,
                    s,
                    if d.is_nan() {
                        f64::NAN
                    } else if d >= cutoff {
                        1.0
                    } else {
                        0.0
                    },
                );
            }
        }
        out
    }
}

/// Block densities and fit of the partition `membership` (one label per
/// node). See the module documentation.
pub fn block_model(g: &Graph, membership: &[i64], opts: BlockModelOptions) -> Result<BlockModel> {
    let n = g.node_count();
    if membership.len() != n {
        return Err(Error::LengthMismatch {
            what: "membership",
            expected: n,
            found: membership.len(),
        });
    }
    let w = if opts.weighted {
        Some(g.edge_weights().ok_or(Error::MissingWeights)?)
    } else {
        None
    };
    let mut seen: HashMap<i64, u32> = HashMap::new();
    let block: Vec<u32> = membership
        .iter()
        .map(|&l| {
            let len = seen.len() as u32;
            *seen.entry(l).or_insert(len)
        })
        .collect();
    let k = seen.len();
    let mut sizes = vec![0usize; k];
    for &b in &block {
        sizes[b as usize] += 1;
    }

    // Tie values a_ij for the ordered pairs with a tie.
    let mut a: HashMap<(u32, u32), f64> = HashMap::new();
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        if s == t && !opts.loops {
            continue;
        }
        let x = w.map_or(1.0, |w| w[e]);
        let mut add = |i: u32, j: u32| {
            let v = a.entry((i, j)).or_insert(0.0);
            *v = if opts.weighted { *v + x } else { 1.0 };
        };
        add(s, t);
        if !g.is_directed() && s != t {
            add(t, s);
        }
    }

    let mut sums = Dense::zeros(k, k);
    let mut pairs = Dense::zeros(k, k);
    for r in 0..k {
        for s in 0..k {
            let (nr, ns) = (sizes[r] as f64, sizes[s] as f64);
            let p = if r == s && !opts.loops {
                nr * (nr - 1.0)
            } else {
                nr * ns
            };
            pairs.set(r, s, p);
        }
    }
    // Sum in a fixed order so results do not depend on hashing.
    let mut entries: Vec<((u32, u32), f64)> = a.into_iter().collect();
    entries.sort_by_key(|e| e.0);
    for &((i, j), x) in &entries {
        let (r, s) = (block[i as usize] as usize, block[j as usize] as usize);
        sums.set(r, s, sums.get(r, s) + x);
    }
    let mut density = Dense::zeros(k, k);
    for r in 0..k {
        for s in 0..k {
            let p = pairs.get(r, s);
            density.set(
                r,
                s,
                if p > 0.0 {
                    sums.get(r, s) / p
                } else {
                    f64::NAN
                },
            );
        }
    }
    let total_pairs: f64 = pairs.as_slice().iter().sum();
    let total_sum: f64 = sums.as_slice().iter().sum();
    let overall = if total_pairs > 0.0 {
        total_sum / total_pairs
    } else {
        f64::NAN
    };

    // Correlation between a_ij and the density of its block, over every
    // counted pair. Pairs without a tie contribute a_ij = 0, so the moments
    // are accumulated per block and then corrected with the tie entries.
    let r_squared = if total_pairs > 0.0 {
        let mut sum_f = 0.0; // sum of fitted values
        let mut sum_ff = 0.0; // sum of squared fitted values
        for r in 0..k {
            for s in 0..k {
                let p = pairs.get(r, s);
                if p > 0.0 {
                    let d = density.get(r, s);
                    sum_f += p * d;
                    sum_ff += p * d * d;
                }
            }
        }
        let mut sum_a = 0.0;
        let mut sum_aa = 0.0;
        let mut sum_af = 0.0;
        for &((i, j), x) in &entries {
            let (r, s) = (block[i as usize] as usize, block[j as usize] as usize);
            let d = density.get(r, s);
            sum_a += x;
            sum_aa += x * x;
            sum_af += x * d;
        }
        let m = total_pairs;
        let cov = sum_af - sum_a * sum_f / m;
        let var_a = sum_aa - sum_a * sum_a / m;
        let var_f = sum_ff - sum_f * sum_f / m;
        // Rounding can leave a tiny non-zero variance where the exact one
        // is zero (all densities equal); treat that as no variance.
        let scale = (sum_aa.abs() + sum_ff.abs()).max(f64::MIN_POSITIVE);
        if var_a <= 1e-12 * scale || var_f <= 1e-12 * scale {
            f64::NAN
        } else {
            (cov * cov) / (var_a * var_f)
        }
    } else {
        f64::NAN
    };

    Ok(BlockModel {
        membership: block,
        sizes,
        sums,
        pairs,
        density,
        overall,
        r_squared,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], w: Option<&[f64]>, directed: bool) -> Graph {
        let s: Vec<i64> = e.iter().map(|x| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        Graph::from_edges(n, &s, &t, w, directed).unwrap()
    }

    #[test]
    fn densities_of_two_cliques_joined_by_one_edge() {
        // Triangles {0,1,2} and {3,4,5} and the edge 2-3.
        let e = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)];
        let g = graph(6, &e, None, false);
        let bm = block_model(&g, &[7, 7, 7, 9, 9, 9], Default::default()).unwrap();
        assert_eq!(bm.sizes, vec![3, 3]);
        assert_eq!(bm.density.get(0, 0), 1.0);
        assert_eq!(bm.density.get(1, 1), 1.0);
        assert_eq!(bm.density.get(0, 1), 1.0 / 9.0);
        assert_eq!(bm.density.get(1, 0), 1.0 / 9.0);
        assert_eq!(bm.overall, 14.0 / 30.0);
        let img = bm.image(bm.overall);
        assert_eq!(img.as_slice(), &[1.0, 0.0, 0.0, 1.0]);
        // R^2 by direct computation over the 30 ordered pairs.
        let blocks = [0, 0, 0, 1, 1, 1];
        let adj = |i: usize, j: usize| {
            e.iter().any(|&(a, b)| {
                (a as usize, b as usize) == (i, j) || (b as usize, a as usize) == (i, j)
            }) as u8 as f64
        };
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        for i in 0..6 {
            for j in 0..6 {
                if i != j {
                    xs.push(adj(i, j));
                    ys.push(bm.density.get(blocks[i], blocks[j]));
                }
            }
        }
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let (mx, my) = (mean(&xs), mean(&ys));
        let cov: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
        let vx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
        let vy: f64 = ys.iter().map(|y| (y - my) * (y - my)).sum();
        let r2 = cov * cov / (vx * vy);
        assert!(
            (bm.r_squared - r2).abs() < 1e-12,
            "{} vs {r2}",
            bm.r_squared
        );
    }

    #[test]
    fn directed_weighted_loops_and_singletons() {
        // 0 -> 1 (weight 2, twice), 1 -> 0 (3), 0 -> 0 loop (5), node 2 alone.
        let g = graph(
            3,
            &[(0, 1), (0, 1), (1, 0), (0, 0)],
            Some(&[2.0, 2.0, 3.0, 5.0]),
            true,
        );
        let bin = block_model(&g, &[0, 0, 1], Default::default()).unwrap();
        assert_eq!(bin.density.get(0, 0), 1.0);
        assert!(bin.density.get(1, 1).is_nan());
        assert_eq!(bin.density.get(0, 1), 0.0);
        assert!(bin.image(0.5).get(1, 1).is_nan());
        let w = block_model(
            &g,
            &[0, 0, 1],
            BlockModelOptions {
                weighted: true,
                loops: true,
            },
        )
        .unwrap();
        // Block (0,0) with loops: 4 pairs, sum 4 + 3 + 5.
        assert_eq!(w.density.get(0, 0), 12.0 / 4.0);
        assert_eq!(w.density.get(1, 1), 0.0);
        assert_eq!(w.pairs.get(1, 1), 1.0);
        assert!(block_model(&g, &[0, 1], Default::default()).is_err());
        // A single block: fitted values have no variance.
        let one = block_model(&g, &[4, 4, 4], Default::default()).unwrap();
        assert!(one.r_squared.is_nan());
    }
}
