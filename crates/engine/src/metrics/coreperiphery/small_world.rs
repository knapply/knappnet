//! The small-world quotient.
//!
//! Watts and Strogatz (1998, `watts1998collective`) call a network a small
//! world when its clustering C is much higher than that of a random
//! network of the same size and density while its characteristic path
//! length L is about as short. The small-world quotient (Uzzi and Spiro
//! 2005, `uzzi2005collaboration`; "small-world-ness" in Humphries and
//! Gurney 2008, `humphries2008network`) combines the two ratios:
//!
//! ```text
//! Q = (C / C_random) / (L / L_random).
//! ```
//!
//! Here C_random and L_random are the means over `reps` simulated G(n, m)
//! networks with the same number of nodes and edges
//! ([`crate::generate::gnm`]), rather than the approximations
//! C_random ~ k / n and L_random ~ ln n / ln k, which fail for sparse or
//! small networks.
//!
//! Conventions:
//! - Everything is computed on the undirected simple reading of the network
//!   (directions ignored, parallel edges once, loops ignored), and weights
//!   are not used; the random networks are undirected and simple, with as
//!   many edges as that reading.
//! - C is either the global transitivity (three times the triangles over
//!   the connected triples) or the mean local clustering coefficient of
//!   Watts and Strogatz with nodes of fewer than two neighbors counted as
//!   0 (NetworkX's `average_clustering`), as chosen; see
//!   [`crate::metrics::descriptive`].
//! - L is the mean distance over pairs of nodes joined by a path (igraph's
//!   `mean_distance(unconnected = TRUE)`), in the observed and in each
//!   random network.
//! - Replicate r uses stream (seed, r), and replicates run in parallel,
//!   so results do not depend on the thread count (up to the order of
//!   floating-point sums of the replicate means, which is fixed).
//! - Undefined parts (no connected triples, no pairs at distance > 0) give
//!   NaN.

use crate::error::{Error, Result};
use crate::generate::gnm;
use crate::graph::Graph;
use crate::metrics::descriptive::{
    average_clustering, mean_distance, transitivity, ClusteringOptions, PathOptions, Undefined,
};
use crate::metrics::topography::{simplify, ParallelWeight};
use crate::parallel;

/// Options for [`small_world`].
#[derive(Debug, Clone, Copy)]
pub struct SmallWorldOptions {
    /// Use the mean local clustering coefficient (true) or global
    /// transitivity (false).
    pub local: bool,
    /// Number of random networks.
    pub reps: usize,
    /// Seed for the random networks.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Result of [`small_world`].
#[derive(Debug, Clone, PartialEq)]
pub struct SmallWorld {
    /// Clustering of the network.
    pub clustering: f64,
    /// Mean distance of the network.
    pub mean_distance: f64,
    /// Mean clustering of the random networks.
    pub random_clustering: f64,
    /// Mean of the random networks' mean distances.
    pub random_mean_distance: f64,
    /// The small-world quotient.
    pub quotient: f64,
}

fn measures(g: &Graph, local: bool) -> Result<(f64, f64)> {
    let copts = ClusteringOptions {
        undefined: if local {
            Undefined::Zero
        } else {
            Undefined::Nan
        },
        threads: 1,
        ..Default::default()
    };
    let c = if local {
        average_clustering(g, copts)?
    } else {
        transitivity(g, copts)
    };
    let popts = PathOptions {
        threads: 1,
        ..Default::default()
    };
    Ok((c, mean_distance(g, popts, true)?))
}

/// Small-world quotient of a network. See the module documentation.
pub fn small_world(g: &Graph, opts: SmallWorldOptions) -> Result<SmallWorld> {
    if opts.reps == 0 {
        return Err(Error::InvalidArgument(
            "the number of random networks must be at least 1".into(),
        ));
    }
    let simple = simplify(g, true, ParallelWeight::Max)?;
    let n = simple.node_count();
    let m = simple.edge_count() as u64;
    let (c, l) = measures(&simple, opts.local)?;
    let nets = gnm(n, m, false, false, opts.reps, opts.seed, opts.threads)?;
    let threads = parallel::threads_for(
        opts.threads,
        opts.reps,
        parallel::search_work(opts.reps, n, 2 * m as usize).saturating_mul(n as u64),
    );
    let parts = parallel::for_each_claimed(opts.reps, threads, Vec::new, |acc, r| {
        let e = &nets[r];
        let s: Vec<i64> = e.from.iter().map(|&x| x as i64).collect();
        let t: Vec<i64> = e.to.iter().map(|&x| x as i64).collect();
        let res = Graph::from_edges(n, &s, &t, None, false).and_then(|h| measures(&h, opts.local));
        acc.push((r, res));
    });
    let mut all: Vec<(usize, Result<(f64, f64)>)> = parts.into_iter().flatten().collect();
    all.sort_by_key(|x| x.0);
    let (mut sc, mut sl) = (0.0, 0.0);
    for (_, res) in all {
        let (cr, lr) = res?;
        sc += cr;
        sl += lr;
    }
    let reps = opts.reps as f64;
    let (cr, lr) = (sc / reps, sl / reps);
    Ok(SmallWorld {
        clustering: c,
        mean_distance: l,
        random_clustering: cr,
        random_mean_distance: lr,
        quotient: (c / cr) / (l / lr),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ring_lattice_is_a_small_world_after_rewiring() {
        let e = crate::generate::small_world(200, 8, 0.05, 1, 3, 1)
            .unwrap()
            .remove(0);
        let s: Vec<i64> = e.from.iter().map(|&x| x as i64).collect();
        let t: Vec<i64> = e.to.iter().map(|&x| x as i64).collect();
        let g = Graph::from_edges(200, &s, &t, None, false).unwrap();
        for local in [false, true] {
            let opts = SmallWorldOptions {
                local,
                reps: 20,
                seed: 1,
                threads: 2,
            };
            let r = small_world(&g, opts).unwrap();
            assert!(r.quotient > 3.0, "{r:?}");
            // Random clustering is near the density, 8 / 199.
            assert!((r.random_clustering - 8.0 / 199.0).abs() < 0.01, "{r:?}");
            let again = small_world(&g, SmallWorldOptions { threads: 1, ..opts }).unwrap();
            assert_eq!(r, again);
        }
        let opts = SmallWorldOptions {
            local: false,
            reps: 0,
            seed: 1,
            threads: 1,
        };
        assert!(small_world(&g, opts).is_err());
    }
}
