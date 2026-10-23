//! Auxiliary statistics for goodness of fit (Lospinoso and Snijders 2019,
//! `lospinoso2019goodness`).
//!
//! The fit of an estimated model is judged by comparing statistics of the
//! observed networks with their distribution over the networks simulated
//! in phase 3 of the estimation. The engine computes, for one network,
//! the four families RSiena's `sienaGOF` offers for one-mode networks:
//!
//! - the in-degree distribution: the number of actors with in-degree
//!   0, 1, ..., n - 1;
//! - the out-degree distribution, likewise;
//! - the triad census in MAN order (16 counts; Holland and Leinhardt 1976,
//!   `holland1976local`);
//! - the geodesic distribution: the number of ordered pairs of distinct
//!   actors at directed distance 1, 2, ..., n - 1, and the number with no
//!   directed path.
//!
//! The statistics are summed over periods by the caller, who also chooses
//! the levels to keep, whether to cumulate them, and computes the
//! Mahalanobis distance and its Monte Carlo p-value (which need a
//! generalized inverse of the simulated covariance matrix).

use super::network::Network;
use crate::graph::Graph;
use crate::metrics::triad::triad_census;

/// Number of auxiliary statistics for a network on `n` actors.
pub fn aux_len(n: usize) -> usize {
    3 * n + 16
}

/// The auxiliary statistics of `net`, in the order listed in the module
/// documentation.
pub fn aux_statistics(net: &Network) -> Vec<f64> {
    let n = net.n();
    let mut out = vec![0.0; aux_len(n)];
    if n == 0 {
        return out;
    }
    for i in 0..n {
        out[net.indeg(i) as usize] += 1.0;
        out[n + net.outdeg(i) as usize] += 1.0;
    }
    let (src, tgt): (Vec<i64>, Vec<i64>) = (0..n)
        .flat_map(|i| net.out_neighbors(i).map(move |j| (i as i64, j as i64)))
        .unzip();
    let g = Graph::from_edges(n, &src, &tgt, None, true).expect("valid indices");
    let census = triad_census(&g, 1);
    out[2 * n..2 * n + 16].copy_from_slice(&census);
    // Breadth-first search from every actor.
    let base = 2 * n + 16;
    let mut dist = vec![usize::MAX; n];
    let mut queue = Vec::with_capacity(n);
    for s in 0..n {
        dist.iter_mut().for_each(|d| *d = usize::MAX);
        dist[s] = 0;
        queue.clear();
        queue.push(s);
        let mut head = 0;
        while head < queue.len() {
            let v = queue[head];
            head += 1;
            for w in net.out_neighbors(v) {
                if dist[w] == usize::MAX {
                    dist[w] = dist[v] + 1;
                    queue.push(w);
                }
            }
        }
        for (t, &d) in dist.iter().enumerate() {
            if t == s {
                continue;
            }
            if d == usize::MAX {
                out[base + n - 1] += 1.0;
            } else {
                out[base + d - 1] += 1.0;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn statistics_of_a_small_network() {
        // 0 -> 1 -> 2, 0 -> 2, 3 isolated.
        let net = Network::from_adjacency(4, &[0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        let s = aux_statistics(&net);
        let n = 4;
        assert_eq!(&s[..n], &[2.0, 1.0, 1.0, 0.0]); // in-degrees 0,1,2,0
        assert_eq!(&s[n..2 * n], &[2.0, 1.0, 1.0, 0.0]); // out-degrees 2,1,0,0
        let census = &s[2 * n..2 * n + 16];
        assert_eq!(census.iter().sum::<f64>(), 4.0);
        assert_eq!(census[8], 1.0); // 030T: the transitive triple 0, 1, 2
        let geo = &s[2 * n + 16..];
        // Distance 1: three pairs; no pair at distance 2 (0 reaches 2
        // directly); unreachable: 12 - 3 = 9.
        assert_eq!(geo[0], 3.0);
        assert_eq!(geo[1], 0.0);
        assert_eq!(geo[n - 1], 9.0);
    }
}
