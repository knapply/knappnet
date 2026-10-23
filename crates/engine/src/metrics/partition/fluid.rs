//! Fluid communities.
//!
//! Pares, Garcia-Gasulla, Vilalta, Moreno, Ayguade, Labarta, Cortes, and
//! Suzumura (2018, `pares2018fluid`) place k communities on k distinct
//! random nodes, each with density 1, where the density of a community is
//! 1 / (its number of nodes). Nodes are then visited in a new random order
//! in every sweep; a node joins the community with the largest total
//! density over the node and its neighbors, keeping its own community when
//! that is among the largest and otherwise choosing uniformly among them.
//! The method stops after a sweep in which no node changes community.
//! Because a community's density at any node is at most 1 and a node that
//! is a community's last member sees density 1 from it, no community ever
//! disappears: the result has exactly k groups.
//!
//! Conventions:
//! - Undirected, connected networks with at least two nodes, as igraph's
//!   `cluster_fluid_communities()` requires (it also requires a simple
//!   network and rejects parallel edges and loops; here parallel edges
//!   count once and loops are ignored, which is the same neighborhood).
//!   Directed networks are rejected rather than symmetrized.
//! - Weights are not used.
//! - 1 <= k <= n (igraph requires k < n).
//! - Total densities are compared exactly as fractions (count / size),
//!   so ties are exact; NetworkX's `asyn_fluidc()` compares them with a
//!   tolerance of 1e-4, which makes no difference unless communities
//!   exceed about 100 nodes.
//! - Random numbers come from [`crate::rng::Stream`] with the seed (the
//!   starting nodes from stream 0, sweep s from stream s + 1), so a seed
//!   reproduces the result on any platform. The method runs on one
//!   thread, and stops after at most 100,000 sweeps.
//! - Groups are labelled 0, 1, ... in order of their first node.

use super::{relabel, require_undirected};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::rng::Stream;

const MAX_SWEEPS: u64 = 100_000;

/// Fluid communities of a connected undirected network with `k`
/// communities. See the module documentation.
pub fn fluid_communities(g: &Graph, k: usize, seed: u64) -> Result<Vec<usize>> {
    require_undirected(g, "Fluid communities")?;
    let n = g.node_count();
    if n < 2 {
        return Err(Error::InvalidArgument(
            "fluid communities need a network with at least two nodes".into(),
        ));
    }
    if k == 0 || k > n {
        return Err(Error::InvalidArgument(format!(
            "the number of communities must be between 1 and {n}, not {k}"
        )));
    }
    let adj = g.simple_undirected(false);
    // Connectedness by breadth-first search from node 0.
    let mut seen = vec![false; n];
    let mut stack = vec![0usize];
    seen[0] = true;
    while let Some(v) = stack.pop() {
        for &u in adj.neighbors(v) {
            if !seen[u as usize] {
                seen[u as usize] = true;
                stack.push(u as usize);
            }
        }
    }
    if seen.iter().any(|&s| !s) {
        return Err(Error::InvalidArgument(
            "fluid communities need a connected network".into(),
        ));
    }
    const NONE: usize = usize::MAX;
    let mut comm = vec![NONE; n];
    let mut size = vec![0u64; k];
    let mut perm: Vec<usize> = (0..n).collect();
    Stream::new(seed, 0).permutation(&mut perm);
    for (c, &v) in perm[..k].iter().enumerate() {
        comm[v] = c;
        size[c] = 1;
    }
    let mut count = vec![0u64; k];
    let mut touched: Vec<usize> = Vec::new();
    let mut best: Vec<usize> = Vec::new();
    for sweep in 0..MAX_SWEEPS {
        let mut rng = Stream::new(seed, sweep + 1);
        rng.permutation(&mut perm);
        let mut changed = false;
        for &v in &perm {
            for u in std::iter::once(v).chain(adj.neighbors(v).iter().map(|&u| u as usize)) {
                let c = comm[u];
                if c != NONE {
                    if count[c] == 0 {
                        touched.push(c);
                    }
                    count[c] += 1;
                }
            }
            // Largest count / size, compared exactly.
            best.clear();
            for &c in &touched {
                match best.first() {
                    None => best.push(c),
                    Some(&b) => {
                        let (lhs, rhs) = (count[c] * size[b], count[b] * size[c]);
                        if lhs > rhs {
                            best.clear();
                            best.push(c);
                        } else if lhs == rhs {
                            best.push(c);
                        }
                    }
                }
            }
            for &c in &touched {
                count[c] = 0;
            }
            touched.clear();
            if best.is_empty() || best.contains(&comm[v]) {
                continue;
            }
            best.sort_unstable();
            let c = best[rng.below(best.len() as u64) as usize];
            if comm[v] != NONE {
                size[comm[v]] -= 1;
            }
            comm[v] = c;
            size[c] += 1;
            changed = true;
        }
        if !changed {
            break;
        }
    }
    Ok(relabel(&comm))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn two_cliques() {
        let mut e = Vec::new();
        for c in 0..2i64 {
            for i in 0..5 {
                for j in (i + 1)..5 {
                    e.push((5 * c + i, 5 * c + j));
                }
            }
        }
        e.push((4, 5));
        let g = und(10, &e);
        let mut found = 0;
        for seed in 0..20 {
            let m = fluid_communities(&g, 2, seed).unwrap();
            assert_eq!(*m.iter().max().unwrap(), 1);
            assert_eq!(m, fluid_communities(&g, 2, seed).unwrap());
            if m == vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1] {
                found += 1;
            }
        }
        assert!(found >= 10, "{found}");
        assert_eq!(
            fluid_communities(&g, 10, 1).unwrap(),
            (0..10).collect::<Vec<_>>()
        );
        assert_eq!(fluid_communities(&g, 1, 1).unwrap(), vec![0; 10]);
    }

    #[test]
    fn requirements() {
        assert!(fluid_communities(&und(4, &[(0, 1), (2, 3)]), 2, 1).is_err());
        assert!(fluid_communities(&und(1, &[]), 1, 1).is_err());
        assert!(fluid_communities(&und(2, &[(0, 1)]), 3, 1).is_err());
        let d = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        assert!(fluid_communities(&d, 1, 1).is_err());
    }
}
