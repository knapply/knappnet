//! Square clustering and generalized degree.
//!
//! [`square_clustering`] is the share of the possible squares (4-cycles)
//! through a node that exist (Lind, González, and Herrmann 2005,
//! `lind2005cycles`), in the form NetworkX 3.7's `square_clustering`
//! computes (`hagberg2008exploring`, after Zhang et al. 2008,
//! `zhang2008clustering`): for node v with neighbour set N (k = |N|),
//!
//! C4(v) = q / (Σ_{u ∈ N} |N(u)| (k − 1) − k (k − 1) − t − q),
//!
//! where t = Σ_{u ∈ N} |N(u) ∩ N| counts the triangles at v twice and q,
//! the number of squares through v, is half of Σ_x p_x (p_x − 1) over
//! every node x ≠ v that is a neighbour or a neighbour's neighbour, p_x =
//! |N ∩ N(x)|; C4 is 0 when k < 2 or the denominator is not positive.
//!
//! [`generalized_degree`] (Zlatić, Garlaschelli, and Caldarelli 2012,
//! `zlatic2012networks`, as NetworkX's `generalized_degree`) counts, for every node, its edges by the number
//! of triangles they lie in: entry k is the number of neighbours w with
//! |N(v) ∩ N(w)| = k.
//!
//! Conventions: self-loops are ignored and parallel edges merged (neighbour
//! sets); weights are ignored. Square clustering of a directed graph uses
//! successor sets, as NetworkX does with `G._adj`; generalized degree is
//! defined for undirected graphs (direction is ignored for directed ones,
//! which NetworkX rejects). Values are computed per node in parallel and
//! do not depend on the thread count.

use crate::error::Result;
use crate::graph::Graph;
use crate::parallel;

/// Distinct out-neighbours (all neighbours when undirected or when
/// `undirected` is set) without the node itself.
fn rows(g: &Graph, undirected: bool) -> Vec<Vec<u32>> {
    let n = g.node_count();
    (0..n)
        .map(|v| {
            let mut r: Vec<u32> = g.out().neighbors(v).to_vec();
            if undirected && g.is_directed() {
                r.extend_from_slice(g.inc().neighbors(v));
                r.sort_unstable();
            }
            r.dedup();
            r.retain(|&u| u as usize != v);
            r
        })
        .collect()
}

/// Square clustering of every node. See the module documentation.
pub fn square_clustering(g: &Graph, threads: usize) -> Result<Vec<f64>> {
    let n = g.node_count();
    let adj = rows(g, false);
    let arcs: usize = adj.iter().map(Vec::len).sum();
    let avg = arcs / n.max(1) + 1;
    let work = (n as u64).saturating_mul((avg * avg * avg) as u64);
    let threads = parallel::threads_for(threads, n, work);
    Ok(parallel::map_items(
        n,
        threads,
        || (vec![false; n], vec![false; n], Vec::new()),
        |(in_n, seen, hop), v| {
            let nv = &adj[v];
            let k = nv.len() as i64;
            if k - 1 <= 0 {
                return 0.0;
            }
            for &u in nv {
                in_n[u as usize] = true;
            }
            let common = |x: usize, in_n: &[bool]| -> i64 {
                adj[x].iter().filter(|&&y| in_n[y as usize]).count() as i64
            };
            let (mut uw, mut tri, mut sq) = (0i64, 0i64, 0i64);
            hop.clear();
            for &u in nv {
                let u = u as usize;
                uw += adj[u].len() as i64 * (k - 1);
                let p = common(u, in_n);
                tri += p;
                sq += p * (p - 1);
                for &x in &adj[u] {
                    let x = x as usize;
                    if x != v && !in_n[x] && !seen[x] {
                        seen[x] = true;
                        hop.push(x);
                    }
                }
            }
            for &x in hop.iter() {
                let p = common(x, in_n);
                sq += p * (p - 1);
                seen[x] = false;
            }
            for &u in nv {
                in_n[u as usize] = false;
            }
            let sq = sq / 2;
            let potential = uw - k * (k - 1) - tri - sq;
            if potential > 0 {
                sq as f64 / potential as f64
            } else {
                0.0
            }
        },
    ))
}

/// Generalized degree of every node: (k, number of neighbours sharing k
/// triangles with it) for every k that occurs, increasing. See the module
/// documentation.
pub fn generalized_degree(g: &Graph, threads: usize) -> Vec<Vec<(u32, u32)>> {
    let n = g.node_count();
    let adj = rows(g, true);
    let arcs: usize = adj.iter().map(Vec::len).sum();
    let work = (arcs as u64).saturating_mul(arcs as u64 / n.max(1) as u64 + 1);
    let threads = parallel::threads_for(threads, n, work);
    parallel::map_items(
        n,
        threads,
        || vec![false; n],
        |in_n, v| {
            for &u in &adj[v] {
                in_n[u as usize] = true;
            }
            let mut ks: Vec<u32> = adj[v]
                .iter()
                .map(|&w| {
                    adj[w as usize]
                        .iter()
                        .filter(|&&y| in_n[y as usize])
                        .count() as u32
                })
                .collect();
            for &u in &adj[v] {
                in_n[u as usize] = false;
            }
            ks.sort_unstable();
            let mut out: Vec<(u32, u32)> = Vec::new();
            for k in ks {
                match out.last_mut() {
                    Some(last) if last.0 == k => last.1 += 1,
                    _ => out.push((k, 1)),
                }
            }
            out
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn squares_of_a_four_cycle_and_a_complete_graph() {
        // In a 4-cycle every node is on its one possible square.
        let c4 = Graph::from_edges(4, &[0, 1, 2, 3], &[1, 2, 3, 0], None, false).unwrap();
        assert_eq!(square_clustering(&c4, 1).unwrap(), vec![1.0; 4]);
        // K4: each node: k = 3, uw = 3 * 3 * 2 = 18, k(k-1) = 6, t = 6,
        // q = 3 (three 4-cycles through it), so 3 / (18 - 6 - 6 - 3) = 1.
        let k4 =
            Graph::from_edges(4, &[0, 0, 0, 1, 1, 2], &[1, 2, 3, 2, 3, 3], None, false).unwrap();
        assert_eq!(square_clustering(&k4, 1).unwrap(), vec![1.0; 4]);
        let gd = generalized_degree(&k4, 1);
        assert_eq!(gd[0], vec![(2, 3)]);
        assert_eq!(generalized_degree(&c4, 1)[0], vec![(0, 2)]);
    }
}
