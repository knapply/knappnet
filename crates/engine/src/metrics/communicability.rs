//! Communicability and communicability betweenness (Estrada and Hatano
//! 2008, `estrada2008communicability`; Estrada, Higham, and Hatano 2009,
//! `estrada2009communicability`).
//!
//! The communicability of nodes p and q is G_pq = (e^A)_pq, the number of
//! walks from p to q weighted by 1/k! for walks of length k. With
//! A = V diag(lambda) V' the eigendecomposition of the symmetric adjacency
//! matrix ([`crate::linalg::symmetric_eigen`]), G_pq = sum over j of
//! V_pj V_qj e^lambda_j.
//!
//! Communicability betweenness of r is
//! omega_r = 1/C sum over p != q, both != r, of (G_pq - G_pq(r)) / G_pq,
//! where G(r) is the communicability of the network with the edges of r
//! removed and C = (n - 1)^2 - (n - 1), the number of ordered pairs of
//! other nodes (for n > 2). Each r needs its own eigendecomposition, so
//! the cost is O(n^4); the nodes are processed in parallel and each value
//! depends only on its node, so results do not depend on the thread count.
//!
//! Conventions:
//! - The binary adjacency matrix of the undirected reading: direction and
//!   weights are ignored, parallel edges count once, and a self-loop puts 1
//!   on the diagonal, as NetworkX's `communicability` and
//!   `communicability_betweenness_centrality` on graphs (NetworkX refuses
//!   directed graphs and multigraphs).
//! - Disconnected networks: e^A is computed per component, so nodes in
//!   different components have communicability exactly 0 (instead of
//!   rounding noise), and pairs in different components contribute 0 to
//!   communicability betweenness. NetworkX divides by those zero entries and
//!   returns NaN or meaningless values there.

use crate::error::Result;
use crate::graph::Graph;
use crate::linalg;
use crate::matrix::Dense;
use crate::parallel;
use crate::union_find::UnionFind;

/// The binary symmetric adjacency (with loops on the diagonal) and the
/// components, each a list of nodes in increasing order.
fn binary(g: &Graph) -> (Vec<Vec<bool>>, Vec<Vec<usize>>) {
    let n = g.node_count();
    let mut a = vec![vec![false; n]; n];
    let mut uf = UnionFind::new(n);
    for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
        let (s, t) = (s as usize, t as usize);
        a[s][t] = true;
        a[t][s] = true;
        uf.union(s, t);
    }
    let mut comps: Vec<Vec<usize>> = Vec::new();
    let mut label = vec![usize::MAX; n];
    for v in 0..n {
        let r = uf.find(v);
        if label[r] == usize::MAX {
            label[r] = comps.len();
            comps.push(Vec::new());
        }
        comps[label[r]].push(v);
    }
    (a, comps)
}

/// e^A for the adjacency restricted to `nodes` without the node `skip`
/// (as a matrix over `nodes`, with zero row and column for `skip`).
fn expm(a: &[Vec<bool>], nodes: &[usize], skip: Option<usize>) -> Result<Dense> {
    let kept: Vec<usize> = (0..nodes.len()).filter(|&i| Some(i) != skip).collect();
    let k = kept.len();
    let mut m = Dense::zeros(k, k);
    for (i, &pi) in kept.iter().enumerate() {
        for (j, &pj) in kept.iter().enumerate() {
            if a[nodes[pi]][nodes[pj]] {
                m.set(i, j, 1.0);
            }
        }
    }
    let (values, vectors) = linalg::symmetric_eigen(&m)?;
    let e: Vec<f64> = values.iter().map(|l| l.exp()).collect();
    let mut out = Dense::zeros(nodes.len(), nodes.len());
    for (i, &pi) in kept.iter().enumerate() {
        for (j, &pj) in kept.iter().enumerate().skip(i) {
            let mut s = 0.0;
            for (l, el) in e.iter().enumerate() {
                s += vectors.get(i, l) * vectors.get(j, l) * el;
            }
            out.set(pi, pj, s);
            out.set(pj, pi, s);
        }
    }
    Ok(out)
}

/// The communicability matrix (n by n). See the module documentation.
pub fn communicability(g: &Graph) -> Result<Dense> {
    let n = g.node_count();
    let (a, comps) = binary(g);
    let mut out = Dense::zeros(n, n);
    for nodes in &comps {
        let e = expm(&a, nodes, None)?;
        for (i, &p) in nodes.iter().enumerate() {
            for (j, &q) in nodes.iter().enumerate() {
                out.set(p, q, e.get(i, j));
            }
        }
    }
    Ok(out)
}

/// Communicability betweenness of every node. See the module
/// documentation.
pub fn communicability_betweenness(g: &Graph, threads: usize) -> Result<Vec<f64>> {
    let n = g.node_count();
    let (a, comps) = binary(g);
    let mut comp_of = vec![0usize; n];
    let mut pos = vec![0usize; n];
    let mut full: Vec<Dense> = Vec::with_capacity(comps.len());
    for (c, nodes) in comps.iter().enumerate() {
        for (i, &v) in nodes.iter().enumerate() {
            comp_of[v] = c;
            pos[v] = i;
        }
        full.push(expm(&a, nodes, None)?);
    }
    let work: u64 = comps.iter().map(|c| (c.len() as u64).pow(3)).sum();
    let threads = parallel::threads_for(threads, n, work.saturating_mul(10));
    let values: Vec<Result<f64>> = parallel::map_items(
        n,
        threads,
        || (),
        |_, r| {
            let c = comp_of[r];
            let nodes = &comps[c];
            let g0 = &full[c];
            let gr = expm(&a, nodes, Some(pos[r]))?;
            let k = nodes.len();
            let mut total = 0.0;
            for p in 0..k {
                if p == pos[r] {
                    continue;
                }
                for q in 0..k {
                    if q == p || q == pos[r] {
                        continue;
                    }
                    let base = g0.get(p, q);
                    total += (base - gr.get(p, q)) / base;
                }
            }
            Ok(total)
        },
    );
    let scale = if n > 2 {
        let m = n as f64 - 1.0;
        1.0 / (m * m - m)
    } else {
        1.0
    };
    values.into_iter().map(|v| v.map(|x| x * scale)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_and_star() {
        // Path 0-1: e^A = [[cosh 1, sinh 1], [sinh 1, cosh 1]].
        let p = Graph::from_edges(2, &[0], &[1], None, false).unwrap();
        let c = communicability(&p).unwrap();
        assert!((c.get(0, 0) - 1f64.cosh()).abs() < 1e-12);
        assert!((c.get(0, 1) - 1f64.sinh()).abs() < 1e-12);
        // In a star every walk between leaves passes the centre, so the
        // centre's betweenness is 1 per pair of leaves; a leaf only lies on
        // longer walks.
        let s = Graph::from_edges(4, &[0, 0, 0], &[1, 2, 3], None, false).unwrap();
        let b = communicability_betweenness(&s, 2).unwrap();
        assert!((b[0] - 1.0).abs() < 1e-12, "{b:?}");
        assert!(b[1] > 0.0 && b[1] < 1.0 && (b[1] - b[2]).abs() < 1e-12);
    }
}
