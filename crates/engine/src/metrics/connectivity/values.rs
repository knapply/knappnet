//! Local and global node and edge connectivity.
//!
//! Definitions (Menger 1927, `menger1927allgemeinen`; Ford and Fulkerson
//! 1956, `ford1956maximal`):
//! - The local edge connectivity lambda(s, t) is the largest number of
//!   edge-disjoint paths from s to t, equal to the fewest edges whose
//!   removal leaves no path from s to t.
//! - The local node connectivity kappa(s, t) is the largest number of
//!   internally node-disjoint paths from s to t; an edge from s to t
//!   counts as one such path. For non-adjacent s and t it equals the
//!   fewest nodes (other than s and t) whose removal leaves no path.
//! - The (global) edge connectivity lambda(G) and node connectivity
//!   kappa(G) are the least local values over ordered pairs (unordered if
//!   undirected); kappa(G) is taken over non-adjacent pairs, and a
//!   complete graph on n nodes has kappa = n - 1 (Harary 1969,
//!   `harary1969graph`, chapter 5). A network that is not connected
//!   (weakly, if directed) has both equal to 0, and so does a single node.
//!
//! Algorithms: each local value is one unit-capacity maximum flow (Even
//! and Tarjan 1975, `even1975network`), on the edge network or the split
//! network of [`super`]. Global node connectivity follows Even's (1975,
//! `even1975algorithm`) reduction as NetworkX implements it: with v a node
//! of least degree (for a directed network, least of in- and out-degree)
//! and K its degree, it is the least of K and kappa over the pairs (v, w)
//! with w not adjacent to v (both orders if directed) and the
//! non-adjacent pairs of v's neighbours (Esfahanian 2013,
//! `esfahanian2013connectivity`). Global edge connectivity of an
//! undirected network is the least of its least degree and lambda(v, w)
//! over the members w of a dominating set D, v being D's first member:
//! when lambda is below the least degree, both sides of a minimum cut hold
//! a node all of whose neighbours are on its side, so D meets both sides
//! (Esfahanian and Hakimi 1984, `esfahanian1984computing`; D is built
//! greedily in node order). For a directed network it is the least of
//! lambda(v_i, v_{i+1}) around the cycle v_0, ..., v_{n-1}, v_0 (a minimum
//! cut separates some consecutive pair; Esfahanian 2013). Every flow stops
//! once it reaches the best value found so far.
//!
//! Conventions (see also [`super`]): simple graph, directions kept. The
//! values equal NetworkX's `node_connectivity`, `edge_connectivity`,
//! `all_pairs_node_connectivity`, and `average_node_connectivity`, except
//! that NetworkX raises an error for the edge connectivity of a directed
//! network with one node (this returns 0), and for an undirected complete
//! graph returns its least degree counting self-loops twice (this ignores
//! loops). A network with no nodes is an error.

use super::{check_pair, Simple};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::matrix::Dense;
use crate::parallel;

/// kappa(s, t) on a prepared split network, stopping at `limit`.
pub(crate) fn kappa(
    net: &mut crate::flow::FlowNetwork,
    n: usize,
    s: usize,
    t: usize,
    limit: Option<f64>,
) -> u64 {
    net.reset();
    net.max_flow(n + s, t, limit)
        .expect("finite unit capacities and distinct terminals") as u64
}

/// lambda(s, t) on a prepared edge network, stopping at `limit`.
pub(crate) fn lambda(
    net: &mut crate::flow::FlowNetwork,
    s: usize,
    t: usize,
    limit: Option<f64>,
) -> u64 {
    net.reset();
    net.max_flow(s, t, limit)
        .expect("finite unit capacities and distinct terminals") as u64
}

/// Local node connectivity kappa(s, t). See the module documentation.
pub fn local_node_connectivity(g: &Graph, s: usize, t: usize) -> Result<u64> {
    check_pair(s, t, g.node_count())?;
    let simple = Simple::new(g);
    let (mut net, _) = simple.node_network(1.0);
    Ok(kappa(&mut net, simple.n, s, t, None))
}

/// Local edge connectivity lambda(s, t). See the module documentation.
pub fn local_edge_connectivity(g: &Graph, s: usize, t: usize) -> Result<u64> {
    check_pair(s, t, g.node_count())?;
    let mut net = Simple::new(g).edge_network();
    Ok(lambda(&mut net, s, t, None))
}

fn nonempty(g: &Graph) -> Result<()> {
    if g.node_count() == 0 {
        Err(Error::InvalidArgument(
            "connectivity is not defined for a network with no nodes".into(),
        ))
    } else {
        Ok(())
    }
}

/// The node v of least degree (in-, out- or the smaller of both if
/// directed), its degree, and the ordered pairs whose local node
/// connectivity decides kappa(G) (Even's reduction).
pub(crate) fn even_pairs(simple: &Simple) -> (usize, usize, Vec<(usize, usize)>) {
    let n = simple.n;
    let degree = |v: usize| {
        if simple.directed {
            simple.out[v].len().min(simple.inc[v].len())
        } else {
            simple.out[v].len()
        }
    };
    let v = (0..n).min_by_key(|&v| degree(v)).unwrap_or(0);
    let k = degree(v);
    let nbrs = simple.neighbors_any(v);
    let mut is_nbr = vec![false; n];
    for &w in &nbrs {
        is_nbr[w as usize] = true;
    }
    let mut pairs = Vec::new();
    for (w, &adjacent) in is_nbr.iter().enumerate() {
        if w != v && !adjacent {
            pairs.push((v, w));
            if simple.directed {
                pairs.push((w, v));
            }
        }
    }
    for (i, &x) in nbrs.iter().enumerate() {
        for (j, &y) in nbrs.iter().enumerate() {
            let (x, y) = (x as usize, y as usize);
            if simple.directed {
                if i != j && !simple.has_arc(x, y) {
                    pairs.push((x, y));
                }
            } else if i < j && !simple.has_arc(x, y) {
                pairs.push((x, y));
            }
        }
    }
    (v, k, pairs)
}

/// Global node connectivity kappa(G). See the module documentation.
pub fn node_connectivity(g: &Graph, threads: usize) -> Result<u64> {
    nonempty(g)?;
    let simple = Simple::new(g);
    if !simple.is_connected() {
        return Ok(0);
    }
    let (_, k, pairs) = even_pairs(&simple);
    let (base, _) = simple.node_network(1.0);
    let n = simple.n;
    let parts = parallel::for_each_item(
        pairs.len(),
        threads,
        || (base.clone(), k as u64),
        |(net, best), i| {
            let (s, t) = pairs[i];
            let f = kappa(net, n, s, t, Some(*best as f64));
            *best = (*best).min(f);
        },
    );
    Ok(parts.into_iter().map(|p| p.1).min().unwrap_or(k as u64))
}

/// The ordered pairs whose local edge connectivity, with the least degree,
/// decides lambda(G): around the node cycle if directed, from the first
/// member of a greedy dominating set to the others if undirected.
pub(crate) fn edge_pairs(simple: &Simple) -> Vec<(usize, usize)> {
    let n = simple.n;
    if simple.directed {
        return (0..n).map(|i| (i, (i + 1) % n)).collect();
    }
    let mut dominated = vec![false; n];
    let mut set = Vec::new();
    for v in 0..n {
        if !dominated[v] {
            set.push(v);
            dominated[v] = true;
            for &w in &simple.out[v] {
                dominated[w as usize] = true;
            }
        }
    }
    set.iter().skip(1).map(|&w| (set[0], w)).collect()
}

/// Global edge connectivity lambda(G). See the module documentation.
pub fn edge_connectivity(g: &Graph, threads: usize) -> Result<u64> {
    nonempty(g)?;
    let simple = Simple::new(g);
    let n = simple.n;
    if !simple.is_connected() || n == 1 {
        return Ok(0);
    }
    let bound = (0..n)
        .map(|v| {
            if simple.directed {
                simple.out[v].len().min(simple.inc[v].len())
            } else {
                simple.out[v].len()
            }
        })
        .min()
        .unwrap_or(0) as u64;
    let pairs = edge_pairs(&simple);
    let base = simple.edge_network();
    let parts = parallel::for_each_item(
        pairs.len(),
        threads,
        || (base.clone(), bound),
        |(net, best), i| {
            let (s, t) = pairs[i];
            let f = lambda(net, s, t, Some(*best as f64));
            *best = (*best).min(f);
        },
    );
    Ok(parts.into_iter().map(|p| p.1).min().unwrap_or(bound))
}

/// kappa(u, v) for every ordered pair, as an n-by-n matrix with zeros on
/// the diagonal (symmetric if undirected). NetworkX's
/// `all_pairs_node_connectivity` holds the same values (it omits the
/// diagonal).
pub fn all_pairs_node_connectivity(g: &Graph, threads: usize) -> Dense {
    let simple = Simple::new(g);
    let n = simple.n;
    let (base, _) = simple.node_network(1.0);
    let directed = simple.directed;
    // Row u: kappa(u, v) for v > u (every v != u if directed).
    let parts = parallel::for_each_item(
        n,
        threads,
        || (base.clone(), Vec::new()),
        |(net, out): &mut (crate::flow::FlowNetwork, Vec<(usize, usize, f64)>), u| {
            for v in 0..n {
                if v == u || (!directed && v < u) {
                    continue;
                }
                let f = kappa(net, n, u, v, None);
                out.push((u, v, f as f64));
            }
        },
    );
    let mut m = Dense::zeros(n, n);
    for (_, list) in parts {
        for (u, v, f) in list {
            m.set(u, v, f);
            if !directed {
                m.set(v, u, f);
            }
        }
    }
    m
}

/// Mean of the off-diagonal entries of a connectivity matrix from
/// [`all_pairs_node_connectivity`] over ordered pairs (unordered if
/// undirected), 0 with fewer than two nodes: NetworkX's
/// `average_node_connectivity`.
pub fn average_node_connectivity(m: &Dense, directed: bool) -> f64 {
    let n = m.rows();
    if n < 2 {
        return 0.0;
    }
    let mut total: u64 = 0;
    let mut pairs: u64 = 0;
    for u in 0..n {
        for v in 0..n {
            if u == v || (!directed && v < u) {
                continue;
            }
            total += m.get(u, v) as u64;
            pairs += 1;
        }
    }
    total as f64 / pairs as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn small_examples() {
        // Two triangles sharing node 2: kappa = 1, lambda = 2.
        let g = graph(5, &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)], false);
        assert_eq!(node_connectivity(&g, 1).unwrap(), 1);
        assert_eq!(edge_connectivity(&g, 2).unwrap(), 2);
        assert_eq!(local_node_connectivity(&g, 0, 4).unwrap(), 1);
        assert_eq!(local_node_connectivity(&g, 0, 1).unwrap(), 2);
        assert_eq!(local_edge_connectivity(&g, 0, 4).unwrap(), 2);
        // Complete graph K5 (with a loop and a repeated edge).
        let mut e = Vec::new();
        for a in 0..5 {
            for b in (a + 1)..5 {
                e.push((a, b));
            }
        }
        e.push((0, 0));
        e.push((0, 1));
        let k5 = graph(5, &e, false);
        assert_eq!(node_connectivity(&k5, 0).unwrap(), 4);
        assert_eq!(edge_connectivity(&k5, 0).unwrap(), 4);
        assert_eq!(
            average_node_connectivity(&all_pairs_node_connectivity(&k5, 2), false),
            4.0
        );
        // Directed cycle: 1 both ways.
        let c = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0)], true);
        assert_eq!(node_connectivity(&c, 1).unwrap(), 1);
        assert_eq!(edge_connectivity(&c, 1).unwrap(), 1);
        let m = all_pairs_node_connectivity(&c, 1);
        assert_eq!(m.get(0, 2), 1.0);
        // Disconnected and trivial networks.
        let d = graph(3, &[(0, 1)], false);
        assert_eq!(node_connectivity(&d, 1).unwrap(), 0);
        assert_eq!(edge_connectivity(&graph(1, &[], true), 1).unwrap(), 0);
        assert!(node_connectivity(&graph(0, &[], false), 1).is_err());
        assert!(local_node_connectivity(&d, 1, 1).is_err());
    }
}
