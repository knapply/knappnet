//! Measures built on all shortest-path distances: center, periphery, and
//! barycenter; the harmonic diameter; the Wiener, Schultz, Gutman, and
//! hyper-Wiener indices; global, local, and pairwise efficiency; and
//! closeness vitality.
//!
//! Sources:
//! - **Center, periphery, barycenter**: the nodes of smallest and largest
//!   eccentricity (Hage and Harary 1995, `hage1995eccentricity`), and the
//!   nodes of smallest total distance to the others, the median of
//!   Smart and Slater (1999, `smart1999center`).
//! - **Harmonic diameter**: n(n - 1) divided by the sum of 1/d over ordered
//!   pairs, the harmonic mean distance of Marchiori and Latora (2000,
//!   `marchiori2000harmony`).
//! - **Wiener index** (Wiener 1947, `wiener1947structural`): the sum of
//!   distances over unordered pairs (ordered pairs in directed networks).
//!   **Schultz index** (Schultz 1989, `schultz1989topological`; in the
//!   degree-distance form of Gutman 1994, `gutman1994selected`): the sum
//!   over unordered pairs of (deg u + deg v) d(u, v). **Gutman index**
//!   (Gutman 1994): the sum over unordered pairs of deg u deg v d(u, v).
//!   **Hyper-Wiener index** (Randic 1993, `randic1993novel`; for graphs
//!   with cycles, Klein, Lukovits, and Gutman 1995, `klein1995hyper`): half
//!   the sum over unordered pairs of d + d^2.
//! - **Efficiency** (Latora and Marchiori 2001, `latora2001efficient`): the
//!   efficiency of a pair is 1/d (0 without a path); global efficiency
//!   averages it over ordered pairs of distinct nodes; local efficiency
//!   averages, over the nodes, the global efficiency of the subgraph
//!   induced by each node's neighbors.
//! - **Closeness vitality** (Koschutzki et al. 2005, `koschutzki2005centrality`,
//!   section 3.6): the Wiener index of the network minus that of the
//!   network without the node.
//!
//! Conventions:
//! - Distances are those of [`crate::metrics::descriptive::distances`]:
//!   `mode` for directed networks, weights as distances or strengths
//!   (positive), parallel edges through the lightest, self-loops ignored.
//!   Directed networks use the distances from each node (NetworkX's
//!   `shortest_path_length(G, source)`).
//! - Disconnected networks: the Wiener family is infinite, as in NetworkX
//!   (a directed network must be strongly connected). Closeness vitality is
//!   then W(G) - W(G - v) with infinities: infinite when only G is
//!   disconnected, minus infinity when only G - v is, and NaN when both are
//!   (NetworkX computes the same differences); it is NaN for a single node.
//!   Center, periphery, and barycenter are defined for (strongly) connected
//!   networks, and are empty otherwise (NetworkX raises); with
//!   [`Reachable::Ignore`] center and periphery instead use igraph's
//!   eccentricity over the nodes each node reaches.
//! - Ties: eccentricities are compared exactly; total distances for the
//!   barycenter are compared with [`Tie`] (sums in different orders can
//!   differ in the last bits; integer lengths are exact).
//! - The Schultz and Gutman indices use vertex degrees (parallel edges
//!   count, a self-loop counts twice) even when distances are weighted,
//!   which is what NetworkX's `schultz_index` and `gutman_index` compute.
//!   They and the hyper-Wiener index are defined for undirected networks.
//! - The hyper-Wiener index: [`HyperWiener::Standard`] is the published
//!   index, (1/2) sum over unordered pairs of d + d^2 (so the path on
//!   three nodes has 5); NetworkX's `hyper_wiener_index` documents and
//!   computes the sum over unordered pairs of d + d^2, twice the published
//!   value ([`HyperWiener::NetworkX`]).
//! - Efficiency: global efficiency is 0 for fewer than two nodes. Local
//!   efficiency uses the undirected reading (distinct neighbors other than
//!   the node), and unweighted distances inside each neighborhood, as
//!   NetworkX's `local_efficiency`; NetworkX also counts a node as its own
//!   neighbor when it has a self-loop, which `loop_is_neighbor` reproduces.
//!   It is NaN for a network without nodes.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::{degree, DegreeMode};
use crate::metrics::descriptive::{distance_rows, eccentricity, path_summary, PathOptions};
use crate::parallel;
use crate::paths::Tie;

/// Per-source totals of one pass over the distance rows.
#[derive(Debug, Clone, Copy, Default)]
struct RowStats {
    reached: usize,
    sum: f64,
    sum_sq: f64,
    sum_inv: f64,
    longest: f64,
    sum_deg: f64,
}

fn row_stats(g: &Graph, opts: PathOptions, deg: Option<&[f64]>) -> Result<Vec<RowStats>> {
    let n = g.node_count();
    let mut out = vec![RowStats::default(); n];
    distance_rows(g, None, opts, |k, row| {
        let mut s = RowStats::default();
        for (v, &d) in row.iter().enumerate() {
            if v == k || d.is_infinite() {
                continue;
            }
            s.reached += 1;
            s.sum += d;
            s.sum_sq += d * d;
            if d != 0.0 {
                s.sum_inv += 1.0 / d;
            }
            s.longest = s.longest.max(d);
            if let Some(deg) = deg {
                s.sum_deg += d * deg[v];
            }
        }
        out[k] = s;
    })?;
    Ok(out)
}

fn connected(stats: &[RowStats]) -> bool {
    let n = stats.len();
    stats.iter().all(|s| s.reached + 1 == n)
}

/// How center and periphery treat nodes that do not reach every other.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Reachable {
    /// Defined only for (strongly) connected networks; empty otherwise
    /// (NetworkX).
    #[default]
    Require,
    /// Eccentricity over the nodes each node reaches (igraph).
    Ignore,
}

impl Reachable {
    /// Decodes the C interface codes: 0 = require, 1 = ignore.
    pub fn from_code(code: i32) -> Result<Reachable> {
        match code {
            0 => Ok(Reachable::Require),
            1 => Ok(Reachable::Ignore),
            _ => Err(Error::InvalidArgument(format!(
                "unknown code {code} for unreachable nodes"
            ))),
        }
    }
}

fn eccentricities(g: &Graph, opts: PathOptions, reach: Reachable) -> Result<Option<Vec<f64>>> {
    match reach {
        Reachable::Ignore => Ok(Some(eccentricity(g, opts)?)),
        Reachable::Require => {
            let stats = row_stats(g, opts, None)?;
            if !connected(&stats) {
                return Ok(None);
            }
            Ok(Some(stats.iter().map(|s| s.longest).collect()))
        }
    }
}

/// Nodes of smallest eccentricity. See the module documentation.
pub fn center(g: &Graph, opts: PathOptions, reach: Reachable) -> Result<Vec<u32>> {
    Ok(match eccentricities(g, opts, reach)? {
        None => Vec::new(),
        Some(e) => {
            let best = e.iter().copied().fold(f64::INFINITY, f64::min);
            (0..e.len() as u32)
                .filter(|&v| e[v as usize] == best)
                .collect()
        }
    })
}

/// Nodes of largest eccentricity. See the module documentation.
pub fn periphery(g: &Graph, opts: PathOptions, reach: Reachable) -> Result<Vec<u32>> {
    Ok(match eccentricities(g, opts, reach)? {
        None => Vec::new(),
        Some(e) => {
            let best = e.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            (0..e.len() as u32)
                .filter(|&v| e[v as usize] == best)
                .collect()
        }
    })
}

/// The barycenter (nodes of smallest total distance to the others) and
/// every node's total distance (infinite when it does not reach every
/// other node). The barycenter is empty unless the network is (strongly)
/// connected.
pub fn barycenter(g: &Graph, opts: PathOptions, tie: Tie) -> Result<(Vec<u32>, Vec<f64>)> {
    let stats = row_stats(g, opts, None)?;
    let n = stats.len();
    let totals: Vec<f64> = stats
        .iter()
        .map(|s| {
            if s.reached + 1 == n {
                s.sum
            } else {
                f64::INFINITY
            }
        })
        .collect();
    if !connected(&stats) {
        return Ok((Vec::new(), totals));
    }
    let best = totals.iter().copied().fold(f64::INFINITY, f64::min);
    let nodes = (0..n as u32)
        .filter(|&v| tie.equal(totals[v as usize], best))
        .collect();
    Ok((nodes, totals))
}

/// Harmonic diameter: n(n - 1) / the sum of 1/d over ordered pairs
/// (infinite when no pair is joined and n > 1, NaN for n <= 1).
pub fn harmonic_diameter(g: &Graph, opts: PathOptions) -> Result<f64> {
    let n = g.node_count() as f64;
    let stats = row_stats(g, opts, None)?;
    let inv: f64 = stats.iter().map(|s| s.sum_inv).sum();
    Ok(if inv != 0.0 {
        n * (n - 1.0) / inv
    } else if n > 1.0 {
        f64::INFINITY
    } else {
        f64::NAN
    })
}

/// Wiener index: the sum of distances over unordered pairs (ordered pairs
/// in directed networks); infinite when the network is not (strongly)
/// connected, NaN without nodes.
pub fn wiener_index(g: &Graph, opts: PathOptions) -> Result<f64> {
    if g.node_count() == 0 {
        return Ok(f64::NAN);
    }
    let s = path_summary(g, None, opts)?;
    let n = g.node_count() as u64;
    if s.pairs != n * (n - 1) {
        return Ok(f64::INFINITY);
    }
    Ok(if g.is_directed() && opts.mode != DegreeMode::All {
        s.total
    } else {
        s.total / 2.0
    })
}

fn require_undirected(g: &Graph, what: &str) -> Result<()> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(format!(
            "the {what} is defined for undirected networks"
        )));
    }
    Ok(())
}

/// Schultz index (degree distance). See the module documentation.
pub fn schultz_index(g: &Graph, opts: PathOptions) -> Result<f64> {
    require_undirected(g, "Schultz index")?;
    let deg = degree(g, DegreeMode::All, true);
    let stats = row_stats(g, opts, Some(&deg))?;
    if g.node_count() == 0 {
        return Ok(f64::NAN);
    }
    if !connected(&stats) {
        return Ok(f64::INFINITY);
    }
    let total: f64 = stats
        .iter()
        .enumerate()
        .map(|(u, s)| deg[u] * s.sum + s.sum_deg)
        .sum();
    Ok(total / 2.0)
}

/// Gutman index. See the module documentation.
pub fn gutman_index(g: &Graph, opts: PathOptions) -> Result<f64> {
    require_undirected(g, "Gutman index")?;
    let deg = degree(g, DegreeMode::All, true);
    let stats = row_stats(g, opts, Some(&deg))?;
    if g.node_count() == 0 {
        return Ok(f64::NAN);
    }
    if !connected(&stats) {
        return Ok(f64::INFINITY);
    }
    let total: f64 = stats
        .iter()
        .enumerate()
        .map(|(u, s)| deg[u] * s.sum_deg)
        .sum();
    Ok(total / 2.0)
}

/// Which hyper-Wiener index. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HyperWiener {
    /// (1/2) sum over unordered pairs of d + d^2 (Randic; Klein et al.).
    #[default]
    Standard,
    /// The sum over unordered pairs of d + d^2 (NetworkX).
    NetworkX,
}

impl HyperWiener {
    /// Decodes the C interface codes: 0 = standard, 1 = NetworkX.
    pub fn from_code(code: i32) -> Result<HyperWiener> {
        match code {
            0 => Ok(HyperWiener::Standard),
            1 => Ok(HyperWiener::NetworkX),
            _ => Err(Error::InvalidArgument(format!(
                "unknown hyper-Wiener convention code {code}"
            ))),
        }
    }
}

/// Hyper-Wiener index. See the module documentation.
pub fn hyper_wiener_index(g: &Graph, opts: PathOptions, which: HyperWiener) -> Result<f64> {
    require_undirected(g, "hyper-Wiener index")?;
    if g.node_count() == 0 {
        return Ok(f64::NAN);
    }
    let stats = row_stats(g, opts, None)?;
    if !connected(&stats) {
        return Ok(f64::INFINITY);
    }
    // Ordered sums count each unordered pair twice.
    let ordered: f64 = stats.iter().map(|s| s.sum + s.sum_sq).sum();
    Ok(match which {
        HyperWiener::Standard => ordered / 4.0,
        HyperWiener::NetworkX => ordered / 2.0,
    })
}

/// Global efficiency: the mean of 1/d over ordered pairs of distinct nodes
/// (0 for fewer than two nodes).
pub fn global_efficiency(g: &Graph, opts: PathOptions) -> Result<f64> {
    let n = g.node_count();
    if n < 2 {
        return Ok(0.0);
    }
    let stats = row_stats(g, opts, None)?;
    let total: f64 = stats.iter().map(|s| s.sum_inv).sum();
    Ok(total / (n as f64 * (n as f64 - 1.0)))
}

/// Local efficiency (unweighted, undirected reading). See the module
/// documentation.
pub fn local_efficiency(g: &Graph, loop_is_neighbor: bool, threads: usize) -> Result<f64> {
    let n = g.node_count();
    if n == 0 {
        return Ok(f64::NAN);
    }
    let s = g.simple_undirected(false);
    let looped: Vec<bool> = {
        let mut l = vec![false; n];
        for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
            if a == b {
                l[a as usize] = true;
            }
        }
        l
    };
    let work: u64 = (0..n)
        .map(|v| (s.degree(v) as u64).pow(2))
        .sum::<u64>()
        .saturating_mul(4);
    let threads = parallel::threads_for(threads, n, work);
    let values = parallel::map_items(
        n,
        threads,
        || (vec![u32::MAX; n], vec![u32::MAX; n]),
        |(member, dist), v| {
            let mut nodes: Vec<u32> = s.neighbors(v).to_vec();
            if loop_is_neighbor && looped[v] {
                let at = nodes.partition_point(|&x| (x as usize) < v);
                nodes.insert(at, v as u32);
            }
            let k = nodes.len();
            if k < 2 {
                return 0.0;
            }
            for (i, &x) in nodes.iter().enumerate() {
                member[x as usize] = i as u32;
            }
            let mut total = 0.0;
            let mut queue: Vec<u32> = Vec::with_capacity(k);
            for &root in &nodes {
                queue.clear();
                queue.push(root);
                dist[root as usize] = 0;
                let mut head = 0;
                while head < queue.len() {
                    let x = queue[head] as usize;
                    head += 1;
                    for &y in s.neighbors(x) {
                        let y = y as usize;
                        if member[y] == u32::MAX || dist[y] != u32::MAX {
                            continue;
                        }
                        dist[y] = dist[x] + 1;
                        total += 1.0 / dist[y] as f64;
                        queue.push(y as u32);
                    }
                }
                for &x in &queue {
                    dist[x as usize] = u32::MAX;
                }
            }
            for &x in &nodes {
                member[x as usize] = u32::MAX;
            }
            total / (k as f64 * (k as f64 - 1.0))
        },
    );
    Ok(values.iter().sum::<f64>() / n as f64)
}

/// Closeness vitality of every node. See the module documentation.
pub fn closeness_vitality(g: &Graph, opts: PathOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    if n <= 1 {
        return Ok(vec![f64::NAN; n]);
    }
    let whole = wiener_index(g, opts)?;
    let halve = !(g.is_directed() && opts.mode != DegreeMode::All);
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let w = g.edge_weights();
    let mut out = Vec::with_capacity(n);
    for v in 0..n as u32 {
        let keep: Vec<usize> = (0..g.edge_count())
            .filter(|&e| src[e] != v && tgt[e] != v)
            .collect();
        let s: Vec<i64> = keep.iter().map(|&e| src[e] as i64).collect();
        let t: Vec<i64> = keep.iter().map(|&e| tgt[e] as i64).collect();
        let kw: Option<Vec<f64>> = w.map(|w| keep.iter().map(|&e| w[e]).collect());
        let h = Graph::from_edges(n, &s, &t, kw.as_deref(), g.is_directed())?;
        let sources: Vec<u32> = (0..n as u32).filter(|&u| u != v).collect();
        let summary = path_summary(&h, Some(&sources), opts)?;
        let others = (n - 1) as u64;
        let after = if summary.pairs == others * (others.saturating_sub(1)) {
            if halve {
                summary.total / 2.0
            } else {
                summary.total
            }
        } else {
            f64::INFINITY
        };
        out.push(whole - after);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(n: usize) -> Graph {
        let s: Vec<i64> = (0..n as i64 - 1).collect();
        let t: Vec<i64> = (1..n as i64).collect();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn indices_of_paths() {
        let p = path(3);
        let o = PathOptions::default();
        assert_eq!(wiener_index(&p, o).unwrap(), 4.0);
        assert_eq!(
            hyper_wiener_index(&p, o, HyperWiener::Standard).unwrap(),
            5.0
        );
        assert_eq!(
            hyper_wiener_index(&p, o, HyperWiener::NetworkX).unwrap(),
            10.0
        );
        // Degrees 1, 2, 1: Schultz = (1+2)1 + (2+1)1 + (1+1)2 = 10;
        // Gutman = 2 + 2 + 2 = 6.
        assert_eq!(schultz_index(&p, o).unwrap(), 10.0);
        assert_eq!(gutman_index(&p, o).unwrap(), 6.0);
        assert_eq!(center(&p, o, Reachable::Require).unwrap(), vec![1]);
        assert_eq!(periphery(&p, o, Reachable::Require).unwrap(), vec![0, 2]);
        assert_eq!(barycenter(&p, o, Tie::Exact).unwrap().0, vec![1]);
        // Pairs: 1, 1, 1/2 each way: 5/2 * 2 = 5; 6 / 5.
        assert!((harmonic_diameter(&p, o).unwrap() - 6.0 / 5.0).abs() < 1e-12);
        assert!((global_efficiency(&p, o).unwrap() - 5.0 / 6.0).abs() < 1e-12);
        // Removing the middle disconnects the rest; removing an end leaves
        // one edge (Wiener 1).
        let v = closeness_vitality(&p, o).unwrap();
        assert_eq!(v, vec![3.0, f64::NEG_INFINITY, 3.0]);
    }

    #[test]
    fn local_efficiency_of_a_triangle_with_a_pendant() {
        // Triangle 0-1-2 and pendant 2-3: node 2's neighbors {0, 1, 3}
        // have one edge (0-1): efficiency (1 + 1) / 6 = 1/3. Nodes 0 and 1
        // have adjacent neighbor pairs (efficiency 1); node 3 has one
        // neighbor (0).
        let g = Graph::from_edges(4, &[0, 1, 0, 2], &[1, 2, 2, 3], None, false).unwrap();
        let e = local_efficiency(&g, false, 2).unwrap();
        assert!((e - (1.0 + 1.0 + 1.0 / 3.0) / 4.0).abs() < 1e-12);
    }
}
