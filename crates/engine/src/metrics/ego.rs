//! Ego networks: neighborhoods and egocentric measures.
//!
//! An ego network is a focal node (ego), the nodes it is tied to (its
//! alters), and the ties among them (Wasserman and Faust 1994,
//! `wasserman1994social`, section 2.3; Perry, Pescosolido, and Borgatti
//! 2018, `perry2018egocentric`). With `order` k the alters are every node
//! within k steps of ego, which generalizes the first-order neighborhood.
//!
//! Measures computed for every ego at once, in parallel:
//! - size: the number of alters;
//! - ties: the number of ties among the alters (ego's own ties excluded);
//! - pairs: the number of pairs of alters, k(k - 1)/2, or k(k - 1) ordered
//!   pairs in a directed network, for k alters;
//! - density: ties / pairs (ego network density as in Marsden 2002,
//!   `marsden2002egocentric`, and UCINET's ego network measures);
//! - components: the number of weak components among the alters once ego
//!   is removed (isolated alters are components of size 1);
//! - ego betweenness: the betweenness of ego within its own ego network
//!   (Everett and Borgatti 2005, `everett2005ego`; Marsden 2002), computed
//!   with Brandes' accumulation (Brandes 2001, `brandes2001faster`) on the
//!   ego network alone.
//!
//! Burt's structural-hole measures, which are also egocentric, live in
//! [`crate::metrics::brokerage`].
//!
//! Conventions:
//! - Alters in directed networks are chosen by `mode`: `Out` takes the
//!   nodes ego reaches along outgoing arcs (successors), `In` the nodes that
//!   reach ego (predecessors), and `All` ignores direction. This is igraph's
//!   `ego()`/`make_ego_graph(mode = ...)`; NetworkX's `ego_graph()` uses
//!   successors, or both directions with `undirected = True`. Undirected
//!   networks ignore the mode.
//! - Distances count steps, so weights are ignored by every function here.
//! - The ego network is the subgraph induced by ego and its alters: every
//!   tie among them counts, whatever its direction and whatever `mode` chose
//!   the alters. Ties keep their direction, so in a directed network an
//!   alter pair joined both ways contributes two ties.
//! - Ties among alters are distinct ties between distinct nodes: parallel
//!   edges count once and self-loops are ignored, the adjacency-matrix
//!   reading that UCINET uses. Ego betweenness is computed on the same
//!   binary ties, so parallel edges are not distinct shortest paths.
//! - Ego betweenness counts each unordered pair of alters once in an
//!   undirected network (as betweenness does) and each ordered pair in a
//!   directed one. It is not normalized; its maximum is the number of pairs.
//! - An ego without alters has size 0, ties 0, pairs 0, components 0, ego
//!   betweenness 0, and density NaN. One alter gives density NaN (no
//!   pairs).
//! - `mindist` in [`neighborhood`] drops nodes closer than `mindist` steps,
//!   as igraph's `mindist`; `mindist = 1` drops ego itself.

use crate::bitset::{for_each_bit, MultiBfs, Step, BATCH};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::closeness::oriented;
use crate::metrics::degree::DegreeMode;
use crate::parallel;
use crate::union_find::UnionFind;

/// Which nodes form a neighborhood.
#[derive(Debug, Clone, Copy)]
pub struct NeighborhoodOptions {
    /// Largest number of steps from ego.
    pub order: usize,
    /// Direction of the steps in directed graphs.
    pub mode: DegreeMode,
    /// Smallest number of steps from ego; 0 keeps ego, 1 drops it.
    pub mindist: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for NeighborhoodOptions {
    fn default() -> Self {
        NeighborhoodOptions {
            order: 1,
            mode: DegreeMode::All,
            mindist: 0,
            threads: 0,
        }
    }
}

/// Options for [`ego_summary`].
#[derive(Debug, Clone, Copy)]
pub struct EgoOptions {
    /// Alters are the nodes within this many steps of ego (at least 1).
    pub order: usize,
    /// Direction of the steps in directed graphs.
    pub mode: DegreeMode,
    /// Whether to compute ego betweenness, the costliest measure.
    pub betweenness: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for EgoOptions {
    fn default() -> Self {
        EgoOptions {
            order: 1,
            mode: DegreeMode::All,
            betweenness: true,
            threads: 0,
        }
    }
}

/// Egocentric measures of every node, in node order. See the module
/// documentation for definitions.
#[derive(Debug, Clone, PartialEq)]
pub struct EgoSummary {
    /// Number of alters.
    pub size: Vec<f64>,
    /// Ties among alters.
    pub ties: Vec<f64>,
    /// Pairs (ordered, if directed) of alters.
    pub pairs: Vec<f64>,
    /// `ties / pairs`, NaN without pairs.
    pub density: Vec<f64>,
    /// Weak components among the alters once ego is removed.
    pub components: Vec<f64>,
    /// Betweenness of ego within its ego network, if requested.
    pub betweenness: Option<Vec<f64>>,
}

fn check_ego(g: &Graph, ego: usize) -> Result<()> {
    if ego >= g.node_count() {
        return Err(Error::InvalidArgument(format!(
            "node index {ego} is out of range for a network with {} nodes",
            g.node_count()
        )));
    }
    Ok(())
}

/// Per-thread breadth-first search state, reused across egos.
struct Bfs {
    /// Steps from the current ego; `u32::MAX` if not reached.
    dist: Vec<u32>,
    /// Reached nodes in order of discovery; doubles as the queue.
    order: Vec<u32>,
}

impl Bfs {
    fn new(n: usize) -> Bfs {
        Bfs {
            dist: vec![u32::MAX; n],
            order: Vec::new(),
        }
    }

    /// Fills `order` with ego and every node within `order_max` steps, in
    /// order of discovery (so in non-decreasing distance).
    fn run(&mut self, g: &Graph, ego: usize, order_max: usize, mode: DegreeMode) {
        self.order.clear();
        self.dist[ego] = 0;
        self.order.push(ego as u32);
        let directed = g.is_directed();
        let limit = order_max.min(u32::MAX as usize - 1) as u32;
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let d = self.dist[v];
            if d >= limit {
                continue;
            }
            let follow_out = !directed || mode != DegreeMode::In;
            let follow_in = directed && mode != DegreeMode::Out;
            if follow_out {
                for &w in g.out().neighbors(v) {
                    if self.dist[w as usize] == u32::MAX {
                        self.dist[w as usize] = d + 1;
                        self.order.push(w);
                    }
                }
            }
            if follow_in {
                for &w in g.inc().neighbors(v) {
                    if self.dist[w as usize] == u32::MAX {
                        self.dist[w as usize] = d + 1;
                        self.order.push(w);
                    }
                }
            }
        }
    }

    fn reset(&mut self) {
        for &v in &self.order {
            self.dist[v as usize] = u32::MAX;
        }
    }
}

/// The nodes within `opts.order` steps of `ego` and at least
/// `opts.mindist` steps away, in ascending node order.
pub fn neighborhood(g: &Graph, ego: usize, opts: NeighborhoodOptions) -> Result<Vec<u32>> {
    check_ego(g, ego)?;
    let mut bfs = Bfs::new(g.node_count());
    bfs.run(g, ego, opts.order, opts.mode);
    let mut out: Vec<u32> = bfs
        .order
        .iter()
        .copied()
        .filter(|&v| bfs.dist[v as usize] as usize >= opts.mindist)
        .collect();
    out.sort_unstable();
    Ok(out)
}

/// [`neighborhood`] for several egos, in parallel, in the order given.
pub fn neighborhoods(
    g: &Graph,
    egos: &[usize],
    opts: NeighborhoodOptions,
) -> Result<Vec<Vec<u32>>> {
    for &e in egos {
        check_ego(g, e)?;
    }
    let n = g.node_count();
    let partials = parallel::for_each_item(
        egos.len(),
        opts.threads,
        || (Bfs::new(n), Vec::new()),
        |(bfs, found): &mut (Bfs, Vec<(usize, Vec<u32>)>), i| {
            bfs.run(g, egos[i], opts.order, opts.mode);
            let mut members: Vec<u32> = bfs
                .order
                .iter()
                .copied()
                .filter(|&v| bfs.dist[v as usize] as usize >= opts.mindist)
                .collect();
            members.sort_unstable();
            bfs.reset();
            found.push((i, members));
        },
    );
    let mut out = vec![Vec::new(); egos.len()];
    for (_, found) in partials {
        for (i, members) in found {
            out[i] = members;
        }
    }
    Ok(out)
}

/// Neighborhood sizes for the given egos (every node if `None`), in
/// parallel. This is igraph's `ego_size()`: with `mindist = 0` ego itself
/// is counted.
pub fn neighborhood_sizes(
    g: &Graph,
    egos: Option<&[usize]>,
    opts: NeighborhoodOptions,
) -> Result<Vec<usize>> {
    let n = g.node_count();
    if let Some(list) = egos {
        for &e in list {
            check_ego(g, e)?;
        }
    }
    let count = egos.map_or(n, |l| l.len());
    let sources: Vec<u32> = match egos {
        Some(l) => l.iter().map(|&e| e as u32).collect(),
        None => (0..n as u32).collect(),
    };
    // Breadth-first search from 64 egos at a time
    // (`crate::bitset::MultiBfs`) over the arcs the mode follows; one
    // batch is one work item, so the result does not depend on the thread
    // count.
    let h = oriented(g, opts.mode)?;
    let adj = h.as_ref();
    let max_level = opts.order.min(u32::MAX as usize - 1) as u32;
    let partials = parallel::for_each_item(
        (count + BATCH - 1) / BATCH,
        opts.threads,
        || (MultiBfs::new(n), Vec::new()),
        |(ms, found): &mut (MultiBfs, Vec<(usize, [usize; BATCH])>), b| {
            let lo = b * BATCH;
            let hi = (lo + BATCH).min(count);
            let mut sizes = [usize::from(opts.mindist == 0); BATCH];
            ms.run(adj, &sources[lo..hi], max_level, |step| {
                if let Step::Reached { level, sources, .. } = step {
                    if level as usize >= opts.mindist {
                        for_each_bit(sources, |i| sizes[i] += 1);
                    }
                }
            });
            found.push((lo, sizes));
        },
    );
    let mut out = vec![0usize; count];
    for (_, found) in partials {
        for (lo, sizes) in found {
            let hi = (lo + BATCH).min(count);
            out[lo..hi].copy_from_slice(&sizes[..hi - lo]);
        }
    }
    Ok(out)
}

/// Per-thread state for [`ego_summary`].
struct EgoWork {
    bfs: Bfs,
    /// Position of each node in the current ego network; `u32::MAX` if
    /// outside it.
    local: Vec<u32>,
    offsets: Vec<usize>,
    adj: Vec<u32>,
    sets: UnionFind,
    hops: Vec<u32>,
    sigma: Vec<f64>,
    delta: Vec<f64>,
    queue: Vec<u32>,
    /// (ego, [size, ties, pairs, density, components, betweenness]).
    found: Vec<(u32, [f64; 6])>,
}

impl EgoWork {
    fn new(n: usize) -> EgoWork {
        EgoWork {
            bfs: Bfs::new(n),
            local: vec![u32::MAX; n],
            offsets: Vec::new(),
            adj: Vec::new(),
            sets: UnionFind::new(0),
            hops: Vec::new(),
            sigma: Vec::new(),
            delta: Vec::new(),
            queue: Vec::new(),
            found: Vec::new(),
        }
    }

    fn measure(&mut self, g: &Graph, ego: usize, opts: &EgoOptions) {
        self.bfs.run(g, ego, opts.order, opts.mode);
        // Ego is local node 0; alters follow in order of discovery.
        let members = &self.bfs.order;
        let s = members.len();
        for (i, &v) in members.iter().enumerate() {
            self.local[v as usize] = i as u32;
        }

        // Binary adjacency of the induced ego network: outgoing arcs (both
        // directions for undirected graphs), distinct, without loops.
        self.offsets.clear();
        self.adj.clear();
        self.offsets.push(0);
        for (a, &v) in members.iter().enumerate() {
            let start = self.adj.len();
            for &w in g.out().neighbors(v as usize) {
                let b = self.local[w as usize];
                if b != u32::MAX && b as usize != a {
                    self.adj.push(b);
                }
            }
            self.adj[start..].sort_unstable();
            let mut keep = start;
            for k in start..self.adj.len() {
                if k == start || self.adj[k] != self.adj[keep - 1] {
                    self.adj[keep] = self.adj[k];
                    keep += 1;
                }
            }
            self.adj.truncate(keep);
            self.offsets.push(self.adj.len());
        }

        let k = (s - 1) as f64;
        let directed = g.is_directed();
        let mut arcs = 0usize;
        for a in 1..s {
            for &b in &self.adj[self.offsets[a]..self.offsets[a + 1]] {
                if b != 0 {
                    arcs += 1;
                }
            }
        }
        let ties = if directed { arcs } else { arcs / 2 } as f64;
        let pairs = if directed {
            k * (k - 1.0)
        } else {
            k * (k - 1.0) / 2.0
        };
        let density = if pairs > 0.0 { ties / pairs } else { f64::NAN };

        // Weak components among alters by union-find (ego excluded).
        self.sets.reset(s);
        let mut components = s.saturating_sub(1);
        for a in 1..s {
            for idx in self.offsets[a]..self.offsets[a + 1] {
                let b = self.adj[idx];
                if b == 0 {
                    continue;
                }
                if self.sets.union(a, b as usize).is_some() {
                    components -= 1;
                }
            }
        }

        let betweenness = if opts.betweenness && s > 2 {
            self.ego_betweenness(s, directed)
        } else {
            0.0
        };

        for &v in &self.bfs.order {
            self.local[v as usize] = u32::MAX;
        }
        self.bfs.reset();
        self.found.push((
            ego as u32,
            [k, ties, pairs, density, components as f64, betweenness],
        ));
    }

    /// Dependency of local node 0 (ego) summed over every other source,
    /// by breadth-first Brandes accumulation over successors.
    fn ego_betweenness(&mut self, s: usize, directed: bool) -> f64 {
        self.hops.clear();
        self.hops.resize(s, u32::MAX);
        self.sigma.clear();
        self.sigma.resize(s, 0.0);
        self.delta.clear();
        self.delta.resize(s, 0.0);
        let mut total = 0.0;
        for src in 1..s {
            self.queue.clear();
            self.queue.push(src as u32);
            self.hops[src] = 0;
            self.sigma[src] = 1.0;
            let mut head = 0;
            while head < self.queue.len() {
                let v = self.queue[head] as usize;
                head += 1;
                let next = self.hops[v] + 1;
                for &w in &self.adj[self.offsets[v]..self.offsets[v + 1]] {
                    let w = w as usize;
                    if self.hops[w] == u32::MAX {
                        self.hops[w] = next;
                        self.queue.push(w as u32);
                    }
                    if self.hops[w] == next {
                        self.sigma[w] += self.sigma[v];
                    }
                }
            }
            for i in (0..self.queue.len()).rev() {
                let v = self.queue[i] as usize;
                let next = self.hops[v] + 1;
                let mut dv = 0.0;
                for &w in &self.adj[self.offsets[v]..self.offsets[v + 1]] {
                    let w = w as usize;
                    if self.hops[w] == next {
                        dv += (1.0 + self.delta[w]) / self.sigma[w];
                    }
                }
                self.delta[v] = self.sigma[v] * dv;
            }
            total += self.delta[0];
            for &v in &self.queue {
                let v = v as usize;
                self.hops[v] = u32::MAX;
                self.sigma[v] = 0.0;
                self.delta[v] = 0.0;
            }
        }
        if directed {
            total
        } else {
            total / 2.0
        }
    }
}

/// Egocentric measures of every node. See the module documentation.
pub fn ego_summary(g: &Graph, opts: EgoOptions) -> Result<EgoSummary> {
    if opts.order == 0 {
        return Err(Error::InvalidArgument(
            "the order of an ego network must be at least 1".into(),
        ));
    }
    let n = g.node_count();
    let partials = parallel::for_each_item(
        n,
        opts.threads,
        || EgoWork::new(n),
        |ws, ego| ws.measure(g, ego, &opts),
    );
    let mut out = EgoSummary {
        size: vec![0.0; n],
        ties: vec![0.0; n],
        pairs: vec![0.0; n],
        density: vec![f64::NAN; n],
        components: vec![0.0; n],
        betweenness: opts.betweenness.then(|| vec![0.0; n]),
    };
    for ws in partials {
        for (ego, m) in ws.found {
            let e = ego as usize;
            out.size[e] = m[0];
            out.ties[e] = m[1];
            out.pairs[e] = m[2];
            out.density[e] = m[3];
            out.components[e] = m[4];
            if let Some(b) = out.betweenness.as_mut() {
                b[e] = m[5];
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(n: usize, s: &[i64], t: &[i64], directed: bool) -> Graph {
        Graph::from_edges(n, s, t, None, directed).unwrap()
    }

    #[test]
    fn neighborhoods_follow_mode_and_order() {
        // 0 -> 1 -> 2 -> 3, 4 -> 0.
        let d = g(5, &[0, 1, 2, 4], &[1, 2, 3, 0], true);
        let o = |order, mode, mindist| NeighborhoodOptions {
            order,
            mode,
            mindist,
            threads: 1,
        };
        assert_eq!(
            neighborhood(&d, 0, o(1, DegreeMode::Out, 0)).unwrap(),
            vec![0, 1]
        );
        assert_eq!(
            neighborhood(&d, 0, o(1, DegreeMode::In, 0)).unwrap(),
            vec![0, 4]
        );
        assert_eq!(
            neighborhood(&d, 0, o(1, DegreeMode::All, 1)).unwrap(),
            vec![1, 4]
        );
        assert_eq!(
            neighborhood(&d, 0, o(2, DegreeMode::Out, 2)).unwrap(),
            vec![2]
        );
        assert_eq!(
            neighborhood_sizes(&d, None, o(2, DegreeMode::All, 0)).unwrap(),
            vec![4, 5, 4, 3, 3]
        );
        let many = neighborhoods(&d, &[3, 0], o(1, DegreeMode::In, 0)).unwrap();
        assert_eq!(many, vec![vec![2, 3], vec![0, 4]]);
        assert!(neighborhood(&d, 5, o(1, DegreeMode::All, 0)).is_err());
    }

    #[test]
    fn kite_like_ego_measures() {
        // Ego 0 with alters 1..=4; ties 1-2 (twice, parallel) and 3-4,
        // plus a loop at 1 and an outside node 5 tied to 4.
        let u = g(
            6,
            &[0, 0, 0, 0, 1, 1, 3, 1, 4],
            &[1, 2, 3, 4, 2, 2, 4, 1, 5],
            false,
        );
        let s = ego_summary(&u, EgoOptions::default()).unwrap();
        assert_eq!(s.size[0], 4.0);
        assert_eq!(s.ties[0], 2.0);
        assert_eq!(s.pairs[0], 6.0);
        assert_eq!(s.density[0], 2.0 / 6.0);
        assert_eq!(s.components[0], 2.0);
        // Ego brokers the four pairs across the two dyads.
        assert_eq!(s.betweenness.as_ref().unwrap()[0], 4.0);
        // Node 5 has one alter: no pairs.
        assert_eq!(s.size[5], 1.0);
        assert!(s.density[5].is_nan());
        assert_eq!(s.components[5], 1.0);
    }

    #[test]
    fn star_ego_betweenness_is_the_number_of_pairs() {
        let star = g(5, &[0, 0, 0, 0], &[1, 2, 3, 4], false);
        let s = ego_summary(&star, EgoOptions::default()).unwrap();
        assert_eq!(s.betweenness.unwrap(), vec![6.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(s.components[0], 4.0);
        let d = g(3, &[1, 0], &[0, 2], true);
        let s = ego_summary(&d, EgoOptions::default()).unwrap();
        // Only 1 -> 0 -> 2 passes through 0.
        assert_eq!(s.betweenness.as_ref().unwrap()[0], 1.0);
        assert_eq!(s.pairs[0], 2.0);
        let out = EgoOptions {
            mode: DegreeMode::Out,
            ..Default::default()
        };
        let s = ego_summary(&d, out).unwrap();
        assert_eq!(s.size, vec![1.0, 1.0, 0.0]);
        assert!(s.density[2].is_nan());
        assert_eq!(s.components[2], 0.0);
    }

    #[test]
    fn rejects_order_zero() {
        let u = g(2, &[0], &[1], false);
        let bad = EgoOptions {
            order: 0,
            ..Default::default()
        };
        assert!(ego_summary(&u, bad).is_err());
    }
}
