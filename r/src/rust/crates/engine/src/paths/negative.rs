//! Shortest paths with negative lengths, negative cycles, and Johnson's
//! all-pairs algorithm.
//!
//! - [`bellman_ford`]: the Bellman-Ford-Moore algorithm (Bellman 1958,
//!   `bellman1958routing`, with Ford's and Moore's label-correcting
//!   formulation) with a first-in first-out queue of nodes whose distance
//!   changed, as described by Cherkassky, Goldberg, and Radzik (1996,
//!   `cherkassky1996shortest`). Nodes enter the queue in passes; a node queued in
//!   pass n (for n nodes) proves a negative cycle, which is then read off
//!   the parent pointers (every cycle of parent pointers built by
//!   relaxations is negative; Cherkassky and Goldberg 1999,
//!   `cherkassky1999negative`).
//! - [`goldberg_radzik`]: Goldberg and Radzik's algorithm
//!   (`goldberg1993heuristic`). Each pass takes the nodes whose distance
//!   changed and have an arc of negative reduced cost
//!   d(u) + len(u, v) - d(v), finds every node reachable from them along
//!   such arcs, orders those nodes topologically by depth-first search, and
//!   scans them in that order. A cycle of arcs with negative reduced cost
//!   found by the depth-first search is a negative cycle. It needs O(nm)
//!   time like Bellman-Ford, and usually far fewer scans.
//! - [`johnson`]: Johnson's algorithm (`johnson1977efficient`): potentials
//!   h from Bellman-Ford from a virtual source joined to every node by arcs
//!   of length 0 make every reduced length len(u, v) + h(u) - h(v)
//!   non-negative, and Dijkstra's algorithm then runs from every source in
//!   parallel. Distances are summed along the resulting shortest-path
//!   trees in the original lengths, so they are as exact as Dijkstra's
//!   (integer lengths are exact).
//!
//! Conventions:
//! - Only negative cycles reachable from the sources are detected, as in
//!   NetworkX (`bellman_ford_*`, `goldberg_radzik`, `find_negative_cycle`).
//!   [`has_negative_cycle`] looks for one anywhere (NetworkX's
//!   `negative_edge_cycle`), as does [`johnson`].
//! - In an undirected network an edge of negative length is a negative
//!   cycle (walk it back and forth), as in NetworkX. A self-loop of
//!   negative length is a negative cycle; other self-loops are ignored.
//! - A negative cycle is returned as a closed walk following the arcs
//!   (first node repeated at the end), as NetworkX's `find_negative_cycle`.
//! - Parents and the order of predecessor lists follow
//!   [`crate::paths::search`]: nodes are taken in order of final distance,
//!   ties by index.

use crate::error::{Error, Result};
use crate::parallel;
use crate::paths::search::{fill_rows, predecessor_lists, SearchOptions, Tree};
use crate::paths::{NodeLists, PathGraph, Tie, NO_NODE};

/// A value, or a negative cycle as a closed walk.
pub type OrCycle<T> = std::result::Result<T, Vec<u32>>;

/// The result of a search that allows negative lengths.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// Shortest paths (see [`Tree`]).
    Tree(Tree),
    /// A negative cycle reachable from the sources, as a closed walk.
    NegativeCycle(Vec<u32>),
}

impl Outcome {
    /// The tree, or an error naming a negative cycle.
    pub fn tree(self) -> Result<Tree> {
        match self {
            Outcome::Tree(t) => Ok(t),
            Outcome::NegativeCycle(c) => Err(negative_cycle_error(&c)),
        }
    }
}

fn negative_cycle_error(cycle: &[u32]) -> Error {
    let shown: Vec<String> = cycle.iter().take(12).map(|v| v.to_string()).collect();
    Error::InvalidArgument(format!(
        "the network has a negative cycle ({}{}), so shortest paths are unbounded",
        shown.join(" -> "),
        if cycle.len() > 12 { " -> ..." } else { "" }
    ))
}

/// State shared by the label-correcting searches.
struct Labels {
    dist: Vec<f64>,
    parent: Vec<u32>,
    parent_arc: Vec<u32>,
    origin: Vec<u32>,
}

impl Labels {
    fn new(n: usize, sources: &[u32]) -> Labels {
        let mut l = Labels {
            dist: vec![f64::INFINITY; n],
            parent: vec![NO_NODE; n],
            parent_arc: vec![NO_NODE; n],
            origin: vec![NO_NODE; n],
        };
        for &s in sources {
            l.dist[s as usize] = 0.0;
            l.origin[s as usize] = s;
        }
        l
    }

    /// A cycle of parent pointers whose arcs have negative total length,
    /// as a closed walk in arc direction, if there is one.
    fn parent_cycle(&self, pg: &PathGraph) -> Option<Vec<u32>> {
        let n = self.parent.len();
        let mut color = vec![0u8; n];
        for start in 0..n {
            if color[start] != 0 {
                continue;
            }
            let mut walk = Vec::new();
            let mut v = start;
            while color[v] == 0 {
                color[v] = 1;
                walk.push(v);
                if self.parent[v] == NO_NODE {
                    break;
                }
                v = self.parent[v] as usize;
            }
            if color[v] == 1 && self.parent[v] != NO_NODE {
                // v is on the current walk: the cycle is v's segment.
                let at = walk.iter().position(|&x| x == v).expect("on walk");
                let mut cycle: Vec<u32> = walk[at..].iter().map(|&x| x as u32).collect();
                let total: f64 = cycle
                    .iter()
                    .map(|&x| pg.len(self.parent_arc[x as usize] as usize))
                    .sum();
                for &x in &walk {
                    color[x] = 2;
                }
                if total < 0.0 {
                    // Parent pointers run against the arcs.
                    cycle.reverse();
                    cycle.push(cycle[0]);
                    return Some(cycle);
                }
                continue;
            }
            for &x in &walk {
                color[x] = 2;
            }
        }
        None
    }

    fn into_tree(self, pg: &PathGraph) -> Tree {
        let n = self.dist.len();
        let mut order: Vec<u32> = (0..n as u32)
            .filter(|&v| self.dist[v as usize].is_finite())
            .collect();
        order.sort_by(|&a, &b| {
            self.dist[a as usize]
                .total_cmp(&self.dist[b as usize])
                .then(a.cmp(&b))
        });
        let parent_edge = self
            .parent_arc
            .iter()
            .map(|&a| {
                if a == NO_NODE {
                    NO_NODE
                } else {
                    pg.edge(a as usize)
                }
            })
            .collect();
        Tree {
            dist: self.dist,
            parent: self.parent,
            parent_edge,
            origin: self.origin,
            order,
        }
    }
}

fn check_sources(pg: &PathGraph, sources: &[u32]) -> Result<()> {
    for &s in sources {
        pg.check_node(s as usize, "source")?;
    }
    Ok(())
}

/// Shortest paths from the nearest of `sources` with lengths of any sign,
/// by the Bellman-Ford-Moore algorithm, or a negative cycle reachable from
/// them. See the module documentation.
pub fn bellman_ford(pg: &PathGraph, sources: &[u32]) -> Result<Outcome> {
    check_sources(pg, sources)?;
    let n = pg.node_count();
    let mut l = Labels::new(n, sources);
    let mut queue = std::collections::VecDeque::new();
    let mut pass = vec![u32::MAX; n];
    for &s in sources {
        if pass[s as usize] == u32::MAX {
            pass[s as usize] = 0;
            queue.push_back(s);
        }
    }
    // Once a node is queued in pass n, look for a cycle of parent
    // pointers after every further n relaxations until one appears.
    let mut suspicious = false;
    let mut since_check = 0usize;
    while let Some(u) = queue.pop_front() {
        let u = u as usize;
        let p = pass[u];
        pass[u] = u32::MAX;
        let du = l.dist[u];
        for a in pg.out_arcs(u) {
            let v = pg.head(a) as usize;
            let alt = du + pg.len(a);
            if alt < l.dist[v] {
                l.dist[v] = alt;
                l.parent[v] = u as u32;
                l.parent_arc[v] = a as u32;
                l.origin[v] = l.origin[u];
                if v == u {
                    // A negative self-loop.
                    return Ok(Outcome::NegativeCycle(vec![u as u32, u as u32]));
                }
                if pass[v] == u32::MAX {
                    pass[v] = p.saturating_add(1);
                    queue.push_back(v as u32);
                    if pass[v] as usize >= n {
                        suspicious = true;
                    }
                }
                since_check += 1;
            }
        }
        if suspicious && since_check >= n.max(1) {
            since_check = 0;
            if let Some(c) = l.parent_cycle(pg) {
                return Ok(Outcome::NegativeCycle(c));
            }
        }
    }
    Ok(Outcome::Tree(l.into_tree(pg)))
}

/// Predecessor lists (every predecessor on a shortest path, see
/// [`crate::paths::search`]) and distances from the nearest of `sources`
/// with lengths of any sign, or a negative cycle.
pub fn bellman_ford_predecessors(
    pg: &PathGraph,
    sources: &[u32],
    tie: Tie,
) -> Result<OrCycle<(NodeLists, Vec<f64>)>> {
    match bellman_ford(pg, sources)? {
        Outcome::NegativeCycle(c) => Ok(Err(c)),
        Outcome::Tree(t) => {
            let lists = tree_predecessors(pg, &t, tie);
            Ok(Ok((lists, t.dist)))
        }
    }
}

/// Predecessor lists (see [`crate::paths::search`]) of a tree of shortest
/// paths, taking nodes in its order.
pub fn tree_predecessors(pg: &PathGraph, t: &Tree, tie: Tie) -> NodeLists {
    let mut rank = vec![u32::MAX; pg.node_count()];
    for (i, &v) in t.order.iter().enumerate() {
        rank[v as usize] = i as u32;
    }
    predecessor_lists(pg, &t.dist, &t.order, &rank, tie, false)
}

/// Shortest paths from `source` with lengths of any sign by Goldberg and
/// Radzik's algorithm, or a negative cycle reachable from it.
pub fn goldberg_radzik(pg: &PathGraph, source: u32) -> Result<Outcome> {
    check_sources(pg, &[source])?;
    let n = pg.node_count();
    let mut l = Labels::new(n, &[source]);
    // Negative self-loops are negative cycles of one arc.
    let mut labeled = vec![source];
    let mut is_labeled = vec![false; n];
    // DFS state: 0 unvisited, 1 on the stack, 2 finished (per pass).
    let mut state = vec![0u8; n];
    let mut stamp = vec![0u32; n];
    let mut topo: Vec<u32> = Vec::new();
    let mut stack: Vec<(u32, usize)> = Vec::new();
    let mut passes = 0usize;
    while !labeled.is_empty() {
        passes += 1;
        if passes > n + 1 {
            if let Some(c) = l.parent_cycle(pg) {
                return Ok(Outcome::NegativeCycle(c));
            }
        }
        let pass = passes as u32;
        topo.clear();
        for &b in &labeled {
            is_labeled[b as usize] = false;
        }
        // Depth-first search along arcs of negative reduced cost from
        // every labeled node that has one.
        for &b in &labeled {
            let b = b as usize;
            if stamp[b] == pass && state[b] != 0 {
                continue;
            }
            let db = l.dist[b];
            let has_negative = pg
                .out_arcs(b)
                .any(|a| db + pg.len(a) < l.dist[pg.head(a) as usize]);
            if !has_negative {
                continue;
            }
            stamp[b] = pass;
            state[b] = 1;
            stack.push((b as u32, pg.out_arcs(b).start));
            while let Some(top) = stack.len().checked_sub(1) {
                let (u, mut next) = stack[top];
                let u = u as usize;
                let end = pg.out_arcs(u).end;
                let du = l.dist[u];
                let mut child = None;
                while next < end {
                    let a = next;
                    next += 1;
                    let v = pg.head(a) as usize;
                    if du + pg.len(a) >= l.dist[v] {
                        continue;
                    }
                    if stamp[v] != pass {
                        stamp[v] = pass;
                        state[v] = 0;
                    }
                    match state[v] {
                        0 => {
                            child = Some(v);
                            break;
                        }
                        1 => {
                            // A cycle of negative reduced costs: the stack
                            // from v to u, closed by the arc u -> v.
                            let at = stack
                                .iter()
                                .position(|&(x, _)| x as usize == v)
                                .expect("on the stack");
                            let mut cycle: Vec<u32> = stack[at..].iter().map(|&(x, _)| x).collect();
                            cycle.push(v as u32);
                            return Ok(Outcome::NegativeCycle(cycle));
                        }
                        _ => {}
                    }
                }
                stack[top].1 = next;
                match child {
                    Some(v) => {
                        state[v] = 1;
                        stack.push((v as u32, pg.out_arcs(v).start));
                    }
                    None => {
                        state[u] = 2;
                        topo.push(u as u32);
                        stack.pop();
                    }
                }
            }
        }
        // Scan in topological order (reverse post-order).
        labeled.clear();
        for &u in topo.iter().rev() {
            let u = u as usize;
            let du = l.dist[u];
            for a in pg.out_arcs(u) {
                let v = pg.head(a) as usize;
                let alt = du + pg.len(a);
                if alt < l.dist[v] {
                    if v == u {
                        return Ok(Outcome::NegativeCycle(vec![u as u32, u as u32]));
                    }
                    l.dist[v] = alt;
                    l.parent[v] = u as u32;
                    l.parent_arc[v] = a as u32;
                    l.origin[v] = l.origin[u];
                    if !is_labeled[v] {
                        is_labeled[v] = true;
                        labeled.push(v as u32);
                    }
                }
            }
        }
        labeled.sort_unstable();
    }
    Ok(Outcome::Tree(l.into_tree(pg)))
}

/// A negative cycle reachable from `source`, as a closed walk, or `None`.
pub fn find_negative_cycle(pg: &PathGraph, source: u32) -> Result<Option<Vec<u32>>> {
    Ok(match bellman_ford(pg, &[source])? {
        Outcome::NegativeCycle(c) => Some(c),
        Outcome::Tree(_) => None,
    })
}

/// A negative cycle anywhere in the network (Bellman-Ford from every node
/// at once), as a closed walk, or `None`.
pub fn has_negative_cycle(pg: &PathGraph) -> Result<Option<Vec<u32>>> {
    let all: Vec<u32> = (0..pg.node_count() as u32).collect();
    Ok(match bellman_ford(pg, &all)? {
        Outcome::NegativeCycle(c) => Some(c),
        Outcome::Tree(_) => None,
    })
}

/// Distances (row-major, one row per source) and parents ([`NO_NODE`] for
/// sources and unreached nodes) from every source in `sources` (every node
/// when `None`) by Johnson's algorithm, or a negative cycle. Rows are
/// independent, so the result does not depend on the thread count.
pub fn johnson(
    pg: &PathGraph,
    sources: Option<&[u32]>,
    threads: usize,
) -> Result<OrCycle<(Vec<f64>, Vec<u32>)>> {
    let n = pg.node_count();
    let all: Vec<u32> = (0..n as u32).collect();
    let sources = sources.unwrap_or(&all);
    check_sources(pg, sources)?;
    let h = match bellman_ford(pg, &all)? {
        Outcome::NegativeCycle(c) => return Ok(Err(c)),
        Outcome::Tree(t) => t.dist,
    };
    let reduced: Vec<f64> = (0..pg.arc_count())
        .map(|a| {
            let (u, v) = (pg.tail(a) as usize, pg.head(a) as usize);
            (pg.len(a) + h[u] - h[v]).max(0.0)
        })
        .collect();
    let rg = pg.with_lengths(reduced);
    let mut dist = vec![f64::INFINITY; sources.len() * n];
    let mut parent = vec![NO_NODE; sources.len() * n];
    let threads = parallel::threads_for(
        threads,
        sources.len(),
        parallel::search_work(sources.len(), n, pg.arc_count()),
    );
    fill_rows(&mut dist, &mut parent, n, threads, |ws, k, drow, prow| {
        // The reduced lengths are non-negative and the sources valid.
        let _ = ws.run(&rg, &sources[k..k + 1], SearchOptions::default());
        for &v in ws.order() {
            let v = v as usize;
            let p = ws.parent(v);
            prow[v] = p;
            drow[v] = if p == NO_NODE {
                0.0
            } else {
                drow[p as usize] + pg.len(ws.parent_arc(v) as usize)
            };
        }
    });
    Ok(Ok((dist, parent)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::paths::search::{all_pairs, single_source};
    use crate::paths::{ArcLengths, PathGraphOptions};
    use crate::rng::Stream;

    fn pg(n: usize, s: &[i64], t: &[i64], len: &[f64], directed: bool) -> PathGraph {
        let g = Graph::from_edges(n, s, t, None, directed).unwrap();
        PathGraph::new(
            &g,
            PathGraphOptions {
                lengths: ArcLengths::PerEdge(len),
                ..Default::default()
            },
        )
        .unwrap()
    }

    /// Exact distances by n - 1 rounds of relaxing every arc, or None if a
    /// round n still relaxes (a negative cycle reachable from s).
    fn naive(p: &PathGraph, s: usize) -> Option<Vec<f64>> {
        let n = p.node_count();
        let mut d = vec![f64::INFINITY; n];
        d[s] = 0.0;
        for round in 0..=n {
            let mut changed = false;
            for a in 0..p.arc_count() {
                let (u, v) = (p.tail(a) as usize, p.head(a) as usize);
                if d[u] + p.len(a) < d[v] {
                    d[v] = d[u] + p.len(a);
                    changed = true;
                }
            }
            if !changed {
                return Some(d);
            }
            if round == n {
                return None;
            }
        }
        unreachable!()
    }

    fn check_cycle(p: &PathGraph, c: &[u32]) {
        assert_eq!(c.first(), c.last());
        assert!(c.len() >= 2);
        let total: f64 = c
            .windows(2)
            .map(|w| p.len(p.arc_between(w[0] as usize, w[1] as usize).unwrap()))
            .sum();
        assert!(total < 0.0, "{c:?} has length {total}");
    }

    #[test]
    fn agree_with_naive_relaxation() {
        for seed in 0..200u64 {
            let mut rng = Stream::new(seed, 3);
            let n = 2 + rng.below(9) as usize;
            let m = rng.below(3 * n as u64) as usize;
            let s: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            let t: Vec<i64> = (0..m).map(|_| rng.below(n as u64) as i64).collect();
            // Mostly positive lengths with some negative ones.
            let len: Vec<f64> = (0..m).map(|_| rng.below(12) as f64 - 2.0).collect();
            let p = pg(n, &s, &t, &len, true);
            for src in 0..n {
                let want = naive(&p, src);
                let bf = bellman_ford(&p, &[src as u32]).unwrap();
                let gr = goldberg_radzik(&p, src as u32).unwrap();
                for out in [bf, gr] {
                    match (&want, out) {
                        (Some(d), Outcome::Tree(t)) => {
                            assert_eq!(&t.dist, d, "seed {seed} source {src}");
                            for (v, &dv) in d.iter().enumerate() {
                                if let Some(path) = t.path_to(v) {
                                    assert_eq!(p.path_length(&path), Some(dv));
                                }
                            }
                        }
                        (None, Outcome::NegativeCycle(c)) => check_cycle(&p, &c),
                        (w, o) => panic!("seed {seed} source {src}: {w:?} vs {o:?}"),
                    }
                }
            }
            match johnson(&p, None, 2).unwrap() {
                Ok((d, _)) => {
                    for src in 0..n {
                        assert_eq!(Some(d[src * n..(src + 1) * n].to_vec()), naive(&p, src));
                    }
                }
                Err(c) => check_cycle(&p, &c),
            }
        }
    }

    #[test]
    fn undirected_negative_edge_is_a_cycle_and_johnson_matches_dijkstra() {
        let p = pg(3, &[0, 1], &[1, 2], &[1.0, -1.0], false);
        let c = find_negative_cycle(&p, 0).unwrap().unwrap();
        check_cycle(&p, &c);
        let p = pg(4, &[0, 1, 2, 0], &[1, 2, 3, 3], &[1.0, 2.0, 0.5, 4.0], true);
        let (d, _) = johnson(&p, None, 1).unwrap().unwrap();
        let (d2, _) = all_pairs(&p, None, f64::INFINITY, 1).unwrap();
        assert_eq!(d, d2);
        let t = single_source(&p, 0, SearchOptions::default()).unwrap();
        assert_eq!(bellman_ford(&p, &[0]).unwrap().tree().unwrap().dist, t.dist);
        assert!(has_negative_cycle(&p).unwrap().is_none());
        let p = pg(2, &[0, 1], &[1, 1], &[1.0, -0.5], true);
        assert_eq!(has_negative_cycle(&p).unwrap(), Some(vec![1, 1]));
    }
}
