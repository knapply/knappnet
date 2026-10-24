//! Single-source, multi-source, and all-pairs shortest paths with
//! parents, and the predecessor lists of the shortest-path DAG.
//!
//! Unit lengths use breadth-first search; other non-negative lengths use
//! Dijkstra's algorithm (`dijkstra1959note`) with the indexed 4-ary heap of
//! [`crate::weights`] (`johnson1975priority`). A multi-source search starts
//! every source at distance 0, so each node gets the distance to its
//! nearest source and the source it was reached from (the cells of a
//! graph Voronoi diagram, `erwig2000graph`).
//!
//! Conventions (see also [`crate::paths`]):
//! - Lengths must be non-negative; zero lengths are allowed, as in
//!   NetworkX's `dijkstra_*` functions.
//! - `cutoff`: only nodes at distance at most `cutoff` are reached (an arc
//!   whose head would be farther is not followed), as NetworkX's `cutoff`.
//! - `target`: the search stops once the target's distance is final. Other
//!   nodes whose distance was final by then are reported; the rest are
//!   reported as unreached, as NetworkX does when given a target.
//! - `reverse`: follow arcs backward, giving distances *to* the sources
//!   (NetworkX's `single_target_*` functions).
//! - The parent of a node is the first node (in settling order) whose arc
//!   reached the node's final distance; the source of a multi-source search
//!   is inherited from the parent. Among tied nodes the smaller index is
//!   settled first, so results do not depend on threads or timing.
//! - Predecessor lists ([`predecessors`]) hold every node `u` with an arc
//!   `u -> v` such that d(u) + len = d(v) under the tie rule, in the order
//!   the predecessors were settled, each once however many parallel arcs
//!   join them, and never `v` itself. With zero-length arcs two nodes at
//!   the same distance can be each other's predecessors, as in NetworkX.

use crate::error::{Error, Result};
use crate::parallel;
use crate::paths::{NodeLists, PathGraph, Tie, NO_NODE};
use crate::weights::IndexedHeap;

/// Limits of a search.
#[derive(Debug, Clone, Copy)]
pub struct SearchOptions {
    /// Largest distance reached (infinity for no limit).
    pub cutoff: f64,
    /// Stop once this node's distance is final.
    pub target: Option<u32>,
    /// Follow arcs backward.
    pub reverse: bool,
}

impl Default for SearchOptions {
    fn default() -> Self {
        SearchOptions {
            cutoff: f64::INFINITY,
            target: None,
            reverse: false,
        }
    }
}

/// The result of a search: for every node its distance (infinity when not
/// reached), parent and the edge to it from the parent ([`NO_NODE`] for
/// sources and unreached nodes), and the source it was reached from; and
/// the reached nodes in the order their distances became final.
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    /// Distance from the nearest source.
    pub dist: Vec<f64>,
    /// Parent on a shortest path, or [`NO_NODE`].
    pub parent: Vec<u32>,
    /// Edge (input index) from the parent, or [`NO_NODE`].
    pub parent_edge: Vec<u32>,
    /// The source each node was reached from, or [`NO_NODE`].
    pub origin: Vec<u32>,
    /// Reached nodes in settling order.
    pub order: Vec<u32>,
}

impl Tree {
    /// The path from the node's source to `v`, or `None` if `v` was not
    /// reached.
    pub fn path_to(&self, v: usize) -> Option<Vec<u32>> {
        if self.dist[v].is_infinite() && self.origin[v] == NO_NODE {
            return None;
        }
        let mut path = vec![v as u32];
        let mut cur = v;
        while self.parent[cur] != NO_NODE {
            cur = self.parent[cur] as usize;
            path.push(cur as u32);
        }
        path.reverse();
        Some(path)
    }
}

/// Reusable scratch space for searches on graphs with `n` nodes.
#[derive(Debug, Clone)]
pub struct Searcher {
    dist: Vec<f64>,
    parent: Vec<u32>,
    parent_arc: Vec<u32>,
    origin: Vec<u32>,
    /// Position in `order` once the node's distance is final.
    rank: Vec<u32>,
    order: Vec<u32>,
    touched: Vec<u32>,
    heap: IndexedHeap,
}

impl Searcher {
    /// Scratch space for graphs with `n` nodes.
    pub fn new(n: usize) -> Searcher {
        Searcher {
            dist: vec![f64::INFINITY; n],
            parent: vec![NO_NODE; n],
            parent_arc: vec![NO_NODE; n],
            origin: vec![NO_NODE; n],
            rank: vec![u32::MAX; n],
            order: Vec::new(),
            touched: Vec::new(),
            heap: IndexedHeap::new(n),
        }
    }

    fn reset(&mut self) {
        for &v in &self.touched {
            let v = v as usize;
            self.dist[v] = f64::INFINITY;
            self.parent[v] = NO_NODE;
            self.parent_arc[v] = NO_NODE;
            self.origin[v] = NO_NODE;
            self.rank[v] = u32::MAX;
        }
        self.touched.clear();
        self.order.clear();
        self.heap.clear();
    }

    fn check(pg: &PathGraph, sources: &[u32], opts: &SearchOptions) -> Result<()> {
        for &s in sources {
            pg.check_node(s as usize, "source")?;
        }
        if let Some(t) = opts.target {
            pg.check_node(t as usize, "target")?;
        }
        if opts.cutoff.is_nan() {
            return Err(Error::InvalidArgument("the cutoff must not be NaN".into()));
        }
        Ok(())
    }

    /// Runs a search from `sources` on `pg`: breadth-first search for unit
    /// lengths, Dijkstra's algorithm otherwise. Lengths must be
    /// non-negative.
    pub fn run(&mut self, pg: &PathGraph, sources: &[u32], opts: SearchOptions) -> Result<()> {
        Self::check(pg, sources, &opts)?;
        if pg.is_unit() {
            self.bfs(pg, sources, opts);
            return Ok(());
        }
        pg.require_nonnegative("Dijkstra's algorithm")?;
        self.dijkstra(pg, sources, opts, |a, d| d + pg.len(a));
        Ok(())
    }

    /// Dijkstra's algorithm where following arc `a` from a node whose label
    /// is final at `d` reaches its other end at `arrival(a, d)`. The labels
    /// are correct when `arrival(a, d) >= d` and `arrival` is non-decreasing
    /// in `d` (the FIFO property of time-dependent networks,
    /// `dreyfus1969appraisal`); with `d + len(a)` this is the ordinary
    /// search. `opts.cutoff` bounds the labels.
    pub fn run_with<F>(
        &mut self,
        pg: &PathGraph,
        sources: &[u32],
        opts: SearchOptions,
        arrival: F,
    ) -> Result<()>
    where
        F: Fn(usize, f64) -> f64,
    {
        Self::check(pg, sources, &opts)?;
        self.dijkstra(pg, sources, opts, arrival);
        Ok(())
    }

    fn bfs(&mut self, pg: &PathGraph, sources: &[u32], opts: SearchOptions) {
        self.reset();
        for &s in sources {
            let s = s as usize;
            if self.rank[s] != u32::MAX {
                continue;
            }
            self.dist[s] = 0.0;
            self.origin[s] = s as u32;
            self.rank[s] = self.order.len() as u32;
            self.order.push(s as u32);
            self.touched.push(s as u32);
        }
        if let Some(t) = opts.target {
            if self.rank[t as usize] != u32::MAX {
                return;
            }
        }
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let next = self.dist[v] + 1.0;
            if next > opts.cutoff {
                continue;
            }
            let (out, inn) = if opts.reverse {
                (0..0, pg.in_arcs(v))
            } else {
                (pg.out_arcs(v), &[][..])
            };
            for a in out.chain(inn.iter().map(|&a| a as usize)) {
                let w = if opts.reverse { pg.tail(a) } else { pg.head(a) } as usize;
                if self.rank[w] != u32::MAX {
                    continue;
                }
                self.dist[w] = next;
                self.parent[w] = v as u32;
                self.parent_arc[w] = a as u32;
                self.origin[w] = self.origin[v];
                self.rank[w] = self.order.len() as u32;
                self.order.push(w as u32);
                self.touched.push(w as u32);
                if opts.target == Some(w as u32) {
                    return;
                }
            }
        }
    }

    fn dijkstra<F>(&mut self, pg: &PathGraph, sources: &[u32], opts: SearchOptions, arrival: F)
    where
        F: Fn(usize, f64) -> f64,
    {
        self.reset();
        for &s in sources {
            let s = s as usize;
            if self.origin[s] != NO_NODE {
                continue;
            }
            self.dist[s] = 0.0;
            self.origin[s] = s as u32;
            self.touched.push(s as u32);
            self.heap.push_or_decrease(s as u32, 0.0);
        }
        while let Some((d, node)) = self.heap.pop() {
            let v = node as usize;
            self.rank[v] = self.order.len() as u32;
            self.order.push(node);
            if opts.target == Some(node) {
                break;
            }
            let (out, inn) = if opts.reverse {
                (0..0, pg.in_arcs(v))
            } else {
                (pg.out_arcs(v), &[][..])
            };
            for a in out.chain(inn.iter().map(|&a| a as usize)) {
                let w = if opts.reverse { pg.tail(a) } else { pg.head(a) } as usize;
                if self.rank[w] != u32::MAX {
                    continue;
                }
                let alt = arrival(a, d);
                if alt > opts.cutoff || alt >= self.dist[w] {
                    continue;
                }
                if self.origin[w] == NO_NODE && self.dist[w].is_infinite() {
                    self.touched.push(w as u32);
                }
                self.dist[w] = alt;
                self.parent[w] = v as u32;
                self.parent_arc[w] = a as u32;
                self.origin[w] = self.origin[v];
                self.heap.push_or_decrease(w as u32, alt);
            }
        }
        // Nodes still waiting when a target stopped the search are not final.
        if !self.heap.entries().is_empty() {
            for &(_, w) in self.heap.entries() {
                let w = w as usize;
                self.dist[w] = f64::INFINITY;
                self.parent[w] = NO_NODE;
                self.parent_arc[w] = NO_NODE;
                self.origin[w] = NO_NODE;
            }
            self.heap.clear();
        }
    }

    /// Whether `v`'s distance is final.
    #[inline]
    pub fn reached(&self, v: usize) -> bool {
        self.rank[v] != u32::MAX
    }

    /// Distance of `v` (infinity when not reached).
    #[inline]
    pub fn dist(&self, v: usize) -> f64 {
        if self.reached(v) {
            self.dist[v]
        } else {
            f64::INFINITY
        }
    }

    /// Parent of `v`, or [`NO_NODE`].
    #[inline]
    pub fn parent(&self, v: usize) -> u32 {
        if self.reached(v) {
            self.parent[v]
        } else {
            NO_NODE
        }
    }

    /// Arc from the parent of `v`, or [`NO_NODE`].
    #[inline]
    pub fn parent_arc(&self, v: usize) -> u32 {
        if self.reached(v) {
            self.parent_arc[v]
        } else {
            NO_NODE
        }
    }

    /// Reached nodes in settling order.
    #[inline]
    pub fn order(&self) -> &[u32] {
        &self.order
    }

    /// Position of `v` in [`Searcher::order`] (`u32::MAX` if not reached).
    #[inline]
    pub fn rank(&self, v: usize) -> u32 {
        self.rank[v]
    }

    /// The path from the source to `v`, if `v` was reached.
    pub fn path_to(&self, v: usize) -> Option<Vec<u32>> {
        if !self.reached(v) {
            return None;
        }
        let mut path = vec![v as u32];
        let mut cur = v;
        while self.parent[cur] != NO_NODE {
            cur = self.parent[cur] as usize;
            path.push(cur as u32);
        }
        path.reverse();
        Some(path)
    }

    /// The result of the last search, with edges for parents.
    pub fn tree(&self, pg: &PathGraph) -> Tree {
        let n = self.dist.len();
        let mut t = Tree {
            dist: vec![f64::INFINITY; n],
            parent: vec![NO_NODE; n],
            parent_edge: vec![NO_NODE; n],
            origin: vec![NO_NODE; n],
            order: self.order.clone(),
        };
        for &v in &self.order {
            let v = v as usize;
            t.dist[v] = self.dist[v];
            t.parent[v] = self.parent[v];
            t.origin[v] = self.origin[v];
            if self.parent_arc[v] != NO_NODE {
                t.parent_edge[v] = pg.edge(self.parent_arc[v] as usize);
            }
        }
        t
    }

    /// Predecessor lists of the last search (see the module
    /// documentation), one list per node.
    pub fn predecessors(&self, pg: &PathGraph, tie: Tie, reverse: bool) -> NodeLists {
        predecessor_lists(pg, &self.dist, &self.order, &self.rank, tie, reverse)
    }
}

/// Predecessor lists given final distances `dist` of the nodes in `order`
/// (the order distances became final; `rank[v]` is v's position in it,
/// `u32::MAX` for nodes not reached). See the module documentation.
pub(crate) fn predecessor_lists(
    pg: &PathGraph,
    dist: &[f64],
    order: &[u32],
    rank: &[u32],
    tie: Tie,
    reverse: bool,
) -> NodeLists {
    let n = pg.node_count();
    let mut lists: Vec<Vec<u32>> = vec![Vec::new(); n];
    for &v in order {
        let v = v as usize;
        let (out, inn) = if reverse {
            (0..0, pg.in_arcs(v))
        } else {
            (pg.out_arcs(v), &[][..])
        };
        for a in out.chain(inn.iter().map(|&a| a as usize)) {
            let w = if reverse { pg.tail(a) } else { pg.head(a) } as usize;
            if w == v || rank[w] == u32::MAX {
                continue;
            }
            let alt = dist[v] + pg.len(a);
            let dw = dist[w];
            if alt == dw || (tie.equal(alt, dw) && rank[v] < rank[w]) {
                let list = &mut lists[w];
                if list.last() != Some(&(v as u32)) {
                    list.push(v as u32);
                }
            }
        }
    }
    let mut out = NodeLists::new();
    for list in &lists {
        out.push(list);
    }
    out
}

/// Shortest paths from one source. See the module documentation.
pub fn single_source(pg: &PathGraph, source: u32, opts: SearchOptions) -> Result<Tree> {
    multi_source(pg, &[source], opts)
}

/// Shortest paths from the nearest of several sources. See the module
/// documentation.
pub fn multi_source(pg: &PathGraph, sources: &[u32], opts: SearchOptions) -> Result<Tree> {
    let mut ws = Searcher::new(pg.node_count());
    ws.run(pg, sources, opts)?;
    Ok(ws.tree(pg))
}

/// Predecessor lists from one search from `sources` (see the module
/// documentation), with the distances.
pub fn predecessors(
    pg: &PathGraph,
    sources: &[u32],
    opts: SearchOptions,
    tie: Tie,
) -> Result<(NodeLists, Vec<f64>)> {
    let mut ws = Searcher::new(pg.node_count());
    ws.run(pg, sources, opts)?;
    let dist = (0..pg.node_count()).map(|v| ws.dist(v)).collect();
    Ok((ws.predecessors(pg, tie, opts.reverse), dist))
}

/// Distances (row-major, one row per source) and parents ([`NO_NODE`] for
/// sources and unreached nodes) from every source in `sources` (every node
/// when `None`), searched in parallel; `cutoff` bounds the distances. Rows
/// are independent, so the result does not depend on the thread count.
pub fn all_pairs(
    pg: &PathGraph,
    sources: Option<&[u32]>,
    cutoff: f64,
    threads: usize,
) -> Result<(Vec<f64>, Vec<u32>)> {
    let n = pg.node_count();
    let all: Vec<u32>;
    let sources = match sources {
        Some(s) => s,
        None => {
            all = (0..n as u32).collect();
            &all
        }
    };
    let opts = SearchOptions {
        cutoff,
        ..Default::default()
    };
    Searcher::check(pg, sources, &opts)?;
    if !pg.is_unit() {
        pg.require_nonnegative("Dijkstra's algorithm")?;
    }
    let mut dist = vec![f64::INFINITY; sources.len() * n];
    let mut parent = vec![NO_NODE; sources.len() * n];
    let threads = parallel::threads_for(
        threads,
        sources.len(),
        parallel::search_work(sources.len(), n, pg.arc_count()),
    );
    fill_rows(
        &mut dist,
        &mut parent,
        n,
        threads,
        |ws: &mut Searcher, k, drow, prow| {
            // Errors were ruled out above.
            let _ = ws.run(pg, &sources[k..k + 1], opts);
            for &v in ws.order() {
                let v = v as usize;
                drow[v] = ws.dist[v];
                prow[v] = ws.parent[v];
            }
        },
    );
    Ok((dist, parent))
}

/// Calls `f(state, k, dist_row, parent_row)` for every row `k` of two
/// row-major matrices with `n` columns, in parallel, each thread with its
/// own [`Searcher`].
pub(crate) fn fill_rows<F>(dist: &mut [f64], parent: &mut [u32], n: usize, threads: usize, f: F)
where
    F: Fn(&mut Searcher, usize, &mut [f64], &mut [u32]) + Sync,
{
    if n == 0 {
        return;
    }
    let rows = dist.len() / n;
    let threads = parallel::resolve_threads(threads, rows);
    if threads <= 1 {
        let mut ws = Searcher::new(n);
        for (k, (d, p)) in dist.chunks_mut(n).zip(parent.chunks_mut(n)).enumerate() {
            f(&mut ws, k, d, p);
        }
        return;
    }
    let queue = std::sync::Mutex::new(dist.chunks_mut(n).zip(parent.chunks_mut(n)).enumerate());
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let (queue, f) = (&queue, &f);
            scope.spawn(move || {
                let mut ws = Searcher::new(n);
                loop {
                    let next = queue.lock().unwrap_or_else(|e| e.into_inner()).next();
                    match next {
                        Some((k, (d, p))) => f(&mut ws, k, d, p),
                        None => break,
                    }
                }
            });
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use crate::metrics::degree::DegreeMode;
    use crate::paths::{ArcLengths, PathGraphOptions};
    use crate::weights::WeightKind;

    fn pg(n: usize, s: &[i64], t: &[i64], w: Option<&[f64]>, directed: bool) -> PathGraph {
        let g = Graph::from_edges(n, s, t, w, directed).unwrap();
        let lengths = if w.is_some() {
            ArcLengths::Weights(WeightKind::Distance)
        } else {
            ArcLengths::Unit
        };
        PathGraph::new(
            &g,
            PathGraphOptions {
                lengths,
                ..Default::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn bfs_and_dijkstra_parents() {
        // Diamond 0-1-3, 0-2-3 plus a tail 3-4.
        let p = pg(5, &[0, 0, 1, 2, 3], &[1, 2, 3, 3, 4], None, false);
        let t = single_source(&p, 0, SearchOptions::default()).unwrap();
        assert_eq!(t.dist, vec![0.0, 1.0, 1.0, 2.0, 3.0]);
        assert_eq!(t.path_to(4).unwrap(), vec![0, 1, 3, 4]);
        let w = [1.0, 1.0, 2.0, 1.0, 0.5];
        let p = pg(5, &[0, 0, 1, 2, 3], &[1, 2, 3, 3, 4], Some(&w), false);
        let t = single_source(&p, 0, SearchOptions::default()).unwrap();
        assert_eq!(t.dist, vec![0.0, 1.0, 1.0, 2.0, 2.5]);
        assert_eq!(t.path_to(4).unwrap(), vec![0, 2, 3, 4]);
        assert_eq!(t.parent_edge[3], 3);
        let (pred, _) = predecessors(&p, &[0], SearchOptions::default(), Tie::Exact).unwrap();
        assert_eq!(pred.get(3), &[2]);
        let w = [1.0, 1.0, 1.0, 1.0, 1.0];
        let p = pg(5, &[0, 0, 1, 2, 3], &[1, 2, 3, 3, 4], Some(&w), false);
        let (pred, _) = predecessors(&p, &[0], SearchOptions::default(), Tie::Exact).unwrap();
        assert_eq!(pred.get(3), &[1, 2]);
        assert!(pred.get(0).is_empty());
    }

    #[test]
    fn cutoff_target_and_reverse() {
        let w = [1.0, 2.0, 3.0];
        let p = pg(4, &[0, 1, 2], &[1, 2, 3], Some(&w), true);
        let opts = SearchOptions {
            cutoff: 3.0,
            ..Default::default()
        };
        let t = single_source(&p, 0, opts).unwrap();
        assert_eq!(t.dist[2], 3.0);
        assert!(t.dist[3].is_infinite());
        let opts = SearchOptions {
            reverse: true,
            ..Default::default()
        };
        let t = single_source(&p, 3, opts).unwrap();
        assert_eq!(t.dist, vec![6.0, 5.0, 3.0, 0.0]);
        assert_eq!(t.path_to(0).unwrap(), vec![3, 2, 1, 0]);
        let opts = SearchOptions {
            target: Some(1),
            ..Default::default()
        };
        let t = single_source(&p, 0, opts).unwrap();
        assert_eq!(t.dist[1], 1.0);
        assert!(t.dist[3].is_infinite());
    }

    #[test]
    fn multi_source_records_the_nearest_source() {
        let p = pg(6, &[0, 1, 2, 3, 4], &[1, 2, 3, 4, 5], None, false);
        let t = multi_source(&p, &[0, 5], SearchOptions::default()).unwrap();
        assert_eq!(t.dist, vec![0.0, 1.0, 2.0, 2.0, 1.0, 0.0]);
        assert_eq!(t.origin, vec![0, 0, 0, 5, 5, 5]);
    }

    #[test]
    fn all_pairs_matches_single_source_for_any_thread_count() {
        let w = [1.0, 2.5, 0.5, 4.0, 1.0, 0.0];
        let p = pg(5, &[0, 1, 2, 3, 0, 4], &[1, 2, 3, 4, 3, 2], Some(&w), true);
        let (d1, p1) = all_pairs(&p, None, f64::INFINITY, 1).unwrap();
        let (d4, p4) = all_pairs(&p, None, f64::INFINITY, 4).unwrap();
        assert_eq!(d1, d4);
        assert_eq!(p1, p4);
        for s in 0..5 {
            let t = single_source(&p, s as u32, SearchOptions::default()).unwrap();
            assert_eq!(&d1[s * 5..s * 5 + 5], &t.dist[..]);
        }
    }

    #[test]
    fn rejects_negative_lengths() {
        let p = pg(2, &[0], &[1], Some(&[-1.0]), true);
        assert!(single_source(&p, 0, SearchOptions::default()).is_err());
        let g = Graph::from_edges(2, &[0], &[1], None, true).unwrap();
        let p = PathGraph::new(
            &g,
            PathGraphOptions {
                mode: DegreeMode::All,
                ..Default::default()
            },
        )
        .unwrap();
        let t = single_source(&p, 1, SearchOptions::default()).unwrap();
        assert_eq!(t.dist, vec![1.0, 0.0]);
    }
}
