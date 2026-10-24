//! Degree-preserving rewiring and the measures built on it: double edge
//! swaps, connected double edge swaps, directed three-edge swaps, random
//! and lattice reference graphs, the small-world coefficients σ and ω, and
//! the rich-club coefficient.
//!
//! # Swaps
//!
//! A double edge swap replaces edges u–v and x–y by u–x and v–y when
//! neither new edge exists, which keeps every degree (Maslov and Sneppen
//! 2002, `maslov2002specificity`). [`SwapGraph`] holds a simple graph with O(1) edge
//! tests, insertions, and removals, uniform neighbour draws, and draws of
//! a node with probability proportional to its degree (a uniform edge
//! end); it is the primitive every function here uses, and it can serve
//! incremental updates.
//!
//! The procedures are those of NetworkX 3.7 (`hagberg2008exploring`),
//! step for step, with draws from the counter-based stream (`seed`):
//!
//! - [`double_edge_swap`] (`double_edge_swap`): draw u and x with
//!   probability proportional to degree (redrawing equal nodes), v and y
//!   uniform neighbours (redrawing when v = y), and swap when x is not a
//!   neighbour of u and y not a neighbour of v; an attempt that passes the
//!   redraws counts toward `max_tries`, and exceeding it is an error.
//! - [`connected_double_edge_swap`] (`connected_double_edge_swap`): the
//!   same swaps in windows (the procedure NetworkX attributes to
//!   Gkantsidis, Mihail, and Zegura 2003): while the window is below 3 each swap is
//!   checked with a path search between u and v and undone if they are
//!   disconnected; larger windows are checked once for connectivity and
//!   undone together; the window grows by one after a success and halves
//!   after a failure. `nswap` attempts are made; the number of swaps kept
//!   is returned.
//! - [`directed_edge_swap`] (`directed_edge_swap`): a three-edge swap of a
//!   directed path a → b → c → d into a → c → b → d (Erdős, Miklós, and
//!   Toroczkai 2010, `erdos2010simple`), keeping in- and out-degrees; a is drawn with
//!   probability proportional to total degree and each next node as a
//!   uniform successor; every draw of a counts toward `max_tries`.
//! - [`random_reference`] and [`lattice_reference`] (`random_reference`,
//!   `lattice_reference`): `niter` times the number of edges rounds, each
//!   of at most ⌊2m / (n − 1)⌋ counted attempts at swapping a–b, c–d into
//!   a–d, c–b; draws with b or d among the other three nodes are redrawn;
//!   with `connectivity` a swap that disconnects a from b is undone. The
//!   lattice version (Sporns and Zwi 2004, `sporns2004small`) keeps a
//!   swap only when D[a, b] + D[c, d] ≥ D[a, c] + D[b, d] for the
//!   distance matrix D (by default the ring distance between node
//!   positions), which is the test NetworkX applies.
//!
//! NetworkX's loops redraw forever when no draw can pass (a star, for
//! example); here a run of 100 (m + 1) redraws in a row ends the round
//! (references) or fails with an error (swaps) instead.
//!
//! # Small-world coefficients and the rich club
//!
//! σ = (C / C_r) / (L / L_r) (Humphries and Gurney 2008,
//! `humphries2008network`), with C the transitivity, L the average
//! shortest path length, and C_r, L_r their means over `nrand` random
//! references; ω = L_r / L − C / C_l (Telesford et al. 2011,
//! `telesford2011ubiquity`), with C the average clustering, L_r from
//! random references with 2 `niter` rounds, and C_l the largest average
//! clustering of the graph and its `nrand` lattice references. As in
//! NetworkX these need a connected graph. The references are drawn from
//! streams (seed, r), so the results do not depend on the thread count.
//!
//! The rich-club coefficient φ(k) = 2 E_k / (N_k (N_k − 1)) for every
//! degree k with N_k > 1 nodes of degree above k, E_k the edges among
//! them (Zhou and Mondragón 2004, `zhou2004rich`), divided, when
//! normalized, by its mean over `n_samples` graphs randomized with Q m
//! double edge swaps (Colizza et al. 2006, `colizza2006detecting`), as
//! NetworkX's `rich_club_coefficient`.
//!
//! Conventions: undirected simple graphs (parallel edges merged; NetworkX
//! rejects multigraphs) except for the directed swap; self-loops are kept
//! by the swaps (a looped node is its own neighbour and its loop counts
//! twice toward its degree, as in NetworkX) and are an error for the rich
//! club, as in NetworkX. Weights are ignored. Results list edges (u, v)
//! with u ≤ v (directed: as arcs), sorted.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::matrix::Dense;
use crate::metrics::descriptive::{
    average_clustering, mean_distance, transitivity, ClusteringOptions, PathOptions, Undefined,
};
use crate::parallel;
use crate::rng::Stream;

/// A simple graph under rewiring. See the module documentation.
#[derive(Debug, Clone)]
pub struct SwapGraph {
    directed: bool,
    /// Out-neighbours (neighbours when undirected; a loop once).
    succ: Vec<Vec<u32>>,
    /// In-neighbours (directed only).
    pred: Vec<Vec<u32>>,
    /// Position of v in succ[u] and of u in pred[v], keyed by arc (u, v).
    at: HashMap<(u32, u32), (u32, u32)>,
    /// Edge ends at the start: a node drawn from here is drawn with
    /// probability proportional to its (total) degree.
    ends: Vec<u32>,
    edges: usize,
}

impl SwapGraph {
    /// The simple graph of `g` (parallel edges merged, loops kept).
    pub fn new(g: &Graph) -> SwapGraph {
        let n = g.node_count();
        let mut s = SwapGraph {
            directed: g.is_directed(),
            succ: vec![Vec::new(); n],
            pred: vec![Vec::new(); n],
            at: HashMap::new(),
            ends: Vec::new(),
            edges: 0,
        };
        for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
            if !s.has_edge(a, b) {
                s.add_edge(a, b);
            }
        }
        for v in 0..n {
            let d = s.degree(v);
            s.ends.extend(std::iter::repeat(v as u32).take(d));
        }
        s
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.succ.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges
    }

    /// Degree of v as NetworkX counts it (in + out for directed graphs; an
    /// undirected loop twice).
    pub fn degree(&self, v: usize) -> usize {
        if self.directed {
            self.succ[v].len() + self.pred[v].len()
        } else {
            self.succ[v].len() + self.succ[v].contains(&(v as u32)) as usize
        }
    }

    /// Out-neighbours (neighbours) of v, in no particular order.
    pub fn neighbors(&self, v: usize) -> &[u32] {
        &self.succ[v]
    }

    /// Whether the edge (arc) u–v exists.
    #[inline]
    pub fn has_edge(&self, u: u32, v: u32) -> bool {
        self.at.contains_key(&(u, v))
    }

    fn push(&mut self, u: u32, v: u32) {
        let i = self.succ[u as usize].len() as u32;
        self.succ[u as usize].push(v);
        let j = if self.directed {
            let j = self.pred[v as usize].len() as u32;
            self.pred[v as usize].push(u);
            j
        } else {
            0
        };
        self.at.insert((u, v), (i, j));
    }

    fn pull(&mut self, u: u32, v: u32) {
        let (i, j) = self.at.remove(&(u, v)).expect("the arc exists");
        let row = &mut self.succ[u as usize];
        row.swap_remove(i as usize);
        if (i as usize) < row.len() {
            let moved = row[i as usize];
            self.at.get_mut(&(u, moved)).expect("moved arc").0 = i;
        }
        if self.directed {
            let col = &mut self.pred[v as usize];
            col.swap_remove(j as usize);
            if (j as usize) < col.len() {
                let moved = col[j as usize];
                self.at.get_mut(&(moved, v)).expect("moved arc").1 = j;
            }
        }
    }

    /// Adds the edge (arc) u–v, which must be absent.
    pub fn add_edge(&mut self, u: u32, v: u32) {
        self.push(u, v);
        if !self.directed && u != v {
            self.push(v, u);
        }
        self.edges += 1;
    }

    /// Removes the edge (arc) u–v, which must exist.
    pub fn remove_edge(&mut self, u: u32, v: u32) {
        self.pull(u, v);
        if !self.directed && u != v {
            self.pull(v, u);
        }
        self.edges -= 1;
    }

    /// A node drawn with probability proportional to its initial degree.
    #[inline]
    pub fn by_degree(&self, rng: &mut Stream) -> u32 {
        self.ends[rng.below(self.ends.len() as u64) as usize]
    }

    /// A uniform out-neighbour (neighbour) of v; v must have one.
    #[inline]
    pub fn random_neighbor(&self, v: u32, rng: &mut Stream) -> u32 {
        let row = &self.succ[v as usize];
        row[rng.below(row.len() as u64) as usize]
    }

    /// Whether v is reachable from u (direction ignored).
    pub fn has_path(&self, u: u32, v: u32, seen: &mut Vec<bool>, queue: &mut Vec<u32>) -> bool {
        if u == v {
            return true;
        }
        seen.clear();
        seen.resize(self.node_count(), false);
        queue.clear();
        seen[u as usize] = true;
        queue.push(u);
        let mut head = 0;
        while head < queue.len() {
            let x = queue[head] as usize;
            head += 1;
            for &w in self.succ[x].iter().chain(&self.pred[x]) {
                if !seen[w as usize] {
                    if w == v {
                        return true;
                    }
                    seen[w as usize] = true;
                    queue.push(w);
                }
            }
        }
        false
    }

    /// Whether the graph is connected (direction ignored).
    pub fn is_connected(&self) -> bool {
        let n = self.node_count();
        if n == 0 {
            return false;
        }
        let mut seen = vec![false; n];
        let mut stack = vec![0u32];
        seen[0] = true;
        let mut count = 1;
        while let Some(x) = stack.pop() {
            for &w in self.succ[x as usize].iter().chain(&self.pred[x as usize]) {
                if !seen[w as usize] {
                    seen[w as usize] = true;
                    count += 1;
                    stack.push(w);
                }
            }
        }
        count == n
    }

    /// The edges, (u, v) with u ≤ v unless directed, sorted.
    pub fn edge_list(&self) -> Vec<(u32, u32)> {
        let mut out: Vec<(u32, u32)> = self
            .at
            .keys()
            .copied()
            .filter(|&(u, v)| self.directed || u <= v)
            .collect();
        out.sort_unstable();
        out
    }

    /// The graph as an engine graph.
    pub fn to_graph(&self) -> Graph {
        let e = self.edge_list();
        let s: Vec<i64> = e.iter().map(|p| p.0 as i64).collect();
        let t: Vec<i64> = e.iter().map(|p| p.1 as i64).collect();
        Graph::from_edges(self.node_count(), &s, &t, None, self.directed).expect("valid edges")
    }

    fn redraw_limit(&self) -> u64 {
        100 * (self.edges as u64 + 1)
    }
}

fn undirected_only(g: &Graph, what: &str) -> Result<()> {
    if g.is_directed() {
        return Err(Error::InvalidArgument(format!(
            "{what} is defined for undirected graphs"
        )));
    }
    Ok(())
}

fn too_many_redraws() -> Error {
    Error::InvalidArgument(
        "no swap can be drawn: every draw was rejected 100 (m + 1) times in a row".into(),
    )
}

/// `nswap` double edge swaps; returns the new edges. See the module
/// documentation.
pub fn double_edge_swap(
    g: &Graph,
    nswap: u64,
    max_tries: u64,
    seed: u64,
) -> Result<Vec<(u32, u32)>> {
    undirected_only(g, "the double edge swap")?;
    let mut s = SwapGraph::new(g);
    swap_in_place(&mut s, nswap, max_tries, &mut Stream::new(seed, 0))?;
    Ok(s.edge_list())
}

fn swap_in_place(s: &mut SwapGraph, nswap: u64, max_tries: u64, rng: &mut Stream) -> Result<()> {
    if nswap > max_tries {
        return Err(Error::InvalidArgument(
            "the number of swaps exceeds the number of tries allowed".into(),
        ));
    }
    if s.node_count() < 4 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than four nodes".into(),
        ));
    }
    if s.edge_count() < 2 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than 2 edges".into(),
        ));
    }
    let (mut n, mut swaps, mut redraws) = (0u64, 0u64, 0u64);
    while swaps < nswap {
        let u = s.by_degree(rng);
        let x = s.by_degree(rng);
        let (v, y) = if u == x {
            (0, 0)
        } else {
            (s.random_neighbor(u, rng), s.random_neighbor(x, rng))
        };
        if u == x || v == y {
            redraws += 1;
            if redraws > s.redraw_limit() {
                return Err(too_many_redraws());
            }
            continue;
        }
        redraws = 0;
        if !s.has_edge(u, x) && !s.has_edge(v, y) {
            s.add_edge(u, x);
            s.add_edge(v, y);
            s.remove_edge(u, v);
            s.remove_edge(x, y);
            swaps += 1;
        }
        if n >= max_tries {
            return Err(Error::InvalidArgument(format!(
                "the maximum number of swap attempts ({n}) was exceeded before {nswap} swaps"
            )));
        }
        n += 1;
    }
    Ok(())
}

/// Connectivity-preserving double edge swaps: `nswap` attempts; returns
/// the new edges and the number of swaps kept. See the module
/// documentation.
pub fn connected_double_edge_swap(
    g: &Graph,
    nswap: u64,
    seed: u64,
) -> Result<(Vec<(u32, u32)>, u64)> {
    undirected_only(g, "the connected double edge swap")?;
    let mut s = SwapGraph::new(g);
    if !s.is_connected() {
        return Err(Error::InvalidArgument("the graph is not connected".into()));
    }
    if s.node_count() < 4 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than four nodes".into(),
        ));
    }
    let rng = &mut Stream::new(seed, 0);
    let (mut seen, mut queue) = (Vec::new(), Vec::new());
    let (mut n, mut swapcount, mut window) = (0u64, 0i64, 1u64);
    let mut redraws = 0u64;
    // One proposal: Some((u, v, x, y)) once the redraws pass.
    let mut draw = |s: &SwapGraph, rng: &mut Stream| -> Result<(u32, u32, u32, u32)> {
        loop {
            let u = s.by_degree(rng);
            let x = s.by_degree(rng);
            if u != x {
                let v = s.random_neighbor(u, rng);
                let y = s.random_neighbor(x, rng);
                if v != y {
                    redraws = 0;
                    return Ok((u, v, x, y));
                }
            }
            redraws += 1;
            if redraws > s.redraw_limit() {
                return Err(too_many_redraws());
            }
        }
    };
    while n < nswap {
        let mut wcount = 0;
        let mut swapped: Vec<(u32, u32, u32, u32)> = Vec::new();
        if window < 3 {
            let mut fail = false;
            while wcount < window && n < nswap {
                let (u, v, x, y) = draw(&s, rng)?;
                let did = !s.has_edge(u, x) && !s.has_edge(v, y);
                if did {
                    s.remove_edge(u, v);
                    s.remove_edge(x, y);
                    s.add_edge(u, x);
                    s.add_edge(v, y);
                    swapcount += 1;
                }
                n += 1;
                if s.has_path(u, v, &mut seen, &mut queue) {
                    wcount += 1;
                } else {
                    s.add_edge(u, v);
                    s.add_edge(x, y);
                    s.remove_edge(u, x);
                    s.remove_edge(v, y);
                    swapcount -= 1;
                    fail = true;
                }
            }
            window = if fail { (window + 1) / 2 } else { window + 1 };
        } else {
            while wcount < window && n < nswap {
                let (u, v, x, y) = draw(&s, rng)?;
                if !s.has_edge(u, x) && !s.has_edge(v, y) {
                    s.remove_edge(u, v);
                    s.remove_edge(x, y);
                    s.add_edge(u, x);
                    s.add_edge(v, y);
                    swapped.push((u, v, x, y));
                    swapcount += 1;
                }
                n += 1;
                wcount += 1;
            }
            if s.is_connected() {
                window += 1;
            } else {
                while let Some((u, v, x, y)) = swapped.pop() {
                    s.add_edge(u, v);
                    s.add_edge(x, y);
                    s.remove_edge(u, x);
                    s.remove_edge(v, y);
                    swapcount -= 1;
                }
                window = (window + 1) / 2;
            }
        }
    }
    Ok((s.edge_list(), swapcount.max(0) as u64))
}

/// `nswap` directed three-edge swaps; returns the new arcs. See the
/// module documentation.
pub fn directed_edge_swap(
    g: &Graph,
    nswap: u64,
    max_tries: u64,
    seed: u64,
) -> Result<Vec<(u32, u32)>> {
    if !g.is_directed() {
        return Err(Error::InvalidArgument(
            "the directed edge swap needs a directed graph".into(),
        ));
    }
    let mut s = SwapGraph::new(g);
    if nswap > max_tries {
        return Err(Error::InvalidArgument(
            "the number of swaps exceeds the number of tries allowed".into(),
        ));
    }
    if s.node_count() < 4 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than four nodes".into(),
        ));
    }
    if s.edge_count() < 3 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than 3 edges".into(),
        ));
    }
    let rng = &mut Stream::new(seed, 0);
    let (mut tries, mut swaps) = (0u64, 0u64);
    while swaps < nswap {
        let a = s.by_degree(rng);
        tries += 1;
        if tries > max_tries {
            return Err(Error::InvalidArgument(format!(
                "the maximum number of swap attempts ({tries}) was exceeded before {nswap} swaps"
            )));
        }
        if s.neighbors(a as usize).is_empty() {
            continue;
        }
        let b = s.random_neighbor(a, rng);
        if a == b || s.neighbors(b as usize).is_empty() {
            continue;
        }
        let c = s.random_neighbor(b, rng);
        if b == c || s.neighbors(c as usize).is_empty() {
            continue;
        }
        let d = s.random_neighbor(c, rng);
        if c == d {
            continue;
        }
        if !s.has_edge(a, c) && !s.has_edge(b, d) && !s.has_edge(c, b) {
            s.add_edge(a, c);
            s.add_edge(c, b);
            s.add_edge(b, d);
            s.remove_edge(a, b);
            s.remove_edge(b, c);
            s.remove_edge(c, d);
            swaps += 1;
        }
    }
    Ok(s.edge_list())
}

/// The ring distance between node positions, NetworkX's default lattice.
fn ring(n: usize, i: usize, j: usize) -> f64 {
    let d = i.abs_diff(j);
    d.min(n - d) as f64
}

/// Which reference graph [`reference`] builds.
#[derive(Clone, Copy)]
enum Reference<'a> {
    Random,
    /// Lattice, with a distance matrix or the ring distance.
    Lattice(Option<&'a Dense>),
}

fn reference(
    g: &Graph,
    niter: u64,
    connectivity: bool,
    kind: Reference,
    rng: &mut Stream,
) -> Result<SwapGraph> {
    undirected_only(g, "a reference graph")?;
    let mut s = SwapGraph::new(g);
    let (nn, ne) = (s.node_count(), s.edge_count());
    if nn < 4 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than four nodes".into(),
        ));
    }
    if ne < 2 {
        return Err(Error::InvalidArgument(
            "the graph has fewer than 2 edges".into(),
        ));
    }
    if let Reference::Lattice(Some(d)) = kind {
        if d.rows() != nn || d.cols() != nn {
            return Err(Error::InvalidArgument(format!(
                "the distance matrix must be {nn} by {nn}"
            )));
        }
    }
    // Whether the lattice test lets a-b, c-d become a-d, c-b.
    let towards_lattice = |a: u32, b: u32, c: u32, d: u32| -> bool {
        let dist = |i: u32, j: u32| match kind {
            Reference::Lattice(Some(m)) => m.get(i as usize, j as usize),
            _ => ring(nn, i as usize, j as usize),
        };
        match kind {
            Reference::Random => true,
            Reference::Lattice(_) => dist(a, b) + dist(c, d) >= dist(a, c) + dist(b, d),
        }
    };
    let ntries = (nn as f64 * ne as f64 / (nn as f64 * (nn as f64 - 1.0) / 2.0)) as u64;
    let (mut seen, mut queue) = (Vec::new(), Vec::new());
    for _ in 0..niter.saturating_mul(ne as u64) {
        let (mut n, mut redraws) = (0u64, 0u64);
        while n < ntries {
            let a = s.by_degree(rng);
            let c = s.by_degree(rng);
            let mut reject = a == c;
            let (mut b, mut d) = (0, 0);
            if !reject {
                b = s.random_neighbor(a, rng);
                d = s.random_neighbor(c, rng);
                reject = b == a || b == c || b == d || d == a || d == c;
            }
            if reject {
                redraws += 1;
                if redraws > s.redraw_limit() {
                    break;
                }
                continue;
            }
            redraws = 0;
            if !s.has_edge(a, d) && !s.has_edge(c, b) && towards_lattice(a, b, c, d) {
                s.add_edge(a, d);
                s.add_edge(c, b);
                s.remove_edge(a, b);
                s.remove_edge(c, d);
                if connectivity && !s.has_path(a, b, &mut seen, &mut queue) {
                    s.remove_edge(a, d);
                    s.remove_edge(c, b);
                    s.add_edge(a, b);
                    s.add_edge(c, d);
                } else {
                    break;
                }
            }
            n += 1;
        }
    }
    Ok(s)
}

/// A random reference graph (NetworkX's `random_reference`); returns its
/// edges. See the module documentation.
pub fn random_reference(
    g: &Graph,
    niter: u64,
    connectivity: bool,
    seed: u64,
) -> Result<Vec<(u32, u32)>> {
    let s = reference(
        g,
        niter,
        connectivity,
        Reference::Random,
        &mut Stream::new(seed, 0),
    )?;
    Ok(s.edge_list())
}

/// A lattice reference graph (NetworkX's `lattice_reference`), towards
/// the distance matrix `d` (the ring distance between node positions when
/// `None`); returns its edges. See the module documentation.
pub fn lattice_reference(
    g: &Graph,
    niter: u64,
    d: Option<&Dense>,
    connectivity: bool,
    seed: u64,
) -> Result<Vec<(u32, u32)>> {
    let s = reference(
        g,
        niter,
        connectivity,
        Reference::Lattice(d),
        &mut Stream::new(seed, 0),
    )?;
    Ok(s.edge_list())
}

/// NetworkX's clustering conventions: 0 / 0 counts as 0. Inside a
/// parallel loop over references the measures run on one thread.
fn clustering_nx() -> ClusteringOptions {
    ClusteringOptions {
        undefined: Undefined::Zero,
        threads: 1,
        ..Default::default()
    }
}

fn paths_nx() -> PathOptions {
    PathOptions {
        threads: 1,
        ..Default::default()
    }
}

fn connected_aspl(g: &Graph) -> Result<f64> {
    let s = SwapGraph::new(g);
    if !s.is_connected() {
        return Err(Error::InvalidArgument(
            "the small-world coefficients need a connected graph".into(),
        ));
    }
    mean_distance(g, paths_nx(), false)
}

/// The small-world coefficient σ. See the module documentation.
pub fn sigma(g: &Graph, niter: u64, nrand: usize, seed: u64, threads: usize) -> Result<f64> {
    undirected_only(g, "sigma")?;
    let l = connected_aspl(g)?;
    let c = transitivity(g, clustering_nx());
    let work = parallel::search_work(nrand * 8, g.node_count(), g.edge_count())
        .saturating_mul(niter.max(1));
    let threads = parallel::threads_for(threads, nrand, work);
    let refs: Vec<Result<(f64, f64)>> = parallel::map_items(
        nrand,
        threads,
        || (),
        |_, r| {
            let s = reference(
                g,
                niter,
                true,
                Reference::Random,
                &mut Stream::new(seed, r as u64),
            )?;
            let h = s.to_graph();
            Ok((
                transitivity(&h, clustering_nx()),
                mean_distance(&h, paths_nx(), false)?,
            ))
        },
    );
    let (mut cr, mut lr) = (0.0, 0.0);
    for x in refs {
        let (a, b) = x?;
        cr += a;
        lr += b;
    }
    let (cr, lr) = (cr / nrand as f64, lr / nrand as f64);
    Ok(c / cr / (l / lr))
}

/// The small-world coefficient ω. See the module documentation.
pub fn omega(g: &Graph, niter: u64, nrand: usize, seed: u64, threads: usize) -> Result<f64> {
    undirected_only(g, "omega")?;
    let l = connected_aspl(g)?;
    let c = average_clustering(g, clustering_nx())?;
    let work = parallel::search_work(nrand * 8, g.node_count(), g.edge_count())
        .saturating_mul(niter.max(1));
    let threads = parallel::threads_for(threads, nrand, work);
    let refs: Vec<Result<(f64, f64)>> = parallel::map_items(
        nrand,
        threads,
        || (),
        |_, r| {
            // The random and lattice references of round r use streams 2r
            // and 2r + 1.
            let rnd = reference(
                g,
                niter * 2,
                true,
                Reference::Random,
                &mut Stream::new(seed, 2 * r as u64),
            )?
            .to_graph();
            let lat = reference(
                g,
                niter,
                true,
                Reference::Lattice(None),
                &mut Stream::new(seed, 2 * r as u64 + 1),
            )?
            .to_graph();
            Ok((
                mean_distance(&rnd, paths_nx(), false)?,
                average_clustering(&lat, clustering_nx())?,
            ))
        },
    );
    let (mut lr, mut cl) = (0.0, c);
    for x in refs {
        let (a, b) = x?;
        lr += a;
        if b > cl {
            cl = b;
        }
    }
    let lr = lr / nrand as f64;
    Ok(lr / l - c / cl)
}

/// φ(k) for k = 0, 1, ... while more than one node has degree above k.
fn rich_club_raw(s: &SwapGraph) -> Vec<f64> {
    let n = s.node_count();
    let deg: Vec<usize> = (0..n).map(|v| s.degree(v)).collect();
    let max = deg.iter().copied().max().unwrap_or(0);
    let mut above = vec![0usize; max + 2];
    for &d in &deg {
        above[d] += 1;
    }
    // above[k] = nodes with degree > k.
    let mut acc = 0;
    for k in (0..=max + 1).rev() {
        let here = above[k];
        above[k] = acc;
        acc += here;
    }
    let mut edges_above = vec![0usize; max + 2];
    for (u, v) in s.edge_list() {
        let m = deg[u as usize].min(deg[v as usize]);
        // The edge joins nodes of degree > k for every k < m.
        edges_above[m] += 1;
    }
    let mut acc = 0;
    for k in (0..=max + 1).rev() {
        let here = edges_above[k];
        edges_above[k] = acc;
        acc += here;
    }
    if s.edge_count() == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for k in 0..=max {
        let nk = above[k];
        if nk <= 1 {
            break;
        }
        out.push(2.0 * edges_above[k] as f64 / (nk as f64 * (nk as f64 - 1.0)));
    }
    out
}

/// The rich-club coefficient φ(k), k = 0, 1, ..., normalized by its mean
/// over `n_samples` graphs randomized with `q` m double edge swaps (stream
/// (seed, sample)) when `normalized`. See the module documentation.
pub fn rich_club_coefficient(
    g: &Graph,
    normalized: bool,
    q: u64,
    n_samples: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<f64>> {
    undirected_only(g, "the rich-club coefficient")?;
    if g.edge_sources()
        .iter()
        .zip(g.edge_targets())
        .any(|(a, b)| a == b)
    {
        return Err(Error::InvalidArgument(
            "the rich-club coefficient is not defined for graphs with self-loops".into(),
        ));
    }
    let s = SwapGraph::new(g);
    let rc = rich_club_raw(&s);
    if !normalized {
        return Ok(rc);
    }
    if n_samples == 0 {
        return Err(Error::InvalidArgument(
            "n_samples must be at least 1".into(),
        ));
    }
    let e = s.edge_count() as u64;
    let work = (q.saturating_mul(e)).saturating_mul(64);
    let threads = parallel::threads_for(threads, n_samples, work);
    let samples: Vec<Result<Vec<f64>>> = parallel::map_items(
        n_samples,
        threads,
        || (),
        |_, r| {
            let mut x = s.clone();
            swap_in_place(&mut x, q * e, q * e * 10, &mut Stream::new(seed, r as u64))?;
            Ok(rich_club_raw(&x))
        },
    );
    let mut avg = vec![0.0; rc.len()];
    for smp in samples {
        let smp = smp?;
        for (k, a) in avg.iter_mut().enumerate() {
            *a += smp.get(k).copied().unwrap_or(0.0) / n_samples as f64;
        }
    }
    if let Some(k) = avg.iter().position(|&a| a == 0.0) {
        return Err(Error::InvalidArgument(format!(
            "rich-club normalization failed at degree {k}: the mean coefficient of the random \
             graphs is zero; try more samples"
        )));
    }
    Ok(rc.iter().zip(&avg).map(|(r, a)| r / a).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn degrees(n: usize, e: &[(u32, u32)]) -> Vec<usize> {
        let mut d = vec![0; n];
        for &(a, b) in e {
            d[a as usize] += 1;
            d[b as usize] += 1;
        }
        d
    }

    #[test]
    fn swaps_keep_degrees() {
        let g = Graph::from_edges(
            8,
            &[0, 1, 2, 3, 4, 5, 6, 7, 0, 2],
            &[1, 2, 3, 4, 5, 6, 7, 0, 4, 6],
            None,
            false,
        )
        .unwrap();
        let before = degrees(8, &SwapGraph::new(&g).edge_list());
        let e = double_edge_swap(&g, 5, 1000, 1).unwrap();
        assert_eq!(degrees(8, &e), before);
        let (e, k) = connected_double_edge_swap(&g, 20, 2).unwrap();
        assert_eq!(degrees(8, &e), before);
        assert!(k <= 20);
        let e = random_reference(&g, 2, true, 3).unwrap();
        assert_eq!(degrees(8, &e), before);
        let e = lattice_reference(&g, 2, None, true, 4).unwrap();
        assert_eq!(degrees(8, &e), before);
        // A complete graph admits no swap: sigma is 1 and omega 0.
        let mut s = Vec::new();
        let mut t = Vec::new();
        for a in 0..6i64 {
            for b in a + 1..6 {
                s.push(a);
                t.push(b);
            }
        }
        let k6 = Graph::from_edges(6, &s, &t, None, false).unwrap();
        assert!((sigma(&k6, 2, 3, 1, 1).unwrap() - 1.0).abs() < 1e-12);
        assert!(omega(&k6, 2, 3, 1, 1).unwrap().abs() < 1e-12);
        // Rich club of a star with a triangle: degrees 4, 3, 2, 1, 1, 1.
        let h =
            Graph::from_edges(6, &[0, 0, 0, 0, 1, 1], &[1, 2, 3, 4, 2, 5], None, false).unwrap();
        let rc = rich_club_coefficient(&h, false, 1, 1, 0, 1).unwrap();
        // k = 0: 6 nodes, 6 edges -> 12 / 30; k = 1: nodes 0, 1, 2 with
        // 3 edges -> 1; k = 2: nodes 0, 1 with 1 edge -> 1.
        assert_eq!(rc, vec![0.4, 1.0, 1.0]);
    }
}
