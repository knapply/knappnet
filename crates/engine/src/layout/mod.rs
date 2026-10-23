//! Node coordinates for drawing a network in the plane.
//!
//! | Function | Method | Cost | Source |
//! |---|---|---|---|
//! | [`stress_layout`] | stress majorization on all pairs (the default up to [`FULL_STRESS_AUTO_MAX`] nodes per component) | O(n^2) memory and time per iteration | Gansner, Koren, and North 2004 (`gansner2004graph`) |
//! | [`stress_layout`] | sparse stress with pivots (the default for larger components) | O(k n + m) per iteration | Ortmann, Klimenta, and Brandes 2017 (`ortmann2017sparse`) |
//! | [`fr_layout`] | force-directed placement, exact or with a Barnes-Hut quadtree | O(n^2) or O(n log n) per iteration | Fruchterman and Reingold 1991 (`fruchterman1991graph`), Barnes and Hut 1986 (`barnes1986hierarchical`) |
//! | [`circle_layout`], [`star_layout`], [`grid_layout`] | fixed positions | O(n) | |
//! | [`tree_layout`] | layered tree over a breadth-first spanning tree | O(n + m) | Reingold and Tilford 1981 (`reingold1981tidier`), aesthetics only |
//! | [`twomode_layout`] | two rows, ordered by the barycenter heuristic | O(m log m) per sweep | Sugiyama, Tagawa, and Toda 1981 (`sugiyama1981methods`) |
//!
//! Stress layouts start from pivot multidimensional scaling (Brandes and
//! Pich 2007, `brandes2007eigensolver`), the initialization recommended by
//! Brandes and Pich (2009, `brandes2009experimental`). The stress of a
//! layout is
//!
//! ```text
//! stress(X) = sum over pairs i < j in the same component of
//!             w_ij (||x_i - x_j|| - d_ij)^2,   w_ij = d_ij^-2,
//! ```
//!
//! where d_ij is the shortest-path length; [`layout_stress`] evaluates it.
//! The weights d^-2 are those of Kamada and Kawai (1989,
//! `kamada1989algorithm`), so this is the energy that their spring model
//! and NetworkX's `kamada_kawai_layout()` minimize.
//!
//! Conventions, shared by every method:
//! - Direction is ignored: a directed network is drawn from its
//!   undirected skeleton.
//! - Self-loops are ignored, and parallel edges count once. With weights,
//!   the shortest of several parallel edges gives the target distance
//!   (stress), and their strengths add up (force-directed attraction).
//! - Weights follow [`crate::weights`]: distances, or strengths converted to
//!   lengths 1/w, and they must be strictly positive. In stress layouts the
//!   target distance between two nodes is their weighted shortest-path
//!   length; in force-directed layouts a tie's attraction is multiplied by
//!   its strength, 1 / length. Without weights every edge has length 1.
//! - Disconnected networks: stress, force-directed, and tree layouts treat
//!   each weak component separately and then pack the components' bounding
//!   boxes in rows, largest first, with a gap of one typical edge length
//!   (the shelf algorithm "next-fit decreasing height" of Coffman, Garey,
//!   Johnson, and Tarjan 1980, `coffman1980performance`). Isolated nodes
//!   are components of one node and end up in rows of their own. Circle,
//!   star, grid, and two-mode layouts place all nodes together.
//! - Coordinates are in units of path length (hops without weights) for
//!   stress layouts and of the ideal edge length 1 for force-directed ones.
//!   Each stress or force-directed component is centered and rotated so its
//!   principal axis is horizontal, and reflected so that its
//!   lowest-numbered node lies left of center; layouts are otherwise only
//!   defined up to rotation and reflection.
//! - Randomness comes from [`crate::rng::Stream`] keyed by the seed and the
//!   component's lowest-numbered node. Every parallel loop computes each
//!   node's value on one thread and every sum over nodes runs in node
//!   order, so a layout is bit-for-bit identical for a given seed on any
//!   number of threads.
//! - Approximations (pivots in sparse stress, the Barnes-Hut quadtree) only
//!   change how positions are optimized. No measure reported elsewhere in
//!   the library depends on a layout.

mod force;
mod simple;
mod stress;

pub use force::{fr_layout, FrApproximation, FrOptions, FR_EXACT_AUTO_MAX};
pub use simple::{circle_layout, grid_layout, star_layout, tree_layout, twomode_layout};
pub use stress::{layout_stress, stress_layout, StressAlgorithm, StressOptions, FULL_STRESS_MAX};

/// Largest component that [`StressAlgorithm::Auto`] lays out with stress
/// majorization on all pairs; larger components use sparse stress.
pub const FULL_STRESS_AUTO_MAX: usize = 2000;

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::Result;
use crate::graph::Graph;
use crate::metrics::topography::{components, ComponentKind};
use crate::parallel;
use crate::weights::{path_lengths, Lengths, PathAlgorithm, WeightKind};

/// Node coordinates, aligned with the graph's node indices.
#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    /// Horizontal coordinates.
    pub x: Vec<f64>,
    /// Vertical coordinates.
    pub y: Vec<f64>,
}

impl Layout {
    fn zeros(n: usize) -> Layout {
        Layout {
            x: vec![0.0; n],
            y: vec![0.0; n],
        }
    }
}

/// How shortest-path lengths are measured inside a [`Part`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Metric {
    /// Every arc has the same length; breadth-first search.
    Unit,
    /// Positive integer lengths up to the given maximum; Dial's algorithm.
    Integer(u32),
    /// Positive real lengths; Dijkstra's algorithm.
    Real,
}

/// One weak component of the undirected skeleton, with local node indices
/// `0..nodes.len()` in increasing order of global index. Self-loops are
/// dropped and parallel arcs merged.
#[derive(Debug, Clone)]
pub(crate) struct Part {
    /// Global node index of each local node, ascending.
    pub nodes: Vec<u32>,
    pub offsets: Vec<usize>,
    /// Local neighbor indices, ascending within each row.
    pub targets: Vec<u32>,
    /// Shortest length among the merged parallel arcs.
    pub lens: Vec<f64>,
    /// Sum of 1 / length over the merged parallel arcs (1 per tie without
    /// weights).
    pub strength: Vec<f64>,
    pub metric: Metric,
}

impl Part {
    #[inline]
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    #[inline]
    pub fn arcs(&self, v: usize) -> std::ops::Range<usize> {
        self.offsets[v]..self.offsets[v + 1]
    }

    #[inline]
    pub fn degree(&self, v: usize) -> usize {
        self.offsets[v + 1] - self.offsets[v]
    }
}

/// Splits `g`, read as undirected, into weak components with arc lengths.
/// Also returns the mean merged arc length (the uniform length, normally 1,
/// if there are no arcs), the typical edge length used for gaps between
/// packed components.
pub(crate) fn parts(g: &Graph, weighted: bool, kind: WeightKind) -> Result<(Vec<Part>, f64)> {
    let n = g.node_count();
    let owned;
    let und = if g.is_directed() {
        let src: Vec<i64> = g.edge_sources().iter().map(|&v| v as i64).collect();
        let tgt: Vec<i64> = g.edge_targets().iter().map(|&v| v as i64).collect();
        owned = Graph::from_edges(n, &src, &tgt, g.edge_weights(), false)?;
        &owned
    } else {
        g
    };
    let lengths = path_lengths(und, weighted, kind, PathAlgorithm::Auto)?;
    let adj = und.out();
    // Uniform lengths keep their common value, so coordinates stay in
    // units of the weights.
    let unit = match (&lengths, adj.all_weights()) {
        (Lengths::Uniform, Some(w)) if weighted && !w.is_empty() => match kind {
            WeightKind::Distance => w[0],
            WeightKind::Strength => 1.0 / w[0],
        },
        _ => 1.0,
    };
    let metric = match &lengths {
        Lengths::Uniform => Metric::Unit,
        Lengths::Integer(_, max) => Metric::Integer(*max),
        Lengths::Real(_) => Metric::Real,
    };
    let arc_len = |a: usize| -> f64 {
        match &lengths {
            Lengths::Uniform => unit,
            Lengths::Integer(v, _) => v[a] as f64,
            Lengths::Real(v) => v[a],
        }
    };

    let comps = components(und, ComponentKind::Weak);
    let mut members: Vec<Vec<u32>> = comps.sizes.iter().map(|&s| Vec::with_capacity(s)).collect();
    let mut local = vec![0u32; n];
    for (v, (l, &c)) in local.iter_mut().zip(&comps.membership).enumerate() {
        let c = c as usize;
        *l = members[c].len() as u32;
        members[c].push(v as u32);
    }

    let mut total_len = 0.0;
    let mut total_arcs = 0usize;
    let mut out = Vec::with_capacity(members.len());
    for nodes in members {
        let mut offsets = Vec::with_capacity(nodes.len() + 1);
        let mut targets = Vec::new();
        let mut lens = Vec::new();
        let mut strength = Vec::new();
        offsets.push(0);
        for &v in &nodes {
            let v = v as usize;
            let range = adj.arc_range(v);
            for (a, &t) in range.clone().zip(adj.neighbors(v)) {
                if t as usize == v {
                    continue;
                }
                let lt = local[t as usize];
                let len = arc_len(a);
                // Rows are sorted by target, so parallel arcs are adjacent.
                if targets.len() > *offsets.last().unwrap() && *targets.last().unwrap() == lt {
                    let last = lens.len() - 1;
                    if len < lens[last] {
                        lens[last] = len;
                    }
                    if weighted {
                        strength[last] += 1.0 / len;
                    }
                } else {
                    targets.push(lt);
                    lens.push(len);
                    strength.push(if weighted { 1.0 / len } else { 1.0 });
                }
            }
            offsets.push(targets.len());
        }
        total_len += lens.iter().sum::<f64>();
        total_arcs += lens.len();
        out.push(Part {
            nodes,
            offsets,
            targets,
            lens,
            strength,
            metric,
        });
    }
    let scale = if total_arcs > 0 {
        total_len / total_arcs as f64
    } else {
        unit
    };
    Ok((out, scale))
}

/// Reusable buffers for single-source shortest paths in a [`Part`].
pub(crate) struct Sssp {
    queue: Vec<u32>,
    idist: Vec<u64>,
    buckets: Vec<Vec<u32>>,
    heap: BinaryHeap<HeapItem>,
}

#[derive(PartialEq)]
struct HeapItem(f64, u32);
impl Eq for HeapItem {}
impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Smallest distance first; ties by node index.
        other
            .0
            .total_cmp(&self.0)
            .then_with(|| other.1.cmp(&self.1))
    }
}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Sssp {
    pub fn new(n: usize) -> Sssp {
        Sssp {
            queue: Vec::with_capacity(n),
            idist: vec![u64::MAX; n],
            buckets: Vec::new(),
            heap: BinaryHeap::new(),
        }
    }

    /// Shortest-path lengths from `s` to every node of `p` into `dist`.
    pub fn run(&mut self, p: &Part, s: usize, dist: &mut [f64]) {
        for d in dist.iter_mut() {
            *d = f64::INFINITY;
        }
        match p.metric {
            Metric::Unit => self.bfs(p, s, dist),
            Metric::Integer(max) => self.dial(p, s, max, dist),
            Metric::Real => self.dijkstra(p, s, dist),
        }
    }

    fn bfs(&mut self, p: &Part, s: usize, dist: &mut [f64]) {
        // Every arc has the same length; hop counts are exact integers.
        let unit = p.lens.first().copied().unwrap_or(1.0);
        self.queue.clear();
        self.queue.push(s as u32);
        self.idist[s] = 0;
        let mut head = 0;
        while head < self.queue.len() {
            let v = self.queue[head] as usize;
            head += 1;
            let next = self.idist[v] + 1;
            for &t in &p.targets[p.arcs(v)] {
                let t = t as usize;
                if self.idist[t] == u64::MAX {
                    self.idist[t] = next;
                    self.queue.push(t as u32);
                }
            }
        }
        for &v in &self.queue {
            let v = v as usize;
            dist[v] = self.idist[v] as f64 * unit;
            self.idist[v] = u64::MAX;
        }
    }

    /// Dial's bucket queue (`dial1969algorithm`) with exact integer sums.
    fn dial(&mut self, p: &Part, s: usize, max: u32, dist: &mut [f64]) {
        let width = max as usize + 1;
        if self.buckets.len() != width {
            self.buckets = vec![Vec::new(); width];
        }
        self.queue.clear();
        self.idist[s] = 0;
        self.buckets[0].push(s as u32);
        let mut pending = 1usize;
        let mut d: u64 = 0;
        while pending > 0 {
            let slot = (d % width as u64) as usize;
            while let Some(v) = self.buckets[slot].pop() {
                pending -= 1;
                let v = v as usize;
                if self.idist[v] != d || dist[v] != f64::INFINITY {
                    continue;
                }
                dist[v] = d as f64;
                self.queue.push(v as u32);
                for a in p.arcs(v) {
                    let t = p.targets[a] as usize;
                    let nd = d + p.lens[a] as u64;
                    if nd < self.idist[t] {
                        self.idist[t] = nd;
                        self.buckets[(nd % width as u64) as usize].push(t as u32);
                        pending += 1;
                    }
                }
            }
            d += 1;
        }
        for &v in &self.queue {
            self.idist[v as usize] = u64::MAX;
        }
    }

    /// Dijkstra's algorithm (`dijkstra1959note`) with a binary heap.
    fn dijkstra(&mut self, p: &Part, s: usize, dist: &mut [f64]) {
        self.heap.clear();
        dist[s] = 0.0;
        self.heap.push(HeapItem(0.0, s as u32));
        // `idist` marks settled nodes here.
        self.queue.clear();
        while let Some(HeapItem(d, v)) = self.heap.pop() {
            let v = v as usize;
            if self.idist[v] == 0 || d > dist[v] {
                continue;
            }
            self.idist[v] = 0;
            self.queue.push(v as u32);
            for a in p.arcs(v) {
                let t = p.targets[a] as usize;
                let nd = d + p.lens[a];
                if nd < dist[t] {
                    dist[t] = nd;
                    self.heap.push(HeapItem(nd, t as u32));
                }
            }
        }
        for &v in &self.queue {
            self.idist[v as usize] = u64::MAX;
        }
    }
}

/// Fills `out[i] = f(i)` in parallel. Each element is computed on exactly
/// one thread, so the result does not depend on the thread count.
pub(crate) fn par_fill<T, F>(out: &mut [T], threads: usize, min_chunk: usize, f: F)
where
    T: Send,
    F: Fn(usize, &mut T) + Sync,
{
    let n = out.len();
    let t = parallel::resolve_threads(threads, n / min_chunk.max(1));
    if t <= 1 {
        for (i, o) in out.iter_mut().enumerate() {
            f(i, o);
        }
        return;
    }
    let chunk = (n + t - 1) / t;
    std::thread::scope(|scope| {
        let f = &f;
        for (c, block) in out.chunks_mut(chunk).enumerate() {
            scope.spawn(move || {
                let start = c * chunk;
                for (k, o) in block.iter_mut().enumerate() {
                    f(start + k, o);
                }
            });
        }
    });
}

/// Centers a component's coordinates, rotates its principal axis to the
/// horizontal, and reflects it so its first (lowest-numbered) node lies at
/// x <= 0 and y >= 0.
pub(crate) fn orient(xy: &mut [[f64; 2]]) {
    let n = xy.len();
    if n == 0 {
        return;
    }
    let (mut mx, mut my) = (0.0, 0.0);
    for p in xy.iter() {
        mx += p[0];
        my += p[1];
    }
    mx /= n as f64;
    my /= n as f64;
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    for p in xy.iter_mut() {
        p[0] -= mx;
        p[1] -= my;
        sxx += p[0] * p[0];
        syy += p[1] * p[1];
        sxy += p[0] * p[1];
    }
    let theta = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let (s, c) = theta.sin_cos();
    for p in xy.iter_mut() {
        let (x, y) = (p[0], p[1]);
        p[0] = x * c + y * s;
        p[1] = -x * s + y * c;
    }
    let scale = xy
        .iter()
        .fold(0.0f64, |m, p| m.max(p[0].abs()).max(p[1].abs()));
    let tiny = 1e-9 * scale;
    let flip_x = xy[0][0] > tiny;
    let flip_y = xy[0][1] < -tiny;
    for p in xy.iter_mut() {
        if flip_x {
            p[0] = -p[0];
        }
        if flip_y {
            p[1] = -p[1];
        }
    }
}

/// Packs laid-out components into rows (next-fit decreasing height,
/// `coffman1980performance`) with a gap of `gap` around each, and centers
/// the result at the origin. `pieces[k]` holds the global node indices of a
/// component and their coordinates.
pub(crate) fn pack(n: usize, pieces: &[(Vec<u32>, Vec<[f64; 2]>)], gap: f64) -> Layout {
    let mut out = Layout::zeros(n);
    if pieces.is_empty() {
        return out;
    }
    struct Box {
        minx: f64,
        maxy: f64,
        w: f64,
        h: f64,
    }
    let boxes: Vec<Box> = pieces
        .iter()
        .map(|(_, xy)| {
            let (mut minx, mut maxx, mut miny, mut maxy) = (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            );
            for p in xy {
                minx = minx.min(p[0]);
                maxx = maxx.max(p[0]);
                miny = miny.min(p[1]);
                maxy = maxy.max(p[1]);
            }
            Box {
                minx,
                maxy,
                w: maxx - minx + gap,
                h: maxy - miny + gap,
            }
        })
        .collect();
    let mut order: Vec<usize> = (0..pieces.len()).collect();
    // Stable: equal heights keep the order of their lowest-numbered nodes.
    order.sort_by(|&a, &b| boxes[b].h.total_cmp(&boxes[a].h));
    let area: f64 = boxes.iter().map(|b| b.w * b.h).sum();
    let widest = boxes.iter().fold(0.0f64, |m, b| m.max(b.w));
    let width = widest.max(area.sqrt());

    let (mut cursor, mut top, mut shelf) = (0.0f64, 0.0f64, 0.0f64);
    let (mut minx, mut maxx, mut miny, mut maxy) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for (k, &i) in order.iter().enumerate() {
        let b = &boxes[i];
        if k == 0 {
            shelf = b.h;
        } else if cursor + b.w > width * (1.0 + 1e-12) {
            top -= shelf;
            cursor = 0.0;
            shelf = b.h;
        }
        let dx = cursor + gap / 2.0 - b.minx;
        let dy = top - gap / 2.0 - b.maxy;
        let (nodes, xy) = &pieces[i];
        for (&v, p) in nodes.iter().zip(xy) {
            let (x, y) = (p[0] + dx, p[1] + dy);
            out.x[v as usize] = x;
            out.y[v as usize] = y;
            minx = minx.min(x);
            maxx = maxx.max(x);
            miny = miny.min(y);
            maxy = maxy.max(y);
        }
        cursor += b.w;
    }
    let (cx, cy) = ((minx + maxx) / 2.0, (miny + maxy) / 2.0);
    for v in 0..n {
        out.x[v] -= cx;
        out.y[v] -= cy;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parts_merge_parallel_arcs_and_drop_loops() {
        let g = Graph::from_edges(
            4,
            &[0, 0, 1, 2],
            &[1, 1, 1, 3],
            Some(&[2.0, 3.0, 1.0, 4.0]),
            true,
        )
        .unwrap();
        let (p, scale) = parts(&g, true, WeightKind::Distance).unwrap();
        assert_eq!(p.len(), 2);
        assert_eq!(p[0].nodes, vec![0, 1]);
        assert_eq!(p[0].targets, vec![1, 0]);
        assert_eq!(p[0].lens, vec![2.0, 2.0]);
        assert!((p[0].strength[0] - (0.5 + 1.0 / 3.0)).abs() < 1e-15);
        assert_eq!(p[1].nodes, vec![2, 3]);
        assert_eq!(scale, (2.0 + 2.0 + 4.0 + 4.0) / 4.0);
    }

    #[test]
    fn shortest_paths_agree_across_algorithms() {
        // A weighted cycle with a chord; lengths are integers, so Dial,
        // Dijkstra, and (for unit lengths) BFS can be compared exactly.
        let s = [0, 1, 2, 3, 4, 0];
        let t = [1, 2, 3, 4, 0, 2];
        let w = [1.0, 2.0, 1.0, 3.0, 1.0, 5.0];
        let g = Graph::from_edges(5, &s, &t, Some(&w), false).unwrap();
        let (mut p, _) = parts(&g, true, WeightKind::Distance).unwrap();
        assert_eq!(p[0].metric, Metric::Integer(5));
        let mut ws = Sssp::new(5);
        let mut a = vec![0.0; 5];
        ws.run(&p[0], 0, &mut a);
        assert_eq!(a, vec![0.0, 1.0, 3.0, 4.0, 1.0]);
        p[0].metric = Metric::Real;
        let mut b = vec![0.0; 5];
        ws.run(&p[0], 0, &mut b);
        assert_eq!(a, b);
        let (q, _) = parts(&g, false, WeightKind::Distance).unwrap();
        ws.run(&q[0], 0, &mut a);
        assert_eq!(a, vec![0.0, 1.0, 1.0, 2.0, 1.0]);
    }

    #[test]
    fn packing_separates_components() {
        let pieces = vec![
            (vec![0, 1], vec![[-1.0, 0.0], [1.0, 0.0]]),
            (vec![2], vec![[0.0, 0.0]]),
            (vec![3], vec![[0.0, 0.0]]),
        ];
        let l = pack(4, &pieces, 1.0);
        // Boxes are disjoint: the isolates are at least the gap away from
        // the pair and from each other.
        let d = |a: usize, b: usize| ((l.x[a] - l.x[b]).powi(2) + (l.y[a] - l.y[b]).powi(2)).sqrt();
        assert!((d(0, 1) - 2.0).abs() < 1e-12);
        assert!(d(2, 3) >= 1.0 - 1e-12);
        assert!(d(0, 2) >= 1.0 - 1e-12 && d(1, 2) >= 1.0 - 1e-12);
    }
}
