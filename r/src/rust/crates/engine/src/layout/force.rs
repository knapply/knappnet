//! Force-directed placement (Fruchterman and Reingold 1991,
//! `fruchterman1991graph`).
//!
//! Nodes repel each other with force k^2 / e and adjacent nodes attract
//! with force e^2 / k, where e is their distance and k = 1 is the ideal
//! edge length (the paper's k = C sqrt(area / n) with C = 1 and area n).
//! Each iteration computes every node's displacement from the current
//! positions and then moves each node along its displacement by at most
//! the current temperature. The temperature starts at one tenth of the
//! width sqrt(n) of the starting square (the paper's W / 10) and decreases
//! linearly to zero over the iterations. The paper's frame walls are not
//! used. Starting positions are uniform in a square of side sqrt(n), drawn
//! from [`crate::rng::Stream`] with the seed.
//!
//! Repulsion between all pairs costs O(n^2) per iteration. With
//! [`FrApproximation::BarnesHut`] it is approximated with a quadtree
//! (Barnes and Hut 1986, `barnes1986hierarchical`): a cell of width s whose
//! center of mass lies at distance r from the node acts as a single body of
//! its total mass when s / r < theta, so each iteration costs O(n log n).
//! The approximation changes only the simulated forces, and so the
//! positions; nothing computed from a network depends on a layout. The
//! quadtree splits cells until each holds one node (or 48 levels, which
//! only coincident nodes reach). Fruchterman and Reingold's own grid
//! variant, which ignores repulsion beyond 2k, is not implemented.
//!
//! Conventions are those of [`crate::layout`]: direction and self-loops are
//! ignored, components are laid out separately and packed, and the layout
//! is reproducible for a given seed on any number of threads. Parallel
//! edges count once without weights. With weights, the attraction between
//! two nodes is multiplied by the strength of their tie: the sum over
//! parallel edges of 1 / length, where length is the weight read as a
//! distance or the reciprocal of a strength (as igraph multiplies
//! attraction by edge weights read as strengths). Two nodes at the same
//! position are pushed apart in a direction fixed by their indices.

use super::{orient, pack, par_fill, parts, Layout, Part};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::rng::Stream;
use crate::weights::WeightKind;

/// How repulsion between all pairs of nodes is computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FrApproximation {
    /// Exact up to [`FR_EXACT_AUTO_MAX`] nodes per component, Barnes-Hut
    /// beyond.
    #[default]
    Auto,
    /// Every pair, O(n^2) per iteration.
    Exact,
    /// A Barnes-Hut quadtree, O(n log n) per iteration.
    BarnesHut,
}

impl FrApproximation {
    /// Decodes the C interface codes: 0 = auto, 1 = exact, 2 = Barnes-Hut.
    pub fn from_code(code: i32) -> Result<FrApproximation> {
        match code {
            0 => Ok(FrApproximation::Auto),
            1 => Ok(FrApproximation::Exact),
            2 => Ok(FrApproximation::BarnesHut),
            _ => Err(Error::InvalidArgument(format!(
                "unknown force approximation code {code}"
            ))),
        }
    }
}

/// Largest component that [`FrApproximation::Auto`] simulates exactly.
pub const FR_EXACT_AUTO_MAX: usize = 1000;

/// Options for [`fr_layout`].
#[derive(Debug, Clone, Copy)]
pub struct FrOptions {
    /// Use edge weights (they must be strictly positive).
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Number of iterations (default 500); 0 uses the default.
    pub iterations: usize,
    /// Exact or Barnes-Hut repulsion.
    pub approximation: FrApproximation,
    /// Barnes-Hut opening angle (default 1); 0 or less uses the default.
    pub theta: f64,
    /// Seed for the starting positions.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for FrOptions {
    fn default() -> FrOptions {
        FrOptions {
            weighted: false,
            weights_as: WeightKind::Distance,
            iterations: 0,
            approximation: FrApproximation::Auto,
            theta: 0.0,
            seed: 0,
            threads: 0,
        }
    }
}

const DEFAULT_ITER: usize = 500;
const DEFAULT_THETA: f64 = 1.0;
/// Distances below this (in units of k) are clamped so forces stay finite.
const MIN_DIST: f64 = 1e-3;
const MAX_DEPTH: usize = 48;

/// A force-directed layout. See the module documentation.
pub fn fr_layout(g: &Graph, opts: &FrOptions) -> Result<Layout> {
    let n = g.node_count();
    let (parts, _) = parts(g, opts.weighted, opts.weights_as)?;
    let iterations = if opts.iterations == 0 {
        DEFAULT_ITER
    } else {
        opts.iterations
    };
    let theta = if opts.theta > 0.0 {
        opts.theta
    } else {
        DEFAULT_THETA
    };
    let mut pieces = Vec::with_capacity(parts.len());
    for p in parts {
        let bh = match opts.approximation {
            FrApproximation::Auto => p.len() > FR_EXACT_AUTO_MAX,
            FrApproximation::Exact => false,
            FrApproximation::BarnesHut => true,
        };
        let threads = if p.len() < 256 { 1 } else { opts.threads };
        let mut xy = simulate(&p, iterations, bh, theta, opts.seed, threads);
        orient(&mut xy);
        pieces.push((p.nodes, xy));
    }
    Ok(pack(n, &pieces, 1.0))
}

/// A unit direction for two coincident nodes, fixed by their indices and
/// opposite for the two of them.
fn tie_direction(i: usize, j: usize) -> (f64, f64) {
    let (a, b) = if i < j { (i, j) } else { (j, i) };
    // The golden angle spreads directions evenly.
    let angle = 2.399_963_229_728_653 * ((a as f64) * 1.618_033_988_749_895 + b as f64);
    let (s, c) = angle.sin_cos();
    if i < j {
        (c, s)
    } else {
        (-c, -s)
    }
}

/// Repulsion on node i from a body of mass `m` at offset (dx, dy) = x_i -
/// x_body: m k^2 / e along the offset, with k = 1.
#[inline]
fn repel(dx: f64, dy: f64, m: f64, i: usize, j: usize) -> (f64, f64) {
    let e2 = dx * dx + dy * dy;
    if e2 < MIN_DIST * MIN_DIST {
        let (ux, uy) = if e2 > 0.0 {
            let e = e2.sqrt();
            (dx / e, dy / e)
        } else {
            tie_direction(i, j)
        };
        return (ux * m / MIN_DIST, uy * m / MIN_DIST);
    }
    (dx * m / e2, dy * m / e2)
}

fn simulate(
    p: &Part,
    iterations: usize,
    barnes_hut: bool,
    theta: f64,
    seed: u64,
    threads: usize,
) -> Vec<[f64; 2]> {
    let c = p.len();
    if c == 1 {
        return vec![[0.0, 0.0]];
    }
    let side = (c as f64).sqrt();
    let mut stream = Stream::new(seed, p.nodes[0] as u64);
    let mut z: Vec<[f64; 2]> = (0..c)
        .map(|_| {
            let x = side * (stream.uniform() - 0.5);
            let y = side * (stream.uniform() - 0.5);
            [x, y]
        })
        .collect();
    let mut disp = vec![[0.0; 2]; c];
    let mut tree = QuadTree::default();
    let t0 = side / 10.0;
    for it in 0..iterations {
        let temp = t0 * (1.0 - it as f64 / iterations as f64);
        if barnes_hut {
            tree.build(&z);
        }
        par_fill(&mut disp, threads, 64, |i, out| {
            let zi = z[i];
            let (mut fx, mut fy) = if barnes_hut {
                tree.force(&z, i, theta)
            } else {
                let (mut fx, mut fy) = (0.0, 0.0);
                for (j, zj) in z.iter().enumerate() {
                    if j != i {
                        let (ax, ay) = repel(zi[0] - zj[0], zi[1] - zj[1], 1.0, i, j);
                        fx += ax;
                        fy += ay;
                    }
                }
                (fx, fy)
            };
            for a in p.arcs(i) {
                let zj = z[p.targets[a] as usize];
                let (dx, dy) = (zi[0] - zj[0], zi[1] - zj[1]);
                // Attraction e^2 / k along -offset/e, times the tie strength.
                let e = (dx * dx + dy * dy).sqrt();
                fx -= dx * e * p.strength[a];
                fy -= dy * e * p.strength[a];
            }
            *out = [fx, fy];
        });
        for (q, d) in z.iter_mut().zip(&disp) {
            let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
            if len > 0.0 {
                let step = len.min(temp) / len;
                q[0] += d[0] * step;
                q[1] += d[1] * step;
            }
        }
    }
    z
}

#[derive(Debug, Clone, Copy)]
struct Cell {
    /// Center of mass.
    cx: f64,
    cy: f64,
    mass: f64,
    /// Cell width.
    width: f64,
    /// Children, or u32::MAX; a leaf has none and owns `points[first..first + count]`.
    child: [u32; 4],
    first: u32,
    count: u32,
}

#[derive(Default)]
struct QuadTree {
    cells: Vec<Cell>,
    points: Vec<u32>,
}

impl QuadTree {
    fn build(&mut self, z: &[[f64; 2]]) {
        self.cells.clear();
        self.points.clear();
        self.points.extend(0..z.len() as u32);
        let (mut x0, mut x1, mut y0, mut y1) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for q in z {
            x0 = x0.min(q[0]);
            x1 = x1.max(q[0]);
            y0 = y0.min(q[1]);
            y1 = y1.max(q[1]);
        }
        let width = (x1 - x0).max(y1 - y0).max(MIN_DIST);
        let mut points = std::mem::take(&mut self.points);
        self.split(z, &mut points, 0, x0, y0, width, 0);
        self.points = points;
    }

    /// Builds the cell for `points[lo..lo + len]` inside the square with
    /// lower-left corner (x0, y0) and side `width`; returns its index.
    #[allow(clippy::too_many_arguments)]
    fn split(
        &mut self,
        z: &[[f64; 2]],
        points: &mut [u32],
        lo: usize,
        x0: f64,
        y0: f64,
        width: f64,
        depth: usize,
    ) -> u32 {
        let idx = self.cells.len() as u32;
        let len = points.len();
        let (mut cx, mut cy) = (0.0, 0.0);
        for &v in points.iter() {
            cx += z[v as usize][0];
            cy += z[v as usize][1];
        }
        self.cells.push(Cell {
            cx: cx / len as f64,
            cy: cy / len as f64,
            mass: len as f64,
            width,
            child: [u32::MAX; 4],
            first: lo as u32,
            count: len as u32,
        });
        if len <= 1 || depth >= MAX_DEPTH {
            return idx;
        }
        let half = width / 2.0;
        let (mx, my) = (x0 + half, y0 + half);
        // Partition into quadrants: [x < mx, y < my], [x >= mx, y < my],
        // [x < mx, y >= my], [x >= mx, y >= my].
        let south = partition(points, |v| z[v as usize][1] < my);
        let (s, nrth) = points.split_at_mut(south);
        let sw = partition(s, |v| z[v as usize][0] < mx);
        let nw = partition(nrth, |v| z[v as usize][0] < mx);
        let bounds = [
            (0, sw, x0, y0),
            (sw, south, mx, y0),
            (south, south + nw, x0, my),
            (south + nw, len, mx, my),
        ];
        let mut child = [u32::MAX; 4];
        for (q, &(a, b, qx, qy)) in bounds.iter().enumerate() {
            if b > a {
                child[q] = self.split(z, &mut points[a..b], lo + a, qx, qy, half, depth + 1);
            }
        }
        self.cells[idx as usize].child = child;
        idx
    }

    /// Approximate total repulsion on node i.
    fn force(&self, z: &[[f64; 2]], i: usize, theta: f64) -> (f64, f64) {
        let (xi, yi) = (z[i][0], z[i][1]);
        let (mut fx, mut fy) = (0.0, 0.0);
        let mut stack = [0u32; 4 * MAX_DEPTH + 4];
        let mut top = 1;
        stack[0] = 0;
        while top > 0 {
            top -= 1;
            let cell = &self.cells[stack[top] as usize];
            let leaf = cell.child.iter().all(|&ch| ch == u32::MAX);
            if leaf {
                let first = cell.first as usize;
                for &j in &self.points[first..first + cell.count as usize] {
                    let j = j as usize;
                    if j != i {
                        let (ax, ay) = repel(xi - z[j][0], yi - z[j][1], 1.0, i, j);
                        fx += ax;
                        fy += ay;
                    }
                }
                continue;
            }
            let (dx, dy) = (xi - cell.cx, yi - cell.cy);
            let r2 = dx * dx + dy * dy;
            if cell.width * cell.width < theta * theta * r2 {
                let (ax, ay) = repel(dx, dy, cell.mass, i, i);
                fx += ax;
                fy += ay;
            } else {
                for &ch in &cell.child {
                    if ch != u32::MAX {
                        stack[top] = ch;
                        top += 1;
                    }
                }
            }
        }
        (fx, fy)
    }
}

/// Moves the elements satisfying `pred` to the front (in place, so the
/// order within each group is deterministic but not preserved); returns
/// how many satisfy it.
fn partition(v: &mut [u32], pred: impl Fn(u32) -> bool) -> usize {
    let mut k = 0;
    for i in 0..v.len() {
        if pred(v[i]) {
            v.swap(i, k);
            k += 1;
        }
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring(n: usize) -> Graph {
        let s: Vec<i64> = (0..n as i64).collect();
        let t: Vec<i64> = (0..n as i64).map(|i| (i + 1) % n as i64).collect();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn barnes_hut_force_is_close_to_exact() {
        let mut stream = Stream::new(3, 0);
        let z: Vec<[f64; 2]> = (0..400)
            .map(|_| [10.0 * stream.uniform(), 10.0 * stream.uniform()])
            .collect();
        let mut tree = QuadTree::default();
        tree.build(&z);
        for i in [0, 17, 399] {
            let (ex, ey) = z.iter().enumerate().filter(|(j, _)| *j != i).fold(
                (0.0, 0.0),
                |(fx, fy), (j, zj)| {
                    let (ax, ay) = repel(z[i][0] - zj[0], z[i][1] - zj[1], 1.0, i, j);
                    (fx + ax, fy + ay)
                },
            );
            let (bx, by) = tree.force(&z, i, 0.5);
            let err = ((ex - bx).powi(2) + (ey - by).powi(2)).sqrt();
            let size = (ex * ex + ey * ey).sqrt();
            assert!(err <= 0.05 * size + 0.05, "{err} vs {size}");
            // theta near 0 opens every cell: exact up to summation order.
            let (cx, cy) = tree.force(&z, i, 1e-9);
            assert!((cx - ex).abs() < 1e-9 * (1.0 + ex.abs()));
            assert!((cy - ey).abs() < 1e-9 * (1.0 + ey.abs()));
        }
    }

    #[test]
    fn a_ring_has_roughly_even_edges() {
        let g = ring(12);
        let l = fr_layout(&g, &FrOptions::default()).unwrap();
        let lens: Vec<f64> = (0..12)
            .map(|i| {
                let j = (i + 1) % 12;
                ((l.x[i] - l.x[j]).powi(2) + (l.y[i] - l.y[j]).powi(2)).sqrt()
            })
            .collect();
        let mean = lens.iter().sum::<f64>() / 12.0;
        assert!(
            lens.iter().all(|e| (e - mean).abs() < 0.2 * mean),
            "{lens:?}"
        );
    }

    #[test]
    fn coincident_nodes_separate() {
        let g = Graph::from_edges(3, &[], &[], None, false).unwrap();
        let z = [[0.0, 0.0]; 3];
        let (ax, ay) = repel(z[0][0] - z[1][0], z[0][1] - z[1][1], 1.0, 0, 1);
        let (bx, by) = repel(z[1][0] - z[0][0], z[1][1] - z[0][1], 1.0, 1, 0);
        assert!((ax + bx).abs() < 1e-12 && (ay + by).abs() < 1e-12);
        assert!(ax.abs() + ay.abs() > 0.0);
        assert!(fr_layout(&g, &FrOptions::default()).is_ok());
    }
}
