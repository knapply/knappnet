//! Stress layouts: stress majorization on all pairs, and sparse stress with
//! pivots for large components.
//!
//! Sources:
//! - Stress majorization (Gansner, Koren, and North 2005,
//!   `gansner2005graph`) minimizes stress(X) = sum_{i<j} w_ij
//!   (||x_i - x_j|| - d_ij)^2 with w_ij = d_ij^-2 (Kamada and Kawai 1989,
//!   `kamada1989algorithm`). Holding every node but i fixed, the majorizing
//!   function of Gansner et al. is minimized by the localized update
//!
//!   ```text
//!   J_i = sum_{j != i} w_ij (x_j + d_ij (x_i - x_j) / ||x_i - x_j||) / sum_{j != i} w_ij.
//!   ```
//!
//!   Applying J to all nodes at once (a Jacobi sweep) lets every node be
//!   updated in parallel, but a Jacobi sweep is not itself a majorization
//!   step and can occasionally increase stress. When it does, the sweep is
//!   replaced by the half step (x + J(x)) / 2, which minimizes a separable
//!   majorizer (it bounds the quadratic part of stress with the Laplacian
//!   inequality L_w <= 2 D_w) and therefore never increases stress. Stress
//!   is thus non-increasing, and iteration stops when an iteration lowers
//!   it by less than a relative `tol` (default 1e-4, the value used by
//!   Gansner et al.) or after `max_iter` iterations (default 500). This
//!   needs every shortest-path length, so time per iteration and memory are
//!   O(n^2) in the component size.
//! - Stress has local minima. Components of at most 100 nodes are therefore
//!   also laid out from 4 uniformly random starts (drawn from the seed) and
//!   the layout with the least stress is kept, with a tolerance of 1e-5.
//!   On the validation fixtures this matches or improves on NetworkX's
//!   `kamada_kawai_layout()`, which minimizes the same energy
//!   (`validation/layout_check.py`).
//! - Sparse stress (Ortmann, Klimenta, and Brandes 2017,
//!   `ortmann2017sparse`) keeps the edge terms and replaces the other pairs
//!   by terms to k pivots (default 100): node i is pulled toward distance
//!   d_ip from pivot p (for pivots that are not i's neighbors) with weight
//!   s / d_ip^2, where s counts the nodes of p's region (the nodes closer
//!   to p than to any other pivot, ties going to the currently smallest
//!   region) whose distance to p is at most d_ip / 2. Following their
//!   Algorithm 1, the initial layout is scaled so the edge lengths sum to
//!   the sum of the edges' target distances. Pivots are a uniform random
//!   sample (one of the strategies Ortmann et al. compare; they found
//!   k-means sampling slightly better, and max/min sampling, used here for
//!   the all-pairs start, gave much worse sparse layouts in our tests).
//!   Ortmann et al. update nodes one at a time. Here nodes are updated in
//!   classes that can be updated at once without one node reading
//!   another's new position (each pivot alone, then a greedy coloring of
//!   the other nodes in node order), a multicolor Gauss-Seidel sweep that
//!   runs in parallel and does not depend on the thread count. The pivot
//!   terms are not symmetric, so the whole layout can keep turning slowly
//!   once its shape has settled; iteration therefore stops when the
//!   movement left after removing the best rigid motion (orthogonal
//!   Procrustes), relative to the sum of distances from the centroid,
//!   falls below `tol` (default 1e-4, as in their Algorithm 1), or after
//!   `max_iter` sweeps (default 100; in our tests stress hardly changes
//!   after the first 50). Pivot weights and distances are stored in single
//!   precision. Time per sweep is O(k n + m) and memory O(k n).
//! - Both start from pivot multidimensional scaling (Brandes and Pich 2007,
//!   `brandes2007eigensolver`), as recommended by Brandes and Pich (2009,
//!   `brandes2009experimental`). For all-pairs stress the 50 pivots are
//!   chosen by the max/min rule of Brandes and Pich: the first is the
//!   lowest-numbered node of largest degree, and each next pivot is the
//!   node farthest from those chosen so far (lowest-numbered on ties). The
//!   two leading eigenvectors of C'C, where C is the double-centered matrix
//!   of squared distances to the pivots, come from power iteration with
//!   deflation (`mises1929praktische`). Nodes with identical distances to
//!   every pivot coincide in that start, so every coordinate is then moved
//!   by a uniform random amount of at most 1% of the typical edge length,
//!   drawn from the seed.
//!
//! Conventions are those of [`crate::layout`]: direction and self-loops are
//! ignored, parallel edges give their shortest length, weights are
//! distances or strengths (lengths 1/w), and components are laid out
//! separately and packed. Coordinates are in units of path length.
//!
//! [`StressAlgorithm::Auto`] uses all pairs for components of up to
//! [`super::FULL_STRESS_AUTO_MAX`] nodes and sparse stress for larger ones.
//! Sparse stress approximates the stress objective, not any reported
//! measure; [`layout_stress`] always evaluates the exact stress.

use super::{orient, pack, par_fill, parts, Layout, Part, Sssp, FULL_STRESS_AUTO_MAX};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;
use crate::rng::Stream;
use crate::weights::WeightKind;

/// Largest component for which [`StressAlgorithm::Full`] is allowed: its
/// distance matrix takes 8 n^2 bytes (800 MB at this size).
pub const FULL_STRESS_MAX: usize = 10_000;

/// Which stress model to optimize.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StressAlgorithm {
    /// All pairs up to [`FULL_STRESS_AUTO_MAX`] nodes per component, sparse
    /// stress beyond.
    #[default]
    Auto,
    /// Stress majorization on all pairs (Gansner, Koren, and North 2005).
    Full,
    /// Sparse stress with pivots (Ortmann, Klimenta, and Brandes 2017).
    Sparse,
}

impl StressAlgorithm {
    /// Decodes the C interface codes: 0 = auto, 1 = all pairs, 2 = sparse.
    pub fn from_code(code: i32) -> Result<StressAlgorithm> {
        match code {
            0 => Ok(StressAlgorithm::Auto),
            1 => Ok(StressAlgorithm::Full),
            2 => Ok(StressAlgorithm::Sparse),
            _ => Err(Error::InvalidArgument(format!(
                "unknown stress algorithm code {code}"
            ))),
        }
    }
}

/// Options for [`stress_layout`].
#[derive(Debug, Clone, Copy)]
pub struct StressOptions {
    /// Use edge weights (they must be strictly positive).
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// All pairs, sparse, or chosen by component size.
    pub algorithm: StressAlgorithm,
    /// Number of pivots for sparse stress (default 100) and for the pivot
    /// MDS start of all-pairs stress (default 50); 0 uses the default.
    pub pivots: usize,
    /// Iteration cap (default 500 for all pairs, 100 for sparse stress);
    /// 0 uses the default.
    pub max_iter: usize,
    /// Relative convergence tolerance; 0 or less uses the default: 1e-5
    /// for all-pairs stress on components of at most 100 nodes, 1e-4
    /// otherwise.
    pub tol: f64,
    /// Seed for the perturbation of the starting layout, the random
    /// restarts of small components, and the sample of sparse-stress
    /// pivots.
    pub seed: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for StressOptions {
    fn default() -> StressOptions {
        StressOptions {
            weighted: false,
            weights_as: WeightKind::Distance,
            algorithm: StressAlgorithm::Auto,
            pivots: 0,
            max_iter: 0,
            tol: 0.0,
            seed: 0,
            threads: 0,
        }
    }
}

const DEFAULT_ITER: usize = 500;
const SPARSE_ITER: usize = 100;
const DEFAULT_TOL: f64 = 1e-4;
/// Tolerance for components small enough to be restarted.
const SMALL_TOL: f64 = 1e-5;
const MDS_PIVOTS: usize = 50;
const SPARSE_PIVOTS: usize = 100;
/// Components of at most this many nodes are also laid out from
/// [`RESTARTS`] random starts, keeping the layout with the least stress.
const RESTART_MAX: usize = 100;
const RESTARTS: usize = 4;
/// Components smaller than this are laid out one per thread instead of
/// with parallel loops inside.
const SMALL_PART: usize = 256;

/// A stress layout. See the module documentation.
pub fn stress_layout(g: &Graph, opts: &StressOptions) -> Result<Layout> {
    let n = g.node_count();
    let (parts, scale) = parts(g, opts.weighted, opts.weights_as)?;
    if opts.algorithm == StressAlgorithm::Full {
        if let Some(p) = parts.iter().find(|p| p.len() > FULL_STRESS_MAX) {
            return Err(Error::InvalidArgument(format!(
                "a component has {} nodes; stress on all pairs is limited to {} nodes \
                 per component (use sparse stress)",
                p.len(),
                FULL_STRESS_MAX
            )));
        }
    }
    let threads = opts.threads;
    let one = |p: &Part, threads: usize| -> Vec<[f64; 2]> {
        let mut xy = layout_part(p, opts, threads);
        orient(&mut xy);
        xy
    };

    // Small components run one per thread, large ones one at a time with
    // parallel loops inside. Each component's result is independent of the
    // thread count either way.
    let small: Vec<usize> = (0..parts.len())
        .filter(|&i| parts[i].len() < SMALL_PART)
        .collect();
    let mut coords: Vec<Vec<[f64; 2]>> = vec![Vec::new(); parts.len()];
    let mut done: Vec<Vec<[f64; 2]>> = vec![Vec::new(); small.len()];
    par_fill(&mut done, threads, 4, |k, out| {
        *out = one(&parts[small[k]], 1);
    });
    for (k, xy) in small.iter().zip(done) {
        coords[*k] = xy;
    }
    for (i, p) in parts.iter().enumerate() {
        if p.len() >= SMALL_PART {
            coords[i] = one(p, threads);
        }
    }
    let pieces: Vec<(Vec<u32>, Vec<[f64; 2]>)> = parts
        .into_iter()
        .zip(coords)
        .map(|(p, xy)| (p.nodes, xy))
        .collect();
    Ok(pack(n, &pieces, scale))
}

fn layout_part(p: &Part, opts: &StressOptions, threads: usize) -> Vec<[f64; 2]> {
    let c = p.len();
    match c {
        0 => return Vec::new(),
        1 => return vec![[0.0, 0.0]],
        2 => {
            let d = p.lens[0];
            return vec![[-d / 2.0, 0.0], [d / 2.0, 0.0]];
        }
        _ => {}
    }
    let full = match opts.algorithm {
        StressAlgorithm::Full => true,
        StressAlgorithm::Sparse => false,
        StressAlgorithm::Auto => c <= FULL_STRESS_AUTO_MAX,
    };
    let max_iter = match (opts.max_iter, full) {
        (0, true) => DEFAULT_ITER,
        (0, false) => SPARSE_ITER,
        (m, _) => m,
    };
    let tol = if opts.tol > 0.0 {
        opts.tol
    } else if full && c <= RESTART_MAX {
        SMALL_TOL
    } else {
        DEFAULT_TOL
    };
    let mut stream = Stream::new(opts.seed, p.nodes[0] as u64);
    let jitter = 0.01 * p.lens.iter().sum::<f64>() / p.lens.len() as f64;
    if full {
        let k = if opts.pivots == 0 {
            MDS_PIVOTS
        } else {
            opts.pivots
        };
        full_stress(p, k.min(c), max_iter, tol, jitter, &mut stream, threads)
    } else {
        let k = if opts.pivots == 0 {
            SPARSE_PIVOTS
        } else {
            opts.pivots
        };
        sparse_stress(p, k.min(c), max_iter, tol, jitter, &mut stream, threads)
    }
}

/// Moves every coordinate by a uniform amount in [-a, a].
fn perturb(xy: &mut [[f64; 2]], a: f64, stream: &mut Stream) {
    for p in xy.iter_mut() {
        p[0] += a * (2.0 * stream.uniform() - 1.0);
        p[1] += a * (2.0 * stream.uniform() - 1.0);
    }
}

/// Max/min pivot selection (Brandes and Pich 2007) over the rows of an
/// all-pairs distance matrix.
fn maxmin_pivots(p: &Part, k: usize, dist: &[Vec<f64>]) -> Vec<usize> {
    let c = p.len();
    let mut first = 0;
    for v in 1..c {
        if p.degree(v) > p.degree(first) {
            first = v;
        }
    }
    let mut pivots = Vec::with_capacity(k);
    let mut nearest = vec![f64::INFINITY; c];
    let mut next = first;
    while pivots.len() < k {
        pivots.push(next);
        for (m, &x) in nearest.iter_mut().zip(&dist[next]) {
            if x < *m {
                *m = x;
            }
        }
        let mut best = 0;
        for v in 1..c {
            if nearest[v] > nearest[best] {
                best = v;
            }
        }
        if nearest[best] == 0.0 {
            break;
        }
        next = best;
    }
    pivots
}

/// Pivot MDS (Brandes and Pich 2007). `d(i, j)` is the distance from node
/// i to pivot j.
fn pivot_mds(
    c: usize,
    k: usize,
    d: impl Fn(usize, usize) -> f64 + Sync,
    threads: usize,
) -> Vec<[f64; 2]> {
    // Double-centered squared distances, C_ij = -1/2 (d_ij^2 - r_i - c_j + g).
    let mut col = vec![0.0; k];
    let mut row = vec![0.0; c];
    for (i, r) in row.iter_mut().enumerate() {
        for (j, cj) in col.iter_mut().enumerate() {
            let x = d(i, j) * d(i, j);
            *r += x;
            *cj += x;
        }
    }
    let grand: f64 = row.iter().sum::<f64>() / (c * k) as f64;
    for r in row.iter_mut() {
        *r /= k as f64;
    }
    for cj in col.iter_mut() {
        *cj /= c as f64;
    }
    let entry = |i: usize, j: usize| -0.5 * (d(i, j) * d(i, j) - row[i] - col[j] + grand);

    // C'C summed over fixed blocks of rows, then in block order, so the sum
    // does not depend on the thread count.
    const BLOCK: usize = 1024;
    let blocks = (c + BLOCK - 1) / BLOCK;
    let mut partial: Vec<Vec<f64>> = vec![Vec::new(); blocks];
    par_fill(&mut partial, threads, 1, |b, out| {
        let mut m = vec![0.0; k * k];
        let mut ci = vec![0.0; k];
        for i in b * BLOCK..((b + 1) * BLOCK).min(c) {
            for (j, x) in ci.iter_mut().enumerate() {
                *x = entry(i, j);
            }
            for a in 0..k {
                let ca = ci[a];
                if ca == 0.0 {
                    continue;
                }
                let r = &mut m[a * k..(a + 1) * k];
                for (rb, &cb) in r.iter_mut().zip(&ci) {
                    *rb += ca * cb;
                }
            }
        }
        *out = m;
    });
    let mut m = vec![0.0; k * k];
    for part in &partial {
        for (a, b) in m.iter_mut().zip(part) {
            *a += b;
        }
    }

    let v1 = power_iteration(&m, k, None, 1);
    let lambda1 = rayleigh(&m, k, &v1);
    let mut m2 = m.clone();
    for a in 0..k {
        for b in 0..k {
            m2[a * k + b] -= lambda1 * v1[a] * v1[b];
        }
    }
    let v2 = power_iteration(&m2, k, Some(&v1), 2);
    let mut xy = vec![[0.0; 2]; c];
    par_fill(&mut xy, threads, 1024, |i, out| {
        let (mut x, mut y) = (0.0, 0.0);
        for j in 0..k {
            let e = entry(i, j);
            x += e * v1[j];
            y += e * v2[j];
        }
        *out = [x, y];
    });
    xy
}

fn rayleigh(m: &[f64], k: usize, v: &[f64]) -> f64 {
    let mut s = 0.0;
    for a in 0..k {
        let mut r = 0.0;
        for b in 0..k {
            r += m[a * k + b] * v[b];
        }
        s += v[a] * r;
    }
    s
}

/// Leading eigenvector of a symmetric positive semidefinite k-by-k matrix
/// by power iteration, kept orthogonal to `ortho` if given. Returns the
/// zero vector if the matrix annihilates the start.
///
/// The start is a fixed random vector (a constant stream chosen by `salt`,
/// independent of the user's seed). The second eigenvector must start from
/// a different vector than the first: when the two leading eigenvalues are
/// equal, the first converges to the projection of its start onto their
/// eigenspace, and that start, made orthogonal to it, would have nothing
/// left in the eigenspace.
fn power_iteration(m: &[f64], k: usize, ortho: Option<&[f64]>, salt: u64) -> Vec<f64> {
    let mut start = Stream::new(0x7069_766f_745f_6d64, salt);
    let mut v: Vec<f64> = (0..k).map(|_| start.uniform() - 0.5).collect();
    let project = |v: &mut [f64]| {
        if let Some(u) = ortho {
            let dot: f64 = v.iter().zip(u).map(|(a, b)| a * b).sum();
            for (a, b) in v.iter_mut().zip(u) {
                *a -= dot * b;
            }
        }
    };
    let normalize = |v: &mut [f64]| -> f64 {
        let norm = v.iter().map(|a| a * a).sum::<f64>().sqrt();
        if norm > 0.0 {
            for a in v.iter_mut() {
                *a /= norm;
            }
        }
        norm
    };
    project(&mut v);
    if normalize(&mut v) == 0.0 {
        return vec![0.0; k];
    }
    let mut w = vec![0.0; k];
    for _ in 0..10_000 {
        for a in 0..k {
            let mut s = 0.0;
            for b in 0..k {
                s += m[a * k + b] * v[b];
            }
            w[a] = s;
        }
        project(&mut w);
        if normalize(&mut w) == 0.0 {
            return vec![0.0; k];
        }
        let change = v
            .iter()
            .zip(&w)
            .fold(0.0f64, |acc, (a, b)| acc.max((a - b).abs()));
        std::mem::swap(&mut v, &mut w);
        if change < 1e-10 {
            break;
        }
    }
    v
}

/// Stress majorization on all pairs.
fn full_stress(
    p: &Part,
    k: usize,
    max_iter: usize,
    tol: f64,
    jitter: f64,
    stream: &mut Stream,
    threads: usize,
) -> Vec<[f64; 2]> {
    let c = p.len();
    let t = crate::linalg::phase::start("stress: all-pairs distances");
    let mut dist: Vec<Vec<f64>> = vec![Vec::new(); c];
    par_fill(&mut dist, threads, 16, |s, row| {
        let mut ws = Sssp::new(c);
        let mut d = vec![0.0; c];
        ws.run(p, s, &mut d);
        *row = d;
    });
    drop(t);
    let t = crate::linalg::phase::start("stress: pivot mds start");
    let pivots = maxmin_pivots(p, k, &dist);
    let mut z = pivot_mds(c, pivots.len(), |i, j| dist[pivots[j]][i], threads);
    drop(t);
    // Reciprocal distances; the update and the stress need only these.
    let mut inv = dist;
    par_fill(&mut inv, threads, 16, |i, row| {
        for (j, x) in row.iter_mut().enumerate() {
            *x = if j == i { 0.0 } else { 1.0 / *x };
        }
    });

    // Scales a start to minimize stress over all multiples of it.
    let rescale = |z: &mut Vec<[f64; 2]>| {
        let mut ratio = vec![[0.0; 2]; c];
        par_fill(&mut ratio, threads, 64, |i, out| {
            let (mut a, mut b) = (0.0, 0.0);
            for j in i + 1..c {
                let e = ((z[i][0] - z[j][0]).powi(2) + (z[i][1] - z[j][1]).powi(2)).sqrt();
                let r = inv[i][j];
                a += e * r;
                b += e * e * r * r;
            }
            *out = [a, b];
        });
        let (a, b) = ratio
            .iter()
            .fold((0.0, 0.0), |(a, b), r| (a + r[0], b + r[1]));
        if b > 0.0 {
            let s = a / b;
            for q in z.iter_mut() {
                q[0] *= s;
                q[1] *= s;
            }
        }
    };

    // One pass: J(z) and each node's share of the stress of z. With
    // r = 1 / d and w = r^2, a term contributes w x_j + r (x_i - x_j) / e
    // to the numerator, w to the denominator, and (e r - 1)^2 to stress.
    //
    // Rows are taken LANES at a time (see [`stress_rows`]): each row's sums
    // still run over j in order, so the values do not depend on the
    // grouping.
    let pass = |z: &[[f64; 2]], j: &mut Vec<[f64; 3]>| -> f64 {
        let mut tiles = vec![[[0.0; 3]; LANES]; (c + LANES - 1) / LANES];
        par_fill(&mut tiles, threads, 16, |g, out| {
            let row = |l: usize| (g * LANES + l).min(c - 1);
            let rows: [&[f64]; LANES] = std::array::from_fn(|l| &inv[row(l)][..]);
            let zi: [[f64; 2]; LANES] = std::array::from_fn(|l| z[row(l)]);
            *out = crate::linalg::simd::dispatch!(stress_rows(z, &rows, &zi));
        });
        for (i, r) in j.iter_mut().enumerate() {
            *r = tiles[i / LANES][i % LANES];
        }
        // Each pair is counted from both ends.
        j.iter().map(|r| r[2]).sum::<f64>() / 2.0
    };

    // Majorization from z; returns the final layout and its stress.
    let majorize = |mut z: Vec<[f64; 2]>| -> (Vec<[f64; 2]>, f64) {
        let _t = crate::linalg::phase::start("stress: majorization");
        let mut jz = vec![[0.0; 3]; c];
        let mut jy = vec![[0.0; 3]; c];
        let mut sz = pass(&z, &mut jz);
        let mut y = vec![[0.0; 2]; c];
        for _ in 0..max_iter {
            for (q, r) in y.iter_mut().zip(&jz) {
                *q = [r[0], r[1]];
            }
            let mut sy = pass(&y, &mut jy);
            if sy > sz {
                // The Jacobi sweep raised stress; the half step cannot.
                for ((q, r), zi) in y.iter_mut().zip(&jz).zip(&z) {
                    *q = [(zi[0] + r[0]) / 2.0, (zi[1] + r[1]) / 2.0];
                }
                sy = pass(&y, &mut jy);
            }
            let converged = sz - sy <= tol * sz;
            std::mem::swap(&mut z, &mut y);
            std::mem::swap(&mut jz, &mut jy);
            sz = sy;
            if converged {
                break;
            }
        }
        (z, sz)
    };

    rescale(&mut z);
    perturb(&mut z, jitter, stream);
    let (mut best, mut best_stress) = majorize(z);
    let restarts = if c <= RESTART_MAX { RESTARTS } else { 0 };
    let side = (c as f64).sqrt();
    for _ in 0..restarts {
        let mut z: Vec<[f64; 2]> = (0..c)
            .map(|_| [side * stream.uniform(), side * stream.uniform()])
            .collect();
        rescale(&mut z);
        let (l, s) = majorize(z);
        if s < best_stress {
            best = l;
            best_stress = s;
        }
    }
    best
}

/// Rows per group in the majorization pass of [`full_stress`].
const LANES: usize = 4;

/// For LANES nodes i at once (positions `zi`, reciprocal distance rows
/// `rows`), the update J(z)_i and the node's share of the stress: with
/// r = 1 / d and w = r^2, each j contributes w x_j + r (x_i - x_j) / e to
/// the numerator, w to the denominator, and (e r - 1)^2 to stress. The
/// rows are independent and each one's sums run over j in order, as for a
/// single row, so the compiler may keep the rows in vector lanes without
/// changing any result.
#[inline(always)]
fn stress_rows(
    z: &[[f64; 2]],
    rows: &[&[f64]; LANES],
    zi: &[[f64; 2]; LANES],
) -> [[f64; 3]; LANES] {
    let c = z.len();
    let rows: [&[f64]; LANES] = std::array::from_fn(|l| &rows[l][..c]);
    let (mut nx, mut ny, mut den, mut s) = ([0.0; LANES], [0.0; LANES], [0.0; LANES], [0.0; LANES]);
    for (k, zj) in z.iter().enumerate() {
        for l in 0..LANES {
            let r = rows[l][k];
            let w = r * r;
            let (dx, dy) = (zi[l][0] - zj[0], zi[l][1] - zj[1]);
            let e = (dx * dx + dy * dy).sqrt();
            // The diagonal has r = 0 and adds nothing.
            nx[l] += w * zj[0];
            ny[l] += w * zj[1];
            if e > 0.0 {
                let f = r / e;
                nx[l] += f * dx;
                ny[l] += f * dy;
            }
            den[l] += w;
            let u = e * r - 1.0;
            s[l] += if r > 0.0 { u * u } else { 0.0 };
        }
    }
    std::array::from_fn(|l| [nx[l] / den[l], ny[l] / den[l], s[l]])
}

/// Sparse stress with pivots (Ortmann, Klimenta, and Brandes 2017).
fn sparse_stress(
    p: &Part,
    k: usize,
    max_iter: usize,
    tol: f64,
    jitter: f64,
    stream: &mut Stream,
    threads: usize,
) -> Vec<[f64; 2]> {
    let c = p.len();
    // Pivots: a uniform random sample, searched in parallel.
    let mut perm = vec![0usize; c];
    stream.permutation(&mut perm);
    let pivots: Vec<usize> = perm[..k].to_vec();
    drop(perm);
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); k];
    par_fill(&mut cols, threads, 1, |t, out| {
        let mut ws = Sssp::new(c);
        let mut d = vec![0.0; c];
        ws.run(p, pivots[t], &mut d);
        *out = d;
    });

    // Regions: each node joins its nearest pivot, ties going to the
    // currently smallest region (then the lowest pivot).
    let mut size = vec![0usize; k];
    let mut region = vec![0usize; c];
    for i in 0..c {
        let mut best = 0;
        for t in 1..k {
            let (dt, db) = (cols[t][i], cols[best][i]);
            if dt < db || (dt == db && size[t] < size[best]) {
                best = t;
            }
        }
        region[i] = best;
        size[best] += 1;
    }
    let mut members: Vec<Vec<f64>> = size.iter().map(|&s| Vec::with_capacity(s)).collect();
    for i in 0..c {
        members[region[i]].push(cols[region[i]][i]);
    }
    for m in members.iter_mut() {
        m.sort_by(|a, b| a.total_cmp(b));
    }
    // Node-major terms (w'_ip, d_ip) in single precision, which halves the
    // memory traffic of each iteration. The weight w'_ip = s / d_ip^2 is
    // zero for i = p and for pivots adjacent to i (those pairs are edge
    // terms).
    let mut pivot_of = vec![usize::MAX; c];
    for (t, &v) in pivots.iter().enumerate() {
        pivot_of[v] = t;
    }
    let mut terms = vec![[0.0f32; 2]; c * k];
    par_fill(&mut terms, threads, 4096, |idx, out| {
        let (i, t) = (idx / k, idx % k);
        let d = cols[t][i];
        let w = if pivots[t] == i || d == 0.0 {
            0.0
        } else {
            let s = members[t].partition_point(|&x| x <= d / 2.0) as f64;
            s / (d * d)
        };
        *out = [w as f32, d as f32];
    });
    drop(members);
    drop(cols);
    for i in 0..c {
        for a in p.arcs(i) {
            let t = pivot_of[p.targets[a] as usize];
            if t != usize::MAX {
                terms[i * k + t][0] = 0.0;
            }
        }
    }

    let mut z = pivot_mds(c, k, |i, t| terms[i * k + t][1] as f64, threads);
    // Scale so the edge lengths sum to the edges' target distances.
    let (mut target, mut drawn) = (0.0, 0.0);
    for i in 0..c {
        for a in p.arcs(i) {
            let j = p.targets[a] as usize;
            target += p.lens[a];
            drawn += ((z[i][0] - z[j][0]).powi(2) + (z[i][1] - z[j][1]).powi(2)).sqrt();
        }
    }
    if drawn > 0.0 {
        let s = target / drawn;
        for q in z.iter_mut() {
            q[0] *= s;
            q[1] *= s;
        }
    }
    perturb(&mut z, jitter, stream);

    // Classes of nodes that can be updated together: each pivot on its
    // own (every node reads every pivot), then a greedy coloring of the
    // other nodes in node order, so no two nodes of a class are adjacent.
    let mut classes: Vec<Vec<u32>> = pivots.iter().map(|&v| vec![v as u32]).collect();
    let first_color = classes.len();
    let mut color = vec![usize::MAX; c];
    let mut used: Vec<usize> = Vec::new();
    for i in 0..c {
        if pivot_of[i] != usize::MAX {
            continue;
        }
        used.clear();
        for &t in &p.targets[p.arcs(i)] {
            let ct = color[t as usize];
            if ct != usize::MAX {
                used.push(ct);
            }
        }
        used.sort_unstable();
        used.dedup();
        let mut col = 0;
        for &u in &used {
            if u == col {
                col += 1;
            } else if u > col {
                break;
            }
        }
        color[i] = col;
        if first_color + col >= classes.len() {
            classes.push(Vec::new());
        }
        classes[first_color + col].push(i as u32);
    }

    let update = |z: &[[f64; 2]], i: usize| -> [f64; 2] {
        let zi = z[i];
        let (mut nx, mut ny, mut den) = (0.0, 0.0, 0.0);
        let mut term = |zj: [f64; 2], d: f64, w: f64| {
            let (dx, dy) = (zi[0] - zj[0], zi[1] - zj[1]);
            let e = (dx * dx + dy * dy).sqrt();
            if e > 0.0 {
                nx += w * (zj[0] + d * dx / e);
                ny += w * (zj[1] + d * dy / e);
            } else {
                nx += w * zj[0];
                ny += w * zj[1];
            }
            den += w;
        };
        for a in p.arcs(i) {
            let d = p.lens[a];
            term(z[p.targets[a] as usize], d, 1.0 / (d * d));
        }
        for (t, wd) in terms[i * k..(i + 1) * k].iter().enumerate() {
            if wd[0] > 0.0 {
                term(z[pivots[t]], wd[1] as f64, wd[0] as f64);
            }
        }
        if den > 0.0 {
            [nx / den, ny / den]
        } else {
            zi
        }
    };

    let mut fresh: Vec<[f64; 2]> = Vec::new();
    let mut old = z.clone();
    for _ in 0..max_iter {
        old.copy_from_slice(&z);
        for class in &classes {
            fresh.clear();
            fresh.resize(class.len(), [0.0; 2]);
            par_fill(&mut fresh, threads, 2048, |m, out| {
                *out = update(&z, class[m] as usize);
            });
            for (&v, q) in class.iter().zip(&fresh) {
                z[v as usize] = *q;
            }
        }
        if relative_change(&old, &z) <= tol {
            break;
        }
    }
    z
}

/// Movement from `old` to `new` that is not a rigid motion, relative to the
/// size of `old`: the sum of distances between `new` and `old` after
/// translating and rotating `old` onto `new` by least squares (orthogonal
/// Procrustes in the plane), divided by the sum of distances of `old` from
/// its centroid. The pivot terms of sparse stress are not symmetric, so the
/// whole layout can keep turning slowly after its shape has settled; that
/// turning does not count as change.
fn relative_change(old: &[[f64; 2]], new: &[[f64; 2]]) -> f64 {
    let n = old.len() as f64;
    let (mut ox, mut oy, mut nx, mut ny) = (0.0, 0.0, 0.0, 0.0);
    for (a, b) in old.iter().zip(new) {
        ox += a[0];
        oy += a[1];
        nx += b[0];
        ny += b[1];
    }
    ox /= n;
    oy /= n;
    nx /= n;
    ny /= n;
    let (mut dot, mut cross) = (0.0, 0.0);
    for (a, b) in old.iter().zip(new) {
        let (a0, a1) = (a[0] - ox, a[1] - oy);
        let (b0, b1) = (b[0] - nx, b[1] - ny);
        dot += a0 * b0 + a1 * b1;
        cross += a0 * b1 - a1 * b0;
    }
    let (sin, cos) = cross.atan2(dot).sin_cos();
    let (mut moved, mut spread) = (0.0, 0.0);
    for (a, b) in old.iter().zip(new) {
        let (a0, a1) = (a[0] - ox, a[1] - oy);
        let (r0, r1) = (cos * a0 - sin * a1, sin * a0 + cos * a1);
        moved += ((r0 - (b[0] - nx)).powi(2) + (r1 - (b[1] - ny)).powi(2)).sqrt();
        spread += (a0 * a0 + a1 * a1).sqrt();
    }
    if spread > 0.0 {
        moved / spread
    } else {
        0.0
    }
}

/// The stress of a layout: the sum over pairs of nodes in the same
/// component of (||x_i - x_j|| - d_ij)^2 / d_ij^2, where d_ij is the
/// shortest-path length in the undirected skeleton (weighted as in
/// [`stress_layout`]). Takes one shortest-path search per node.
pub fn layout_stress(
    g: &Graph,
    layout: &Layout,
    weighted: bool,
    weights_as: WeightKind,
    threads: usize,
) -> Result<f64> {
    let n = g.node_count();
    if layout.x.len() != n || layout.y.len() != n {
        return Err(Error::LengthMismatch {
            what: "layout",
            expected: n,
            found: layout.x.len().min(layout.y.len()),
        });
    }
    let (parts, _) = parts(g, weighted, weights_as)?;
    let mut total = 0.0;
    for p in &parts {
        let c = p.len();
        let t = parallel::resolve_threads(threads, c / 16);
        let mut sums = vec![0.0; c];
        par_fill(&mut sums, t, 16, |i, out| {
            let mut ws = Sssp::new(c);
            let mut d = vec![0.0; c];
            ws.run(p, i, &mut d);
            let vi = p.nodes[i] as usize;
            let mut s = 0.0;
            for (&vj, &dj) in p.nodes[i + 1..].iter().zip(&d[i + 1..]) {
                let vj = vj as usize;
                let e = ((layout.x[vi] - layout.x[vj]).powi(2)
                    + (layout.y[vi] - layout.y[vj]).powi(2))
                .sqrt();
                s += (e - dj) * (e - dj) / (dj * dj);
            }
            *out = s;
        });
        total += sums.iter().sum::<f64>();
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, edges: &[(i64, i64)]) -> Graph {
        let s: Vec<i64> = edges.iter().map(|e| e.0).collect();
        let t: Vec<i64> = edges.iter().map(|e| e.1).collect();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn a_path_is_drawn_straight_with_unit_edges() {
        let g = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        let l = stress_layout(&g, &StressOptions::default()).unwrap();
        for i in 0..4 {
            let e = ((l.x[i + 1] - l.x[i]).powi(2) + (l.y[i + 1] - l.y[i]).powi(2)).sqrt();
            assert!((e - 1.0).abs() < 1e-3, "{e}");
        }
        assert!(layout_stress(&g, &l, false, WeightKind::Distance, 1).unwrap() < 1e-5);
        // Principal axis horizontal, node 0 on the left.
        assert!(l.x[0] < l.x[4]);
        assert!(l.y.iter().all(|y| y.abs() < 1e-2));
    }

    #[test]
    fn a_square_is_a_square() {
        let g = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        let opts = StressOptions {
            tol: 1e-12,
            ..Default::default()
        };
        let l = stress_layout(&g, &opts).unwrap();
        let d = |a: usize, b: usize| ((l.x[a] - l.x[b]).powi(2) + (l.y[a] - l.y[b]).powi(2)).sqrt();
        // The stress optimum for C4 (targets 1 on sides, 2 on diagonals) is a
        // square whose side s minimizes 4 (s - 1)^2 + 2 (s sqrt2 - 2)^2 / 4,
        // so s = (8 + 2 sqrt2) / 10.
        let s = (8.0 + 2.0 * 2f64.sqrt()) / 10.0;
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 0)] {
            assert!((d(a, b) - s).abs() < 1e-3, "{} vs {s}", d(a, b));
        }
        assert!((d(0, 2) - s * 2f64.sqrt()).abs() < 1e-3);
    }

    #[test]
    fn sparse_stress_is_close_to_full_stress_on_a_grid() {
        let mut e = Vec::new();
        for r in 0..6i64 {
            for c in 0..6i64 {
                let v = r * 6 + c;
                if c < 5 {
                    e.push((v, v + 1));
                }
                if r < 5 {
                    e.push((v, v + 6));
                }
            }
        }
        let g = graph(36, &e);
        let full = stress_layout(&g, &StressOptions::default()).unwrap();
        let sf = layout_stress(&g, &full, false, WeightKind::Distance, 1).unwrap();
        // The grid drawn as a square grid has stress 29.10 (computed
        // independently); stress majorization does much better.
        let mut grid = Layout {
            x: vec![0.0; 36],
            y: vec![0.0; 36],
        };
        for v in 0..36 {
            grid.x[v] = (v % 6) as f64;
            grid.y[v] = (v / 6) as f64;
        }
        let sg = layout_stress(&g, &grid, false, WeightKind::Distance, 1).unwrap();
        assert!((sg - 29.104).abs() < 1e-3, "{sg}");
        assert!(sf <= 0.5 * sg, "{sf} vs {sg}");
        // With every node a pivot, each region is a single node, every
        // weight is 1 / d^2, and sparse stress is the full model.
        let sparse = |k: usize| {
            let opts = StressOptions {
                algorithm: StressAlgorithm::Sparse,
                pivots: k,
                ..Default::default()
            };
            let l = stress_layout(&g, &opts).unwrap();
            layout_stress(&g, &l, false, WeightKind::Distance, 1).unwrap()
        };
        let s36 = sparse(36);
        assert!((s36 - sf).abs() <= 1e-3 * sf, "{s36} vs {sf}");
        // A third of the nodes as pivots is a coarse approximation on so
        // small a graph, but stays close.
        let s12 = sparse(12);
        assert!(s12 <= 1.3 * sf, "{s12} vs {sf}");
    }

    #[test]
    fn weights_scale_the_drawing() {
        let s = [0, 1, 2];
        let t = [1, 2, 0];
        let g = Graph::from_edges(3, &s, &t, Some(&[2.0, 2.0, 2.0]), false).unwrap();
        let opts = StressOptions {
            weighted: true,
            ..Default::default()
        };
        let l = stress_layout(&g, &opts).unwrap();
        let e = ((l.x[0] - l.x[1]).powi(2) + (l.y[0] - l.y[1]).powi(2)).sqrt();
        assert!((e - 2.0).abs() < 1e-6, "{e}");
        let strength = StressOptions {
            weighted: true,
            weights_as: WeightKind::Strength,
            ..Default::default()
        };
        let l = stress_layout(&g, &strength).unwrap();
        let e = ((l.x[0] - l.x[1]).powi(2) + (l.y[0] - l.y[1]).powi(2)).sqrt();
        assert!((e - 0.5).abs() < 1e-6, "{e}");
    }
}
