//! Structural equivalence and CONCOR.
//!
//! Two nodes are structurally equivalent when they have the same ties to
//! and from every other node (Lorrain and White 1971,
//! `lorrain1971structural`). In practice nodes are compared by how far
//! apart their *tie profiles* are.
//!
//! # Tie profiles
//!
//! Let a_ik be the tie from i to k: 1 or 0 in a binary network (parallel
//! edges count once), or the total weight of the edges from i to k in a
//! weighted network. The profile of node i is its row and its column of
//! the adjacency matrix, (a_ik, a_ki) for every k. Nodes i and j are
//! compared over k other than i and j (Wasserman and Faust 1994,
//! `wasserman1994social`, sec. 9.2), so neither the diagonal nor the pair's
//! own cells a_ij and a_ji enter the comparison:
//!
//! - **Euclidean distance** (Burt 1976, `burt1976positions`):
//!   d_ij = sqrt(sum_{k != i, j} [(a_ik - a_jk)^2 + (a_ki - a_kj)^2]).
//! - **Hamming distance**: the same sum with absolute differences, which in
//!   a binary network counts the ties that differ.
//! - **Correlation** (Wasserman and Faust 1994, sec. 9.2): the Pearson
//!   correlation of the 2(n - 2) values (a_ik, a_ki) and (a_jk, a_kj),
//!   k != i, j, with means taken over those values.
//!
//! These are exactly the conventions of `sna::sedist()` with its defaults
//! (`diag = FALSE`): sna replaces the diagonal with missing values, stacks
//! each node's row and column, and compares profiles over the entries
//! present in both, which drops k = i and k = j. Following sna, the
//! diagonal of the result is 0 for the distances, and for correlation it
//! is 1, or 0 for a node whose profile is constant. A correlation that is
//! undefined because a profile has no variance over the compared entries
//! (or because there are none, as when n <= 2) is reported as 0, as sna
//! does.
//!
//! An undirected network is the directed network with every edge
//! reciprocated, so each profile holds every tie twice: distances are
//! those of the symmetric adjacency matrix, as in sna, and correlation is
//! unaffected. Self-loops are ignored. Weights may be any finite values.
//!
//! The pairwise sums run over the union of the two nodes' neighbor lists,
//! so all pairs cost O(nm) rather than O(n^3), and each term is evaluated
//! directly rather than by expanding squares, so identical profiles give
//! exactly zero distance. For binary profiles of dense networks the same
//! sums are counted with word-wide Boolean operations on bit rows
//! (O(n^3 / 64), exact, and identical to the sparse sums); every measure
//! is symmetric in the two nodes, so each pair is evaluated once.
//! The iterated correlations of CONCOR are tiled products over centred
//! columns, each sum still taken in order (O(b^3) for a block of b nodes).
//!
//! # CONCOR
//!
//! CONCOR (Breiger, Boorman, and Arabie 1975, `breiger1975algorithm`;
//! Wasserman and Faust 1994, sec. 9.4) repeatedly correlates the columns of
//! a correlation matrix. Starting from the matrix of profile correlations
//! above, the iteration M_{t+1} = cor(M_t) drives every entry towards +1 or
//! -1, which splits the nodes into two blocks by sign. Each block is then
//! split again, starting from the submatrix of the initial correlations
//! for its members, up to the requested depth, giving at most 2^depth
//! blocks. Conventions:
//!
//! - The initial correlations are the structural-equivalence correlations
//!   defined above (so the diagonal and the pair's own cells are excluded,
//!   as in sna and the concorR package).
//! - Nodes whose profile is constant (for example isolates, ignoring
//!   self-loops) cannot be correlated with any other node. They are set
//!   aside before the first split and placed in blocks of their own, one
//!   block per distinct constant value; such nodes are structurally
//!   equivalent to one another.
//! - Iteration stops when every entry has absolute value at least
//!   1 - `tolerance` (default 1e-7), when an iteration changes no entry by
//!   more than [`CONCOR_STALL`] (a fixed point that is not +1 or -1), or
//!   after `max_iter` iterations (default 100). The last two count as not
//!   converged; the split then uses the signs reached. Later correlations use all entries of each column, with the mean
//!   taken over the block. A column with no variance has correlation 0 with
//!   every other column.
//! - The split is by the sign of the converged correlations with the
//!   block's lowest-numbered node: entries greater than zero go with that
//!   node, entries of zero or less form the other block. A block that does
//!   not divide (all entries positive) is kept whole, and blocks of one
//!   node are not split.
//! - CONCOR can be numerically ill-conditioned. Highly symmetric networks
//!   (a directed cycle, for example) may have fixed points whose entries
//!   are not +1 or -1; stopping there, rather than iterating until rounding
//!   error grows, keeps the result reproducible. A node whose correlation
//!   with the reference node decays towards zero is assigned by rounding
//!   error once it is small enough. [`Concor`] counts splits that did not
//!   converge and splits in which, after the first iteration, some node's
//!   correlation with the reference node fell below [`CONCOR_NEAR_ZERO`] in
//!   absolute value; the latter means the partition may differ between
//!   implementations or platforms and should not be relied on.
//! - Blocks are numbered from 1 in the order of a depth-first walk of the
//!   splits, with the side holding the lower-numbered node first, followed
//!   by the blocks of constant-profile nodes in order of their first
//!   member.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::linalg::simd;
use crate::matrix::Dense;
use crate::parallel;

/// How to compare tie profiles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquivalenceMethod {
    /// Euclidean distance (Burt 1976).
    Euclidean,
    /// Pearson correlation (Wasserman and Faust 1994).
    Correlation,
    /// Hamming (city-block) distance.
    Hamming,
}

impl EquivalenceMethod {
    /// Decodes the codes used by the C interface: 0 = Euclidean,
    /// 1 = correlation, 2 = Hamming.
    pub fn from_code(code: i32) -> Result<EquivalenceMethod> {
        match code {
            0 => Ok(EquivalenceMethod::Euclidean),
            1 => Ok(EquivalenceMethod::Correlation),
            2 => Ok(EquivalenceMethod::Hamming),
            _ => Err(Error::InvalidArgument(format!(
                "unknown structural equivalence method code {code}"
            ))),
        }
    }
}

/// Options for [`structural_equivalence`].
#[derive(Debug, Clone, Copy)]
pub struct EquivalenceOptions {
    /// How to compare profiles.
    pub method: EquivalenceMethod,
    /// Use total edge weights instead of binary ties.
    pub weighted: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for EquivalenceOptions {
    fn default() -> Self {
        EquivalenceOptions {
            method: EquivalenceMethod::Euclidean,
            weighted: false,
            threads: 0,
        }
    }
}

/// The n-by-n matrix of profile distances or correlations.
///
/// See the module documentation for conventions.
pub fn structural_equivalence(g: &Graph, opts: EquivalenceOptions) -> Result<Dense> {
    let p = Profiles::new(g, opts.weighted)?;
    Ok(p.compare(opts.method, opts.threads))
}

/// Correlations with the reference node smaller than this in absolute
/// value during the iteration mark a split as numerically fragile.
pub const CONCOR_NEAR_ZERO: f64 = 1e-6;

/// An iteration that changes no entry by more than this has reached a fixed
/// point and stops.
pub const CONCOR_STALL: f64 = 1e-12;

/// Options for [`concor`].
#[derive(Debug, Clone, Copy)]
pub struct ConcorOptions {
    /// Number of rounds of splitting (at least 1).
    pub depth: usize,
    /// Use total edge weights instead of binary ties.
    pub weighted: bool,
    /// Iteration stops when every |entry| >= 1 - tolerance.
    pub tolerance: f64,
    /// Maximum number of iterations per split.
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for ConcorOptions {
    fn default() -> Self {
        ConcorOptions {
            depth: 1,
            weighted: false,
            tolerance: 1e-7,
            max_iter: 100,
            threads: 0,
        }
    }
}

/// Result of [`concor`].
#[derive(Debug, Clone, PartialEq)]
pub struct Concor {
    /// Block of each node, numbered from 1.
    pub membership: Vec<i64>,
    /// Number of splits that reached `max_iter` without converging.
    pub unconverged: usize,
    /// Number of splits in which a correlation with the reference node came
    /// within [`CONCOR_NEAR_ZERO`] of zero, so rounding may have decided
    /// the side of some node.
    pub near_zero: usize,
}

/// CONCOR block partition.
///
/// See the module documentation for conventions.
pub fn concor(g: &Graph, opts: ConcorOptions) -> Result<Concor> {
    if opts.depth == 0 {
        return Err(Error::InvalidArgument(
            "CONCOR depth must be at least 1".into(),
        ));
    }
    if !(opts.tolerance > 0.0 && opts.tolerance < 1.0) {
        return Err(Error::InvalidArgument(
            "CONCOR tolerance must be between 0 and 1".into(),
        ));
    }
    let p = Profiles::new(g, opts.weighted)?;
    let n = p.n;

    // Set aside constant profiles, grouped by their value.
    let mut constant: Vec<(f64, Vec<usize>)> = Vec::new();
    let mut active = Vec::new();
    for v in 0..n {
        match p.constant_value(v) {
            Some(c) => match constant.iter_mut().find(|(x, _)| *x == c) {
                Some((_, members)) => members.push(v),
                None => constant.push((c, vec![v])),
            },
            None => active.push(v),
        }
    }

    let init = p.compare(EquivalenceMethod::Correlation, opts.threads);
    let mut blocks = vec![active];
    let mut unconverged = 0;
    let mut near_zero = 0;
    for _ in 0..opts.depth {
        let mut next = Vec::with_capacity(blocks.len() * 2);
        for block in blocks {
            if block.len() < 2 {
                if !block.is_empty() {
                    next.push(block);
                }
                continue;
            }
            let it = iterate(&init, &block, opts);
            unconverged += !it.converged as usize;
            near_zero += it.near_zero as usize;
            let m = it.matrix;
            let mut first = Vec::new();
            let mut second = Vec::new();
            for (c, &v) in block.iter().enumerate() {
                if m.get(0, c) > 0.0 {
                    first.push(v);
                } else {
                    second.push(v);
                }
            }
            for part in [first, second] {
                if !part.is_empty() {
                    next.push(part);
                }
            }
        }
        blocks = next;
    }

    let mut membership = vec![0i64; n];
    let mut label = 0;
    for block in blocks.iter().chain(constant.iter().map(|(_, b)| b)) {
        label += 1;
        for &v in block {
            membership[v] = label;
        }
    }
    Ok(Concor {
        membership,
        unconverged,
        near_zero,
    })
}

struct Iteration {
    matrix: Dense,
    converged: bool,
    near_zero: bool,
}

/// Iterates correlations on the submatrix of `init` for `block`.
fn iterate(init: &Dense, block: &[usize], opts: ConcorOptions) -> Iteration {
    let b = block.len();
    let mut m = Dense::zeros(b, b);
    for (r, &i) in block.iter().enumerate() {
        for (c, &j) in block.iter().enumerate() {
            m.set(r, c, init.get(i, j));
        }
    }
    let cutoff = 1.0 - opts.tolerance;
    let converged = |m: &Dense| m.as_slice().iter().all(|x| x.abs() >= cutoff);
    let near = |m: &Dense| m.row(0).iter().any(|x| x.abs() < CONCOR_NEAR_ZERO);
    // Exact zeros in the starting correlations are data, not rounding, so
    // only iterates are checked for entries near zero.
    let mut near_zero = false;
    for _ in 0..opts.max_iter {
        if converged(&m) {
            break;
        }
        let next = {
            let _t = crate::linalg::phase::start("concor: column correlations");
            column_correlations(&m, opts.threads)
        };
        let step = next
            .as_slice()
            .iter()
            .zip(m.as_slice())
            .fold(0.0f64, |acc, (a, b)| acc.max((a - b).abs()));
        m = next;
        near_zero |= near(&m);
        if step <= CONCOR_STALL {
            break;
        }
    }
    Iteration {
        converged: converged(&m),
        near_zero,
        matrix: m,
    }
}

/// Pearson correlations between the columns of a square matrix, computed
/// with centred columns (two passes) for accuracy. Columns without variance
/// have correlation 0 with every column, including themselves.
fn column_correlations(m: &Dense, threads: usize) -> Dense {
    let b = m.rows();
    let threads = parallel::threads_for(
        threads,
        b,
        (b as u64).saturating_mul(b as u64).saturating_mul(b as u64) / 2,
    );
    // Centred columns, stored as rows for contiguous access, GROUP columns
    // per task. Each task reads the rows of `m` in order, so every column's
    // mean and sum of squares add their terms in row order, as one sum per
    // column would.
    let mut centred_data = vec![0.0; b * b];
    let mut norm = vec![0.0; b];
    if b > 0 {
        let start: f64 = std::iter::empty::<f64>().sum();
        let mut norms = vec![0.0; b];
        let mut both: Vec<(&mut [f64], &mut [f64])> = centred_data
            .chunks_mut(GROUP * b)
            .zip(norms.chunks_mut(GROUP))
            .collect();
        parallel::fill_chunks(
            &mut both,
            1,
            threads,
            Vec::new,
            |mean: &mut Vec<f64>, g, part| {
                let (rows, ns) = (&mut *part[0].0, &mut *part[0].1);
                let c0 = g * GROUP;
                let w = ns.len();
                mean.clear();
                mean.resize(w, start);
                for r in 0..b {
                    for (acc, &x) in mean.iter_mut().zip(&m.row(r)[c0..c0 + w]) {
                        *acc += x;
                    }
                }
                mean.iter_mut().for_each(|x| *x /= b as f64);
                ns.iter_mut().for_each(|x| *x = 0.0);
                for r in 0..b {
                    for (c, &x) in m.row(r)[c0..c0 + w].iter().enumerate() {
                        let y = x - mean[c];
                        rows[c * b + r] = y;
                        ns[c] += y * y;
                    }
                }
                ns.iter_mut().for_each(|x| *x = x.sqrt());
            },
        );
        drop(both);
        norm = norms;
    }
    let centred = Dense::from_row_major(b, b, centred_data).expect("square");
    // Correlation is symmetric and each dot product is the same sequence
    // of exactly rounded products whichever column comes first, so only
    // j >= i is computed and mirrored. The products are tiled: each block
    // of RB rows of `centred` is packed so that the RB values for one k are
    // adjacent, and each pass over k updates RB x CB running sums, one per
    // pair, from one packed load and CB scalars. Every sum still adds its
    // terms in order of k, starting from the value the standard library's
    // sum starts from, so the result is that of one sum per pair; the
    // compiler can keep the RB sums of a column in vector registers since
    // they are independent.
    //
    // Work is split into square tasks of GROUP rows by GROUP columns (on or
    // above the diagonal), so that the packed rows and the columns a task
    // reads stay in cache while it computes GROUP^2 sums from them; a pass
    // over every column per block of rows would instead be limited by
    // memory bandwidth once several threads run. Which task computes a
    // pair does not change its sum.
    let zero: f64 = std::iter::empty::<f64>().sum();
    let blocks = (b + RB - 1) / RB;
    let mut packed = vec![0.0; blocks * b * RB];
    for (blk, pack) in packed.chunks_exact_mut(b * RB).enumerate() {
        let i0 = blk * RB;
        let i1 = (i0 + RB).min(b);
        for r in 0..RB {
            let row = centred.row((i0 + r).min(i1 - 1));
            for (k, &x) in row.iter().enumerate() {
                pack[k * RB + r] = x;
            }
        }
    }
    let groups = (b + GROUP - 1) / GROUP;
    let tasks: Vec<(usize, usize)> = (0..groups)
        .flat_map(|gi| (gi..groups).map(move |gj| (gi, gj)))
        .collect();
    let tiles = parallel::map_items(
        tasks.len(),
        threads,
        || (),
        |_, t| {
            let (gi, gj) = tasks[t];
            let (r0, r1) = (gi * GROUP, ((gi + 1) * GROUP).min(b));
            let (c0, c1) = (gj * GROUP, ((gj + 1) * GROUP).min(b));
            let mut tile = vec![0.0; GROUP * GROUP];
            for i0 in (r0..r1).step_by(RB) {
                let i1 = (i0 + RB).min(r1);
                let pack = &packed[(i0 / RB) * b * RB..(i0 / RB + 1) * b * RB];
                let mut j = c0.max(i0);
                while j < c1 {
                    let cols: [&[f64]; CB] =
                        std::array::from_fn(|c| centred.row((j + c).min(b - 1)));
                    let acc = simd::dispatch!(correlation_tile(pack, &cols, b, zero));
                    for (r, i) in (i0..i1).enumerate() {
                        for (c, sums) in acc.iter().enumerate() {
                            let jj = j + c;
                            if jj < c1 && jj >= i && norm[i] > 0.0 && norm[jj] > 0.0 {
                                tile[(i - r0) * GROUP + (jj - c0)] =
                                    (sums[r] / (norm[i] * norm[jj])).clamp(-1.0, 1.0);
                            }
                        }
                    }
                    j += CB;
                }
            }
            tile
        },
    );
    // Each row of the result reads its entries from the tiles on and above
    // the diagonal, the row's own for j >= i and the mirrored ones below.
    let mut index = vec![usize::MAX; groups * groups];
    for (t, &(gi, gj)) in tasks.iter().enumerate() {
        index[gi * groups + gj] = t;
    }
    let mut out = vec![0.0; b * b];
    parallel::fill_chunks(
        &mut out,
        b,
        threads,
        || (),
        |_, i, row| {
            let gi = i / GROUP;
            for (j, x) in row.iter_mut().enumerate() {
                let gj = j / GROUP;
                *x = if j >= i {
                    tiles[index[gi * groups + gj]][(i - gi * GROUP) * GROUP + (j - gj * GROUP)]
                } else {
                    tiles[index[gj * groups + gi]][(j - gj * GROUP) * GROUP + (i - gi * GROUP)]
                };
            }
        },
    );
    Dense::from_row_major(b, b, out).expect("square")
}

/// Rows and columns per task in [`column_correlations`] (a multiple of
/// `RB` and `CB`).
const GROUP: usize = 64;
/// Rows per tile in [`column_correlations`].
const RB: usize = 8;
/// Columns per tile in [`column_correlations`].
const CB: usize = 4;

/// The RB x CB dot products of the packed rows (`pack[k * RB + r]`) with
/// the columns `cols`, each summed in order of k from `zero`.
#[inline(always)]
fn correlation_tile(pack: &[f64], cols: &[&[f64]; CB], b: usize, zero: f64) -> [[f64; RB]; CB] {
    let mut acc = [[zero; RB]; CB];
    let (y0, y1, y2, y3) = (&cols[0][..b], &cols[1][..b], &cols[2][..b], &cols[3][..b]);
    for (k, p) in pack[..b * RB].chunks_exact(RB).enumerate() {
        let y = [y0[k], y1[k], y2[k], y3[k]];
        for (sums, &yc) in acc.iter_mut().zip(&y) {
            for (s, &x) in sums.iter_mut().zip(p) {
                *s += x * yc;
            }
        }
    }
    acc
}

/// Distinct ties without self-loops: out[v] holds (k, a_vk), inc[v] holds
/// (k, a_kv), each sorted by k with nonzero values only.
struct Profiles {
    n: usize,
    out: Vec<Vec<(u32, f64)>>,
    inc: Vec<Vec<(u32, f64)>>,
}

impl Profiles {
    fn new(g: &Graph, weighted: bool) -> Result<Profiles> {
        if weighted && !g.is_weighted() {
            return Err(Error::MissingWeights);
        }
        let n = g.node_count();
        let build = |csr: &crate::graph::Csr| -> Vec<Vec<(u32, f64)>> {
            (0..n)
                .map(|v| {
                    let mut row: Vec<(u32, f64)> = Vec::new();
                    let w = csr.weights(v);
                    for (k, &t) in csr.neighbors(v).iter().enumerate() {
                        if t as usize == v {
                            continue;
                        }
                        let x = match (weighted, w) {
                            (true, Some(w)) => w[k],
                            _ => 1.0,
                        };
                        match row.last_mut() {
                            Some(last) if last.0 == t => {
                                if weighted {
                                    last.1 += x;
                                }
                            }
                            _ => row.push((t, x)),
                        }
                    }
                    row.retain(|&(_, x)| x != 0.0);
                    row
                })
                .collect()
        };
        let out = build(g.out());
        let inc = if g.is_directed() {
            build(g.inc())
        } else {
            out.clone()
        };
        Ok(Profiles { n, out, inc })
    }

    /// The common value of every entry of v's profile, if it is constant.
    fn constant_value(&self, v: usize) -> Option<f64> {
        let others = self.n - 1;
        let row = &self.out[v];
        let col = &self.inc[v];
        if row.is_empty() && col.is_empty() {
            return Some(0.0);
        }
        if row.len() == others && col.len() == others {
            let c = row[0].1;
            if row.iter().chain(col).all(|&(_, x)| x == c) {
                return Some(c);
            }
        }
        None
    }

    /// The matrix of [`Profiles::pair`] over all pairs. Every method is
    /// symmetric bit for bit (swapping i and j swaps the two arguments of
    /// each term, and every term and the final combination are symmetric
    /// in them), so only i < j is evaluated and mirrored. Binary profiles
    /// of networks dense enough are compared as bit rows (see
    /// [`BitProfiles`]), which gives the same values.
    fn compare(&self, method: EquivalenceMethod, threads: usize) -> Dense {
        let n = self.n;
        let bits = BitProfiles::new(self);
        let rows = parallel::for_each_item(n, threads, Vec::new, |out, i| {
            let mut row = vec![0.0; n - i];
            row[0] = self.pair(method, i, i);
            for (d, slot) in row.iter_mut().enumerate().skip(1) {
                let j = i + d;
                *slot = match &bits {
                    Some(b) => b.pair(method, i, j),
                    None => self.pair(method, i, j),
                };
            }
            out.push((i, row));
        });
        let mut m = vec![0.0; n * n];
        for part in rows {
            for (i, row) in part {
                for (d, x) in row.into_iter().enumerate() {
                    let j = i + d;
                    m[i * n + j] = x;
                    m[j * n + i] = x;
                }
            }
        }
        Dense::from_row_major(n, n, m).expect("square")
    }

    fn pair(&self, method: EquivalenceMethod, i: usize, j: usize) -> f64 {
        if i == j {
            return match method {
                EquivalenceMethod::Correlation => self.self_correlation(i),
                _ => 0.0,
            };
        }
        let skip = [i as u32, j as u32];
        match method {
            EquivalenceMethod::Euclidean => {
                let mut s = 0.0;
                let mut f = |x: f64, y: f64| s += (x - y) * (x - y);
                merge(&self.out[i], &self.out[j], skip, &mut f);
                merge(&self.inc[i], &self.inc[j], skip, &mut f);
                s.sqrt()
            }
            EquivalenceMethod::Hamming => {
                let mut s = 0.0;
                let mut f = |x: f64, y: f64| s += (x - y).abs();
                merge(&self.out[i], &self.out[j], skip, &mut f);
                merge(&self.inc[i], &self.inc[j], skip, &mut f);
                s
            }
            EquivalenceMethod::Correlation => {
                let mut st = Moments::default();
                let mut f = |x: f64, y: f64| st.add(x, y);
                merge(&self.out[i], &self.out[j], skip, &mut f);
                merge(&self.inc[i], &self.inc[j], skip, &mut f);
                // Entries absent from both lists are zeros in both profiles.
                st.pearson(2 * (self.n - 2))
            }
        }
    }

    /// sna's diagonal: the profile correlated with itself over k != i.
    fn self_correlation(&self, i: usize) -> f64 {
        if self.n < 2 {
            return 0.0;
        }
        let count = 2 * (self.n - 1);
        let mut st = Moments::default();
        for &(_, x) in self.out[i].iter().chain(&self.inc[i]) {
            st.add(x, x);
        }
        st.pearson(count)
    }
}

/// Binary profiles as rows of bits: for each node, its out-row and then
/// its in-row, `words` 64-bit words each.
///
/// When every profile entry is 0 or 1, each sum that [`Profiles::pair`]
/// accumulates is a count, and counts are exact in floating point (below
/// 2^53) whatever the order of the additions. They are popcounts: over
/// k other than i and j, the sum of x y is the number of common bits of the
/// rows of i and j (the bit of a node in its own rows is never set, so i
/// and j drop out by themselves), the sum of x (and of x^2) is the degree
/// of i less its bits at j, and the sum of (x - y)^2 and of |x - y| is the
/// number of differing bits less those at i and j. The values are those of
/// the sparse merge exactly (Kepner and Gilbert 2011, `kepner2011graph`,
/// on graph algorithms as Boolean matrix products).
struct BitProfiles {
    n: usize,
    words: usize,
    bits: Vec<u64>,
    degree: Vec<f64>,
}

impl BitProfiles {
    /// Bit rows for binary profiles, or `None` when some entry is not 1 or
    /// when the network is too sparse for words to beat merging lists.
    fn new(p: &Profiles) -> Option<BitProfiles> {
        let n = p.n;
        let words = (n + 63) / 64;
        let entries: usize = p.out.iter().chain(&p.inc).map(|r| r.len()).sum();
        let binary = p.out.iter().chain(&p.inc).flatten().all(|&(_, x)| x == 1.0);
        // A merge step per entry of the two profiles against two words per
        // 64 columns, with a memory cap of 256 MiB.
        let worthwhile = n >= 2 && 2 * entries >= 2 * words * n && n * words <= 1 << 24;
        if !binary || !worthwhile {
            return None;
        }
        let mut bits = vec![0u64; 2 * words * n];
        for v in 0..n {
            for (half, row) in [&p.out[v], &p.inc[v]].into_iter().enumerate() {
                let base = (2 * v + half) * words;
                for &(k, _) in row {
                    bits[base + k as usize / 64] |= 1u64 << (k % 64);
                }
            }
        }
        let degree = (0..n)
            .map(|v| (p.out[v].len() + p.inc[v].len()) as f64)
            .collect();
        Some(BitProfiles {
            n,
            words,
            bits,
            degree,
        })
    }

    fn row(&self, v: usize, half: usize) -> &[u64] {
        let base = (2 * v + half) * self.words;
        &self.bits[base..base + self.words]
    }

    fn bit(&self, v: usize, half: usize, k: usize) -> f64 {
        ((self.row(v, half)[k / 64] >> (k % 64)) & 1) as f64
    }

    /// [`Profiles::pair`] for i != j.
    fn pair(&self, method: EquivalenceMethod, i: usize, j: usize) -> f64 {
        let count = |op: fn(u64, u64) -> u64| -> f64 {
            let mut c = 0u64;
            for half in 0..2 {
                for (a, b) in self.row(i, half).iter().zip(self.row(j, half)) {
                    c += op(*a, *b).count_ones() as u64;
                }
            }
            c as f64
        };
        let at_j = self.bit(i, 0, j) + self.bit(i, 1, j);
        let at_i = self.bit(j, 0, i) + self.bit(j, 1, i);
        match method {
            EquivalenceMethod::Euclidean => (count(|a, b| a ^ b) - at_j - at_i).sqrt(),
            EquivalenceMethod::Hamming => count(|a, b| a ^ b) - at_j - at_i,
            EquivalenceMethod::Correlation => {
                let sx = self.degree[i] - at_j;
                let sy = self.degree[j] - at_i;
                let st = Moments {
                    sx,
                    sy,
                    sxx: sx,
                    syy: sy,
                    sxy: count(|a, b| a & b),
                };
                st.pearson(2 * (self.n - 2))
            }
        }
    }
}

/// Calls `f(x, y)` for every k in the union of two sorted sparse rows,
/// skipping the indices in `skip`; a missing entry is 0.
fn merge(a: &[(u32, f64)], b: &[(u32, f64)], skip: [u32; 2], f: &mut impl FnMut(f64, f64)) {
    let (mut p, mut q) = (0, 0);
    while p < a.len() || q < b.len() {
        let ka = a.get(p).map_or(u32::MAX, |e| e.0);
        let kb = b.get(q).map_or(u32::MAX, |e| e.0);
        let k = ka.min(kb);
        let x = if ka == k {
            p += 1;
            a[p - 1].1
        } else {
            0.0
        };
        let y = if kb == k {
            q += 1;
            b[q - 1].1
        } else {
            0.0
        };
        if k != skip[0] && k != skip[1] {
            f(x, y);
        }
    }
}

/// Sums for a Pearson correlation over nonzero entries; the remaining
/// entries are zeros in both profiles.
#[derive(Default)]
struct Moments {
    sx: f64,
    sy: f64,
    sxx: f64,
    syy: f64,
    sxy: f64,
}

impl Moments {
    fn add(&mut self, x: f64, y: f64) {
        self.sx += x;
        self.sy += y;
        self.sxx += x * x;
        self.syy += y * y;
        self.sxy += x * y;
    }

    /// Correlation over `count` entries; 0 when undefined. With integer
    /// ties every intermediate quantity is an exact integer.
    fn pearson(&self, count: usize) -> f64 {
        if count < 2 {
            return 0.0;
        }
        let c = count as f64;
        let vx = c * self.sxx - self.sx * self.sx;
        let vy = c * self.syy - self.sy * self.sy;
        if vx <= 0.0 || vy <= 0.0 {
            return 0.0;
        }
        let cov = c * self.sxy - self.sx * self.sy;
        (cov / (vx.sqrt() * vy.sqrt())).clamp(-1.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dir(n: usize, e: &[(i64, i64)]) -> Graph {
        let s: Vec<i64> = e.iter().map(|x| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        Graph::from_edges(n, &s, &t, None, true).unwrap()
    }

    /// Direct evaluation of the definitions on a dense matrix.
    fn dense_reference(a: &[Vec<f64>], method: EquivalenceMethod) -> Vec<Vec<f64>> {
        let n = a.len();
        let mut out = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let ks: Vec<usize> = (0..n).filter(|&k| k != i && k != j).collect();
                let x: Vec<f64> = ks.iter().flat_map(|&k| [a[i][k], a[k][i]]).collect();
                let y: Vec<f64> = ks.iter().flat_map(|&k| [a[j][k], a[k][j]]).collect();
                out[i][j] = match method {
                    EquivalenceMethod::Euclidean => x
                        .iter()
                        .zip(&y)
                        .map(|(p, q)| (p - q) * (p - q))
                        .sum::<f64>()
                        .sqrt(),
                    EquivalenceMethod::Hamming => {
                        x.iter().zip(&y).map(|(p, q)| (p - q).abs()).sum()
                    }
                    EquivalenceMethod::Correlation => {
                        let len = x.len() as f64;
                        if len < 2.0 {
                            0.0
                        } else {
                            let mx = x.iter().sum::<f64>() / len;
                            let my = y.iter().sum::<f64>() / len;
                            let sxy: f64 = x.iter().zip(&y).map(|(p, q)| (p - mx) * (q - my)).sum();
                            let sxx: f64 = x.iter().map(|p| (p - mx) * (p - mx)).sum();
                            let syy: f64 = y.iter().map(|q| (q - my) * (q - my)).sum();
                            if sxx == 0.0 || syy == 0.0 {
                                0.0
                            } else {
                                sxy / (sxx * syy).sqrt()
                            }
                        }
                    }
                };
            }
        }
        out
    }

    #[test]
    fn matches_dense_definition() {
        let e = [
            (0, 1),
            (1, 2),
            (2, 0),
            (3, 0),
            (3, 1),
            (4, 4),
            (1, 3),
            (1, 2),
        ];
        let g = dir(5, &e);
        let mut a = vec![vec![0.0; 5]; 5];
        for &(s, t) in &e {
            if s != t {
                a[s as usize][t as usize] = 1.0;
            }
        }
        for method in [
            EquivalenceMethod::Euclidean,
            EquivalenceMethod::Hamming,
            EquivalenceMethod::Correlation,
        ] {
            let m = structural_equivalence(
                &g,
                EquivalenceOptions {
                    method,
                    threads: 2,
                    ..Default::default()
                },
            )
            .unwrap();
            let r = dense_reference(&a, method);
            for (i, ri) in r.iter().enumerate() {
                for (j, &rij) in ri.iter().enumerate() {
                    let want = if i == j && method == EquivalenceMethod::Correlation {
                        // sna's diagonal: 1, or 0 for a constant profile.
                        if i == 4 {
                            0.0
                        } else {
                            1.0
                        }
                    } else {
                        rij
                    };
                    assert!((m.get(i, j) - want).abs() < 1e-12, "{method:?} {i} {j}");
                }
            }
        }
    }

    #[test]
    fn equivalent_nodes_have_zero_distance() {
        // 1 and 2 both receive from 0 and send to 3.
        let g = dir(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        let m = structural_equivalence(&g, EquivalenceOptions::default()).unwrap();
        assert_eq!(m.get(1, 2), 0.0);
        let c = structural_equivalence(
            &g,
            EquivalenceOptions {
                method: EquivalenceMethod::Correlation,
                ..Default::default()
            },
        )
        .unwrap();
        assert!((c.get(1, 2) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn bit_rows_match_the_sparse_merge() {
        let mut r = crate::rng::Stream::new(11, 0);
        for (n, p, directed) in [
            (3usize, 0.9, true),
            (40, 0.5, true),
            (64, 0.3, false),
            (65, 0.6, true),
            (130, 0.4, true),
            (130, 0.2, false),
        ] {
            let (mut s, mut t) = (Vec::new(), Vec::new());
            for i in 0..n as i64 {
                for j in 0..n as i64 {
                    if (directed || i < j) && r.uniform() < p {
                        s.push(i);
                        t.push(j);
                    }
                }
            }
            let g = Graph::from_edges(n, &s, &t, None, directed).unwrap();
            let prof = Profiles::new(&g, false).unwrap();
            let bits = BitProfiles::new(&prof).expect("dense enough for bit rows");
            for method in [
                EquivalenceMethod::Euclidean,
                EquivalenceMethod::Hamming,
                EquivalenceMethod::Correlation,
            ] {
                for i in 0..n {
                    for j in 0..n {
                        if i != j {
                            let (a, b) = (bits.pair(method, i, j), prof.pair(method, i, j));
                            assert_eq!(a.to_bits(), b.to_bits(), "{n} {method:?} {i} {j}");
                            let c = prof.pair(method, j, i);
                            assert_eq!(b.to_bits(), c.to_bits(), "symmetry {n} {i} {j}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn column_correlations_match_one_sum_per_pair() {
        let mut r = crate::rng::Stream::new(5, 0);
        for b in [1usize, 2, 3, 4, 5, 7, 9, 30] {
            let mut m = Dense::zeros(b, b);
            for i in 0..b {
                for j in 0..b {
                    // Some constant columns, which correlate 0 with all.
                    let x = if j % 4 == 3 { 0.5 } else { r.uniform() };
                    m.set(i, j, x);
                }
            }
            for threads in [1, 2, 3] {
                let fast = column_correlations(&m, threads);
                for i in 0..b {
                    for j in 0..b {
                        let col = |c: usize| -> (Vec<f64>, f64) {
                            let mean = (0..b).map(|r| m.get(r, c)).sum::<f64>() / b as f64;
                            let x: Vec<f64> = (0..b).map(|r| m.get(r, c) - mean).collect();
                            let ss: f64 = x.iter().map(|v| v * v).sum();
                            (x, ss.sqrt())
                        };
                        let ((xi, ni), (xj, nj)) = (col(i), col(j));
                        let want = if ni > 0.0 && nj > 0.0 {
                            let dot: f64 = xi.iter().zip(&xj).map(|(x, y)| x * y).sum();
                            (dot / (ni * nj)).clamp(-1.0, 1.0)
                        } else {
                            0.0
                        };
                        assert_eq!(fast.get(i, j).to_bits(), want.to_bits(), "{b} {i} {j}");
                    }
                }
            }
        }
    }

    #[test]
    fn concor_separates_two_positions() {
        // Complete bipartite K_{3,3}: the two sides are the positions.
        let mut e = Vec::new();
        for a in 0..3 {
            for b in 3..6 {
                e.push((a, b));
            }
        }
        let s: Vec<i64> = e.iter().map(|x: &(i64, i64)| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        let g = Graph::from_edges(7, &s, &t, None, false).unwrap();
        let r = concor(&g, ConcorOptions::default()).unwrap();
        // Node 6 is an isolate and forms its own block.
        assert_eq!(r.membership, vec![1, 1, 1, 2, 2, 2, 3]);
        assert_eq!((r.unconverged, r.near_zero), (0, 0));
        assert!(concor(
            &g,
            ConcorOptions {
                depth: 0,
                ..Default::default()
            }
        )
        .is_err());
    }
}
