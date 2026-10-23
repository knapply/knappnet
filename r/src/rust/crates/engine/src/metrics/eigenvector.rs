//! Eigenvector centrality.
//!
//! A node's score is proportional to the sum of its neighbors' scores: the
//! entries of the principal eigenvector of the adjacency matrix (Bonacich
//! 1972, `bonacich1972factoring`; see also Bonacich 1987,
//! `bonacich1987power`). For a nonnegative matrix whose graph is strongly
//! connected, the Perron-Frobenius theorem guarantees that this eigenvector
//! is unique up to scale and positive.
//!
//! # Algorithms
//!
//! Symmetric matrices (undirected graphs, or mode `All`) use the Lanczos
//! method (Lanczos 1950, `lanczos1950iteration`) with thick restarts (Wu and
//! Simon 2000, `wu2000thick`; the symmetric case of Stewart's Krylov-Schur
//! method, `stewart2002krylov`), full reorthogonalization by classical
//! Gram-Schmidt with the DGKS correction (Daniel, Gragg, Kaufman, and
//! Stewart 1976, `daniel1976reorthogonalization`), and the small projected
//! eigenproblems solved by the cyclic Jacobi method (Jacobi 1846,
//! `jacobi1846leichtes`; Golub and Van Loan 2013, `golub2013matrix`,
//! section 8.5). Each connected component is solved separately, starting
//! from the all-ones vector. A Krylov method started from a vector b finds
//! the component of b in the leading eigenspace, so when that eigenvalue is
//! repeated the answer is the same as the limit of power iteration from b.
//! The Lanczos method needs far fewer matrix-vector products than power
//! iteration (its error falls with the square root of the relative gap
//! between the two largest eigenvalues rather than with the gap itself),
//! and it is not slowed down by bipartite or periodic graphs, since it seeks
//! the algebraically largest eigenvalue rather than the largest in modulus.
//!
//! Directed matrices use power iteration (von Mises and Pollaczek-Geiringer
//! 1929, `mises1929praktische`) on A' + I, where A' is the transposed (mode
//! `In`) or plain (mode `Out`) adjacency matrix. The shift by the identity
//! leaves the eigenvectors unchanged and makes the dominant eigenvalue
//! strictly largest in modulus even for periodic graphs, where plain power
//! iteration on A oscillates; NetworkX's `eigenvector_centrality` uses the
//! same shift. The iteration starts from the all-ones vector and rescales to
//! a largest entry of 1 after every step. Convergence is geometric with
//! ratio |lambda_2 + 1| / (lambda_1 + 1). Power iteration is kept for
//! directed graphs because its limit is defined on every nonnegative matrix,
//! including reducible ones, whose leading eigenvalue can be defective.
//!
//! Both solvers stop under the same rule: when one more step of power
//! iteration on A' + I would change no score, with the largest score scaled
//! to 1, by more than `tol` (default 1e-12). For the Lanczos method this is
//! checked on the residual: max_i |(A'x)_i - lambda x_i| <= tol (lambda + 1)
//! max_i |x_i|, which is that change to first order. If this does not happen
//! within `max_iter` matrix-vector products (default 10,000; per component
//! for symmetric matrices) an error is returned rather than an inaccurate
//! answer. Matrix-vector products and the vector operations of the Lanczos
//! method run in parallel on large graphs; sums are taken over fixed blocks
//! of nodes in a fixed order, so the result does not depend on the thread
//! count.
//!
//! For symmetric matrices, a component is skipped (scores 0) when Collatz's
//! bound shows that its largest eigenvalue is below one already found: the
//! largest eigenvalue of a nonnegative matrix is at most the largest row
//! sum, and more generally lies between the minimum and maximum of
//! (A'x)_i / x_i for any positive x (Collatz 1942,
//! `collatz1942einschliessungssatz`). Components are visited in decreasing
//! order of a lower bound on their largest eigenvalue (the largest of the
//! mean row sum and the square root of the largest sum of squares in a row),
//! so most components are never solved.
//!
//! # Conventions
//!
//! - Directed graphs: with mode `In` (igraph with `directed = TRUE`,
//!   NetworkX) a node's score is the sum of the scores of nodes with arcs
//!   *to* it, the left eigenvector of A. With `Out` (sna's `evcent`) it sums
//!   the nodes it points to, the right eigenvector. `All` ignores direction,
//!   so two opposite arcs become two parallel edges (igraph's R default,
//!   `directed = FALSE`). Undirected graphs ignore the mode.
//! - Parallel edges add up: their count, or their weights.
//! - Self-loops ([`LoopWeight`]): igraph puts twice the number of loops on
//!   the diagonal of an undirected graph (each loop is two edge endpoints)
//!   and NetworkX once; in directed graphs both count a loop once. sna
//!   ignores the diagonal.
//! - Weights are used as given (larger means a stronger tie) and must be
//!   nonnegative.
//! - Scaling ([`EigenScale`]): `Max` makes the largest score 1 (igraph);
//!   `Unit` gives the vector Euclidean length 1 (NetworkX, sna, UCINET).
//!   Scores are nonnegative; rounding can leave entries of the order of
//!   1e-16 below zero in the Lanczos solution, and these are set to zero.
//! - Graphs without cycles (including graphs without edges): the only
//!   eigenvalue is 0. Following igraph, nodes whose score feeds no other
//!   node's score (nodes without outgoing arcs for mode `In`, without
//!   incoming arcs for `Out`, and every node in a graph without edges) score
//!   1 before scaling, every other node 0, and the eigenvalue is 0.
//! - Disconnected graphs (not strongly connected, for directed graphs) have
//!   no unique positive eigenvector, and [`Eigenvector::connected`] is false
//!   so callers can warn. igraph returns scores without checking; NetworkX's
//!   `eigenvector_centrality_numpy` refuses. Here:
//!   - Symmetric matrices: components whose largest eigenvalue is below the
//!     overall largest score 0. If several components attain it (to a
//!     relative 1e-9), each contributes its unit Perron vector u times
//!     sum(u), which is the limit of power iteration from equal scores.
//!   - Directed matrices: the result is the limit of the iteration, which
//!     is zero outside the strongly connected components that carry the
//!     principal eigenvalue and those they reach (`In`) or that reach them
//!     (`Out`). When several components share the principal eigenvalue and
//!     one reaches another, the matrix is defective, the iteration converges
//!     only sublinearly, and an error is returned.

use crate::error::{Error, Result};
use crate::graph::{Csr, Graph};
use crate::metrics::closeness::oriented;
use crate::metrics::degree::DegreeMode;
use crate::metrics::topography::{components, ComponentKind};
use crate::parallel;

/// How eigenvector scores are scaled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EigenScale {
    /// Largest score 1 (igraph).
    #[default]
    Max,
    /// Euclidean norm 1 (NetworkX, sna, UCINET).
    Unit,
}

impl EigenScale {
    /// Decodes the C interface codes: 0 = max, 1 = unit.
    pub fn from_code(code: i32) -> Result<EigenScale> {
        match code {
            0 => Ok(EigenScale::Max),
            1 => Ok(EigenScale::Unit),
            _ => Err(Error::InvalidArgument(format!(
                "unknown eigenvector scale code {code}"
            ))),
        }
    }
}

/// What a self-loop contributes to the diagonal of an undirected graph's
/// adjacency matrix. Directed loops always count once unless ignored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LoopWeight {
    /// Twice its count or weight (igraph).
    #[default]
    Twice,
    /// Once (NetworkX).
    Once,
    /// Not at all (sna, UCINET).
    Ignore,
}

impl LoopWeight {
    /// Decodes the C interface codes: 0 = twice, 1 = once, 2 = ignore.
    pub fn from_code(code: i32) -> Result<LoopWeight> {
        match code {
            0 => Ok(LoopWeight::Twice),
            1 => Ok(LoopWeight::Once),
            2 => Ok(LoopWeight::Ignore),
            _ => Err(Error::InvalidArgument(format!(
                "unknown self-loop weight code {code}"
            ))),
        }
    }

    /// The diagonal entry a loop of weight 1 contributes.
    pub(crate) fn factor(self, directed: bool) -> f64 {
        match (directed, self) {
            (_, LoopWeight::Ignore) => 0.0,
            (false, LoopWeight::Twice) => 2.0,
            _ => 1.0,
        }
    }
}

/// Default convergence tolerance (largest change in any score, with the
/// largest score scaled to 1).
pub const DEFAULT_TOL: f64 = 1e-12;
/// Default iteration cap.
pub const DEFAULT_MAX_ITER: usize = 10_000;

/// Options for [`eigenvector_centrality`].
#[derive(Debug, Clone, Copy)]
pub struct EigenvectorOptions {
    /// Direction for directed graphs (`In` is the usual definition).
    pub mode: DegreeMode,
    /// Use edge weights as given (they must be nonnegative).
    pub weighted: bool,
    /// Scaling of the result.
    pub scale: EigenScale,
    /// Self-loops in undirected graphs.
    pub loops: LoopWeight,
    /// Convergence tolerance.
    pub tol: f64,
    /// Cap on matrix-vector products (per component, for symmetric
    /// matrices).
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for EigenvectorOptions {
    fn default() -> Self {
        EigenvectorOptions {
            mode: DegreeMode::In,
            weighted: false,
            scale: EigenScale::Max,
            loops: LoopWeight::Twice,
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
            threads: 0,
        }
    }
}

/// Result of [`eigenvector_centrality`].
#[derive(Debug, Clone, PartialEq)]
pub struct Eigenvector {
    /// Scores in node order.
    pub values: Vec<f64>,
    /// The principal eigenvalue (0 for graphs without cycles).
    pub eigenvalue: f64,
    /// Matrix-vector products used, summed over the components solved (0
    /// when no iteration was needed).
    pub iterations: usize,
    /// Whether the graph is connected (strongly, if directed and the mode
    /// is not `All`). When false the scores depend on more than the
    /// network; see the module documentation.
    pub connected: bool,
}

/// Arcs whose number exceeds this use several threads per step.
pub(crate) const PARALLEL_ARCS: usize = 1 << 16;

/// Eigenvector centrality. See the module documentation for conventions.
pub fn eigenvector_centrality(g: &Graph, opts: EigenvectorOptions) -> Result<Eigenvector> {
    if opts.tol.is_nan() || opts.tol <= 0.0 || opts.max_iter == 0 {
        return Err(Error::InvalidArgument(
            "the tolerance and the iteration cap must be positive".into(),
        ));
    }
    let n = g.node_count();
    if n == 0 {
        return Ok(Eigenvector {
            values: Vec::new(),
            eigenvalue: f64::NAN,
            iterations: 0,
            connected: true,
        });
    }
    if opts.weighted {
        check_nonnegative(g, "eigenvector centrality")?;
    }
    // Rows of the iteration matrix: x_i <- sum over row(i) of w x_j.
    // `In` sums over arcs arriving at i; `Out` over arcs leaving i.
    let mode = if g.is_directed() {
        opts.mode
    } else {
        DegreeMode::All
    };
    let h = oriented(
        g,
        if mode == DegreeMode::All {
            mode
        } else {
            DegreeMode::Out
        },
    )?;
    let h = h.as_ref();
    let rows = if h.is_directed() && mode == DegreeMode::In {
        h.inc()
    } else {
        h.out()
    };
    let op = SparseOp::from_rows(
        rows,
        n,
        opts.weighted,
        opts.loops.factor(h.is_directed()),
        opts.threads,
    );

    let kind = if h.is_directed() {
        ComponentKind::Strong
    } else {
        ComponentKind::Weak
    };
    let connected = components(h, kind).count() == 1;

    // Acyclic case: every eigenvalue is 0. A graph has a cycle with
    // positive weight exactly when some strong component has a positive
    // arc inside it.
    let comp = components(h, ComponentKind::Strong).membership;
    let cyclic = (0..n).any(|v| op.row(v).any(|(j, w)| comp[j] == comp[v] && w > 0.0));
    if !cyclic {
        // Nodes with no positive entries in their column of the matrix
        // feed no other score: igraph's sinks, for mode `In`.
        let mut sink = vec![1.0; n];
        for v in 0..n {
            for (j, w) in op.row(v) {
                if w > 0.0 {
                    sink[j] = 0.0;
                }
            }
        }
        scale(&mut sink, opts.scale);
        return Ok(Eigenvector {
            values: sink,
            eigenvalue: 0.0,
            iterations: 0,
            connected,
        });
    }

    let (mut x, eigenvalue, iterations) = if h.is_directed() {
        power_iteration(&op, opts)?
    } else {
        symmetric_components(h, &op, opts)?
    };
    scale(&mut x, opts.scale);
    Ok(Eigenvector {
        values: x,
        eigenvalue,
        iterations,
        connected,
    })
}

/// Refuses negative edge weights.
pub(crate) fn check_nonnegative(g: &Graph, what: &str) -> Result<()> {
    let w = g.edge_weights().ok_or(Error::MissingWeights)?;
    if let Some(edge) = w.iter().position(|&x| x < 0.0) {
        return Err(Error::InvalidArgument(format!(
            "edge {edge} has weight {}; {what} requires nonnegative weights",
            w[edge]
        )));
    }
    Ok(())
}

/// Power iteration on A' + I for directed matrices, iterated as one block.
/// Returns the scores (largest 1 before the final combination), the
/// eigenvalue, and the iteration count.
fn power_iteration(op: &SparseOp, opts: EigenvectorOptions) -> Result<(Vec<f64>, f64, usize)> {
    let n = op.dim();
    let mut x = vec![1.0; n];
    let mut y = vec![0.0; n];
    for iter in 1..=opts.max_iter {
        op.apply(&x, &mut y);
        let mut top = 0.0f64;
        for v in 0..n {
            y[v] += x[v];
            top = top.max(y[v]);
        }
        let mut change = 0.0f64;
        for v in 0..n {
            y[v] /= top;
            change = change.max((y[v] - x[v]).abs());
        }
        std::mem::swap(&mut x, &mut y);
        if change <= opts.tol {
            // Eigenvalue from sum(A'x) / sum(x), then the scores rescaled
            // by sum(x) / sum(x^2) (a positive constant).
            op.apply(&x, &mut y);
            let (mut s0, mut s1, mut s2) = (0.0f64, 0.0f64, 0.0f64);
            for v in 0..n {
                s0 += y[v];
                s1 += x[v];
                s2 += x[v] * x[v];
            }
            let lambda = if s1 > 0.0 { s0 / s1 } else { 0.0 };
            if s2 > 0.0 {
                let w = s1 / s2;
                for xv in x.iter_mut() {
                    *xv *= w;
                }
            }
            return Ok((x, lambda, iter));
        }
    }
    Err(Error::InvalidArgument(format!(
        "eigenvector centrality did not converge to tolerance {} within {} iterations. \
         Power iteration converges slowly when the largest eigenvalue is close to the \
         next one, and not at all to a single answer when several strongly connected \
         components share it; a larger iteration cap may help in the first case",
        opts.tol, opts.max_iter
    )))
}

/// Relative difference within which two components' largest eigenvalues
/// are treated as equal.
const TIE: f64 = 1e-9;

/// The Lanczos solution of a symmetric matrix, one connected component at a
/// time. Returns unscaled scores, the eigenvalue, and the number of
/// matrix-vector products.
fn symmetric_components(
    h: &Graph,
    op: &SparseOp,
    opts: EigenvectorOptions,
) -> Result<(Vec<f64>, f64, usize)> {
    let n = op.dim();
    let comps = components(h, ComponentKind::Weak);
    let k = comps.count();
    // Bounds on each component's largest eigenvalue: the largest row sum
    // above (Collatz), and below the mean row sum (the Rayleigh quotient of
    // the all-ones vector) and the square root of the largest sum of
    // squares in a row (lambda^2 >= (A^2)_vv for symmetric A).
    let mut upper = vec![0.0f64; k];
    let mut sum = vec![0.0f64; k];
    let mut square = vec![0.0f64; k];
    for v in 0..n {
        let c = comps.membership[v] as usize;
        let (s, q) = op.row_sums(v);
        upper[c] = upper[c].max(s);
        sum[c] += s;
        square[c] = square[c].max(q);
    }
    let lower: Vec<f64> = (0..k)
        .map(|c| (sum[c] / comps.sizes[c] as f64).max(square[c].sqrt()))
        .collect();
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&a, &b| lower[b].total_cmp(&lower[a]).then(a.cmp(&b)));

    // Nodes of each component, in node order.
    let mut start = vec![0usize; k + 1];
    for &c in &comps.membership {
        start[c as usize + 1] += 1;
    }
    for c in 0..k {
        start[c + 1] += start[c];
    }
    let mut members = vec![0u32; n];
    let mut cursor = start.clone();
    for (v, &c) in comps.membership.iter().enumerate() {
        members[cursor[c as usize]] = v as u32;
        cursor[c as usize] += 1;
    }

    let mut best = 0.0f64;
    let mut solved: Vec<(usize, f64, Vec<f64>)> = Vec::new();
    let mut matvecs = 0usize;
    let mut local = vec![0u32; n];
    for &c in &order {
        if upper[c] <= 0.0 || upper[c] < best * (1.0 - TIE) {
            continue;
        }
        let nodes = &members[start[c]..start[c + 1]];
        let sub;
        let a = if k == 1 {
            op
        } else {
            sub = op.restrict(nodes, &mut local);
            &sub
        };
        let ones = vec![1.0; nodes.len()];
        let res = symmetric_leading(
            a.dim(),
            &|x: &[f64], y: &mut [f64]| a.apply(x, y),
            &ones,
            &Stop {
                tol: opts.tol,
                shift: 1.0,
                max_matvecs: opts.max_iter,
            },
            a.threads,
        );
        let Some(res) = res else {
            return Err(Error::InvalidArgument(format!(
                "eigenvector centrality did not converge to tolerance {} within {} \
                 matrix-vector products; a larger iteration cap may help",
                opts.tol, opts.max_iter
            )));
        };
        matvecs += res.matvecs;
        best = best.max(res.value);
        solved.push((c, res.value, res.vector));
    }
    // Each component that attains the largest eigenvalue contributes its
    // unit Perron vector u times sum(u), the limit of power iteration from
    // equal scores; the others contribute nothing.
    let mut x = vec![0.0; n];
    for (c, value, u) in solved {
        if value < best * (1.0 - TIE) {
            continue;
        }
        let total: f64 = u.iter().sum();
        for (&v, &uv) in members[start[c]..start[c + 1]].iter().zip(&u) {
            x[v as usize] = uv * total;
        }
    }
    Ok((x, best, matvecs))
}

fn scale(x: &mut [f64], how: EigenScale) {
    let d = match how {
        EigenScale::Max => x.iter().copied().fold(0.0f64, f64::max),
        EigenScale::Unit => x.iter().map(|v| v * v).sum::<f64>().sqrt(),
    };
    if d > 0.0 {
        for v in x {
            *v /= d;
        }
    }
}

// ---- Sparse matrices ------------------------------------------------------

/// A sparse matrix in row form, with entries ready for multiplication: arc
/// weights (or 1) with self-loops multiplied by their diagonal factor.
/// Shared with the spectral measures.
pub(crate) struct SparseOp {
    offsets: Vec<usize>,
    cols: Vec<u32>,
    /// `None` when every entry is 1.
    vals: Option<Vec<f64>>,
    /// Threads used by [`SparseOp::apply`] and by the solvers working on
    /// this matrix (1 for small matrices).
    pub(crate) threads: usize,
}

impl SparseOp {
    /// The `dim` by `dim` matrix whose row v lists `rows.neighbors(v)`,
    /// with weights if `weighted` and self-loops multiplied by
    /// `loop_factor`.
    pub(crate) fn from_rows(
        rows: &Csr,
        dim: usize,
        weighted: bool,
        loop_factor: f64,
        threads: usize,
    ) -> SparseOp {
        let mut offsets = Vec::with_capacity(dim + 1);
        offsets.push(0);
        let mut cols = Vec::with_capacity(rows.arc_count());
        let w = if weighted { rows.all_weights() } else { None };
        let ones = w.is_none() && loop_factor == 1.0;
        let mut vals = if ones {
            None
        } else {
            Some(Vec::with_capacity(rows.arc_count()))
        };
        for v in 0..dim {
            for (k, &j) in rows.neighbors(v).iter().enumerate() {
                cols.push(j);
                if let Some(vals) = vals.as_mut() {
                    let x = w.map_or(1.0, |w| w[rows.arc_range(v).start + k]);
                    vals.push(if j as usize == v { x * loop_factor } else { x });
                }
            }
            offsets.push(cols.len());
        }
        let threads = if cols.len() >= PARALLEL_ARCS {
            parallel::resolve_threads(threads, dim)
        } else {
            1
        };
        SparseOp {
            offsets,
            cols,
            vals,
            threads,
        }
    }

    /// Number of rows (and columns).
    pub(crate) fn dim(&self) -> usize {
        self.offsets.len() - 1
    }

    /// The entries of row v as (column, value).
    pub(crate) fn row(&self, v: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let r = self.offsets[v]..self.offsets[v + 1];
        let vals = self.vals.as_ref().map(|x| &x[r.clone()]);
        self.cols[r]
            .iter()
            .enumerate()
            .map(move |(k, &j)| (j as usize, vals.map_or(1.0, |x| x[k])))
    }

    /// Sum of row v, and the sum of squares of its entries with parallel
    /// entries (which are adjacent, since rows are sorted) merged first.
    fn row_sums(&self, v: usize) -> (f64, f64) {
        let (mut s, mut q) = (0.0f64, 0.0f64);
        let mut run: Option<(usize, f64)> = None;
        for (j, w) in self.row(v) {
            s += w;
            run = match run {
                Some((c, acc)) if c == j => Some((c, acc + w)),
                Some((_, acc)) => {
                    q += acc * acc;
                    Some((j, w))
                }
                None => Some((j, w)),
            };
        }
        if let Some((_, acc)) = run {
            q += acc * acc;
        }
        (s, q)
    }

    /// The submatrix on `nodes` (rows and columns), which must be closed
    /// under the matrix's nonzero pattern (a union of components).
    /// `local` is scratch space of the full dimension.
    fn restrict(&self, nodes: &[u32], local: &mut [u32]) -> SparseOp {
        for (i, &v) in nodes.iter().enumerate() {
            local[v as usize] = i as u32;
        }
        self.submatrix(nodes, local)
    }

    /// The rows `rows`, with column j renumbered `col_index[j]` and dropped
    /// when that is `u32::MAX`.
    pub(crate) fn submatrix(&self, rows: &[u32], col_index: &[u32]) -> SparseOp {
        let mut offsets = Vec::with_capacity(rows.len() + 1);
        offsets.push(0);
        let mut cols = Vec::new();
        let mut vals = self.vals.as_ref().map(|_| Vec::new());
        for &v in rows {
            let r = self.offsets[v as usize]..self.offsets[v as usize + 1];
            for k in r {
                let c = col_index[self.cols[k] as usize];
                if c == u32::MAX {
                    continue;
                }
                cols.push(c);
                if let (Some(dst), Some(src)) = (vals.as_mut(), self.vals.as_ref()) {
                    dst.push(src[k]);
                }
            }
            offsets.push(cols.len());
        }
        let threads = if cols.len() >= PARALLEL_ARCS {
            self.threads.max(1)
        } else {
            1
        };
        SparseOp {
            offsets,
            cols,
            vals,
            threads,
        }
    }

    /// The transpose, for a matrix with `ncol` columns.
    pub(crate) fn transpose(&self, ncol: usize) -> SparseOp {
        let mut offsets = vec![0usize; ncol + 1];
        for &j in &self.cols {
            offsets[j as usize + 1] += 1;
        }
        for c in 0..ncol {
            offsets[c + 1] += offsets[c];
        }
        let mut cursor = offsets.clone();
        let mut cols = vec![0u32; self.cols.len()];
        let mut vals = self.vals.as_ref().map(|_| vec![0.0; self.cols.len()]);
        for v in 0..self.dim() {
            for k in self.offsets[v]..self.offsets[v + 1] {
                let j = self.cols[k] as usize;
                let slot = cursor[j];
                cols[slot] = v as u32;
                if let (Some(dst), Some(src)) = (vals.as_mut(), self.vals.as_ref()) {
                    dst[slot] = src[k];
                }
                cursor[j] += 1;
            }
        }
        let threads = self.threads;
        SparseOp {
            offsets,
            cols,
            vals,
            threads,
        }
    }

    /// y = A x.
    pub(crate) fn apply(&self, x: &[f64], y: &mut [f64]) {
        let n = self.dim();
        let row = |v: usize| -> f64 {
            let r = self.offsets[v]..self.offsets[v + 1];
            let mut acc = 0.0;
            match &self.vals {
                None => {
                    for &j in &self.cols[r] {
                        acc += x[j as usize];
                    }
                }
                Some(vals) => {
                    for (&j, &w) in self.cols[r.clone()].iter().zip(&vals[r]) {
                        acc += w * x[j as usize];
                    }
                }
            }
            acc
        };
        if self.threads <= 1 {
            for (v, yv) in y.iter_mut().enumerate() {
                *yv = row(v);
            }
            return;
        }
        // Contiguous blocks, one value per row, so the result does not
        // depend on the thread count.
        let threads = self.threads;
        let block = (n + threads - 1) / threads;
        let parts = parallel::for_each_item(
            threads,
            threads,
            Vec::new,
            |out: &mut Vec<(usize, Vec<f64>)>, b| {
                let lo = (b * block).min(n);
                let hi = ((b + 1) * block).min(n);
                out.push((lo, (lo..hi).map(row).collect()));
            },
        );
        for part in parts {
            for (lo, vals) in part {
                y[lo..lo + vals.len()].copy_from_slice(&vals);
            }
        }
    }
}

// ---- The Lanczos method ---------------------------------------------------

/// When the Lanczos method stops: max_i |r_i| <= tol (lambda + shift)
/// max_i |x_i| for the residual r = Ax - lambda x, or after `max_matvecs`
/// matrix-vector products without success.
pub(crate) struct Stop {
    pub(crate) tol: f64,
    pub(crate) shift: f64,
    pub(crate) max_matvecs: usize,
}

/// A converged eigenpair.
pub(crate) struct Leading {
    /// The eigenvalue (Rayleigh quotient of `vector`).
    pub(crate) value: f64,
    /// Unit eigenvector, signed so its entries sum to a nonnegative number,
    /// with negative entries (rounding error, for nonnegative matrices) set
    /// to zero.
    pub(crate) vector: Vec<f64>,
    /// Matrix-vector products used.
    pub(crate) matvecs: usize,
}

/// Largest Krylov basis kept in memory.
const MAX_BASIS: usize = 20;
/// Ritz vectors kept at a thick restart.
const KEEP: usize = 8;
/// Basis memory is capped at about this many doubles; larger problems use a
/// smaller basis.
const MAX_BASIS_DOUBLES: usize = 1 << 27;
/// Rows per block for sums over vectors (fixed, so results do not depend
/// on the thread count).
const BLOCK: usize = 4096;
/// Lanczos steps between convergence checks within a cycle.
const CHECK_EVERY: usize = 4;

/// The algebraically largest eigenvalue of the symmetric matrix applied by
/// `apply`, and its eigenvector, by the thick-restart Lanczos method
/// started from `start`. Returns `None` if `stop.max_matvecs` products do
/// not suffice.
pub(crate) fn symmetric_leading(
    dim: usize,
    apply: &dyn Fn(&[f64], &mut [f64]),
    start: &[f64],
    stop: &Stop,
    threads: usize,
) -> Option<Leading> {
    debug_assert!(dim > 0 && start.len() == dim);
    let m = dim
        .min(MAX_BASIS)
        .min((MAX_BASIS_DOUBLES / dim).max(8))
        .max(1);
    let keep = KEEP.min(m.saturating_sub(2)).max(1);
    let blocks = Blocks::new(dim, threads);
    let mut basis = Basis {
        dim,
        data: vec![0.0; dim * m],
    };
    let mut hm = vec![0.0f64; m * m];
    let mut v = start.to_vec();
    let norm = blocks.dot(&v, &v).sqrt();
    if norm == 0.0 || !norm.is_finite() {
        return None;
    }
    v.iter_mut().for_each(|x| *x /= norm);
    basis.set_col(0, &v);
    let mut w = vec![0.0; dim];
    let mut y = vec![0.0; dim];
    let mut matvecs = 0usize;
    let mut j0 = 0usize;
    // The residual estimate must fall below this before the true residual
    // is checked; it is tightened after a failed check.
    let mut threshold = f64::INFINITY;
    loop {
        let mut j = j0;
        loop {
            if matvecs >= stop.max_matvecs {
                return None;
            }
            apply(&v, &mut w);
            matvecs += 1;
            let (coef, before, beta) = basis.orthogonalize(j + 1, &mut w, &blocks);
            for (i, &c) in coef.iter().enumerate() {
                hm[i * m + j] = c;
                hm[j * m + i] = c;
            }
            let size = j + 1;
            let breakdown = beta <= 1e-12 * before || beta.is_nan();
            let end = breakdown || size == m;
            if end || (size - j0) % CHECK_EVERY == 0 {
                let (theta, s) = symmetric_eigen(&hm, m, size);
                let top = argmax(&theta);
                let estimate = if breakdown {
                    0.0
                } else {
                    beta * s[(size - 1) * size + top].abs()
                };
                let limit = stop.tol * (theta[top].max(0.0) + stop.shift);
                if estimate <= limit.min(threshold) {
                    // Candidate: check the true residual.
                    let coeffs: Vec<f64> = (0..size).map(|l| s[l * size + top]).collect();
                    let x = basis.combine(size, &coeffs, &blocks);
                    if matvecs >= stop.max_matvecs {
                        return None;
                    }
                    apply(&x, &mut y);
                    matvecs += 1;
                    let xx = blocks.dot(&x, &x);
                    let value = blocks.dot(&x, &y) / xx;
                    let (mut rinf, mut xinf) = (0.0f64, 0.0f64);
                    for i in 0..dim {
                        rinf = rinf.max((y[i] - value * x[i]).abs());
                        xinf = xinf.max(x[i].abs());
                    }
                    let need = stop.tol * (value.max(0.0) + stop.shift) * xinf;
                    if rinf <= need {
                        return Some(finish_vector(x, value, matvecs, &blocks));
                    }
                    threshold = if rinf > 0.0 {
                        estimate.min(limit) * (need / rinf) * 0.5
                    } else {
                        estimate.min(limit) * 0.5
                    };
                    if breakdown {
                        // The Krylov space is invariant but the residual is
                        // not small enough: restart from the candidate.
                        let nx = xx.sqrt();
                        v = x.iter().map(|xi| xi / nx).collect();
                        basis.set_col(0, &v);
                        hm.iter_mut().for_each(|e| *e = 0.0);
                        j0 = 0;
                        break;
                    }
                }
                if end {
                    if breakdown {
                        // An invariant subspace that holds the answer passes
                        // the check above; anything else is rounding noise.
                        return None;
                    }
                    // Thick restart: keep the leading Ritz vectors, then
                    // continue from the residual direction.
                    let mut idx: Vec<usize> = (0..size).collect();
                    idx.sort_by(|&a, &b| theta[b].total_cmp(&theta[a]).then(a.cmp(&b)));
                    let p = keep.min(size - 1);
                    let mut svec = vec![0.0; size * p];
                    for q in 0..p {
                        for l in 0..size {
                            svec[l * p + q] = s[l * size + idx[q]];
                        }
                    }
                    basis.rotate(size, &svec, p, &blocks);
                    hm.iter_mut().for_each(|e| *e = 0.0);
                    for q in 0..p {
                        hm[q * m + q] = theta[idx[q]];
                        let c = beta * s[(size - 1) * size + idx[q]];
                        hm[q * m + p] = c;
                        hm[p * m + q] = c;
                    }
                    v = w.iter().map(|x| x / beta).collect();
                    basis.set_col(p, &v);
                    j0 = p;
                    break;
                }
            }
            v = w.iter().map(|x| x / beta).collect();
            basis.set_col(j + 1, &v);
            hm[(j + 1) * m + j] = beta;
            hm[j * m + j + 1] = beta;
            j += 1;
        }
    }
}

fn finish_vector(mut x: Vec<f64>, value: f64, matvecs: usize, blocks: &Blocks) -> Leading {
    let total: f64 = blocks.sum(&x);
    if total < 0.0 {
        x.iter_mut().for_each(|v| *v = -*v);
    }
    for v in x.iter_mut() {
        if *v < 0.0 {
            *v = 0.0;
        }
    }
    let norm = blocks.dot(&x, &x).sqrt();
    if norm > 0.0 {
        x.iter_mut().for_each(|v| *v /= norm);
    }
    Leading {
        value,
        vector: x,
        matvecs,
    }
}

fn argmax(x: &[f64]) -> usize {
    let mut best = 0;
    for (i, &v) in x.iter().enumerate() {
        if v > x[best] {
            best = i;
        }
    }
    best
}

/// Fixed blocks of rows for reductions, processed in parallel when there
/// are several and more than one thread.
struct Blocks {
    dim: usize,
    count: usize,
    threads: usize,
}

impl Blocks {
    fn new(dim: usize, threads: usize) -> Blocks {
        let count = (dim + BLOCK - 1) / BLOCK;
        Blocks {
            dim,
            count,
            threads: threads.max(1).min(count.max(1)),
        }
    }

    fn range(&self, b: usize) -> std::ops::Range<usize> {
        (b * BLOCK)..((b + 1) * BLOCK).min(self.dim)
    }

    /// `f` on every block, results in block order.
    fn map<T: Send>(&self, f: impl Fn(std::ops::Range<usize>) -> T + Sync) -> Vec<T> {
        if self.threads <= 1 {
            return (0..self.count).map(|b| f(self.range(b))).collect();
        }
        let parts = parallel::for_each_item(
            self.count,
            self.threads,
            Vec::new,
            |out: &mut Vec<(usize, T)>, b| out.push((b, f(self.range(b)))),
        );
        let mut slots: Vec<Option<T>> = (0..self.count).map(|_| None).collect();
        for part in parts {
            for (b, t) in part {
                slots[b] = Some(t);
            }
        }
        slots
            .into_iter()
            .map(|t| t.expect("every block is computed"))
            .collect()
    }

    fn dot(&self, a: &[f64], b: &[f64]) -> f64 {
        self.map(|r| {
            let mut s = 0.0;
            for i in r {
                s += a[i] * b[i];
            }
            s
        })
        .into_iter()
        .sum()
    }

    fn sum(&self, a: &[f64]) -> f64 {
        self.map(|r| a[r].iter().sum::<f64>()).into_iter().sum()
    }
}

/// Krylov basis vectors, each stored contiguously: entry i of vector l at
/// `data[l * dim + i]`.
struct Basis {
    dim: usize,
    data: Vec<f64>,
}

impl Basis {
    fn col(&self, l: usize) -> &[f64] {
        &self.data[l * self.dim..(l + 1) * self.dim]
    }

    fn set_col(&mut self, l: usize, v: &[f64]) {
        let dim = self.dim;
        self.data[l * dim..(l + 1) * dim].copy_from_slice(v);
    }

    /// Orthogonalizes `w` against the first `cols` basis vectors. In exact
    /// arithmetic the Lanczos vector A v_j has components only along v_j
    /// and v_(j-1), so those are removed first; one pass of classical
    /// Gram-Schmidt against every vector then restores orthogonality lost
    /// to rounding, and is repeated when it cancels most of what is left
    /// (the DGKS criterion). Returns the coefficients and the norm of `w`
    /// before and after.
    fn orthogonalize(&self, cols: usize, w: &mut [f64], blocks: &Blocks) -> (Vec<f64>, f64, f64) {
        let before = blocks.dot(w, w).sqrt();
        let lo = cols.saturating_sub(2);
        let mut coef = vec![0.0; cols];
        let local = self.project(lo..cols, w, blocks);
        let mut mid = self.subtract(lo..cols, &local, w, blocks);
        for (c, x) in coef[lo..].iter_mut().zip(&local) {
            *c += x;
        }
        loop {
            let again = self.project(0..cols, w, blocks);
            let after = self.subtract(0..cols, &again, w, blocks);
            for (c, x) in coef.iter_mut().zip(&again) {
                *c += x;
            }
            if after >= std::f64::consts::FRAC_1_SQRT_2 * mid || after == 0.0 {
                return (coef, before, after);
            }
            if mid == before {
                // Second full pass already done.
                return (coef, before, after);
            }
            mid = before;
        }
    }

    /// V^T w over the vectors in `range`.
    fn project(&self, range: std::ops::Range<usize>, w: &[f64], blocks: &Blocks) -> Vec<f64> {
        let k = range.len();
        let parts = blocks.map(|r| {
            let wr = &w[r.clone()];
            range
                .clone()
                .map(|l| {
                    let v = &self.col(l)[r.clone()];
                    let mut s = 0.0;
                    for (a, b) in v.iter().zip(wr) {
                        s += a * b;
                    }
                    s
                })
                .collect::<Vec<f64>>()
        });
        let mut total = vec![0.0f64; k];
        for p in parts {
            for (t, x) in total.iter_mut().zip(&p) {
                *t += x;
            }
        }
        total
    }

    /// w -= V h over the vectors in `range`; returns the new ||w||.
    fn subtract(
        &self,
        range: std::ops::Range<usize>,
        h: &[f64],
        w: &mut [f64],
        blocks: &Blocks,
    ) -> f64 {
        let src: &[f64] = w;
        let parts = blocks.map(|r| {
            let mut out = src[r.clone()].to_vec();
            for (l, &c) in range.clone().zip(h) {
                let v = &self.col(l)[r.clone()];
                for (o, &x) in out.iter_mut().zip(v) {
                    *o -= c * x;
                }
            }
            let sq: f64 = out.iter().map(|x| x * x).sum();
            (out, sq)
        });
        let mut total = 0.0f64;
        let mut lo = 0;
        for (vals, sq) in parts {
            w[lo..lo + vals.len()].copy_from_slice(&vals);
            lo += vals.len();
            total += sq;
        }
        total.sqrt()
    }

    /// V s for coefficients `s` over the first `cols` vectors.
    fn combine(&self, cols: usize, s: &[f64], blocks: &Blocks) -> Vec<f64> {
        let parts = blocks.map(|r| {
            let mut out = vec![0.0; r.len()];
            for (l, &c) in s.iter().enumerate().take(cols) {
                let v = &self.col(l)[r.clone()];
                for (o, &x) in out.iter_mut().zip(v) {
                    *o += c * x;
                }
            }
            out
        });
        parts.concat()
    }

    /// Replaces the first `p` vectors by V S, for a `cols` by `p` matrix S
    /// stored by rows.
    fn rotate(&mut self, cols: usize, s: &[f64], p: usize, blocks: &Blocks) {
        // Block by block, so each basis vector is read once.
        let parts = blocks.map(|r| {
            let len = r.len();
            let mut out = vec![0.0; p * len];
            for l in 0..cols {
                let v = &self.col(l)[r.clone()];
                for q in 0..p {
                    let c = s[l * p + q];
                    for (o, &x) in out[q * len..(q + 1) * len].iter_mut().zip(v) {
                        *o += c * x;
                    }
                }
            }
            out
        });
        let dim = self.dim;
        for (b, part) in parts.into_iter().enumerate() {
            let r = blocks.range(b);
            let len = r.len();
            for q in 0..p {
                self.data[q * dim + r.start..q * dim + r.end]
                    .copy_from_slice(&part[q * len..(q + 1) * len]);
            }
        }
    }
}

/// Eigenvalues and eigenvectors of the leading `size` by `size` block of the
/// symmetric matrix `a` (row stride `stride`), by the cyclic Jacobi method.
/// Eigenvector q is column q of the returned row-major `size` by `size`
/// matrix.
pub(crate) fn symmetric_eigen(a: &[f64], stride: usize, size: usize) -> (Vec<f64>, Vec<f64>) {
    let n = size;
    let mut b = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            // Symmetrize against rounding in the projected matrix.
            b[i * n + j] = 0.5 * (a[i * stride + j] + a[j * stride + i]);
        }
    }
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _sweep in 0..100 {
        let mut off = 0.0f64;
        let mut diag = 0.0f64;
        for p in 0..n {
            diag += b[p * n + p] * b[p * n + p];
            for q in (p + 1)..n {
                off += b[p * n + q] * b[p * n + q];
            }
        }
        // Stop when the off-diagonal part is below rounding level relative
        // to the whole matrix (Golub and Van Loan 2013, section 8.5.3).
        if off == 0.0 || off <= 1e-32 * (diag + 2.0 * off) {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = b[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let app = b[p * n + p];
                let aqq = b[q * n + q];
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta.abs() > 1e150 {
                    0.5 / theta
                } else {
                    let sign = if theta < 0.0 { -1.0 } else { 1.0 };
                    sign / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    if k == p || k == q {
                        continue;
                    }
                    let akp = b[k * n + p];
                    let akq = b[k * n + q];
                    let np = c * akp - s * akq;
                    let nq = s * akp + c * akq;
                    b[k * n + p] = np;
                    b[p * n + k] = np;
                    b[k * n + q] = nq;
                    b[q * n + k] = nq;
                }
                b[p * n + p] = app - t * apq;
                b[q * n + q] = aqq + t * apq;
                b[p * n + q] = 0.0;
                b[q * n + p] = 0.0;
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let values = (0..n).map(|i| b[i * n + i]).collect();
    (values, v)
}
#[cfg(test)]
mod tests {
    use super::*;

    fn run(g: &Graph, opts: EigenvectorOptions) -> Eigenvector {
        eigenvector_centrality(g, opts).unwrap()
    }

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)
    }

    #[test]
    fn star_closed_form() {
        // Star with k leaves: eigenvalue sqrt(k), hub sqrt(k) times a leaf.
        let k = 4.0f64;
        let g = Graph::from_edges(5, &[0, 0, 0, 0], &[1, 2, 3, 4], None, false).unwrap();
        let e = run(&g, Default::default());
        assert!((e.eigenvalue - k.sqrt()).abs() < 1e-10);
        assert!(close(&e.values, &[1.0, 0.5, 0.5, 0.5, 0.5], 1e-10));
        assert!(e.connected);
        let u = run(
            &g,
            EigenvectorOptions {
                scale: EigenScale::Unit,
                ..Default::default()
            },
        );
        let norm: f64 = u.values.iter().map(|v| v * v).sum();
        assert!((norm - 1.0).abs() < 1e-12);
    }

    #[test]
    fn bipartite_graph_converges() {
        // A path is bipartite; plain power iteration on A would oscillate.
        let g = Graph::from_edges(4, &[0, 1, 2], &[1, 2, 3], None, false).unwrap();
        let e = run(&g, Default::default());
        let phi = (1.0 + 5f64.sqrt()) / 2.0;
        assert!((e.eigenvalue - phi).abs() < 1e-10);
        assert!(close(&e.values, &[1.0 / phi, 1.0, 1.0, 1.0 / phi], 1e-10));
    }

    #[test]
    fn directed_modes() {
        // 0 -> 1, 1 -> 2, 2 -> 0, 2 -> 1: in-scores favor 1.
        let g = Graph::from_edges(3, &[0, 1, 2, 2], &[1, 2, 0, 1], None, true).unwrap();
        let ins = run(&g, Default::default());
        let outs = run(
            &g,
            EigenvectorOptions {
                mode: DegreeMode::Out,
                ..Default::default()
            },
        );
        // Check A^T x = lambda x and A x = lambda x directly.
        let a = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0]];
        for (i, a_i) in a.iter().enumerate() {
            let col: f64 = (0..3).map(|j| a[j][i] * ins.values[j]).sum();
            let row: f64 = (0..3).map(|j| a_i[j] * outs.values[j]).sum();
            assert!((col - ins.eigenvalue * ins.values[i]).abs() < 1e-10);
            assert!((row - outs.eigenvalue * outs.values[i]).abs() < 1e-10);
        }
        assert!(ins.values[1] > ins.values[0]);
    }

    #[test]
    fn loops_and_weights() {
        // Undirected loop at 0 counts twice by default (A_00 = 2).
        let g = Graph::from_edges(2, &[0, 0], &[0, 1], Some(&[1.0, 3.0]), false).unwrap();
        let twice = run(&g, Default::default());
        // Binary A = [[2, 1], [1, 0]]: eigenvalue 1 + sqrt(2).
        assert!((twice.eigenvalue - (1.0 + 2f64.sqrt())).abs() < 1e-10);
        let once = run(
            &g,
            EigenvectorOptions {
                loops: LoopWeight::Once,
                weighted: true,
                ..Default::default()
            },
        );
        // Weighted A = [[1, 3], [3, 0]]: eigenvalue (1 + sqrt(37)) / 2.
        assert!((once.eigenvalue - (1.0 + 37f64.sqrt()) / 2.0).abs() < 1e-10);
        let neg = Graph::from_edges(2, &[0], &[1], Some(&[-1.0]), false).unwrap();
        let o = EigenvectorOptions {
            weighted: true,
            ..Default::default()
        };
        assert!(eigenvector_centrality(&neg, o).is_err());
    }

    #[test]
    fn acyclic_and_empty_graphs() {
        let empty = Graph::from_edges(3, &[], &[], None, false).unwrap();
        let e = run(&empty, Default::default());
        assert_eq!(e.values, vec![1.0, 1.0, 1.0]);
        assert_eq!(e.eigenvalue, 0.0);
        assert!(!e.connected);
        // 0 -> 1 -> 2: with in-ties the sink 2 scores 1.
        let dag = Graph::from_edges(3, &[0, 1], &[1, 2], None, true).unwrap();
        assert_eq!(run(&dag, Default::default()).values, vec![0.0, 0.0, 1.0]);
        let out = EigenvectorOptions {
            mode: DegreeMode::Out,
            ..Default::default()
        };
        assert_eq!(run(&dag, out).values, vec![1.0, 0.0, 0.0]);
    }

    #[test]
    fn disconnected_graph_is_flagged() {
        // Triangle plus an edge: the triangle (eigenvalue 2) dominates.
        let g = Graph::from_edges(5, &[0, 1, 2, 3], &[1, 2, 0, 4], None, false).unwrap();
        let e = run(&g, Default::default());
        assert!(!e.connected);
        assert!((e.eigenvalue - 2.0).abs() < 1e-10);
        assert!(close(&e.values, &[1.0, 1.0, 1.0, 0.0, 0.0], 1e-9));
    }

    #[test]
    fn thread_count_does_not_change_results() {
        let n = 3000usize;
        let mut s = Vec::new();
        let mut t = Vec::new();
        for i in 0..n {
            for d in [1, 7, 31, 101, 997, 1999] {
                s.push(i as i64);
                t.push(((i * d + 3) % n) as i64);
            }
        }
        let g = Graph::from_edges(n, &s, &t, None, true).unwrap();
        let one = run(
            &g,
            EigenvectorOptions {
                threads: 1,
                ..Default::default()
            },
        );
        let four = run(
            &g,
            EigenvectorOptions {
                threads: 4,
                ..Default::default()
            },
        );
        assert_eq!(one, four);
    }
}
