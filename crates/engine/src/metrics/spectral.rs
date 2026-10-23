//! Spectral and walk-based centralities: PageRank, Katz centrality,
//! Bonacich's power centrality, and Kleinberg's hubs and authorities.
//!
//! Sources:
//! - PageRank (Brin and Page 1998, `brin1998anatomy`): the stationary
//!   distribution of a random walk that follows an out-arc with probability
//!   `damping` (chosen in proportion to arc weights) and otherwise jumps to
//!   a node drawn from the personalization vector p (uniform by default).
//!   A walker at a node without out-arcs (a dangling node) jumps instead;
//!   the treatment of dangling nodes follows Langville and Meyer (2004,
//!   `langville2004deeper`). Computed by power iteration from the uniform
//!   distribution until the scores, which sum to 1, change by at most `tol`
//!   in total (L1 norm); the error is then at most damping / (1 - damping)
//!   times that change.
//! - Katz centrality (Katz 1953, `katz1953new`) in the form of NetworkX's
//!   `katz_centrality` and igraph's `alpha_centrality` (Bonacich and Lloyd
//!   2001, `bonacich2001eigenvector`): x = alpha A' x + beta, so x =
//!   beta (1 + alpha A' 1 + alpha^2 A'^2 1 + ...), counting walks ending at
//!   each node with weight alpha^length. Katz's own index omits the term
//!   for walks of length 0 (x / beta - 1 for constant beta). The series
//!   converges only when alpha < 1 / lambda_max, the reciprocal of the
//!   spectral radius of A; larger values are an error. Solved by the
//!   fixed-point iteration x <- alpha A' x + beta.
//! - Bonacich's power centrality (Bonacich 1987, `bonacich1987power`):
//!   c(alpha, beta) = alpha (I - beta A)^-1 A 1, with alpha chosen so that
//!   the sum of squared scores equals the number of nodes, as in sna's
//!   `bonpow()` and igraph's `power_centrality()`. beta may be negative
//!   (power from being connected to weak alters). When |beta| lambda_max < 1
//!   the system is solved by the fixed-point iteration y <- A 1 + beta A y;
//!   otherwise the series diverges but the linear system usually still has
//!   a solution, which is found by Gaussian elimination with partial
//!   pivoting on the dense matrix (up to 5,000 nodes). A singular system is
//!   an error.
//! - Hubs and authorities (Kleinberg 1999, `kleinberg1999authoritative`):
//!   authority scores are the leading eigenvector of A^T A and hub scores of
//!   A A^T; equivalently, the leading right and left singular vectors of A.
//!   Kleinberg's iteration alternates a <- A^T h and h <- A a from h = 1.
//!   Here the authority vector is computed with the Lanczos method of
//!   [`crate::metrics::eigenvector`] on A^T A (never formed), one block of
//!   the bipartite hub-authority graph at a time, and the hub vector is
//!   A a.
//!
//! Spectral radii are bounded with Collatz's inequalities (Collatz 1942,
//! `collatz1942einschliessungssatz`): for an irreducible nonnegative matrix
//! and any positive x, min_i (A x)_i / x_i <= lambda_max <= max_i
//! (A x)_i / x_i. The spectral radius of a reducible matrix is the largest
//! over its strongly connected components, and the bounds are tightened by
//! power iteration on each component until they decide whether alpha (or
//! |beta|) lambda_max < 1.
//!
//! # Conventions
//!
//! | | igraph | NetworkX | sna |
//! |---|---|---|---|
//! | PageRank, undirected self-loop | adds 2 to the node's degree (two arc ends) | adds 1 | (none) |
//! | PageRank, dangling nodes | jump by the personalization vector | jump by the personalization vector (or a separate `dangling` vector) | (none) |
//! | PageRank, parallel edges | add up | add up | (none) |
//! | Katz | `alpha_centrality(alpha, exo)`: unnormalized, in-ties, self-loops ignored by default | `katz_centrality`: in-ties, unit Euclidean norm by default, self-loops count once | (none) |
//! | Bonacich power | `power_centrality`: out-ties, diagonal ignored unless `loops = TRUE`, parallel edges add up | (none) | `bonpow`: out-ties, adjacency matrix (parallel edges collapse), diagonal ignored unless `diag = TRUE` |
//! | Hubs and authorities scale | maximum 1 | scores sum to 1 | (none) |
//!
//! Further conventions:
//! - Weights are used as given (larger means a stronger tie) and must be
//!   nonnegative; a node whose out-arcs all have weight 0 is dangling.
//! - Directed graphs: PageRank, Katz, and authorities follow arcs into a
//!   node; Bonacich power and hubs follow arcs out of it. PageRank with
//!   `directed = false` and the other measures with mode `All` ignore
//!   direction, so opposite arcs become two parallel edges.
//! - Self-loops ([`LoopWeight`]): in undirected graphs a loop adds twice
//!   its weight to the diagonal (igraph) or once (NetworkX), or nothing; in
//!   directed graphs a loop counts once unless ignored.
//! - Disconnected graphs: PageRank is defined on every graph. Katz and
//!   Bonacich scores are sums over walks and need no connectivity. Hubs and
//!   authorities, like eigenvector centrality, are unique only when the
//!   largest singular value of A is simple; when several blocks of the
//!   hub-authority graph share it (to a relative 1e-9), the result is the
//!   limit of Kleinberg's iteration from h = 1: each such block contributes
//!   its unit authority vector v times v^T (A^T 1).
//! - Undefined results are NaN: hubs and authorities of a graph without
//!   positive arcs, and Bonacich scores when A 1 = 0 (no ties), since alpha
//!   is then 0 / 0.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::linalg;
use crate::matrix::Dense;
use crate::metrics::closeness::oriented;
use crate::metrics::degree::DegreeMode;
use crate::metrics::eigenvector::{
    check_nonnegative, symmetric_leading, LoopWeight, SparseOp, Stop,
};
use crate::metrics::topography::{components, ComponentKind};

/// Default tolerance of the iterative solvers.
pub const DEFAULT_TOL: f64 = 1e-12;
/// Default iteration cap.
pub const DEFAULT_MAX_ITER: usize = 10_000;

/// Where a PageRank walker at a node without out-arcs jumps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Dangling {
    /// By the personalization vector (igraph, NetworkX).
    #[default]
    Personalized,
    /// To any node with equal probability, whatever the personalization.
    Uniform,
}

impl Dangling {
    /// Decodes the C interface codes: 0 = personalized, 1 = uniform.
    pub fn from_code(code: i32) -> Result<Dangling> {
        match code {
            0 => Ok(Dangling::Personalized),
            1 => Ok(Dangling::Uniform),
            _ => Err(Error::InvalidArgument(format!(
                "unknown dangling-node code {code}"
            ))),
        }
    }
}

/// Options for [`pagerank`].
#[derive(Debug, Clone)]
pub struct PagerankOptions {
    /// Probability of following an arc (0 <= damping <= 1).
    pub damping: f64,
    /// Follow arc directions in directed graphs.
    pub directed: bool,
    /// Use edge weights as given.
    pub weighted: bool,
    /// Jump distribution, one nonnegative value per node, rescaled to sum
    /// to 1; `None` is uniform.
    pub personalization: Option<Vec<f64>>,
    /// Treatment of dangling nodes.
    pub dangling: Dangling,
    /// Self-loops in undirected graphs.
    pub loops: LoopWeight,
    /// Convergence tolerance on the L1 change of the scores.
    pub tol: f64,
    /// Iteration cap.
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for PagerankOptions {
    fn default() -> Self {
        PagerankOptions {
            damping: 0.85,
            directed: true,
            weighted: false,
            personalization: None,
            dangling: Dangling::Personalized,
            loops: LoopWeight::Twice,
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
            threads: 0,
        }
    }
}

fn check_iteration(tol: f64, max_iter: usize) -> Result<()> {
    if tol.is_nan() || tol <= 0.0 || max_iter == 0 {
        return Err(Error::InvalidArgument(
            "the tolerance and the iteration cap must be positive".into(),
        ));
    }
    Ok(())
}

/// The matrix whose row i lists the arcs a measure follows *into* i (mode
/// `In`), out of i (`Out`), or both ways (`All`, and undirected graphs).
fn matrix(
    g: &Graph,
    mode: DegreeMode,
    weighted: bool,
    loops: LoopWeight,
    threads: usize,
) -> Result<SparseOp> {
    let n = g.node_count();
    let mode = if g.is_directed() {
        mode
    } else {
        DegreeMode::All
    };
    let h = oriented(
        g,
        if mode == DegreeMode::All {
            DegreeMode::All
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
    Ok(SparseOp::from_rows(
        rows,
        n,
        weighted,
        loops.factor(h.is_directed()),
        threads,
    ))
}

/// PageRank. See the module documentation.
pub fn pagerank(g: &Graph, opts: &PagerankOptions) -> Result<Vec<f64>> {
    check_iteration(opts.tol, opts.max_iter)?;
    if !(0.0..=1.0).contains(&opts.damping) {
        return Err(Error::InvalidArgument(format!(
            "damping is {}; it must be between 0 and 1",
            opts.damping
        )));
    }
    let n = g.node_count();
    if n == 0 {
        return Ok(Vec::new());
    }
    if opts.weighted {
        check_nonnegative(g, "PageRank")?;
    }
    let p = match &opts.personalization {
        None => vec![1.0 / n as f64; n],
        Some(p) => probabilities(p, n, "personalization")?,
    };
    let jump = match opts.dangling {
        Dangling::Personalized => p.clone(),
        Dangling::Uniform => vec![1.0 / n as f64; n],
    };
    let mode = if opts.directed {
        DegreeMode::In
    } else {
        DegreeMode::All
    };
    let op = matrix(g, mode, opts.weighted, opts.loops, opts.threads)?;
    // Out-strength of j: the sum of column j of the in-arc matrix.
    let mut out = vec![0.0f64; n];
    for i in 0..n {
        for (j, w) in op.row(i) {
            out[j] += w;
        }
    }
    let d = opts.damping;
    let mut x = vec![1.0 / n as f64; n];
    let mut z = vec![0.0; n];
    let mut y = vec![0.0; n];
    for _ in 0..opts.max_iter {
        let mut stuck = 0.0;
        for j in 0..n {
            if out[j] > 0.0 {
                z[j] = x[j] / out[j];
            } else {
                z[j] = 0.0;
                stuck += x[j];
            }
        }
        op.apply(&z, &mut y);
        let mut total = 0.0;
        for i in 0..n {
            y[i] = d * (y[i] + stuck * jump[i]) + (1.0 - d) * p[i];
            total += y[i];
        }
        let mut change = 0.0;
        for i in 0..n {
            y[i] /= total;
            change += (y[i] - x[i]).abs();
        }
        std::mem::swap(&mut x, &mut y);
        if change <= opts.tol {
            return Ok(x);
        }
    }
    Err(Error::InvalidArgument(format!(
        "PageRank did not converge to tolerance {} within {} iterations",
        opts.tol, opts.max_iter
    )))
}

/// Checks a vector of nonnegative weights and rescales it to sum to 1.
fn probabilities(p: &[f64], n: usize, what: &str) -> Result<Vec<f64>> {
    if p.len() != n {
        return Err(Error::LengthMismatch {
            what: "personalization",
            expected: n,
            found: p.len(),
        });
    }
    if p.iter().any(|&v| !v.is_finite() || v < 0.0) {
        return Err(Error::InvalidArgument(format!(
            "the {what} vector must hold finite, nonnegative values"
        )));
    }
    let total: f64 = p.iter().sum();
    if total <= 0.0 {
        return Err(Error::InvalidArgument(format!(
            "the {what} vector must have a positive sum"
        )));
    }
    Ok(p.iter().map(|v| v / total).collect())
}

/// Options for [`katz`].
#[derive(Debug, Clone)]
pub struct KatzOptions {
    /// Attenuation factor; must satisfy 0 <= alpha < 1 / lambda_max.
    pub alpha: f64,
    /// Exogenous term: one value for every node, or one per node.
    pub beta: Vec<f64>,
    /// Direction: `In` (the usual definition), `Out`, or `All`.
    pub mode: DegreeMode,
    /// Use edge weights as given.
    pub weighted: bool,
    /// Self-loops.
    pub loops: LoopWeight,
    /// Scale to unit Euclidean norm (NetworkX's default).
    pub normalized: bool,
    /// Convergence tolerance (relative change of the largest score).
    pub tol: f64,
    /// Iteration cap.
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for KatzOptions {
    fn default() -> Self {
        KatzOptions {
            alpha: 0.1,
            beta: vec![1.0],
            mode: DegreeMode::In,
            weighted: false,
            loops: LoopWeight::Once,
            normalized: true,
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
            threads: 0,
        }
    }
}

/// Expands a constant or per-node vector.
fn per_node(v: &[f64], n: usize, what: &'static str) -> Result<Vec<f64>> {
    let out = match v.len() {
        1 => vec![v[0]; n],
        len if len == n => v.to_vec(),
        len => {
            return Err(Error::LengthMismatch {
                what,
                expected: n,
                found: len,
            })
        }
    };
    if out.iter().any(|x| !x.is_finite()) {
        return Err(Error::InvalidArgument(format!("{what} must be finite")));
    }
    Ok(out)
}

/// Katz centrality. See the module documentation.
pub fn katz(g: &Graph, opts: &KatzOptions) -> Result<Vec<f64>> {
    check_iteration(opts.tol, opts.max_iter)?;
    let n = g.node_count();
    if !opts.alpha.is_finite() || opts.alpha < 0.0 {
        return Err(Error::InvalidArgument(format!(
            "alpha is {}; it must be a nonnegative number below 1 / lambda_max",
            opts.alpha
        )));
    }
    let beta = per_node(&opts.beta, n, "beta")?;
    if n == 0 {
        return Ok(Vec::new());
    }
    if opts.weighted {
        check_nonnegative(g, "Katz centrality")?;
    }
    let op = matrix(g, opts.mode, opts.weighted, opts.loops, opts.threads)?;
    let q = if opts.alpha > 0.0 {
        opts.alpha * check_attenuation(g, &op, opts.alpha, "alpha", opts.max_iter)?
    } else {
        0.0
    };
    let mut x = solve_series(
        &op,
        &beta,
        opts.alpha,
        q,
        opts.tol,
        opts.max_iter,
        "Katz centrality",
    )?;
    if opts.normalized {
        let norm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm > 0.0 {
            x.iter_mut().for_each(|v| *v /= norm);
        }
    }
    Ok(x)
}

/// Solves x = b + t A x by fixed-point iteration, for |t| lambda_max at
/// most q, with q below 1. Stops when the change, times q / (1 - q), is at
/// most `tol` times the largest |x_i|; for a matrix whose norm is below q
/// this bounds the error.
fn solve_series(
    op: &SparseOp,
    b: &[f64],
    t: f64,
    q: f64,
    tol: f64,
    max_iter: usize,
    what: &str,
) -> Result<Vec<f64>> {
    let n = b.len();
    let mut x = b.to_vec();
    if t == 0.0 {
        return Ok(x);
    }
    let factor = (q / (1.0 - q)).max(1.0);
    let mut y = vec![0.0; n];
    for _ in 0..max_iter {
        op.apply(&x, &mut y);
        let (mut change, mut top) = (0.0f64, 0.0f64);
        for i in 0..n {
            y[i] = b[i] + t * y[i];
            change = change.max((y[i] - x[i]).abs());
            top = top.max(y[i].abs());
        }
        std::mem::swap(&mut x, &mut y);
        if change * factor <= tol * top || top == 0.0 {
            return Ok(x);
        }
    }
    Err(Error::InvalidArgument(format!(
        "{what} did not converge to tolerance {tol} within {max_iter} iterations"
    )))
}

/// Checks that `t` times the spectral radius of `op` is below 1 and
/// returns an upper bound on the spectral radius (below 1 / t). The radius
/// is the largest over strongly connected components; each is bounded with
/// Collatz's inequalities, tightened by power iteration on A + I.
fn check_attenuation(g: &Graph, op: &SparseOp, t: f64, name: &str, max_iter: usize) -> Result<f64> {
    let n = op.dim();
    let limit = 1.0 / t;
    // Strong components of the matrix's pattern (the transpose and the
    // symmetrized pattern have the same components as the graph).
    let comp = if g.is_directed() {
        let h = pattern_graph(op)?;
        components(&h, ComponentKind::Strong)
    } else {
        components(g, ComponentKind::Weak)
    };
    let k = comp.count();
    let mut start = vec![0usize; k + 1];
    for &c in &comp.membership {
        start[c as usize + 1] += 1;
    }
    for c in 0..k {
        start[c + 1] += start[c];
    }
    let mut members = vec![0u32; n];
    let mut cursor = start.clone();
    for (v, &c) in comp.membership.iter().enumerate() {
        members[cursor[c as usize]] = v as u32;
        cursor[c as usize] += 1;
    }
    let mut index = vec![u32::MAX; n];
    let mut radius = 0.0f64;
    for c in 0..k {
        let nodes = &members[start[c]..start[c + 1]];
        if nodes.len() == 1 {
            let v = nodes[0] as usize;
            let diag: f64 = op.row(v).filter(|&(j, _)| j == v).map(|(_, w)| w).sum();
            if diag >= limit {
                return Err(too_large(name, t, diag));
            }
            radius = radius.max(diag);
            continue;
        }
        for (i, &v) in nodes.iter().enumerate() {
            index[v as usize] = i as u32;
        }
        let sub = op.submatrix(nodes, &index);
        for &v in nodes {
            index[v as usize] = u32::MAX;
        }
        let m = sub.dim();
        let mut x = vec![1.0; m];
        let mut y = vec![0.0; m];
        let mut decided = false;
        for _ in 0..max_iter {
            sub.apply(&x, &mut y);
            let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
            let mut top = 0.0f64;
            for i in 0..m {
                let r = y[i] / x[i];
                lo = lo.min(r);
                hi = hi.max(r);
                y[i] += x[i];
                top = top.max(y[i]);
            }
            if lo >= limit {
                return Err(too_large(name, t, lo));
            }
            if hi < limit {
                radius = radius.max(hi);
                decided = true;
                break;
            }
            for i in 0..m {
                x[i] = y[i] / top;
            }
        }
        if !decided {
            return Err(Error::InvalidArgument(format!(
                "{name} = {t} is too close to 1 / lambda_max to decide whether the \
                 series converges; choose a smaller value"
            )));
        }
    }
    Ok(radius)
}

fn too_large(name: &str, t: f64, lambda: f64) -> Error {
    Error::InvalidArgument(format!(
        "{name} = {t} is not below 1 / lambda_max: the largest eigenvalue of the \
         adjacency matrix is at least {lambda}, so 1 / lambda_max is at most {}",
        1.0 / lambda
    ))
}

/// A directed graph with the nonzero pattern of `op` (an arc i -> j for
/// each positive entry (i, j)).
fn pattern_graph(op: &SparseOp) -> Result<Graph> {
    let n = op.dim();
    let (mut s, mut t) = (Vec::new(), Vec::new());
    for i in 0..n {
        for (j, w) in op.row(i) {
            if w > 0.0 {
                s.push(i as i64);
                t.push(j as i64);
            }
        }
    }
    Graph::from_edges(n, &s, &t, None, true)
}

/// Options for [`bonacich_power`].
#[derive(Debug, Clone, Copy)]
pub struct BonacichOptions {
    /// The attenuation (beta, igraph's and sna's `exponent`); may be
    /// negative.
    pub exponent: f64,
    /// Direction: `Out` (Bonacich, igraph, sna), `In`, or `All`.
    pub mode: DegreeMode,
    /// Use edge weights as given.
    pub weighted: bool,
    /// Self-loops (the diagonal). igraph and sna ignore it by default.
    pub loops: LoopWeight,
    /// Divide by the sum of the scores instead of scaling the sum of
    /// squares to n.
    pub rescale: bool,
    /// Tolerance of the iterative solver.
    pub tol: f64,
    /// Iteration cap.
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for BonacichOptions {
    fn default() -> Self {
        BonacichOptions {
            exponent: 1.0,
            mode: DegreeMode::Out,
            weighted: false,
            loops: LoopWeight::Ignore,
            rescale: false,
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
            threads: 0,
        }
    }
}

/// Largest network solved by dense elimination when the power series
/// diverges.
const DENSE_MAX: usize = 5000;

/// Bonacich's power centrality. See the module documentation.
pub fn bonacich_power(g: &Graph, opts: BonacichOptions) -> Result<Vec<f64>> {
    check_iteration(opts.tol, opts.max_iter)?;
    if !opts.exponent.is_finite() {
        return Err(Error::InvalidArgument("the exponent must be finite".into()));
    }
    let n = g.node_count();
    if n == 0 {
        return Ok(Vec::new());
    }
    if opts.weighted {
        check_nonnegative(g, "power centrality")?;
    }
    // Row i of A lists the arcs out of i for mode `Out` (A[i][j] = 1 for
    // an arc i -> j), into i for `In`.
    let op = matrix(g, opts.mode, opts.weighted, opts.loops, opts.threads)?;
    let ones = vec![1.0; n];
    let mut b = vec![0.0; n];
    op.apply(&ones, &mut b);
    let beta = opts.exponent;
    let series = if beta == 0.0 {
        Some(0.0)
    } else {
        match check_attenuation(g, &op, beta.abs(), "|exponent|", opts.max_iter) {
            Ok(r) => Some(r * beta.abs()),
            Err(_) => None,
        }
    };
    let mut y = match series {
        Some(q) => solve_series(
            &op,
            &b,
            beta,
            q,
            opts.tol,
            opts.max_iter,
            "power centrality",
        )?,
        None => {
            if n > DENSE_MAX {
                return Err(Error::InvalidArgument(format!(
                    "|exponent| x lambda_max >= 1, so power centrality needs a dense \
                     linear solve, which is limited to {DENSE_MAX} nodes; use a smaller \
                     |exponent|"
                )));
            }
            dense_solve(&op, beta, &b)?
        }
    };
    let sq: f64 = y.iter().map(|v| v * v).sum();
    let alpha = (n as f64 / sq).sqrt();
    if opts.rescale {
        let total: f64 = y.iter().sum();
        y.iter_mut().for_each(|v| *v /= total);
    } else {
        y.iter_mut().for_each(|v| *v *= alpha);
    }
    Ok(y)
}

/// Solves (I - t A) y = b by Gaussian elimination with partial pivoting
/// ([`crate::linalg::solve`]).
fn dense_solve(op: &SparseOp, t: f64, b: &[f64]) -> Result<Vec<f64>> {
    let n = op.dim();
    let mut a = Dense::identity(n);
    for i in 0..n {
        for (j, w) in op.row(i) {
            a.set(i, j, a.get(i, j) - t * w);
        }
    }
    linalg::solve(&a, b, 1e-12).ok_or_else(|| {
        Error::InvalidArgument(
            "I - exponent x A is singular (exponent is the reciprocal of an \
             eigenvalue), so power centrality is undefined; choose another exponent"
                .into(),
        )
    })
}

/// How hub and authority scores are scaled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HitsScale {
    /// Largest score 1 (igraph).
    #[default]
    Max,
    /// Scores sum to 1 (NetworkX).
    Sum,
    /// Unit Euclidean norm (igraph with `scale = FALSE`).
    Unit,
}

impl HitsScale {
    /// Decodes the C interface codes: 0 = max, 1 = sum, 2 = unit.
    pub fn from_code(code: i32) -> Result<HitsScale> {
        match code {
            0 => Ok(HitsScale::Max),
            1 => Ok(HitsScale::Sum),
            2 => Ok(HitsScale::Unit),
            _ => Err(Error::InvalidArgument(format!(
                "unknown hub and authority scale code {code}"
            ))),
        }
    }
}

/// Options for [`hits`].
#[derive(Debug, Clone, Copy)]
pub struct HitsOptions {
    /// Use edge weights as given.
    pub weighted: bool,
    /// Self-loops in undirected graphs.
    pub loops: LoopWeight,
    /// Scaling of both vectors.
    pub scale: HitsScale,
    /// Convergence tolerance (relative residual of the eigenpair).
    pub tol: f64,
    /// Cap on matrix-vector products per block.
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for HitsOptions {
    fn default() -> Self {
        HitsOptions {
            weighted: false,
            loops: LoopWeight::Twice,
            scale: HitsScale::Max,
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
            threads: 0,
        }
    }
}

/// Result of [`hits`].
#[derive(Debug, Clone, PartialEq)]
pub struct Hits {
    /// Hub scores.
    pub hubs: Vec<f64>,
    /// Authority scores.
    pub authorities: Vec<f64>,
    /// The largest eigenvalue of A^T A (the squared largest singular value).
    pub value: f64,
    /// Whether the leading singular value is simple, so the scores do not
    /// depend on the starting vector.
    pub unique: bool,
}

/// Relative difference within which two blocks' singular values tie.
const TIE: f64 = 1e-9;

/// Hubs and authorities. See the module documentation.
pub fn hits(g: &Graph, opts: HitsOptions) -> Result<Hits> {
    check_iteration(opts.tol, opts.max_iter)?;
    let n = g.node_count();
    if opts.weighted {
        check_nonnegative(g, "hub and authority scores")?;
    }
    // Rows of `a` are arcs out of each node.
    let a = matrix(g, DegreeMode::Out, opts.weighted, opts.loops, opts.threads)?;
    let at = a.transpose(n);
    // Blocks of the bipartite graph with hub copies 0..n and authority
    // copies n..2n, joined by every positive entry.
    let mut parent: Vec<usize> = (0..2 * n).collect();
    fn find(p: &mut [usize], mut v: usize) -> usize {
        while p[v] != v {
            p[v] = p[p[v]];
            v = p[v];
        }
        v
    }
    let mut any = false;
    for i in 0..n {
        for (j, w) in a.row(i) {
            if w > 0.0 {
                any = true;
                let (x, y) = (find(&mut parent, i), find(&mut parent, n + j));
                if x != y {
                    parent[x.max(y)] = x.min(y);
                }
            }
        }
    }
    if !any {
        return Ok(Hits {
            hubs: vec![f64::NAN; n],
            authorities: vec![f64::NAN; n],
            value: 0.0,
            unique: false,
        });
    }
    // Group authorities (nodes with a positive in-entry) and hubs by block
    // root, in node order.
    let has_in: Vec<bool> = (0..n).map(|j| at.row(j).any(|(_, w)| w > 0.0)).collect();
    let has_out: Vec<bool> = (0..n).map(|i| a.row(i).any(|(_, w)| w > 0.0)).collect();
    let mut block_of = vec![usize::MAX; 2 * n];
    let mut blocks: Vec<(Vec<u32>, Vec<u32>)> = Vec::new();
    for v in 0..2 * n {
        let present = if v < n { has_out[v] } else { has_in[v - n] };
        if !present {
            continue;
        }
        let r = find(&mut parent, v);
        if block_of[r] == usize::MAX {
            block_of[r] = blocks.len();
            blocks.push((Vec::new(), Vec::new()));
        }
        let b = &mut blocks[block_of[r]];
        if v < n {
            b.0.push(v as u32);
        } else {
            b.1.push((v - n) as u32);
        }
    }
    // Bounds on each block's sigma^2: above, (max row sum)(max column sum)
    // of A; below, the largest column sum of squares (||A e_j||^2).
    let row_sum: Vec<f64> = (0..n).map(|i| a.row(i).map(|(_, w)| w).sum()).collect();
    let col_sum: Vec<f64> = (0..n).map(|j| at.row(j).map(|(_, w)| w).sum()).collect();
    let col_sq: Vec<f64> = (0..n).map(|j| merged_squares(&at, j)).collect();
    let k = blocks.len();
    let mut upper = vec![0.0f64; k];
    let mut lower = vec![0.0f64; k];
    for (b, (hubs, auths)) in blocks.iter().enumerate() {
        let r = hubs
            .iter()
            .map(|&i| row_sum[i as usize])
            .fold(0.0, f64::max);
        let c = auths
            .iter()
            .map(|&j| col_sum[j as usize])
            .fold(0.0, f64::max);
        upper[b] = r * c;
        lower[b] = auths
            .iter()
            .map(|&j| col_sq[j as usize])
            .fold(0.0, f64::max);
    }
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&x, &y| lower[y].total_cmp(&lower[x]).then(x.cmp(&y)));

    let mut hub_index = vec![u32::MAX; n];
    let mut auth_index = vec![u32::MAX; n];
    let mut best = 0.0f64;
    let mut solved: Vec<(usize, f64, Vec<f64>, Vec<f64>)> = Vec::new();
    for &b in &order {
        if upper[b] < best * (1.0 - TIE) {
            continue;
        }
        let (hubs, auths) = &blocks[b];
        for (i, &v) in hubs.iter().enumerate() {
            hub_index[v as usize] = i as u32;
        }
        for (i, &v) in auths.iter().enumerate() {
            auth_index[v as usize] = i as u32;
        }
        let ab = a.submatrix(hubs, &auth_index);
        let abt = at.submatrix(auths, &hub_index);
        for &v in hubs {
            hub_index[v as usize] = u32::MAX;
        }
        for &v in auths {
            auth_index[v as usize] = u32::MAX;
        }
        // Kleinberg's first authority vector, A^T 1, restricted.
        let ones = vec![1.0; hubs.len()];
        let mut s = vec![0.0; auths.len()];
        abt.apply(&ones, &mut s);
        let apply = |x: &[f64], y: &mut [f64]| {
            let mut t = vec![0.0; ab.dim()];
            ab.apply(x, &mut t);
            abt.apply(&t, y);
        };
        let threads = ab.threads.max(abt.threads);
        let res = symmetric_leading(
            auths.len(),
            &apply,
            &s,
            &Stop {
                tol: opts.tol,
                shift: 0.0,
                max_matvecs: opts.max_iter,
            },
            threads,
        )
        .ok_or_else(|| {
            Error::InvalidArgument(format!(
                "hub and authority scores did not converge to tolerance {} within {} \
                 matrix-vector products; a larger iteration cap may help",
                opts.tol, opts.max_iter
            ))
        })?;
        best = best.max(res.value);
        solved.push((b, res.value, res.vector, s));
    }
    let mut auth = vec![0.0; n];
    let mut tied = 0usize;
    for (b, value, v, s) in &solved {
        if *value < best * (1.0 - TIE) {
            continue;
        }
        tied += 1;
        let weight: f64 = v.iter().zip(s).map(|(x, y)| x * y).sum();
        for (&j, &vj) in blocks[*b].1.iter().zip(v) {
            auth[j as usize] = vj * weight;
        }
    }
    let mut hubs = vec![0.0; n];
    a.apply(&auth, &mut hubs);
    rescale(&mut auth, opts.scale);
    rescale(&mut hubs, opts.scale);
    Ok(Hits {
        hubs,
        authorities: auth,
        value: best,
        unique: tied == 1,
    })
}

/// Sum of squares of row v with parallel entries merged.
fn merged_squares(op: &SparseOp, v: usize) -> f64 {
    let mut q = 0.0;
    let mut run: Option<(usize, f64)> = None;
    for (j, w) in op.row(v) {
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
    q
}

fn rescale(x: &mut [f64], how: HitsScale) {
    for v in x.iter_mut() {
        if *v < 0.0 {
            *v = 0.0;
        }
    }
    let d = match how {
        HitsScale::Max => x.iter().copied().fold(0.0f64, f64::max),
        HitsScale::Sum => x.iter().sum(),
        HitsScale::Unit => x.iter().map(|v| v * v).sum::<f64>().sqrt(),
    };
    if d > 0.0 {
        x.iter_mut().for_each(|v| *v /= d);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)
    }

    // 0 -> 1, 0 -> 2, 1 -> 2, 2 -> 0, 3 -> 2, 2 -> 4; node 4 is dangling.
    fn small() -> Graph {
        Graph::from_edges(5, &[0, 0, 1, 2, 3, 2], &[1, 2, 2, 0, 2, 4], None, true).unwrap()
    }

    #[test]
    fn pagerank_hand_values() {
        // Dense solutions of the defining linear system.
        let g = small();
        let x = pagerank(&g, &PagerankOptions::default()).unwrap();
        let want = [
            0.2142011097,
            0.1574496602,
            0.3477339318,
            0.0664141886,
            0.2142011097,
        ];
        assert!(close(&x, &want, 1e-9));
        let p = PagerankOptions {
            personalization: Some(vec![1.0, 0.0, 0.0, 0.0, 0.0]),
            ..Default::default()
        };
        let want = [0.3928645968, 0.1669674536, 0.3088897892, 0.0, 0.1312781604];
        assert!(close(&pagerank(&g, &p).unwrap(), &want, 1e-9));
        let u = PagerankOptions {
            dangling: Dangling::Uniform,
            ..p
        };
        let want = [
            0.3166510816,
            0.1629073936,
            0.3254597594,
            0.0283306839,
            0.1666510816,
        ];
        assert!(close(&pagerank(&g, &u).unwrap(), &want, 1e-9));
        let bad = PagerankOptions {
            damping: 1.5,
            ..Default::default()
        };
        assert!(pagerank(&g, &bad).is_err());
    }

    #[test]
    fn katz_counts_walks_and_checks_alpha() {
        // 0 -> 1 -> 2: x = 1 + alpha A^T x gives 1, 1 + a, 1 + a + a^2.
        let dag = Graph::from_edges(3, &[0, 1], &[1, 2], None, true).unwrap();
        let o = KatzOptions {
            alpha: 10.0,
            normalized: false,
            ..Default::default()
        };
        assert_eq!(katz(&dag, &o).unwrap(), vec![1.0, 11.0, 111.0]);
        // lambda_max of `small` is about 1.3247 (the plastic number).
        let g = small();
        let ok = KatzOptions {
            alpha: 0.2,
            normalized: false,
            loops: LoopWeight::Ignore,
            ..Default::default()
        };
        let want = [1.3445378151, 1.2689075630, 1.7226890756, 1.0, 1.3445378151];
        assert!(close(&katz(&g, &ok).unwrap(), &want, 1e-9));
        let too_big = KatzOptions {
            alpha: 0.76,
            ..ok.clone()
        };
        assert!(katz(&g, &too_big).is_err());
    }

    #[test]
    fn power_centrality_hand_values() {
        // Star around 0 with an extra edge 1 - 4.
        let g = Graph::from_edges(5, &[0, 0, 0, 1], &[1, 2, 3, 4], None, false).unwrap();
        let run = |exponent, rescale| {
            bonacich_power(
                &g,
                BonacichOptions {
                    exponent,
                    rescale,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let want = [
            1.4696735043,
            1.1238679739,
            0.7780624435,
            0.7780624435,
            0.6051596782,
        ];
        assert!(close(&run(0.5, false), &want, 1e-9));
        let want = [
            1.8605210188,
            -0.6201736729,
            -0.6201736729,
            -0.6201736729,
            0.6201736729,
        ];
        assert!(close(&run(-0.5, false), &want, 1e-9));
        let s: f64 = run(0.5, true).iter().sum();
        assert!((s - 1.0).abs() < 1e-12);
        // exponent 1 is beyond 1 / lambda_max (about 1.9), so the dense
        // solver is used; check that (I - A) y is proportional to A 1.
        let y = run(1.0, false);
        let a = [
            [0., 1., 1., 1., 0.],
            [1., 0., 0., 0., 1.],
            [1., 0., 0., 0., 0.],
            [1., 0., 0., 0., 0.],
            [0., 1., 0., 0., 0.],
        ];
        let a1: Vec<f64> = a.iter().map(|r| r.iter().sum()).collect();
        let lhs: Vec<f64> = (0..5)
            .map(|i| y[i] - (0..5).map(|j| a[i][j] * y[j]).sum::<f64>())
            .collect();
        let k = lhs[0] / a1[0];
        assert!((0..5).all(|i| (lhs[i] - k * a1[i]).abs() < 1e-9));
        // A directed cycle has eigenvalue 1, so exponent 1 is singular.
        let cyc = Graph::from_edges(3, &[0, 1, 2], &[1, 2, 0], None, true).unwrap();
        assert!(bonacich_power(&cyc, BonacichOptions::default()).is_err());
    }

    #[test]
    fn hits_hand_values_and_ties() {
        let r = hits(&small(), HitsOptions::default()).unwrap();
        let h = 0.5f64.sqrt();
        assert!(close(&r.hubs, &[1.0, h, 0.0, h, 0.0], 1e-9));
        assert!(close(
            &r.authorities,
            &[0.0, 2f64.sqrt() - 1.0, 1.0, 0.0, 0.0],
            1e-9
        ));
        assert!(r.unique);
        // A directed 4-cycle: A^T A = I, every block ties, and Kleinberg's
        // iteration from equal scores stays put.
        let cyc = Graph::from_edges(4, &[0, 1, 2, 3], &[1, 2, 3, 0], None, true).unwrap();
        let r = hits(&cyc, HitsOptions::default()).unwrap();
        assert!(!r.unique);
        assert_eq!(r.hubs, vec![1.0; 4]);
        let empty = Graph::from_edges(3, &[], &[], None, true).unwrap();
        let r = hits(&empty, HitsOptions::default()).unwrap();
        assert!(r.hubs.iter().all(|v| v.is_nan()));
    }
}
