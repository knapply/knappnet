//! The quadratic assignment procedure (QAP) for dyadic data.
//!
//! QAP tests an association between n-by-n matrices (networks, or any
//! dyadic data such as distances or similarities) against a reference
//! distribution formed by relabeling the nodes of one matrix: the rows and
//! columns are permuted simultaneously, which keeps the dependence among
//! dyads that share a node intact (Hubert and Schultz 1976,
//! `hubert1976quadratic`; Hubert 1987, `hubert1987assignment`;
//! Krackhardt 1987, `krackhardt1987qap`).
//!
//! Two procedures are provided:
//!
//! - [`qap_cor`]: the Pearson correlation of two matrices over their dyads,
//!   with a QAP reference distribution formed by permuting `x`. This is the
//!   test sna computes with `qaptest(..., gcor)`.
//! - [`mrqap`]: multiple regression QAP. Ordinary least squares on the
//!   vectorized dyads gives the coefficients, standard errors,
//!   t-statistics, and R-squared, which match `lm()` on the vectorized
//!   dyads and sna's `netlm()` exactly. Permutation p-values come from one
//!   of two methods:
//!   - **Y-permutation** (Krackhardt 1988, `krackhardt1988predicting`):
//!     the response matrix is permuted and the full model refitted. This is
//!     sna's `nullhyp = "qapy"`.
//!   - **Double semi-partialling**, DSP (Dekker, Krackhardt, and Snijders
//!     2007, `dekker2007sensitivity`): for each coefficient k, the matrix of
//!     column k is replaced by its residuals from a regression on every
//!     other column (the intercept included); that residual matrix is
//!     permuted and the model with y, the other columns, and the permuted
//!     residuals is refitted. The statistic is the t-statistic of the
//!     permuted residual column. Dekker et al. show that DSP with this
//!     pivotal t-statistic is robust to collinearity and skewness, which is
//!     why it is the default. This is sna's `nullhyp = "qap"` (equivalently
//!     `"qapspp"`), including its treatment of the intercept, which is
//!     partialled and permuted like any other column. As in sna, a model
//!     with a single column (one predictor and no intercept, or an
//!     intercept alone) uses Y-permutation, since there is nothing to
//!     partial out.
//!
//! Conventions:
//! - **Dyads.** In [`QapMode::Digraph`] mode every ordered pair (i, j) with
//!   i != j is a dyad. In [`QapMode::Graph`] mode (undirected or symmetric
//!   data) each unordered pair is used once, read from the lower triangle
//!   (cells with row > column), which is sna's convention; the upper
//!   triangle is ignored, and permuted matrices are also read as symmetric
//!   from their lower triangle. For symmetric data this is the same as
//!   using the upper triangle. The diagonal (self-dyads) is excluded unless
//!   `diag` is set; under simultaneous permutation of rows and columns
//!   diagonal cells map only to diagonal cells.
//! - **Missing values** (NaN, which is how R stores `NA`) are handled by
//!   listwise deletion of dyads: each fit or correlation uses the dyads
//!   observed in every matrix that enters it. A matrix's missing cells move
//!   with it when it is permuted, so a permuted fit may use a different set
//!   of dyads than the observed one. In DSP the residuals of column k come
//!   from a regression on the other columns only (y does not enter it), so
//!   they are defined on the dyads where every predictor is observed. These
//!   rules reproduce sna (`netlm()` drops incomplete rows in every fit;
//!   `gcor()` uses `use = "complete.obs"`). With missing values in y only,
//!   DSP also matches sna. With a missing value in a predictor, sna's DSP
//!   is not well defined (it assigns a residual vector shorter than the
//!   cells it fills, so R recycles it), and results are not comparable.
//! - **Permutations.** Monte Carlo draws use [`crate::rng::Stream`]: draw i
//!   (0-based) is the permutation produced by stream i of the seed, so
//!   results are identical for any thread count. Draws are uniform over all
//!   n! permutations and may repeat, and the identity can occur. The
//!   exhaustive mode instead evaluates every one of the n! permutations
//!   exactly once (n <= 10), which includes the identity (the observed
//!   labeling); its p-values are the exact permutation p-values.
//! - **p-values** are reported in sna's three forms, as proportions of the
//!   reference distribution: Pr(stat <= observed), Pr(stat >= observed),
//!   and Pr(|stat| >= |observed|). The reference distribution is the set of
//!   draws only; the observed statistic is not added to it, which is sna's
//!   convention. Counts are returned as well, so the form that includes the
//!   observed statistic, (count + 1) / (draws + 1), which is a valid
//!   p-value for any number of draws (Phipson and Smyth 2010,
//!   `phipson2010permutation`), can be computed from them. Two statistics
//!   are treated as equal when they differ by at most 1e-10 times
//!   max(1, |observed|), so that permutations giving the same statistic in
//!   exact arithmetic (automorphisms, binary data) count as ties despite
//!   rounding. A draw whose statistic is undefined (NaN: for example, a
//!   permuted design that is rank deficient, or fewer dyads than
//!   coefficients after deletion) is left out of the reference
//!   distribution and of the denominator; the number of valid draws is
//!   reported.
//! - **Least squares** uses Householder QR without pivoting. A design whose
//!   column loses all but a fraction 1e-7 of its norm to the preceding
//!   columns (sna's default `tol`) is rank deficient, which is an error for
//!   the observed fit. Standard errors use the residual variance
//!   rss / (N - p), with N the number of dyads used and p the number of
//!   columns. R-squared, adjusted R-squared, and the F statistic follow
//!   `summary.lm()`: with an intercept the model sum of squares is taken
//!   about the mean of the fitted values, and without one about zero.

use crate::error::{Error, Result};
use crate::matrix::Dense;
use crate::parallel;
use crate::rng::Stream;

/// Relative tolerance for treating a permuted statistic as tied with the
/// observed one.
const TIE_EPS: f64 = 1e-10;

/// Largest n for exhaustive enumeration (10! = 3,628,800 permutations).
pub const MAX_EXHAUSTIVE_N: usize = 10;

/// Which dyads of an n-by-n matrix are used.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QapMode {
    /// Every ordered pair: directed data.
    #[default]
    Digraph,
    /// Each unordered pair once, read from the lower triangle: undirected
    /// or symmetric data.
    Graph,
}

impl QapMode {
    /// Decodes the C ABI code: 0 = digraph, 1 = graph.
    pub fn from_code(code: i32) -> Result<QapMode> {
        match code {
            0 => Ok(QapMode::Digraph),
            1 => Ok(QapMode::Graph),
            _ => Err(Error::InvalidArgument(format!(
                "unknown QAP mode code {code}"
            ))),
        }
    }
}

/// How the permutation reference distribution is formed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Permutations {
    /// `reps` uniformly random permutations; draw i uses stream i of `seed`.
    Random {
        /// Number of draws.
        reps: usize,
        /// Seed of the counter-based generator.
        seed: u64,
    },
    /// Every one of the n! permutations exactly once (n <= 10).
    Exhaustive,
}

impl Default for Permutations {
    fn default() -> Self {
        Permutations::Random {
            reps: 1000,
            seed: 0,
        }
    }
}

/// Options for [`qap_cor`].
#[derive(Debug, Clone, Copy, Default)]
pub struct QapCorOptions {
    /// Which dyads are used.
    pub mode: QapMode,
    /// Whether self-dyads (the diagonal) are used.
    pub diag: bool,
    /// How the reference distribution is formed.
    pub permutations: Permutations,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
    /// Keep every draw's statistic in the result.
    pub keep_draws: bool,
}

/// Counts of a reference distribution relative to an observed statistic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Tail {
    /// Draws with statistic <= observed.
    pub le: u64,
    /// Draws with statistic >= observed.
    pub ge: u64,
    /// Draws with |statistic| >= |observed|.
    pub abs_ge: u64,
    /// Draws with a defined (non-NaN) statistic.
    pub valid: u64,
}

impl Tail {
    fn add(&mut self, stat: f64, observed: f64) {
        if stat.is_nan() {
            return;
        }
        let eps = TIE_EPS * observed.abs().max(1.0);
        self.valid += 1;
        if stat <= observed + eps {
            self.le += 1;
        }
        if stat >= observed - eps {
            self.ge += 1;
        }
        if stat.abs() >= observed.abs() - eps {
            self.abs_ge += 1;
        }
    }

    fn merge(&mut self, other: &Tail) {
        self.le += other.le;
        self.ge += other.ge;
        self.abs_ge += other.abs_ge;
        self.valid += other.valid;
    }

    fn ratio(count: u64, valid: u64) -> f64 {
        if valid == 0 {
            f64::NAN
        } else {
            count as f64 / valid as f64
        }
    }

    /// Pr(stat <= observed) over the valid draws.
    pub fn p_le(&self) -> f64 {
        Tail::ratio(self.le, self.valid)
    }

    /// Pr(stat >= observed) over the valid draws.
    pub fn p_ge(&self) -> f64 {
        Tail::ratio(self.ge, self.valid)
    }

    /// Pr(|stat| >= |observed|) over the valid draws.
    pub fn p_abs(&self) -> f64 {
        Tail::ratio(self.abs_ge, self.valid)
    }
}

/// Result of [`qap_cor`].
#[derive(Debug, Clone)]
pub struct QapCorResult {
    /// Observed Pearson correlation.
    pub statistic: f64,
    /// Number of dyads used for the observed correlation.
    pub dyads: usize,
    /// Number of permutations evaluated.
    pub draws: u64,
    /// Counts relative to the observed correlation.
    pub tail: Tail,
    /// Each draw's correlation, in draw order, when `keep_draws` is set;
    /// otherwise empty.
    pub distribution: Vec<f64>,
}

/// Method used for MRQAP p-values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MrqapMethod {
    /// Double semi-partialling (Dekker, Krackhardt, and Snijders 2007).
    #[default]
    DoubleSemiPartialling,
    /// Permutation of the response (Krackhardt 1988).
    YPermutation,
}

impl MrqapMethod {
    /// Decodes the C ABI code: 0 = DSP, 1 = Y-permutation.
    pub fn from_code(code: i32) -> Result<MrqapMethod> {
        match code {
            0 => Ok(MrqapMethod::DoubleSemiPartialling),
            1 => Ok(MrqapMethod::YPermutation),
            _ => Err(Error::InvalidArgument(format!(
                "unknown MRQAP method code {code}"
            ))),
        }
    }
}

/// Test statistic for MRQAP p-values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MrqapStatistic {
    /// The t-statistic (pivotal; recommended by Dekker et al. 2007).
    #[default]
    T,
    /// The raw coefficient.
    Beta,
}

impl MrqapStatistic {
    /// Decodes the C ABI code: 0 = t, 1 = beta.
    pub fn from_code(code: i32) -> Result<MrqapStatistic> {
        match code {
            0 => Ok(MrqapStatistic::T),
            1 => Ok(MrqapStatistic::Beta),
            _ => Err(Error::InvalidArgument(format!(
                "unknown MRQAP statistic code {code}"
            ))),
        }
    }
}

/// Options for [`mrqap`].
#[derive(Debug, Clone, Copy)]
pub struct MrqapOptions {
    /// Which dyads are used.
    pub mode: QapMode,
    /// Whether self-dyads (the diagonal) are used.
    pub diag: bool,
    /// Whether to include an intercept (a matrix of ones, first column).
    pub intercept: bool,
    /// How p-values are computed.
    pub method: MrqapMethod,
    /// Which statistic the permutation test uses.
    pub statistic: MrqapStatistic,
    /// How the reference distribution is formed.
    pub permutations: Permutations,
    /// Relative tolerance for detecting a rank-deficient design.
    pub tol: f64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
    /// Keep every draw's statistics in the result.
    pub keep_draws: bool,
}

impl Default for MrqapOptions {
    fn default() -> Self {
        MrqapOptions {
            mode: QapMode::Digraph,
            diag: false,
            intercept: true,
            method: MrqapMethod::DoubleSemiPartialling,
            statistic: MrqapStatistic::T,
            permutations: Permutations::default(),
            tol: 1e-7,
            threads: 0,
            keep_draws: false,
        }
    }
}

/// Result of [`mrqap`]. Vectors have one entry per column of the design:
/// the intercept first (if any), then the predictors in order.
#[derive(Debug, Clone)]
pub struct MrqapResult {
    /// OLS coefficients.
    pub coefficients: Vec<f64>,
    /// Standard errors.
    pub std_errors: Vec<f64>,
    /// t-statistics.
    pub t: Vec<f64>,
    /// Counts of the reference distribution of the chosen statistic.
    pub tails: Vec<Tail>,
    /// Number of dyads used in the observed fit.
    pub dyads: usize,
    /// Residual degrees of freedom (dyads minus columns).
    pub df_residual: usize,
    /// Residual sum of squares.
    pub rss: f64,
    /// Model sum of squares (about the mean with an intercept).
    pub mss: f64,
    /// R-squared.
    pub r_squared: f64,
    /// Adjusted R-squared.
    pub adj_r_squared: f64,
    /// Residual standard error.
    pub sigma: f64,
    /// F statistic for all non-intercept coefficients.
    pub f_statistic: f64,
    /// Method actually used (Y-permutation when the design has one column).
    pub method: MrqapMethod,
    /// Number of permutations evaluated.
    pub draws: u64,
    /// Each draw's statistics, row-major (draw by column), when
    /// `keep_draws` is set; otherwise empty. NaN marks an undefined draw.
    pub distribution: Vec<f64>,
}

// ---- Dyads ------------------------------------------------------------

/// The cells (row, column) used as dyads.
fn dyad_cells(n: usize, mode: QapMode, diag: bool) -> Vec<(usize, usize)> {
    let mut cells = Vec::new();
    for c in 0..n {
        for r in 0..n {
            let keep = match mode {
                QapMode::Digraph => r != c || diag,
                QapMode::Graph => r > c || (diag && r == c),
            };
            if keep {
                cells.push((r, c));
            }
        }
    }
    cells
}

/// Value of a cell, reading graph-mode matrices from the lower triangle.
#[inline]
fn cell(m: &Dense, mode: QapMode, r: usize, c: usize) -> f64 {
    match mode {
        QapMode::Graph if r < c => m.get(c, r),
        _ => m.get(r, c),
    }
}

/// Values of `m` at `cells` after relabeling node v as `perm[v]`
/// (the permuted matrix is m[perm, perm], as in sna's `rmperm()`).
fn gather(
    m: &Dense,
    mode: QapMode,
    cells: &[(usize, usize)],
    perm: Option<&[usize]>,
    out: &mut [f64],
) {
    match perm {
        None => {
            for (o, &(r, c)) in out.iter_mut().zip(cells) {
                *o = cell(m, mode, r, c);
            }
        }
        Some(p) => {
            for (o, &(r, c)) in out.iter_mut().zip(cells) {
                *o = cell(m, mode, p[r], p[c]);
            }
        }
    }
}

fn check_square(what: &str, m: &Dense, n: usize) -> Result<()> {
    if m.rows() != n || m.cols() != n {
        return Err(Error::InvalidArgument(format!(
            "{what} is {}x{}, but a {n}x{n} matrix is required",
            m.rows(),
            m.cols()
        )));
    }
    Ok(())
}

fn check_permutations(p: Permutations, n: usize) -> Result<u64> {
    match p {
        Permutations::Random { reps, .. } => Ok(reps as u64),
        Permutations::Exhaustive => {
            if n > MAX_EXHAUSTIVE_N {
                return Err(Error::InvalidArgument(format!(
                    "exhaustive enumeration needs n <= {MAX_EXHAUSTIVE_N}; n is {n}"
                )));
            }
            Ok((1..=n as u64).product())
        }
    }
}

/// Permutation number `k` (0-based) of `0..perm.len()` in lexicographic
/// order, decoded from the factorial number system.
fn unrank(mut k: u64, perm: &mut [usize], pool: &mut Vec<usize>) {
    let n = perm.len();
    pool.clear();
    pool.extend(0..n);
    let mut fact: u64 = (1..n as u64).product();
    for (i, slot) in perm.iter_mut().enumerate() {
        let d = (k / fact) as usize;
        k %= fact;
        *slot = pool.remove(d);
        if i + 1 < n {
            fact /= (n - 1 - i) as u64;
        }
    }
}

/// Fills `perm` with permutation `item` of the reference distribution.
fn draw_permutation(p: Permutations, item: usize, perm: &mut [usize], pool: &mut Vec<usize>) {
    match p {
        Permutations::Random { seed, .. } => Stream::new(seed, item as u64).permutation(perm),
        Permutations::Exhaustive => unrank(item as u64, perm, pool),
    }
}

// ---- Correlation --------------------------------------------------------

/// Pearson correlation over the pairs where both values are observed.
/// Returns (r, pairs used); r is NaN if either is constant.
fn pearson(a: &[f64], b: &[f64]) -> (f64, usize) {
    let (mut n, mut sa, mut sb) = (0usize, 0.0, 0.0);
    for (&x, &y) in a.iter().zip(b) {
        if !x.is_nan() && !y.is_nan() {
            n += 1;
            sa += x;
            sb += y;
        }
    }
    if n < 2 {
        return (f64::NAN, n);
    }
    let (ma, mb) = (sa / n as f64, sb / n as f64);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (&x, &y) in a.iter().zip(b) {
        if !x.is_nan() && !y.is_nan() {
            let (dx, dy) = (x - ma, y - mb);
            sab += dx * dy;
            saa += dx * dx;
            sbb += dy * dy;
        }
    }
    if saa == 0.0 || sbb == 0.0 {
        return (f64::NAN, n);
    }
    (sab / (saa.sqrt() * sbb.sqrt()), n)
}

/// QAP test of the correlation between two n-by-n matrices.
///
/// See the module documentation for conventions.
pub fn qap_cor(y: &Dense, x: &Dense, opts: QapCorOptions) -> Result<QapCorResult> {
    let n = y.rows();
    check_square("y", y, n)?;
    check_square("x", x, n)?;
    let draws = check_permutations(opts.permutations, n)?;
    let cells = dyad_cells(n, opts.mode, opts.diag);
    let mut yv = vec![0.0; cells.len()];
    let mut xv = vec![0.0; cells.len()];
    gather(y, opts.mode, &cells, None, &mut yv);
    gather(x, opts.mode, &cells, None, &mut xv);
    let (observed, dyads) = pearson(&yv, &xv);
    if observed.is_nan() {
        return Err(Error::InvalidArgument(
            "the correlation is undefined: fewer than two complete dyads, or a matrix is \
             constant over them"
                .into(),
        ));
    }

    let keep = opts.keep_draws;
    struct State {
        tail: Tail,
        kept: Vec<(usize, f64)>,
        perm: Vec<usize>,
        pool: Vec<usize>,
        xp: Vec<f64>,
    }
    let states = parallel::for_each_item(
        draws as usize,
        opts.threads,
        || State {
            tail: Tail::default(),
            kept: Vec::new(),
            perm: vec![0; n],
            pool: Vec::with_capacity(n),
            xp: vec![0.0; cells.len()],
        },
        |s, item| {
            draw_permutation(opts.permutations, item, &mut s.perm, &mut s.pool);
            gather(x, opts.mode, &cells, Some(&s.perm), &mut s.xp);
            let (r, _) = pearson(&yv, &s.xp);
            s.tail.add(r, observed);
            if keep {
                s.kept.push((item, r));
            }
        },
    );

    let mut tail = Tail::default();
    let mut distribution = if keep {
        vec![0.0; draws as usize]
    } else {
        Vec::new()
    };
    for s in &states {
        tail.merge(&s.tail);
        for &(i, r) in &s.kept {
            distribution[i] = r;
        }
    }
    Ok(QapCorResult {
        statistic: observed,
        dyads,
        draws,
        tail,
        distribution,
    })
}

// ---- Least squares ------------------------------------------------------

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// An ordinary least squares fit.
struct Fit {
    coef: Vec<f64>,
    se: Vec<f64>,
    /// Residuals at every row; NaN at rows dropped for missing values.
    residuals: Vec<f64>,
    /// Fitted values at the rows used.
    fitted_used: Vec<f64>,
    rows: usize,
    rss: f64,
}

impl Fit {
    fn stat(&self, k: usize, which: MrqapStatistic) -> f64 {
        match which {
            MrqapStatistic::T => self.coef[k] / self.se[k],
            MrqapStatistic::Beta => self.coef[k],
        }
    }
}

/// Reusable buffers for [`ols`].
#[derive(Default)]
struct Workspace {
    a: Vec<f64>,
    used: Vec<usize>,
    diag: Vec<f64>,
    rinv: Vec<f64>,
}

/// Least squares of `y` on `cols` (all of the same length) by Householder
/// QR, after deleting rows with a missing value in `y` or any column.
/// Returns `None` if the design is rank deficient (a column keeps at most a
/// fraction `tol` of its norm after the preceding columns are projected
/// out) or has no residual degrees of freedom.
fn ols(cols: &[&[f64]], y: &[f64], tol: f64, ws: &mut Workspace) -> Option<Fit> {
    let p = cols.len();
    let len = y.len();
    ws.used.clear();
    for i in 0..len {
        if !y[i].is_nan() && cols.iter().all(|c| !c[i].is_nan()) {
            ws.used.push(i);
        }
    }
    let rows = ws.used.len();
    if rows <= p {
        return None;
    }
    // Column-major [X | y], rows by (p + 1).
    ws.a.clear();
    for c in cols.iter().copied().chain(std::iter::once(y)) {
        ws.a.extend(ws.used.iter().map(|&i| c[i]));
    }
    let a = &mut ws.a;
    ws.diag.clear();
    ws.diag.resize(p, 0.0);
    let diag = &mut ws.diag;
    for j in 0..p {
        let (left, right) = a.split_at_mut((j + 1) * rows);
        let col = &mut left[j * rows..];
        let orig = dot(col, col).sqrt();
        let alpha = dot(&col[j..], &col[j..]).sqrt();
        if orig == 0.0 || alpha <= tol * orig {
            return None;
        }
        // Reflector mapping col[j..] to s e1 with s = -sign(col[j]) alpha.
        let s = if col[j] >= 0.0 { -alpha } else { alpha };
        col[j] -= s;
        let v = &col[j..];
        let vv = dot(v, v);
        for k in 0..(p - j) {
            let target = &mut right[k * rows + j..(k + 1) * rows];
            let f = 2.0 * dot(v, target) / vv;
            for (t, vi) in target.iter_mut().zip(v) {
                *t -= f * vi;
            }
        }
        diag[j] = s;
    }
    // R[i][j] is a[j * rows + i] above the diagonal and diag[i] on it.
    let r = |i: usize, j: usize| if i == j { diag[i] } else { a[j * rows + i] };
    let qty = &a[p * rows..(p + 1) * rows];
    let rss: f64 = qty[p..].iter().map(|v| v * v).sum();
    let mut coef = vec![0.0; p];
    for i in (0..p).rev() {
        let mut s = qty[i];
        for (j, c) in coef.iter().enumerate().skip(i + 1) {
            s -= r(i, j) * c;
        }
        coef[i] = s / r(i, i);
    }
    // diag((X'X)^-1) = squared row norms of R^-1, which is upper triangular.
    ws.rinv.clear();
    ws.rinv.resize(p * p, 0.0);
    let rinv = &mut ws.rinv;
    for j in 0..p {
        rinv[j * p + j] = 1.0 / r(j, j);
        for i in (0..j).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s += r(i, k) * rinv[k * p + j];
            }
            rinv[i * p + j] = -s / r(i, i);
        }
    }
    let sigma2 = rss / (rows - p) as f64;
    let se = (0..p)
        .map(|i| {
            let v: f64 = (i..p).map(|j| rinv[i * p + j] * rinv[i * p + j]).sum();
            (v * sigma2).sqrt()
        })
        .collect();
    let mut residuals = vec![f64::NAN; len];
    let mut fitted_used = Vec::with_capacity(rows);
    for &i in &ws.used {
        let f: f64 = cols.iter().zip(&coef).map(|(c, b)| c[i] * b).sum();
        fitted_used.push(f);
        residuals[i] = y[i] - f;
    }
    Some(Fit {
        coef,
        se,
        residuals,
        fitted_used,
        rows,
        rss,
    })
}

// ---- MRQAP --------------------------------------------------------------

/// Multiple regression QAP of `y` on the matrices in `x`.
///
/// See the module documentation for conventions.
pub fn mrqap(y: &Dense, x: &[Dense], opts: MrqapOptions) -> Result<MrqapResult> {
    let n = y.rows();
    check_square("y", y, n)?;
    for (i, m) in x.iter().enumerate() {
        check_square(&format!("predictor {}", i + 1), m, n)?;
    }
    if !(opts.tol >= 0.0 && opts.tol < 1.0) {
        return Err(Error::InvalidArgument("tol must be in [0, 1)".into()));
    }
    let draws = check_permutations(opts.permutations, n)?;
    let cells = dyad_cells(n, opts.mode, opts.diag);
    let len = cells.len();

    // Design columns (intercept first), vectorized over the dyads.
    let vectorize = |m: &Dense| {
        let mut v = vec![0.0; len];
        gather(m, opts.mode, &cells, None, &mut v);
        v
    };
    let yv = vectorize(y);
    let mut xv: Vec<Vec<f64>> = Vec::new();
    if opts.intercept {
        xv.push(vec![1.0; len]);
    }
    xv.extend(x.iter().map(vectorize));
    let p = xv.len();
    if p == 0 {
        return Err(Error::InvalidArgument("the model has no columns".into()));
    }
    let xs: Vec<&[f64]> = xv.iter().map(|v| v.as_slice()).collect();

    let mut ws = Workspace::default();
    let fit = ols(&xs, &yv, opts.tol, &mut ws).ok_or_else(|| {
        Error::InvalidArgument(
            "the design is rank deficient, or there are no more complete dyads than \
             coefficients"
                .into(),
        )
    })?;
    let t: Vec<f64> = (0..p).map(|k| fit.coef[k] / fit.se[k]).collect();
    let observed: Vec<f64> = (0..p).map(|k| fit.stat(k, opts.statistic)).collect();

    // Summary statistics as in summary.lm().
    let df_residual = fit.rows - p;
    let rdf = df_residual as f64;
    let mss: f64 = if opts.intercept {
        let m = fit.fitted_used.iter().sum::<f64>() / fit.rows as f64;
        fit.fitted_used.iter().map(|f| (f - m) * (f - m)).sum()
    } else {
        fit.fitted_used.iter().map(|f| f * f).sum()
    };
    let df_int = if opts.intercept { 1.0 } else { 0.0 };
    let r_squared = mss / (mss + fit.rss);
    let adj_r_squared = 1.0 - (1.0 - r_squared) * ((fit.rows as f64 - df_int) / rdf);
    let f_statistic = if p as f64 > df_int {
        (mss / (p as f64 - df_int)) / (fit.rss / rdf)
    } else {
        f64::NAN
    };
    let sigma = (fit.rss / rdf).sqrt();

    // With one column there is nothing to partial out.
    let method = if p == 1 {
        MrqapMethod::YPermutation
    } else {
        opts.method
    };

    // For DSP, the residual matrix of each column on the others. Residuals
    // exist where every predictor is observed and are NaN elsewhere.
    let mut residual_mats: Vec<Dense> = Vec::new();
    if method == MrqapMethod::DoubleSemiPartialling {
        for k in 0..p {
            let others: Vec<&[f64]> = (0..p).filter(|&j| j != k).map(|j| xs[j]).collect();
            let part = ols(&others, xs[k], opts.tol, &mut ws).ok_or_else(|| {
                Error::InvalidArgument(format!(
                    "column {} cannot be partialled on the other columns",
                    k + 1
                ))
            })?;
            let mut m = Dense::from_row_major(n, n, vec![f64::NAN; n * n])?;
            for (&(r, c), &e) in cells.iter().zip(&part.residuals) {
                m.set(r, c, e);
            }
            residual_mats.push(m);
        }
    }

    let keep = opts.keep_draws;
    struct State {
        tails: Vec<Tail>,
        kept: Vec<(usize, Vec<f64>)>,
        perm: Vec<usize>,
        pool: Vec<usize>,
        buf: Vec<f64>,
        ws: Workspace,
    }
    let states = parallel::for_each_item(
        draws as usize,
        opts.threads,
        || State {
            tails: vec![Tail::default(); p],
            kept: Vec::new(),
            perm: vec![0; n],
            pool: Vec::with_capacity(n),
            buf: vec![0.0; len],
            ws: Workspace::default(),
        },
        |s, item| {
            draw_permutation(opts.permutations, item, &mut s.perm, &mut s.pool);
            let mut stats = vec![f64::NAN; p];
            match method {
                MrqapMethod::YPermutation => {
                    gather(y, opts.mode, &cells, Some(&s.perm), &mut s.buf);
                    if let Some(f) = ols(&xs, &s.buf, opts.tol, &mut s.ws) {
                        for (k, st) in stats.iter_mut().enumerate() {
                            *st = f.stat(k, opts.statistic);
                        }
                    }
                }
                MrqapMethod::DoubleSemiPartialling => {
                    for (k, st) in stats.iter_mut().enumerate() {
                        let rm = &residual_mats[k];
                        gather(rm, opts.mode, &cells, Some(&s.perm), &mut s.buf);
                        let mut design: Vec<&[f64]> =
                            (0..p).filter(|&j| j != k).map(|j| xs[j]).collect();
                        design.push(&s.buf);
                        if let Some(f) = ols(&design, &yv, opts.tol, &mut s.ws) {
                            *st = f.stat(p - 1, opts.statistic);
                        }
                    }
                }
            }
            for (tail, (&st, &obs)) in s.tails.iter_mut().zip(stats.iter().zip(&observed)) {
                tail.add(st, obs);
            }
            if keep {
                s.kept.push((item, stats));
            }
        },
    );

    let mut tails = vec![Tail::default(); p];
    let mut distribution = if keep {
        vec![0.0; draws as usize * p]
    } else {
        Vec::new()
    };
    for s in &states {
        for (t, st) in tails.iter_mut().zip(&s.tails) {
            t.merge(st);
        }
        for (i, stats) in &s.kept {
            distribution[i * p..(i + 1) * p].copy_from_slice(stats);
        }
    }

    Ok(MrqapResult {
        coefficients: fit.coef,
        std_errors: fit.se,
        t,
        tails,
        dyads: fit.rows,
        df_residual,
        rss: fit.rss,
        mss,
        r_squared,
        adj_r_squared,
        sigma,
        f_statistic,
        method,
        draws,
        distribution,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat(n: usize, f: impl Fn(usize, usize) -> f64) -> Dense {
        let mut m = Dense::zeros(n, n);
        for r in 0..n {
            for c in 0..n {
                m.set(r, c, f(r, c));
            }
        }
        m
    }

    #[test]
    fn dyads_follow_mode_and_diag() {
        assert_eq!(dyad_cells(3, QapMode::Digraph, false).len(), 6);
        assert_eq!(dyad_cells(3, QapMode::Digraph, true).len(), 9);
        assert_eq!(
            dyad_cells(3, QapMode::Graph, false),
            vec![(1, 0), (2, 0), (2, 1)]
        );
        assert_eq!(dyad_cells(3, QapMode::Graph, true).len(), 6);
    }

    #[test]
    fn unrank_enumerates_every_permutation_once() {
        let mut seen = std::collections::HashSet::new();
        let (mut perm, mut pool) = (vec![0; 4], Vec::new());
        for k in 0..24 {
            unrank(k, &mut perm, &mut pool);
            assert!(seen.insert(perm.clone()));
        }
        unrank(0, &mut perm, &mut pool);
        assert_eq!(perm, vec![0, 1, 2, 3]);
        unrank(23, &mut perm, &mut pool);
        assert_eq!(perm, vec![3, 2, 1, 0]);
    }

    #[test]
    fn ols_matches_hand_computation() {
        // y = 1 + 2x exactly, then perturbed: simple regression by hand.
        let x = [0.0, 1.0, 2.0, 3.0];
        let y = [1.0, 3.5, 4.5, 7.0];
        let ones = [1.0; 4];
        let f = ols(&[&ones, &x], &y, 1e-7, &mut Workspace::default()).unwrap();
        // slope = Sxy / Sxx = 9.5 / 5, intercept = ybar - slope * xbar.
        assert!((f.coef[1] - 1.9).abs() < 1e-12);
        assert!((f.coef[0] - (4.0 - 1.9 * 1.5)).abs() < 1e-12);
        let rss: f64 = (0..4).map(|i| f.residuals[i].powi(2)).sum();
        assert!((f.rss - rss).abs() < 1e-12);
        // se(slope) = sqrt(rss / (n - 2) / Sxx).
        assert!((f.se[1] - (rss / 2.0 / 5.0).sqrt()).abs() < 1e-12);
        // A rank-deficient design is rejected.
        assert!(ols(&[&ones, &ones], &y, 1e-7, &mut Workspace::default()).is_none());
    }

    #[test]
    fn identical_matrices_correlate_perfectly() {
        let a = mat(6, |r, c| ((r * 7 + c * 3) % 5) as f64);
        let res = qap_cor(
            &a,
            &a,
            QapCorOptions {
                permutations: Permutations::Exhaustive,
                ..Default::default()
            },
        )
        .unwrap();
        assert!((res.statistic - 1.0).abs() < 1e-12);
        assert_eq!(res.draws, 720);
        // Only relabelings that preserve the matrix reach r = 1, including
        // the identity, so Pr(r >= 1) is at least 1/720.
        assert!(res.tail.ge >= 1);
        assert_eq!(res.tail.le, 720);
    }

    #[test]
    fn results_do_not_depend_on_threads() {
        let y = mat(9, |r, c| ((r * 31 + c * 17) % 11) as f64);
        let x1 = mat(9, |r, c| ((r * 13 + c * 5) % 7) as f64);
        let x2 = mat(9, |r, c| (r as f64 - c as f64).abs());
        for method in [
            MrqapMethod::DoubleSemiPartialling,
            MrqapMethod::YPermutation,
        ] {
            let run = |threads| {
                mrqap(
                    &y,
                    &[x1.clone(), x2.clone()],
                    MrqapOptions {
                        method,
                        permutations: Permutations::Random {
                            reps: 200,
                            seed: 42,
                        },
                        threads,
                        keep_draws: true,
                        ..Default::default()
                    },
                )
                .unwrap()
            };
            let (a, b) = (run(1), run(4));
            assert_eq!(a.tails, b.tails);
            assert_eq!(a.distribution, b.distribution);
        }
        let run = |threads| {
            qap_cor(
                &y,
                &x1,
                QapCorOptions {
                    permutations: Permutations::Random { reps: 300, seed: 7 },
                    threads,
                    keep_draws: true,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        assert_eq!(run(1).distribution, run(3).distribution);
    }

    #[test]
    fn one_column_falls_back_to_y_permutation() {
        let y = mat(5, |r, c| (r + 2 * c) as f64);
        let res = mrqap(&y, &[], MrqapOptions::default()).unwrap();
        assert_eq!(res.method, MrqapMethod::YPermutation);
        assert!(res.f_statistic.is_nan());
    }

    #[test]
    fn graph_mode_reads_lower_triangle() {
        let y = mat(4, |r, c| if r > c { (r * c) as f64 } else { 99.0 });
        let x = mat(4, |r, c| if r > c { (r + c) as f64 } else { -5.0 });
        let a = qap_cor(
            &y,
            &x,
            QapCorOptions {
                mode: QapMode::Graph,
                permutations: Permutations::Exhaustive,
                ..Default::default()
            },
        )
        .unwrap();
        let ys = mat(4, |r, c| (r.max(c) * r.min(c)) as f64);
        let xs = mat(4, |r, c| (r + c) as f64);
        let b = qap_cor(
            &ys,
            &xs,
            QapCorOptions {
                mode: QapMode::Graph,
                permutations: Permutations::Exhaustive,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(a.statistic, b.statistic);
        assert_eq!(a.tail, b.tail);
        assert_eq!(a.dyads, 6);
    }
}
