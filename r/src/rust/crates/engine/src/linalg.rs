//! Dense linear algebra on [`Dense`] matrices: products, Cholesky and
//! Gauss-Jordan inverses, Gaussian elimination, the symmetric
//! eigendecomposition, and the Moore-Penrose inverse.
//!
//! These routines serve the small systems of the estimation methods (the
//! p-by-p matrices of model fits), the n-by-n problems of
//! classical scaling, power centrality, and the leading eigenvector
//! method, and the projected eigenproblems of the Lanczos solver. They are
//! O(n^3) algorithms with O(n^2) memory, not a replacement for LAPACK on
//! large matrices. The product is blocked for cache reuse, and the
//! eigendecomposition works on the transpose of the EISPACK working matrix
//! so that its inner loops run over contiguous memory; neither change
//! alters any arithmetic operation or its order, so results do not depend
//! on these choices.
//!
//! Sources:
//! - Cholesky factorization and the triangular solves: Golub and Van Loan
//!   2013, `golub2013matrix`, Algorithm 4.2.1 and section 3.1.
//! - Gaussian elimination and Gauss-Jordan inversion with partial
//!   pivoting: Golub and Van Loan 2013, section 3.4.
//! - Symmetric eigendecomposition: Householder reduction to tridiagonal
//!   form followed by the implicit QL algorithm, the procedures `tred2`
//!   (Martin, Reinsch, and Wilkinson 1968, `martin1968householder`) and
//!   `tql2` (Bowdler, Martin, Reinsch, and Wilkinson 1968, `bowdler1968ql`)
//!   as published in EISPACK (Smith et al. 1976, `smith1976matrix`, which
//!   is in the public domain); Golub and Van Loan 2013, sections 8.3 and
//!   8.2.
//! - Moore-Penrose inverse through the eigendecomposition: Golub and Van
//!   Loan 2013, section 5.5.
//! - Conjugate gradients for sparse symmetric positive semidefinite
//!   systems given as a matrix-vector product, with a diagonal (Jacobi)
//!   preconditioner: Hestenes and Stiefel 1952, `hestenes1952methods`;
//!   Golub and Van Loan 2013, section 11.5.
//!
//! Numerical properties:
//! - [`cholesky`] succeeds exactly when every pivot is positive and finite,
//!   so it doubles as the test for positive definiteness to working
//!   precision; no symmetry check is made, and only the lower triangle of
//!   the input is read.
//! - [`invert`] and [`solve`] treat a pivot as zero when its magnitude is
//!   at most `tol` times the largest magnitude in the input matrix, and
//!   then report singularity instead of dividing by rounding noise.
//!   Callers pick `tol` (1e-12 or 1e-13 in this crate) and the message.
//! - [`symmetric_eigen`] returns eigenvalues in increasing order with
//!   orthonormal eigenvectors; the residual ||A v - lambda v|| is of order
//!   n times machine precision times ||A||. The sign of each eigenvector
//!   is a by-product of the rotations, so callers that need a fixed sign
//!   choose one. The QL iteration fails only when 60 iterations on one
//!   eigenvalue do not reduce the subdiagonal to rounding level.
//! - [`pseudo_inverse`] symmetrizes its input first and drops eigenvalues
//!   at most `tol` times the largest in magnitude, returning the rank it
//!   kept.
//! - [`conjugate_gradient`] stops when the residual norm is at most `tol`
//!   times the norm of the right-hand side. A singular positive
//!   semidefinite system (a graph Laplacian) is solved when the right-hand
//!   side is orthogonal to the null space, starting from zero; the
//!   solution is then the one orthogonal to the null space in the
//!   preconditioned inner product.

use crate::error::{Error, Result};
use crate::matrix::Dense;

// ---- Instruction sets ------------------------------------------------------

/// Kernels compiled for wider vector instructions, chosen at run time.
///
/// `dispatch!(kernel(args))` calls `kernel` through a copy compiled with
/// AVX2 enabled when the processor supports it (x86-64), and through the
/// baseline copy otherwise. The kernel must be `#[inline(always)]` so that
/// its body is compiled inside each copy. Only element-wise loops are
/// vectorized: the compiler never reassociates floating-point sums or
/// fuses a multiply and an add without being told to, so every copy
/// performs the same IEEE operations in the same order and the results are
/// identical bit for bit on every processor.
pub(crate) mod simd {
    /// Calls `$f(args)` through an AVX2 build when available.
    macro_rules! dispatch {
        ($f:ident ( $($arg:expr),* $(,)? )) => {{
            #[cfg(target_arch = "x86_64")]
            {
                #[target_feature(enable = "avx2")]
                unsafe fn wide<T>(f: impl FnOnce() -> T) -> T {
                    f()
                }
                if std::is_x86_feature_detected!("avx2") {
                    // SAFETY: the processor supports AVX2.
                    unsafe { wide(|| $f($($arg),*)) }
                } else {
                    $f($($arg),*)
                }
            }
            #[cfg(not(target_arch = "x86_64"))]
            {
                $f($($arg),*)
            }
        }};
    }
    pub(crate) use dispatch;
}

// ---- Products ---------------------------------------------------------------

/// The transpose of `a`.
pub fn transpose(a: &Dense) -> Dense {
    let mut out = Dense::zeros(a.cols(), a.rows());
    for r in 0..a.rows() {
        for c in 0..a.cols() {
            out.set(c, r, a.get(r, c));
        }
    }
    out
}

/// The product A B. Panics unless the inner dimensions agree.
///
/// The loops are blocked over columns of B and over the inner dimension
/// so that a panel of B stays in cache while every row of A passes over
/// it. Each entry still accumulates its terms in increasing order of the
/// inner index, so the result is that of the plain row-by-row product.
pub fn matmul(a: &Dense, b: &Dense) -> Dense {
    assert_eq!(a.cols(), b.rows(), "matmul: inner dimensions differ");
    const KB: usize = 64;
    const CB: usize = 512;
    let (rows, inner, cols) = (a.rows(), a.cols(), b.cols());
    let (av, bv) = (a.as_slice(), b.as_slice());
    let mut out = vec![0.0; rows * cols];
    let mut c0 = 0;
    while c0 < cols {
        let c1 = (c0 + CB).min(cols);
        let mut k0 = 0;
        while k0 < inner {
            let k1 = (k0 + KB).min(inner);
            for r in 0..rows {
                let orow = &mut out[r * cols + c0..r * cols + c1];
                for k in k0..k1 {
                    let v = av[r * inner + k];
                    if v == 0.0 {
                        continue;
                    }
                    let brow = &bv[k * cols + c0..k * cols + c1];
                    for (o, x) in orow.iter_mut().zip(brow) {
                        *o += v * x;
                    }
                }
            }
            k0 = k1;
        }
        c0 = c1;
    }
    Dense::from_row_major(rows, cols, out).expect("dimensions")
}

/// The product A x.
pub fn matvec(a: &Dense, x: &[f64]) -> Vec<f64> {
    (0..a.rows())
        .map(|r| a.row(r).iter().zip(x).map(|(u, v)| u * v).sum())
        .collect()
}

/// The quadratic form x' A x.
pub fn quadratic_form(a: &Dense, x: &[f64]) -> f64 {
    matvec(a, x).iter().zip(x).map(|(u, v)| u * v).sum()
}

// ---- Cholesky -----------------------------------------------------------------

/// The lower Cholesky factor L with A = L L', or `None` if `a` is not
/// positive definite to working precision.
pub fn cholesky(a: &Dense) -> Option<Dense> {
    let p = a.rows();
    let mut l = Dense::zeros(p, p);
    for j in 0..p {
        let mut d = a.get(j, j);
        for k in 0..j {
            d -= l.get(j, k) * l.get(j, k);
        }
        if d <= 0.0 || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        l.set(j, j, d);
        for i in j + 1..p {
            let mut s = a.get(i, j);
            for k in 0..j {
                s -= l.get(i, k) * l.get(j, k);
            }
            l.set(i, j, s / d);
        }
    }
    Some(l)
}

/// Solves A x = b for symmetric positive definite A by Cholesky
/// factorization, or `None` if A is not positive definite.
pub fn solve_spd(a: &Dense, b: &[f64]) -> Option<Vec<f64>> {
    let p = a.rows();
    let l = cholesky(a)?;
    let mut z = vec![0.0; p];
    for i in 0..p {
        let mut s = b[i];
        for (k, zk) in z.iter().enumerate().take(i) {
            s -= l.get(i, k) * zk;
        }
        z[i] = s / l.get(i, i);
    }
    let mut x = vec![0.0; p];
    for i in (0..p).rev() {
        let mut s = z[i];
        for (k, xk) in x.iter().enumerate().skip(i + 1) {
            s -= l.get(k, i) * xk;
        }
        x[i] = s / l.get(i, i);
    }
    Some(x)
}

/// The inverse of a symmetric positive definite matrix by Cholesky
/// factorization (A^{-1} = L^{-T} L^{-1}), or `None` if `a` is not
/// positive definite.
pub fn cholesky_inverse(a: &Dense) -> Option<Dense> {
    let p = a.rows();
    let l = cholesky(a)?;
    let mut linv = Dense::zeros(p, p);
    for c in 0..p {
        for r in c..p {
            let mut s = if r == c { 1.0 } else { 0.0 };
            for k in c..r {
                s -= l.get(r, k) * linv.get(k, c);
            }
            linv.set(r, c, s / l.get(r, r));
        }
    }
    let mut inv = Dense::zeros(p, p);
    for r in 0..p {
        for c in 0..=r {
            let mut s = 0.0;
            for k in r..p {
                s += linv.get(k, r) * linv.get(k, c);
            }
            inv.set(r, c, s);
            inv.set(c, r, s);
        }
    }
    Some(inv)
}

/// The inverse of a symmetric positive definite matrix on several threads,
/// or `None` if `a` is not positive definite (to working precision).
///
/// The Cholesky factor is computed column by column, the rows below the
/// diagonal of each column in parallel (each entry a dot product of two
/// row prefixes, so memory is read contiguously); then every column of the
/// inverse is found by a forward and a backward triangular solve, columns
/// in parallel. Every entry is computed by the same operations in the same
/// order whatever the number of threads, so results do not depend on it;
/// they can differ from [`cholesky_inverse`] in the last bits, since the
/// operations are grouped differently. O(n^3) work in all, about twice
/// that of [`cholesky_inverse`], spread over the threads.
pub fn spd_inverse_parallel(a: &Dense, threads: usize) -> Option<Dense> {
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
    let n = a.rows();
    assert_eq!(
        a.cols(),
        n,
        "spd_inverse_parallel: the matrix must be square"
    );
    if n == 0 {
        return Some(Dense::zeros(0, 0));
    }
    let threads = crate::parallel::threads_for(threads, n, (n as u64).pow(3) / 2);
    // The lower factor, row-major, shared between the members.
    let l: Vec<AtomicU64> = (0..n * n).map(|_| AtomicU64::new(0)).collect();
    let s: Vec<AtomicU64> = (0..n).map(|_| AtomicU64::new(0)).collect();
    let failed = AtomicBool::new(false);
    let get = |i: usize, j: usize| f64::from_bits(l[i * n + j].load(Ordering::Relaxed));
    crate::parallel::crew(threads, |crew| {
        let size = crew.size();
        let me = crew.index();
        for j in 0..n {
            // Rows j.. are dealt out in turn so that the work stays even as
            // the columns shrink.
            let first = j + (me + size - j % size) % size;
            let mut i = first;
            while i < n {
                let mut acc = a.get(i, j);
                for k in 0..j {
                    acc -= get(i, k) * get(j, k);
                }
                s[i].store(acc.to_bits(), Ordering::Relaxed);
                i += size;
            }
            crew.wait();
            let d = f64::from_bits(s[j].load(Ordering::Relaxed));
            if !(d > 0.0) || !d.is_finite() {
                failed.store(true, Ordering::Relaxed);
                return;
            }
            let d = d.sqrt();
            let mut i = first;
            while i < n {
                let x = f64::from_bits(s[i].load(Ordering::Relaxed));
                let v = if i == j { d } else { x / d };
                l[i * n + j].store(v.to_bits(), Ordering::Relaxed);
                i += size;
            }
            crew.wait();
        }
    });
    if failed.load(Ordering::Relaxed) {
        return None;
    }
    let l: Vec<f64> = l
        .into_iter()
        .map(|x| f64::from_bits(x.into_inner()))
        .collect();
    // U = L' row-major, for contiguous backward solves.
    let mut u = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            u[j * n + i] = l[i * n + j];
        }
    }
    let mut inv = vec![0.0; n * n];
    // Row c of the (symmetric) inverse is the solution of L L' x = e_c.
    crate::parallel::fill_chunks(
        &mut inv,
        n,
        threads,
        || vec![0.0; n],
        |y, c, row| {
            y.iter_mut().for_each(|v| *v = 0.0);
            for i in c..n {
                let li = &l[i * n..i * n + i];
                let mut acc = if i == c { 1.0 } else { 0.0 };
                for k in c..i {
                    acc -= li[k] * y[k];
                }
                y[i] = acc / l[i * n + i];
            }
            for i in (0..n).rev() {
                let ui = &u[i * n..(i + 1) * n];
                let mut acc = y[i];
                for k in i + 1..n {
                    acc -= ui[k] * row[k];
                }
                row[i] = acc / u[i * n + i];
            }
        },
    );
    Dense::from_row_major(n, n, inv).ok()
}

// ---- Gaussian elimination -----------------------------------------------------

/// The largest magnitude in `a`, floored at the smallest positive double
/// so that a relative pivot threshold is never zero.
fn magnitude(a: &Dense) -> f64 {
    a.as_slice()
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE)
}

/// The row at or below `c` holding the largest magnitude in column `c`
/// (the first such row on ties).
fn pivot_row(m: &Dense, c: usize) -> usize {
    let mut pivot = c;
    for r in c + 1..m.rows() {
        if m.get(r, c).abs() > m.get(pivot, c).abs() {
            pivot = r;
        }
    }
    pivot
}

fn swap_rows(m: &mut Dense, a: usize, b: usize) {
    if a != b {
        for k in 0..m.cols() {
            let (x, y) = (m.get(a, k), m.get(b, k));
            m.set(a, k, y);
            m.set(b, k, x);
        }
    }
}

/// The inverse of a square matrix by Gauss-Jordan elimination with partial
/// pivoting, or `None` if a pivot is at most `tol` times the largest
/// magnitude in `a` (singular to working precision).
pub fn invert(a: &Dense, tol: f64) -> Option<Dense> {
    let p = a.rows();
    assert_eq!(a.cols(), p, "invert: the matrix must be square");
    let mut m = a.clone();
    let mut inv = Dense::identity(p);
    let threshold = tol * magnitude(a);
    for c in 0..p {
        let pivot = pivot_row(&m, c);
        if m.get(pivot, c).abs() <= threshold {
            return None;
        }
        swap_rows(&mut m, c, pivot);
        swap_rows(&mut inv, c, pivot);
        let d = m.get(c, c);
        for k in 0..p {
            m.set(c, k, m.get(c, k) / d);
            inv.set(c, k, inv.get(c, k) / d);
        }
        for r in 0..p {
            if r == c {
                continue;
            }
            let f = m.get(r, c);
            if f == 0.0 {
                continue;
            }
            for k in 0..p {
                m.set(r, k, m.get(r, k) - f * m.get(c, k));
                inv.set(r, k, inv.get(r, k) - f * inv.get(c, k));
            }
        }
    }
    Some(inv)
}

/// The inverse of a symmetric matrix: by Cholesky factorization when it is
/// positive definite, otherwise by [`invert`] with threshold `tol`.
pub fn invert_symmetric(a: &Dense, tol: f64) -> Option<Dense> {
    cholesky_inverse(a).or_else(|| invert(a, tol))
}

/// Solves A x = b for square A by Gaussian elimination with partial
/// pivoting, or `None` if a pivot is at most `tol` times the largest
/// magnitude in `a`.
pub fn solve(a: &Dense, b: &[f64], tol: f64) -> Option<Vec<f64>> {
    let n = a.rows();
    assert_eq!(a.cols(), n, "solve: the matrix must be square");
    let mut m = a.clone();
    let mut x = b.to_vec();
    let threshold = tol * magnitude(a);
    for k in 0..n {
        let p = pivot_row(&m, k);
        if m.get(p, k).abs() <= threshold {
            return None;
        }
        swap_rows(&mut m, k, p);
        x.swap(k, p);
        let pivot = m.get(k, k);
        for i in k + 1..n {
            let f = m.get(i, k) / pivot;
            if f == 0.0 {
                continue;
            }
            m.set(i, k, 0.0);
            for j in k + 1..n {
                m.set(i, j, m.get(i, j) - f * m.get(k, j));
            }
            x[i] -= f * x[k];
        }
    }
    for k in (0..n).rev() {
        let mut s = x[k];
        for (j, xj) in x.iter().enumerate().skip(k + 1) {
            s -= m.get(k, j) * xj;
        }
        x[k] = s / m.get(k, k);
    }
    Some(x)
}

// ---- Symmetric eigendecomposition ------------------------------------------------

/// Eigenvalues in increasing order and the matching orthonormal
/// eigenvectors, as the columns of a square matrix (`vectors.get(i, j)`
/// is component i of eigenvector j), of the symmetric matrix `a`.
///
/// Only the symmetric part is meaningful: the reduction reads the whole
/// matrix but assumes `a.get(i, j) == a.get(j, i)`. Fails if the QL
/// iteration does not converge (see the module documentation).
pub fn symmetric_eigen(a: &Dense) -> Result<(Vec<f64>, Dense)> {
    let n = a.rows();
    assert_eq!(a.cols(), n, "symmetric_eigen: the matrix must be square");
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    if n == 0 {
        return Ok((d, a.clone()));
    }
    // The reduction and the QL sweeps work on the transpose W = V' of the
    // EISPACK working matrix V, so that the columns of V they traverse
    // (Householder vectors, eigenvectors under rotation) are contiguous
    // rows of W. Every arithmetic operation, and its order, is that of
    // EISPACK; only the memory layout differs.
    let mut w = transpose(a).as_slice().to_vec();
    tred2(&mut w, &mut d, &mut e);
    tql2(&mut w, &mut d, &mut e)?;
    let w = Dense::from_row_major(n, n, w).expect("square");
    Ok((d, transpose(&w)))
}

/// Householder reduction to tridiagonal form, accumulating the
/// transformations (EISPACK `tred2`) in `w`, which holds the transpose of
/// the EISPACK working matrix V (row-major, so V(r, c) is `w[c * n + r]`):
/// on return `d` holds the diagonal and `e[1..]` the subdiagonal.
#[allow(clippy::needless_range_loop)]
fn tred2(w: &mut [f64], d: &mut [f64], e: &mut [f64]) {
    let n = d.len();
    for j in 0..n {
        d[j] = w[j * n + n - 1];
    }
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for k in 0..i {
            scale += d[k].abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = w[j * n + i - 1];
                w[j * n + i] = 0.0;
                w[i * n + j] = 0.0;
            }
        } else {
            for k in 0..i {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for j in 0..i {
                e[j] = 0.0;
            }
            for j in 0..i {
                f = d[j];
                w[i * n + j] = f;
                let col = &w[j * n..j * n + i];
                g = e[j] + col[j] * f;
                for k in j + 1..i {
                    g += col[k] * d[k];
                    e[k] += col[k] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                let col = &mut w[j * n..j * n + n];
                for k in j..i {
                    col[k] -= f * e[k] + g * d[k];
                }
                d[j] = col[i - 1];
                col[i] = 0.0;
            }
        }
        d[i] = h;
    }
    for i in 0..n - 1 {
        w[i * n + n - 1] = w[i * n + i];
        w[i * n + i] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            let (head, tail) = w.split_at_mut((i + 1) * n);
            let hv = &tail[..=i];
            for k in 0..=i {
                d[k] = hv[k] / h;
            }
            for j in 0..=i {
                let col = &mut head[j * n..j * n + i + 1];
                let mut g = 0.0;
                for k in 0..=i {
                    g += hv[k] * col[k];
                }
                for k in 0..=i {
                    col[k] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            w[(i + 1) * n + k] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = w[j * n + n - 1];
        w[j * n + n - 1] = 0.0;
    }
    w[(n - 1) * n + n - 1] = 1.0;
    e[0] = 0.0;
}

/// Implicit QL iterations on the tridiagonal matrix (d, e), accumulating
/// the rotations in `w` (the transpose of the EISPACK matrix, as in
/// [`tred2`]), then sorting eigenvalues ascending (EISPACK `tql2`).
#[allow(clippy::needless_range_loop)]
fn tql2(w: &mut [f64], d: &mut [f64], e: &mut [f64]) -> Result<()> {
    let n = d.len();
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    let mut tst1: f64 = 0.0;
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        // e[n - 1] is 0, so m < n.
        if m > l {
            let mut iter = 0;
            loop {
                iter += 1;
                if iter > 60 {
                    return Err(Error::InvalidArgument(
                        "the symmetric eigenvalue iteration did not converge".into(),
                    ));
                }
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for i in l + 2..n {
                    d[i] -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    let (lo, hi) = w.split_at_mut((i + 1) * n);
                    let vi = &mut lo[i * n..];
                    let vi1 = &mut hi[..n];
                    for (a, b) in vi.iter_mut().zip(vi1.iter_mut()) {
                        let hk = *b;
                        *b = s * *a + c * hk;
                        *a = c * *a - s * hk;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
    // Selection sort into ascending order, moving eigenvectors along.
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for j in i + 1..n {
            if d[j] < p {
                k = j;
                p = d[j];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            let (lo, hi) = w.split_at_mut(k * n);
            lo[i * n..i * n + n].swap_with_slice(&mut hi[..n]);
        }
    }
    Ok(())
}

/// The result of [`conjugate_gradient`].
#[derive(Debug, Clone, PartialEq)]
pub struct CgSolution {
    /// The approximate solution.
    pub x: Vec<f64>,
    /// Iterations used.
    pub iterations: usize,
    /// Final relative residual ||b - A x|| / ||b||.
    pub residual: f64,
}

/// Solves A x = b by preconditioned conjugate gradients, where
/// `apply(v, out)` writes A v to `out` for a symmetric positive
/// semidefinite A, and `diag` holds the diagonal used as preconditioner
/// (entries that are not positive are treated as 1). Stops when the
/// relative residual is at most `tol` or after `max_iter` iterations.
pub fn conjugate_gradient<F>(
    apply: F,
    diag: &[f64],
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> CgSolution
where
    F: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let dot = |a: &[f64], c: &[f64]| a.iter().zip(c).map(|(x, y)| x * y).sum::<f64>();
    let inv: Vec<f64> = diag
        .iter()
        .map(|&d| if d > 0.0 { 1.0 / d } else { 1.0 })
        .collect();
    let mut x = vec![0.0; n];
    let bnorm = dot(b, b).sqrt();
    if bnorm == 0.0 {
        return CgSolution {
            x,
            iterations: 0,
            residual: 0.0,
        };
    }
    let mut r = b.to_vec();
    let mut z: Vec<f64> = r.iter().zip(&inv).map(|(r, m)| r * m).collect();
    let mut p = z.clone();
    let mut q = vec![0.0; n];
    let mut rz = dot(&r, &z);
    let mut iterations = 0;
    let mut residual = 1.0;
    while iterations < max_iter {
        apply(&p, &mut q);
        let pq = dot(&p, &q);
        if pq <= 0.0 || !pq.is_finite() {
            break;
        }
        let alpha = rz / pq;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        iterations += 1;
        residual = dot(&r, &r).sqrt() / bnorm;
        if residual <= tol {
            break;
        }
        for i in 0..n {
            z[i] = r[i] * inv[i];
        }
        let rz_next = dot(&r, &z);
        let beta = rz_next / rz;
        rz = rz_next;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    CgSolution {
        x,
        iterations,
        residual,
    }
}

/// The Moore-Penrose inverse of a symmetric matrix through its
/// eigendecomposition, with the rank used. Eigenvalues whose magnitude is
/// at most `tol` times the largest are treated as zero. The input is
/// symmetrized (averaged with its transpose) first.
pub fn pseudo_inverse(a: &Dense, tol: f64) -> Result<(Dense, usize)> {
    let p = a.rows();
    let mut sym = a.clone();
    for i in 0..p {
        for j in 0..i {
            let v = 0.5 * (a.get(i, j) + a.get(j, i));
            sym.set(i, j, v);
            sym.set(j, i, v);
        }
    }
    let (values, vectors) = symmetric_eigen(&sym)?;
    let largest = values.iter().fold(0.0f64, |s, v| s.max(v.abs()));
    let mut inv = Dense::zeros(p, p);
    let mut rank = 0;
    for (k, &lambda) in values.iter().enumerate() {
        if lambda.abs() <= tol * largest || lambda == 0.0 {
            continue;
        }
        rank += 1;
        let w = 1.0 / lambda;
        for i in 0..p {
            let vi = vectors.get(i, k);
            if vi == 0.0 {
                continue;
            }
            for j in 0..p {
                inv.set(i, j, inv.get(i, j) + w * vi * vectors.get(j, k));
            }
        }
    }
    Ok((inv, rank))
}

/// Coarse phase timers for the benchmark harness
/// (`examples/bench_models.rs`).
///
/// Long-running algorithms open a [`phase::start`] timer around their main
/// phases (for example local moves and aggregation in Louvain, or drawing
/// permutations and computing statistics in QAP). Timing is off unless
/// [`phase::enable`] switched it on, and then costs one relaxed atomic load
/// per timer; when on, each closed timer adds its elapsed time to a global
/// table, summed over threads (so parallel phases report CPU time).
#[doc(hidden)]
pub mod phase {
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::Mutex;
    use std::time::Instant;

    static ENABLED: AtomicBool = AtomicBool::new(false);
    static TOTALS: Mutex<Vec<(&'static str, f64, u64)>> = Mutex::new(Vec::new());

    /// Switches phase timing on or off.
    pub fn enable(on: bool) {
        ENABLED.store(on, Ordering::Relaxed);
    }

    /// Returns and clears the totals: (phase, seconds, number of timers).
    pub fn take() -> Vec<(&'static str, f64, u64)> {
        std::mem::take(&mut *TOTALS.lock().unwrap_or_else(|e| e.into_inner()))
    }

    /// A running timer; it records its phase when dropped.
    pub struct Timer(Option<(&'static str, Instant)>);

    /// Starts timing `name` if timing is enabled.
    #[inline]
    pub fn start(name: &'static str) -> Timer {
        if ENABLED.load(Ordering::Relaxed) {
            Timer(Some((name, Instant::now())))
        } else {
            Timer(None)
        }
    }

    impl Drop for Timer {
        fn drop(&mut self) {
            if let Some((name, t0)) = self.0.take() {
                let dt = t0.elapsed().as_secs_f64();
                let mut tot = TOTALS.lock().unwrap_or_else(|e| e.into_inner());
                match tot.iter_mut().find(|x| x.0 == name) {
                    Some(x) => {
                        x.1 += dt;
                        x.2 += 1;
                    }
                    None => tot.push((name, dt, 1)),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    fn square(data: Vec<f64>) -> Dense {
        let n = (data.len() as f64).sqrt() as usize;
        Dense::from_row_major(n, n, data).unwrap()
    }

    fn assert_identity(m: &Dense, tol: f64) {
        for r in 0..m.rows() {
            for c in 0..m.cols() {
                let want = if r == c { 1.0 } else { 0.0 };
                assert!(
                    (m.get(r, c) - want).abs() < tol,
                    "({r}, {c}) = {}",
                    m.get(r, c)
                );
            }
        }
    }

    #[test]
    fn products_and_transpose() {
        let a = Dense::from_row_major(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        let t = transpose(&a);
        assert_eq!((t.rows(), t.cols()), (3, 2));
        assert_eq!(t.row(1), &[2.0, 5.0]);
        let ata = matmul(&t, &a);
        assert_eq!(ata.row(0), &[17.0, 22.0, 27.0]);
        assert_eq!(matvec(&a, &[1.0, 0.0, -1.0]), vec![-2.0, -2.0]);
        assert_eq!(
            quadratic_form(&square(vec![2.0, 0.0, 0.0, 3.0]), &[1.0, 2.0]),
            14.0
        );
    }

    #[test]
    fn cholesky_solves_and_inverts_spd_matrices() {
        let a = square(vec![4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0]);
        let l = cholesky(&a).unwrap();
        let back = matmul(&l, &transpose(&l));
        for r in 0..3 {
            for c in 0..3 {
                assert!((back.get(r, c) - a.get(r, c)).abs() < 1e-14);
            }
        }
        assert_identity(&matmul(&a, &cholesky_inverse(&a).unwrap()), 1e-12);
        let x = solve_spd(&a, &[1.0, 2.0, 3.0]).unwrap();
        let b = matvec(&a, &x);
        assert!((b[0] - 1.0).abs() < 1e-12 && (b[2] - 3.0).abs() < 1e-12);
        // Indefinite: no factor.
        assert!(cholesky(&square(vec![0.0, 1.0, 1.0, 0.0])).is_none());
        assert!(solve_spd(&square(vec![1.0, 2.0, 2.0, 1.0]), &[1.0, 1.0]).is_none());
    }

    #[test]
    fn gauss_jordan_and_elimination_agree_with_hand_values() {
        let a = square(vec![4.0, 7.0, 2.0, 6.0]);
        let inv = invert(&a, 1e-13).unwrap();
        for (x, y) in inv.as_slice().iter().zip(&[0.6, -0.7, -0.2, 0.4]) {
            assert!((x - y).abs() < 1e-12);
        }
        assert_identity(&matmul(&a, &inv), 1e-12);
        // Indefinite but invertible: the fallback of invert_symmetric.
        let b = square(vec![0.0, 1.0, 1.0, 0.0]);
        let inv = invert_symmetric(&b, 1e-12).unwrap();
        assert!((inv.get(0, 1) - 1.0).abs() < 1e-12 && inv.get(0, 0).abs() < 1e-12);
        let singular = square(vec![1.0, 2.0, 2.0, 4.0]);
        assert!(invert(&singular, 1e-13).is_none());
        assert!(invert_symmetric(&singular, 1e-12).is_none());
        assert!(solve(&singular, &[1.0, 1.0], 1e-12).is_none());
        // Elimination needs the pivot swap on this matrix.
        let c = square(vec![0.0, 2.0, 1.0, 1.0, 3.0, 0.0, 2.0, 0.0, 1.0]);
        let x = solve(&c, &[5.0, 7.0, 4.0], 1e-12).unwrap();
        let back = matvec(&c, &x);
        for (u, v) in back.iter().zip(&[5.0, 7.0, 4.0]) {
            assert!((u - v).abs() < 1e-12);
        }
        assert!(invert(&Dense::zeros(0, 0), 1e-12).is_some());
    }

    fn check_eigen(a: &Dense) {
        let n = a.rows();
        let (vals, vecs) = symmetric_eigen(a).unwrap();
        let scale = a.as_slice().iter().fold(1.0f64, |m, x| m.max(x.abs()));
        assert!(vals.windows(2).all(|w| w[0] <= w[1]));
        for (j, &lambda) in vals.iter().enumerate() {
            for i in 0..n {
                let av: f64 = (0..n).map(|k| a.get(i, k) * vecs.get(k, j)).sum();
                assert!(
                    (av - lambda * vecs.get(i, j)).abs() < 1e-10 * scale * n as f64,
                    "residual at ({i}, {j})"
                );
            }
            for k in 0..n {
                let dot: f64 = (0..n).map(|i| vecs.get(i, j) * vecs.get(i, k)).sum();
                let expect = if j == k { 1.0 } else { 0.0 };
                assert!((dot - expect).abs() < 1e-10, "orthogonality ({j}, {k})");
            }
        }
    }

    #[test]
    fn eigen_decomposition_of_random_symmetric_matrices() {
        let mut rng = Stream::new(42, 0);
        for &n in &[1usize, 2, 3, 5, 10, 31] {
            let mut a = Dense::zeros(n, n);
            for i in 0..n {
                for j in 0..=i {
                    let x = rng.uniform() * 2.0 - 1.0;
                    a.set(i, j, x);
                    a.set(j, i, x);
                }
            }
            check_eigen(&a);
        }
        for seed in 0..6 {
            let mut rng = Stream::new(seed, 1);
            let n = 2 + seed as usize * 3;
            let mut a = Dense::zeros(n, n);
            for i in 0..n {
                for j in 0..=i {
                    let x = rng.below(2) as f64;
                    a.set(i, j, x);
                    a.set(j, i, x);
                }
            }
            check_eigen(&a);
        }
        // Repeated eigenvalues (complete graph), zero and empty matrices.
        let n = 6;
        let k: Vec<f64> = (0..n * n)
            .map(|x| if x / n == x % n { 0.0 } else { 1.0 })
            .collect();
        check_eigen(&square(k));
        check_eigen(&Dense::zeros(4, 4));
        let (vals, vecs) = symmetric_eigen(&Dense::zeros(0, 0)).unwrap();
        assert!(vals.is_empty() && vecs.rows() == 0);
    }

    #[test]
    fn parallel_spd_inverse_matches_cholesky_inverse() {
        let mut rng = crate::rng::Stream::new(4, 0);
        for &n in &[1usize, 2, 7, 40, 97] {
            let b: Vec<f64> = (0..n * n).map(|_| rng.uniform() - 0.5).collect();
            let b = square(b);
            let mut a = matmul(&b, &transpose(&b));
            for i in 0..n {
                a.set(i, i, a.get(i, i) + 0.5);
            }
            let want = cholesky_inverse(&a).unwrap();
            let one = spd_inverse_parallel(&a, 1).unwrap();
            let four = spd_inverse_parallel(&a, 4).unwrap();
            assert_eq!(one, four);
            for i in 0..n {
                for j in 0..n {
                    assert!(
                        (one.get(i, j) - want.get(i, j)).abs()
                            < 1e-9 * (1.0 + want.get(i, j).abs())
                    );
                }
            }
        }
        let mut not_pd = Dense::identity(3);
        not_pd.set(1, 1, -1.0);
        assert!(spd_inverse_parallel(&not_pd, 2).is_none());
    }

    #[test]
    fn conjugate_gradient_solves_laplacian_systems() {
        // Laplacian of the path 0-1-2-3 with conductances 1, 2, 4; inject
        // 1 at node 0 and extract it at node 3: potential drops 1, 1/2, 1/4.
        let c = [1.0, 2.0, 4.0];
        let apply = |v: &[f64], out: &mut [f64]| {
            out.iter_mut().for_each(|x| *x = 0.0);
            for (i, &w) in c.iter().enumerate() {
                let d = w * (v[i] - v[i + 1]);
                out[i] += d;
                out[i + 1] -= d;
            }
        };
        let diag = [1.0, 3.0, 6.0, 4.0];
        let s = conjugate_gradient(apply, &diag, &[1.0, 0.0, 0.0, -1.0], 1e-14, 100);
        assert!(s.residual <= 1e-14);
        let x = &s.x;
        assert!((x[0] - x[1] - 1.0).abs() < 1e-12);
        assert!((x[1] - x[2] - 0.5).abs() < 1e-12);
        assert!((x[2] - x[3] - 0.25).abs() < 1e-12);
    }

    #[test]
    fn known_spectra() {
        let (vals, _) = symmetric_eigen(&square(vec![2.0, 1.0, 1.0, 2.0])).unwrap();
        assert!((vals[0] - 1.0).abs() < 1e-14 && (vals[1] - 3.0).abs() < 1e-14);
        // Path on 4 nodes: 2 cos(k pi / 5).
        let p = square(vec![
            0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
        ]);
        let (vals, _) = symmetric_eigen(&p).unwrap();
        for (k, v) in vals.iter().enumerate() {
            let want = 2.0 * ((4 - k) as f64 * std::f64::consts::PI / 5.0).cos();
            assert!((v - want).abs() < 1e-13);
        }
        let (vals, vecs) = symmetric_eigen(&square(vec![5.0])).unwrap();
        assert_eq!((vals[0], vecs.get(0, 0)), (5.0, 1.0));
    }

    #[test]
    fn pseudo_inverse_of_a_rank_one_matrix() {
        // v v' with v = (1, 2): the pseudo-inverse is v v' / |v|^4.
        let s = square(vec![1.0, 2.0, 2.0, 4.0]);
        let (pinv, rank) = pseudo_inverse(&s, 1e-10).unwrap();
        assert_eq!(rank, 1);
        for (x, y) in pinv
            .as_slice()
            .iter()
            .zip(&[1.0 / 25.0, 2.0 / 25.0, 2.0 / 25.0, 4.0 / 25.0])
        {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
        // Full rank: agrees with the inverse.
        let a = square(vec![4.0, 1.0, 1.0, 3.0]);
        let (pinv, rank) = pseudo_inverse(&a, 1e-10).unwrap();
        assert_eq!(rank, 2);
        assert_identity(&matmul(&a, &pinv), 1e-12);
    }
}
