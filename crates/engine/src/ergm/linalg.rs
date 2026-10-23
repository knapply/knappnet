//! Small dense linear algebra for the p-by-p systems of estimation: a
//! symmetric positive definite solve by Cholesky factorization (Golub and
//! Van Loan 2013, `golub2013matrix`, Algorithm 4.2.1) with a pivoted
//! Gauss-Jordan fallback for matrices that are not positive definite.

use crate::error::{Error, Result};
use crate::matrix::Dense;

/// Inverts a symmetric matrix. Positive definite matrices use Cholesky;
/// others fall back to Gauss-Jordan elimination with partial pivoting.
/// Fails when the matrix is singular to working precision.
pub fn invert_symmetric(a: &Dense) -> Result<Dense> {
    if let Some(inv) = cholesky_inverse(a) {
        return Ok(inv);
    }
    gauss_jordan_inverse(a)
}

/// The inverse by Cholesky, or `None` if `a` is not positive definite.
pub fn cholesky_inverse(a: &Dense) -> Option<Dense> {
    let p = a.rows();
    let l = cholesky(a)?;
    // Invert L, then A^{-1} = L^{-T} L^{-1}.
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

/// The lower Cholesky factor, or `None` if `a` is not positive definite.
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

/// Solves A x = b for symmetric positive definite A, or `None`.
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

fn gauss_jordan_inverse(a: &Dense) -> Result<Dense> {
    let p = a.rows();
    let mut m = a.clone();
    let mut inv = Dense::zeros(p, p);
    for i in 0..p {
        inv.set(i, i, 1.0);
    }
    let scale = a
        .as_slice()
        .iter()
        .fold(0.0f64, |acc, v| acc.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    for c in 0..p {
        let pivot = (c..p)
            .max_by(|&x, &y| m.get(x, c).abs().partial_cmp(&m.get(y, c).abs()).unwrap())
            .unwrap();
        if m.get(pivot, c).abs() <= 1e-12 * scale {
            return Err(Error::InvalidArgument(
                "the information matrix is singular; the model may have collinear \
                 statistics or a statistic that is constant over the sample"
                    .into(),
            ));
        }
        if pivot != c {
            for k in 0..p {
                let (x, y) = (m.get(c, k), m.get(pivot, k));
                m.set(c, k, y);
                m.set(pivot, k, x);
                let (x, y) = (inv.get(c, k), inv.get(pivot, k));
                inv.set(c, k, y);
                inv.set(pivot, k, x);
            }
        }
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
    Ok(inv)
}

/// y = A x.
pub fn matvec(a: &Dense, x: &[f64]) -> Vec<f64> {
    (0..a.rows())
        .map(|r| a.row(r).iter().zip(x).map(|(u, v)| u * v).sum())
        .collect()
}

/// x' A x.
pub fn quadratic_form(a: &Dense, x: &[f64]) -> f64 {
    matvec(a, x).iter().zip(x).map(|(u, v)| u * v).sum()
}

/// A B.
pub fn matmul(a: &Dense, b: &Dense) -> Dense {
    let mut out = Dense::zeros(a.rows(), b.cols());
    for r in 0..a.rows() {
        for k in 0..a.cols() {
            let v = a.get(r, k);
            if v == 0.0 {
                continue;
            }
            for c in 0..b.cols() {
                out.set(r, c, out.get(r, c) + v * b.get(k, c));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_round_trips() {
        let a =
            Dense::from_row_major(3, 3, vec![4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0]).unwrap();
        let inv = invert_symmetric(&a).unwrap();
        let id = matmul(&a, &inv);
        for r in 0..3 {
            for c in 0..3 {
                let want = if r == c { 1.0 } else { 0.0 };
                assert!((id.get(r, c) - want).abs() < 1e-12);
            }
        }
        // Indefinite but invertible: the fallback path.
        let b = Dense::from_row_major(2, 2, vec![0.0, 1.0, 1.0, 0.0]).unwrap();
        let inv = invert_symmetric(&b).unwrap();
        assert!((inv.get(0, 1) - 1.0).abs() < 1e-12 && inv.get(0, 0).abs() < 1e-12);
        let s = Dense::from_row_major(2, 2, vec![1.0, 2.0, 2.0, 4.0]).unwrap();
        assert!(invert_symmetric(&s).is_err());
        let x = solve_spd(&a, &[1.0, 2.0, 3.0]).unwrap();
        let back = matvec(&a, &x);
        assert!((back[0] - 1.0).abs() < 1e-12 && (back[2] - 3.0).abs() < 1e-12);
    }
}
