//! Small dense linear algebra for the estimation algorithm: inversion by
//! Gauss-Jordan elimination with partial pivoting, and the Moore-Penrose
//! inverse of a symmetric matrix through its eigendecomposition.
//!
//! Matrices are row-major `p` by `p` slices. These routines are for the
//! handful of parameters of a SAOM, not for large systems.

use crate::error::{Error, Result};
use crate::metrics::scaling::symmetric_eigen;

/// The inverse of a square row-major matrix, or an error when it is
/// singular to working precision.
pub fn invert(a: &[f64], p: usize) -> Result<Vec<f64>> {
    debug_assert_eq!(a.len(), p * p);
    let mut m = a.to_vec();
    let mut inv = vec![0.0; p * p];
    for i in 0..p {
        inv[i * p + i] = 1.0;
    }
    let scale = a
        .iter()
        .fold(0.0f64, |s, v| s.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    for c in 0..p {
        // Partial pivoting.
        let mut pivot = c;
        for r in c + 1..p {
            if m[r * p + c].abs() > m[pivot * p + c].abs() {
                pivot = r;
            }
        }
        let v = m[pivot * p + c];
        if v.abs() <= 1e-13 * scale {
            return Err(Error::InvalidArgument(
                "a matrix in the estimation algorithm is singular (the derivative or \
                 covariance matrix of the statistics); some effects may be perfectly \
                 collinear or have statistics that never vary"
                    .into(),
            ));
        }
        if pivot != c {
            for j in 0..p {
                m.swap(c * p + j, pivot * p + j);
                inv.swap(c * p + j, pivot * p + j);
            }
        }
        let inv_v = 1.0 / v;
        for j in 0..p {
            m[c * p + j] *= inv_v;
            inv[c * p + j] *= inv_v;
        }
        for r in 0..p {
            if r == c {
                continue;
            }
            let f = m[r * p + c];
            if f == 0.0 {
                continue;
            }
            for j in 0..p {
                m[r * p + j] -= f * m[c * p + j];
                inv[r * p + j] -= f * inv[c * p + j];
            }
        }
    }
    Ok(inv)
}

/// The Moore-Penrose inverse of a symmetric row-major matrix. Eigenvalues
/// below `tol` times the largest are treated as zero, and the rank used is
/// returned with the inverse.
pub fn pseudo_inverse(a: &[f64], p: usize, tol: f64) -> Result<(Vec<f64>, usize)> {
    debug_assert_eq!(a.len(), p * p);
    let mut sym = a.to_vec();
    for i in 0..p {
        for j in 0..i {
            let v = 0.5 * (sym[i * p + j] + sym[j * p + i]);
            sym[i * p + j] = v;
            sym[j * p + i] = v;
        }
    }
    let (values, vectors) = symmetric_eigen(sym, p)?;
    let largest = values.iter().fold(0.0f64, |s, v| s.max(v.abs()));
    let mut inv = vec![0.0; p * p];
    let mut rank = 0;
    for (k, &lambda) in values.iter().enumerate() {
        if lambda.abs() <= tol * largest || lambda == 0.0 {
            continue;
        }
        rank += 1;
        let w = 1.0 / lambda;
        for i in 0..p {
            let vi = vectors[i * p + k];
            if vi == 0.0 {
                continue;
            }
            for j in 0..p {
                inv[i * p + j] += w * vi * vectors[j * p + k];
            }
        }
    }
    Ok((inv, rank))
}

/// `a * b` for row-major `p` by `p` matrices.
pub fn multiply(a: &[f64], b: &[f64], p: usize) -> Vec<f64> {
    let mut out = vec![0.0; p * p];
    for i in 0..p {
        for k in 0..p {
            let aik = a[i * p + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..p {
                out[i * p + j] += aik * b[k * p + j];
            }
        }
    }
    out
}

/// `a * v` for a row-major `p` by `p` matrix and a vector.
pub fn multiply_vector(a: &[f64], v: &[f64], p: usize) -> Vec<f64> {
    (0..p)
        .map(|i| (0..p).map(|j| a[i * p + j] * v[j]).sum())
        .collect()
}

/// The transpose of a row-major `p` by `p` matrix.
pub fn transpose(a: &[f64], p: usize) -> Vec<f64> {
    let mut out = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..p {
            out[j * p + i] = a[i * p + j];
        }
    }
    out
}

/// `v' a v` for a row-major `p` by `p` matrix.
pub fn quadratic_form(a: &[f64], v: &[f64], p: usize) -> f64 {
    let av = multiply_vector(a, v, p);
    v.iter().zip(&av).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_and_pseudo_inverse_of_known_matrices() {
        let a = [4.0, 7.0, 2.0, 6.0];
        let inv = invert(&a, 2).unwrap();
        let expect = [0.6, -0.7, -0.2, 0.4];
        for (x, y) in inv.iter().zip(&expect) {
            assert!((x - y).abs() < 1e-12);
        }
        assert!(invert(&[1.0, 2.0, 2.0, 4.0], 2).is_err());
        // Rank-one symmetric matrix: pseudo-inverse is v v' / |v|^4.
        let s = [1.0, 2.0, 2.0, 4.0];
        let (pinv, rank) = pseudo_inverse(&s, 2, 1e-10).unwrap();
        assert_eq!(rank, 1);
        for (x, y) in pinv
            .iter()
            .zip(&[1.0 / 25.0, 2.0 / 25.0, 2.0 / 25.0, 4.0 / 25.0])
        {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
        let prod = multiply(&a, &inv, 2);
        for i in 0..2 {
            for j in 0..2 {
                let e = if i == j { 1.0 } else { 0.0 };
                assert!((prod[i * 2 + j] - e).abs() < 1e-12);
            }
        }
        assert_eq!(
            transpose(&[1.0, 2.0, 3.0, 4.0], 2),
            vec![1.0, 3.0, 2.0, 4.0]
        );
        assert_eq!(quadratic_form(&[2.0, 0.0, 0.0, 3.0], &[1.0, 2.0], 2), 14.0);
    }
}
