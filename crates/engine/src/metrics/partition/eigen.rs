//! Eigenvalues and eigenvectors of a dense real symmetric matrix.
//!
//! The matrix is reduced to tridiagonal form by Householder reflections
//! and the tridiagonal matrix is diagonalized by the QL algorithm with
//! implicit Wilkinson shifts, accumulating the transformations (Golub and
//! Van Loan 2013, `golub2013matrix`, sections 8.3 and 8.2; the procedures
//! `tred2` and `tql2` of Bowdler, Martin, Reinsch, and Wilkinson 1968,
//! `bowdler1968ql`, and Martin, Reinsch, and Wilkinson 1968,
//! `martin1968householder`, follow the same plan). Cost O(n^3), memory
//! O(n^2). Used for the leading eigenvector of modularity matrices and the
//! spectral steps of core-periphery fitting.

/// Eigenvalues in increasing order and the matching unit eigenvectors as
/// columns of a row-major n-by-n matrix (`vectors[i * n + j]` is component
/// i of eigenvector j). `a` is row-major and must be symmetric.
pub(crate) fn symmetric_eigen(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(a.len(), n * n);
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut z = a.to_vec();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    tridiagonalize(&mut z, n, &mut d, &mut e);
    ql_implicit(&mut d, &mut e, &mut z, n);
    // Sort eigenpairs by eigenvalue.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| d[i].total_cmp(&d[j]));
    let values: Vec<f64> = order.iter().map(|&i| d[i]).collect();
    let mut vectors = vec![0.0; n * n];
    for (new, &old) in order.iter().enumerate() {
        for r in 0..n {
            vectors[r * n + new] = z[r * n + old];
        }
    }
    (values, vectors)
}

/// Householder reduction of the symmetric matrix in `z` (row-major) to
/// tridiagonal form: on return `d` holds the diagonal, `e[i]` the
/// subdiagonal element (i, i - 1) (with `e[0] = 0`), and `z` the orthogonal
/// matrix Q with A = Q T Q'.
fn tridiagonalize(z: &mut [f64], n: usize, d: &mut [f64], e: &mut [f64]) {
    let at = |i: usize, j: usize| i * n + j;
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        if l > 0 {
            let scale: f64 = (0..=l).map(|k| z[at(i, k)].abs()).sum();
            if scale == 0.0 {
                e[i] = z[at(i, l)];
            } else {
                for k in 0..=l {
                    z[at(i, k)] /= scale;
                    h += z[at(i, k)] * z[at(i, k)];
                }
                let f = z[at(i, l)];
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                z[at(i, l)] = f - g;
                let mut f = 0.0;
                for j in 0..=l {
                    z[at(j, i)] = z[at(i, j)] / h;
                    let mut g = 0.0;
                    for k in 0..=j {
                        g += z[at(j, k)] * z[at(i, k)];
                    }
                    for k in (j + 1)..=l {
                        g += z[at(k, j)] * z[at(i, k)];
                    }
                    e[j] = g / h;
                    f += e[j] * z[at(i, j)];
                }
                let hh = f / (h + h);
                for j in 0..=l {
                    let f = z[at(i, j)];
                    let g = e[j] - hh * f;
                    e[j] = g;
                    for k in 0..=j {
                        z[at(j, k)] -= f * e[k] + g * z[at(i, k)];
                    }
                }
            }
        } else {
            e[i] = z[at(i, l)];
        }
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    // Accumulate the transformations.
    for i in 0..n {
        if d[i] != 0.0 {
            for j in 0..i {
                let mut g = 0.0;
                for k in 0..i {
                    g += z[at(i, k)] * z[at(k, j)];
                }
                for k in 0..i {
                    z[at(k, j)] -= g * z[at(k, i)];
                }
            }
        }
        d[i] = z[at(i, i)];
        z[at(i, i)] = 1.0;
        for j in 0..i {
            z[at(j, i)] = 0.0;
            z[at(i, j)] = 0.0;
        }
    }
}

/// QL iteration with implicit shifts on the tridiagonal matrix (d, e),
/// applying the rotations to the columns of `z`.
fn ql_implicit(d: &mut [f64], e: &mut [f64], z: &mut [f64], n: usize) {
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    for l in 0..n {
        let mut iter = 0;
        loop {
            // Find a negligible subdiagonal element.
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            assert!(iter <= 60, "QL iteration did not converge");
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r.abs() } else { -r.abs() });
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut i = m;
            let mut early = false;
            while i > l {
                i -= 1;
                let f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    early = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                for k in 0..n {
                    let f = z[k * n + i + 1];
                    z[k * n + i + 1] = s * z[k * n + i] + c * f;
                    z[k * n + i] = c * z[k * n + i] - s * f;
                }
            }
            if early {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    fn check(a: &[f64], n: usize) {
        let (vals, vecs) = symmetric_eigen(a, n);
        let scale = a.iter().fold(1.0f64, |m, x| m.max(x.abs()));
        for j in 0..n {
            for i in 0..n {
                let av: f64 = (0..n).map(|k| a[i * n + k] * vecs[k * n + j]).sum();
                assert!(
                    (av - vals[j] * vecs[i * n + j]).abs() < 1e-10 * scale * n as f64,
                    "residual"
                );
            }
            for k in 0..n {
                let dot: f64 = (0..n).map(|i| vecs[i * n + j] * vecs[i * n + k]).sum();
                let want = if j == k { 1.0 } else { 0.0 };
                assert!((dot - want).abs() < 1e-10);
            }
        }
        assert!(vals.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn known_spectra() {
        let (vals, _) = symmetric_eigen(&[2.0, 1.0, 1.0, 2.0], 2);
        assert!((vals[0] - 1.0).abs() < 1e-14 && (vals[1] - 3.0).abs() < 1e-14);
        // Path on 4 nodes: 2 cos(k pi / 5).
        let p = [
            0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
        ];
        let (vals, _) = symmetric_eigen(&p, 4);
        for (k, v) in vals.iter().enumerate() {
            let want = 2.0 * ((4 - k) as f64 * std::f64::consts::PI / 5.0).cos();
            assert!((v - want).abs() < 1e-13);
        }
        let (vals, vecs) = symmetric_eigen(&[5.0], 1);
        assert_eq!((vals[0], vecs[0]), (5.0, 1.0));
    }

    #[test]
    fn random_matrices() {
        for seed in 0..10 {
            let mut rng = Stream::new(seed, 0);
            let n = 2 + seed as usize * 3;
            let mut a = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..=i {
                    let x = if seed % 2 == 0 {
                        rng.uniform() - 0.5
                    } else {
                        (rng.below(2)) as f64
                    };
                    a[i * n + j] = x;
                    a[j * n + i] = x;
                }
            }
            check(&a, n);
        }
        // Repeated eigenvalues (complete graph) and zero matrices.
        let n = 6;
        let k: Vec<f64> = (0..n * n)
            .map(|x| if x / n == x % n { 0.0 } else { 1.0 })
            .collect();
        check(&k, n);
        check(&vec![0.0; 16], 4);
    }
}
