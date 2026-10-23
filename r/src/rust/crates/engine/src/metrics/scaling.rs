//! Classical (metric) multidimensional scaling and stress.
//!
//! Classical scaling (Torgerson 1952, `torgerson1952multidimensional`;
//! Gower 1966, `gower1966some`) places n points in k dimensions so that
//! their Euclidean distances approximate a matrix of dissimilarities D. The
//! squared dissimilarities are double-centred,
//! B = -1/2 J D^(2) J with J = I - 11'/n,
//! B is decomposed as V L V' with eigenvalues l_1 >= l_2 >= ... >= l_n, and
//! the coordinates on dimension j are v_j sqrt(l_j). When D holds the
//! Euclidean distances of some configuration, B is its matrix of centred
//! inner products and the configuration is recovered exactly (up to
//! rotation, reflection, and translation) from the positive eigenvalues.
//!
//! The eigendecomposition is [`crate::linalg::symmetric_eigen`]
//! (Householder reduction and implicit QL, after EISPACK's `tred2` and
//! `tql2`).
//!
//! Conventions, matching R's `stats::cmdscale()`:
//! - D must be square, symmetric, finite, non-negative, with a zero
//!   diagonal. Any dissimilarity can be scaled, Euclidean or not; for
//!   non-Euclidean D some eigenvalues are negative.
//! - k is between 1 and n - 1. Only dimensions with a positive eigenvalue
//!   get coordinates, so fewer than k columns can be returned. An
//!   eigenvalue counts as positive when it exceeds 1e-12 times the largest
//!   one: rounding leaves values of that order where the exact eigenvalue
//!   is zero (cmdscale keeps any value above 0, so it can return extra
//!   dimensions made of rounding error).
//! - All n eigenvalues are reported, in decreasing order.
//! - Goodness of fit is reported as cmdscale's pair (g_1, g_2), with
//!   g_i = sum_{j <= k} l_j / sum_j T_i(l_j), T_1(l) = |l|, T_2(l) =
//!   max(l, 0).
//! - The sign of each dimension is arbitrary in the method; here it is
//!   fixed so that the first coordinate whose magnitude exceeds 1e-8 times
//!   the largest magnitude in that dimension is positive, so results are
//!   reproducible. Dimensions with equal eigenvalues are not unique: any
//!   rotation within them fits equally well.
//! - Stress is Kruskal's (1964, `kruskal1964multidimensional`) stress
//!   formula 1, sqrt(sum (d_ij - delta_ij)^2 / sum d_ij^2), over pairs i < j,
//!   where d are the distances between the points and delta the
//!   dissimilarities themselves (metric stress: the disparities of
//!   Kruskal's monotone regression are replaced by the data).

use crate::error::{Error, Result};
use crate::linalg::symmetric_eigen;
use crate::matrix::Dense;

/// Result of [`classical_mds`].
#[derive(Debug, Clone, PartialEq)]
pub struct Mds {
    /// Coordinates, one row per point and one column per dimension with a
    /// positive eigenvalue among the first k.
    pub points: Dense,
    /// All n eigenvalues of the double-centred matrix, in decreasing order.
    pub eigenvalues: Vec<f64>,
    /// cmdscale's goodness of fit (g_1, g_2).
    pub gof: [f64; 2],
}

/// Classical multidimensional scaling of the dissimilarity matrix `d` into
/// at most `k` dimensions. See the module documentation.
pub fn classical_mds(d: &Dense, k: usize) -> Result<Mds> {
    let n = d.rows();
    if d.cols() != n {
        return Err(Error::InvalidArgument(format!(
            "the dissimilarity matrix must be square, not {} by {}",
            n,
            d.cols()
        )));
    }
    if n < 2 {
        return Err(Error::InvalidArgument(
            "classical scaling needs at least two points".into(),
        ));
    }
    if k < 1 || k > n - 1 {
        return Err(Error::InvalidArgument(format!(
            "the number of dimensions must be between 1 and {} (one less than the number of points), not {k}",
            n - 1
        )));
    }
    check_dissimilarities(d)?;

    // Double centring of the squared dissimilarities.
    let sq: Vec<f64> = d.as_slice().iter().map(|x| x * x).collect();
    let row_mean: Vec<f64> = (0..n)
        .map(|i| sq[i * n..(i + 1) * n].iter().sum::<f64>() / n as f64)
        .collect();
    let grand = row_mean.iter().sum::<f64>() / n as f64;
    let mut b = Dense::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            b.set(
                i,
                j,
                -0.5 * (sq[i * n + j] - row_mean[i] - row_mean[j] + grand),
            );
        }
    }
    let (values, vectors) = symmetric_eigen(&b)?;
    // symmetric_eigen sorts ascending; reverse to decreasing order.
    let order: Vec<usize> = (0..n).rev().collect();
    let eigenvalues: Vec<f64> = order.iter().map(|&j| values[j]).collect();
    // Rounding leaves eigenvalues of order n * eps * l_1 where the exact
    // ones are zero; such dimensions carry no information.
    let floor = 1e-12 * eigenvalues[0].abs().max(f64::MIN_POSITIVE);
    let dims = eigenvalues[..k].iter().take_while(|&&l| l > floor).count();
    let mut points = Dense::zeros(n, dims);
    for (c, &j) in order.iter().take(dims).enumerate() {
        let scale = values[j].sqrt();
        let col: Vec<f64> = (0..n).map(|i| vectors.get(i, j) * scale).collect();
        let big = col.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let flip = col
            .iter()
            .find(|x| x.abs() > 1e-8 * big)
            .is_some_and(|&x| x < 0.0);
        for (i, &x) in col.iter().enumerate() {
            points.set(i, c, if flip { -x } else { x });
        }
    }
    let top: f64 = eigenvalues[..k].iter().sum();
    let abs_total: f64 = eigenvalues.iter().map(|l| l.abs()).sum();
    let pos_total: f64 = eigenvalues.iter().map(|l| l.max(0.0)).sum();
    Ok(Mds {
        points,
        eigenvalues,
        gof: [top / abs_total, top / pos_total],
    })
}

/// Kruskal's stress formula 1 of a configuration (one row per point)
/// against the dissimilarities `d`. See the module documentation.
pub fn stress(d: &Dense, points: &Dense) -> Result<f64> {
    let n = d.rows();
    if d.cols() != n || points.rows() != n {
        return Err(Error::InvalidArgument(format!(
            "the configuration has {} points but the dissimilarity matrix is {} by {}",
            points.rows(),
            n,
            d.cols()
        )));
    }
    check_dissimilarities(d)?;
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..n {
        for j in i + 1..n {
            let dist = points
                .row(i)
                .iter()
                .zip(points.row(j))
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            let diff = dist - d.get(i, j);
            num += diff * diff;
            den += dist * dist;
        }
    }
    Ok(if den > 0.0 {
        (num / den).sqrt()
    } else {
        f64::NAN
    })
}

fn check_dissimilarities(d: &Dense) -> Result<()> {
    let n = d.rows();
    let scale = d.as_slice().iter().fold(0.0f64, |m, x| m.max(x.abs()));
    for i in 0..n {
        for j in 0..n {
            let x = d.get(i, j);
            if !x.is_finite() {
                return Err(Error::InvalidArgument(format!(
                    "dissimilarity ({}, {}) is missing or not finite",
                    i + 1,
                    j + 1
                )));
            }
            if x < 0.0 {
                return Err(Error::InvalidArgument(format!(
                    "dissimilarity ({}, {}) is negative ({x})",
                    i + 1,
                    j + 1
                )));
            }
            if i == j && x != 0.0 {
                return Err(Error::InvalidArgument(format!(
                    "the diagonal of the dissimilarity matrix must be 0, but entry ({}, {}) is {x}",
                    i + 1,
                    i + 1
                )));
            }
            if (x - d.get(j, i)).abs() > 1e-12 * scale {
                return Err(Error::InvalidArgument(format!(
                    "the dissimilarity matrix is not symmetric: entries ({}, {}) and ({}, {}) differ",
                    i + 1,
                    j + 1,
                    j + 1,
                    i + 1
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn euclid(pts: &[[f64; 2]]) -> Dense {
        let n = pts.len();
        let mut d = Dense::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let dx = pts[i][0] - pts[j][0];
                let dy = pts[i][1] - pts[j][1];
                d.set(i, j, (dx * dx + dy * dy).sqrt());
            }
        }
        d
    }

    #[test]
    fn recovers_a_planar_configuration() {
        let pts = [[0.0, 0.0], [3.0, 0.0], [3.0, 1.0], [0.5, 2.0], [-1.0, 1.0]];
        let d = euclid(&pts);
        let m = classical_mds(&d, 2).unwrap();
        assert_eq!(m.points.cols(), 2);
        // Distances are reproduced exactly; the remaining eigenvalues vanish.
        let back = Dense::from_row_major(5, 2, m.points.as_slice().to_vec()).unwrap();
        let rows: Vec<[f64; 2]> = (0..5).map(|i| [back.get(i, 0), back.get(i, 1)]).collect();
        let d2 = euclid(&rows);
        for i in 0..5 {
            for j in 0..5 {
                assert!((d.get(i, j) - d2.get(i, j)).abs() < 1e-10);
            }
        }
        assert!(m.eigenvalues[2..].iter().all(|l| l.abs() < 1e-10));
        assert!((m.gof[0] - 1.0).abs() < 1e-10 && (m.gof[1] - 1.0).abs() < 1e-10);
        assert!(stress(&d, &m.points).unwrap() < 1e-10);
        // Coordinates are centred.
        for c in 0..2 {
            let s: f64 = (0..5).map(|i| m.points.get(i, c)).sum();
            assert!(s.abs() < 1e-10);
        }
        // One dimension: fewer columns than points, GOF below 1.
        let m1 = classical_mds(&d, 1).unwrap();
        assert_eq!(m1.points.cols(), 1);
        assert!(m1.gof[0] < 1.0);
        assert!(stress(&d, &m1.points).unwrap() > 0.0);
    }

    #[test]
    fn rejects_bad_input_and_drops_non_positive_dimensions() {
        let d = Dense::from_row_major(2, 2, vec![0.0, 1.0, 2.0, 0.0]).unwrap();
        assert!(classical_mds(&d, 1).is_err());
        let d = Dense::from_row_major(2, 2, vec![1.0, 1.0, 1.0, 0.0]).unwrap();
        assert!(classical_mds(&d, 1).is_err());
        let d = Dense::from_row_major(2, 2, vec![0.0, 1.0, 1.0, 0.0]).unwrap();
        assert!(classical_mds(&d, 2).is_err());
        let m = classical_mds(&d, 1).unwrap();
        assert!((m.points.get(0, 0).abs() - 0.5).abs() < 1e-12);
        assert!(m.points.get(0, 0) > 0.0);
        // Three collinear points asked for two dimensions: one is returned.
        let pts = [[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]];
        let m = classical_mds(&euclid(&pts), 2).unwrap();
        assert_eq!(m.points.cols(), 1);
    }
}
