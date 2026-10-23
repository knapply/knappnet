//! Distances and similarities between two partitions of the same nodes.
//!
//! With n nodes, n_ij the number of nodes in group i of the first partition
//! and group j of the second, a_i and b_j the group sizes, and entropies
//! and mutual information in natural units,
//!
//! ```text
//! H(A)   = -sum_i (a_i/n) ln(a_i/n)
//! I(A,B) =  sum_ij (n_ij/n) ln(n n_ij / (a_i b_j))
//! ```
//!
//! the measures are:
//! - Variation of information (Meila 2007, `meila2007comparing`):
//!   VI = H(A) + H(B) - 2 I(A,B), a metric on partitions.
//! - Normalized mutual information (Danon, Diaz-Guilera, Duch, and Arenas
//!   2005, `danon2005comparing`): NMI = 2 I(A,B) / (H(A) + H(B)).
//! - Split-join distance (van Dongen 2000, `vandongen2000performance`):
//!   the projection distance of A from B is n - sum_i max_j n_ij, and the
//!   split-join distance is the sum of the two projection distances.
//! - Rand index (Rand 1971, `rand1971objective`): the fraction of the
//!   n(n-1)/2 node pairs that both partitions put together or both put
//!   apart.
//! - Adjusted Rand index (Hubert and Arabie 1985, `hubert1985comparing`):
//!   the Rand index corrected for chance under the hypergeometric model,
//!
//!   ```text
//!   ARI = (S - S_a S_b / C) / ((S_a + S_b)/2 - S_a S_b / C),
//!   ```
//!
//!   with S = sum_ij C(n_ij, 2), S_a = sum_i C(a_i, 2),
//!   S_b = sum_j C(b_j, 2), and C = C(n, 2).
//! - Adjusted mutual information (Vinh, Epps, and Bailey 2010,
//!   `vinh2010information`): AMI = (I - E(I)) / (mean(H(A), H(B)) - E(I)),
//!   with E(I) the expected mutual information under the hypergeometric
//!   model of their equation (24a) and the arithmetic mean of the two
//!   entropies (the normalization that matches NMI above).
//!
//! Conventions (matching igraph's `compare()` for the measures it offers):
//! - Labels are arbitrary integers; only equality matters.
//! - Both partitions must have the same length. With no nodes every
//!   measure is NaN.
//! - NMI is 1 when both partitions have zero entropy (each one group).
//! - The Rand and adjusted Rand indices are undefined for fewer than two
//!   nodes and are NaN there (igraph raises an error).
//! - The adjusted Rand index is 1 when its denominator vanishes, which
//!   happens exactly when both partitions are a single group or both are
//!   all singletons (and so identical). The adjusted mutual information is
//!   1 in the same two cases, and NaN if its denominator vanishes in any
//!   other.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::metrics::community::dense_labels;

/// A measure computed by [`compare_partitions`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompareMethod {
    /// Variation of information (natural logarithms).
    VariationOfInformation,
    /// Normalized mutual information of Danon et al.
    NormalizedMutualInformation,
    /// Split-join distance (sum of the two projection distances).
    SplitJoin,
    /// Rand index.
    Rand,
    /// Adjusted Rand index of Hubert and Arabie.
    AdjustedRand,
    /// Adjusted mutual information of Vinh et al. (arithmetic mean).
    AdjustedMutualInformation,
}

impl CompareMethod {
    /// The method for an integer code: 0 = VI, 1 = NMI, 2 = split-join,
    /// 3 = Rand, 4 = adjusted Rand, 5 = adjusted mutual information.
    pub fn from_code(code: i32) -> Result<CompareMethod> {
        Ok(match code {
            0 => CompareMethod::VariationOfInformation,
            1 => CompareMethod::NormalizedMutualInformation,
            2 => CompareMethod::SplitJoin,
            3 => CompareMethod::Rand,
            4 => CompareMethod::AdjustedRand,
            5 => CompareMethod::AdjustedMutualInformation,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown partition comparison method code {code}"
                )))
            }
        })
    }
}

/// Contingency table of two partitions: the cells n_ij (sparse) and the
/// marginal sizes.
struct Table {
    n: usize,
    cells: Vec<usize>,
    a: Vec<usize>,
    b: Vec<usize>,
    /// For each cell, its row and column.
    index: Vec<(usize, usize)>,
}

fn table(x: &[i64], y: &[i64]) -> Result<Table> {
    if x.len() != y.len() {
        return Err(Error::LengthMismatch {
            what: "second membership",
            expected: x.len(),
            found: y.len(),
        });
    }
    let (xa, ka) = dense_labels(x);
    let (yb, kb) = dense_labels(y);
    let mut a = vec![0usize; ka];
    let mut b = vec![0usize; kb];
    let mut pos: HashMap<(usize, usize), usize> = HashMap::new();
    let mut cells = Vec::new();
    let mut index = Vec::new();
    for (&i, &j) in xa.iter().zip(&yb) {
        a[i] += 1;
        b[j] += 1;
        let next = cells.len();
        let c = *pos.entry((i, j)).or_insert(next);
        if c == cells.len() {
            cells.push(0);
            index.push((i, j));
        }
        cells[c] += 1;
    }
    Ok(Table {
        n: x.len(),
        cells,
        a,
        b,
        index,
    })
}

fn entropy(sizes: &[usize], n: f64) -> f64 {
    -sizes
        .iter()
        .filter(|&&s| s > 0)
        .map(|&s| {
            let p = s as f64 / n;
            p * p.ln()
        })
        .sum::<f64>()
}

fn mutual_information(t: &Table) -> f64 {
    let n = t.n as f64;
    t.cells
        .iter()
        .zip(&t.index)
        .map(|(&c, &(i, j))| {
            let c = c as f64;
            c / n * (n * c / (t.a[i] as f64 * t.b[j] as f64)).ln()
        })
        .sum()
}

fn pairs(k: usize) -> f64 {
    let k = k as f64;
    k * (k - 1.0) / 2.0
}

/// Expected mutual information of two random partitions with the given
/// group sizes (Vinh, Epps, and Bailey 2010, equation 24a).
fn expected_mutual_information(a: &[usize], b: &[usize], n: usize) -> f64 {
    // ln k! for k = 0..=n.
    let mut lf = vec![0.0f64; n + 1];
    for k in 1..=n {
        lf[k] = lf[k - 1] + (k as f64).ln();
    }
    let nf = n as f64;
    let mut e = 0.0;
    for &ai in a {
        for &bj in b {
            let lo = (ai + bj).saturating_sub(n).max(1);
            let hi = ai.min(bj);
            let fixed = lf[ai] + lf[bj] + lf[n - ai] + lf[n - bj] - lf[n];
            for nij in lo..=hi {
                let x = nij as f64;
                let log_p = fixed - lf[nij] - lf[ai - nij] - lf[bj - nij] - lf[n + nij - ai - bj];
                e += x / nf * (nf * x / (ai as f64 * bj as f64)).ln() * log_p.exp();
            }
        }
    }
    e
}

/// Both partitions are one group, or both are all singletons.
fn trivially_identical(t: &Table) -> bool {
    (t.a.len() == 1 && t.b.len() == 1) || (t.a.len() == t.n && t.b.len() == t.n)
}

/// The two projection distances of van Dongen (2000): of the first
/// partition from the second and of the second from the first.
pub fn split_join_distance(x: &[i64], y: &[i64]) -> Result<(usize, usize)> {
    let t = table(x, y)?;
    let mut row_max = vec![0usize; t.a.len()];
    let mut col_max = vec![0usize; t.b.len()];
    for (&c, &(i, j)) in t.cells.iter().zip(&t.index) {
        row_max[i] = row_max[i].max(c);
        col_max[j] = col_max[j].max(c);
    }
    Ok((
        t.n - row_max.iter().sum::<usize>(),
        t.n - col_max.iter().sum::<usize>(),
    ))
}

/// Compares two partitions (one label per node in each). See the module
/// documentation for the measures and conventions.
pub fn compare_partitions(x: &[i64], y: &[i64], method: CompareMethod) -> Result<f64> {
    let t = table(x, y)?;
    if t.n == 0 {
        return Ok(f64::NAN);
    }
    let n = t.n as f64;
    Ok(match method {
        CompareMethod::VariationOfInformation => {
            let v = entropy(&t.a, n) + entropy(&t.b, n) - 2.0 * mutual_information(&t);
            v.max(0.0)
        }
        CompareMethod::NormalizedMutualInformation => {
            let h = entropy(&t.a, n) + entropy(&t.b, n);
            if h == 0.0 {
                1.0
            } else {
                (2.0 * mutual_information(&t) / h).clamp(0.0, 1.0)
            }
        }
        CompareMethod::SplitJoin => {
            let (p, q) = split_join_distance(x, y)?;
            (p + q) as f64
        }
        CompareMethod::Rand | CompareMethod::AdjustedRand => {
            if t.n < 2 {
                return Ok(f64::NAN);
            }
            let s: f64 = t.cells.iter().map(|&c| pairs(c)).sum();
            let sa: f64 = t.a.iter().map(|&c| pairs(c)).sum();
            let sb: f64 = t.b.iter().map(|&c| pairs(c)).sum();
            let c = pairs(t.n);
            if method == CompareMethod::Rand {
                // Pairs together in both plus pairs apart in both.
                (c + 2.0 * s - sa - sb) / c
            } else {
                let expected = sa * sb / c;
                let denom = (sa + sb) / 2.0 - expected;
                if denom == 0.0 {
                    1.0
                } else {
                    (s - expected) / denom
                }
            }
        }
        CompareMethod::AdjustedMutualInformation => {
            if trivially_identical(&t) {
                return Ok(1.0);
            }
            let i = mutual_information(&t);
            let e = expected_mutual_information(&t.a, &t.b, t.n);
            let mean_h = (entropy(&t.a, n) + entropy(&t.b, n)) / 2.0;
            let denom = mean_h - e;
            if denom == 0.0 {
                f64::NAN
            } else {
                (i - e) / denom
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use CompareMethod::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn identical_and_relabelled() {
        let x = [0, 0, 1, 1, 2];
        let y = [7, 7, 3, 3, 5];
        assert!(close(
            compare_partitions(&x, &y, VariationOfInformation).unwrap(),
            0.0
        ));
        assert!(close(
            compare_partitions(&x, &y, NormalizedMutualInformation).unwrap(),
            1.0
        ));
        assert_eq!(compare_partitions(&x, &y, SplitJoin).unwrap(), 0.0);
        assert!(close(compare_partitions(&x, &y, Rand).unwrap(), 1.0));
        assert!(close(
            compare_partitions(&x, &y, AdjustedRand).unwrap(),
            1.0
        ));
        assert!(close(
            compare_partitions(&x, &y, AdjustedMutualInformation).unwrap(),
            1.0
        ));
    }

    #[test]
    fn hand_values() {
        // {0,1,2}{3} against {0,1}{2,3}: n_ij = [[2,1],[0,1]].
        let x = [0, 0, 0, 1];
        let y = [0, 0, 1, 1];
        // Pairs: together in both {0,1}; apart in both {0,3},{1,3}: 3 of 6.
        assert!(close(compare_partitions(&x, &y, Rand).unwrap(), 0.5));
        // S = 1, S_a = 3, S_b = 2, C = 6: (1 - 1) / (2.5 - 1) = 0.
        assert!(close(
            compare_partitions(&x, &y, AdjustedRand).unwrap(),
            0.0
        ));
        // Projection distances: 4 - (2 + 1) and 4 - (2 + 1).
        assert_eq!(split_join_distance(&x, &y).unwrap(), (1, 1));
        let ha = -(0.75f64 * 0.75f64.ln() + 0.25 * 0.25f64.ln());
        let hb = 2f64.ln();
        let i = 0.5 * (4.0 * 2.0 / 6.0f64).ln()
            + 0.25 * (4.0 * 1.0 / 6.0f64).ln()
            + 0.25 * (4.0 * 1.0 / 2.0f64).ln();
        let vi = compare_partitions(&x, &y, VariationOfInformation).unwrap();
        assert!(close(vi, ha + hb - 2.0 * i));
        let nmi = compare_partitions(&x, &y, NormalizedMutualInformation).unwrap();
        assert!(close(nmi, 2.0 * i / (ha + hb)));
    }

    #[test]
    fn edge_cases() {
        assert!(compare_partitions(&[], &[], Rand).unwrap().is_nan());
        assert!(compare_partitions(&[1], &[2], Rand).unwrap().is_nan());
        assert_eq!(
            compare_partitions(&[1], &[2], NormalizedMutualInformation).unwrap(),
            1.0
        );
        assert!(compare_partitions(&[0, 1], &[0], Rand).is_err());
        // Both one group: ARI and AMI are 1 by convention.
        let one = [3, 3, 3];
        assert_eq!(compare_partitions(&one, &one, AdjustedRand).unwrap(), 1.0);
        assert_eq!(
            compare_partitions(&one, &one, AdjustedMutualInformation).unwrap(),
            1.0
        );
        // One group against all singletons: NMI 0, ARI 0.
        let single = [0, 1, 2];
        assert!(close(
            compare_partitions(&one, &single, NormalizedMutualInformation).unwrap(),
            0.0
        ));
        assert!(close(
            compare_partitions(&one, &single, AdjustedRand).unwrap(),
            0.0
        ));
    }

    #[test]
    fn expected_mutual_information_by_enumeration() {
        // E[I] over every labelling with sizes (2, 2) and (3, 1) of 4 nodes:
        // fix the first partition {0,1}{2,3}; the singleton of the second
        // is equally likely to be any node, and every choice gives the
        // same table [[2,0],[1,1]] up to symmetry.
        let i = 0.5 * (4.0 * 2.0 / 6.0f64).ln()
            + 0.25 * (4.0 * 1.0 / 6.0f64).ln()
            + 0.25 * (4.0 * 1.0 / 2.0f64).ln();
        let e = expected_mutual_information(&[2, 2], &[3, 1], 4);
        assert!(close(e, i));
    }
}
