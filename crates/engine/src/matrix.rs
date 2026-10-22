//! A small dense, row-major matrix of doubles, for dyadic data (QAP),
//! equivalence and correlation matrices, and other n-by-n results.

use crate::error::{Error, Result};

/// A dense row-major matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Dense {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Dense {
    /// A matrix of zeros.
    pub fn zeros(rows: usize, cols: usize) -> Dense {
        Dense {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// Wraps row-major data, checking its length.
    pub fn from_row_major(rows: usize, cols: usize, data: Vec<f64>) -> Result<Dense> {
        if data.len() != rows * cols {
            return Err(Error::LengthMismatch {
                what: "matrix data",
                expected: rows * cols,
                found: data.len(),
            });
        }
        Ok(Dense { rows, cols, data })
    }

    /// Wraps column-major data (as R stores matrices), checking its length.
    pub fn from_col_major(rows: usize, cols: usize, data: &[f64]) -> Result<Dense> {
        if data.len() != rows * cols {
            return Err(Error::LengthMismatch {
                what: "matrix data",
                expected: rows * cols,
                found: data.len(),
            });
        }
        let mut m = Dense::zeros(rows, cols);
        for c in 0..cols {
            for r in 0..rows {
                m.data[r * cols + c] = data[c * rows + r];
            }
        }
        Ok(m)
    }

    /// Copies the matrix out in column-major order.
    pub fn to_col_major(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.data.len()];
        for r in 0..self.rows {
            for c in 0..self.cols {
                out[c * self.rows + r] = self.data[r * self.cols + c];
            }
        }
        out
    }

    /// Number of rows.
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Element (r, c).
    #[inline]
    pub fn get(&self, r: usize, c: usize) -> f64 {
        self.data[r * self.cols + c]
    }

    /// Sets element (r, c).
    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: f64) {
        self.data[r * self.cols + c] = v;
    }

    /// Row `r` as a slice.
    #[inline]
    pub fn row(&self, r: usize) -> &[f64] {
        &self.data[r * self.cols..(r + 1) * self.cols]
    }

    /// The underlying row-major data.
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layouts_round_trip() {
        let m = Dense::from_col_major(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        assert_eq!(m.row(0), &[1.0, 3.0, 5.0]);
        assert_eq!(m.get(1, 2), 6.0);
        assert_eq!(m.to_col_major(), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert!(Dense::from_row_major(2, 2, vec![1.0]).is_err());
    }
}
