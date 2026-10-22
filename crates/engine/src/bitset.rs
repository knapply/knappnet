//! Dense adjacency as a bit matrix, for binary (unweighted) networks.
//!
//! An unweighted simple network needs only one bit per ordered pair to say
//! whether an arc exists. When at least 1 in 32 pairs is connected, a bit
//! matrix takes no more memory than 32-bit adjacency lists and answers
//! "is there an arc from a to b?" with a single load instead of a search.
//! Algorithms that test membership heavily (triad census, cliques,
//! clustering) use it for dense networks.

/// A square bit matrix with rows padded to whole 64-bit words.
#[derive(Debug, Clone)]
pub struct BitMatrix {
    words_per_row: usize,
    bits: Vec<u64>,
}

/// Largest bit matrix built automatically, in bytes.
pub const MAX_BYTES: usize = 256 << 20;

impl BitMatrix {
    /// An empty n-by-n matrix.
    pub fn new(n: usize) -> BitMatrix {
        let words_per_row = n.saturating_add(63) / 64;
        BitMatrix {
            words_per_row,
            bits: vec![0; words_per_row * n],
        }
    }

    /// Whether a bit matrix is the better representation for `n` nodes and
    /// `arcs` distinct arcs: at least as compact as adjacency lists, and no
    /// larger than [`MAX_BYTES`].
    pub fn worthwhile(n: usize, arcs: usize) -> bool {
        let bytes = (n.saturating_add(63) / 64)
            .saturating_mul(8)
            .saturating_mul(n);
        bytes <= MAX_BYTES && bytes <= arcs.saturating_mul(4)
    }

    /// Sets the bit for (row, col).
    #[inline]
    pub fn set(&mut self, row: usize, col: usize) {
        self.bits[row * self.words_per_row + col / 64] |= 1u64 << (col % 64);
    }

    /// Reads the bit for (row, col).
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> bool {
        self.bits[row * self.words_per_row + col / 64] >> (col % 64) & 1 == 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_and_get() {
        let mut m = BitMatrix::new(130);
        m.set(0, 129);
        m.set(129, 0);
        m.set(64, 64);
        assert!(m.get(0, 129) && m.get(129, 0) && m.get(64, 64));
        assert!(!m.get(0, 128) && !m.get(1, 129));
        assert!(BitMatrix::worthwhile(100, 1000));
        assert!(!BitMatrix::worthwhile(100_000, 1_000_000));
    }
}
