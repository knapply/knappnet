//! Dense binary network state for ministep simulation.
//!
//! SAOM networks are small (tens to a few thousand actors) and change one
//! tie at a time, so the state is a dense 0/1 matrix stored twice, by row
//! and by column, with the out- and in-degrees kept up to date. Two-path
//! counts from or to one actor are then contiguous row scans.

/// A directed binary network on `n` actors without self-ties.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Network {
    n: usize,
    /// `rows[i * n + j]` = x_ij.
    rows: Vec<u8>,
    /// `cols[j * n + i]` = x_ij (the transpose, for column scans).
    cols: Vec<u8>,
    outdeg: Vec<u32>,
    indeg: Vec<u32>,
}

impl Network {
    /// An empty network.
    pub fn empty(n: usize) -> Network {
        Network {
            n,
            rows: vec![0; n * n],
            cols: vec![0; n * n],
            outdeg: vec![0; n],
            indeg: vec![0; n],
        }
    }

    /// A network from a row-major 0/1 adjacency matrix (`adj[i * n + j]`).
    /// Diagonal entries are ignored; any non-zero value is a tie.
    pub fn from_adjacency(n: usize, adj: &[u8]) -> Network {
        let mut net = Network::empty(n);
        for i in 0..n {
            for j in 0..n {
                if i != j && adj[i * n + j] != 0 {
                    net.set(i, j, true);
                }
            }
        }
        net
    }

    /// Number of actors.
    #[inline]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Whether the tie i -> j is present.
    #[inline]
    pub fn has(&self, i: usize, j: usize) -> bool {
        self.rows[i * self.n + j] != 0
    }

    /// Row i of the adjacency matrix (out-ties of i).
    #[inline]
    pub fn row(&self, i: usize) -> &[u8] {
        &self.rows[i * self.n..(i + 1) * self.n]
    }

    /// Column j of the adjacency matrix (in-ties of j), contiguous.
    #[inline]
    pub fn col(&self, j: usize) -> &[u8] {
        &self.cols[j * self.n..(j + 1) * self.n]
    }

    /// Out-degree of i.
    #[inline]
    pub fn outdeg(&self, i: usize) -> u32 {
        self.outdeg[i]
    }

    /// In-degree of i.
    #[inline]
    pub fn indeg(&self, i: usize) -> u32 {
        self.indeg[i]
    }

    /// Number of ties.
    pub fn tie_count(&self) -> usize {
        self.outdeg.iter().map(|&d| d as usize).sum()
    }

    /// Sets x_ij. Self-ties are ignored.
    pub fn set(&mut self, i: usize, j: usize, present: bool) {
        if i == j {
            return;
        }
        let n = self.n;
        let was = self.rows[i * n + j] != 0;
        if was == present {
            return;
        }
        let v = present as u8;
        self.rows[i * n + j] = v;
        self.cols[j * n + i] = v;
        if present {
            self.outdeg[i] += 1;
            self.indeg[j] += 1;
        } else {
            self.outdeg[i] -= 1;
            self.indeg[j] -= 1;
        }
    }

    /// Flips x_ij.
    pub fn toggle(&mut self, i: usize, j: usize) {
        let present = self.has(i, j);
        self.set(i, j, !present);
    }

    /// Number of ordered pairs whose tie differs between the two networks.
    pub fn distance(&self, other: &Network) -> usize {
        self.rows
            .iter()
            .zip(&other.rows)
            .filter(|(a, b)| a != b)
            .count()
    }

    /// The adjacency matrix in row-major order.
    pub fn adjacency(&self) -> &[u8] {
        &self.rows
    }

    /// Out-neighbors of i, ascending.
    pub fn out_neighbors(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.row(i)
            .iter()
            .enumerate()
            .filter(|(_, &x)| x != 0)
            .map(|(j, _)| j)
    }

    /// In-neighbors of j, ascending.
    pub fn in_neighbors(&self, j: usize) -> impl Iterator<Item = usize> + '_ {
        self.col(j)
            .iter()
            .enumerate()
            .filter(|(_, &x)| x != 0)
            .map(|(i, _)| i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degrees_and_transpose_follow_toggles() {
        let mut net = Network::empty(4);
        net.set(0, 1, true);
        net.set(2, 1, true);
        net.set(1, 1, true); // ignored
        assert_eq!(net.outdeg(0), 1);
        assert_eq!(net.indeg(1), 2);
        assert_eq!(net.col(1), &[1, 0, 1, 0]);
        assert_eq!(net.in_neighbors(1).collect::<Vec<_>>(), vec![0, 2]);
        net.toggle(0, 1);
        assert!(!net.has(0, 1));
        assert_eq!(net.indeg(1), 1);
        assert_eq!(net.tie_count(), 1);
        let other = Network::from_adjacency(4, &[0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(other.tie_count(), 1);
        assert_eq!(net.distance(&other), 2);
    }
}
