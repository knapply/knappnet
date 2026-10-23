//! Disjoint sets (union-find) over node indices.
//!
//! Tarjan 1975, `tarjan1975efficiency`. Union by index (the smaller root
//! absorbs the larger, so the representative of a set is its smallest
//! member ever chosen as a root) and path halving keep trees shallow;
//! roots are therefore deterministic given the sequence of unions, which
//! the component, blockmodel, and fragmentation code relies on when it
//! relabels sets in order of their roots.

/// A forest of disjoint sets over `0..n`.
#[derive(Debug, Clone)]
pub struct UnionFind {
    parent: Vec<u32>,
}

impl UnionFind {
    /// `n` singleton sets.
    pub fn new(n: usize) -> UnionFind {
        UnionFind {
            parent: (0..n as u32).collect(),
        }
    }

    /// Resets to `n` singleton sets, reusing the allocation.
    pub fn reset(&mut self, n: usize) {
        self.parent.clear();
        self.parent.extend(0..n as u32);
    }

    /// Number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    /// Whether there are no elements.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// The representative of `v`'s set, with path halving.
    #[inline]
    pub fn find(&mut self, mut v: usize) -> usize {
        while self.parent[v] as usize != v {
            let p = self.parent[v] as usize;
            self.parent[v] = self.parent[p];
            v = self.parent[v] as usize;
        }
        v
    }

    /// Whether `a` and `b` are in the same set.
    #[inline]
    pub fn same(&mut self, a: usize, b: usize) -> bool {
        self.find(a) == self.find(b)
    }

    /// Merges the sets of `a` and `b`. Returns `(kept, absorbed)`, the
    /// smaller and larger of the two roots, when they were distinct, and
    /// `None` when the elements were already together.
    #[inline]
    pub fn union(&mut self, a: usize, b: usize) -> Option<(usize, usize)> {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return None;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.parent[hi] = lo as u32;
        Some((lo, hi))
    }

    /// Whether `v` is the representative of its set.
    #[inline]
    pub fn is_root(&self, v: usize) -> bool {
        self.parent[v] as usize == v
    }

    /// The representative of every element, fully compressed.
    pub fn roots(&mut self) -> Vec<u32> {
        (0..self.len()).map(|v| self.find(v) as u32).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unions_keep_the_smaller_root() {
        let mut s = UnionFind::new(6);
        assert_eq!(s.union(4, 5), Some((4, 5)));
        assert_eq!(s.union(5, 1), Some((1, 4)));
        assert_eq!(s.union(1, 4), None);
        assert!(s.same(5, 1) && !s.same(0, 1));
        assert!(s.is_root(1) && !s.is_root(4));
        assert_eq!(s.roots(), vec![0, 1, 2, 3, 1, 1]);
        s.reset(2);
        assert_eq!(s.len(), 2);
        assert_eq!(s.roots(), vec![0, 1]);
        assert!(UnionFind::new(0).is_empty());
    }
}
