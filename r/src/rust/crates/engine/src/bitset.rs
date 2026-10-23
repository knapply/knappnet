//! Bit-parallel structures for binary (unweighted) networks: dense
//! adjacency as a bit matrix, and breadth-first search from many sources at
//! once.
//!
//! An unweighted simple network needs only one bit per ordered pair to say
//! whether an arc exists. When at least 1 in 32 pairs is connected, a bit
//! matrix takes no more memory than 32-bit adjacency lists and answers
//! "is there an arc from a to b?" with a single load instead of a search.
//! Algorithms that test membership heavily (triad census, cliques,
//! clustering) use it for dense networks.
//!
//! Breadth-first search from every node (the all-pairs measures: closeness,
//! harmonic centrality, eccentricity, mean distance, reach) repeats the
//! same work n times. [`MultiBfs`] runs up to 64 searches together, one
//! bit per source: a word per node records which sources have reached it,
//! and one pass over the frontier's arcs advances every search by one
//! level (Then et al. 2014, `then2014more`). Each arc is then scanned once
//! per level in which its tail is on some search's frontier, instead of
//! once per search, and the per-arc work is a word-wide OR. The result of
//! every search is the same as a plain breadth-first search.
//!
//! Searches follow a [`View`] of a graph: its out-arcs, its in-arcs, or
//! both (direction ignored), read from the graph's own adjacency lists, so
//! no reversed or symmetrized copy of the graph is built.

use crate::graph::{Csr, Graph};
use crate::metrics::degree::DegreeMode;

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

    /// The words of row `row`.
    #[inline]
    pub fn row(&self, row: usize) -> &[u64] {
        &self.bits[row * self.words_per_row..(row + 1) * self.words_per_row]
    }
}

/// The arcs a search follows: one or two adjacency lists per node, read in
/// place from a [`Graph`]. `push` lists the arcs leaving a node in the
/// chosen direction, `pull` the arcs entering it (the same lists when
/// direction is ignored).
#[derive(Debug, Clone, Copy)]
pub struct View<'a> {
    push: (&'a Csr, Option<&'a Csr>),
    pull: (&'a Csr, Option<&'a Csr>),
}

impl<'a> View<'a> {
    /// The arcs of `g` that paths in `mode` follow: out-arcs (`Out`),
    /// in-arcs (`In`), or both (`All`). Undirected graphs ignore `mode`.
    pub fn new(g: &'a Graph, mode: DegreeMode) -> View<'a> {
        let (out, inc) = (g.out(), g.inc());
        if !g.is_directed() {
            return View {
                push: (out, None),
                pull: (out, None),
            };
        }
        match mode {
            DegreeMode::Out => View {
                push: (out, None),
                pull: (inc, None),
            },
            DegreeMode::In => View {
                push: (inc, None),
                pull: (out, None),
            },
            DegreeMode::All => View {
                push: (out, Some(inc)),
                pull: (out, Some(inc)),
            },
        }
    }

    /// Number of nodes.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.push.0.node_count()
    }

    /// Number of arcs followed, counting an arc once per list it is in.
    #[inline]
    pub fn arc_count(&self) -> usize {
        self.push.0.arc_count() + self.push.1.map_or(0, |c| c.arc_count())
    }

    /// Number of arcs leaving `v`.
    #[inline]
    pub fn degree(&self, v: usize) -> usize {
        self.push.0.degree(v) + self.push.1.map_or(0, |c| c.degree(v))
    }

    /// Heads of the arcs leaving `v`: one or two slices.
    #[inline]
    pub fn heads(&self, v: usize) -> (&'a [u32], &'a [u32]) {
        (
            self.push.0.neighbors(v),
            self.push.1.map_or(&[][..], |c| c.neighbors(v)),
        )
    }

    /// Tails of the arcs entering `v`: one or two slices.
    #[inline]
    pub fn tails(&self, v: usize) -> (&'a [u32], &'a [u32]) {
        (
            self.pull.0.neighbors(v),
            self.pull.1.map_or(&[][..], |c| c.neighbors(v)),
        )
    }
}

/// Largest number of sources one [`MultiBfs::run`] advances together.
pub const BATCH: usize = 64;

/// A level is expanded bottom-up when its frontier has more than
/// 1/`PULL_SHARE` of all arcs.
const PULL_SHARE: usize = 8;

/// What [`MultiBfs::run`] reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    /// `node` was first reached at `level` (at least 1) by the searches
    /// from the sources whose bits are set in `sources`.
    Reached {
        /// The level, that is, the distance from those sources.
        level: u32,
        /// The node reached.
        node: u32,
        /// The sources that reached it at this level, one bit each.
        sources: u64,
    },
    /// Every node of `level` has been reported (not sent for a level that
    /// reaches nothing).
    LevelDone(u32),
}

/// Breadth-first search from up to [`BATCH`] sources at once, one bit per
/// source. Scratch space is reused across calls.
#[derive(Debug, Clone)]
pub struct MultiBfs {
    /// Sources that have reached each node.
    seen: Vec<u64>,
    /// Sources for which each node is on the current level.
    frontier: Vec<u64>,
    /// Sources whose current frontier has an arc into each node.
    next: Vec<u64>,
    /// Nodes on the current level for at least one source.
    active: Vec<u32>,
    /// Nodes reached by at least one source.
    reached: Vec<u32>,
    /// Nodes whose `next` set is not empty, in order of first push.
    touched: Vec<u32>,
}

impl MultiBfs {
    /// Scratch space for a network with `n` nodes.
    pub fn new(n: usize) -> MultiBfs {
        MultiBfs {
            seen: vec![0; n],
            frontier: vec![0; n],
            next: vec![0; n],
            active: Vec::new(),
            reached: Vec::new(),
            touched: Vec::new(),
        }
    }

    /// Searches along the arcs of `view` from `sources` (at most [`BATCH`]
    /// of them; source i is bit i), stopping after `max_level` levels, and
    /// reports every [`Step`] to `on`. A source is at level 0 of its own
    /// search and is not reported for it. Within a level, nodes are
    /// reported in increasing index order.
    ///
    /// Each level is expanded top-down (every frontier node ORs its set
    /// into the heads of its arcs) or, when the frontier's arcs are a large
    /// share of all arcs, bottom-up (every node not yet reached by all
    /// sources ORs the sets of the tails of its arcs, stopping once all
    /// sources have reached it), the direction-optimizing scheme of
    /// Beamer, Asanovic, and Patterson (2012, `beamer2012direction`)
    /// applied to source sets. Both give the same sets.
    pub fn run(&mut self, view: View, sources: &[u32], max_level: u32, mut on: impl FnMut(Step)) {
        assert!(sources.len() <= BATCH, "at most {BATCH} sources per batch");
        let n = view.node_count();
        let arcs = view.arc_count();
        let all = if sources.len() == BATCH {
            u64::MAX
        } else {
            (1u64 << sources.len()) - 1
        };
        self.active.clear();
        self.reached.clear();
        for (i, &s) in sources.iter().enumerate() {
            let s = s as usize;
            if self.seen[s] == 0 {
                self.active.push(s as u32);
                self.reached.push(s as u32);
            }
            self.seen[s] |= 1 << i;
            self.frontier[s] |= 1 << i;
        }
        let mut level = 0;
        while !self.active.is_empty() && level < max_level {
            level += 1;
            let frontier_arcs: usize = self.active.iter().map(|&v| view.degree(v as usize)).sum();
            let bottom_up = frontier_arcs > arcs / PULL_SHARE;
            self.touched.clear();
            if bottom_up {
                for w in 0..n {
                    let seen = self.seen[w];
                    if seen == all {
                        continue;
                    }
                    let (a, b) = view.tails(w);
                    let mut acc = 0;
                    'pull: for list in [a, b] {
                        for &v in list {
                            acc |= self.frontier[v as usize];
                            if acc | seen == all {
                                break 'pull;
                            }
                        }
                    }
                    self.next[w] = acc;
                }
            } else {
                // `next` is clean at this point.
                for &v in &self.active {
                    let f = self.frontier[v as usize];
                    let (a, b) = view.heads(v as usize);
                    for list in [a, b] {
                        for &w in list {
                            let x = &mut self.next[w as usize];
                            if *x == 0 {
                                self.touched.push(w);
                            }
                            *x |= f;
                        }
                    }
                }
            }
            // The frontier sets of this level have been read; clear them
            // so that only the next level's sets remain.
            for &v in &self.active {
                self.frontier[v as usize] = 0;
            }
            self.active.clear();
            // Visit the nodes with a non-empty `next` in increasing index
            // order: sort the list when it is short, scan every node when
            // it is long (or unknown, after a bottom-up step).
            if !bottom_up && self.touched.len() < n / 32 {
                self.touched.sort_unstable();
                for i in 0..self.touched.len() {
                    let w = self.touched[i] as usize;
                    self.settle(w, level, &mut on);
                }
            } else {
                for w in 0..n {
                    if self.next[w] != 0 {
                        self.settle(w, level, &mut on);
                    }
                }
            }
            if !self.active.is_empty() {
                on(Step::LevelDone(level));
            }
        }
        // `next` is cleared by every settle, and frontiers after every
        // level but the last, so only these need resetting.
        for &v in &self.reached {
            self.seen[v as usize] = 0;
            self.frontier[v as usize] = 0;
        }
    }

    /// Moves the sources pushed to `w` into its frontier, reporting those
    /// that reach it for the first time, and clears `next[w]`.
    #[inline]
    fn settle(&mut self, w: usize, level: u32, on: &mut impl FnMut(Step)) {
        let pushed = self.next[w];
        self.next[w] = 0;
        let new = pushed & !self.seen[w];
        if new != 0 {
            if self.seen[w] == 0 {
                self.reached.push(w as u32);
            }
            self.seen[w] |= new;
            self.frontier[w] = new;
            self.active.push(w as u32);
            on(Step::Reached {
                level,
                node: w as u32,
                sources: new,
            });
        }
    }

    /// Nodes reached by at least one source of the last run, sources
    /// included, in order of first discovery.
    pub fn reached(&self) -> &[u32] {
        &self.reached
    }
}

/// Calls `f(lane)` for every set bit of `bits`, lowest first.
#[inline]
pub fn for_each_bit(mut bits: u64, mut f: impl FnMut(usize)) {
    while bits != 0 {
        f(bits.trailing_zeros() as usize);
        bits &= bits - 1;
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
        assert_eq!(m.row(0).len(), 3);
        assert_eq!(m.row(0)[2], 1 << 1);
    }

    fn plain_bfs(adj: &Csr, s: usize) -> Vec<u32> {
        let mut d = vec![u32::MAX; adj.node_count()];
        d[s] = 0;
        let mut q = vec![s as u32];
        let mut head = 0;
        while head < q.len() {
            let v = q[head] as usize;
            head += 1;
            for &w in adj.neighbors(v) {
                if d[w as usize] == u32::MAX {
                    d[w as usize] = d[v] + 1;
                    q.push(w);
                }
            }
        }
        d
    }

    #[test]
    fn matches_plain_search_and_resets() {
        // Two components, a directed cycle, and an isolate.
        let g = Graph::from_edges(
            9,
            &[0, 1, 2, 3, 4, 5, 6, 7],
            &[1, 2, 3, 0, 5, 6, 7, 4],
            None,
            true,
        )
        .unwrap();
        let adj = g.out();
        let mut ms = MultiBfs::new(9);
        for sources in [vec![0u32, 4, 8, 2], (0..9).collect::<Vec<u32>>()] {
            let mut levels = vec![vec![u32::MAX; 9]; sources.len()];
            let mut done = Vec::new();
            ms.run(
                View::new(&g, DegreeMode::Out),
                &sources,
                u32::MAX,
                |step| match step {
                    Step::Reached {
                        level,
                        node,
                        sources,
                    } => for_each_bit(sources, |i| {
                        assert_eq!(levels[i][node as usize], u32::MAX);
                        levels[i][node as usize] = level;
                    }),
                    Step::LevelDone(level) => done.push(level),
                },
            );
            assert_eq!(done, vec![1, 2, 3]);
            for (i, &s) in sources.iter().enumerate() {
                let mut expected = plain_bfs(adj, s as usize);
                expected[s as usize] = u32::MAX;
                assert_eq!(levels[i], expected, "source {s}");
            }
            assert!(ms.seen.iter().all(|&x| x == 0));
            assert!(ms.next.iter().all(|&x| x == 0));
            assert!(ms.frontier.iter().all(|&x| x == 0));
        }
    }

    #[test]
    fn bounded_levels_and_repeated_sources() {
        // Path 0 - 1 - 2 - 3.
        let g = Graph::from_edges(4, &[0, 1, 2], &[1, 2, 3], None, false).unwrap();
        let mut ms = MultiBfs::new(4);
        let mut seen = Vec::new();
        ms.run(View::new(&g, DegreeMode::Out), &[0, 0, 3], 2, |step| {
            if let Step::Reached {
                level,
                node,
                sources,
            } = step
            {
                seen.push((level, node, sources));
            }
        });
        assert_eq!(
            seen,
            vec![(1, 1, 0b011), (1, 2, 0b100), (2, 1, 0b100), (2, 2, 0b011)]
        );
        assert_eq!(ms.reached(), &[0, 3, 1, 2]);
        assert!(ms.seen.iter().all(|&x| x == 0) && ms.next.iter().all(|&x| x == 0));
    }
}
