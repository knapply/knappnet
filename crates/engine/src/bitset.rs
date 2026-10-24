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

    /// Makes this an empty n-by-n matrix, reusing its storage.
    pub fn reset(&mut self, n: usize) {
        self.words_per_row = n.saturating_add(63) / 64;
        self.bits.clear();
        self.bits.resize(self.words_per_row * n, 0);
    }

    /// Number of 64-bit words in each row.
    #[inline]
    pub fn words_per_row(&self) -> usize {
        self.words_per_row
    }
}

/// Shortest-path levels and counts from one source in a small dense
/// network held as bit rows: the forward phase of Brandes' (2001,
/// `brandes2001faster`) accumulation, with every level found by word
/// operations on rows instead of a scan of adjacency lists.
///
/// A level is expanded top-down (the next level is the OR of the frontier's
/// out-rows, minus the nodes already reached) or, when the frontier is
/// large next to the nodes not yet reached, bottom-up (each node not yet
/// reached tests its in-row against the previous level), as in Beamer,
/// Asanovic, and Patterson (2012, `beamer2012direction`). The number of
/// shortest paths to a node at level l is the sum of the counts of its
/// in-neighbours at level l - 1, read from its in-row ANDed with that
/// level; while every count involved is 1 (always at level 2) the sum is a
/// population count. Counts are exact integers in `f64` up to 2^53, as in
/// the list-based search.
///
/// [`RowPaths::run_for`] counts paths only to the nodes a caller asks for
/// and to the nodes on their shortest paths, found by one pass back up the
/// levels; callers that need a few counts (ego betweenness needs those to
/// the nodes beyond ego) skip most of the population counts.
///
/// A search costs O(n w) word operations for n nodes and w words per row,
/// plus one addition per arc between consecutive levels beyond the second
/// on the shortest paths counted.
#[derive(Debug, Clone, Default)]
pub struct RowPaths {
    /// Level of every node; `u32::MAX` if not reached.
    level: Vec<u32>,
    /// Number of shortest paths from the source; 0 if not reached or not
    /// counted.
    sigma: Vec<f64>,
    /// Reached nodes by level, each level in increasing node order.
    order: Vec<u32>,
    /// Start of each level in `order`.
    starts: Vec<usize>,
    /// The nodes of each level as bits, `words` words per level.
    levels: Vec<u64>,
    /// The nodes of each level whose counts are needed.
    need: Vec<u64>,
    seen: Vec<u64>,
    words: usize,
}

impl RowPaths {
    /// Empty scratch space; it grows to fit the first search.
    pub fn new() -> RowPaths {
        RowPaths::default()
    }

    /// Searches from `source` in the network on nodes `0..n` whose arcs
    /// are the set bits of `out` (row a has bit b for an arc a -> b) and
    /// whose transpose is `inc`, stopping after `max_level` levels, and
    /// counts the shortest paths to every node reached. Both matrices are
    /// n by n; an undirected network passes the same matrix twice.
    pub fn run(
        &mut self,
        out: &BitMatrix,
        inc: &BitMatrix,
        n: usize,
        source: usize,
        max_level: u32,
    ) {
        self.search(out, inc, n, source, max_level);
        let total = self.levels.len();
        self.need.clear();
        self.need.extend_from_slice(&self.levels[..total]);
        self.count(inc);
    }

    /// As [`RowPaths::run`], but counts shortest paths only to the nodes
    /// at each level l whose bit is set in `wanted(l, w)` (word w of a node
    /// mask) and to every node on one of their shortest paths. Other
    /// counts read as 0.
    pub fn run_for(
        &mut self,
        out: &BitMatrix,
        inc: &BitMatrix,
        n: usize,
        source: usize,
        max_level: u32,
        wanted: impl Fn(usize, usize) -> u64,
    ) {
        self.search(out, inc, n, source, max_level);
        let words = self.words;
        let depth = self.depth();
        self.need.clear();
        self.need.resize(self.levels.len(), 0);
        for l in (1..depth).rev() {
            for w in 0..words {
                self.need[l * words + w] |= self.levels[l * words + w] & wanted(l, w);
            }
            if l < 2 {
                continue;
            }
            // Every in-neighbour one level up of a needed node is needed.
            let (above, here) = self.need.split_at_mut(l * words);
            let above = &mut above[(l - 1) * words..];
            let prev = &self.levels[(l - 1) * words..l * words];
            for (w, &bits) in here[..words].iter().enumerate() {
                for_each_bit(bits, |bit| {
                    for ((x, r), p) in above.iter_mut().zip(inc.row(w * 64 + bit)).zip(prev) {
                        *x |= r & p;
                    }
                });
            }
        }
        self.count(inc);
    }

    /// Finds the levels of a search without counting paths.
    fn search(
        &mut self,
        out: &BitMatrix,
        inc: &BitMatrix,
        n: usize,
        source: usize,
        max_level: u32,
    ) {
        let words = out.words_per_row();
        debug_assert_eq!(inc.words_per_row(), words);
        self.words = words;
        if self.level.len() < n {
            self.level.resize(n, u32::MAX);
            self.sigma.resize(n, 0.0);
        }
        for &v in &self.order {
            self.level[v as usize] = u32::MAX;
            self.sigma[v as usize] = 0.0;
        }
        self.order.clear();
        self.starts.clear();
        self.levels.clear();
        self.levels.resize(words, 0);
        self.seen.clear();
        self.seen.resize(words, 0);
        // Padding bits past node n - 1 count as reached, so they are never
        // candidates of a bottom-up level.
        if n % 64 != 0 {
            self.seen[words - 1] = !0u64 << (n % 64);
        }
        self.level[source] = 0;
        self.order.push(source as u32);
        self.starts.push(0);
        self.seen[source / 64] |= 1 << (source % 64);
        self.levels[source / 64] |= 1 << (source % 64);
        let mut unseen = n - 1;
        let mut depth = 0u32;
        while depth < max_level && unseen > 0 {
            let first = *self.starts.last().expect("a level");
            let width = self.order.len() - first;
            if width == 0 {
                break;
            }
            depth += 1;
            let start = self.order.len();
            self.starts.push(start);
            let base = self.levels.len();
            self.levels.resize(base + words, 0);
            let (done, next) = self.levels.split_at_mut(base);
            let frontier = &done[base - words..];
            if 2 * width < unseen {
                // Top-down.
                for &v in &self.order[first..] {
                    for (x, r) in next.iter_mut().zip(out.row(v as usize)) {
                        *x |= r;
                    }
                }
                for (x, s) in next.iter_mut().zip(&self.seen) {
                    *x &= !s;
                }
            } else {
                // Bottom-up.
                for (w, (x, &s)) in next.iter_mut().zip(&self.seen).enumerate() {
                    let mut hit = 0u64;
                    for_each_bit(!s, |bit| {
                        let row = inc.row(w * 64 + bit);
                        if row.iter().zip(frontier).any(|(a, b)| a & b != 0) {
                            hit |= 1 << bit;
                        }
                    });
                    *x = hit;
                }
            }
            for (w, (&x, s)) in next.iter().zip(self.seen.iter_mut()).enumerate() {
                let order = &mut self.order;
                for_each_bit(x, |bit| order.push((w * 64 + bit) as u32));
                *s |= x;
            }
            for &j in &self.order[start..] {
                self.level[j as usize] = depth;
            }
            unseen -= self.order.len() - start;
        }
        // Close the last level, unless it was empty (then drop its bits;
        // its start already closes the level before it).
        if self.starts.last() != Some(&self.order.len()) {
            self.starts.push(self.order.len());
        } else {
            self.levels.truncate(self.levels.len() - words);
        }
    }

    /// Counts the shortest paths to the nodes set in `need`, level by
    /// level; `need` must hold every needed node's in-neighbours one level
    /// up.
    fn count(&mut self, inc: &BitMatrix) {
        let words = self.words;
        self.sigma[self.order[0] as usize] = 1.0;
        let mut unit = true;
        for l in 1..self.depth() {
            unit = count_level(
                inc,
                &self.need[l * words..(l + 1) * words],
                &self.levels[(l - 1) * words..l * words],
                &mut self.sigma,
                unit,
            );
        }
    }

    /// Level of `v` in the last search; `u32::MAX` if it was not reached.
    #[inline]
    pub fn level(&self, v: usize) -> u32 {
        self.level[v]
    }

    /// Number of shortest paths from the source to `v` in the last search
    /// (1 for the source, 0 if `v` was not reached or not counted).
    #[inline]
    pub fn sigma(&self, v: usize) -> f64 {
        self.sigma[v]
    }

    /// Number of levels reached by the last search, the source's included.
    #[inline]
    pub fn depth(&self) -> usize {
        self.starts.len() - 1
    }

    /// The nodes at level `l` of the last search as bits (`l` less than
    /// [`RowPaths::depth`]).
    #[inline]
    pub fn level_bits(&self, l: usize) -> &[u64] {
        &self.levels[l * self.words..(l + 1) * self.words]
    }

    /// Nodes at level `l` of the last search, in increasing order (empty
    /// past the last level).
    #[inline]
    pub fn nodes_at(&self, l: usize) -> &[u32] {
        if l + 1 >= self.starts.len() {
            return &[];
        }
        &self.order[self.starts[l]..self.starts[l + 1]]
    }
}

/// Sets `sigma[j]` for every node j set in `nodes` to the sum of `sigma`
/// over the nodes set in both `inc.row(j)` and `prev` (a population count
/// when every count in `prev` that is read is 1, `unit`), and returns
/// whether every new count is 1. Population counts dominate the cost of
/// [`RowPaths`], so on x86-64 processors with the POPCNT instruction this
/// runs through a copy compiled to use it; both copies compute the same
/// integers.
fn count_level(
    inc: &BitMatrix,
    nodes: &[u64],
    prev: &[u64],
    sigma: &mut [f64],
    unit: bool,
) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        #[target_feature(enable = "popcnt")]
        unsafe fn fast(
            inc: &BitMatrix,
            nodes: &[u64],
            prev: &[u64],
            sigma: &mut [f64],
            unit: bool,
        ) -> bool {
            count_level_body(inc, nodes, prev, sigma, unit)
        }
        if std::is_x86_feature_detected!("popcnt") {
            // SAFETY: the processor supports POPCNT.
            return unsafe { fast(inc, nodes, prev, sigma, unit) };
        }
    }
    count_level_body(inc, nodes, prev, sigma, unit)
}

#[inline(always)]
fn count_level_body(
    inc: &BitMatrix,
    nodes: &[u64],
    prev: &[u64],
    sigma: &mut [f64],
    unit: bool,
) -> bool {
    let mut all_unit = true;
    for (w, &bits) in nodes.iter().enumerate() {
        for_each_bit(bits, |bit| {
            let j = w * 64 + bit;
            let row = inc.row(j);
            let c = if unit {
                row.iter()
                    .zip(prev)
                    .map(|(a, b)| (a & b).count_ones())
                    .sum::<u32>() as f64
            } else {
                let mut total = 0.0;
                for (v, (a, b)) in row.iter().zip(prev).enumerate() {
                    for_each_bit(a & b, |bit| total += sigma[v * 64 + bit]);
                }
                total
            };
            all_unit &= c == 1.0;
            sigma[j] = c;
        });
    }
    all_unit
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

/// Whether searching from [`BATCH`] sources at a time is expected to beat
/// one plain search per source when searches stop after `max_level` levels.
/// A batch saves work only where the searches of its sources overlap; with
/// average out-degree d, a search of radius r reaches about d^r nodes, so
/// the searches of a batch are likely to overlap once [`BATCH`] d^r reaches
/// n. Below that, plain searches touch less memory per node reached.
pub fn batches_pay(view: View, max_level: u32) -> bool {
    let n = view.node_count();
    if n == 0 {
        return false;
    }
    let d = view.arc_count() as f64 / n as f64;
    let ball = d.powf(max_level.min(64) as f64).min(n as f64);
    BATCH as f64 * ball >= n as f64
}

/// Scratch space for plain breadth-first searches from one source along a
/// [`View`], reused across sources.
#[derive(Debug, Clone)]
pub struct SingleBfs {
    /// Level of every node reached by the last search; `u32::MAX` if not.
    level: Vec<u32>,
    /// Nodes reached by the last search, in order of discovery.
    order: Vec<u32>,
}

impl SingleBfs {
    /// Scratch space for a network with `n` nodes.
    pub fn new(n: usize) -> SingleBfs {
        SingleBfs {
            level: vec![u32::MAX; n],
            order: Vec::new(),
        }
    }

    /// Searches from `source` along the arcs of `view`, stopping after
    /// `max_level` levels, and returns the nodes reached (the source
    /// first) in order of discovery, which is non-decreasing level.
    pub fn run(&mut self, view: View, source: u32, max_level: u32) -> &[u32] {
        for &v in &self.order {
            self.level[v as usize] = u32::MAX;
        }
        self.order.clear();
        self.level[source as usize] = 0;
        self.order.push(source);
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let next = self.level[v] + 1;
            if next > max_level {
                break;
            }
            let (a, b) = view.heads(v);
            for list in [a, b] {
                for &w in list {
                    if self.level[w as usize] == u32::MAX {
                        self.level[w as usize] = next;
                        self.order.push(w);
                    }
                }
            }
        }
        &self.order
    }

    /// Level of `v` in the last search; `u32::MAX` if it was not reached.
    #[inline]
    pub fn level(&self, v: u32) -> u32 {
        self.level[v as usize]
    }
}

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
    /// reported in an order fixed by the graph and the sources alone.
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
            // Visit the nodes with a non-empty `next`: in the order they
            // were first pushed to after a top-down step, by a scan of
            // every node after a bottom-up step.
            if !bottom_up {
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
