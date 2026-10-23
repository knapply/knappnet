//! Random graphs: G(n, p), G(n, m), the configuration model, U|MAN,
//! Watts-Strogatz small worlds, Barabasi-Albert preferential attachment,
//! and stochastic block models.
//!
//! Sources:
//! - G(n, p), each possible edge present independently with probability p
//!   (Gilbert 1959, `gilbert1959random`), and G(n, m), a graph drawn
//!   uniformly from those with exactly m edges (Erdos and Renyi 1959,
//!   `erdos1959random`).
//! - G(n, p) and the blocks of a stochastic block model are sampled by
//!   geometric skipping: the gap to the next present edge in a fixed order
//!   of candidate pairs is drawn directly, so the cost is O(n + m) rather
//!   than O(n^2) (Batagelj and Brandes 2005, `batagelj2005efficient`).
//! - G(n, m) draws m distinct pairs with Floyd's algorithm for sampling
//!   without replacement (Bentley and Floyd 1987, `bentley1987sample`).
//! - The configuration model matches the "stubs" (edge ends) of a degree
//!   sequence uniformly at random (Bollobas 1980, `bollobas1980probabilistic`;
//!   Molloy and Reed 1995, `molloy1995critical`; Newman, Strogatz, and Watts
//!   2001, `newman2001random`).
//! - U|MAN draws a directed graph uniformly from those with given numbers
//!   of mutual, asymmetric, and null dyads (Holland and Leinhardt 1976,
//!   `holland1976local`).
//! - Small worlds follow the rewiring procedure of Watts and Strogatz
//!   (1998, `watts1998collective`).
//! - Preferential attachment follows Barabasi and Albert (1999,
//!   `barabasi1999emergence`), with the degree-weighted node list of
//!   Batagelj and Brandes (2005) for linear time.
//! - Stochastic block models (Holland, Laskey, and Leinhardt 1983,
//!   `holland1983stochastic`) join a node of block r and a node of block s
//!   independently with probability P[r][s].
//!
//! Conventions:
//! - Nodes are `0..n`. Every generator returns an edge list in increasing
//!   order of (from, to), with from < to in undirected graphs (from <= to
//!   when loops are allowed), except the configuration model multigraph,
//!   whose edges are sorted the same way and may repeat.
//! - Randomness comes from [`crate::rng::Stream`]: replicate r of a batch
//!   uses stream (seed, r), so a seed reproduces every replicate on any
//!   platform and with any number of threads.
//! - Unless stated otherwise graphs are simple: no self-loops, no parallel
//!   edges. G(n, p) and G(n, m) can include loops on request (n loops become
//!   possible, as in igraph's `sample_gnp(loops = TRUE)`).
//! - The configuration model can return the multigraph (loops and parallel
//!   edges kept, so every node has exactly its degree, a loop counting 2 in
//!   undirected graphs) or a simple graph. The simple variant repeats the
//!   whole matching until it contains no loop and no parallel edge, which
//!   gives every simple graph with the degree sequence the same
//!   probability (each corresponds to the same number of stub matchings,
//!   the product of the factorials of the degrees, and of in- and
//!   out-degrees for directed graphs). Sequences for which simple matchings
//!   are rare stop with an error after `max_tries` attempts rather than
//!   running indefinitely; erasing loops and parallel edges instead would
//!   change the degrees and is not offered.
//! - Watts-Strogatz starts from a ring in which every node is joined to its
//!   k nearest neighbors (k / 2 on each side, k even). Each lattice edge
//!   (u, u + j) is considered once, lap by lap for j = 1, ..., k / 2 and
//!   node by node within a lap, and with probability p its far end is
//!   moved to a node drawn uniformly from those that are neither u nor
//!   already adjacent to u (Watts and Strogatz 1998, as NetworkX's
//!   `watts_strogatz_graph`). A node already adjacent to every other is
//!   left alone. The graph stays simple and keeps n k / 2 edges.
//! - Preferential attachment starts from a star on m + 1 nodes (node 0 the
//!   center), as NetworkX's `barabasi_albert_graph`; each later node joins
//!   m distinct existing nodes, each drawn with probability proportional to
//!   its degree (repeated draws of a node already chosen are discarded).
//!   The graph has m (n - m) edges and is undirected and simple.
//! - Stochastic block models take block sizes (nodes are numbered block by
//!   block) and a k-by-k matrix of probabilities, which must be symmetric
//!   for undirected graphs.

use std::collections::HashSet;

use crate::error::{Error, Result};
use crate::parallel;
use crate::rng::Stream;

/// An edge list on nodes `0..n`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeList {
    /// Number of nodes.
    pub n: usize,
    /// Whether edges are directed.
    pub directed: bool,
    /// Edge sources.
    pub from: Vec<u32>,
    /// Edge targets.
    pub to: Vec<u32>,
}

impl EdgeList {
    fn new(n: usize, directed: bool) -> EdgeList {
        EdgeList {
            n,
            directed,
            from: Vec::new(),
            to: Vec::new(),
        }
    }

    fn push(&mut self, a: u64, b: u64) {
        self.from.push(a as u32);
        self.to.push(b as u32);
    }

    /// Number of edges.
    pub fn len(&self) -> usize {
        self.from.len()
    }

    /// Whether there are no edges.
    pub fn is_empty(&self) -> bool {
        self.from.is_empty()
    }

    fn sort(&mut self) {
        let mut pairs: Vec<(u32, u32)> = self
            .from
            .iter()
            .copied()
            .zip(self.to.iter().copied())
            .collect();
        pairs.sort_unstable();
        self.from = pairs.iter().map(|p| p.0).collect();
        self.to = pairs.iter().map(|p| p.1).collect();
    }
}

fn check_n(n: usize) -> Result<()> {
    if n > u32::MAX as usize {
        return Err(Error::TooLarge {
            what: "node count",
            value: n,
        });
    }
    Ok(())
}

fn check_probability(p: f64, what: &str) -> Result<()> {
    if !(0.0..=1.0).contains(&p) {
        return Err(Error::InvalidArgument(format!(
            "{what} must be a probability between 0 and 1, not {p}"
        )));
    }
    Ok(())
}

/// Number of candidate pairs among n nodes.
fn pair_count(n: u64, directed: bool, loops: bool) -> u64 {
    match (directed, loops) {
        (true, false) => n * n.saturating_sub(1),
        (true, true) => n * n,
        (false, false) => n * n.saturating_sub(1) / 2,
        (false, true) => n * (n + 1) / 2,
    }
}

/// Integer square root (largest r with r * r <= x).
fn isqrt(x: u64) -> u64 {
    let mut r = (x as f64).sqrt() as u64;
    while r > 0 && r.checked_mul(r).map_or(true, |s| s > x) {
        r -= 1;
    }
    while (r + 1).checked_mul(r + 1).map_or(false, |s| s <= x) {
        r += 1;
    }
    r
}

/// The k-th candidate pair, in increasing order of (row, column), where
/// the candidates are: directed, all ordered pairs (without the diagonal
/// unless `loops`); undirected, pairs (a, b) with a < b (a <= b with
/// `loops`), ordered by b then a. The undirected pair is returned with the
/// smaller node first.
fn decode(k: u64, n: u64, directed: bool, loops: bool) -> (u64, u64) {
    match (directed, loops) {
        (true, true) => (k / n, k % n),
        (true, false) => {
            let row = k / (n - 1);
            let col = k % (n - 1);
            (row, if col >= row { col + 1 } else { col })
        }
        (false, false) => {
            // Column b holds pairs (0, b), ..., (b - 1, b): b (b - 1) / 2
            // pairs come before it.
            let mut b = (1 + isqrt(1 + 8 * k)) / 2;
            while b * (b - 1) / 2 > k {
                b -= 1;
            }
            while (b + 1) * b / 2 <= k {
                b += 1;
            }
            (k - b * (b - 1) / 2, b)
        }
        (false, true) => {
            // Column b holds (0, b), ..., (b, b): b (b + 1) / 2 before it.
            let mut b = (isqrt(1 + 8 * k) - 1) / 2;
            while b * (b + 1) / 2 > k {
                b -= 1;
            }
            while (b + 1) * (b + 2) / 2 <= k {
                b += 1;
            }
            (k - b * (b + 1) / 2, b)
        }
    }
}

/// Calls `emit(k)` for every index k in `0..count` kept independently with
/// probability p, in increasing order, by geometric skipping (Batagelj and
/// Brandes 2005).
fn bernoulli_indices(count: u64, p: f64, rng: &mut Stream, mut emit: impl FnMut(u64)) {
    if count == 0 || p <= 0.0 {
        return;
    }
    if p >= 1.0 {
        (0..count).for_each(emit);
        return;
    }
    let log_q = (-p).ln_1p();
    let mut k: u64 = 0;
    let mut first = true;
    loop {
        // 1 - u lies in (0, 1], so the logarithm is finite and <= 0.
        let u = rng.uniform();
        let skip = ((-u).ln_1p() / log_q).floor();
        let step = if first { skip } else { skip + 1.0 };
        if !(step < (count - k) as f64) {
            return;
        }
        k += step as u64;
        if k >= count {
            return;
        }
        first = false;
        emit(k);
    }
}

/// One G(n, p) graph from `rng`.
fn gnp_one(n: usize, p: f64, directed: bool, loops: bool, rng: &mut Stream) -> EdgeList {
    let mut out = EdgeList::new(n, directed);
    let nn = n as u64;
    bernoulli_indices(pair_count(nn, directed, loops), p, rng, |k| {
        let (a, b) = decode(k, nn, directed, loops);
        out.push(a, b);
    });
    // Undirected candidates are ordered by column; sort by row.
    if !directed {
        out.sort();
    }
    out
}

/// One G(n, m) graph from `rng`.
fn gnm_one(n: usize, m: u64, directed: bool, loops: bool, rng: &mut Stream) -> EdgeList {
    let nn = n as u64;
    let total = pair_count(nn, directed, loops);
    let mut chosen: Vec<u64> = if m.saturating_mul(2) <= total {
        // Floyd's algorithm: m distinct indices, each m-subset equally
        // likely.
        let mut set: HashSet<u64> = HashSet::with_capacity(m as usize);
        for j in (total - m)..total {
            let t = rng.below(j + 1);
            if !set.insert(t) {
                set.insert(j);
            }
        }
        set.into_iter().collect()
    } else {
        // Dense: choose the total - m pairs to leave out instead.
        let mut skip: HashSet<u64> = HashSet::with_capacity((total - m) as usize);
        for j in m..total {
            let t = rng.below(j + 1);
            if !skip.insert(t) {
                skip.insert(j);
            }
        }
        (0..total).filter(|k| !skip.contains(k)).collect()
    };
    chosen.sort_unstable();
    let mut out = EdgeList::new(n, directed);
    for k in chosen {
        let (a, b) = decode(k, nn, directed, loops);
        out.push(a, b);
    }
    out.sort();
    out
}

/// Runs `one(rng)` for replicates `0..reps`, replicate r with stream
/// (seed, r), in parallel; results are in replicate order.
fn replicate<F>(reps: usize, seed: u64, threads: usize, one: F) -> Result<Vec<EdgeList>>
where
    F: Fn(&mut Stream) -> Result<EdgeList> + Sync,
{
    let parts = parallel::for_each_item(reps, threads, Vec::new, |acc, r| {
        let mut rng = Stream::new(seed, r as u64);
        acc.push((r, one(&mut rng)));
    });
    let mut all: Vec<(usize, Result<EdgeList>)> = parts.into_iter().flatten().collect();
    all.sort_unstable_by_key(|x| x.0);
    all.into_iter().map(|x| x.1).collect()
}

/// `reps` independent G(n, p) graphs (Gilbert 1959). See the module
/// documentation.
pub fn gnp(
    n: usize,
    p: f64,
    directed: bool,
    loops: bool,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    check_n(n)?;
    check_probability(p, "p")?;
    replicate(reps, seed, threads, |rng| {
        Ok(gnp_one(n, p, directed, loops, rng))
    })
}

/// `reps` independent G(n, m) graphs (Erdos and Renyi 1959). See the
/// module documentation.
pub fn gnm(
    n: usize,
    m: u64,
    directed: bool,
    loops: bool,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    check_n(n)?;
    let total = pair_count(n as u64, directed, loops);
    if m > total {
        return Err(Error::InvalidArgument(format!(
            "a {} graph on {n} nodes has at most {total} possible edges, so {m} edges \
             cannot be placed",
            if directed { "directed" } else { "undirected" }
        )));
    }
    replicate(reps, seed, threads, |rng| {
        Ok(gnm_one(n, m, directed, loops, rng))
    })
}

/// Degree sequences for [`configuration`]: one sequence for undirected
/// graphs, out- and in-degrees for directed graphs.
#[derive(Debug, Clone)]
pub enum Degrees {
    /// Degree of every node of an undirected graph (a loop counts 2).
    Undirected(Vec<u64>),
    /// Out- and in-degree of every node of a directed graph.
    Directed(Vec<u64>, Vec<u64>),
}

fn stubs(degrees: &[u64]) -> Result<Vec<u32>> {
    let total: u64 = degrees.iter().sum();
    if total > (u32::MAX as u64) {
        return Err(Error::TooLarge {
            what: "sum of degrees",
            value: total as usize,
        });
    }
    let mut out = Vec::with_capacity(total as usize);
    for (v, &d) in degrees.iter().enumerate() {
        out.extend(std::iter::repeat(v as u32).take(d as usize));
    }
    Ok(out)
}

fn shuffle(x: &mut [u32], rng: &mut Stream) {
    for i in (1..x.len()).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        x.swap(i, j);
    }
}

fn is_simple(e: &EdgeList) -> bool {
    // Edges are sorted, so repeats are adjacent.
    e.from.iter().zip(&e.to).all(|(a, b)| a != b)
        && e.from
            .windows(2)
            .zip(e.to.windows(2))
            .all(|(f, t)| f[0] != f[1] || t[0] != t[1])
}

fn configuration_one(degrees: &Degrees, rng: &mut Stream) -> Result<EdgeList> {
    match degrees {
        Degrees::Undirected(d) => {
            let mut s = stubs(d)?;
            shuffle(&mut s, rng);
            let mut out = EdgeList::new(d.len(), false);
            for pair in s.chunks_exact(2) {
                let (a, b) = (pair[0].min(pair[1]), pair[0].max(pair[1]));
                out.push(a as u64, b as u64);
            }
            out.sort();
            Ok(out)
        }
        Degrees::Directed(dout, din) => {
            let out_stubs = stubs(dout)?;
            let mut in_stubs = stubs(din)?;
            shuffle(&mut in_stubs, rng);
            let mut out = EdgeList::new(dout.len(), true);
            for (&a, &b) in out_stubs.iter().zip(&in_stubs) {
                out.push(a as u64, b as u64);
            }
            out.sort();
            Ok(out)
        }
    }
}

/// `reps` configuration-model graphs with the given degrees. With
/// `simple`, each replicate repeats the matching (up to `max_tries` times)
/// until it has no loop or parallel edge. See the module documentation.
pub fn configuration(
    degrees: &Degrees,
    simple: bool,
    max_tries: usize,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    match degrees {
        Degrees::Undirected(d) => {
            check_n(d.len())?;
            let total: u64 = d.iter().sum();
            if total % 2 == 1 {
                return Err(Error::InvalidArgument(format!(
                    "the degrees sum to {total}, an odd number, so no undirected graph has them"
                )));
            }
        }
        Degrees::Directed(o, i) => {
            check_n(o.len())?;
            if o.len() != i.len() {
                return Err(Error::LengthMismatch {
                    what: "in-degrees",
                    expected: o.len(),
                    found: i.len(),
                });
            }
            let (so, si): (u64, u64) = (o.iter().sum(), i.iter().sum());
            if so != si {
                return Err(Error::InvalidArgument(format!(
                    "out-degrees sum to {so} but in-degrees sum to {si}; they must be equal"
                )));
            }
        }
    }
    if simple && max_tries == 0 {
        return Err(Error::InvalidArgument(
            "max_tries must be at least 1".into(),
        ));
    }
    replicate(reps, seed, threads, |rng| {
        if !simple {
            return configuration_one(degrees, rng);
        }
        for _ in 0..max_tries {
            let e = configuration_one(degrees, rng)?;
            if is_simple(&e) {
                return Ok(e);
            }
        }
        Err(Error::InvalidArgument(format!(
            "no simple graph was found in {max_tries} attempts; simple graphs are rare \
             (or impossible) for this degree sequence in the configuration model"
        )))
    })
}

fn uman_one(n: usize, mutual: u64, asym: u64, rng: &mut Stream) -> EdgeList {
    // Choose the mutual + asym non-null dyads uniformly among the
    // n (n - 1) / 2, then which of them are mutual, then a direction for
    // each asymmetric one.
    let chosen = gnm_one(n, mutual + asym, false, false, rng);
    let k = chosen.len();
    let mut order: Vec<usize> = (0..k).collect();
    // Partial Fisher-Yates: the first `mutual` positions are a uniform
    // random subset.
    for i in 0..(mutual as usize).min(k) {
        let j = i + rng.below((k - i) as u64) as usize;
        order.swap(i, j);
    }
    let mut is_mutual = vec![false; k];
    for &i in &order[..mutual as usize] {
        is_mutual[i] = true;
    }
    let mut out = EdgeList::new(n, true);
    for i in 0..k {
        let (a, b) = (chosen.from[i] as u64, chosen.to[i] as u64);
        if is_mutual[i] {
            out.push(a, b);
            out.push(b, a);
        } else if rng.next_u32() & 1 == 0 {
            out.push(a, b);
        } else {
            out.push(b, a);
        }
    }
    out.sort();
    out
}

/// `reps` directed graphs drawn uniformly from those with `mutual` mutual
/// and `asym` asymmetric dyads, the rest null (U|MAN; Holland and
/// Leinhardt 1976).
pub fn uman(
    n: usize,
    mutual: u64,
    asym: u64,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    check_n(n)?;
    let dyads = pair_count(n as u64, false, false);
    if mutual.saturating_add(asym) > dyads {
        return Err(Error::InvalidArgument(format!(
            "{n} nodes have {dyads} dyads, fewer than the {mutual} mutual plus {asym} \
             asymmetric dyads requested"
        )));
    }
    replicate(reps, seed, threads, |rng| {
        Ok(uman_one(n, mutual, asym, rng))
    })
}

fn small_world_one(n: usize, k: usize, p: f64, rng: &mut Stream) -> EdgeList {
    let mut adj: Vec<HashSet<u32>> = vec![HashSet::new(); n];
    let half = k / 2;
    for u in 0..n {
        for j in 1..=half {
            let v = (u + j) % n;
            adj[u].insert(v as u32);
            adj[v].insert(u as u32);
        }
    }
    for j in 1..=half {
        for u in 0..n {
            let v = ((u + j) % n) as u32;
            if rng.uniform() >= p {
                continue;
            }
            if adj[u].len() >= n - 1 {
                continue;
            }
            // The lattice edge may already have been moved away from u by
            // an earlier rewiring of v; only edges still present move.
            if !adj[u].contains(&v) {
                continue;
            }
            let w = loop {
                let w = rng.below(n as u64) as u32;
                if w as usize != u && !adj[u].contains(&w) {
                    break w;
                }
            };
            adj[u].remove(&v);
            adj[v as usize].remove(&(u as u32));
            adj[u].insert(w);
            adj[w as usize].insert(u as u32);
        }
    }
    let mut out = EdgeList::new(n, false);
    for (u, nbrs) in adj.iter().enumerate() {
        let mut row: Vec<u32> = nbrs.iter().copied().filter(|&v| v as usize > u).collect();
        row.sort_unstable();
        for v in row {
            out.push(u as u64, v as u64);
        }
    }
    out
}

/// `reps` Watts-Strogatz small-world graphs: a ring of n nodes each joined
/// to its k nearest neighbors, with each edge rewired with probability p.
/// See the module documentation.
pub fn small_world(
    n: usize,
    k: usize,
    p: f64,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    check_n(n)?;
    check_probability(p, "p")?;
    if k % 2 == 1 || k >= n {
        return Err(Error::InvalidArgument(format!(
            "k (the number of nearest neighbors) must be even and less than the number \
             of nodes, not {k} with {n} nodes"
        )));
    }
    replicate(reps, seed, threads, |rng| Ok(small_world_one(n, k, p, rng)))
}

fn preferential_one(n: usize, m: usize, rng: &mut Stream) -> EdgeList {
    let mut out = EdgeList::new(n, false);
    // Every node appears in `pool` once per edge end, so a uniform draw
    // from the pool picks a node with probability proportional to degree.
    let mut pool: Vec<u32> = Vec::with_capacity(2 * m * n);
    for leaf in 1..=m {
        out.push(0, leaf as u64);
        pool.push(0);
        pool.push(leaf as u32);
    }
    let mut targets: Vec<u32> = Vec::with_capacity(m);
    for source in (m + 1)..n {
        targets.clear();
        while targets.len() < m {
            let t = pool[rng.below(pool.len() as u64) as usize];
            if !targets.contains(&t) {
                targets.push(t);
            }
        }
        for &t in &targets {
            out.push(t as u64, source as u64);
            pool.push(t);
            pool.push(source as u32);
        }
    }
    out.sort();
    out
}

/// `reps` Barabasi-Albert preferential attachment graphs with n nodes,
/// each new node bringing m edges. See the module documentation.
pub fn preferential(
    n: usize,
    m: usize,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    check_n(n)?;
    if m == 0 || m >= n {
        return Err(Error::InvalidArgument(format!(
            "m (edges per new node) must be at least 1 and less than the number of nodes, \
             not {m} with {n} nodes"
        )));
    }
    replicate(reps, seed, threads, |rng| Ok(preferential_one(n, m, rng)))
}

fn sbm_one(
    sizes: &[usize],
    probs: &[f64],
    directed: bool,
    loops: bool,
    rng: &mut Stream,
) -> EdgeList {
    let k = sizes.len();
    let mut start = vec![0u64; k + 1];
    for r in 0..k {
        start[r + 1] = start[r] + sizes[r] as u64;
    }
    let n = start[k] as usize;
    let mut out = EdgeList::new(n, directed);
    for r in 0..k {
        for s in 0..k {
            if !directed && s < r {
                continue;
            }
            let p = probs[r * k + s];
            let (nr, ns) = (sizes[r] as u64, sizes[s] as u64);
            if r == s {
                bernoulli_indices(pair_count(nr, directed, loops), p, rng, |idx| {
                    let (a, b) = decode(idx, nr, directed, loops);
                    out.push(start[r] + a, start[r] + b);
                });
            } else {
                bernoulli_indices(nr * ns, p, rng, |idx| {
                    out.push(start[r] + idx / ns, start[s] + idx % ns);
                });
            }
        }
    }
    out.sort();
    out
}

/// `reps` stochastic block model graphs. `sizes` gives the number of
/// nodes in each block and `probs` the k-by-k matrix of tie probabilities
/// in row-major order. See the module documentation.
pub fn sbm(
    sizes: &[usize],
    probs: &[f64],
    directed: bool,
    loops: bool,
    reps: usize,
    seed: u64,
    threads: usize,
) -> Result<Vec<EdgeList>> {
    let k = sizes.len();
    if probs.len() != k * k {
        return Err(Error::LengthMismatch {
            what: "block probability matrix",
            expected: k * k,
            found: probs.len(),
        });
    }
    let n: usize = sizes.iter().sum();
    check_n(n)?;
    for &p in probs {
        check_probability(p, "every block probability")?;
    }
    if !directed {
        for r in 0..k {
            for s in 0..r {
                if probs[r * k + s] != probs[s * k + r] {
                    return Err(Error::InvalidArgument(
                        "the block probability matrix of an undirected graph must be symmetric"
                            .into(),
                    ));
                }
            }
        }
    }
    replicate(reps, seed, threads, |rng| {
        Ok(sbm_one(sizes, probs, directed, loops, rng))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one(v: Result<Vec<EdgeList>>) -> EdgeList {
        v.unwrap().pop().unwrap()
    }

    #[test]
    fn decode_enumerates_every_pair_once_in_order() {
        for &(directed, loops) in &[(true, false), (true, true), (false, false), (false, true)] {
            let n = 7u64;
            let total = pair_count(n, directed, loops);
            let pairs: Vec<(u64, u64)> =
                (0..total).map(|k| decode(k, n, directed, loops)).collect();
            let mut expect = Vec::new();
            for a in 0..n {
                for b in 0..n {
                    let ok = match (directed, loops) {
                        (true, false) => a != b,
                        (true, true) => true,
                        (false, false) => a < b,
                        (false, true) => a <= b,
                    };
                    if ok {
                        expect.push((a, b));
                    }
                }
            }
            let mut sorted = pairs.clone();
            sorted.sort_unstable();
            assert_eq!(sorted, expect);
            if directed {
                assert_eq!(pairs, expect);
            }
        }
        // Large indices decode exactly.
        let n = 3_000_000u64;
        let k = pair_count(n, false, false) - 1;
        assert_eq!(decode(k, n, false, false), (n - 2, n - 1));
        let k = pair_count(n, false, true) - 1;
        assert_eq!(decode(k, n, false, true), (n - 1, n - 1));
    }

    #[test]
    fn gnp_extremes_and_reproducibility() {
        assert!(one(gnp(10, 0.0, false, false, 1, 1, 1)).is_empty());
        let full = one(gnp(6, 1.0, false, false, 1, 1, 1));
        assert_eq!(full.len(), 15);
        assert_eq!(one(gnp(5, 1.0, true, true, 1, 1, 1)).len(), 25);
        let a = gnp(200, 0.05, true, false, 4, 9, 1).unwrap();
        let b = gnp(200, 0.05, true, false, 4, 9, 3).unwrap();
        assert_eq!(a, b);
        assert_ne!(a[0], a[1]);
        for e in &a {
            assert!(is_simple(e));
        }
        assert!(gnp(3, 1.5, false, false, 1, 1, 1).is_err());
    }

    #[test]
    fn gnp_edge_count_is_binomial() {
        // 2000 replicates of G(30, 0.2): the mean edge count should be
        // 87 = 435 * 0.2 within a few standard errors.
        let reps = gnp(30, 0.2, false, false, 2000, 3, 0).unwrap();
        let mean = reps.iter().map(|e| e.len() as f64).sum::<f64>() / 2000.0;
        let se = (435.0f64 * 0.2 * 0.8 / 2000.0).sqrt();
        assert!((mean - 87.0).abs() < 4.0 * se, "{mean}");
    }

    #[test]
    fn gnm_has_exactly_m_distinct_edges() {
        for &(directed, loops) in &[(true, false), (true, true), (false, false), (false, true)] {
            let total = pair_count(9, directed, loops);
            for &m in &[0u64, 1, 10, total - 5, total] {
                let e = one(gnm(9, m, directed, loops, 1, 5, 1));
                assert_eq!(e.len() as u64, m);
                let mut seen = HashSet::new();
                for (&a, &b) in e.from.iter().zip(&e.to) {
                    assert!(loops || a != b);
                    assert!(directed || a <= b);
                    assert!(seen.insert((a, b)));
                }
            }
        }
        assert!(gnm(4, 7, false, false, 1, 1, 1).is_err());
    }

    #[test]
    fn gnm_is_uniform_over_pairs() {
        // Every pair of 5 nodes appears in about m / 10 of the graphs.
        let reps = gnm(5, 3, false, false, 20_000, 2, 0).unwrap();
        let mut count = [[0u32; 5]; 5];
        for e in &reps {
            for (&a, &b) in e.from.iter().zip(&e.to) {
                count[a as usize][b as usize] += 1;
            }
        }
        for a in 0..5 {
            for b in (a + 1)..5 {
                assert!((5600..6400).contains(&count[a][b]), "{}", count[a][b]);
            }
        }
    }

    #[test]
    fn configuration_keeps_degrees() {
        let d = vec![3u64, 2, 2, 2, 1, 0, 4];
        for e in configuration(&Degrees::Undirected(d.clone()), false, 1, 50, 1, 0).unwrap() {
            let mut deg = vec![0u64; d.len()];
            for (&a, &b) in e.from.iter().zip(&e.to) {
                deg[a as usize] += 1;
                deg[b as usize] += 1;
            }
            assert_eq!(deg, d);
        }
        for e in configuration(&Degrees::Undirected(d.clone()), true, 10_000, 50, 1, 0).unwrap() {
            assert!(is_simple(&e));
        }
        let (o, i) = (vec![2u64, 1, 0, 1], vec![0u64, 1, 2, 1]);
        for e in configuration(
            &Degrees::Directed(o.clone(), i.clone()),
            true,
            1000,
            20,
            1,
            0,
        )
        .unwrap()
        {
            assert!(is_simple(&e));
            let mut dout = vec![0u64; 4];
            let mut din = vec![0u64; 4];
            for (&a, &b) in e.from.iter().zip(&e.to) {
                dout[a as usize] += 1;
                din[b as usize] += 1;
            }
            assert_eq!((dout, din), (o.clone(), i.clone()));
        }
        assert!(configuration(&Degrees::Undirected(vec![1, 1, 1]), false, 1, 1, 1, 1).is_err());
        // A node of degree 4 among 3 nodes cannot be simple.
        assert!(configuration(&Degrees::Undirected(vec![4, 2, 2]), true, 50, 1, 1, 1).is_err());
    }

    #[test]
    fn simple_configuration_is_uniform() {
        // Degrees (2, 2, 2, 1, 1) have 3 simple realizations up to the
        // path/triangle shapes; count labelled graphs exactly: a triangle
        // on {0, 1, 2} plus the edge 3-4, or a path 3-a-b-c-4 through all
        // of 0, 1, 2 (3! = 6 orders, each its own graph). 7 graphs, each
        // with probability 1/7.
        let reps = configuration(
            &Degrees::Undirected(vec![2, 2, 2, 1, 1]),
            true,
            1000,
            14_000,
            4,
            0,
        )
        .unwrap();
        let mut counts = std::collections::HashMap::new();
        for e in reps {
            *counts.entry((e.from, e.to)).or_insert(0u32) += 1;
        }
        assert_eq!(counts.len(), 7);
        for &c in counts.values() {
            assert!((1800..2200).contains(&c), "{c}");
        }
    }

    #[test]
    fn uman_has_exact_dyad_census() {
        for e in uman(12, 7, 15, 30, 2, 0).unwrap() {
            let set: HashSet<(u32, u32)> =
                e.from.iter().copied().zip(e.to.iter().copied()).collect();
            assert_eq!(set.len(), e.len());
            let mutual = set
                .iter()
                .filter(|&&(a, b)| a < b && set.contains(&(b, a)))
                .count();
            let asym = set.iter().filter(|&&(a, b)| !set.contains(&(b, a))).count();
            assert_eq!((mutual, asym), (7, 15));
            assert!(e.from.iter().zip(&e.to).all(|(a, b)| a != b));
        }
        assert!(uman(4, 4, 3, 1, 1, 1).is_err());
    }

    #[test]
    fn small_world_keeps_edge_count() {
        let ring = one(small_world(10, 4, 0.0, 1, 1, 1));
        assert_eq!(ring.len(), 20);
        assert!(ring.from.iter().zip(&ring.to).all(|(&a, &b)| {
            let d = (b as i64 - a as i64).rem_euclid(10);
            d == 1 || d == 2 || d == 8 || d == 9
        }));
        for e in small_world(30, 6, 0.3, 20, 3, 0).unwrap() {
            assert_eq!(e.len(), 90);
            assert!(is_simple(&e));
        }
        let full = one(small_world(5, 4, 1.0, 1, 1, 1));
        assert_eq!(full.len(), 10);
        assert!(small_world(10, 3, 0.1, 1, 1, 1).is_err());
        assert!(small_world(4, 4, 0.1, 1, 1, 1).is_err());
    }

    #[test]
    fn preferential_attachment_shape() {
        for e in preferential(50, 3, 10, 7, 0).unwrap() {
            assert_eq!(e.len(), 3 * 47);
            assert!(is_simple(&e));
            let mut deg = vec![0usize; 50];
            for (&a, &b) in e.from.iter().zip(&e.to) {
                deg[a as usize] += 1;
                deg[b as usize] += 1;
            }
            assert!(deg.iter().all(|&d| d >= 1));
            assert!(deg[4..].iter().all(|&d| d >= 3));
        }
        assert!(preferential(3, 3, 1, 1, 1).is_err());
    }

    #[test]
    fn sbm_blocks() {
        let e = one(sbm(&[3, 4], &[1.0, 0.0, 0.0, 1.0], false, false, 1, 1, 1));
        assert_eq!(e.len(), 3 + 6);
        assert!(e
            .from
            .iter()
            .zip(&e.to)
            .all(|(&a, &b)| (a < 3) == (b < 3) && a < b));
        let d = one(sbm(&[2, 2], &[0.0, 1.0, 0.0, 0.0], true, false, 1, 1, 1));
        assert_eq!(d.from, vec![0, 0, 1, 1]);
        assert_eq!(d.to, vec![2, 3, 2, 3]);
        assert!(sbm(&[2, 2], &[0.0, 1.0, 0.0, 0.0], false, false, 1, 1, 1).is_err());
        assert!(sbm(&[2], &[0.5, 0.5], false, false, 1, 1, 1).is_err());
    }
}
