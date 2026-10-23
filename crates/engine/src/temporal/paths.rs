//! Time-respecting paths: earliest arrival, latest departure, fastest and
//! shortest walks, temporal reachability, and temporal closeness.
//!
//! The walk model, the window, and the algorithms are described in the
//! documentation of the parent module ([`crate::temporal`]); sources are
//! Wu et al. (2014, `wu2014path`), Kempe, Kleinberg, and Kumar (2002,
//! `kempe2002connectivity`), Holme (2005, `holme2005network`), Holme and
//! Saramäki (2012, `holme2012temporal`), and Pan and Saramäki (2011,
//! `pan2011path`).
//!
//! Values reported for a walk from i to j:
//!
//! | Criterion | Value | i = j | No walk |
//! |---|---|---|---|
//! | Earliest arrival | arrival time of the last contact | window start | +inf |
//! | Latest departure | departure time of the first contact | window end | -inf |
//! | Fastest | arrival of the last minus departure of the first | 0 | +inf |
//! | Shortest | number of contacts | 0 | +inf |

use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

use super::{Stream, TemporalNetwork};
use crate::error::{Error, Result};
use crate::parallel;

/// What a time-respecting walk optimizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Criterion {
    /// The earliest arrival at the destination.
    EarliestArrival,
    /// The latest departure from the origin that still arrives in the window.
    LatestDeparture,
    /// The smallest elapsed time from departure to arrival.
    Fastest,
    /// The fewest contacts.
    Shortest,
}

impl Criterion {
    /// Decodes the codes used by the C interface: 0 = earliest arrival,
    /// 1 = latest departure, 2 = fastest, 3 = shortest.
    pub fn from_code(code: i32) -> Result<Criterion> {
        match code {
            0 => Ok(Criterion::EarliestArrival),
            1 => Ok(Criterion::LatestDeparture),
            2 => Ok(Criterion::Fastest),
            3 => Ok(Criterion::Shortest),
            _ => Err(Error::InvalidArgument(format!(
                "unknown temporal path criterion code {code}"
            ))),
        }
    }
}

/// Whether walks leave the roots or arrive at them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Walks from the root to every node.
    Out,
    /// Walks from every node to the root.
    In,
}

impl Direction {
    /// Decodes the codes used by the C interface: 1 = out, 2 = in (as for
    /// degree modes).
    pub fn from_code(code: i32) -> Result<Direction> {
        match code {
            1 => Ok(Direction::Out),
            2 => Ok(Direction::In),
            _ => Err(Error::InvalidArgument(format!(
                "unknown temporal direction code {code}; use 1 (out) or 2 (in)"
            ))),
        }
    }
}

/// How long traversing a contact takes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Traversal {
    /// The same time for every contact (zero by default).
    Constant(f64),
    /// The contact's own duration, end minus start.
    Duration,
}

impl Default for Traversal {
    fn default() -> Traversal {
        Traversal::Constant(0.0)
    }
}

/// Options for time-respecting walks.
#[derive(Debug, Clone, Copy)]
pub struct PathOptions {
    /// Start of the window; `None` uses the earliest start time.
    pub start: Option<f64>,
    /// End of the window; `None` uses the latest arrival of any contact.
    pub end: Option<f64>,
    /// Smallest admissible waiting time at an intermediate node.
    pub min_wait: f64,
    /// Largest admissible waiting time (`f64::INFINITY` for no limit).
    pub max_wait: f64,
    /// Require waiting times strictly greater than `min_wait`.
    pub strict: bool,
    /// Traversal time of a contact.
    pub traversal: Traversal,
    /// Worker threads; 0 uses all cores.
    pub threads: usize,
}

impl Default for PathOptions {
    fn default() -> PathOptions {
        PathOptions {
            start: None,
            end: None,
            min_wait: 0.0,
            max_wait: f64::INFINITY,
            strict: false,
            traversal: Traversal::default(),
            threads: 0,
        }
    }
}

/// The window and waiting rule after validation, in the orientation of
/// the stream being swept (negated and swapped for reversed streams).
#[derive(Debug, Clone, Copy)]
struct Rule {
    lo: f64,
    hi: f64,
    min_wait: f64,
    max_wait: f64,
    strict: bool,
}

impl Rule {
    #[inline]
    fn usable(&self, wait: f64) -> bool {
        if self.strict {
            wait > self.min_wait
        } else {
            wait >= self.min_wait
        }
    }

    #[inline]
    fn expired(&self, wait: f64) -> bool {
        wait > self.max_wait
    }

    /// Whether an arrival can be extended at the very time it arrives,
    /// which lets contacts with equal times chain.
    #[inline]
    fn zero_wait(&self) -> bool {
        !self.strict && self.min_wait == 0.0
    }

    fn reversed(&self) -> Rule {
        Rule {
            lo: -self.hi,
            hi: -self.lo,
            ..*self
        }
    }
}

/// Validates options and resolves the window to finite bounds.
fn resolve(net: &TemporalNetwork, opts: &PathOptions) -> Result<Rule> {
    if let Traversal::Constant(x) = opts.traversal {
        if !(x.is_finite() && x >= 0.0) {
            return Err(Error::InvalidArgument(
                "the traversal time must be finite and non-negative".into(),
            ));
        }
    }
    if !(opts.min_wait.is_finite() && opts.min_wait >= 0.0) {
        return Err(Error::InvalidArgument(
            "the minimum waiting time must be finite and non-negative".into(),
        ));
    }
    if opts.max_wait.is_nan() || opts.max_wait < opts.min_wait {
        return Err(Error::InvalidArgument(
            "the maximum waiting time must be at least the minimum waiting time".into(),
        ));
    }
    let lambda = |c: usize| match opts.traversal {
        Traversal::Constant(x) => x,
        Traversal::Duration => net.end[c] - net.start[c],
    };
    let lo = match opts.start {
        Some(t) => t,
        None => net.first_time().unwrap_or(0.0),
    };
    let hi = match opts.end {
        Some(t) => t,
        None => (0..net.contact_count())
            .map(|c| net.start[c] + lambda(c))
            .reduce(f64::max)
            .unwrap_or(lo),
    };
    if lo.is_nan() || hi.is_nan() || lo == f64::INFINITY || hi == f64::NEG_INFINITY {
        return Err(Error::InvalidArgument(
            "the window must have a start and end that are numbers".into(),
        ));
    }
    if lo > hi {
        return Err(Error::InvalidArgument(
            "the window starts after it ends".into(),
        ));
    }
    Ok(Rule {
        lo,
        hi,
        min_wait: opts.min_wait,
        max_wait: opts.max_wait,
        strict: opts.strict,
    })
}

/// What labels arrivals carry. Keys are maximized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// Reachability only (key always 0); value = earliest arrival.
    Reach,
    /// Key = departure time of the walk's first contact; value = arrival
    /// minus key.
    Start,
    /// Key = minus the number of contacts; value = minus key.
    Hops,
}

impl Kind {
    #[inline]
    fn fresh(self, dep: f64) -> f64 {
        match self {
            Kind::Reach => 0.0,
            Kind::Start => dep,
            Kind::Hops => -1.0,
        }
    }

    #[inline]
    fn extend(self, key: f64) -> f64 {
        match self {
            Kind::Hops => key - 1.0,
            _ => key,
        }
    }

    #[inline]
    fn value(self, arr: f64, key: f64) -> f64 {
        match self {
            Kind::Reach => arr,
            Kind::Start => arr - key,
            Kind::Hops => -key,
        }
    }
}

/// An arrival waiting until it may be extended.
#[derive(Debug, Clone, Copy)]
struct Pending {
    arr: f64,
    node: u32,
    key: f64,
}

impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Pending {}
impl PartialOrd for Pending {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Pending {
    // Reversed, so the standard max-heap pops the earliest arrival first.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .arr
            .total_cmp(&self.arr)
            .then(other.node.cmp(&self.node))
    }
}

/// Arrivals not yet usable, in order of arrival time. When the stream's
/// arrivals are non-decreasing a FIFO queue suffices.
#[derive(Debug)]
enum PendingQueue {
    Fifo(VecDeque<Pending>),
    Heap(BinaryHeap<Pending>),
}

impl PendingQueue {
    fn clear(&mut self) {
        match self {
            PendingQueue::Fifo(q) => q.clear(),
            PendingQueue::Heap(h) => h.clear(),
        }
    }

    fn push(&mut self, p: Pending) {
        match self {
            PendingQueue::Fifo(q) => q.push_back(p),
            PendingQueue::Heap(h) => h.push(p),
        }
    }

    fn peek(&self) -> Option<&Pending> {
        match self {
            PendingQueue::Fifo(q) => q.front(),
            PendingQueue::Heap(h) => h.peek(),
        }
    }

    fn pop(&mut self) -> Option<Pending> {
        match self {
            PendingQueue::Fifo(q) => q.pop_front(),
            PendingQueue::Heap(h) => h.pop(),
        }
    }
}

/// Per-thread scratch space for sweeps from one root at a time.
struct Workspace {
    /// Best usable key per node without a waiting limit.
    best: Vec<f64>,
    /// Usable (arrival, key) pairs per node with a waiting limit: arrivals
    /// increase from front to back and keys decrease, so the front is the
    /// best unexpired label (a monotone queue).
    window: Vec<VecDeque<(f64, f64)>>,
    touched: Vec<u32>,
    is_touched: Vec<bool>,
    pending: PendingQueue,
    /// Nodes whose label improved through a zero-wait arrival inside the
    /// current group of equal departure times.
    group: BinaryHeap<(OrderedKey, u32)>,
    stack: Vec<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderedKey(f64);
impl Eq for OrderedKey {}
impl PartialOrd for OrderedKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrderedKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl Workspace {
    fn new(n: usize, limited: bool, fifo: bool) -> Workspace {
        Workspace {
            best: vec![f64::NEG_INFINITY; n],
            window: if limited {
                (0..n).map(|_| VecDeque::new()).collect()
            } else {
                Vec::new()
            },
            touched: Vec::new(),
            is_touched: vec![false; n],
            pending: if fifo {
                PendingQueue::Fifo(VecDeque::new())
            } else {
                PendingQueue::Heap(BinaryHeap::new())
            },
            group: BinaryHeap::new(),
            stack: Vec::new(),
        }
    }

    fn reset(&mut self) {
        for &v in &self.touched {
            let v = v as usize;
            self.best[v] = f64::NEG_INFINITY;
            if let Some(w) = self.window.get_mut(v) {
                w.clear();
            }
            self.is_touched[v] = false;
        }
        self.touched.clear();
        self.pending.clear();
        self.group.clear();
        self.stack.clear();
    }

    fn touch(&mut self, v: usize) {
        if !self.is_touched[v] {
            self.is_touched[v] = true;
            self.touched.push(v as u32);
        }
    }

    /// The best label usable at `v` at time `t`, if any.
    #[inline]
    fn query(&mut self, v: usize, t: f64, rule: &Rule) -> Option<f64> {
        if self.window.is_empty() {
            let b = self.best[v];
            (b > f64::NEG_INFINITY).then_some(b)
        } else {
            let w = &mut self.window[v];
            while let Some(&(a, _)) = w.front() {
                if rule.expired(t - a) {
                    w.pop_front();
                } else {
                    break;
                }
            }
            w.front().map(|&(_, k)| k)
        }
    }

    /// Makes an arrival usable. Arrivals are inserted in non-decreasing
    /// order of arrival time. Returns whether the node's best label at
    /// that time improved.
    #[inline]
    fn insert(&mut self, v: usize, arr: f64, key: f64) -> bool {
        self.touch(v);
        if self.window.is_empty() {
            if key > self.best[v] {
                self.best[v] = key;
                true
            } else {
                false
            }
        } else {
            let w = &mut self.window[v];
            let before = w.front().map(|&(_, k)| k);
            while let Some(&(_, k)) = w.back() {
                if k <= key {
                    w.pop_back();
                } else {
                    break;
                }
            }
            w.push_back((arr, key));
            before.map_or(true, |b| key > b)
        }
    }
}

/// Positions in `stream` of the arcs in `range` whose tail is `v` (the
/// stream is sorted by tail within equal departure times).
fn tail_range(stream: &Stream, range: (usize, usize), v: u32) -> (usize, usize) {
    let slice = &stream.from[range.0..range.1];
    let a = slice.partition_point(|&x| x < v);
    let b = slice.partition_point(|&x| x <= v);
    (range.0 + a, range.0 + b)
}

/// One sweep from `root` over a stream, writing the value for every node
/// into `out` (initialized here). For [`Kind::Start`], `pairs` (when given)
/// receives every (node, start, arrival) produced, from which the
/// fastest-path profile follows.
fn sweep(
    stream: &Stream,
    root: u32,
    rule: &Rule,
    kind: Kind,
    ws: &mut Workspace,
    out: &mut [f64],
    mut pairs: Option<&mut Vec<(u32, f64, f64)>>,
) {
    ws.reset();
    out.iter_mut().for_each(|x| *x = f64::INFINITY);
    out[root as usize] = match kind {
        Kind::Reach => rule.lo,
        _ => 0.0,
    };
    let k = stream.len();
    let mut i = stream.dep.partition_point(|&d| d < rule.lo);
    while i < k {
        let t = stream.dep[i];
        if t > rule.hi {
            break;
        }
        let mut j = i + 1;
        while j < k && stream.dep[j] == t {
            j += 1;
        }
        // Arrivals that have waited long enough become usable.
        while let Some(p) = ws.pending.peek() {
            if rule.usable(t - p.arr) {
                let p = ws.pending.pop().expect("peeked");
                ws.insert(p.node as usize, p.arr, p.key);
            } else {
                break;
            }
        }
        let chain = rule.zero_wait();
        for a in i..j {
            relax(
                stream, a, root, rule, kind, ws, out, &mut pairs, chain, None,
            );
        }
        if chain {
            // Labels that improved within the group are pushed again
            // through the group's arcs, best label first, so each node is
            // expanded once per improvement (labels never improve along
            // an arc).
            while let Some((OrderedKey(key), v)) = ws.group.pop() {
                if ws.query(v as usize, t, rule) != Some(key) {
                    continue;
                }
                let (a0, a1) = tail_range(stream, (i, j), v);
                for a in a0..a1 {
                    relax(
                        stream,
                        a,
                        root,
                        rule,
                        kind,
                        ws,
                        out,
                        &mut pairs,
                        chain,
                        Some(key),
                    );
                }
            }
        }
        i = j;
    }
}

/// Extends the best label usable at the tail of arc `a` (or `forced`)
/// along the arc.
#[allow(clippy::too_many_arguments)]
#[inline]
fn relax(
    stream: &Stream,
    a: usize,
    root: u32,
    rule: &Rule,
    kind: Kind,
    ws: &mut Workspace,
    out: &mut [f64],
    pairs: &mut Option<&mut Vec<(u32, f64, f64)>>,
    chain: bool,
    forced: Option<f64>,
) {
    let (u, v, t, arr) = (stream.from[a], stream.to[a], stream.dep[a], stream.arr[a]);
    if arr > rule.hi || v == root {
        return;
    }
    let key = if u == root {
        kind.fresh(t)
    } else {
        match forced.or_else(|| ws.query(u as usize, t, rule)) {
            Some(k) => kind.extend(k),
            None => return,
        }
    };
    let value = kind.value(arr, key);
    let slot = &mut out[v as usize];
    if value < *slot {
        *slot = value;
    }
    if let Some(p) = pairs.as_deref_mut() {
        p.push((v, key, arr));
    }
    if chain && arr == t {
        if ws.insert(v as usize, arr, key) {
            ws.group.push((OrderedKey(key), v));
        }
    } else {
        ws.pending.push(Pending { arr, node: v, key });
    }
}

/// Earliest arrival without a waiting limit: the earliest arrival at a
/// node is always its most extensible one, so one time per node suffices
/// and no queue is needed. O(m) per root.
fn sweep_earliest(stream: &Stream, root: u32, rule: &Rule, ws: &mut Workspace, out: &mut [f64]) {
    ws.stack.clear();
    out.iter_mut().for_each(|x| *x = f64::INFINITY);
    let r = root as usize;
    let k = stream.len();
    let chain = rule.zero_wait();
    let usable = |out: &[f64], u: usize, t: f64| u == r || rule.usable(t - out[u]);
    let mut i = stream.dep.partition_point(|&d| d < rule.lo);
    while i < k {
        let t = stream.dep[i];
        if t > rule.hi {
            break;
        }
        let mut j = i + 1;
        while j < k && stream.dep[j] == t {
            j += 1;
        }
        for a in i..j {
            let (u, v, arr) = (
                stream.from[a] as usize,
                stream.to[a] as usize,
                stream.arr[a],
            );
            if arr <= rule.hi && usable(out, u, t) && arr < out[v] {
                out[v] = arr;
                if chain && arr == t {
                    ws.stack.push(v as u32);
                }
            }
        }
        while let Some(w) = ws.stack.pop() {
            let (a0, a1) = tail_range(stream, (i, j), w);
            for a in a0..a1 {
                let (v, arr) = (stream.to[a] as usize, stream.arr[a]);
                if arr <= rule.hi && arr < out[v] {
                    out[v] = arr;
                    if arr == t {
                        ws.stack.push(v as u32);
                    }
                }
            }
        }
        i = j;
    }
    out[r] = rule.lo;
}

/// A stream and the rule in its orientation, ready for sweeps.
struct Prepared {
    stream: Stream,
    rule: Rule,
    fifo: bool,
}

fn prepare(
    net: &TemporalNetwork,
    direction: Direction,
    opts: &PathOptions,
) -> Result<(Prepared, Rule)> {
    let rule = resolve(net, opts)?;
    let stream = net.stream(direction, opts.traversal);
    // Arrivals within a group of equal departures can be queued in any
    // order, so a FIFO queue is valid only when they are all equal there
    // and non-decreasing along the stream (always so for a constant
    // traversal time, barring rounding).
    let fifo = (1..stream.len()).all(|i| {
        let (a, b) = (stream.arr[i - 1], stream.arr[i]);
        if stream.dep[i - 1] == stream.dep[i] {
            a == b
        } else {
            a <= b
        }
    });
    let oriented = match direction {
        Direction::Out => rule,
        Direction::In => rule.reversed(),
    };
    Ok((
        Prepared {
            stream,
            rule: oriented,
            fifo,
        },
        rule,
    ))
}

fn check_roots(net: &TemporalNetwork, roots: &[usize]) -> Result<()> {
    if let Some(&r) = roots.iter().find(|&&r| r >= net.node_count()) {
        return Err(Error::InvalidArgument(format!(
            "node index {r} is out of range; valid indices are 0..{}",
            net.node_count()
        )));
    }
    Ok(())
}

/// Runs one sweep per root in parallel; returns rows in root order.
fn sweep_all(p: &Prepared, roots: &[usize], n: usize, kind: Kind, threads: usize) -> Vec<f64> {
    let limited = p.rule.max_wait.is_finite();
    let fast = kind == Kind::Reach && !limited;
    let states = parallel::for_each_item(
        roots.len(),
        threads,
        || {
            (
                Workspace::new(n, limited, p.fifo),
                Vec::<(usize, Vec<f64>)>::new(),
            )
        },
        |(ws, rows), r| {
            let mut row = vec![0.0; n];
            if fast {
                sweep_earliest(&p.stream, roots[r] as u32, &p.rule, ws, &mut row);
            } else {
                sweep(
                    &p.stream,
                    roots[r] as u32,
                    &p.rule,
                    kind,
                    ws,
                    &mut row,
                    None,
                );
            }
            rows.push((r, row));
        },
    );
    let mut out = vec![0.0; roots.len() * n];
    for (_, rows) in states {
        for (r, row) in rows {
            out[r * n..(r + 1) * n].copy_from_slice(&row);
        }
    }
    out
}

/// Optimal time-respecting walks between `roots` and every node.
///
/// Returns a row-major matrix with one row per root and one column per
/// node. With [`Direction::Out`] entry (r, j) describes walks from
/// `roots[r]` to j; with [`Direction::In`], walks from j to `roots[r]`.
/// Each row costs one pass over the contacts, except that earliest
/// arrival into roots and latest departure out of roots are not
/// single-pass problems (Wu et al. 2014): they are computed with one pass
/// per node of the network and cost n passes whatever the number of
/// roots.
pub fn temporal_paths(
    net: &TemporalNetwork,
    criterion: Criterion,
    direction: Direction,
    roots: &[usize],
    opts: &PathOptions,
) -> Result<Vec<f64>> {
    check_roots(net, roots)?;
    let n = net.node_count();
    let natural = match criterion {
        Criterion::EarliestArrival => Direction::Out,
        Criterion::LatestDeparture => Direction::In,
        _ => direction,
    };
    let kind = match criterion {
        Criterion::EarliestArrival | Criterion::LatestDeparture => Kind::Reach,
        Criterion::Fastest => Kind::Start,
        Criterion::Shortest => Kind::Hops,
    };
    let (p, _) = prepare(net, natural, opts)?;
    let negate = natural == Direction::In && kind == Kind::Reach;
    if natural == direction {
        let mut out = sweep_all(&p, roots, n, kind, opts.threads);
        if negate {
            out.iter_mut().for_each(|x| *x = -*x);
        }
        return Ok(out);
    }
    // Transpose a full all-pairs computation in the natural direction.
    let all: Vec<usize> = (0..n).collect();
    let mut full = sweep_all(&p, &all, n, kind, opts.threads);
    if negate {
        full.iter_mut().for_each(|x| *x = -*x);
    }
    let mut out = vec![0.0; roots.len() * n];
    for (r, &root) in roots.iter().enumerate() {
        for j in 0..n {
            out[r * n + j] = full[j * n + root];
        }
    }
    Ok(out)
}

/// For every node, the number of other nodes it reaches by time-respecting
/// walks in the window ([`Direction::Out`], the size of its set of
/// influence) or that reach it ([`Direction::In`], its source set), after
/// Holme (2005). Costs one pass over the contacts per node.
pub fn temporal_reach_counts(
    net: &TemporalNetwork,
    direction: Direction,
    opts: &PathOptions,
) -> Result<Vec<f64>> {
    let n = net.node_count();
    let (p, _) = prepare(net, direction, opts)?;
    let limited = p.rule.max_wait.is_finite();
    let states = parallel::for_each_item(
        n,
        opts.threads,
        || {
            (
                Workspace::new(n, limited, p.fifo),
                vec![0.0; n],
                Vec::<(usize, f64)>::new(),
            )
        },
        |(ws, row, counts), r| {
            if limited {
                sweep(&p.stream, r as u32, &p.rule, Kind::Reach, ws, row, None);
            } else {
                sweep_earliest(&p.stream, r as u32, &p.rule, ws, row);
            }
            let c = row.iter().filter(|x| x.is_finite()).count() - 1;
            counts.push((r, c as f64));
        },
    );
    let mut out = vec![0.0; n];
    for (_, _, counts) in states {
        for (r, c) in counts {
            out[r] = c;
        }
    }
    Ok(out)
}

/// Pan and Saramäki's (2011) temporal closeness of every node, with the
/// average temporal distances it is built from.
#[derive(Debug, Clone)]
pub struct TemporalCloseness {
    /// C_i = (1 / (n - 1)) sum_{j != i} 1 / τ_ij; NaN for n < 2 or an
    /// empty window.
    pub closeness: Vec<f64>,
    /// τ_ij row-major (i the origin), +inf where no walk exists, 0 on the
    /// diagonal; empty unless requested.
    pub distances: Vec<f64>,
}

/// The fastest-path profile of one target: the walks (departure t_k,
/// arrival a_k) such that leaving at t_k arrives at a_k, and no walk
/// leaving at or after t_k arrives earlier, sorted by departure. Both
/// coordinates increase strictly.
fn profile(pairs: &mut [(f64, f64)]) -> Vec<(f64, f64)> {
    pairs.sort_unstable_by(|x, y| y.0.total_cmp(&x.0).then(x.1.total_cmp(&y.1)));
    let mut front = Vec::new();
    let mut best = f64::INFINITY;
    for &(s, a) in pairs.iter() {
        if a < best {
            front.push((s, a));
            best = a;
        }
    }
    front.reverse();
    front
}

/// Average temporal distance over starting times in [lo, hi] with Pan and
/// Saramäki's periodic boundary condition (their Eq. 1): a start after
/// the last walk waits for the first walk of the next period.
fn average_latency(front: &[(f64, f64)], lo: f64, hi: f64) -> f64 {
    if front.is_empty() {
        return f64::INFINITY;
    }
    let span = hi - lo;
    let mut total = 0.0;
    let mut prev = lo;
    for &(t, a) in front {
        let gap = t - prev;
        total += gap * (gap / 2.0 + (a - t));
        prev = t;
    }
    let (t1, a1) = front[0];
    let gap = hi - prev;
    total += gap * (gap / 2.0 + (t1 - lo) + (a1 - t1));
    total / span
}

/// Temporal closeness after Pan and Saramäki (2011, Eqs. 1 and 2):
/// C_i = (1 / (n - 1)) sum_{j != i} 1 / τ_ij, where τ_ij is the time from
/// a starting time t to the earliest arrival at j of a walk from i leaving
/// at or after t, averaged over t in the window, with walks repeating
/// periodically beyond its end. Walks follow `opts`; Pan and Saramäki's
/// own model is `strict` waiting with [`Traversal::Duration`]. The window
/// must have positive length.
pub fn temporal_closeness(
    net: &TemporalNetwork,
    opts: &PathOptions,
    with_distances: bool,
) -> Result<TemporalCloseness> {
    let n = net.node_count();
    let (p, rule) = prepare(net, Direction::Out, opts)?;
    if n < 2 || rule.hi <= rule.lo {
        return Ok(TemporalCloseness {
            closeness: vec![f64::NAN; n],
            distances: if with_distances {
                let mut d = vec![f64::NAN; n * n];
                (0..n).for_each(|i| d[i * n + i] = 0.0);
                d
            } else {
                Vec::new()
            },
        });
    }
    let limited = p.rule.max_wait.is_finite();
    let states = parallel::for_each_item(
        n,
        opts.threads,
        || {
            (
                Workspace::new(n, limited, p.fifo),
                vec![0.0; n],
                Vec::new(),
                Vec::<(usize, f64, Vec<f64>)>::new(),
            )
        },
        |(ws, row, pairs, results), r| {
            pairs.clear();
            sweep(
                &p.stream,
                r as u32,
                &p.rule,
                Kind::Start,
                ws,
                row,
                Some(pairs),
            );
            pairs.sort_unstable_by_key(|x: &(u32, f64, f64)| x.0);
            let mut tau = vec![f64::INFINITY; n];
            tau[r] = 0.0;
            let mut a = 0;
            while a < pairs.len() {
                let v = pairs[a].0;
                let mut b = a;
                while b < pairs.len() && pairs[b].0 == v {
                    b += 1;
                }
                let mut sa: Vec<(f64, f64)> = pairs[a..b].iter().map(|x| (x.1, x.2)).collect();
                tau[v as usize] = average_latency(&profile(&mut sa), rule.lo, rule.hi);
                a = b;
            }
            let sum: f64 = (0..n).filter(|&j| j != r).map(|j| 1.0 / tau[j]).sum();
            results.push((
                r,
                sum / (n - 1) as f64,
                if with_distances { tau } else { Vec::new() },
            ));
        },
    );
    let mut closeness = vec![0.0; n];
    let mut distances = if with_distances {
        vec![0.0; n * n]
    } else {
        Vec::new()
    };
    for (_, _, _, results) in states {
        for (r, c, tau) in results {
            closeness[r] = c;
            if with_distances {
                distances[r * n..(r + 1) * n].copy_from_slice(&tau);
            }
        }
    }
    Ok(TemporalCloseness {
        closeness,
        distances,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn net(
        src: &[i64],
        tgt: &[i64],
        start: &[f64],
        end: Option<&[f64]>,
        d: bool,
    ) -> TemporalNetwork {
        let n = src
            .iter()
            .chain(tgt)
            .copied()
            .max()
            .map_or(0, |x| x as usize + 1);
        TemporalNetwork::new(n, src, tgt, start, end, d).unwrap()
    }

    fn row(
        net: &TemporalNetwork,
        c: Criterion,
        dir: Direction,
        root: usize,
        o: &PathOptions,
    ) -> Vec<f64> {
        temporal_paths(net, c, dir, &[root], o).unwrap()
    }

    #[test]
    fn earliest_arrival_respects_time_order() {
        // 0 -> 1 at 1, 1 -> 2 at 2, 2 -> 3 at 1 (too early), 0 -> 3 at 5.
        let g = net(
            &[0, 1, 2, 0],
            &[1, 2, 3, 3],
            &[1.0, 2.0, 1.0, 5.0],
            None,
            true,
        );
        let o = PathOptions::default();
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![1.0, 1.0, 2.0, 5.0]);
        let sh = row(&g, Criterion::Shortest, Direction::Out, 0, &o);
        assert_eq!(sh, vec![0.0, 1.0, 2.0, 1.0]);
        let fa = row(&g, Criterion::Fastest, Direction::Out, 0, &o);
        assert_eq!(fa, vec![0.0, 0.0, 1.0, 0.0]);
        // Latest departure from each node to 3 within [1, 5].
        let ld = row(&g, Criterion::LatestDeparture, Direction::In, 3, &o);
        assert_eq!(ld, vec![5.0, f64::NEG_INFINITY, 1.0, 5.0]);
    }

    #[test]
    fn equal_times_chain_unless_strict() {
        // Listed so that the second hop precedes the first within the time.
        let g = net(&[1, 0], &[2, 1], &[3.0, 3.0], None, true);
        let mut o = PathOptions::default();
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![3.0, 3.0, 3.0]);
        let sh = row(&g, Criterion::Shortest, Direction::Out, 0, &o);
        assert_eq!(sh, vec![0.0, 1.0, 2.0]);
        o.strict = true;
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![3.0, 3.0, f64::INFINITY]);
        let sh = row(&g, Criterion::Shortest, Direction::Out, 0, &o);
        assert_eq!(sh, vec![0.0, 1.0, f64::INFINITY]);
    }

    #[test]
    fn waiting_limits_and_traversal_times() {
        // 0 -> 1 at 0, 1 -> 2 at 5, 1 -> 2 at 2.
        let g = net(
            &[0, 1, 1],
            &[1, 2, 2],
            &[0.0, 5.0, 2.0],
            Some(&[1.0, 5.0, 2.0]),
            true,
        );
        let mut o = PathOptions {
            max_wait: 3.0,
            ..Default::default()
        };
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![0.0, 0.0, 2.0]);
        o.min_wait = 3.0;
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![0.0, 0.0, f64::INFINITY]);
        o.min_wait = 0.0;
        o.max_wait = f64::INFINITY;
        o.traversal = Traversal::Duration;
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![0.0, 1.0, 2.0]);
        o.traversal = Traversal::Constant(2.5);
        let ea = row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o);
        assert_eq!(ea, vec![0.0, 2.5, 7.5]);
    }

    #[test]
    fn fastest_prefers_late_starts() {
        // Leave at 0 and arrive at 10, or leave at 8 and arrive at 9.
        let g = net(
            &[0, 1, 0, 1],
            &[1, 2, 1, 2],
            &[0.0, 10.0, 8.0, 9.0],
            None,
            true,
        );
        let o = PathOptions::default();
        assert_eq!(row(&g, Criterion::Fastest, Direction::Out, 0, &o)[2], 1.0);
        assert_eq!(row(&g, Criterion::Fastest, Direction::In, 2, &o)[0], 1.0);
        assert_eq!(
            row(&g, Criterion::EarliestArrival, Direction::Out, 0, &o)[2],
            9.0
        );
        let ea_in = row(&g, Criterion::EarliestArrival, Direction::In, 2, &o);
        assert_eq!(ea_in[0], 9.0);
        let ld_out = row(&g, Criterion::LatestDeparture, Direction::Out, 0, &o);
        assert_eq!(ld_out[2], 8.0);
    }

    #[test]
    fn reach_counts_and_undirected_contacts() {
        let g = net(&[0, 1], &[1, 2], &[2.0, 1.0], None, false);
        let o = PathOptions::default();
        assert_eq!(
            temporal_reach_counts(&g, Direction::Out, &o).unwrap(),
            vec![1.0, 2.0, 2.0]
        );
        assert_eq!(
            temporal_reach_counts(&g, Direction::In, &o).unwrap(),
            vec![2.0, 2.0, 1.0]
        );
    }

    #[test]
    fn closeness_matches_pan_saramaki_single_path() {
        // One contact 0 -> 1 at time 4 in the window [0, 10]: the average
        // latency is T/2 + Δt with the periodic condition, here 5 + 0.
        let g = net(&[0], &[1], &[4.0], None, true);
        let o = PathOptions {
            start: Some(0.0),
            end: Some(10.0),
            ..Default::default()
        };
        let c = temporal_closeness(&g, &o, true).unwrap();
        assert_eq!(c.distances, vec![0.0, 5.0, f64::INFINITY, 0.0]);
        assert_eq!(c.closeness, vec![0.2, 0.0]);
    }

    #[test]
    fn rejects_bad_options() {
        let g = net(&[0], &[1], &[0.0], None, true);
        let bad = [
            PathOptions {
                min_wait: -1.0,
                ..Default::default()
            },
            PathOptions {
                min_wait: 2.0,
                max_wait: 1.0,
                ..Default::default()
            },
            PathOptions {
                start: Some(2.0),
                end: Some(1.0),
                ..Default::default()
            },
            PathOptions {
                traversal: Traversal::Constant(f64::NAN),
                ..Default::default()
            },
        ];
        for o in bad {
            assert!(temporal_paths(&g, Criterion::Fastest, Direction::Out, &[0], &o).is_err());
        }
        assert!(temporal_paths(
            &g,
            Criterion::Fastest,
            Direction::Out,
            &[5],
            &PathOptions::default()
        )
        .is_err());
    }
}
