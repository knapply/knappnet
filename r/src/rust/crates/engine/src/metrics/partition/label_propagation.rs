//! Label propagation communities.
//!
//! Raghavan, Albert, and Kumara (2007, `raghavan2007near`) give every node
//! its own label and let nodes, visited asynchronously in a random order,
//! adopt the label that most of their neighbors carry (ties broken
//! uniformly at random), until every node carries a label of maximum
//! frequency among its neighbors. The variants, after those igraph
//! documents for `cluster_label_prop()`:
//!
//! - `Dominance` (igraph's default): each sweep visits the nodes in a new
//!   random order and every node draws uniformly among its dominant
//!   labels; the method stops after a sweep at whose end every node
//!   carries a dominant label (the stopping rule of Raghavan et al.).
//! - `Retention`: a node keeps its label when it is dominant and otherwise
//!   draws uniformly among the dominant labels; the method stops after a
//!   sweep that changes no label.
//! - `Fast`: the fast label propagation of Traag and Subelj (2023,
//!   `traag2023large`). Nodes wait in a queue, initially every node in a
//!   random order; the node at its head draws uniformly among its dominant
//!   labels, and when its label changes, the nodes that read its label,
//!   are not queued, and carry a different label join the end of the
//!   queue. The method stops when the queue is empty.
//! - `SemiSynchronous`: the semi-synchronous method of Cordasco and
//!   Gargano (2010, `cordasco2010community`). The nodes are coloured
//!   greedily (largest degree first, ties by node index, smallest free
//!   colour), and each round updates the colour classes in colour order;
//!   the nodes of a class are not adjacent, so updating them together is
//!   the same as one at a time. A node takes its unique dominant label,
//!   and among several keeps its own if dominant and otherwise takes the
//!   largest ("Prec-Max"). Rounds repeat until every node carries a
//!   dominant label. It uses no random numbers, converges on every
//!   undirected network (Cordasco and Gargano 2010), and matches
//!   NetworkX's `label_propagation_communities()`.
//!
//! Conventions:
//! - The labels a node reads: with `mode` `Out`, labels propagate along
//!   edge directions, so a node reads its in-neighbors; with `In`, against
//!   them (it reads its out-neighbors); with `All`, directions are ignored.
//!   On undirected networks `mode` has no effect. As igraph's
//!   documentation advises, `All` is the natural choice for directed
//!   networks; the semi-synchronous variant needs an undirected network.
//! - Each neighbor counts with the weight of the edge (1 without weights;
//!   weights must be non-negative); parallel edges add up. Self-loops are
//!   ignored, so a node never reads its own label (NetworkX counts a loop
//!   as one neighbor, and igraph's results on networks with loops are
//!   consistent with the same rule, which it does not document).
//! - `initial` labels: non-negative values are labels, negative values
//!   mark unlabeled nodes, which count for nothing until they are
//!   labeled. Without `initial` every node starts with its own label.
//!   `fixed` nodes (which must be labeled) never change. A node that reads
//!   no labeled neighbor keeps its label. Nodes that are still unlabeled
//!   at the end (no labeled node reaches them) each form a group of their
//!   own.
//! - Random numbers come from [`crate::rng::Stream`] with the seed (sweep
//!   s uses stream s; the fast variant uses stream 0), so a seed
//!   reproduces the result on any platform. The method runs on one thread.
//! - A node carries a dominant label when no label has a larger total
//!   weight among the neighbors it reads (exact comparison of the totals).
//! - The randomized variants stop after at most 100,000 sweeps (or
//!   100,000 n queue steps) even if the stopping rule has not been met,
//!   which can happen when ties make labels oscillate.
//! - Groups are labelled 0, 1, ... in order of their first node.

use std::collections::VecDeque;

use super::relabel;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::rng::Stream;

/// The update and stopping rules of [`label_propagation`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LpaVariant {
    /// Draw among dominant labels; stop when every label is dominant.
    Dominance,
    /// Keep a dominant label; stop when a sweep changes nothing.
    Retention,
    /// Queue-based fast label propagation (Traag and Subelj 2023).
    Fast,
    /// Semi-synchronous propagation (Cordasco and Gargano 2010).
    SemiSynchronous,
}

impl LpaVariant {
    /// 0 = dominance, 1 = retention, 2 = fast, 3 = semi-synchronous.
    pub fn from_code(code: i32) -> Result<LpaVariant> {
        Ok(match code {
            0 => LpaVariant::Dominance,
            1 => LpaVariant::Retention,
            2 => LpaVariant::Fast,
            3 => LpaVariant::SemiSynchronous,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown label propagation variant code {code}"
                )))
            }
        })
    }
}

/// Direction in which labels propagate on a directed network.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LpaMode {
    /// Along edge directions (a node reads its in-neighbors).
    Out,
    /// Against edge directions (a node reads its out-neighbors).
    In,
    /// Directions ignored.
    All,
}

impl LpaMode {
    /// 1 = out, 2 = in, 3 = all (as for `DegreeMode`).
    pub fn from_code(code: i32) -> Result<LpaMode> {
        Ok(match code {
            1 => LpaMode::Out,
            2 => LpaMode::In,
            3 => LpaMode::All,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown label propagation mode code {code}"
                )))
            }
        })
    }
}

/// Options for [`label_propagation`].
#[derive(Debug, Clone, Copy)]
pub struct LabelPropagationOptions<'a> {
    /// Update and stopping rules.
    pub variant: LpaVariant,
    /// Direction of propagation on directed networks.
    pub mode: LpaMode,
    /// Count neighbors by edge weight.
    pub weighted: bool,
    /// Initial labels (negative = unlabeled); `None` gives every node its
    /// own label.
    pub initial: Option<&'a [i64]>,
    /// Nodes whose label never changes.
    pub fixed: Option<&'a [bool]>,
    /// Seed of the random stream.
    pub seed: u64,
}

const UNLABELED: usize = usize::MAX;
const MAX_SWEEPS: usize = 100_000;

/// Who reads whom: `reads[v]` are the (neighbor, weight) pairs whose labels
/// v counts, and `readers[v]` the nodes that count v's label.
struct Links {
    reads: Vec<Vec<(usize, f64)>>,
    readers: Vec<Vec<usize>>,
}

fn links(g: &Graph, mode: LpaMode, weighted: bool) -> Result<Links> {
    crate::weights::additive_weights(g, weighted, "label propagation")?;
    let n = g.node_count();
    let use_w = weighted && g.is_weighted();
    let mut reads: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut readers: Vec<Vec<usize>> = vec![Vec::new(); n];
    let add = |reads: &mut Vec<Vec<(usize, f64)>>, csr: &crate::graph::Csr, v: usize| {
        let w = csr.weights(v);
        for (i, &u) in csr.neighbors(v).iter().enumerate() {
            let u = u as usize;
            if u != v {
                let x = if use_w { w.map_or(1.0, |w| w[i]) } else { 1.0 };
                reads[v].push((u, x));
            }
        }
    };
    for v in 0..n {
        if !g.is_directed() {
            add(&mut reads, g.out(), v);
        } else {
            if mode != LpaMode::In {
                add(&mut reads, g.inc(), v);
            }
            if mode != LpaMode::Out {
                add(&mut reads, g.out(), v);
            }
        }
    }
    for (v, row) in reads.iter().enumerate() {
        for &(u, _) in row {
            readers[u].push(v);
        }
    }
    for r in readers.iter_mut() {
        r.sort_unstable();
        r.dedup();
    }
    Ok(Links { reads, readers })
}

/// Scratch space for counting labels.
struct Counter {
    count: Vec<f64>,
    seen: Vec<bool>,
    touched: Vec<usize>,
    dom: Vec<usize>,
}

impl Counter {
    fn new(n: usize) -> Counter {
        Counter {
            count: vec![0.0; n],
            seen: vec![false; n],
            touched: Vec::new(),
            dom: Vec::new(),
        }
    }

    /// The dominant labels among the labeled nodes v reads (by total
    /// weight), sorted ascending, into `self.dom`; empty when v reads no
    /// labeled node.
    fn dominant(&mut self, v: usize, reads: &[Vec<(usize, f64)>], label: &[usize]) {
        self.dom.clear();
        for &(u, x) in &reads[v] {
            let l = label[u];
            if l == UNLABELED {
                continue;
            }
            if !self.seen[l] {
                self.seen[l] = true;
                self.touched.push(l);
            }
            self.count[l] += x;
        }
        let mut best = f64::NEG_INFINITY;
        for &l in &self.touched {
            best = best.max(self.count[l]);
        }
        for &l in &self.touched {
            if self.count[l] == best {
                self.dom.push(l);
            }
            self.count[l] = 0.0;
            self.seen[l] = false;
        }
        self.dom.sort_unstable();
        self.touched.clear();
    }

    /// Whether v's label is dominant (or v reads no labeled node).
    fn satisfied(&mut self, v: usize, reads: &[Vec<(usize, f64)>], label: &[usize]) -> bool {
        self.dominant(v, reads, label);
        self.dom.is_empty() || self.dom.binary_search(&label[v]).is_ok()
    }
}

/// Label propagation communities. See the module documentation.
pub fn label_propagation(g: &Graph, opts: &LabelPropagationOptions) -> Result<Vec<usize>> {
    let n = g.node_count();
    if opts.variant == LpaVariant::SemiSynchronous && g.is_directed() {
        return Err(Error::InvalidArgument(
            "semi-synchronous label propagation is defined for undirected networks only".into(),
        ));
    }
    let mut label: Vec<usize> = match opts.initial {
        None => (0..n).collect(),
        Some(init) => {
            if init.len() != n {
                return Err(Error::LengthMismatch {
                    what: "initial labels",
                    expected: n,
                    found: init.len(),
                });
            }
            let labeled: Vec<i64> = init.iter().copied().filter(|&l| l >= 0).collect();
            let (dense, _) = crate::metrics::community::dense_labels(&labeled);
            let mut it = dense.into_iter();
            init.iter()
                .map(|&l| {
                    if l >= 0 {
                        it.next().unwrap()
                    } else {
                        UNLABELED
                    }
                })
                .collect()
        }
    };
    let fixed: Vec<bool> = match opts.fixed {
        None => vec![false; n],
        Some(f) => {
            if opts.initial.is_none() && f.iter().any(|&x| x) {
                return Err(Error::InvalidArgument(
                    "fixed labels need initial labels".into(),
                ));
            }
            if f.len() != n {
                return Err(Error::LengthMismatch {
                    what: "fixed",
                    expected: n,
                    found: f.len(),
                });
            }
            if let Some(v) = (0..n).find(|&v| f[v] && label[v] == UNLABELED) {
                return Err(Error::InvalidArgument(format!(
                    "node {v} is fixed but has no initial label"
                )));
            }
            f.to_vec()
        }
    };
    let lk = links(g, opts.mode, opts.weighted)?;
    let mut ct = Counter::new(n);
    let mut order: Vec<usize> = (0..n).collect();
    match opts.variant {
        LpaVariant::Dominance | LpaVariant::Retention => {
            let retain = opts.variant == LpaVariant::Retention;
            for sweep in 0..MAX_SWEEPS {
                let mut rng = Stream::new(opts.seed, sweep as u64);
                rng.permutation(&mut order);
                let mut changed = false;
                for &v in &order {
                    if fixed[v] {
                        continue;
                    }
                    ct.dominant(v, &lk.reads, &label);
                    let dom = &ct.dom;
                    if dom.is_empty() || (retain && dom.binary_search(&label[v]).is_ok()) {
                        continue;
                    }
                    let l = dom[rng.below(dom.len() as u64) as usize];
                    changed |= l != label[v];
                    label[v] = l;
                }
                let done = if retain {
                    !changed
                } else {
                    (0..n).all(|v| fixed[v] || ct.satisfied(v, &lk.reads, &label))
                };
                if done {
                    break;
                }
            }
        }
        LpaVariant::Fast => {
            let mut rng = Stream::new(opts.seed, 0);
            rng.permutation(&mut order);
            let mut queued = vec![false; n];
            let mut queue: VecDeque<usize> = VecDeque::new();
            for &v in &order {
                if !fixed[v] {
                    queue.push_back(v);
                    queued[v] = true;
                }
            }
            let mut steps = 0usize;
            let limit = MAX_SWEEPS.saturating_mul(n.max(1));
            while let Some(v) = queue.pop_front() {
                queued[v] = false;
                steps += 1;
                if steps > limit {
                    break;
                }
                ct.dominant(v, &lk.reads, &label);
                if ct.dom.is_empty() {
                    continue;
                }
                let l = ct.dom[rng.below(ct.dom.len() as u64) as usize];
                if l != label[v] {
                    label[v] = l;
                    for &u in &lk.readers[v] {
                        if !queued[u] && !fixed[u] && label[u] != l {
                            queued[u] = true;
                            queue.push_back(u);
                        }
                    }
                }
            }
        }
        LpaVariant::SemiSynchronous => {
            let classes = greedy_colouring(n, &lk.reads);
            loop {
                let stable = (0..n).all(|v| fixed[v] || ct.satisfied(v, &lk.reads, &label));
                if stable {
                    break;
                }
                for class in &classes {
                    // The nodes of a class are not adjacent, so updating
                    // them one at a time equals updating them together.
                    for &v in class {
                        if fixed[v] {
                            continue;
                        }
                        ct.dominant(v, &lk.reads, &label);
                        let dom = &ct.dom;
                        if dom.len() == 1 {
                            label[v] = dom[0];
                        } else if !dom.is_empty() && dom.binary_search(&label[v]).is_err() {
                            label[v] = *dom.last().unwrap();
                        }
                    }
                }
            }
        }
    }
    // Unlabeled nodes form groups of their own.
    let mut next = n;
    let keyed: Vec<usize> = label
        .iter()
        .map(|&l| {
            if l == UNLABELED {
                next += 1;
                next
            } else {
                l
            }
        })
        .collect();
    Ok(relabel(&keyed))
}

/// Greedy colouring, largest degree first (ties by node index), each node
/// taking the smallest colour none of its neighbours has. Returns the
/// colour classes in colour order.
fn greedy_colouring(n: usize, reads: &[Vec<(usize, f64)>]) -> Vec<Vec<usize>> {
    let simple: Vec<Vec<usize>> = reads
        .iter()
        .map(|r| {
            let mut s: Vec<usize> = r.iter().map(|&(u, _)| u).collect();
            s.sort_unstable();
            s.dedup();
            s
        })
        .collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| simple[b].len().cmp(&simple[a].len()).then(a.cmp(&b)));
    let mut colour = vec![usize::MAX; n];
    let mut classes: Vec<Vec<usize>> = Vec::new();
    let mut used: Vec<bool> = Vec::new();
    for &v in &order {
        used.clear();
        used.resize(simple[v].len() + 1, false);
        for &u in &simple[v] {
            if colour[u] < used.len() {
                used[colour[u]] = true;
            }
        }
        let c = used.iter().position(|&x| !x).unwrap();
        colour[v] = c;
        if c == classes.len() {
            classes.push(Vec::new());
        }
        classes[c].push(v);
    }
    classes
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    fn opts(variant: LpaVariant, seed: u64) -> LabelPropagationOptions<'static> {
        LabelPropagationOptions {
            variant,
            mode: LpaMode::All,
            weighted: false,
            initial: None,
            fixed: None,
            seed,
        }
    }

    fn is_fixed_point(g: &Graph, memb: &[usize]) -> bool {
        let lk = links(g, LpaMode::All, false).unwrap();
        let mut ct = Counter::new(g.node_count());
        (0..g.node_count()).all(|v| ct.satisfied(v, &lk.reads, memb))
    }

    #[test]
    fn two_cliques_and_fixed_points() {
        let mut e = Vec::new();
        for c in 0..2i64 {
            for i in 0..5 {
                for j in (i + 1)..5 {
                    e.push((5 * c + i, 5 * c + j));
                }
            }
        }
        e.push((4, 5));
        let g = und(10, &e);
        for variant in [
            LpaVariant::Dominance,
            LpaVariant::Retention,
            LpaVariant::Fast,
            LpaVariant::SemiSynchronous,
        ] {
            for seed in 0..10 {
                let m = label_propagation(&g, &opts(variant, seed)).unwrap();
                assert!(is_fixed_point(&g, &m), "{variant:?} {seed}");
                assert_eq!(m, label_propagation(&g, &opts(variant, seed)).unwrap());
            }
        }
    }

    #[test]
    fn semi_synchronous_by_hand() {
        // Path 0-1-2: colouring 1 -> 0, then 0 and 2 -> 1. Labels (0,1,2):
        // node 1 reads {0, 2}, a tie, takes the largest (2); then 0 and 2
        // read 2. Stable: everyone has label 2.
        let g = und(3, &[(0, 1), (1, 2)]);
        let m = label_propagation(&g, &opts(LpaVariant::SemiSynchronous, 0)).unwrap();
        assert_eq!(m, vec![0, 0, 0]);
    }

    #[test]
    fn initial_fixed_and_directed() {
        // Path 0-1-2-3-4 with ends fixed to different labels.
        let g = und(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        let init = [0, -1, -1, -1, 1];
        let fixed = [true, false, false, false, true];
        let mut o = opts(LpaVariant::Retention, 3);
        o.initial = Some(&init);
        o.fixed = Some(&fixed);
        let m = label_propagation(&g, &o).unwrap();
        assert_eq!((m[0], m[4]), (0, 1));
        assert!(m[1] == m[0] && m[3] == m[4]);
        // Directed out-star from 0: labels flow out from the hub.
        let d = Graph::from_edges(4, &[0, 0, 0], &[1, 2, 3], None, true).unwrap();
        let mut o = opts(LpaVariant::Dominance, 1);
        o.mode = LpaMode::Out;
        assert_eq!(label_propagation(&d, &o).unwrap(), vec![0, 0, 0, 0]);
        o.mode = LpaMode::In;
        // Against the edges the leaves read nothing; the hub draws a leaf.
        let m = label_propagation(&d, &o).unwrap();
        assert!(m[1] != m[2] && m[2] != m[3] && m[1] != m[3]);
        assert!(m[0] == m[1] || m[0] == m[2] || m[0] == m[3]);
        // An unlabeled node that no label reaches is a group of its own.
        let g = und(3, &[(0, 1)]);
        let init = [0, -1, -1];
        let mut o = opts(LpaVariant::Fast, 0);
        o.initial = Some(&init);
        assert_eq!(label_propagation(&g, &o).unwrap(), vec![0, 0, 1]);
        let bad = [false, true, false];
        o.fixed = Some(&bad);
        assert!(label_propagation(&g, &o).is_err());
        let mut s = opts(LpaVariant::SemiSynchronous, 0);
        s.mode = LpaMode::All;
        assert!(label_propagation(&d, &s).is_err());
    }
}
