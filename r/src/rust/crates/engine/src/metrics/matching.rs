//! Matchings: maximum-weight matching (Edmonds' blossom algorithm),
//! minimum-weight maximum-cardinality matching, greedy maximal matching,
//! and the matching tests.
//!
//! [`max_weight_matching`] finds a matching of largest total weight, or,
//! with `max_cardinality`, of largest total weight among the matchings of
//! largest size, by Edmonds' primal-dual blossom algorithm (Edmonds 1965,
//! `edmonds1965paths` and `edmonds1965maximum`) in the O(n³)
//! organization described by Galil (1986, `galil1986efficient`): stages
//! of alternating-tree growth from the free vertices, shrinking odd cycles
//! into blossoms, expanding blossoms whose dual reaches zero, and dual
//! adjustments of the four kinds; every stage either augments the matching
//! or proves it optimal. NetworkX 3.7's `max_weight_matching`
//! (`hagberg2008exploring`) implements the same algorithm.
//! [`min_weight_matching`] is NetworkX's `min_weight_matching`: among the
//! matchings of largest size, one of least total weight (NetworkX reaches
//! it by maximizing 1 + max w − w with `maxcardinality`).
//!
//! [`maximal_matching`] adds edges greedily in the graph's edge order,
//! skipping an edge if either end is already matched (and skipping
//! self-loops), as NetworkX's `maximal_matching` does over `G.edges()`;
//! a binding that passes edges in that order gets NetworkX's matching.
//!
//! Conventions:
//! - Undirected graphs only (direction is ignored for directed ones, which
//!   NetworkX rejects).
//! - Self-loops never belong to a matching.
//! - Parallel edges: the heaviest (lightest, for the minimum) copy of an
//!   edge stands for it (NetworkX rejects multigraphs).
//! - Weights may be any finite numbers; without weights each edge weighs
//!   one. Edges of non-positive weight never improve a maximum-weight
//!   matching, so they appear only when `max_cardinality` needs them.
//! - Exactness: the matching is optimal; its total weight equals
//!   NetworkX's (and so does its size with `max_cardinality`). When
//!   several matchings are optimal the one returned may differ from
//!   NetworkX's. With integer weights all dual values are
//!   multiples of 1/2 and the arithmetic is exact.
//! - Results list matched edges as indices into the graph's edges,
//!   increasing.

use crate::error::{Error, Result};
use crate::graph::Graph;

const NONE: usize = usize::MAX;

/// Candidate pairs, their weights, and the graph edges they stand for.
type Candidates = (Vec<(usize, usize)>, Vec<f64>, Vec<usize>);

/// The candidate edges: distinct non-loop pairs, each with its best copy.
fn candidates(g: &Graph, weighted: bool, prefer_max: bool) -> Result<Candidates> {
    let w = if weighted {
        Some(g.edge_weights().ok_or(Error::MissingWeights)?)
    } else {
        None
    };
    let mut best: Vec<((usize, usize), f64, usize)> = Vec::new();
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        if s == t {
            continue;
        }
        let key = (s.min(t) as usize, s.max(t) as usize);
        best.push((key, w.map_or(1.0, |w| w[e]), e));
    }
    // Sort by pair, then best weight first (ties: first edge).
    best.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| {
                let o = a.1.total_cmp(&b.1);
                if prefer_max {
                    o.reverse()
                } else {
                    o
                }
            })
            .then(a.2.cmp(&b.2))
    });
    best.dedup_by(|b, a| a.0 == b.0);
    Ok((
        best.iter().map(|x| x.0).collect(),
        best.iter().map(|x| x.1).collect(),
        best.iter().map(|x| x.2).collect(),
    ))
}

/// Edmonds' blossom algorithm on `n` vertices and edges `edges` with
/// weights `wt`. Returns `mate[v]`: the vertex matched to v or `NONE`.
struct Blossom<'a> {
    n: usize,
    edges: &'a [(usize, usize)],
    wt: &'a [f64],
    endpoint: Vec<usize>,
    neighbend: Vec<Vec<usize>>,
    mate: Vec<usize>,
    label: Vec<u8>,
    labelend: Vec<usize>,
    inblossom: Vec<usize>,
    blossomparent: Vec<usize>,
    blossomchilds: Vec<Vec<usize>>,
    blossombase: Vec<usize>,
    blossomendps: Vec<Vec<usize>>,
    bestedge: Vec<usize>,
    blossombestedges: Vec<Option<Vec<usize>>>,
    unusedblossoms: Vec<usize>,
    dualvar: Vec<f64>,
    allowedge: Vec<bool>,
    queue: Vec<usize>,
}

impl<'a> Blossom<'a> {
    fn new(n: usize, edges: &'a [(usize, usize)], wt: &'a [f64]) -> Blossom<'a> {
        let m = edges.len();
        let maxweight = wt.iter().cloned().fold(0.0f64, f64::max);
        let mut endpoint = vec![0; 2 * m];
        let mut neighbend = vec![Vec::new(); n];
        for (k, &(i, j)) in edges.iter().enumerate() {
            endpoint[2 * k] = i;
            endpoint[2 * k + 1] = j;
            neighbend[i].push(2 * k + 1);
            neighbend[j].push(2 * k);
        }
        let mut dualvar = vec![maxweight; n];
        dualvar.extend(std::iter::repeat(0.0).take(n));
        Blossom {
            n,
            edges,
            wt,
            endpoint,
            neighbend,
            mate: vec![NONE; n],
            label: vec![0; 2 * n],
            labelend: vec![NONE; 2 * n],
            inblossom: (0..n).collect(),
            blossomparent: vec![NONE; 2 * n],
            blossomchilds: vec![Vec::new(); 2 * n],
            blossombase: (0..n).chain(std::iter::repeat(NONE).take(n)).collect(),
            blossomendps: vec![Vec::new(); 2 * n],
            bestedge: vec![NONE; 2 * n],
            blossombestedges: vec![None; 2 * n],
            unusedblossoms: (n..2 * n).collect(),
            dualvar,
            allowedge: vec![false; m],
            queue: Vec::new(),
        }
    }

    #[inline]
    fn slack(&self, k: usize) -> f64 {
        let (i, j) = self.edges[k];
        self.dualvar[i] + self.dualvar[j] - 2.0 * self.wt[k]
    }

    fn leaves(&self, b: usize, out: &mut Vec<usize>) {
        if b < self.n {
            out.push(b);
        } else {
            for &t in &self.blossomchilds[b] {
                self.leaves(t, out);
            }
        }
    }

    fn leaves_of(&self, b: usize) -> Vec<usize> {
        let mut v = Vec::new();
        self.leaves(b, &mut v);
        v
    }

    fn assign_label(&mut self, w: usize, t: u8, p: usize) {
        let b = self.inblossom[w];
        self.label[w] = t;
        self.label[b] = t;
        self.labelend[w] = p;
        self.labelend[b] = p;
        self.bestedge[w] = NONE;
        self.bestedge[b] = NONE;
        if t == 1 {
            let l = self.leaves_of(b);
            self.queue.extend(l);
        } else if t == 2 {
            let base = self.blossombase[b];
            let mb = self.mate[base];
            self.assign_label(self.endpoint[mb], 1, mb ^ 1);
        }
    }

    /// Traces back from v and w to find a new blossom's base, or `NONE` if
    /// they lie in different trees (an augmenting path).
    fn scan_blossom(&mut self, mut v: usize, mut w: usize) -> usize {
        let mut path = Vec::new();
        let mut base = NONE;
        while v != NONE || w != NONE {
            let mut b = self.inblossom[v];
            if self.label[b] & 4 != 0 {
                base = self.blossombase[b];
                break;
            }
            path.push(b);
            self.label[b] = 5;
            if self.labelend[b] == NONE {
                v = NONE;
            } else {
                v = self.endpoint[self.labelend[b]];
                b = self.inblossom[v];
                v = self.endpoint[self.labelend[b]];
            }
            if w != NONE {
                std::mem::swap(&mut v, &mut w);
            }
        }
        for b in path {
            self.label[b] = 1;
        }
        base
    }

    fn add_blossom(&mut self, base: usize, k: usize) {
        let (mut v, mut w) = self.edges[k];
        let bb = self.inblossom[base];
        let mut bv = self.inblossom[v];
        let mut bw = self.inblossom[w];
        let b = self.unusedblossoms.pop().expect("a free blossom slot");
        self.blossombase[b] = base;
        self.blossomparent[b] = NONE;
        self.blossomparent[bb] = b;
        let mut path = Vec::new();
        let mut endps = Vec::new();
        while bv != bb {
            self.blossomparent[bv] = b;
            path.push(bv);
            endps.push(self.labelend[bv]);
            v = self.endpoint[self.labelend[bv]];
            bv = self.inblossom[v];
        }
        path.push(bb);
        path.reverse();
        endps.reverse();
        endps.push(2 * k);
        while bw != bb {
            self.blossomparent[bw] = b;
            path.push(bw);
            endps.push(self.labelend[bw] ^ 1);
            w = self.endpoint[self.labelend[bw]];
            bw = self.inblossom[w];
        }
        self.label[b] = 1;
        self.labelend[b] = self.labelend[bb];
        self.dualvar[b] = 0.0;
        self.blossomchilds[b] = path.clone();
        self.blossomendps[b] = endps;
        for v in self.leaves_of(b) {
            if self.label[self.inblossom[v]] == 2 {
                self.queue.push(v);
            }
            self.inblossom[v] = b;
        }
        // The least-slack edges from the new blossom to each other
        // S-blossom.
        let mut bestedgeto = vec![NONE; 2 * self.n];
        for &bv in &path {
            let lists: Vec<Vec<usize>> = match self.blossombestedges[bv].take() {
                Some(l) => vec![l],
                None => self
                    .leaves_of(bv)
                    .into_iter()
                    .map(|v| self.neighbend[v].iter().map(|&p| p / 2).collect())
                    .collect(),
            };
            for list in lists {
                for k in list {
                    let (i, j) = self.edges[k];
                    // The end outside the new blossom (if any).
                    let j = if self.inblossom[j] == b { i } else { j };
                    let bj = self.inblossom[j];
                    if bj != b
                        && self.label[bj] == 1
                        && (bestedgeto[bj] == NONE || self.slack(k) < self.slack(bestedgeto[bj]))
                    {
                        bestedgeto[bj] = k;
                    }
                }
            }
            self.bestedge[bv] = NONE;
        }
        let best: Vec<usize> = bestedgeto.into_iter().filter(|&k| k != NONE).collect();
        self.bestedge[b] = NONE;
        for &k in &best {
            if self.bestedge[b] == NONE || self.slack(k) < self.slack(self.bestedge[b]) {
                self.bestedge[b] = k;
            }
        }
        self.blossombestedges[b] = Some(best);
    }

    #[inline]
    fn at(list: &[usize], j: isize) -> usize {
        let len = list.len() as isize;
        list[(((j % len) + len) % len) as usize]
    }

    fn expand_blossom(&mut self, b: usize, endstage: bool) {
        let childs = self.blossomchilds[b].clone();
        for &s in &childs {
            self.blossomparent[s] = NONE;
            if s < self.n {
                self.inblossom[s] = s;
            } else if endstage && self.dualvar[s] == 0.0 {
                self.expand_blossom(s, endstage);
            } else {
                for v in self.leaves_of(s) {
                    self.inblossom[v] = s;
                }
            }
        }
        if !endstage && self.label[b] == 2 {
            // Relabel the children along the even-length path from the
            // child that holds the T-label's entry to the base.
            let endps = self.blossomendps[b].clone();
            let entrychild = self.inblossom[self.endpoint[self.labelend[b] ^ 1]];
            let mut j = childs.iter().position(|&c| c == entrychild).expect("child") as isize;
            let (jstep, endptrick): (isize, usize) = if j & 1 != 0 {
                j -= childs.len() as isize;
                (1, 0)
            } else {
                (-1, 1)
            };
            let mut p = self.labelend[b];
            while j != 0 {
                let a = self.endpoint[p ^ 1];
                self.label[a] = 0;
                let q = Self::at(&endps, j - endptrick as isize) ^ endptrick ^ 1;
                let c = self.endpoint[q];
                self.label[c] = 0;
                self.assign_label(self.endpoint[p ^ 1], 2, p);
                self.allowedge[Self::at(&endps, j - endptrick as isize) / 2] = true;
                j += jstep;
                p = Self::at(&endps, j - endptrick as isize) ^ endptrick;
                self.allowedge[p / 2] = true;
                j += jstep;
            }
            let bv = Self::at(&childs, j);
            let a = self.endpoint[p ^ 1];
            self.label[a] = 2;
            self.label[bv] = 2;
            self.labelend[a] = p;
            self.labelend[bv] = p;
            self.bestedge[bv] = NONE;
            j += jstep;
            while Self::at(&childs, j) != entrychild {
                let bv = Self::at(&childs, j);
                if self.label[bv] == 1 {
                    j += jstep;
                    continue;
                }
                let found = self.leaves_of(bv).into_iter().find(|&v| self.label[v] != 0);
                if let Some(v) = found {
                    self.label[v] = 0;
                    let mb = self.mate[self.blossombase[bv]];
                    let e = self.endpoint[mb];
                    self.label[e] = 0;
                    self.assign_label(v, 2, self.labelend[v]);
                }
                j += jstep;
            }
        }
        self.label[b] = 0;
        self.labelend[b] = NONE;
        self.blossomchilds[b].clear();
        self.blossomendps[b].clear();
        self.blossombase[b] = NONE;
        self.blossombestedges[b] = None;
        self.bestedge[b] = NONE;
        self.unusedblossoms.push(b);
    }

    fn augment_blossom(&mut self, b: usize, v: usize) {
        let mut t = v;
        while self.blossomparent[t] != b {
            t = self.blossomparent[t];
        }
        if t >= self.n {
            self.augment_blossom(t, v);
        }
        let childs = self.blossomchilds[b].clone();
        let endps = self.blossomendps[b].clone();
        let i = childs.iter().position(|&c| c == t).expect("child");
        let mut j = i as isize;
        let (jstep, endptrick): (isize, usize) = if i & 1 != 0 {
            j -= childs.len() as isize;
            (1, 0)
        } else {
            (-1, 1)
        };
        while j != 0 {
            j += jstep;
            let mut t = Self::at(&childs, j);
            let p = Self::at(&endps, j - endptrick as isize) ^ endptrick;
            if t >= self.n {
                self.augment_blossom(t, self.endpoint[p]);
            }
            j += jstep;
            t = Self::at(&childs, j);
            if t >= self.n {
                self.augment_blossom(t, self.endpoint[p ^ 1]);
            }
            let (a, c) = (self.endpoint[p], self.endpoint[p ^ 1]);
            self.mate[a] = p ^ 1;
            self.mate[c] = p;
        }
        let mut rc = childs[i..].to_vec();
        rc.extend_from_slice(&childs[..i]);
        let mut re = endps[i..].to_vec();
        re.extend_from_slice(&endps[..i]);
        self.blossomchilds[b] = rc;
        self.blossomendps[b] = re;
        self.blossombase[b] = self.blossombase[self.blossomchilds[b][0]];
    }

    fn augment_matching(&mut self, k: usize) {
        let (v, w) = self.edges[k];
        for (s0, p0) in [(v, 2 * k + 1), (w, 2 * k)] {
            let (mut s, mut p) = (s0, p0);
            loop {
                let bs = self.inblossom[s];
                if bs >= self.n {
                    self.augment_blossom(bs, s);
                }
                self.mate[s] = p;
                if self.labelend[bs] == NONE {
                    break;
                }
                let t = self.endpoint[self.labelend[bs]];
                let bt = self.inblossom[t];
                s = self.endpoint[self.labelend[bt]];
                let j = self.endpoint[self.labelend[bt] ^ 1];
                if bt >= self.n {
                    self.augment_blossom(bt, j);
                }
                self.mate[j] = self.labelend[bt];
                p = self.labelend[bt] ^ 1;
            }
        }
    }

    fn solve(&mut self, max_cardinality: bool) {
        let n = self.n;
        for _stage in 0..n {
            self.label.iter_mut().for_each(|x| *x = 0);
            self.bestedge.iter_mut().for_each(|x| *x = NONE);
            for b in n..2 * n {
                self.blossombestedges[b] = None;
            }
            self.allowedge.iter_mut().for_each(|x| *x = false);
            self.queue.clear();
            for v in 0..n {
                if self.mate[v] == NONE && self.label[self.inblossom[v]] == 0 {
                    self.assign_label(v, 1, NONE);
                }
            }
            let mut augmented = false;
            loop {
                while let Some(v) = self.queue.pop() {
                    let ends = self.neighbend[v].clone();
                    for p in ends {
                        let k = p / 2;
                        let w = self.endpoint[p];
                        if self.inblossom[v] == self.inblossom[w] {
                            continue;
                        }
                        let mut kslack = 0.0;
                        if !self.allowedge[k] {
                            kslack = self.slack(k);
                            if kslack <= 0.0 {
                                self.allowedge[k] = true;
                            }
                        }
                        if self.allowedge[k] {
                            let bw = self.inblossom[w];
                            if self.label[bw] == 0 {
                                self.assign_label(w, 2, p ^ 1);
                            } else if self.label[bw] == 1 {
                                let base = self.scan_blossom(v, w);
                                if base != NONE {
                                    self.add_blossom(base, k);
                                } else {
                                    self.augment_matching(k);
                                    augmented = true;
                                    break;
                                }
                            } else if self.label[w] == 0 {
                                self.label[w] = 2;
                                self.labelend[w] = p ^ 1;
                            }
                        } else if self.label[self.inblossom[w]] == 1 {
                            let b = self.inblossom[v];
                            if self.bestedge[b] == NONE || kslack < self.slack(self.bestedge[b]) {
                                self.bestedge[b] = k;
                            }
                        } else if self.label[w] == 0
                            && (self.bestedge[w] == NONE || kslack < self.slack(self.bestedge[w]))
                        {
                            self.bestedge[w] = k;
                        }
                    }
                    if augmented {
                        break;
                    }
                }
                if augmented {
                    break;
                }
                // Dual adjustment.
                let mut deltatype = 0u8;
                let mut delta = 0.0;
                let mut deltaedge = NONE;
                let mut deltablossom = NONE;
                if !max_cardinality {
                    deltatype = 1;
                    delta = self.dualvar[..n]
                        .iter()
                        .cloned()
                        .fold(f64::INFINITY, f64::min);
                }
                for v in 0..n {
                    if self.label[self.inblossom[v]] == 0 && self.bestedge[v] != NONE {
                        let d = self.slack(self.bestedge[v]);
                        if deltatype == 0 || d < delta {
                            delta = d;
                            deltatype = 2;
                            deltaedge = self.bestedge[v];
                        }
                    }
                }
                for b in 0..2 * n {
                    if self.blossomparent[b] == NONE
                        && self.label[b] == 1
                        && self.bestedge[b] != NONE
                    {
                        let d = self.slack(self.bestedge[b]) / 2.0;
                        if deltatype == 0 || d < delta {
                            delta = d;
                            deltatype = 3;
                            deltaedge = self.bestedge[b];
                        }
                    }
                }
                for b in n..2 * n {
                    if self.blossombase[b] != NONE
                        && self.blossomparent[b] == NONE
                        && self.label[b] == 2
                        && (deltatype == 0 || self.dualvar[b] < delta)
                    {
                        delta = self.dualvar[b];
                        deltatype = 4;
                        deltablossom = b;
                    }
                }
                if deltatype == 0 {
                    // No further progress with max_cardinality: the
                    // matching is of largest size; finish with a final
                    // vertex dual adjustment.
                    deltatype = 1;
                    delta = self.dualvar[..n]
                        .iter()
                        .cloned()
                        .fold(f64::INFINITY, f64::min)
                        .max(0.0);
                }
                for v in 0..n {
                    match self.label[self.inblossom[v]] {
                        1 => self.dualvar[v] -= delta,
                        2 => self.dualvar[v] += delta,
                        _ => {}
                    }
                }
                for b in n..2 * n {
                    if self.blossombase[b] != NONE && self.blossomparent[b] == NONE {
                        match self.label[b] {
                            1 => self.dualvar[b] += delta,
                            2 => self.dualvar[b] -= delta,
                            _ => {}
                        }
                    }
                }
                match deltatype {
                    1 => break,
                    2 | 3 => {
                        self.allowedge[deltaedge] = true;
                        let (mut i, j) = self.edges[deltaedge];
                        if deltatype == 2 && self.label[self.inblossom[i]] == 0 {
                            i = j;
                        }
                        self.queue.push(i);
                    }
                    _ => self.expand_blossom(deltablossom, false),
                }
            }
            if !augmented {
                break;
            }
            for b in n..2 * n {
                if self.blossomparent[b] == NONE
                    && self.blossombase[b] != NONE
                    && self.label[b] == 1
                    && self.dualvar[b] == 0.0
                {
                    self.expand_blossom(b, true);
                }
            }
        }
    }
}

fn run(n: usize, edges: &[(usize, usize)], wt: &[f64], max_cardinality: bool) -> Vec<usize> {
    if n == 0 || edges.is_empty() {
        return vec![NONE; n];
    }
    let mut b = Blossom::new(n, edges, wt);
    b.solve(max_cardinality);
    (0..n)
        .map(|v| {
            let p = b.mate[v];
            if p == NONE {
                NONE
            } else {
                b.endpoint[p]
            }
        })
        .collect()
}

/// Matched edges (graph edge indices, increasing) from the mate array.
fn matched(mate: &[usize], pairs: &[(usize, usize)], ids: &[usize]) -> Vec<usize> {
    let mut out: Vec<usize> = pairs
        .iter()
        .zip(ids)
        .filter(|((i, j), _)| mate[*i] == *j)
        .map(|(_, &e)| e)
        .collect();
    out.sort_unstable();
    out
}

/// A maximum-weight matching (with `max_cardinality`: of largest weight
/// among the largest matchings). See the module documentation.
pub fn max_weight_matching(g: &Graph, weighted: bool, max_cardinality: bool) -> Result<Vec<usize>> {
    let (pairs, wt, ids) = candidates(g, weighted, true)?;
    if wt.iter().any(|w| !w.is_finite()) {
        return Err(Error::InvalidArgument("weights must be finite".into()));
    }
    let mate = run(g.node_count(), &pairs, &wt, max_cardinality);
    Ok(matched(&mate, &pairs, &ids))
}

/// A matching of least total weight among the matchings of largest size
/// (NetworkX's `min_weight_matching`). See the module documentation.
pub fn min_weight_matching(g: &Graph, weighted: bool) -> Result<Vec<usize>> {
    let (pairs, wt, ids) = candidates(g, weighted, false)?;
    let top = wt.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let inv: Vec<f64> = wt.iter().map(|&w| 1.0 + top - w).collect();
    let mate = run(g.node_count(), &pairs, &inv, true);
    Ok(matched(&mate, &pairs, &ids))
}

/// A maximal matching by greedy scan of the edges in graph order. See the
/// module documentation.
pub fn maximal_matching(g: &Graph) -> Vec<usize> {
    let mut used = vec![false; g.node_count()];
    let mut out = Vec::new();
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        let (s, t) = (s as usize, t as usize);
        if s != t && !used[s] && !used[t] {
            used[s] = true;
            used[t] = true;
            out.push(e);
        }
    }
    out
}

/// What a set of node pairs is with respect to a graph, as NetworkX's
/// `is_matching`, `is_maximal_matching`, and `is_perfect_matching` test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchingCheck {
    /// Every pair is an edge of the graph (not a loop) and no node is
    /// repeated.
    pub matching: bool,
    /// A matching that no edge of the graph can extend.
    pub maximal: bool,
    /// A matching that covers every node.
    pub perfect: bool,
}

/// Tests the pairs `(us[i], vs[i])` (a set of edges: a pair given twice,
/// in either orientation, repeats its nodes and is not a matching, as in
/// NetworkX; in a directed graph a pair must be an arc u -> v). A node
/// outside the graph is an error, as in NetworkX.
pub fn check_matching(g: &Graph, us: &[u32], vs: &[u32]) -> Result<MatchingCheck> {
    let n = g.node_count();
    if us.len() != vs.len() {
        return Err(Error::LengthMismatch {
            what: "second nodes of the pairs",
            expected: us.len(),
            found: vs.len(),
        });
    }
    // An arc a -> b in a directed graph, as NetworkX's has_edge.
    let adjacent = |a: usize, b: usize| g.out().neighbors(a).binary_search(&(b as u32)).is_ok();
    let mut used = vec![false; n];
    let mut ok = true;
    for (&a, &b) in us.iter().zip(vs) {
        let (a, b) = (a as usize, b as usize);
        if a >= n || b >= n {
            return Err(Error::InvalidArgument(format!(
                "pair ({a}, {b}) names a node outside 0..{n}"
            )));
        }
        if a == b || !adjacent(a, b) || used[a] || used[b] {
            ok = false;
        }
        used[a] = true;
        used[b] = true;
    }
    let maximal = ok
        && !g
            .edge_sources()
            .iter()
            .zip(g.edge_targets())
            .any(|(&s, &t)| s != t && !used[s as usize] && !used[t as usize]);
    let perfect = ok && used.iter().all(|&u| u);
    Ok(MatchingCheck {
        matching: ok,
        maximal,
        perfect,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wgraph(n: usize, e: &[(i64, i64, f64)]) -> Graph {
        let s: Vec<i64> = e.iter().map(|x| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        let w: Vec<f64> = e.iter().map(|x| x.2).collect();
        Graph::from_edges(n, &s, &t, Some(&w), false).unwrap()
    }

    fn weight(g: &Graph, m: &[usize]) -> f64 {
        m.iter().map(|&e| g.edge_weights().unwrap()[e]).sum()
    }

    #[test]
    fn classic_cases() {
        // Path 0-1-2-3 with a heavy middle: weight 3 beats 1 + 1.
        let g = wgraph(4, &[(0, 1, 1.0), (1, 2, 3.0), (2, 3, 1.0)]);
        assert_eq!(max_weight_matching(&g, true, false).unwrap(), vec![1]);
        // With max cardinality, both outer edges.
        assert_eq!(max_weight_matching(&g, true, true).unwrap(), vec![0, 2]);
        // A triangle 0-1-2 with pendants 2-3 and 0-4: 11 at best.
        let g = wgraph(
            5,
            &[
                (0, 1, 6.0),
                (1, 2, 6.0),
                (2, 0, 6.0),
                (2, 3, 5.0),
                (0, 4, 5.0),
            ],
        );
        let m = max_weight_matching(&g, true, false).unwrap();
        assert_eq!(weight(&g, &m), 11.0);
        // An S-blossom used for augmentation: 0-1 and 2-3 (15).
        let g = wgraph(4, &[(0, 1, 8.0), (0, 2, 9.0), (1, 2, 10.0), (2, 3, 7.0)]);
        assert_eq!(max_weight_matching(&g, true, false).unwrap(), vec![0, 3]);
        // An S-blossom relabelled as a T-blossom: 0-5, 1-2, 3-4 (17).
        let g = wgraph(
            6,
            &[
                (0, 1, 9.0),
                (0, 2, 8.0),
                (1, 2, 10.0),
                (0, 3, 5.0),
                (3, 4, 4.0),
                (0, 5, 3.0),
            ],
        );
        let m = max_weight_matching(&g, true, false).unwrap();
        assert_eq!(m, vec![2, 4, 5]);
        // Minimum weight among perfect matchings of a 4-cycle.
        let g = wgraph(4, &[(0, 1, 1.0), (1, 2, 5.0), (2, 3, 1.0), (3, 0, 5.0)]);
        assert_eq!(min_weight_matching(&g, true).unwrap(), vec![0, 2]);
        let c = check_matching(&g, &[0, 2], &[1, 3]).unwrap();
        assert!(c.matching && c.maximal && c.perfect);
        assert!(!check_matching(&g, &[0], &[2]).unwrap().matching);
        assert_eq!(maximal_matching(&g), vec![0, 2]);
    }
}
