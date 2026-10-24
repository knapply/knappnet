//! Graph isomorphism, subgraph isomorphism, and monomorphism (VF2++),
//! quick necessary conditions for isomorphism, and tree isomorphism.
//!
//! # VF2++
//!
//! [`Matcher`] enumerates the mappings of a smaller graph SG into a full
//! graph FG by the VF2++ algorithm (Jüttner and Madarasi 2018,
//! `juttner2018vf2`), the refinement of VF2 (Cordella, Foggia, Sansone,
//! and Vento 2004, `cordella2004subgraph`) that NetworkX 3.7's `vf2pp_*`
//! functions implement (`hagberg2008exploring`):
//!
//! - [`Problem::Isomorphism`]: a bijection between the nodes that
//!   preserves adjacency and its absence;
//! - [`Problem::InducedSubgraph`]: an injection of SG into FG such that SG
//!   is isomorphic to the subgraph of FG induced by the image (NetworkX's
//!   `vf2pp_subgraph_*`);
//! - [`Problem::Monomorphism`]: an injection that maps every edge of SG
//!   onto an edge of FG (SG isomorphic to a subgraph of FG, not
//!   necessarily induced; `vf2pp_*monomorphism`).
//!
//! Nodes are visited in VF2++'s matching order (breadth-first from the
//! node with the rarest label and largest degree, preferring nodes with
//! many neighbours already ordered); a candidate pair is accepted after
//! the consistency rules (edges to mapped nodes) and VF2++'s cutting rules
//! (per label, the numbers of neighbours in the frontier and outside it).
//! Every rule is a necessary condition, so the mappings produced are
//! exactly the mappings of the problem, each once; only the order in which
//! they come differs from NetworkX's (which depends on Python set order).
//!
//! Conventions, as in NetworkX's `vf2pp_*`:
//! - Node labels, when given, must be equal for mapped nodes.
//! - Parallel edges: edge multiplicities must be equal (isomorphism and
//!   induced subgraph isomorphism) or no larger in SG (monomorphism), as
//!   for NetworkX multigraphs; simple graphs have multiplicity 1.
//! - Self-loops are counted like other edges: a node maps only onto a node
//!   with the same number of loops (at least as many for monomorphism).
//! - Direction: both graphs directed (arcs must match in direction) or
//!   both undirected.
//! - Weights are ignored.
//! - A mapping is returned as, for every node of SG, the node of FG it
//!   maps to (NetworkX returns the inverse, a dict from FG nodes to SG
//!   nodes).
//! - Graphs without nodes: NetworkX reports no mapping when either graph
//!   has no nodes (so two empty graphs are "not isomorphic");
//!   `empty_matches` gives the mathematical answer instead (one empty
//!   mapping when SG has no nodes).
//!
//! # Quick checks
//!
//! [`could_be_isomorphic`] compares invariants that isomorphic graphs
//! share, as NetworkX's `could_be_isomorphic`: the sorted sequence of the
//! node degrees, then of (degree, triangles), then of (degree, triangles,
//! number of maximal cliques containing the node). Degrees count a
//! self-loop twice (undirected) or once in each direction, and parallel
//! edges each time; triangles and cliques are those of the simple
//! undirected graph and are defined for undirected graphs only (NetworkX
//! raises for directed ones).
//!
//! # Tree isomorphism
//!
//! [`rooted_tree_isomorphism`] decides whether two rooted trees are
//! isomorphic with the roots mapped onto each other by the level-by-level
//! canonical labelling of Aho, Hopcroft, and Ullman (1974,
//! `aho1974design`), in O(n log n), and returns an isomorphism.
//! [`tree_isomorphism`] does the same for free trees, rooted at their
//! centres, as NetworkX's `tree_isomorphism` does. Direction is ignored; a
//! graph that is not a tree (connected with n − 1 edges, parallel edges
//! and loops counted) is an error, as is a graph without nodes. The
//! isomorphism returned is valid but may differ from NetworkX's when the
//! trees have several.

use std::collections::VecDeque;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::cliques::max_cliques;
use crate::metrics::descriptive::triangles;

const NONE: u32 = u32::MAX;

/// A neighbour and the multiplicity of the edges to it.
type Arc = (u32, u32);

/// The mapping problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Problem {
    /// Graph isomorphism.
    Isomorphism,
    /// SG isomorphic to an induced subgraph of FG.
    InducedSubgraph,
    /// SG isomorphic to a subgraph (not necessarily induced) of FG.
    Monomorphism,
}

impl Problem {
    /// Decodes the C interface codes: 0 = isomorphism, 1 = induced
    /// subgraph isomorphism, 2 = monomorphism.
    pub fn from_code(code: i32) -> Result<Problem> {
        match code {
            0 => Ok(Problem::Isomorphism),
            1 => Ok(Problem::InducedSubgraph),
            2 => Ok(Problem::Monomorphism),
            _ => Err(Error::InvalidArgument(format!(
                "unknown mapping problem code {code}"
            ))),
        }
    }
}

/// One graph's adjacency with multiplicities.
#[derive(Debug, Clone)]
struct Side {
    n: usize,
    directed: bool,
    /// Distinct out-neighbours (all neighbours if undirected), sorted, with
    /// edge multiplicities; a loop appears in its own node's row.
    out: Vec<Vec<(u32, u32)>>,
    /// Distinct in-neighbours (directed graphs only).
    inc: Vec<Vec<(u32, u32)>>,
    loops: Vec<u32>,
    /// (in-degree, out-degree), or (degree, 0) when undirected, with loops
    /// counted as in NetworkX.
    deg: Vec<(u32, u32)>,
    label: Vec<u32>,
}

fn grouped(mut pairs: Vec<u32>) -> Vec<(u32, u32)> {
    pairs.sort_unstable();
    let mut out: Vec<(u32, u32)> = Vec::new();
    for v in pairs {
        match out.last_mut() {
            Some(last) if last.0 == v => last.1 += 1,
            _ => out.push((v, 1)),
        }
    }
    out
}

impl Side {
    fn new(g: &Graph, label: Vec<u32>) -> Side {
        let n = g.node_count();
        let directed = g.is_directed();
        let out: Vec<Vec<(u32, u32)>> = (0..n)
            .map(|v| grouped(g.out().neighbors(v).to_vec()))
            .collect();
        let inc: Vec<Vec<(u32, u32)>> = if directed {
            (0..n)
                .map(|v| grouped(g.inc().neighbors(v).to_vec()))
                .collect()
        } else {
            vec![Vec::new(); n]
        };
        let count = |row: &[(u32, u32)], v: usize| {
            row.binary_search_by_key(&(v as u32), |p| p.0)
                .map_or(0, |i| row[i].1)
        };
        let loops: Vec<u32> = (0..n).map(|v| count(&out[v], v)).collect();
        let deg = (0..n)
            .map(|v| {
                let o: u32 = out[v].iter().map(|p| p.1).sum();
                if directed {
                    let i: u32 = inc[v].iter().map(|p| p.1).sum();
                    (i, o)
                } else {
                    (o + loops[v], 0)
                }
            })
            .collect();
        Side {
            n,
            directed,
            out,
            inc,
            loops,
            deg,
            label,
        }
    }

    #[inline]
    fn count(&self, a: usize, b: usize) -> u32 {
        let row = &self.out[a];
        row.binary_search_by_key(&(b as u32), |p| p.0)
            .map_or(0, |i| row[i].1)
    }

    /// Distinct neighbours in either direction (undirected view).
    fn neighbours(&self, v: usize) -> Vec<u32> {
        let mut x: Vec<u32> = self.out[v].iter().map(|p| p.0).collect();
        x.extend(self.inc[v].iter().map(|p| p.0));
        x.sort_unstable();
        x.dedup();
        x
    }
}

/// Mapping state of one side: the partner of every node and how many
/// mapped nodes point to it (`t_out`: it is their successor) and from it
/// (`t_in`).
#[derive(Debug, Clone)]
struct State {
    core: Vec<u32>,
    t_out: Vec<u32>,
    t_in: Vec<u32>,
}

impl State {
    fn new(n: usize) -> State {
        State {
            core: vec![NONE; n],
            t_out: vec![0; n],
            t_in: vec![0; n],
        }
    }

    fn add(&mut self, s: &Side, v: usize, partner: u32) {
        self.core[v] = partner;
        for &(w, _) in &s.out[v] {
            self.t_out[w as usize] += 1;
        }
        for &(w, _) in &s.inc[v] {
            self.t_in[w as usize] += 1;
        }
    }

    fn remove(&mut self, s: &Side, v: usize) {
        self.core[v] = NONE;
        for &(w, _) in &s.out[v] {
            self.t_out[w as usize] -= 1;
        }
        for &(w, _) in &s.inc[v] {
            self.t_in[w as usize] -= 1;
        }
    }

    /// Frontier memberships of an unmapped node as bits: 1 = a successor
    /// of a mapped node, 2 = a predecessor of one (undirected graphs use
    /// bit 1 only), 4 = neither (T-tilde); 0 for a mapped node.
    #[inline]
    fn classes(&self, v: usize) -> u8 {
        if self.core[v] != NONE {
            return 0;
        }
        let b = (self.t_out[v] > 0) as u8 | ((self.t_in[v] > 0) as u8) << 1;
        if b == 0 {
            4
        } else {
            b
        }
    }
}

struct Frame {
    cands: Vec<u32>,
    pos: usize,
}

/// Enumerates the mappings of SG into FG (see the module documentation).
pub struct Matcher {
    fg: Side,
    sg: Side,
    problem: Problem,
    order: Vec<u32>,
    by_label: Vec<Vec<u32>>,
    sf: State,
    ss: State,
    stack: Vec<Frame>,
    started: bool,
    done: bool,
    empty_mapping: bool,
    labels: usize,
    scratch: Vec<[u32; 8]>,
    touched: Vec<u32>,
    steps: u64,
}

impl Matcher {
    /// A matcher for SG into FG. Labels, when given, hold one value per
    /// node of each graph. With `empty_matches`, an SG without nodes has
    /// one (empty) mapping; otherwise, as in NetworkX, a graph without
    /// nodes on either side gives none.
    pub fn new(
        fg: &Graph,
        sg: &Graph,
        fg_labels: Option<&[i64]>,
        sg_labels: Option<&[i64]>,
        problem: Problem,
        empty_matches: bool,
    ) -> Result<Matcher> {
        if fg.is_directed() != sg.is_directed() {
            return Err(Error::InvalidArgument(
                "both graphs must be directed or both undirected".into(),
            ));
        }
        let (nf, ns) = (fg.node_count(), sg.node_count());
        let check = |l: Option<&[i64]>, n: usize, what: &'static str| -> Result<()> {
            match l {
                Some(l) if l.len() != n => Err(Error::LengthMismatch {
                    what,
                    expected: n,
                    found: l.len(),
                }),
                _ => Ok(()),
            }
        };
        check(fg_labels, nf, "labels of the full graph")?;
        check(sg_labels, ns, "labels of the smaller graph")?;
        if fg_labels.is_some() != sg_labels.is_some() {
            return Err(Error::InvalidArgument(
                "labels must be given for both graphs or for neither".into(),
            ));
        }
        // Dense label ids shared by both graphs.
        let mut ids: Vec<i64> = fg_labels
            .unwrap_or(&[])
            .iter()
            .chain(sg_labels.unwrap_or(&[]))
            .copied()
            .collect();
        ids.sort_unstable();
        ids.dedup();
        let dense = |l: Option<&[i64]>, n: usize| -> Vec<u32> {
            match l {
                Some(l) => l
                    .iter()
                    .map(|x| ids.binary_search(x).unwrap() as u32)
                    .collect(),
                None => vec![0; n],
            }
        };
        let labels = ids.len().max(1);
        let fgs = Side::new(fg, dense(fg_labels, nf));
        let sgs = Side::new(sg, dense(sg_labels, ns));
        let mut by_label = vec![Vec::new(); labels];
        for v in 0..nf {
            by_label[fgs.label[v] as usize].push(v as u32);
        }
        let mut m = Matcher {
            sf: State::new(nf),
            ss: State::new(ns),
            fg: fgs,
            sg: sgs,
            problem,
            order: Vec::new(),
            by_label,
            stack: Vec::new(),
            started: false,
            done: false,
            empty_mapping: false,
            labels,
            scratch: vec![[0; 8]; labels],
            touched: Vec::new(),
            steps: 0,
        };
        if ns == 0 || nf == 0 {
            m.done = true;
            m.empty_mapping = empty_matches && ns == 0;
            return Ok(m);
        }
        if !m.sizes_fit() {
            m.done = true;
            return Ok(m);
        }
        m.order = m.matching_order();
        Ok(m)
    }

    /// Number of candidate pairs tested so far.
    pub fn steps(&self) -> u64 {
        self.steps
    }

    #[inline]
    fn fits(&self, s: u32, f: u32) -> bool {
        match self.problem {
            Problem::Isomorphism => s == f,
            _ => s <= f,
        }
    }

    #[inline]
    fn edges_fit(&self, s: u32, f: u32) -> bool {
        match self.problem {
            Problem::Monomorphism => s <= f,
            _ => s == f,
        }
    }

    /// Node and label counts that rule out every mapping.
    fn sizes_fit(&self) -> bool {
        let (ns, nf) = (self.sg.n as u32, self.fg.n as u32);
        if !self.fits(ns, nf) {
            return false;
        }
        let mut cs = vec![0u32; self.labels];
        let mut cf = vec![0u32; self.labels];
        for &l in &self.sg.label {
            cs[l as usize] += 1;
        }
        for &l in &self.fg.label {
            cf[l as usize] += 1;
        }
        (0..self.labels).all(|l| cs[l] == 0 || (cf[l] > 0 && self.fits(cs[l], cf[l])))
    }

    /// VF2++'s matching order of the nodes of SG.
    fn matching_order(&self) -> Vec<u32> {
        let s = &self.sg;
        let n = s.n;
        let nbrs: Vec<Vec<u32>> = (0..n).map(|v| s.neighbours(v)).collect();
        let degree: Vec<usize> = (0..n)
            .map(|v| nbrs[v].len() + nbrs[v].binary_search(&(v as u32)).is_ok() as usize)
            .collect();
        let mut rarity: Vec<i64> = self.by_label.iter().map(|x| x.len() as i64).collect();
        let mut used = vec![0usize; n];
        let mut ordered = vec![false; n];
        let mut order = Vec::with_capacity(n);
        let mut seen = vec![false; n];
        while order.len() < n {
            let root = (0..n)
                .filter(|&v| !ordered[v])
                .min_by(|&a, &b| {
                    rarity[s.label[a] as usize]
                        .cmp(&rarity[s.label[b] as usize])
                        .then(degree[b].cmp(&degree[a]))
                        .then(a.cmp(&b))
                })
                .expect("a node is left");
            let mut layer = vec![root];
            seen[root] = true;
            while !layer.is_empty() {
                let mut next = Vec::new();
                for &v in &layer {
                    for &w in &nbrs[v] {
                        if !seen[w as usize] {
                            seen[w as usize] = true;
                            next.push(w as usize);
                        }
                    }
                }
                let mut rest = layer;
                while !rest.is_empty() {
                    let (i, _) = rest
                        .iter()
                        .enumerate()
                        .min_by(|(_, &a), (_, &b)| {
                            used[b]
                                .cmp(&used[a])
                                .then(degree[b].cmp(&degree[a]))
                                .then(rarity[s.label[a] as usize].cmp(&rarity[s.label[b] as usize]))
                                .then(a.cmp(&b))
                        })
                        .expect("the layer is not empty");
                    let v = rest.swap_remove(i);
                    order.push(v as u32);
                    ordered[v] = true;
                    for &w in &nbrs[v] {
                        used[w as usize] += 1;
                    }
                    rarity[s.label[v] as usize] -= 1;
                }
                layer = next;
            }
        }
        order
    }

    /// Whether node v of FG passes the per-node tests for u of SG.
    #[inline]
    fn node_fits(&self, u: usize, v: usize) -> bool {
        let (s, f) = (&self.sg, &self.fg);
        s.label[u] == f.label[v]
            && self.sf.core[v] == NONE
            && self.edges_fit(s.loops[u], f.loops[v])
            && self.fits(s.deg[u].0, f.deg[v].0)
            && self.fits(s.deg[u].1, f.deg[v].1)
    }

    fn candidates(&self, u: usize) -> Vec<u32> {
        let (s, f) = (&self.sg, &self.fg);
        // A mapped successor n of u: the image of u must precede n's image.
        let from_out = s.out[u]
            .iter()
            .find(|&&(n, _)| n as usize != u && self.ss.core[n as usize] != NONE);
        let from_in = s.inc[u]
            .iter()
            .find(|&&(n, _)| n as usize != u && self.ss.core[n as usize] != NONE);
        let pool: Vec<u32> = if let Some(&(n, _)) = from_out {
            let m = self.ss.core[n as usize] as usize;
            let row = if f.directed { &f.inc[m] } else { &f.out[m] };
            row.iter().map(|p| p.0).collect()
        } else if let Some(&(n, _)) = from_in {
            let m = self.ss.core[n as usize] as usize;
            f.out[m].iter().map(|p| p.0).collect()
        } else {
            let pool = &self.by_label[s.label[u] as usize];
            if self.problem == Problem::Monomorphism {
                pool.clone()
            } else {
                pool.iter()
                    .copied()
                    .filter(|&v| self.sf.classes(v as usize) == 4)
                    .collect()
            }
        };
        pool.into_iter()
            .filter(|&v| self.node_fits(u, v as usize))
            .collect()
    }

    fn feasible(&mut self, u: usize, v: usize) -> bool {
        self.steps += 1;
        let (s, f) = (&self.sg, &self.fg);
        let induced = self.problem != Problem::Monomorphism;
        // Consistency with the mapped nodes.
        for &(n, c) in &s.out[u] {
            let m = self.ss.core[n as usize];
            if n as usize != u && m != NONE && !self.edges_fit(c, f.count(v, m as usize)) {
                return false;
            }
        }
        if s.directed {
            for &(n, c) in &s.inc[u] {
                let m = self.ss.core[n as usize];
                if n as usize != u && m != NONE && !self.edges_fit(c, f.count(m as usize, v)) {
                    return false;
                }
            }
        }
        if induced {
            for &(w, c) in &f.out[v] {
                let m = self.sf.core[w as usize];
                if w as usize != v && m != NONE && s.count(u, m as usize) != c {
                    return false;
                }
            }
            if f.directed {
                for &(w, c) in &f.inc[v] {
                    let m = self.sf.core[w as usize];
                    if w as usize != v && m != NONE && s.count(m as usize, u) != c {
                        return false;
                    }
                }
            }
        }
        // Cutting rules: per label, the neighbours of u in each part of
        // the frontier must fit among those of v. Slots 0..4 count u's
        // neighbours (0 = all, 1 = successors of mapped nodes, 2 =
        // predecessors of mapped nodes, 3 = neither), slots 4..8 v's. The
        // two frontier parts may overlap, so each is a necessary condition
        // for monomorphisms too; the last applies to induced problems only
        // (an edge of FG to a mapped node needs one in SG).
        let rows: [(&[Arc], &[Arc]); 2] = [(&s.out[u], &f.out[v]), (&s.inc[u], &f.inc[v])];
        let parts = if induced { 4 } else { 3 };
        for (ru, rv) in rows {
            for (row, own, base, state, side) in [(ru, u, 0, &self.ss, s), (rv, v, 4, &self.sf, f)]
            {
                for &(x, _) in row {
                    if x as usize == own {
                        continue;
                    }
                    let l = side.label[x as usize] as usize;
                    if self.scratch[l] == [0; 8] {
                        self.touched.push(l as u32);
                    }
                    let k = &mut self.scratch[l];
                    k[base] += 1;
                    let c = state.classes(x as usize);
                    for (bit, slot) in [(1u8, 1), (2, 2), (4, 3)] {
                        if c & bit != 0 {
                            k[base + slot] += 1;
                        }
                    }
                }
            }
            let mut ok = true;
            for &l in &self.touched {
                let k = &self.scratch[l as usize];
                if !(0..parts).all(|c| self.fits(k[c], k[4 + c])) {
                    ok = false;
                }
            }
            for &l in &self.touched {
                self.scratch[l as usize] = [0; 8];
            }
            self.touched.clear();
            if !ok {
                return false;
            }
        }
        true
    }

    fn push(&mut self, u: usize, v: usize) {
        let (s, f) = (&self.sg, &self.fg);
        self.ss.add(s, u, v as u32);
        self.sf.add(f, v, u as u32);
    }

    fn pop(&mut self, u: usize) {
        let v = self.ss.core[u] as usize;
        let (s, f) = (&self.sg, &self.fg);
        self.ss.remove(s, u);
        self.sf.remove(f, v);
    }

    /// The next mapping (for every node of SG, its node of FG), or `None`
    /// when there are no more.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Option<Vec<u32>> {
        if self.done {
            if self.empty_mapping {
                self.empty_mapping = false;
                return Some(Vec::new());
            }
            return None;
        }
        if !self.started {
            self.started = true;
            let cands = self.candidates(self.order[0] as usize);
            self.stack.push(Frame { cands, pos: 0 });
        }
        let ns = self.sg.n;
        while !self.stack.is_empty() {
            let d = self.stack.len() - 1;
            let top = &mut self.stack[d];
            if top.pos == top.cands.len() {
                self.stack.pop();
                if d > 0 {
                    self.pop(self.order[d - 1] as usize);
                }
                continue;
            }
            let v = top.cands[top.pos] as usize;
            top.pos += 1;
            let u = self.order[d] as usize;
            if !self.feasible(u, v) {
                continue;
            }
            if d + 1 == ns {
                let mut map = self.ss.core.clone();
                map[u] = v as u32;
                return Some(map);
            }
            self.push(u, v);
            let cands = self.candidates(self.order[d + 1] as usize);
            self.stack.push(Frame { cands, pos: 0 });
        }
        self.done = true;
        None
    }
}

/// Number of mappings of SG into FG, counting at most `limit`.
pub fn count_mappings(m: &mut Matcher, limit: u64) -> u64 {
    let mut k = 0;
    while k < limit && m.next().is_some() {
        k += 1;
    }
    k
}

// ---- Quick checks ------------------------------------------------------------

/// Whether `g1` and `g2` pass NetworkX's `could_be_isomorphic` checks:
/// equal node counts and, cumulatively, equal sorted sequences of degrees
/// (`degree`), of (degree, triangles) (`triangles`), and of (degree,
/// triangles, maximal cliques per node) (`cliques`). See the module
/// documentation.
pub fn could_be_isomorphic(
    g1: &Graph,
    g2: &Graph,
    degree: bool,
    use_triangles: bool,
    cliques: bool,
    threads: usize,
) -> Result<bool> {
    if (use_triangles || cliques) && (g1.is_directed() || g2.is_directed()) {
        return Err(Error::InvalidArgument(
            "triangle and clique checks are defined for undirected graphs only".into(),
        ));
    }
    if g1.node_count() != g2.node_count() {
        return Ok(false);
    }
    let mut props: [Vec<Vec<u64>>; 2] = [
        vec![Vec::new(); g1.node_count()],
        vec![Vec::new(); g2.node_count()],
    ];
    let consistent = |props: &[Vec<Vec<u64>>; 2]| {
        let mut a = props[0].clone();
        let mut b = props[1].clone();
        a.sort_unstable();
        b.sort_unstable();
        a == b
    };
    let graphs = [g1, g2];
    if degree {
        for (p, g) in props.iter_mut().zip(graphs) {
            for (v, row) in p.iter_mut().enumerate() {
                let loops = g
                    .out()
                    .neighbors(v)
                    .iter()
                    .filter(|&&w| w as usize == v)
                    .count();
                let d = if g.is_directed() {
                    g.out().degree(v) + g.inc().degree(v)
                } else {
                    g.out().degree(v) + loops
                };
                row.push(d as u64);
            }
        }
        if !consistent(&props) {
            return Ok(false);
        }
    }
    if use_triangles {
        for (p, g) in props.iter_mut().zip(graphs) {
            for (row, t) in p.iter_mut().zip(triangles(g, threads)) {
                row.push(t);
            }
        }
        if !consistent(&props) {
            return Ok(false);
        }
    }
    if cliques {
        for (p, g) in props.iter_mut().zip(graphs) {
            let mut count = vec![0u64; g.node_count()];
            for c in max_cliques(g, 1, threads) {
                for v in c {
                    count[v as usize] += 1;
                }
            }
            for (row, c) in p.iter_mut().zip(count) {
                row.push(c);
            }
        }
        if !consistent(&props) {
            return Ok(false);
        }
    }
    Ok(true)
}

// ---- Tree isomorphism ----------------------------------------------------------

/// Neighbour lists of a tree (direction ignored), or an error if `g` is
/// not a tree.
fn tree_adjacency(g: &Graph, what: &str) -> Result<Vec<Vec<u32>>> {
    let n = g.node_count();
    if n == 0 {
        return Err(Error::InvalidArgument(format!(
            "{what} has no nodes, so it is not a tree"
        )));
    }
    let not_tree = || Error::InvalidArgument(format!("{what} is not a tree"));
    if g.edge_count() != n - 1 {
        return Err(not_tree());
    }
    let mut adj = vec![Vec::new(); n];
    for (&s, &t) in g.edge_sources().iter().zip(g.edge_targets()) {
        if s == t {
            return Err(not_tree());
        }
        adj[s as usize].push(t);
        adj[t as usize].push(s);
    }
    // Connected with n - 1 edges.
    let mut seen = vec![false; n];
    let mut stack = vec![0usize];
    seen[0] = true;
    let mut count = 1;
    while let Some(v) = stack.pop() {
        for &w in &adj[v] {
            if !seen[w as usize] {
                seen[w as usize] = true;
                count += 1;
                stack.push(w as usize);
            }
        }
    }
    if count != n {
        return Err(not_tree());
    }
    Ok(adj)
}

/// Parent, depth, and breadth-first order of a tree rooted at `root`.
fn rooted(adj: &[Vec<u32>], root: usize) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let n = adj.len();
    let mut parent = vec![NONE; n];
    let mut depth = vec![0u32; n];
    let mut order = vec![root as u32];
    let mut seen = vec![false; n];
    seen[root] = true;
    let mut head = 0;
    while head < order.len() {
        let v = order[head] as usize;
        head += 1;
        for &w in &adj[v] {
            if !seen[w as usize] {
                seen[w as usize] = true;
                parent[w as usize] = v as u32;
                depth[w as usize] = depth[v] + 1;
                order.push(w);
            }
        }
    }
    (parent, depth, order)
}

fn rooted_iso(a1: &[Vec<u32>], r1: usize, a2: &[Vec<u32>], r2: usize) -> Option<Vec<u32>> {
    let n = a1.len();
    if a2.len() != n {
        return None;
    }
    let t = [rooted(a1, r1), rooted(a2, r2)];
    let h = t[0].1.iter().max().copied().unwrap_or(0);
    if t[1].1.iter().max().copied().unwrap_or(0) != h {
        return None;
    }
    // Canonical labels level by level, bottom up, shared by both trees: a
    // node's label is the rank of the sorted list of its children's labels
    // among the nodes of its level in either tree.
    let mut label: [Vec<u32>; 2] = [vec![0; n], vec![0; n]];
    let mut children: [Vec<Vec<u32>>; 2] = [vec![Vec::new(); n], vec![Vec::new(); n]];
    for side in 0..2 {
        for &v in &t[side].2[1..] {
            let p = t[side].0[v as usize] as usize;
            children[side][p].push(v);
        }
    }
    let mut levels: [Vec<Vec<u32>>; 2] = [
        vec![Vec::new(); h as usize + 1],
        vec![Vec::new(); h as usize + 1],
    ];
    for side in 0..2 {
        for &v in &t[side].2 {
            levels[side][t[side].1[v as usize] as usize].push(v);
        }
    }
    for d in (0..=h as usize).rev() {
        if levels[0][d].len() != levels[1][d].len() {
            return None;
        }
        let mut keys: Vec<(Vec<u32>, usize, u32)> = Vec::new();
        for side in 0..2 {
            for &v in &levels[side][d] {
                let mut k: Vec<u32> = children[side][v as usize]
                    .iter()
                    .map(|&c| label[side][c as usize])
                    .collect();
                k.sort_unstable();
                keys.push((k, side, v));
            }
        }
        keys.sort();
        let mut current = 0u32;
        for i in 0..keys.len() {
            if i > 0 && keys[i].0 != keys[i - 1].0 {
                current += 1;
            }
            label[keys[i].1][keys[i].2 as usize] = current;
        }
        // Equal multisets of labels on both sides, or no isomorphism.
        let mut l0: Vec<u32> = levels[0][d].iter().map(|&v| label[0][v as usize]).collect();
        let mut l1: Vec<u32> = levels[1][d].iter().map(|&v| label[1][v as usize]).collect();
        l0.sort_unstable();
        l1.sort_unstable();
        if l0 != l1 {
            return None;
        }
    }
    if label[0][r1] != label[1][r2] {
        return None;
    }
    let mut map = vec![NONE; n];
    let mut stack = vec![(r1, r2)];
    while let Some((v, w)) = stack.pop() {
        map[v] = w as u32;
        let mut c1 = children[0][v].clone();
        let mut c2 = children[1][w].clone();
        c1.sort_by_key(|&c| (label[0][c as usize], c));
        c2.sort_by_key(|&c| (label[1][c as usize], c));
        for (a, b) in c1.into_iter().zip(c2) {
            stack.push((a as usize, b as usize));
        }
    }
    Some(map)
}

/// An isomorphism of two rooted trees mapping `r1` onto `r2` (for every
/// node of `t1`, its image in `t2`), or `None`. See the module
/// documentation.
pub fn rooted_tree_isomorphism(
    t1: &Graph,
    r1: usize,
    t2: &Graph,
    r2: usize,
) -> Result<Option<Vec<u32>>> {
    let a1 = tree_adjacency(t1, "t1")?;
    let a2 = tree_adjacency(t2, "t2")?;
    if r1 >= a1.len() || r2 >= a2.len() {
        return Err(Error::InvalidArgument(
            "each root must be a node of its tree".into(),
        ));
    }
    Ok(rooted_iso(&a1, r1, &a2, r2))
}

/// The centre (one or two nodes) of a tree: the nodes left after leaves
/// are removed repeatedly.
fn centre(adj: &[Vec<u32>]) -> Vec<usize> {
    let n = adj.len();
    let mut deg: Vec<usize> = adj.iter().map(Vec::len).collect();
    let mut leaves: VecDeque<usize> = (0..n).filter(|&v| deg[v] <= 1).collect();
    let mut left = n;
    let mut removed = vec![false; n];
    while left > 2 {
        let k = leaves.len();
        for _ in 0..k {
            let v = leaves.pop_front().expect("leaves remain");
            removed[v] = true;
            left -= 1;
            for &w in &adj[v] {
                let w = w as usize;
                if !removed[w] {
                    deg[w] -= 1;
                    if deg[w] == 1 {
                        leaves.push_back(w);
                    }
                }
            }
        }
    }
    (0..n).filter(|&v| !removed[v]).collect()
}

/// An isomorphism of two free trees (for every node of `t1`, its image in
/// `t2`), or `None`. See the module documentation.
pub fn tree_isomorphism(t1: &Graph, t2: &Graph) -> Result<Option<Vec<u32>>> {
    let a1 = tree_adjacency(t1, "t1")?;
    let a2 = tree_adjacency(t2, "t2")?;
    if a1.len() != a2.len() {
        return Ok(None);
    }
    let mut d1: Vec<usize> = a1.iter().map(Vec::len).collect();
    let mut d2: Vec<usize> = a2.iter().map(Vec::len).collect();
    d1.sort_unstable();
    d2.sort_unstable();
    if d1 != d2 {
        return Ok(None);
    }
    let (c1, c2) = (centre(&a1), centre(&a2));
    if c1.len() != c2.len() {
        return Ok(None);
    }
    if let Some(m) = rooted_iso(&a1, c1[0], &a2, c2[0]) {
        return Ok(Some(m));
    }
    if c1.len() == 2 {
        return Ok(rooted_iso(&a1, c1[0], &a2, c2[1]));
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    fn all(m: &mut Matcher) -> Vec<Vec<u32>> {
        let mut out = Vec::new();
        while let Some(x) = m.next() {
            out.push(x);
        }
        out.sort();
        out
    }

    #[test]
    fn automorphisms_and_subgraphs_of_small_graphs() {
        let c4 = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0)], false);
        let mut m = Matcher::new(&c4, &c4, None, None, Problem::Isomorphism, false).unwrap();
        assert_eq!(all(&mut m).len(), 8);
        let p3 = graph(3, &[(0, 1), (1, 2)], false);
        // Induced paths on 3 nodes in a 4-cycle: 4 choices of the middle,
        // 2 orientations each.
        let mut m = Matcher::new(&c4, &p3, None, None, Problem::InducedSubgraph, false).unwrap();
        assert_eq!(all(&mut m).len(), 8);
        let k4 = graph(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], false);
        // In K4 no path on 3 nodes is induced, but there are 24
        // monomorphisms.
        let mut m = Matcher::new(&k4, &p3, None, None, Problem::InducedSubgraph, false).unwrap();
        assert_eq!(all(&mut m).len(), 0);
        let mut m = Matcher::new(&k4, &p3, None, None, Problem::Monomorphism, false).unwrap();
        assert_eq!(all(&mut m).len(), 24);
        // Labels restrict mappings.
        let mut m = Matcher::new(
            &c4,
            &c4,
            Some(&[0, 1, 0, 1]),
            Some(&[1, 0, 1, 0]),
            Problem::Isomorphism,
            false,
        )
        .unwrap();
        assert_eq!(all(&mut m).len(), 4);
        // Empty graphs: NetworkX's convention and the mathematical one.
        let e = graph(0, &[], false);
        let mut m = Matcher::new(&e, &e, None, None, Problem::Isomorphism, false).unwrap();
        assert!(m.next().is_none());
        let mut m = Matcher::new(&e, &e, None, None, Problem::Isomorphism, true).unwrap();
        assert_eq!(m.next(), Some(vec![]));
        assert!(m.next().is_none());
    }

    #[test]
    fn directed_cycles() {
        let c3 = graph(3, &[(0, 1), (1, 2), (2, 0)], true);
        let r3 = graph(3, &[(1, 0), (2, 1), (0, 2)], true);
        let mut m = Matcher::new(&c3, &r3, None, None, Problem::Isomorphism, false).unwrap();
        assert_eq!(all(&mut m).len(), 3);
        let t = graph(3, &[(0, 1), (1, 2), (0, 2)], true);
        let mut m = Matcher::new(&c3, &t, None, None, Problem::Isomorphism, false).unwrap();
        assert!(m.next().is_none());
    }

    #[test]
    fn trees() {
        let a = graph(5, &[(0, 1), (1, 2), (1, 3), (3, 4)], false);
        let b = graph(5, &[(4, 3), (3, 2), (3, 0), (0, 1)], false);
        let m = tree_isomorphism(&a, &b).unwrap().unwrap();
        for (&s, &t) in a.edge_sources().iter().zip(a.edge_targets()) {
            let (x, y) = (m[s as usize], m[t as usize]);
            assert!(b
                .edge_sources()
                .iter()
                .zip(b.edge_targets())
                .any(|(&p, &q)| (p == x && q == y) || (p == y && q == x)));
        }
        let star = graph(5, &[(0, 1), (0, 2), (0, 3), (0, 4)], false);
        assert!(tree_isomorphism(&a, &star).unwrap().is_none());
        assert!(rooted_tree_isomorphism(&a, 0, &b, 4).unwrap().is_some());
        assert!(rooted_tree_isomorphism(&a, 0, &b, 3).unwrap().is_none());
        let cyc = graph(3, &[(0, 1), (1, 2), (2, 0)], false);
        assert!(tree_isomorphism(&cyc, &cyc).is_err());
    }
}
