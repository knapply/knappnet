//! Regular and automorphic equivalence.
//!
//! # Regular equivalence (CATREGE)
//!
//! Nodes i and j are regularly equivalent when, for every node k tied to i,
//! some node regularly equivalent to k is tied to j in the same way, and
//! vice versa (White and Reitz 1983, `white1983graph`; Everett and Borgatti
//! 1994, `everett1994regular`). Equivalences with this property are closed
//! under union, so a largest one exists: the maximal regular equivalence.
//! It is found by CATREGE (Borgatti and Everett 1993, `borgatti1993two`):
//! start with every node in one class (or with a given partition), and
//! repeatedly split each class by the set of (tie type, class of the
//! alter) pairs its members have, until no class splits. The tie type
//! between i and k is the pair (a_ik, a_ki), so the graph and its transpose
//! are treated as one multiplex relation: a mutual tie is a different type
//! from a tie in one direction, as in sna's `redist(method = "catrege")`
//! (`butts2008sna`).
//!
//! Conventions:
//! - Without weights a_ik is 1 when there is an edge from i to k, so
//!   parallel edges count once. With weights a_ik is the total weight of
//!   the edges from i to k, and different values are different categories
//!   of tie (CATREGE treats values as categories, not magnitudes); weights
//!   may be any finite numbers other than 0, which means no tie.
//! - An undirected edge is a tie in both directions, so every tie in an
//!   undirected network has the same type.
//! - Self-loops are ignored (sna's default `diag = FALSE`).
//! - A seed partition, if given, is refined; nodes in different seed
//!   classes are never equivalent. Only equality of seed labels matters.
//! - Classes are numbered from 0 in order of their first node.
//! - In an undirected network without a seed partition, every node with at
//!   least one tie ends in the same class as every other such node, and
//!   isolated nodes form a second class: the maximal regular equivalence
//!   is degenerate there, as Borgatti and Everett note.
//!
//! # Automorphic equivalence
//!
//! Nodes are automorphically equivalent when some automorphism of the
//! network (a relabeling of the nodes that maps the network onto itself)
//! maps one onto the other (Everett 1985, `everett1985role`; Borgatti and
//! Everett 1992, `borgatti1992notions`). The classes are the orbits of the
//! automorphism group. They are computed exactly with an
//! individualization-refinement search (McKay and Piperno 2014,
//! `mckay2014practical`): colour refinement gives an equitable partition
//! that every automorphism respects; for two nodes in the same cell, a
//! backtracking search individualizes them, refines both colourings, and
//! tries the candidates for one undecided node at a time until it finds an
//! automorphism mapping one node onto the other or proves that none
//! exists. Every automorphism found merges whole orbits.
//!
//! Conventions:
//! - Defined on binary ties: weights are ignored and parallel edges count
//!   once. Direction is respected: an automorphism must preserve every
//!   arc. Self-loops are kept as a property of their node (a node with a
//!   loop can only map to a node with a loop).
//! - Orbits are numbered from 0 in order of their first node.
//! - The search is exponential in the worst case (highly symmetric
//!   networks such as strongly regular graphs), so it stops with an error
//!   once it has refined `max_steps` colourings, rather than run for an
//!   unbounded time; results are never approximate.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::union_find::UnionFind;

// ---- Regular equivalence -----------------------------------------------------

/// Options for [`regular_equivalence`].
#[derive(Debug, Clone, Copy, Default)]
pub struct RegularOptions {
    /// Treat total edge weights as tie categories.
    pub weighted: bool,
}

/// The maximal regular equivalence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegularEquivalence {
    /// Class of every node, numbered from 0 in order of first node.
    pub membership: Vec<u32>,
    /// Number of refinement rounds that split at least one class.
    pub rounds: usize,
}

/// The maximal regular equivalence by CATREGE, refining `seed` (one label
/// per node) when given. See the module documentation.
pub fn regular_equivalence(
    g: &Graph,
    seed: Option<&[i64]>,
    opts: RegularOptions,
) -> Result<RegularEquivalence> {
    let n = g.node_count();
    if let Some(s) = seed {
        if s.len() != n {
            return Err(Error::LengthMismatch {
                what: "seed partition",
                expected: n,
                found: s.len(),
            });
        }
    }
    let ties = tie_types(g, opts.weighted)?;
    let mut class: Vec<u32> = match seed {
        Some(s) => first_appearance(s.iter().copied()),
        None => vec![0; n],
    };
    let mut count = distinct(&class);
    let mut rounds = 0;
    loop {
        let mut sig: HashMap<Signature, u32> = HashMap::new();
        let mut next = Vec::with_capacity(n);
        for (i, row) in ties.iter().enumerate() {
            let mut set: Vec<(u32, u32, u32)> = row
                .iter()
                .map(|&(k, out, inc)| (out, inc, class[k as usize]))
                .collect();
            set.sort_unstable();
            set.dedup();
            let len = sig.len() as u32;
            next.push(*sig.entry((class[i], set)).or_insert(len));
        }
        let new_count = sig.len();
        if new_count == count {
            break;
        }
        class = next;
        count = new_count;
        rounds += 1;
    }
    Ok(RegularEquivalence {
        membership: first_appearance(class.iter().map(|&c| c as i64)),
        rounds,
    })
}

/// A node's class and the set of (tie out, tie in, alter's class) it has.
type Signature = (u32, Vec<(u32, u32, u32)>);

/// For every node i, its alters k (k != i, tied in either direction) with
/// the category codes of a_ik and a_ki (0 = no tie).
fn tie_types(g: &Graph, weighted: bool) -> Result<Vec<Vec<(u32, u32, u32)>>> {
    let n = g.node_count();
    let w = if weighted {
        Some(g.edge_weights().ok_or(Error::MissingWeights)?)
    } else {
        None
    };
    let mut value: HashMap<(u32, u32), f64> = HashMap::new();
    for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
        if s == t {
            continue;
        }
        let x = w.map_or(1.0, |w| w[e]);
        let mut add = |a: u32, b: u32| {
            let v = value.entry((a, b)).or_insert(0.0);
            *v = if weighted { *v + x } else { 1.0 };
        };
        add(s, t);
        if !g.is_directed() {
            add(t, s);
        }
    }
    // Category codes: 0 for no tie, then each distinct value in order.
    let mut values: Vec<f64> = value.values().copied().filter(|&x| x != 0.0).collect();
    values.sort_by(|a, b| a.total_cmp(b));
    values.dedup();
    let code = |x: f64| -> u32 {
        if x == 0.0 {
            0
        } else {
            values.binary_search_by(|v| v.total_cmp(&x)).unwrap() as u32 + 1
        }
    };
    let mut rows: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); n];
    let mut pairs: Vec<(u32, u32)> = value
        .iter()
        .filter(|(_, &x)| x != 0.0)
        .map(|(&p, _)| p)
        .collect();
    pairs.sort_unstable();
    for &(a, b) in &pairs {
        let out = code(value[&(a, b)]);
        let back = value.get(&(b, a)).copied().unwrap_or(0.0);
        let inc = code(back);
        rows[a as usize].push((b, out, inc));
        // The reverse direction is listed from b's side when b -> a has no
        // tie of its own.
        if back == 0.0 {
            rows[b as usize].push((a, 0, out));
        }
    }
    Ok(rows)
}

fn first_appearance(labels: impl Iterator<Item = i64>) -> Vec<u32> {
    let mut seen: HashMap<i64, u32> = HashMap::new();
    labels
        .map(|l| {
            let len = seen.len() as u32;
            *seen.entry(l).or_insert(len)
        })
        .collect()
}

fn distinct(labels: &[u32]) -> usize {
    let mut v = labels.to_vec();
    v.sort_unstable();
    v.dedup();
    v.len()
}

// ---- Automorphic equivalence -------------------------------------------------

/// Default cap on the number of colourings the automorphism search refines.
pub const DEFAULT_MAX_STEPS: usize = 1_000_000;

/// Orbits of the automorphism group of the network's binary ties, numbered
/// from 0 in order of first node. Fails if the search would refine more
/// than `max_steps` colourings. See the module documentation.
pub fn automorphic_equivalence(g: &Graph, max_steps: usize) -> Result<Vec<u32>> {
    let s = Simple::new(g);
    let n = s.n;
    let mut search = Search {
        s: &s,
        steps: 0,
        max_steps,
    };
    let base = search.refine(s.loops.iter().map(|&l| l as u32).collect())?;
    let mut orbit = UnionFind::new(n);
    // Representatives of the distinct orbits found so far in each cell.
    let mut reps: HashMap<u32, Vec<usize>> = HashMap::new();
    for b in 0..n {
        let cell = reps.entry(base[b]).or_default();
        let mut placed = false;
        for &a in cell.iter() {
            if orbit.same(a, b) {
                placed = true;
                break;
            }
            if let Some(perm) = search.map_node(&base, a, b)? {
                for (v, &image) in perm.iter().enumerate() {
                    orbit.union(v, image);
                }
                placed = true;
                break;
            }
        }
        if !placed {
            cell.push(b);
        }
    }
    let roots: Vec<i64> = (0..n).map(|v| orbit.find(v) as i64).collect();
    Ok(first_appearance(roots.into_iter()))
}

/// Binary ties without parallel edges: for each node, its distinct
/// neighbors other than itself with the arc type (1 = out only, 2 = in
/// only, 3 = both ways), and whether it has a self-loop.
struct Simple {
    n: usize,
    adj: Vec<Vec<(u32, u8)>>,
    loops: Vec<bool>,
}

impl Simple {
    fn new(g: &Graph) -> Simple {
        let n = g.node_count();
        let mut loops = vec![false; n];
        let mut adj: Vec<Vec<(u32, u8)>> = vec![Vec::new(); n];
        for (v, row) in adj.iter_mut().enumerate() {
            for &u in g.out().neighbors(v) {
                if u as usize == v {
                    loops[v] = true;
                } else {
                    row.push((u, 1));
                }
            }
            if g.is_directed() {
                for &u in g.inc().neighbors(v) {
                    if u as usize != v {
                        row.push((u, 2));
                    }
                }
            }
            row.sort_unstable();
            let mut merged: Vec<(u32, u8)> = Vec::with_capacity(row.len());
            for &(u, t) in row.iter() {
                match merged.last_mut() {
                    Some(last) if last.0 == u => last.1 |= t,
                    _ => merged.push((u, t)),
                }
            }
            if !g.is_directed() {
                for x in merged.iter_mut() {
                    x.1 = 3;
                }
            }
            *row = merged;
        }
        Simple { n, adj, loops }
    }

    fn arc(&self, v: usize, u: u32) -> u8 {
        match self.adj[v].binary_search_by_key(&u, |x| x.0) {
            Ok(k) => self.adj[v][k].1,
            Err(_) => 0,
        }
    }
}

struct Search<'a> {
    s: &'a Simple,
    steps: usize,
    max_steps: usize,
}

impl Search<'_> {
    /// Colour refinement to the coarsest equitable partition finer than
    /// `colors`. New colours are ranks of sorted signatures, so isomorphic
    /// coloured networks receive identical colours.
    fn refine(&mut self, mut colors: Vec<u32>) -> Result<Vec<u32>> {
        self.steps += 1;
        if self.steps > self.max_steps {
            return Err(Error::InvalidArgument(format!(
                "the automorphism search refined more than {} colourings without \
                 finishing; the network is too symmetric (or too large) for an exact \
                 answer within that limit",
                self.max_steps
            )));
        }
        let mut count = distinct(&colors);
        loop {
            let sigs: Vec<(u32, Vec<(u8, u32)>)> = (0..self.s.n)
                .map(|v| {
                    let mut nb: Vec<(u8, u32)> = self.s.adj[v]
                        .iter()
                        .map(|&(u, t)| (t, colors[u as usize]))
                        .collect();
                    nb.sort_unstable();
                    (colors[v], nb)
                })
                .collect();
            let mut sorted: Vec<&(u32, Vec<(u8, u32)>)> = sigs.iter().collect();
            sorted.sort();
            sorted.dedup();
            let next: Vec<u32> = sigs
                .iter()
                .map(|x| sorted.binary_search(&x).unwrap() as u32)
                .collect();
            let new_count = sorted.len();
            colors = next;
            if new_count == count {
                return Ok(colors);
            }
            count = new_count;
        }
    }

    /// Gives `v` a colour of its own (one more than the largest) and
    /// refines.
    fn individualize(&mut self, colors: &[u32], v: usize) -> Result<Vec<u32>> {
        let mut c = colors.to_vec();
        c[v] = colors.iter().copied().max().unwrap_or(0) + 1;
        self.refine(c)
    }

    /// An automorphism (as the image of every node) mapping `a` to `b`, if
    /// one exists. `base` is the equitable partition of the network.
    fn map_node(&mut self, base: &[u32], a: usize, b: usize) -> Result<Option<Vec<usize>>> {
        let left = self.individualize(base, a)?;
        let right = self.individualize(base, b)?;
        self.extend(left, right)
    }

    fn extend(&mut self, left: Vec<u32>, right: Vec<u32>) -> Result<Option<Vec<usize>>> {
        let n = self.s.n;
        let mut lh = left.clone();
        let mut rh = right.clone();
        lh.sort_unstable();
        rh.sort_unstable();
        if lh != rh {
            return Ok(None);
        }
        // Size of every colour class.
        let mut size: HashMap<u32, usize> = HashMap::new();
        for &c in &left {
            *size.entry(c).or_insert(0) += 1;
        }
        let pick = (0..n)
            .filter(|&v| size[&left[v]] > 1)
            .min_by_key(|&v| (size[&left[v]], left[v], v));
        match pick {
            None => {
                // Discrete: the colours define the only candidate map.
                let mut by_color = vec![0usize; n];
                for (v, &c) in right.iter().enumerate() {
                    by_color[c as usize] = v;
                }
                let perm: Vec<usize> = left.iter().map(|&c| by_color[c as usize]).collect();
                Ok(self.is_automorphism(&perm).then_some(perm))
            }
            Some(x) => {
                let c = left[x];
                let l = self.individualize(&left, x)?;
                for y in (0..n).filter(|&y| right[y] == c) {
                    let r = self.individualize(&right, y)?;
                    if let Some(perm) = self.extend(l.clone(), r)? {
                        return Ok(Some(perm));
                    }
                }
                Ok(None)
            }
        }
    }

    fn is_automorphism(&self, perm: &[usize]) -> bool {
        let s = self.s;
        (0..s.n).all(|v| {
            s.loops[v] == s.loops[perm[v]]
                && s.adj[v].len() == s.adj[perm[v]].len()
                && s.adj[v]
                    .iter()
                    .all(|&(u, t)| s.arc(perm[v], perm[u as usize] as u32) == t)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let s: Vec<i64> = e.iter().map(|x| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    #[test]
    fn regular_equivalence_of_a_directed_tree() {
        // 0 -> 1, 0 -> 2, 1 -> 3, 2 -> 4: roots, middles, and leaves.
        let g = graph(5, &[(0, 1), (0, 2), (1, 3), (2, 4)], true);
        let r = regular_equivalence(&g, None, Default::default()).unwrap();
        assert_eq!(r.membership, vec![0, 1, 1, 2, 2]);
        // Unbalanced: 0 -> 1 -> 3 and 0 -> 2. Round 1 separates 1 (in and
        // out ties) from the leaves 2 and 3; round 2 separates the leaves,
        // whose senders are now in different classes.
        let g = graph(4, &[(0, 1), (0, 2), (1, 3)], true);
        let r = regular_equivalence(&g, None, Default::default()).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2, 3]);
        assert!(r.rounds >= 2);
    }

    #[test]
    fn undirected_regular_equivalence_is_degenerate() {
        let g = graph(5, &[(0, 1), (1, 2), (2, 3)], false);
        let r = regular_equivalence(&g, None, Default::default()).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 0, 1]);
        // A seed partition is refined: marking the end node 0 of the path
        // 0-1-2-3 separates every position along it.
        let r = regular_equivalence(&g, Some(&[5, 7, 7, 7, 7]), Default::default()).unwrap();
        assert_eq!(r.membership, vec![0, 1, 2, 3, 4]);
        // Marking both ends keeps the symmetric positions together.
        let r = regular_equivalence(&g, Some(&[5, 7, 7, 5, 7]), Default::default()).unwrap();
        assert_eq!(r.membership, vec![0, 1, 1, 0, 2]);
        assert!(regular_equivalence(&g, Some(&[0, 1]), Default::default()).is_err());
    }

    #[test]
    fn mutual_ties_are_their_own_type() {
        // 0 <-> 1 versus 2 -> 3 and 3 -> 2? 0 and 1 have one mutual tie;
        // 2 -> 3 only: different types.
        let g = graph(4, &[(0, 1), (1, 0), (2, 3)], true);
        let r = regular_equivalence(&g, None, Default::default()).unwrap();
        assert_eq!(r.membership, vec![0, 0, 1, 2]);
        // Values as categories.
        let w = Graph::from_edges(4, &[0, 2], &[1, 3], Some(&[1.0, 2.0]), false).unwrap();
        let r = regular_equivalence(&w, None, RegularOptions { weighted: true }).unwrap();
        assert_eq!(r.membership, vec![0, 0, 1, 1]);
        let r = regular_equivalence(&w, None, RegularOptions::default()).unwrap();
        assert_eq!(r.membership, vec![0, 0, 0, 0]);
    }

    #[test]
    fn automorphism_orbits() {
        // Path 0-1-2-3-4: orbits {0,4}, {1,3}, {2}.
        let g = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4)], false);
        assert_eq!(
            automorphic_equivalence(&g, DEFAULT_MAX_STEPS).unwrap(),
            vec![0, 1, 2, 1, 0]
        );
        // Directed path: no symmetry.
        let d = graph(3, &[(0, 1), (1, 2)], true);
        assert_eq!(
            automorphic_equivalence(&d, DEFAULT_MAX_STEPS).unwrap(),
            vec![0, 1, 2]
        );
        // Cycle of 6 with a loop at 0: reflections fixing 0 remain.
        let c = graph(
            6,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 0)],
            false,
        );
        assert_eq!(
            automorphic_equivalence(&c, DEFAULT_MAX_STEPS).unwrap(),
            vec![0, 1, 2, 3, 2, 1]
        );
        // Two disjoint triangles: one orbit; plus an isolate.
        let t = graph(7, &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)], false);
        assert_eq!(
            automorphic_equivalence(&t, DEFAULT_MAX_STEPS).unwrap(),
            vec![0, 0, 0, 0, 0, 0, 1]
        );
        // Colour refinement alone cannot separate these: a 6-cycle and two
        // triangles are both 2-regular, and their union has two orbits.
        let u = graph(
            12,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 0),
                (6, 7),
                (7, 8),
                (8, 6),
                (9, 10),
                (10, 11),
                (11, 9),
            ],
            false,
        );
        assert_eq!(
            automorphic_equivalence(&u, DEFAULT_MAX_STEPS).unwrap(),
            vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
        );
        // The step cap is enforced.
        assert!(automorphic_equivalence(&u, 2).is_err());
    }
}
