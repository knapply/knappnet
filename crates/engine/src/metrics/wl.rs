//! Colour refinement and Weisfeiler-Lehman graph hashes.
//!
//! # Colour refinement
//!
//! One round of colour refinement gives every node a new colour determined
//! by its current colour and the collection of (tie type, colour) pairs of
//! its neighbours (Weisfeiler and Leman 1968, `weisfeiler1968reduction`;
//! the 1-dimensional algorithm). [`Refiner`] runs rounds over any notion of
//! neighbourhood a caller supplies, with the collection read as a multiset
//! (Weisfeiler-Lehman, and the subtree kernel of Shervashidze et al. 2011,
//! `shervashidze2011weisfeiler`) or as a set (CATREGE's regular
//! equivalence, `borgatti1993two`, in [`crate::metrics::regular`]). Colours
//! are numbered in order of first appearance in a dictionary the caller
//! may keep across rounds and across graphs, so colours of different
//! graphs are comparable, as graph kernels and incremental updates need.
//! [`wl_colors`] runs the multiset rounds on a graph.
//!
//! # Weisfeiler-Lehman graph and subgraph hashes
//!
//! [`wl_graph_hash`] and [`wl_subgraph_hashes`] reproduce, character for
//! character, NetworkX 3.7's `weisfeiler_lehman_graph_hash` and
//! `weisfeiler_lehman_subgraph_hashes` (`hagberg2008exploring`), which
//! follow Shervashidze et al. (2011): every node starts with a label (its
//! degree, `in_out` for directed graphs, or a caller's string), and each
//! iteration replaces a node's label by the BLAKE2b hash
//! ([`crate::blake2b`]) of its label followed by its neighbours' labels
//! sorted and concatenated (successors prefixed `s_`, then predecessors
//! prefixed `p_`, for directed graphs; an edge label, when given, replaces
//! the prefix). The graph hash is the hash of the Python text of the tuple
//! of (label, count) pairs of every iteration, each iteration's pairs
//! sorted by label. Without node or edge labels the degree labels count as
//! the first iteration, as in NetworkX.
//!
//! Conventions:
//! - Self-loops: a looped node is its own neighbour (once), and its loop
//!   adds 2 to its undirected degree, or 1 to each of its in- and
//!   out-degrees, as in a NetworkX `Graph` or `DiGraph`.
//! - Parallel edges are merged (NetworkX does not hash multigraphs);
//!   with edge labels they are an error, since a pair of nodes would have
//!   two labels.
//! - Direction is respected in directed graphs.
//! - Weights are ignored; edge labels are strings the caller supplies
//!   (`str()` of the attribute in Python).
//! - Labels must be ASCII (NetworkX encodes them as ASCII and fails
//!   otherwise).
//! - The hashes are the same for every thread count.

use std::collections::HashMap;

use crate::bitset::NeighborSets;
use crate::blake2b;
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;

// ---- Colour refinement ------------------------------------------------------

/// How a node's neighbour colours are collected in a refinement round.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Aggregate {
    /// Each (tie type, colour) pair counts once (CATREGE).
    Set,
    /// Each neighbour counts (Weisfeiler-Lehman).
    Multiset,
}

/// A node's colour and its sorted collection of (tie type, colour) pairs.
type Signature = (u32, Vec<(u64, u32)>);

/// Colour refinement rounds with a dictionary of signatures. See the
/// module documentation.
#[derive(Debug, Clone, Default)]
pub struct Refiner {
    dict: HashMap<Signature, u32>,
    buf: Vec<(u64, u32)>,
}

impl Refiner {
    /// A refiner with an empty dictionary.
    pub fn new() -> Refiner {
        Refiner::default()
    }

    /// Forgets every colour.
    pub fn clear(&mut self) {
        self.dict.clear();
    }

    /// Number of colours in the dictionary.
    pub fn len(&self) -> usize {
        self.dict.len()
    }

    /// Whether the dictionary is empty.
    pub fn is_empty(&self) -> bool {
        self.dict.is_empty()
    }

    /// One round: `row(i, buf)` pushes a (tie type, neighbour) pair onto
    /// `buf` for every neighbour of node i (the same neighbour may appear
    /// with several tie types). Returns the new colour of every node: the
    /// dictionary's number for (colour of i, collection of (tie type,
    /// colour of neighbour)), new signatures numbered in node order after
    /// those already in the dictionary.
    pub fn round(
        &mut self,
        class: &[u32],
        aggregate: Aggregate,
        mut row: impl FnMut(usize, &mut Vec<(u64, u32)>),
    ) -> Vec<u32> {
        let mut next = Vec::with_capacity(class.len());
        for (i, &c) in class.iter().enumerate() {
            self.buf.clear();
            row(i, &mut self.buf);
            for p in self.buf.iter_mut() {
                p.1 = class[p.1 as usize];
            }
            self.buf.sort_unstable();
            if aggregate == Aggregate::Set {
                self.buf.dedup();
            }
            let len = self.dict.len() as u32;
            if let Some(&k) = self.dict.get(&(c, self.buf.clone())) {
                next.push(k);
            } else {
                self.dict.insert((c, self.buf.clone()), len);
                next.push(len);
            }
        }
        next
    }
}

/// Weisfeiler-Lehman colours of `g` after each of `rounds` rounds
/// (multiset aggregation; for directed graphs successors and predecessors
/// are different tie types; parallel edges count once and a looped node
/// is its own neighbour), starting from `initial` labels (only equality
/// matters) or from one colour for all. Colours are numbered in order of
/// first appearance, continuing across rounds through `refiner`, so the
/// colours of graphs refined with the same refiner are comparable.
pub fn wl_colors(
    g: &Graph,
    initial: Option<&[i64]>,
    rounds: usize,
    refiner: &mut Refiner,
) -> Result<Vec<Vec<u32>>> {
    let n = g.node_count();
    let mut class: Vec<u32> = match initial {
        Some(s) if s.len() != n => {
            return Err(Error::LengthMismatch {
                what: "initial labels",
                expected: n,
                found: s.len(),
            })
        }
        Some(s) => {
            let mut seen: HashMap<i64, u32> = HashMap::new();
            s.iter()
                .map(|&x| {
                    let k = seen.len() as u32;
                    *seen.entry(x).or_insert(k)
                })
                .collect()
        }
        None => vec![0; n],
    };
    let adj = Directed::new(g, None)?;
    let mut out = Vec::with_capacity(rounds);
    for _ in 0..rounds {
        class = refiner.round(&class, Aggregate::Multiset, |i, buf| {
            buf.extend(adj.succ(i).iter().map(|&(v, _)| (1, v)));
            if g.is_directed() {
                buf.extend(adj.pred(i).iter().map(|&(v, _)| (2, v)));
            }
        });
        out.push(class.clone());
    }
    Ok(out)
}

// ---- Weisfeiler-Lehman hashes (NetworkX) ------------------------------------

/// Distinct successors and predecessors of every node with the edge each
/// pair came from (undirected graphs: neighbours, `pred` empty).
struct Directed {
    succ: Vec<Vec<(u32, u32)>>,
    pred: Vec<Vec<(u32, u32)>>,
}

impl Directed {
    /// With `labelled`, parallel edges are an error.
    fn new(g: &Graph, labelled: Option<&[String]>) -> Result<Directed> {
        let n = g.node_count();
        let mut succ: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
        let mut pred: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
        for (e, (&s, &t)) in g.edge_sources().iter().zip(g.edge_targets()).enumerate() {
            succ[s as usize].push((t, e as u32));
            if g.is_directed() {
                pred[t as usize].push((s, e as u32));
            } else if s != t {
                succ[t as usize].push((s, e as u32));
            }
        }
        for rows in [&mut succ, &mut pred] {
            for r in rows.iter_mut() {
                r.sort_unstable();
                let before = r.len();
                r.dedup_by_key(|p| p.0);
                if labelled.is_some() && r.len() != before {
                    return Err(Error::InvalidArgument(
                        "edge labels need a network without parallel edges".into(),
                    ));
                }
            }
        }
        Ok(Directed { succ, pred })
    }

    fn succ(&self, v: usize) -> &[(u32, u32)] {
        &self.succ[v]
    }

    fn pred(&self, v: usize) -> &[(u32, u32)] {
        &self.pred[v]
    }
}

/// Options for [`wl_graph_hash`] and [`wl_subgraph_hashes`].
#[derive(Debug, Clone, Copy)]
pub struct WlOptions<'a> {
    /// Number of iterations, at least 1 (NetworkX's default is 3).
    pub iterations: usize,
    /// Bytes per BLAKE2b digest, 1 to 64 (NetworkX's default is 16).
    pub digest_size: usize,
    /// Initial node labels (`str()` of the node attribute), one per node.
    pub node_labels: Option<&'a [String]>,
    /// Edge labels (`str()` of the edge attribute), one per edge.
    pub edge_labels: Option<&'a [String]>,
    /// Subgraph hashes only: start every node's list with the hash of its
    /// initial label.
    pub include_initial_labels: bool,
    /// Threads (0 = all cores).
    pub threads: usize,
}

impl Default for WlOptions<'_> {
    fn default() -> Self {
        WlOptions {
            iterations: 3,
            digest_size: 16,
            node_labels: None,
            edge_labels: None,
            include_initial_labels: false,
            threads: 0,
        }
    }
}

struct Hasher<'a> {
    adj: Directed,
    directed: bool,
    opts: WlOptions<'a>,
}

impl Hasher<'_> {
    fn new<'a>(g: &Graph, opts: WlOptions<'a>) -> Result<(Hasher<'a>, Vec<String>, usize)> {
        let n = g.node_count();
        if opts.iterations == 0 {
            return Err(Error::InvalidArgument(
                "the Weisfeiler-Lehman algorithm needs at least one iteration".into(),
            ));
        }
        if !(1..=64).contains(&opts.digest_size) {
            return Err(Error::InvalidArgument(
                "digest_size must be between 1 and 64 bytes".into(),
            ));
        }
        if let Some(l) = opts.node_labels {
            if l.len() != n {
                return Err(Error::LengthMismatch {
                    what: "node labels",
                    expected: n,
                    found: l.len(),
                });
            }
        }
        if let Some(l) = opts.edge_labels {
            if l.len() != g.edge_count() {
                return Err(Error::LengthMismatch {
                    what: "edge labels",
                    expected: g.edge_count(),
                    found: l.len(),
                });
            }
        }
        for l in opts
            .node_labels
            .into_iter()
            .chain(opts.edge_labels)
            .flatten()
        {
            if !l.is_ascii() {
                return Err(Error::InvalidArgument(format!(
                    "label {l:?} is not ASCII; Weisfeiler-Lehman hashes need ASCII labels"
                )));
            }
        }
        let adj = Directed::new(g, opts.edge_labels)?;
        let mut iterations = opts.iterations;
        let labels: Vec<String> = if let Some(l) = opts.node_labels {
            l.to_vec()
        } else if opts.edge_labels.is_some() {
            vec![String::new(); n]
        } else {
            iterations -= 1;
            if g.is_directed() {
                (0..n)
                    .map(|v| format!("{}_{}", adj.pred(v).len(), adj.succ(v).len()))
                    .collect()
            } else {
                let sets = NeighborSets::new(g, false);
                (0..n).map(|v| sets.simple_degree(v).to_string()).collect()
            }
        };
        let h = Hasher {
            adj,
            directed: g.is_directed(),
            opts,
        };
        Ok((h, labels, iterations))
    }

    /// The text a node's next label is the hash of.
    fn aggregate(&self, v: usize, labels: &[String]) -> String {
        let edge = self.opts.edge_labels;
        let part = |rows: &[(u32, u32)], prefix: &str| -> String {
            let mut items: Vec<String> = rows
                .iter()
                .map(|&(u, e)| {
                    let p = edge.map_or(prefix, |l| l[e as usize].as_str());
                    let mut s = String::with_capacity(p.len() + labels[u as usize].len());
                    s.push_str(p);
                    s.push_str(&labels[u as usize]);
                    s
                })
                .collect();
            items.sort_unstable();
            items.concat()
        };
        let mut s = labels[v].clone();
        if self.directed {
            s.push_str(&part(self.adj.succ(v), "s_"));
            s.push_str(&part(self.adj.pred(v), "p_"));
        } else {
            s.push_str(&part(self.adj.succ(v), ""));
        }
        s
    }

    fn step(&self, labels: &[String]) -> Vec<String> {
        let n = labels.len();
        // A hash costs some thousand steps, plus the sort of the
        // neighbour labels.
        let arcs: usize = self
            .adj
            .succ
            .iter()
            .chain(&self.adj.pred)
            .map(Vec::len)
            .sum();
        let work = (n as u64 * 2048).saturating_add(arcs as u64 * 256);
        let threads = parallel::threads_for(self.opts.threads, n, work);
        parallel::map_items(
            n,
            threads,
            || (),
            |_, v| blake2b::hex(self.aggregate(v, labels).as_bytes(), self.opts.digest_size),
        )
    }
}

/// The Weisfeiler-Lehman graph hash, as NetworkX computes it. See the
/// module documentation.
pub fn wl_graph_hash(g: &Graph, opts: WlOptions) -> Result<String> {
    let (h, mut labels, iterations) = Hasher::new(g, opts)?;
    let mut items: Vec<String> = Vec::new();
    for _ in 0..iterations {
        labels = h.step(&labels);
        let mut sorted: Vec<&String> = labels.iter().collect();
        sorted.sort_unstable();
        let mut i = 0;
        while i < sorted.len() {
            let mut j = i;
            while j < sorted.len() && sorted[j] == sorted[i] {
                j += 1;
            }
            items.push(format!("('{}', {})", sorted[i], j - i));
            i = j;
        }
    }
    // Python's str() of a tuple: a one-element tuple keeps its comma.
    let text = match items.len() {
        0 => "()".to_string(),
        1 => format!("({},)", items[0]),
        _ => format!("({})", items.join(", ")),
    };
    Ok(blake2b::hex(text.as_bytes(), opts.digest_size))
}

/// The Weisfeiler-Lehman subgraph hashes of every node (one list per
/// node, oldest first), as NetworkX computes them. See the module
/// documentation.
pub fn wl_subgraph_hashes(g: &Graph, opts: WlOptions) -> Result<Vec<Vec<String>>> {
    let (h, mut labels, iterations) = Hasher::new(g, opts)?;
    let size = opts.digest_size;
    let mut out: Vec<Vec<String>> = vec![Vec::new(); labels.len()];
    if opts.include_initial_labels {
        for (o, l) in out.iter_mut().zip(&labels) {
            o.push(blake2b::hex(l.as_bytes(), size));
        }
    }
    if iterations < opts.iterations {
        // Degree labels count as the first iteration.
        for (o, l) in out.iter_mut().zip(&labels) {
            o.push(blake2b::hex(l.as_bytes(), size));
        }
    }
    for _ in 0..iterations {
        labels = h.step(&labels);
        for (o, l) in out.iter_mut().zip(&labels) {
            o.push(l.clone());
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hashes_respect_isomorphism() {
        // A path and a star on three nodes are isomorphic; a triangle is
        // not. The ends of the path have equal subgraph hashes.
        let path = Graph::from_edges(3, &[0, 1], &[1, 2], None, false).unwrap();
        let star = Graph::from_edges(3, &[1, 1], &[0, 2], None, false).unwrap();
        let o = WlOptions {
            threads: 1,
            ..Default::default()
        };
        assert_eq!(
            wl_graph_hash(&path, o).unwrap(),
            wl_graph_hash(&star, o).unwrap()
        );
        let tri = Graph::from_edges(3, &[0, 1, 2], &[1, 2, 0], None, false).unwrap();
        assert_ne!(
            wl_graph_hash(&path, o).unwrap(),
            wl_graph_hash(&tri, o).unwrap()
        );
        let sub = wl_subgraph_hashes(&path, o).unwrap();
        assert_eq!(sub[0], sub[2]);
        assert_eq!(sub[0].len(), 3);
    }

    #[test]
    fn refinement_of_a_path() {
        // Path 0-1-2-3-4: rounds split ends, then neighbours of ends.
        let g = Graph::from_edges(5, &[0, 1, 2, 3], &[1, 2, 3, 4], None, false).unwrap();
        let mut r = Refiner::new();
        let c = wl_colors(&g, None, 3, &mut r).unwrap();
        assert_eq!(c[0], [0, 1, 1, 1, 0]);
        assert_eq!(c[1], [2, 3, 4, 3, 2]);
        assert_eq!(c[2], [5, 6, 7, 6, 5]);
    }
}
