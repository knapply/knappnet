//! Triad census.
//!
//! Counts all n-choose-3 node triples by their isomorphism class among the 16
//! directed triad types of Holland and Leinhardt (1976, `holland1976local`),
//! using the subquadratic algorithm of Batagelj and Mrvar (2001,
//! `batagelj2001subquadratic`). Only connected triads are enumerated; triads
//! with exactly one dyad are counted with a closed-form expression, and the
//! empty triad 003 is obtained by subtraction from n-choose-3.
//!
//! Conventions:
//! - The census is defined on simple directed graphs. Self-loops are ignored
//!   and parallel edges collapse to a single arc, as in igraph, NetworkX,
//!   and sna.
//! - An undirected graph is treated as a directed graph in which every edge
//!   is mutual, as in igraph. Only types 003, 102, 201, and 300 can then be
//!   non-zero.
//! - Weights are ignored: the census is defined on binary ties. Dense
//!   networks use a bit matrix for constant-time arc tests
//!   ([`crate::bitset`]); sparse networks use sorted adjacency lists.
//!
//! Output order is the standard MAN labeling:
//! 003, 012, 102, 021D, 021U, 021C, 111D, 111U,
//! 030T, 030C, 201, 120D, 120U, 120C, 210, 300.

use crate::bitset::BitMatrix;
use crate::graph::Graph;
use crate::parallel;

/// Triad type labels in the order returned by [`triad_census`].
pub const TRIAD_NAMES: [&str; 16] = [
    "003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T", "030C", "201", "120D",
    "120U", "120C", "210", "300",
];

/// Maps the 6-bit arc code of an ordered triple (v, u, w) to a triad type
/// (1-based index into `TRIAD_NAMES`). Bits: 1 = v->u, 2 = u->v, 4 = v->w,
/// 8 = w->v, 16 = u->w, 32 = w->u. Table from Batagelj and Mrvar (2001).
const TRICODES: [u8; 64] = [
    1, 2, 2, 3, 2, 4, 6, 8, 2, 6, 5, 7, 3, 8, 7, 11, 2, 6, 4, 8, 5, 9, 9, 13, 6, 10, 9, 14, 7, 14,
    12, 15, 2, 5, 6, 7, 6, 9, 10, 14, 4, 9, 9, 12, 8, 13, 14, 15, 3, 7, 8, 11, 7, 12, 14, 15, 8,
    14, 13, 15, 11, 15, 15, 16,
];

/// Simplified adjacency: sorted, de-duplicated, loop-free.
struct Simple {
    /// Out-neighbors of each node.
    succ: Vec<Vec<u32>>,
    /// Union of in- and out-neighbors of each node.
    nbrs: Vec<Vec<u32>>,
    /// Arc membership as bits, for dense networks.
    bits: Option<BitMatrix>,
}

impl Simple {
    fn new(g: &Graph, use_bits: Option<bool>) -> Simple {
        let n = g.node_count();
        let clean = |v: usize, xs: &[u32]| -> Vec<u32> {
            let mut out: Vec<u32> = xs.iter().copied().filter(|&x| x as usize != v).collect();
            out.dedup(); // input rows are sorted
            out
        };
        let succ: Vec<Vec<u32>> = (0..n).map(|v| clean(v, g.out().neighbors(v))).collect();
        let nbrs = if g.is_directed() {
            (0..n)
                .map(|v| {
                    let pred = clean(v, g.inc().neighbors(v));
                    merge_union(&succ[v], &pred)
                })
                .collect()
        } else {
            succ.clone()
        };
        let arcs: usize = succ.iter().map(Vec::len).sum();
        let bits = use_bits
            .unwrap_or_else(|| BitMatrix::worthwhile(n, arcs))
            .then(|| {
                let mut m = BitMatrix::new(n);
                for (v, row) in succ.iter().enumerate() {
                    for &w in row {
                        m.set(v, w as usize);
                    }
                }
                m
            });
        Simple { succ, nbrs, bits }
    }

    #[inline]
    fn arc(&self, a: usize, b: u32) -> bool {
        match &self.bits {
            Some(m) => m.get(a, b as usize),
            None => self.succ[a].binary_search(&b).is_ok(),
        }
    }

    #[inline]
    fn adjacent(&self, a: usize, b: u32) -> bool {
        match &self.bits {
            Some(m) => m.get(a, b as usize) || m.get(b as usize, a),
            None => self.nbrs[a].binary_search(&b).is_ok(),
        }
    }
}

fn merge_union(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Returns the 16 triad counts in MAN order. Counts are exact; the 003 count
/// is computed in 128-bit arithmetic and returned as f64, which is exact up
/// to 2^53 and correctly rounded beyond that.
pub fn triad_census(g: &Graph, threads: usize) -> [f64; 16] {
    census_impl(g, threads, None)
}

fn census_impl(g: &Graph, threads: usize, use_bits: Option<bool>) -> [f64; 16] {
    let n = g.node_count();
    let s = Simple::new(g, use_bits);

    let partials = parallel::for_each_item(
        n,
        threads,
        || ([0u64; 16], Vec::<u32>::new()),
        |(census, union), v| {
            for &u in &s.nbrs[v] {
                let u_idx = u as usize;
                if u_idx <= v {
                    continue;
                }
                // S = N(u) ∪ N(v) \ {u, v}
                union.clear();
                union.extend(
                    merge_union(&s.nbrs[v], &s.nbrs[u_idx])
                        .into_iter()
                        .filter(|&x| x as usize != v && x != u),
                );
                let dyad = if s.arc(v, u) && s.arc(u_idx, v as u32) {
                    2
                } else {
                    1
                };
                census[dyad] += (n - union.len() - 2) as u64;
                for &w in union.iter() {
                    let w_idx = w as usize;
                    if u_idx < w_idx || (v < w_idx && w_idx < u_idx && !s.adjacent(v, w)) {
                        let code = tricode(&s, v, u_idx, w_idx);
                        census[TRICODES[code] as usize - 1] += 1;
                    }
                }
            }
        },
    );

    let mut counts = [0u64; 16];
    for (part, _) in partials {
        for (c, p) in counts.iter_mut().zip(part) {
            *c += p;
        }
    }

    let n128 = n as u128;
    let total = if n < 3 {
        0
    } else {
        n128 * (n128 - 1) * (n128 - 2) / 6
    };
    let connected: u128 = counts[1..].iter().map(|&c| c as u128).sum();
    let mut out = [0f64; 16];
    out[0] = (total - connected) as f64;
    for i in 1..16 {
        out[i] = counts[i] as f64;
    }
    out
}

#[inline]
fn tricode(s: &Simple, v: usize, u: usize, w: usize) -> usize {
    let (vu, uv) = (s.arc(v, u as u32), s.arc(u, v as u32));
    let (vw, wv) = (s.arc(v, w as u32), s.arc(w, v as u32));
    let (uw, wu) = (s.arc(u, w as u32), s.arc(w, u as u32));
    (vu as usize)
        | (uv as usize) << 1
        | (vw as usize) << 2
        | (wv as usize) << 3
        | (uw as usize) << 4
        | (wu as usize) << 5
}

#[cfg(test)]
mod tests {
    use super::*;

    fn census(n: usize, s: &[i64], t: &[i64], directed: bool) -> [f64; 16] {
        triad_census(&Graph::from_edges(n, s, t, None, directed).unwrap(), 1)
    }

    fn idx(name: &str) -> usize {
        TRIAD_NAMES.iter().position(|&x| x == name).unwrap()
    }

    #[test]
    fn single_triads() {
        let cases: &[(&str, &[i64], &[i64])] = &[
            ("003", &[], &[]),
            ("012", &[0], &[1]),
            ("102", &[0, 1], &[1, 0]),
            ("021D", &[0, 0], &[1, 2]),
            ("021U", &[1, 2], &[0, 0]),
            ("021C", &[0, 1], &[1, 2]),
            ("030T", &[0, 0, 1], &[1, 2, 2]),
            ("030C", &[0, 1, 2], &[1, 2, 0]),
            ("300", &[0, 1, 0, 2, 1, 2], &[1, 0, 2, 0, 2, 1]),
        ];
        for (name, s, t) in cases {
            let c = census(3, s, t, true);
            assert_eq!(c[idx(name)], 1.0, "{name}: {c:?}");
            assert_eq!(c.iter().sum::<f64>(), 1.0, "{name}");
        }
    }

    #[test]
    fn sums_to_n_choose_3_and_ignores_loops_and_multi_edges() {
        let s = [0, 0, 1, 2, 3, 3, 4, 4, 5];
        let t = [1, 1, 2, 0, 3, 4, 0, 5, 1];
        let c = census(7, &s, &t, true);
        assert_eq!(c.iter().sum::<f64>(), 35.0);
        let simple = census(7, &[0, 1, 2, 3, 4, 4, 5], &[1, 2, 0, 4, 0, 5, 1], true);
        assert_eq!(c, simple);
    }

    #[test]
    fn bit_matrix_and_lists_agree() {
        // A dense random digraph, counted both ways.
        let mut state = 12345u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % 40
        };
        let (mut s, mut t) = (Vec::new(), Vec::new());
        for _ in 0..600 {
            s.push(next() as i64);
            t.push(next() as i64);
        }
        let g = Graph::from_edges(40, &s, &t, None, true).unwrap();
        assert_eq!(
            census_impl(&g, 1, Some(true)),
            census_impl(&g, 1, Some(false))
        );
    }

    #[test]
    fn undirected_edges_are_mutual() {
        let c = census(4, &[0, 1, 2], &[1, 2, 0], false);
        assert_eq!(c[idx("300")], 1.0);
        assert_eq!(c[idx("102")], 3.0);
        assert_eq!(c.iter().sum::<f64>(), 4.0);
    }
}
