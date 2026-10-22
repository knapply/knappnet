//! Brokerage: Burt's structural-hole measures and Gould-Fernandez roles.
//!
//! # Structural holes
//!
//! Burt (1992, `burt1992structural`) measures how far a node's contacts are
//! redundant with one another. Both measures start from the *mutual tie*
//! z_ij = a_ij + a_ji, where a_ij is the total weight of the edges from i to
//! j (the number of such edges in a binary network), and from the
//! *proportional tie strength*
//!
//! ```text
//! p_ij = z_ij / sum_k z_ik,
//! ```
//!
//! the share of i's network time and energy invested in j.
//!
//! - **Constraint** (Burt 1992, ch. 2; Burt 2004, `burt2004structural`):
//!   C_i = sum_j (p_ij + sum_q p_iq p_qj)^2, where j ranges over i's
//!   contacts and q over i's contacts other than j. Low constraint means
//!   many structural holes around i.
//! - **Effective size** (Burt 1992, ch. 2): E_i = sum_j (1 - sum_q p_iq
//!   m_jq), with j and q as above and m_jq = z_jq / max_k z_jk, j's tie to q
//!   relative to j's strongest tie. In a binary undirected network without
//!   parallel edges this reduces to n - 2t/n, where n is the number of
//!   contacts and t the number of ties among them (Borgatti 1997,
//!   `borgatti1997structural`); the unit tests check this identity.
//! - **Efficiency** (Burt 1992, ch. 2): effective size divided by the number
//!   of contacts.
//!
//! Conventions, matching igraph's `constraint()` and NetworkX's
//! `constraint()` and `effective_size()`:
//! - A node's contacts are its in- and out-neighbors, other than itself.
//!   Direction only matters through z_ij = a_ij + a_ji, so a mutual pair of
//!   arcs is a stronger tie than a single arc, and an undirected network
//!   gives the same results as the directed network with every edge
//!   reciprocated.
//! - Self-loops are ignored. (NetworkX includes them in the row sums of
//!   z, so its values differ on networks with loops; igraph ignores them.)
//! - Parallel edges add up: in a binary network each edge has weight 1, so
//!   two parallel edges make a tie twice as strong as one, as in igraph and
//!   NetworkX.
//! - With weights, they are read as tie strengths (volume of contact) and
//!   used as given; they must be strictly positive. Whether the network
//!   declares its weights as distances or strengths does not change this.
//! - The sum over j covers i's contacts only, as in igraph and NetworkX.
//!   (Burt's formula is sometimes written with j over all nodes; for a
//!   non-contact j the term p_ij is zero but indirect terms may not be.)
//! - All three measures are undefined, and returned as NaN, for a node
//!   without contacts (an isolate, or a node with only self-loops). igraph
//!   returns 0 for constraint on a node with only self-loops.
//!
//! # Gould-Fernandez brokerage
//!
//! Gould and Fernandez (1989, `gould1989structures`): node j brokers the
//! ordered pair (i, k) when i -> j -> k and there is no arc i -> k, with i,
//! j, k distinct. Given a group for each node, each brokered pair is one of
//! five roles of j:
//!
//! | Role | Groups of (i, j, k) | sna name |
//! |---|---|---|
//! | coordinator | A, A, A | `w_I` |
//! | itinerant (consultant) | A, B, A | `w_O` |
//! | representative | A, A, B | `b_IO` |
//! | gatekeeper | A, B, B | `b_OI` |
//! | liaison | A, B, C | `b_O` |
//!
//! Conventions, matching the raw counts of `sna::brokerage()` (`raw.nli`):
//! - The measure is defined on binary ties: weights are ignored and
//!   parallel edges count once.
//! - Self-loops are ignored.
//! - An undirected edge is a pair of mutual arcs, so in an undirected
//!   network each unordered pair {i, k} that j brokers is counted twice,
//!   once in each direction, as sna does with a symmetric adjacency matrix.
//! - Pairs are ordered: with i <-> j <-> k and no tie between i and k, j
//!   brokers both (i, k) and (k, i).

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::matrix::Dense;
use crate::parallel;

/// Options for [`constraint`], [`effective_size`], and [`efficiency`].
#[derive(Debug, Clone, Copy, Default)]
pub struct StructuralHoleOptions {
    /// Use edge weights as tie strengths (they must be strictly positive).
    pub weighted: bool,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Burt's constraint of every node, in node order. NaN for nodes without
/// contacts.
///
/// See the module documentation for conventions.
pub fn constraint(g: &Graph, opts: StructuralHoleOptions) -> Result<Vec<f64>> {
    let ties = MutualTies::new(g, opts.weighted)?;
    Ok(per_node(&ties, opts.threads, |ties, acc, touched, i| {
        ties.constraint(acc, touched, i)
    }))
}

/// Burt's effective size of every node, in node order. NaN for nodes
/// without contacts.
///
/// See the module documentation for conventions.
pub fn effective_size(g: &Graph, opts: StructuralHoleOptions) -> Result<Vec<f64>> {
    let ties = MutualTies::new(g, opts.weighted)?;
    Ok(per_node(&ties, opts.threads, |ties, acc, touched, i| {
        ties.effective_size(acc, touched, i)
    }))
}

/// Burt's efficiency (effective size divided by the number of contacts) of
/// every node, in node order. NaN for nodes without contacts.
///
/// See the module documentation for conventions.
pub fn efficiency(g: &Graph, opts: StructuralHoleOptions) -> Result<Vec<f64>> {
    let ties = MutualTies::new(g, opts.weighted)?;
    Ok(per_node(&ties, opts.threads, |ties, acc, touched, i| {
        let contacts = ties.contacts(i).len();
        ties.effective_size(acc, touched, i) / contacts as f64
    }))
}

/// Runs `f` for every node in parallel with a per-thread scatter buffer.
fn per_node<F>(ties: &MutualTies, threads: usize, f: F) -> Vec<f64>
where
    F: Fn(&MutualTies, &mut [f64], &mut Vec<u32>, usize) -> f64 + Sync,
{
    let n = ties.n();
    let parts = parallel::for_each_item(
        n,
        threads,
        || (vec![0.0; n], Vec::new(), Vec::new()),
        |(acc, touched, out): &mut (Vec<f64>, Vec<u32>, Vec<(u32, f64)>), i| {
            let v = if ties.contacts(i).is_empty() {
                f64::NAN
            } else {
                f(ties, acc, touched, i)
            };
            out.push((i as u32, v));
        },
    );
    let mut result = vec![0.0; n];
    for (_, _, out) in parts {
        for (i, v) in out {
            result[i as usize] = v;
        }
    }
    result
}

/// Mutual ties z_ij = a_ij + a_ji in CSR form, without self-loops, with
/// each node's row total and maximum.
struct MutualTies {
    offsets: Vec<usize>,
    targets: Vec<u32>,
    z: Vec<f64>,
    total: Vec<f64>,
    max: Vec<f64>,
}

impl MutualTies {
    fn new(g: &Graph, weighted: bool) -> Result<MutualTies> {
        if weighted {
            let w = g.edge_weights().ok_or(Error::MissingWeights)?;
            if let Some(edge) = w.iter().position(|&x| x <= 0.0) {
                return Err(Error::InvalidArgument(format!(
                    "edge {edge} has weight {}; structural-hole measures read weights \
                     as tie strengths and require them to be strictly positive",
                    w[edge]
                )));
            }
        }
        let n = g.node_count();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut targets = Vec::new();
        let mut z = Vec::new();
        offsets.push(0);
        let mut row: Vec<(u32, f64)> = Vec::new();
        for v in 0..n {
            row.clear();
            push_arcs(
                &mut row,
                g.out().neighbors(v),
                g.out().weights(v),
                v,
                weighted,
            );
            if g.is_directed() {
                push_arcs(
                    &mut row,
                    g.inc().neighbors(v),
                    g.inc().weights(v),
                    v,
                    weighted,
                );
                // Rows are each sorted; the merged row needs one more sort.
                row.sort_by_key(|&(t, _)| t);
            }
            // In an undirected graph z_ij = 2 a_ij; the factor cancels in
            // every ratio below, so a_ij is stored instead.
            for &(t, w) in row.iter() {
                if targets.len() > offsets[v] && *targets.last().unwrap() == t {
                    *z.last_mut().unwrap() += w;
                } else {
                    targets.push(t);
                    z.push(w);
                }
            }
            offsets.push(targets.len());
        }
        let mut total = vec![0.0; n];
        let mut max = vec![0.0; n];
        for v in 0..n {
            for &w in &z[offsets[v]..offsets[v + 1]] {
                total[v] += w;
                if w > max[v] {
                    max[v] = w;
                }
            }
        }
        Ok(MutualTies {
            offsets,
            targets,
            z,
            total,
            max,
        })
    }

    #[inline]
    fn n(&self) -> usize {
        self.total.len()
    }

    #[inline]
    fn contacts(&self, v: usize) -> &[u32] {
        &self.targets[self.offsets[v]..self.offsets[v + 1]]
    }

    #[inline]
    fn ties(&self, v: usize) -> &[f64] {
        &self.z[self.offsets[v]..self.offsets[v + 1]]
    }

    /// C_i, using `acc` (all zeros on entry and exit) as a scatter buffer.
    fn constraint(&self, acc: &mut [f64], touched: &mut Vec<u32>, i: usize) -> f64 {
        touched.clear();
        let ti = self.total[i];
        for (&q, &ziq) in self.contacts(i).iter().zip(self.ties(i)) {
            let piq = ziq / ti;
            acc[q as usize] += piq;
            touched.push(q);
            let tq = self.total[q as usize];
            for (&j, &zqj) in self.contacts(q as usize).iter().zip(self.ties(q as usize)) {
                if j as usize != i {
                    acc[j as usize] += piq * (zqj / tq);
                    touched.push(j);
                }
            }
        }
        let mut c = 0.0;
        for &j in self.contacts(i) {
            let x = acc[j as usize];
            c += x * x;
        }
        for &v in touched.iter() {
            acc[v as usize] = 0.0;
        }
        c
    }

    /// E_i, using `acc` (all zeros on entry and exit) as a scatter buffer.
    fn effective_size(&self, acc: &mut [f64], touched: &mut Vec<u32>, i: usize) -> f64 {
        touched.clear();
        let ti = self.total[i];
        for (&q, &ziq) in self.contacts(i).iter().zip(self.ties(i)) {
            let piq = ziq / ti;
            for (&j, &zqj) in self.contacts(q as usize).iter().zip(self.ties(q as usize)) {
                if j as usize != i {
                    // m_jq = z_jq / max_k z_jk, and z is symmetric.
                    acc[j as usize] += piq * (zqj / self.max[j as usize]);
                    touched.push(j);
                }
            }
        }
        let mut e = 0.0;
        for &j in self.contacts(i) {
            e += 1.0 - acc[j as usize];
        }
        for &v in touched.iter() {
            acc[v as usize] = 0.0;
        }
        e
    }
}

fn push_arcs(
    row: &mut Vec<(u32, f64)>,
    nbrs: &[u32],
    weights: Option<&[f64]>,
    v: usize,
    weighted: bool,
) {
    for (k, &t) in nbrs.iter().enumerate() {
        if t as usize == v {
            continue;
        }
        let w = match (weighted, weights) {
            (true, Some(w)) => w[k],
            _ => 1.0,
        };
        row.push((t, w));
    }
}

/// Names of the brokerage roles, in the column order of [`brokerage`].
pub const BROKERAGE_ROLES: [&str; 6] = [
    "coordinator",
    "itinerant",
    "representative",
    "gatekeeper",
    "liaison",
    "total",
];

/// Gould-Fernandez brokerage counts.
///
/// `groups` holds one group code per node; only equality of codes matters.
/// Returns an n-by-6 matrix whose columns are the roles in
/// [`BROKERAGE_ROLES`] order (coordinator, itinerant, representative,
/// gatekeeper, liaison, total), the column order of `sna::brokerage()`.
/// `threads == 0` uses all available cores.
///
/// See the module documentation for conventions.
pub fn brokerage(g: &Graph, groups: &[i64], threads: usize) -> Result<Dense> {
    let n = g.node_count();
    if groups.len() != n {
        return Err(Error::LengthMismatch {
            what: "groups",
            expected: n,
            found: groups.len(),
        });
    }
    let parts = parallel::for_each_item(
        n,
        threads,
        || (Vec::new(), Vec::new(), Vec::new()),
        |(ins, outs, rows): &mut BrokerageScratch, j| {
            distinct_without(ins, g.inc().neighbors(j), j);
            distinct_without(outs, g.out().neighbors(j), j);
            let mut counts = [0u64; 5];
            let gj = groups[j];
            for &i in ins.iter() {
                let gi = groups[i as usize];
                let from_i = g.out().neighbors(i as usize);
                for &k in outs.iter() {
                    if k == i || from_i.binary_search(&k).is_ok() {
                        continue;
                    }
                    let gk = groups[k as usize];
                    let role = if gi == gj {
                        if gj == gk {
                            0 // coordinator
                        } else {
                            2 // representative
                        }
                    } else if gj == gk {
                        3 // gatekeeper
                    } else if gi == gk {
                        1 // itinerant
                    } else {
                        4 // liaison
                    };
                    counts[role] += 1;
                }
            }
            rows.push((j as u32, counts));
        },
    );
    let mut out = Dense::zeros(n, 6);
    for (_, _, rows) in parts {
        for (j, counts) in rows {
            let mut total = 0u64;
            for (r, &c) in counts.iter().enumerate() {
                out.set(j as usize, r, c as f64);
                total += c;
            }
            out.set(j as usize, 5, total as f64);
        }
    }
    Ok(out)
}

/// Per-thread state for [`brokerage`]: a node's distinct in-neighbors and
/// out-neighbors, and the finished rows of counts.
type BrokerageScratch = (Vec<u32>, Vec<u32>, Vec<(u32, [u64; 5])>);

/// Copies a sorted neighbor list without repeats and without `v`.
fn distinct_without(dst: &mut Vec<u32>, src: &[u32], v: usize) {
    dst.clear();
    for &t in src {
        if t as usize != v && dst.last() != Some(&t) {
            dst.push(t);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let s: Vec<i64> = e.iter().map(|x| x.0).collect();
        let t: Vec<i64> = e.iter().map(|x| x.1).collect();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    fn close(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) {
            assert!(
                (x.is_nan() && y.is_nan()) || (x - y).abs() < 1e-12,
                "{a:?} vs {b:?}"
            );
        }
    }

    #[test]
    fn star_and_triangle_by_hand() {
        // Star with centre 0 and three leaves: the centre's contacts are
        // unconnected, so p = 1/3 each and C = 3 (1/3)^2 = 1/3. Each leaf
        // has one contact: C = 1.
        let g = und(4, &[(0, 1), (0, 2), (0, 3)]);
        let o = StructuralHoleOptions::default();
        close(&constraint(&g, o).unwrap(), &[1.0 / 3.0, 1.0, 1.0, 1.0]);
        close(&effective_size(&g, o).unwrap(), &[3.0, 1.0, 1.0, 1.0]);
        close(&efficiency(&g, o).unwrap(), &[1.0, 1.0, 1.0, 1.0]);

        // Triangle: p = 1/2, each term (1/2 + 1/4)^2, C = 2 * 9/16 = 9/8.
        // Effective size 2 - 2*1/2 = 1 (Borgatti 1997).
        let g = und(3, &[(0, 1), (1, 2), (0, 2)]);
        close(&constraint(&g, o).unwrap(), &[1.125; 3]);
        close(&effective_size(&g, o).unwrap(), &[1.0; 3]);
        close(&efficiency(&g, o).unwrap(), &[0.5; 3]);
    }

    #[test]
    fn isolates_and_loops() {
        let g = und(3, &[(0, 1), (2, 2)]);
        let o = StructuralHoleOptions::default();
        let c = constraint(&g, o).unwrap();
        assert_eq!(&c[..2], &[1.0, 1.0]);
        assert!(c[2].is_nan());
        assert!(effective_size(&g, o).unwrap()[2].is_nan());
        assert!(efficiency(&g, o).unwrap()[2].is_nan());
    }

    #[test]
    fn borgatti_simplification_holds_on_simple_graphs() {
        // Kite-like graph; compare with n - 2t/n computed directly.
        let e = [
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 3),
            (2, 3),
            (3, 4),
            (4, 5),
            (2, 5),
        ];
        let g = und(6, &e);
        let es = effective_size(&g, StructuralHoleOptions::default()).unwrap();
        let adj = |a: i64, b: i64| e.iter().any(|&(x, y)| (x, y) == (a, b) || (x, y) == (b, a));
        for v in 0..6i64 {
            let nb: Vec<i64> = (0..6).filter(|&u| u != v && adj(u, v)).collect();
            let mut t = 0;
            for (a, &x) in nb.iter().enumerate() {
                for &y in &nb[a + 1..] {
                    t += adj(x, y) as usize;
                }
            }
            let n = nb.len() as f64;
            assert!((es[v as usize] - (n - 2.0 * t as f64 / n)).abs() < 1e-12);
        }
    }

    #[test]
    fn directed_equals_reciprocated_and_weights_scale() {
        let s = [0, 1, 2, 3];
        let t = [1, 2, 0, 0];
        let d = Graph::from_edges(4, &s, &t, None, true).unwrap();
        let u = Graph::from_edges(4, &s, &t, None, false).unwrap();
        let o = StructuralHoleOptions::default();
        close(&constraint(&d, o).unwrap(), &constraint(&u, o).unwrap());
        // Multiplying every weight by a constant changes nothing.
        let w1 = Graph::from_edges(4, &s, &t, Some(&[1.0, 2.0, 3.0, 4.0]), true).unwrap();
        let w2 = Graph::from_edges(4, &s, &t, Some(&[2.5, 5.0, 7.5, 10.0]), true).unwrap();
        let ow = StructuralHoleOptions {
            weighted: true,
            ..Default::default()
        };
        close(&constraint(&w1, ow).unwrap(), &constraint(&w2, ow).unwrap());
        close(
            &effective_size(&w1, ow).unwrap(),
            &effective_size(&w2, ow).unwrap(),
        );
        let bad = Graph::from_edges(2, &[0], &[1], Some(&[0.0]), true).unwrap();
        assert!(constraint(&bad, ow).is_err());
    }

    #[test]
    fn brokerage_roles_by_hand() {
        // 0 -> 1 -> 2 with no 0 -> 2: node 1 brokers (0, 2).
        let g = Graph::from_edges(3, &[0, 1], &[1, 2], None, true).unwrap();
        let cases: [([i64; 3], usize); 5] = [
            ([0, 0, 0], 0), // coordinator
            ([0, 1, 0], 1), // itinerant
            ([0, 0, 1], 2), // representative
            ([0, 1, 1], 3), // gatekeeper
            ([0, 1, 2], 4), // liaison
        ];
        for (groups, role) in cases {
            let b = brokerage(&g, &groups, 1).unwrap();
            for r in 0..5 {
                assert_eq!(b.get(1, r), (r == role) as u8 as f64);
            }
            assert_eq!(b.get(1, 5), 1.0);
            assert_eq!(b.get(0, 5) + b.get(2, 5), 0.0);
        }
        // A closing arc removes the brokerage; parallel arcs and loops do
        // not change it.
        let g = Graph::from_edges(3, &[0, 1, 0], &[1, 2, 2], None, true).unwrap();
        assert_eq!(brokerage(&g, &[0, 0, 0], 1).unwrap().get(1, 5), 0.0);
        let g = Graph::from_edges(3, &[0, 0, 1, 1], &[1, 1, 2, 1], None, true).unwrap();
        assert_eq!(brokerage(&g, &[0, 0, 0], 1).unwrap().get(1, 5), 1.0);
        // Undirected: both directions count.
        let g = und(3, &[(0, 1), (1, 2)]);
        assert_eq!(brokerage(&g, &[0, 0, 0], 1).unwrap().get(1, 0), 2.0);
        assert!(brokerage(&g, &[0, 0], 1).is_err());
    }
}
