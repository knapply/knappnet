//! Maximal k-plexes.
//!
//! A k-plex is a set S of nodes in which every member is adjacent to at
//! least |S| - k of the other members, that is, has at most k - 1
//! non-neighbors among them (Seidman and Foster 1978, `seidman1978graph`).
//! 1-plexes are cliques. A k-plex is maximal when no
//! other node can join it.
//!
//! Algorithm: being a k-plex is hereditary (every subset of a k-plex is
//! one), so maximal k-plexes can be listed with the set-enumeration scheme
//! of Bron and Kerbosch (1973, `bron1973algorithm`) without pivoting: a
//! growing k-plex R, candidates P that could join it, and excluded nodes
//! X that could join it but whose sets were already listed. R is reported
//! when P and X are empty. Heredity guarantees that a node unable to join
//! R can never join a superset of R, which makes the scheme exact (each
//! maximal k-plex is reported exactly once). Two exact reductions make it
//! practical:
//! - Every member of a k-plex with at least s nodes has at least s - k
//!   neighbors in it, so nodes outside the (s - k)-core of the network can
//!   neither belong to nor extend a k-plex of at least s = `min_size`
//!   nodes, and are removed first (Seidman and Foster 1978; Seidman 1983,
//!   `seidman1983network`).
//! - A k-plex with at least 2k - 1 nodes has diameter at most 2 (Seidman
//!   and Foster 1978): two non-adjacent members each have at least |S| - k
//!   neighbors among the other |S| - 2 members, so at least |S| - 2k + 2
//!   in common. So when `min_size` >= 2k - 1, the
//!   candidates and excluded nodes of a search started at node v are
//!   limited to nodes within two steps of v.
//!
//! Conventions:
//! - The undirected simple reading of the network is used (see
//!   [`crate::metrics::factions`]).
//! - `min_size` drops smaller k-plexes. Any set of at most k nodes is a
//!   k-plex, so small values report many trivial sets; values below
//!   2k - 1 also allow disconnected k-plexes and disable the second
//!   reduction.
//! - The number of maximal k-plexes can grow exponentially with the size
//!   of the network. The search counts its steps (one per set R examined)
//!   and stops with an error after `max_steps`, rather than running for an
//!   unbounded time. Whether that limit is reached does not depend on the
//!   number of threads.
//! - Searches started at different nodes run in parallel. Members are
//!   listed in increasing order and k-plexes in lexicographic order.

use std::sync::atomic::{AtomicU64, Ordering};

use super::{adjacent, canonical_sets, simple_rows};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;

/// Options for [`k_plexes`].
#[derive(Debug, Clone, Copy)]
pub struct KPlexOptions {
    /// Each member may have at most k - 1 non-neighbors among the other
    /// members (k >= 1).
    pub k: usize,
    /// Smallest k-plex reported.
    pub min_size: usize,
    /// Largest number of search steps before giving up with an error.
    pub max_steps: u64,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

struct Search<'a> {
    rows: &'a [Vec<u32>],
    k: usize,
    min_size: usize,
    steps: &'a AtomicU64,
    max_steps: u64,
}

impl Search<'_> {
    /// Whether `r` plus `v` is a k-plex, given `missing[i]`, the number of
    /// non-neighbors of r[i] among the other members of r.
    fn fits(&self, r: &[u32], missing: &[usize], v: u32) -> bool {
        let mut own = 0;
        for (i, &u) in r.iter().enumerate() {
            if !adjacent(self.rows, u, v) {
                own += 1;
                if own >= self.k || missing[i] + 1 >= self.k {
                    return false;
                }
            }
        }
        true
    }

    fn run(
        &self,
        r: &mut Vec<u32>,
        missing: &mut Vec<usize>,
        p: Vec<u32>,
        x: Vec<u32>,
        out: &mut Vec<Vec<u32>>,
    ) -> Result<()> {
        if self.steps.fetch_add(1, Ordering::Relaxed) >= self.max_steps {
            return Err(Error::InvalidArgument(format!(
                "the search for k-plexes exceeded {} steps; raise the step limit, \
                 increase the minimum size, or use a smaller k",
                self.max_steps
            )));
        }
        if p.is_empty() {
            if x.is_empty() && r.len() >= self.min_size {
                out.push(r.clone());
            }
            return Ok(());
        }
        if r.len() + p.len() < self.min_size {
            return Ok(());
        }
        let mut x = x;
        for (i, &v) in p.iter().enumerate() {
            // Add v to R.
            let mut own = 0;
            for (j, &u) in r.iter().enumerate() {
                if !adjacent(self.rows, u, v) {
                    missing[j] += 1;
                    own += 1;
                }
            }
            r.push(v);
            missing.push(own);
            let p2: Vec<u32> = p[i + 1..]
                .iter()
                .copied()
                .filter(|&w| self.fits(r, missing, w))
                .collect();
            let x2: Vec<u32> = x
                .iter()
                .copied()
                .filter(|&w| self.fits(r, missing, w))
                .collect();
            let res = self.run(r, missing, p2, x2, out);
            // Remove v again.
            r.pop();
            missing.pop();
            for (j, &u) in r.iter().enumerate() {
                if !adjacent(self.rows, u, v) {
                    missing[j] -= 1;
                }
            }
            res?;
            x.push(v);
        }
        Ok(())
    }
}

/// Removes nodes until every remaining node has at least `need` remaining
/// neighbors; returns which nodes remain.
fn core_filter(rows: &[Vec<u32>], need: usize) -> Vec<bool> {
    let n = rows.len();
    let mut alive = vec![true; n];
    let mut deg: Vec<usize> = rows.iter().map(|r| r.len()).collect();
    let mut stack: Vec<usize> = (0..n).filter(|&v| deg[v] < need).collect();
    for &v in &stack {
        alive[v] = false;
    }
    while let Some(v) = stack.pop() {
        for &w in &rows[v] {
            let w = w as usize;
            if alive[w] {
                deg[w] -= 1;
                if deg[w] < need {
                    alive[w] = false;
                    stack.push(w);
                }
            }
        }
    }
    alive
}

/// Maximal k-plexes of `g` with at least `min_size` members. See the
/// module documentation.
pub fn k_plexes(g: &Graph, opts: KPlexOptions) -> Result<Vec<Vec<u32>>> {
    let k = opts.k;
    if k == 0 {
        return Err(Error::InvalidArgument("k must be at least 1".into()));
    }
    let full = simple_rows(g);
    let n = full.len();
    let alive = core_filter(&full, opts.min_size.saturating_sub(k));
    let rows: Vec<Vec<u32>> = full
        .iter()
        .enumerate()
        .map(|(v, r)| {
            if alive[v] {
                r.iter().copied().filter(|&w| alive[w as usize]).collect()
            } else {
                Vec::new()
            }
        })
        .collect();
    let local = opts.min_size >= 2 * k - 1;
    let steps = AtomicU64::new(0);
    let search = Search {
        rows: &rows,
        k,
        min_size: opts.min_size,
        steps: &steps,
        max_steps: opts.max_steps,
    };
    let parts = parallel::for_each_item(
        n,
        opts.threads,
        || (Vec::new(), None::<Error>),
        |(out, err): &mut (Vec<Vec<u32>>, Option<Error>), v| {
            if err.is_some() || !alive[v] || steps.load(Ordering::Relaxed) >= opts.max_steps {
                if alive[v] && err.is_none() {
                    *err = Some(Error::InvalidArgument(format!(
                        "the search for k-plexes exceeded {} steps; raise the step limit, \
                         increase the minimum size, or use a smaller k",
                        opts.max_steps
                    )));
                }
                return;
            }
            // Candidates: later nodes; excluded: earlier nodes; both
            // limited to two steps from v when that is exact.
            let pool: Vec<u32> = if local {
                let mut near: Vec<u32> = rows[v].clone();
                for &w in &rows[v] {
                    near.extend_from_slice(&rows[w as usize]);
                }
                near.sort_unstable();
                near.dedup();
                near.retain(|&w| w as usize != v);
                near
            } else {
                (0..n as u32)
                    .filter(|&w| w as usize != v && alive[w as usize])
                    .collect()
            };
            let mut r = vec![v as u32];
            let mut missing = vec![0usize];
            let fits = |w: u32| search.fits(&r, &missing, w);
            let p: Vec<u32> = pool
                .iter()
                .copied()
                .filter(|&w| w as usize > v && fits(w))
                .collect();
            let x: Vec<u32> = pool
                .iter()
                .copied()
                .filter(|&w| (w as usize) < v && fits(w))
                .collect();
            if let Err(e) = search.run(&mut r, &mut missing, p, x, out) {
                *err = Some(e);
            }
        },
    );
    let mut sets = Vec::new();
    for (found, err) in parts {
        if let Some(e) = err {
            return Err(e);
        }
        sets.extend(found);
    }
    Ok(canonical_sets(sets))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::cliques::max_cliques;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    fn opts(k: usize, min_size: usize) -> KPlexOptions {
        KPlexOptions {
            k,
            min_size,
            max_steps: 10_000_000,
            threads: 2,
        }
    }

    /// Every subset of the nodes, checked directly.
    fn brute(g: &Graph, k: usize, min_size: usize) -> Vec<Vec<u32>> {
        let rows = simple_rows(g);
        let n = rows.len();
        let is_plex = |mask: u32| {
            let members: Vec<u32> = (0..n as u32).filter(|&v| mask >> v & 1 == 1).collect();
            members.iter().all(|&a| {
                members
                    .iter()
                    .filter(|&&b| b != a && !adjacent(&rows, a, b))
                    .count()
                    < k
            })
        };
        let plexes: Vec<u32> = (1..(1u32 << n)).filter(|&m| is_plex(m)).collect();
        let mut out = Vec::new();
        for &m in &plexes {
            let maximal = (0..n).all(|v| m >> v & 1 == 1 || !is_plex(m | 1 << v));
            if maximal && m.count_ones() as usize >= min_size {
                out.push((0..n as u32).filter(|&v| m >> v & 1 == 1).collect());
            }
        }
        canonical_sets(out)
    }

    #[test]
    fn one_plexes_are_cliques() {
        let g = und(6, &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2), (4, 5)]);
        assert_eq!(k_plexes(&g, opts(1, 1)).unwrap(), max_cliques(&g, 1, 1));
    }

    #[test]
    fn agrees_with_brute_force() {
        let g = und(
            9,
            &[
                (0, 1),
                (0, 2),
                (1, 2),
                (1, 3),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (4, 6),
                (6, 7),
                (7, 8),
                (2, 7),
            ],
        );
        for k in 1..=3 {
            for min_size in 0..=6 {
                assert_eq!(
                    k_plexes(&g, opts(k, min_size)).unwrap(),
                    brute(&g, k, min_size),
                    "k = {k}, min_size = {min_size}"
                );
            }
        }
    }

    #[test]
    fn step_limit_is_an_error() {
        let g = und(8, &[]);
        let mut o = opts(3, 0);
        o.max_steps = 10;
        assert!(k_plexes(&g, o).is_err());
        assert!(k_plexes(&g, opts(0, 1)).is_err());
    }
}
