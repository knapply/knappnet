//! Goodness-of-fit statistics (Hunter, Goodreau, and Handcock 2008,
//! `hunter2008goodness`).
//!
//! A fitted model is checked by simulating networks from it and comparing
//! distributions of structural features that were not necessarily in the
//! model with those of the observed network:
//!
//! - the degree distribution (in- and out-degree for directed networks):
//!   the number of nodes with each degree 0, ..., n - 1;
//! - edgewise shared partners: the number of ties whose endpoints have
//!   exactly k common partners, k = 0, ..., n - 2 (outgoing two-path
//!   partners i -> k -> j for a directed tie, as in [`super::terms`]);
//! - dyadwise shared partners: the same over every dyad (ordered pairs in
//!   directed networks);
//! - geodesic distances: the number of pairs at each shortest-path length
//!   1, ..., n - 1, then the number of unreachable pairs (unordered pairs
//!   for undirected networks, ordered otherwise), which is `ergm`'s
//!   `distance` statistic;
//! - the triad census for directed networks, in the order of
//!   [`crate::metrics::triad::TRIAD_NAMES`].
//!
//! These are the statistics of `ergm`'s `gof()` defaults. Each is computed
//! from a tie list by the same code for observed and simulated networks.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::triad::triad_census;

use super::state::NetState;

/// A goodness-of-fit statistic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GofStat {
    /// Degree distribution (undirected).
    Degree,
    /// In-degree distribution (directed).
    InDegree,
    /// Out-degree distribution (directed).
    OutDegree,
    /// Edgewise shared partners.
    Esp,
    /// Dyadwise shared partners.
    Dsp,
    /// Geodesic distance distribution.
    Distance,
    /// Triad census (directed).
    TriadCensus,
}

impl GofStat {
    /// Decodes the C interface codes: 0 = degree, 1 = in-degree, 2 =
    /// out-degree, 3 = esp, 4 = dsp, 5 = distance, 6 = triad census.
    pub fn from_code(code: i32) -> Result<GofStat> {
        Ok(match code {
            0 => GofStat::Degree,
            1 => GofStat::InDegree,
            2 => GofStat::OutDegree,
            3 => GofStat::Esp,
            4 => GofStat::Dsp,
            5 => GofStat::Distance,
            6 => GofStat::TriadCensus,
            _ => {
                return Err(Error::InvalidArgument(format!(
                    "unknown goodness-of-fit statistic code {code}"
                )))
            }
        })
    }

    /// Number of values the statistic has for `n` nodes.
    pub fn len(&self, n: usize) -> usize {
        match self {
            GofStat::Degree | GofStat::InDegree | GofStat::OutDegree => n,
            GofStat::Esp | GofStat::Dsp => n.saturating_sub(1),
            GofStat::Distance => n,
            GofStat::TriadCensus => 16,
        }
    }
}

/// Computes `which` for the network with the given ties.
pub fn gof_stat(n: usize, directed: bool, ties: &[(u32, u32)], which: GofStat) -> Result<Vec<f64>> {
    let len = which.len(n);
    let mut out = vec![0.0; len];
    match which {
        GofStat::Degree | GofStat::InDegree | GofStat::OutDegree => {
            if which == GofStat::Degree && directed {
                return Err(Error::InvalidArgument(
                    "use the in- or out-degree distribution for a directed network".into(),
                ));
            }
            if which != GofStat::Degree && !directed {
                return Err(Error::InvalidArgument(
                    "in- and out-degree distributions are for directed networks".into(),
                ));
            }
            let mut deg = vec![0usize; n];
            for &(a, b) in ties {
                match which {
                    GofStat::Degree => {
                        deg[a as usize] += 1;
                        deg[b as usize] += 1;
                    }
                    GofStat::OutDegree => deg[a as usize] += 1,
                    _ => deg[b as usize] += 1,
                }
            }
            for d in deg {
                if d < len {
                    out[d] += 1.0;
                }
            }
        }
        GofStat::Esp | GofStat::Dsp => {
            let mut y = NetState::empty(n, directed);
            for &(a, b) in ties {
                if !y.has(a as usize, b as usize) {
                    y.toggle(a as usize, b as usize);
                }
            }
            let sp = |a: usize, b: usize| {
                if directed {
                    y.two_paths(a, b)
                } else {
                    y.shared_partners_undirected(a, b)
                }
            };
            if which == GofStat::Esp {
                for &(a, b) in y.ties() {
                    let k = sp(a as usize, b as usize);
                    if k < len {
                        out[k] += 1.0;
                    }
                }
            } else {
                for (a, b) in y.dyads() {
                    let k = sp(a as usize, b as usize);
                    if k < len {
                        out[k] += 1.0;
                    }
                }
            }
        }
        GofStat::Distance => {
            // Breadth-first search from every node; the last cell counts
            // unreachable pairs.
            let mut adj = vec![Vec::new(); n];
            for &(a, b) in ties {
                adj[a as usize].push(b as usize);
                if !directed {
                    adj[b as usize].push(a as usize);
                }
            }
            let mut dist = vec![usize::MAX; n];
            let mut queue = Vec::with_capacity(n);
            for s in 0..n {
                dist.iter_mut().for_each(|d| *d = usize::MAX);
                dist[s] = 0;
                queue.clear();
                queue.push(s);
                let mut head = 0;
                while head < queue.len() {
                    let v = queue[head];
                    head += 1;
                    for &w in &adj[v] {
                        if dist[w] == usize::MAX {
                            dist[w] = dist[v] + 1;
                            queue.push(w);
                        }
                    }
                }
                for (t, &dt) in dist.iter().enumerate() {
                    if t == s || (!directed && t < s) {
                        continue;
                    }
                    match dt {
                        usize::MAX => out[len - 1] += 1.0,
                        d if d >= 1 && d - 1 < len - 1 => out[d - 1] += 1.0,
                        _ => {}
                    }
                }
            }
        }
        GofStat::TriadCensus => {
            if !directed {
                return Err(Error::InvalidArgument(
                    "the triad census goodness-of-fit statistic is for directed networks".into(),
                ));
            }
            let src: Vec<i64> = ties.iter().map(|t| t.0 as i64).collect();
            let tgt: Vec<i64> = ties.iter().map(|t| t.1 as i64).collect();
            let g = Graph::from_edges(n, &src, &tgt, None, true)?;
            out.copy_from_slice(&triad_census(&g, 1));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn statistics_of_a_small_network() {
        // Path 0-1-2-3 plus the tie 0-2: one triangle.
        let ties = [(0, 1), (1, 2), (2, 3), (0, 2)];
        let deg = gof_stat(4, false, &ties, GofStat::Degree).unwrap();
        assert_eq!(deg, vec![0.0, 1.0, 2.0, 1.0]);
        let esp = gof_stat(4, false, &ties, GofStat::Esp).unwrap();
        // Ties 0-1, 1-2, 0-2 have one shared partner; 2-3 has none.
        assert_eq!(esp, vec![1.0, 3.0, 0.0]);
        let dsp = gof_stat(4, false, &ties, GofStat::Dsp).unwrap();
        // Dyads: 0-1 (1), 0-2 (1), 0-3 (1: via 2), 1-2 (1), 1-3 (1), 2-3 (0).
        assert_eq!(dsp, vec![1.0, 5.0, 0.0]);
        let dist = gof_stat(4, false, &ties, GofStat::Distance).unwrap();
        // Distances: four pairs at 1, two at 2 (0-3 and 1-3), none further.
        assert_eq!(dist, vec![4.0, 2.0, 0.0, 0.0]);
    }

    #[test]
    fn directed_statistics() {
        let ties = [(0, 1), (1, 2), (0, 2), (2, 0)];
        let indeg = gof_stat(3, true, &ties, GofStat::InDegree).unwrap();
        assert_eq!(indeg, vec![0.0, 2.0, 1.0]);
        let esp = gof_stat(3, true, &ties, GofStat::Esp).unwrap();
        // 0->2 has the two-path 0->1->2; the others have none.
        assert_eq!(esp, vec![3.0, 1.0]);
        let tc = gof_stat(3, true, &ties, GofStat::TriadCensus).unwrap();
        assert_eq!(tc.iter().sum::<f64>(), 1.0);
        let dist = gof_stat(3, true, &ties, GofStat::Distance).unwrap();
        // Ordered pairs: 0->1, 1->2, 0->2, 2->0 at 1; 1->0 at 2, 2->1 at 2.
        assert_eq!(dist, vec![4.0, 2.0, 0.0]);
        assert!(gof_stat(3, true, &ties, GofStat::Degree).is_err());
    }
}
