//! Temporal networks: contacts with times, time-respecting paths, slicing
//! into windows, and event-sequence statistics.
//!
//! A temporal network is a set of nodes `0..n` and a list of *contacts*
//! (events), each with endpoints u and v, a start time s, and an end time
//! e >= s. A contact with e = s is instantaneous (an e-mail, a login, a
//! packet); one with e > s lasts over an interval (a phone call, a
//! session). This is the contact-sequence and interval-graph representation
//! of Holme and Saramäki (2012, `holme2012temporal`). Times are plain
//! `f64` values in whatever unit the caller uses; the R package passes
//! seconds for dates and date-times.
//!
//! # When a contact is active
//!
//! Intervals are half-open, as in networkDynamic's spells: a contact with
//! e > s is active at time t when s <= t < e, so back-to-back contacts
//! [1, 2) and [2, 3) never overlap. An instantaneous contact (e = s) is
//! active only at t = s. A contact is active in a window [a, b) (a < b)
//! when s < b and e > a, or, if instantaneous, when a <= s < b. See
//! [`windows`].
//!
//! # Time-respecting paths
//!
//! Paths follow the model of Wu et al. (2014, `wu2014path`): every contact
//! is a temporal arc (u, v, t, λ) that departs u at its start time t and
//! arrives at v at t + λ, where the traversal time λ is a constant (zero by
//! default) or the contact's own duration e - s (as for a flight). A
//! time-respecting walk is a sequence of arcs, each departing from the node
//! where the previous one arrived, such that the waiting time
//! w = t_{i+1} - (t_i + λ_i) at every intermediate node satisfies
//!
//! - w >= `min_wait` (non-decreasing times when `min_wait` = 0, the
//!   definition of Kempe, Kleinberg, and Kumar 2002, `kempe2002connectivity`,
//!   and of Holme and Saramäki 2012), or w > `min_wait` with `strict` (with
//!   `min_wait` = 0, strictly increasing times, as in Pan and Saramäki 2011,
//!   `pan2011path`); and
//! - w <= `max_wait` (no limit by default).
//!
//! Waiting constraints apply between consecutive contacts only: the walk
//! may leave its first node at any time in the window. Within a window
//! [t_start, t_end], the first contact departs at or after t_start and the
//! last arrives at or before t_end (the closed interval of Wu et al.).
//!
//! Four optimal walks from a root (Wu et al. 2014, Definition 1):
//! earliest arrival, latest departure (to the root), fastest (smallest
//! arrival minus departure), and shortest, here the fewest contacts (Wu et
//! al. minimize total traversal time; with a cost of 1 per contact the
//! same pair-list technique counts hops). Without `max_wait`, an optimal
//! walk can always be shortened to a path that repeats no node. With
//! `max_wait`, waiting can be capped only by moving on, so an optimal walk
//! may need to revisit a node; results are defined over walks.
//!
//! Conventions:
//! - Undirected contacts can be traversed in both directions; directed
//!   ones from u to v only.
//! - Self-loop contacts never move anything and are ignored by every path
//!   computation (in particular they do not reset the waiting clock).
//! - Parallel contacts (same endpoints, any times) are separate contacts.
//! - Node presence intervals are not consulted by path computations: a
//!   contact between two nodes is usable whenever the contact itself is.
//! - Weights play no part in temporal paths.
//! - Interval contacts are traversed from their start time (with λ as
//!   above); walks that enter a long contact part-way through (interval
//!   graphs) are not modelled. Split long contacts into events at the
//!   resolution needed.
//!
//! # Algorithms
//!
//! Arcs are sorted once by departure time (and by tail within equal
//! times) into an *edge stream*, and every computation is one pass over it
//! per root, following Wu et al. (2014, Algorithms 1, 2, 4, and the
//! shortest-path variant): latest departure is the earliest arrival of the
//! time-reversed stream, in which arc (u, v, t, λ) becomes (v, u, -(t + λ),
//! λ). Each node keeps the arrivals that can still be extended, labelled
//! with the best value of the walk that produced them (latest start for
//! fastest walks, fewest hops for shortest walks); an arc departing at t
//! extends the best arrival whose waiting time is admissible. Arrivals
//! become usable once their waiting time reaches `min_wait` (a heap keyed
//! by arrival time) and expire after `max_wait` (a monotone queue keeps the
//! best unexpired label, the sliding-window technique), so waiting
//! constraints cost O(log m) per arrival. With zero traversal time and
//! non-strict waiting, contacts with the same time can chain; such a group
//! is closed by re-scanning the arcs of every node first reached within the
//! group. Roots run in parallel. Total work is O(m log m) per root, and
//! O(m) per root for earliest arrival without `max_wait`.
//!
//! Temporal closeness is Pan and Saramäki's (2011, Eq. 2)
//! C_i = (1 / (n - 1)) sum_j 1 / τ_ij, where τ_ij is the latency from i to
//! j averaged over all starting times in the window with their pair-wise
//! periodic boundary condition (Eq. 1). It is computed exactly from the
//! fastest-path profile, the Pareto set of (departure, arrival) pairs.
//! Reachability, the set of influence and source set, and the reachability
//! ratio follow Holme (2005, `holme2005network`) and Holme and Saramäki
//! (2012). Burstiness and memory coefficients follow Goh and Barabási
//! (2008, `goh2008burstiness`); see [`events`].

pub mod events;
pub mod paths;
pub mod windows;

use std::sync::OnceLock;

use crate::error::{Error, Result};

pub use paths::{
    temporal_closeness, temporal_paths, temporal_reach_counts, Criterion, Direction, PathOptions,
    TemporalCloseness, Traversal,
};

/// A temporal network: contacts with start and end times.
#[derive(Debug)]
pub struct TemporalNetwork {
    n: usize,
    directed: bool,
    src: Vec<u32>,
    tgt: Vec<u32>,
    start: Vec<f64>,
    end: Vec<f64>,
    /// Arcs used by path computations: one per directed non-loop contact,
    /// two per undirected non-loop contact.
    arc_contact: Vec<u32>,
    arc_from: Vec<u32>,
    arc_to: Vec<u32>,
    /// Arc orders for the edge streams, built on first use.
    by_start: OnceLock<Vec<u32>>,
    by_start_desc: OnceLock<Vec<u32>>,
    by_end_desc: OnceLock<Vec<u32>>,
}

impl TemporalNetwork {
    /// Builds a temporal network from contacts with 0-based endpoints.
    ///
    /// `end` is `None` for instantaneous contacts; otherwise it has one
    /// value per contact, at least the contact's start. Times must be
    /// finite.
    pub fn new(
        n: usize,
        sources: &[i64],
        targets: &[i64],
        start: &[f64],
        end: Option<&[f64]>,
        directed: bool,
    ) -> Result<TemporalNetwork> {
        if n > u32::MAX as usize {
            return Err(Error::TooLarge {
                what: "node count",
                value: n,
            });
        }
        let m = sources.len();
        for (what, len) in [("targets", targets.len()), ("start times", start.len())] {
            if len != m {
                return Err(Error::LengthMismatch {
                    what,
                    expected: m,
                    found: len,
                });
            }
        }
        if let Some(e) = end {
            if e.len() != m {
                return Err(Error::LengthMismatch {
                    what: "end times",
                    expected: m,
                    found: e.len(),
                });
            }
        }
        if m > (u32::MAX as usize) / 2 {
            return Err(Error::TooLarge {
                what: "contact count",
                value: m,
            });
        }
        let check = |edge: usize, node: i64| -> Result<u32> {
            if node < 0 || node as u64 >= n as u64 {
                Err(Error::NodeOutOfRange { edge, node, n })
            } else {
                Ok(node as u32)
            }
        };
        let mut src = Vec::with_capacity(m);
        let mut tgt = Vec::with_capacity(m);
        for e in 0..m {
            src.push(check(e, sources[e])?);
            tgt.push(check(e, targets[e])?);
            if !start[e].is_finite() {
                return Err(Error::InvalidArgument(format!(
                    "edge {e} has a missing or non-finite start time"
                )));
            }
        }
        let end = match end {
            None => start.to_vec(),
            Some(e) => {
                for c in 0..m {
                    if !e[c].is_finite() {
                        return Err(Error::InvalidArgument(format!(
                            "edge {c} has a missing or non-finite end time"
                        )));
                    }
                    if e[c] < start[c] {
                        return Err(Error::InvalidArgument(format!(
                            "edge {c} ends before it starts"
                        )));
                    }
                }
                e.to_vec()
            }
        };

        let mut arc_contact = Vec::with_capacity(if directed { m } else { 2 * m });
        let mut arc_from = Vec::with_capacity(arc_contact.capacity());
        let mut arc_to = Vec::with_capacity(arc_contact.capacity());
        for c in 0..m {
            let (u, v) = (src[c], tgt[c]);
            if u == v {
                continue;
            }
            arc_contact.push(c as u32);
            arc_from.push(u);
            arc_to.push(v);
            if !directed {
                arc_contact.push(c as u32);
                arc_from.push(v);
                arc_to.push(u);
            }
        }

        Ok(TemporalNetwork {
            n,
            directed,
            src,
            tgt,
            start: start.to_vec(),
            end,
            arc_contact,
            arc_from,
            arc_to,
            by_start: OnceLock::new(),
            by_start_desc: OnceLock::new(),
            by_end_desc: OnceLock::new(),
        })
    }

    /// Number of nodes.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.n
    }

    /// Number of contacts as supplied (self-loops included).
    #[inline]
    pub fn contact_count(&self) -> usize {
        self.src.len()
    }

    /// Whether contacts are directed.
    #[inline]
    pub fn is_directed(&self) -> bool {
        self.directed
    }

    /// Contact sources in input order.
    pub fn sources(&self) -> &[u32] {
        &self.src
    }

    /// Contact targets in input order.
    pub fn targets(&self) -> &[u32] {
        &self.tgt
    }

    /// Start times in input order.
    pub fn start_times(&self) -> &[f64] {
        &self.start
    }

    /// End times in input order (equal to the start for instantaneous
    /// contacts).
    pub fn end_times(&self) -> &[f64] {
        &self.end
    }

    /// Earliest start time, or `None` without contacts.
    pub fn first_time(&self) -> Option<f64> {
        self.start.iter().copied().reduce(f64::min)
    }

    /// Latest end time, or `None` without contacts.
    pub fn last_time(&self) -> Option<f64> {
        self.end.iter().copied().reduce(f64::max)
    }

    fn order_by_start(&self) -> &[u32] {
        self.by_start.get_or_init(|| {
            let mut idx: Vec<u32> = (0..self.arc_from.len() as u32).collect();
            idx.sort_unstable_by(|&a, &b| {
                let (a, b) = (a as usize, b as usize);
                let (ca, cb) = (self.arc_contact[a] as usize, self.arc_contact[b] as usize);
                self.start[ca]
                    .total_cmp(&self.start[cb])
                    .then(self.arc_from[a].cmp(&self.arc_from[b]))
                    .then(a.cmp(&b))
            });
            idx
        })
    }

    /// Arcs in the order of the time-reversed stream: decreasing `key`
    /// (start or end time), then increasing head (the tail after reversal).
    fn order_desc(&self, by_end: bool) -> &[u32] {
        let cell = if by_end {
            &self.by_end_desc
        } else {
            &self.by_start_desc
        };
        cell.get_or_init(|| {
            let key = if by_end { &self.end } else { &self.start };
            let mut idx: Vec<u32> = (0..self.arc_from.len() as u32).collect();
            idx.sort_unstable_by(|&a, &b| {
                let (a, b) = (a as usize, b as usize);
                let (ca, cb) = (self.arc_contact[a] as usize, self.arc_contact[b] as usize);
                key[cb]
                    .total_cmp(&key[ca])
                    .then(self.arc_to[a].cmp(&self.arc_to[b]))
                    .then(a.cmp(&b))
            });
            idx
        })
    }

    /// The edge stream for one direction and traversal time: arcs as
    /// (tail, head, departure, arrival), sorted by departure and then tail.
    pub(crate) fn stream(&self, direction: Direction, traversal: Traversal) -> Stream {
        let lambda = |c: usize| match traversal {
            Traversal::Constant(x) => x,
            Traversal::Duration => self.end[c] - self.start[c],
        };
        let order: &[u32] = match (direction, traversal) {
            (Direction::Out, _) => self.order_by_start(),
            (Direction::In, Traversal::Constant(_)) => self.order_desc(false),
            (Direction::In, Traversal::Duration) => self.order_desc(true),
        };
        let k = order.len();
        let mut s = Stream {
            from: Vec::with_capacity(k),
            to: Vec::with_capacity(k),
            dep: Vec::with_capacity(k),
            arr: Vec::with_capacity(k),
        };
        for &a in order {
            let a = a as usize;
            let c = self.arc_contact[a] as usize;
            let (t, arrive) = (self.start[c], self.start[c] + lambda(c));
            match direction {
                Direction::Out => {
                    s.from.push(self.arc_from[a]);
                    s.to.push(self.arc_to[a]);
                    s.dep.push(t);
                    s.arr.push(arrive);
                }
                Direction::In => {
                    s.from.push(self.arc_to[a]);
                    s.to.push(self.arc_from[a]);
                    s.dep.push(-arrive);
                    s.arr.push(-t);
                }
            }
        }
        // Adding a constant traversal time can round distinct start times
        // to the same departure, which would break the tail order within
        // equal departures that the sweep relies on. Re-sort if so.
        let sorted = (1..k).all(|i| {
            s.dep[i - 1] < s.dep[i] || (s.dep[i - 1] == s.dep[i] && s.from[i - 1] <= s.from[i])
        });
        if !sorted {
            s.sort();
        }
        s
    }
}

/// Arcs sorted by departure time, then tail.
#[derive(Debug, Clone, Default)]
pub(crate) struct Stream {
    pub(crate) from: Vec<u32>,
    pub(crate) to: Vec<u32>,
    pub(crate) dep: Vec<f64>,
    pub(crate) arr: Vec<f64>,
}

impl Stream {
    fn sort(&mut self) {
        let mut idx: Vec<usize> = (0..self.from.len()).collect();
        idx.sort_by(|&a, &b| {
            self.dep[a]
                .total_cmp(&self.dep[b])
                .then(self.from[a].cmp(&self.from[b]))
        });
        let from = idx.iter().map(|&i| self.from[i]).collect();
        let to = idx.iter().map(|&i| self.to[i]).collect();
        let dep = idx.iter().map(|&i| self.dep[i]).collect();
        let arr = idx.iter().map(|&i| self.arr[i]).collect();
        *self = Stream { from, to, dep, arr };
    }

    pub(crate) fn len(&self) -> usize {
        self.from.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_contacts() {
        assert!(matches!(
            TemporalNetwork::new(2, &[0], &[2], &[0.0], None, true),
            Err(Error::NodeOutOfRange { .. })
        ));
        assert!(TemporalNetwork::new(2, &[0], &[1], &[f64::NAN], None, true).is_err());
        let e = TemporalNetwork::new(2, &[0], &[1], &[2.0], Some(&[1.0]), true).unwrap_err();
        assert_eq!(e.to_string(), "edge 0 ends before it starts");
        assert!(TemporalNetwork::new(2, &[0], &[1, 1], &[0.0], None, true).is_err());
    }

    #[test]
    fn streams_are_sorted_and_reversed() {
        let net = TemporalNetwork::new(
            3,
            &[0, 1, 2, 1],
            &[1, 2, 2, 0],
            &[3.0, 1.0, 0.0, 1.0],
            None,
            true,
        )
        .unwrap();
        let s = net.stream(Direction::Out, Traversal::Constant(0.0));
        // The loop 2 -> 2 is dropped.
        assert_eq!(s.dep, vec![1.0, 1.0, 3.0]);
        assert_eq!(s.from, vec![1, 1, 0]);
        let r = net.stream(Direction::In, Traversal::Constant(1.0));
        assert_eq!(r.dep, vec![-4.0, -2.0, -2.0]);
        assert_eq!(r.from, vec![1, 0, 2]);
        assert_eq!(r.to, vec![0, 1, 1]);
        assert_eq!(r.arr, vec![-3.0, -1.0, -1.0]);
    }

    #[test]
    fn undirected_contacts_give_two_arcs() {
        let net = TemporalNetwork::new(2, &[0], &[1], &[5.0], Some(&[7.0]), false).unwrap();
        let s = net.stream(Direction::Out, Traversal::Duration);
        assert_eq!(s.len(), 2);
        assert_eq!(s.arr, vec![7.0, 7.0]);
        assert_eq!(net.first_time(), Some(5.0));
        assert_eq!(net.last_time(), Some(7.0));
    }
}
