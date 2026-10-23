//! Slicing a temporal network into snapshots, and aggregating contacts
//! into a weighted static network.
//!
//! Activity follows the half-open spells of networkDynamic
//! (`butts2023networkdynamic`): a contact [s, e) with e > s is active at
//! time t when s <= t < e, and an instantaneous contact (e = s) only at
//! t = s. A contact is active in a window [a, b), a < b, when it is active
//! at some time in the window: s < b and e > a for an interval contact,
//! a <= s < b for an instantaneous one. Back-to-back contacts [1, 2) and
//! [2, 3) therefore never overlap, consecutive windows [0, 10) and
//! [10, 20) of a tumbling series never share an instantaneous contact, and
//! an interval contact belongs to every window it overlaps.
//!
//! A series of windows whose starts and ends are both non-decreasing (as
//! in sliding and tumbling series) is computed in one sweep over the
//! contacts sorted by start time, with a heap of active contacts keyed by
//! end time, in O(m log m + total output); other series cost one scan of
//! the contacts that start before each window ends.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::TemporalNetwork;
use crate::error::{Error, Result};

/// Whether a contact [s, e] is active at time t.
#[inline]
pub fn active_at(s: f64, e: f64, t: f64) -> bool {
    s <= t && (t < e || s == t)
}

/// Whether a contact [s, e] is active in the window [a, b).
#[inline]
pub fn active_in(s: f64, e: f64, a: f64, b: f64) -> bool {
    s < b && (e > a || s >= a)
}

/// Indices (in input order) of the contacts active at time `t`.
pub fn contacts_at(net: &TemporalNetwork, t: f64) -> Result<Vec<usize>> {
    if t.is_nan() {
        return Err(Error::InvalidArgument(
            "the time must not be missing".into(),
        ));
    }
    Ok((0..net.contact_count())
        .filter(|&c| active_at(net.start[c], net.end[c], t))
        .collect())
}

/// Contacts active in each window: window k holds
/// `indices[offsets[k]..offsets[k + 1]]`, in input order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct WindowContacts {
    /// Start of each window's indices; one more entry than windows.
    pub offsets: Vec<usize>,
    /// Contact indices.
    pub indices: Vec<usize>,
}

/// Contacts active in each window [from[k], to[k]).
pub fn contacts_in_windows(
    net: &TemporalNetwork,
    from: &[f64],
    to: &[f64],
) -> Result<WindowContacts> {
    if from.len() != to.len() {
        return Err(Error::LengthMismatch {
            what: "window ends",
            expected: from.len(),
            found: to.len(),
        });
    }
    for k in 0..from.len() {
        if from[k].is_nan() || to[k].is_nan() || from[k] >= to[k] {
            return Err(Error::InvalidArgument(format!(
                "window {k} must start before it ends"
            )));
        }
    }
    let m = net.contact_count();
    let mut order: Vec<u32> = (0..m as u32).collect();
    order.sort_unstable_by(|&x, &y| {
        net.start[x as usize]
            .total_cmp(&net.start[y as usize])
            .then(x.cmp(&y))
    });
    let monotone = (1..from.len()).all(|k| from[k - 1] <= from[k] && to[k - 1] <= to[k]);
    let mut out = WindowContacts {
        offsets: Vec::with_capacity(from.len() + 1),
        indices: Vec::new(),
    };
    out.offsets.push(0);
    if !monotone {
        for k in 0..from.len() {
            let stop = order.partition_point(|&c| net.start[c as usize] < to[k]);
            let first = out.indices.len();
            out.indices.extend(
                order[..stop]
                    .iter()
                    .map(|&c| c as usize)
                    .filter(|&c| active_in(net.start[c], net.end[c], from[k], to[k])),
            );
            out.indices[first..].sort_unstable();
            out.offsets.push(out.indices.len());
        }
        return Ok(out);
    }
    // Active contacts are kept in a set with O(1) removal, and in a heap
    // ordered by when they stop being active: an interval contact leaves
    // once the window start reaches its end, an instantaneous one once
    // the start passes it, so at equal end times intervals leave first.
    let mut active: Vec<u32> = Vec::new();
    let mut slot = vec![u32::MAX; m];
    let mut heap: BinaryHeap<Reverse<(Key, u32)>> = BinaryHeap::new();
    let mut next = 0;
    for k in 0..from.len() {
        let (a, b) = (from[k], to[k]);
        while next < m && net.start[order[next] as usize] < b {
            let c = order[next];
            let (s, e) = (net.start[c as usize], net.end[c as usize]);
            slot[c as usize] = active.len() as u32;
            active.push(c);
            heap.push(Reverse((Key(e, s == e), c)));
            next += 1;
        }
        while let Some(&Reverse((_, c))) = heap.peek() {
            let (s, e) = (net.start[c as usize], net.end[c as usize]);
            if active_in(s, e, a, b) {
                break;
            }
            heap.pop();
            let i = slot[c as usize] as usize;
            let last = *active.last().expect("non-empty");
            active.swap_remove(i);
            if last != c {
                slot[last as usize] = i as u32;
            }
        }
        let first = out.indices.len();
        out.indices.extend(active.iter().map(|&c| c as usize));
        out.indices[first..].sort_unstable();
        out.offsets.push(out.indices.len());
    }
    Ok(out)
}

/// Heap key: end time, then instantaneous contacts after intervals.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(f64, bool);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1))
    }
}

/// Windows of width `width` starting at `start`, `start + step`, ..., up
/// to the last one that starts no later than `end` (sliding when
/// `step < width`, tumbling when equal, sampling with gaps when larger).
/// With `end` the last time in the data, every contact falls in some
/// window when `step <= width`.
pub fn regular_windows(
    start: f64,
    end: f64,
    width: f64,
    step: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if !(start.is_finite() && end.is_finite()) || start > end {
        return Err(Error::InvalidArgument(
            "the series must have a finite start no later than its end".into(),
        ));
    }
    if !(width.is_finite() && width > 0.0 && step.is_finite() && step > 0.0) {
        return Err(Error::InvalidArgument(
            "window width and step must be positive and finite".into(),
        ));
    }
    let count = ((end - start) / step).floor() as u64 + 1;
    if count > 100_000_000 {
        return Err(Error::TooLarge {
            what: "number of windows",
            value: count as usize,
        });
    }
    let mut a = Vec::with_capacity(count as usize);
    let mut b = Vec::with_capacity(count as usize);
    for k in 0..count {
        // Multiplying rather than accumulating keeps window starts exact
        // for integer and dyadic steps.
        let s = start + k as f64 * step;
        if s > end {
            break;
        }
        a.push(s);
        b.push(s + width);
    }
    Ok((a, b))
}

/// Contacts grouped by node pair.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Aggregate {
    /// Source of each pair (for undirected networks, as in the pair's
    /// first contact).
    pub from: Vec<u32>,
    /// Target of each pair.
    pub to: Vec<u32>,
    /// Number of contacts.
    pub count: Vec<f64>,
    /// Total duration, clipped to the window when one is given.
    pub duration: Vec<f64>,
    /// Earliest start.
    pub first: Vec<f64>,
    /// Latest end.
    pub last: Vec<f64>,
    /// For each contact considered, the index of its pair (`usize::MAX`
    /// for contacts outside the window).
    pub pair_of: Vec<usize>,
}

/// Groups contacts by node pair (ordered pairs when directed, unordered
/// otherwise), in order of each pair's first contact. With `window`
/// [a, b), only contacts active in it count and durations are clipped to
/// it, so an instantaneous contact contributes a count but no duration.
pub fn aggregate(net: &TemporalNetwork, window: Option<(f64, f64)>) -> Result<Aggregate> {
    if let Some((a, b)) = window {
        if a.is_nan() || b.is_nan() || a >= b {
            return Err(Error::InvalidArgument(
                "the window must start before it ends".into(),
            ));
        }
    }
    let m = net.contact_count();
    let n = net.node_count() as u64;
    let key = |c: usize| -> u64 {
        let (u, v) = (net.src[c] as u64, net.tgt[c] as u64);
        let (u, v) = if net.directed || u <= v {
            (u, v)
        } else {
            (v, u)
        };
        u * n + v
    };
    let keep = |c: usize| match window {
        None => true,
        Some((a, b)) => active_in(net.start[c], net.end[c], a, b),
    };
    let mut order: Vec<(u64, u32)> = (0..m)
        .filter(|&c| keep(c))
        .map(|c| (key(c), c as u32))
        .collect();
    order.sort_unstable();
    let mut agg = Aggregate {
        pair_of: vec![usize::MAX; m],
        ..Default::default()
    };
    // Groups in key order, then renumbered by first contact.
    let mut groups: Vec<(u32, usize, usize)> = Vec::new();
    let mut i = 0;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && order[j].0 == order[i].0 {
            j += 1;
        }
        groups.push((order[i].1, i, j));
        i = j;
    }
    groups.sort_unstable_by_key(|g| g.0);
    for (gid, &(first_contact, i, j)) in groups.iter().enumerate() {
        let c0 = first_contact as usize;
        agg.from.push(net.src[c0]);
        agg.to.push(net.tgt[c0]);
        let (mut count, mut dur) = (0.0, 0.0);
        let (mut first, mut last) = (f64::INFINITY, f64::NEG_INFINITY);
        for &(_, c) in &order[i..j] {
            let c = c as usize;
            let (s, e) = (net.start[c], net.end[c]);
            count += 1.0;
            dur += match window {
                None => e - s,
                Some((a, b)) => (e.min(b) - s.max(a)).max(0.0),
            };
            first = first.min(s);
            last = last.max(e);
            agg.pair_of[c] = gid;
        }
        agg.count.push(count);
        agg.duration.push(dur);
        agg.first.push(first);
        agg.last.push(last);
    }
    Ok(agg)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn net() -> TemporalNetwork {
        // Contacts: [0,2), [2,2], [2,5), [5,5], [1,1].
        TemporalNetwork::new(
            3,
            &[0, 1, 0, 2, 1],
            &[1, 0, 2, 1, 0],
            &[0.0, 2.0, 2.0, 5.0, 1.0],
            Some(&[2.0, 2.0, 5.0, 5.0, 1.0]),
            false,
        )
        .unwrap()
    }

    #[test]
    fn half_open_activity() {
        let g = net();
        assert_eq!(contacts_at(&g, 2.0).unwrap(), vec![1, 2]);
        assert_eq!(contacts_at(&g, 1.0).unwrap(), vec![0, 4]);
        assert_eq!(contacts_at(&g, 5.0).unwrap(), vec![3]);
    }

    #[test]
    fn sweep_matches_direct_scan() {
        let g = net();
        let (a, b) = regular_windows(0.0, 6.0, 2.0, 1.0).unwrap();
        assert_eq!(a, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let w = contacts_in_windows(&g, &a, &b).unwrap();
        for k in 0..a.len() {
            let want: Vec<usize> = (0..5)
                .filter(|&c| active_in(g.start[c], g.end[c], a[k], b[k]))
                .collect();
            assert_eq!(&w.indices[w.offsets[k]..w.offsets[k + 1]], &want[..]);
        }
        assert_eq!(&w.indices[w.offsets[0]..w.offsets[1]], &[0, 4]);
        assert_eq!(&w.indices[w.offsets[2]..w.offsets[3]], &[1, 2]);
        // Unsorted windows take the direct path.
        let w2 = contacts_in_windows(&g, &[2.0, 0.0], &[4.0, 2.0]).unwrap();
        assert_eq!(w2.indices, vec![1, 2, 0, 4]);
        assert!(contacts_in_windows(&g, &[1.0], &[1.0]).is_err());
    }

    #[test]
    fn aggregation_groups_pairs() {
        let g = net();
        let agg = aggregate(&g, None).unwrap();
        assert_eq!(agg.from, vec![0, 0, 2]);
        assert_eq!(agg.to, vec![1, 2, 1]);
        assert_eq!(agg.count, vec![3.0, 1.0, 1.0]);
        assert_eq!(agg.duration, vec![2.0, 3.0, 0.0]);
        assert_eq!(agg.first, vec![0.0, 2.0, 5.0]);
        assert_eq!(agg.pair_of, vec![0, 0, 1, 2, 0]);
        let w = aggregate(&g, Some((1.0, 3.0))).unwrap();
        assert_eq!(w.count, vec![3.0, 1.0]);
        assert_eq!(w.duration, vec![1.0, 1.0]);
    }
}
