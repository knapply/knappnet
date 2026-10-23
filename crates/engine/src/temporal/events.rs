//! Inter-event times, burstiness, and memory.
//!
//! The events of a node are the start times of its contacts (sent,
//! received, or both; a self-loop is one event), and the events of a node
//! pair are the start times of the contacts between them (ordered pairs in
//! directed networks, unordered otherwise, grouped as in
//! [`super::windows::aggregate`]). Sorted, consecutive events give the
//! inter-event times τ_1, ..., τ_k; simultaneous events give τ = 0.
//!
//! Burstiness and memory follow Goh and Barabási (2008,
//! `goh2008burstiness`, Eqs. 1 and 4):
//!
//! - B = (σ - μ) / (σ + μ), where μ and σ are the mean and standard
//!   deviation of the inter-event times. B is -1 for perfectly regular
//!   events, about 0 for a Poisson process, and approaches 1 for very
//!   bursty ones. σ is the standard deviation of the observed
//!   distribution (divisor k).
//! - M = (1 / (k - 1)) sum_{i=1}^{k-1} (τ_i - μ_1)(τ_{i+1} - μ_2) / (σ_1 σ_2),
//!   where μ_1, σ_1 (μ_2, σ_2) are the mean and standard deviation of
//!   τ_1..τ_{k-1} (τ_2..τ_k). With divisor k - 1 for σ_1 and σ_2, as here,
//!   M is the Pearson correlation of consecutive inter-event times.
//!
//! B is undefined (NaN) with fewer than two inter-event times or when all
//! are zero; M with fewer than three, or when either set it correlates is
//! constant. Neither measure corrects for the number of events (Kim and Jo
//! 2016 give a finite-size correction for B); compare nodes with similar
//! activity, or few events, with care. Weights play no part.

use super::windows::aggregate;
use super::TemporalNetwork;
use crate::error::Result;
use crate::metrics::degree::DegreeMode;

/// Whose events are counted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventScope {
    /// Events of each node, with a [`DegreeMode`] for directed networks.
    Node(DegreeMode),
    /// Events of each node pair.
    Pair,
}

/// Inter-event times grouped by node or pair.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InterEvents {
    /// Number of groups: nodes, or pairs (as in [`super::windows::aggregate`]).
    pub groups: usize,
    /// Group of each inter-event time, non-decreasing.
    pub group: Vec<u32>,
    /// Inter-event times, in time order within each group.
    pub gap: Vec<f64>,
    /// Number of events in each group.
    pub events: Vec<f64>,
}

/// Inter-event times of every node or pair.
pub fn inter_event_times(net: &TemporalNetwork, scope: EventScope) -> Result<InterEvents> {
    let m = net.contact_count();
    let mut ev: Vec<(u32, f64)> = Vec::new();
    let groups = match scope {
        EventScope::Node(mode) => {
            let (out, inc) = match (net.is_directed(), mode) {
                (false, _) | (true, DegreeMode::All) => (true, true),
                (true, DegreeMode::Out) => (true, false),
                (true, DegreeMode::In) => (false, true),
            };
            ev.reserve(2 * m);
            for c in 0..m {
                let (u, v, t) = (net.src[c], net.tgt[c], net.start[c]);
                if out || (u == v && inc) {
                    ev.push((u, t));
                }
                if inc && u != v {
                    ev.push((v, t));
                }
            }
            net.node_count()
        }
        EventScope::Pair => {
            let agg = aggregate(net, None)?;
            ev.reserve(m);
            for c in 0..m {
                ev.push((agg.pair_of[c] as u32, net.start[c]));
            }
            agg.from.len()
        }
    };
    ev.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let mut out = InterEvents {
        groups,
        events: vec![0.0; groups],
        ..Default::default()
    };
    for (i, &(g, t)) in ev.iter().enumerate() {
        out.events[g as usize] += 1.0;
        if i > 0 && ev[i - 1].0 == g {
            out.group.push(g);
            out.gap.push(t - ev[i - 1].1);
        }
    }
    Ok(out)
}

/// Burstiness B and memory M of one sequence of inter-event times.
pub fn burstiness_memory(gaps: &[f64]) -> (f64, f64) {
    let k = gaps.len();
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let sd = |x: &[f64], mu: f64| {
        (x.iter().map(|v| (v - mu) * (v - mu)).sum::<f64>() / x.len() as f64).sqrt()
    };
    let b = if k < 2 {
        f64::NAN
    } else {
        let mu = mean(gaps);
        let sigma = sd(gaps, mu);
        if sigma + mu > 0.0 {
            (sigma - mu) / (sigma + mu)
        } else {
            f64::NAN
        }
    };
    let m = if k < 3 {
        f64::NAN
    } else {
        let (x, y) = (&gaps[..k - 1], &gaps[1..]);
        let (m1, m2) = (mean(x), mean(y));
        let (s1, s2) = (sd(x, m1), sd(y, m2));
        if s1 > 0.0 && s2 > 0.0 {
            let cov: f64 = x
                .iter()
                .zip(y)
                .map(|(a, b)| (a - m1) * (b - m2))
                .sum::<f64>()
                / (k - 1) as f64;
            (cov / (s1 * s2)).clamp(-1.0, 1.0)
        } else {
            f64::NAN
        }
    };
    (b, m)
}

/// Event counts, burstiness, and memory of every node or pair.
#[derive(Debug, Clone, PartialEq)]
pub struct Burstiness {
    /// Number of events in each group.
    pub events: Vec<f64>,
    /// Goh and Barabási's B.
    pub burstiness: Vec<f64>,
    /// Goh and Barabási's M.
    pub memory: Vec<f64>,
}

/// Burstiness and memory of every node or pair.
pub fn burstiness(net: &TemporalNetwork, scope: EventScope) -> Result<Burstiness> {
    let ie = inter_event_times(net, scope)?;
    let mut out = Burstiness {
        events: ie.events.clone(),
        burstiness: vec![f64::NAN; ie.groups],
        memory: vec![f64::NAN; ie.groups],
    };
    let mut i = 0;
    while i < ie.gap.len() {
        let g = ie.group[i];
        let mut j = i;
        while j < ie.gap.len() && ie.group[j] == g {
            j += 1;
        }
        let (b, m) = burstiness_memory(&ie.gap[i..j]);
        out.burstiness[g as usize] = b;
        out.memory[g as usize] = m;
        i = j;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_and_bursty_sequences() {
        let (b, m) = burstiness_memory(&[2.0, 2.0, 2.0]);
        assert_eq!(b, -1.0);
        assert!(m.is_nan());
        // Gaps 1, 3: mean 2, sd 1, B = -1/3.
        let (b, _) = burstiness_memory(&[1.0, 3.0]);
        assert!((b + 1.0 / 3.0).abs() < 1e-15);
        // Alternating gaps are perfectly anti-correlated.
        let (_, m) = burstiness_memory(&[1.0, 5.0, 1.0, 5.0]);
        assert!((m + 1.0).abs() < 1e-12);
        assert!(burstiness_memory(&[0.0, 0.0]).0.is_nan());
    }

    #[test]
    fn events_by_node_and_pair() {
        // 0 -> 1 at 0, 1 -> 0 at 4, 0 -> 2 at 1, 2 -> 2 at 3.
        let g = TemporalNetwork::new(
            3,
            &[0, 1, 0, 2],
            &[1, 0, 2, 2],
            &[0.0, 4.0, 1.0, 3.0],
            None,
            true,
        )
        .unwrap();
        let ie = inter_event_times(&g, EventScope::Node(DegreeMode::All)).unwrap();
        assert_eq!(ie.events, vec![3.0, 2.0, 2.0]);
        assert_eq!(ie.group, vec![0, 0, 1, 2]);
        assert_eq!(ie.gap, vec![1.0, 3.0, 4.0, 2.0]);
        let out = inter_event_times(&g, EventScope::Node(DegreeMode::Out)).unwrap();
        assert_eq!(out.events, vec![2.0, 1.0, 1.0]);
        let pairs = inter_event_times(&g, EventScope::Pair).unwrap();
        assert_eq!(pairs.groups, 4);
        assert!(pairs.gap.is_empty());
    }
}
