//! Closeness and harmonic centrality.
//!
//! Closeness is the reciprocal of a node's farness, the sum of its
//! shortest-path lengths to the other nodes (Sabidussi 1966,
//! `sabidussi1966centrality`; Freeman 1978, `freeman1978centrality`).
//! Harmonic centrality sums the reciprocals of those lengths instead, so an
//! unreachable node contributes 0 and the measure needs no special case for
//! disconnected networks (Marchiori and Latora 2000,
//! `marchiori2000harmony`; Rochat 2009, `rochat2009closeness`; Boldi and
//! Vigna 2014, `boldi2014axioms`).
//!
//! Every node is a source of one single-source search, run in parallel.
//! Without weights the searches run 64 at a time with the bit-parallel
//! breadth-first search of [`crate::bitset::MultiBfs`] (Then et al. 2014,
//! `then2014more`); the sum of reciprocal lengths is then formed level by
//! level, as (nodes at that distance) / distance.
//!
//! Conventions:
//! - Mode: `Out` uses lengths of paths leaving the node (igraph's default),
//!   `In` lengths of paths arriving at it (NetworkX's choice for directed
//!   graphs), and `All` ignores direction. Undirected graphs ignore the mode.
//! - Self-loops never lie on a shortest path. Parallel edges do not change
//!   path lengths; with weights the lightest of them counts.
//! - Without weights every edge has length 1 and the search is breadth-first
//!   search. With weights they are distances or strengths (lengths 1/w) and
//!   must be positive; the cheapest exact algorithm is chosen as described in
//!   [`crate::weights`]. When every edge has the same length L, breadth-first
//!   hop counts are multiplied by L, so closeness scales by 1/L and harmonic
//!   centrality by 1/L.
//! - Let n be the number of nodes, r the number of other nodes the node
//!   reaches, and S the sum of the lengths to them. Disconnected networks are
//!   handled by one of four rules ([`Unreachable`]):
//!
//!   | Rule | Raw | Normalized | r = 0 | Matches |
//!   |---|---|---|---|---|
//!   | `Ignore` | 1/S | r/S | NaN | igraph `closeness()` |
//!   | `WassermanFaust` | (normalized only) | (r/S)(r/(n-1)) | 0 | NetworkX `closeness_centrality(wf_improved = True)` |
//!   | `Zero` | 1/S, or 0 if r < n - 1 | (n-1)/S, or 0 if r < n - 1 | 0 | sna `closeness()` |
//!   | `MaxPlusOne` | 1/S' | (n-1)/S' | (n-1)/S' | UCINET Freeman closeness |
//!
//!   `Ignore` averages over reachable nodes only, so a node that reaches one
//!   neighbor scores as high as a hub. `WassermanFaust` scales that average
//!   by the share of nodes reached (Wasserman and Faust 1994, p. 201,
//!   `wasserman1994social`). In `MaxPlusOne`, S' replaces every missing
//!   length with one more than the longest finite length between any two
//!   distinct nodes; it is defined for unweighted networks only, since
//!   UCINET computes closeness on binary ties. All four agree on strongly
//!   connected networks. A single node has closeness NaN under every rule
//!   except `WassermanFaust`, where it is 0, as in NetworkX.
//! - Harmonic centrality is the sum of 1/d over the other nodes (0 for
//!   unreachable ones); normalized, it is divided by n - 1, as in igraph.
//!   sna (`closeness(cmode = "suminvdir")`) and UCINET ("reciprocal
//!   closeness") report the normalized form.

use std::borrow::Cow;

use crate::bitset::{for_each_bit, MultiBfs, Step, View, BATCH};
use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;
use crate::parallel;
use crate::weights::{path_lengths, IndexedHeap, Lengths, PathAlgorithm, WeightKind};

/// How closeness treats nodes that cannot be reached. See the module
/// documentation for the formulas.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Unreachable {
    /// Average over reachable nodes only (igraph).
    #[default]
    Ignore,
    /// Scale by the share of nodes reached (Wasserman and Faust; NetworkX).
    WassermanFaust,
    /// Closeness 0 unless the node reaches every other node (sna).
    Zero,
    /// Missing lengths become the longest observed length plus 1 (UCINET).
    MaxPlusOne,
}

impl Unreachable {
    /// Decodes the C interface codes: 0 = ignore, 1 = Wasserman-Faust,
    /// 2 = zero, 3 = max plus one.
    pub fn from_code(code: i32) -> Result<Unreachable> {
        match code {
            0 => Ok(Unreachable::Ignore),
            1 => Ok(Unreachable::WassermanFaust),
            2 => Ok(Unreachable::Zero),
            3 => Ok(Unreachable::MaxPlusOne),
            _ => Err(Error::InvalidArgument(format!(
                "unknown code {code} for unreachable nodes"
            ))),
        }
    }
}

/// Options for [`closeness`] and [`harmonic_centrality`].
#[derive(Debug, Clone, Copy)]
pub struct ClosenessOptions {
    /// Direction of the paths in directed graphs.
    pub mode: DegreeMode,
    /// Use edge weights (they must be strictly positive).
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Shortest-path algorithm; `Auto` picks the cheapest exact one.
    pub algorithm: PathAlgorithm,
    /// Normalize (see the module documentation).
    pub normalized: bool,
    /// Treatment of unreachable nodes (closeness only).
    pub unreachable: Unreachable,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for ClosenessOptions {
    fn default() -> Self {
        ClosenessOptions {
            mode: DegreeMode::Out,
            weighted: false,
            weights_as: WeightKind::Distance,
            algorithm: PathAlgorithm::Auto,
            normalized: false,
            unreachable: Unreachable::Ignore,
            threads: 0,
        }
    }
}

/// Closeness of every node, in node order.
///
/// See the module documentation for conventions.
pub fn closeness(g: &Graph, opts: ClosenessOptions) -> Result<Vec<f64>> {
    if opts.unreachable == Unreachable::WassermanFaust && !opts.normalized {
        return Err(Error::InvalidArgument(
            "the Wasserman-Faust closeness index is defined only in normalized form".into(),
        ));
    }
    if opts.unreachable == Unreachable::MaxPlusOne && opts.weighted {
        return Err(Error::InvalidArgument(
            "closeness with unreachable lengths set to the longest length plus 1 is \
             defined only for unweighted networks"
                .into(),
        ));
    }
    let n = g.node_count();
    let profiles = distance_profiles(g, &opts)?;
    let longest = profiles.iter().fold(0.0f64, |m, p| m.max(p.longest));
    let others = n.saturating_sub(1) as f64;

    let value = |p: &Profile| -> f64 {
        let r = p.reached as f64;
        match opts.unreachable {
            Unreachable::Ignore => {
                if p.reached == 0 {
                    f64::NAN
                } else if opts.normalized {
                    r / p.total
                } else {
                    1.0 / p.total
                }
            }
            Unreachable::WassermanFaust => {
                if p.reached == 0 {
                    0.0
                } else {
                    (r / p.total) * (r / others)
                }
            }
            Unreachable::Zero => {
                if n < 2 {
                    f64::NAN
                } else if p.reached < n - 1 {
                    0.0
                } else if opts.normalized {
                    others / p.total
                } else {
                    1.0 / p.total
                }
            }
            Unreachable::MaxPlusOne => {
                if n < 2 {
                    return f64::NAN;
                }
                let total = p.total + (others - r) * (longest + 1.0);
                if opts.normalized {
                    others / total
                } else {
                    1.0 / total
                }
            }
        }
    };
    Ok(profiles.iter().map(value).collect())
}

/// Harmonic centrality of every node, in node order. `unreachable` is
/// ignored: unreachable nodes always contribute 0.
///
/// See the module documentation for conventions.
pub fn harmonic_centrality(g: &Graph, opts: ClosenessOptions) -> Result<Vec<f64>> {
    let n = g.node_count();
    let profiles = distance_profiles(g, &opts)?;
    let scale = if opts.normalized && n > 1 {
        1.0 / (n - 1) as f64
    } else {
        1.0
    };
    Ok(profiles.iter().map(|p| p.harmonic * scale).collect())
}

/// What one single-source search records.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Profile {
    /// Other nodes reached.
    pub reached: usize,
    /// Sum of shortest-path lengths to them.
    pub total: f64,
    /// Sum of reciprocal lengths.
    pub harmonic: f64,
    /// Longest finite length.
    pub longest: f64,
}

/// The graph whose outgoing paths are the paths `mode` asks for.
pub(crate) fn oriented(g: &Graph, mode: DegreeMode) -> Result<Cow<'_, Graph>> {
    if !g.is_directed() || mode == DegreeMode::Out {
        return Ok(Cow::Borrowed(g));
    }
    let src: Vec<i64> = g.edge_sources().iter().map(|&v| v as i64).collect();
    let tgt: Vec<i64> = g.edge_targets().iter().map(|&v| v as i64).collect();
    let w = g.edge_weights();
    let h = match mode {
        DegreeMode::In => Graph::from_edges(g.node_count(), &tgt, &src, w, true)?,
        _ => Graph::from_edges(g.node_count(), &src, &tgt, w, false)?,
    };
    Ok(Cow::Owned(h))
}

/// One shortest-path search from every node, in parallel.
pub(crate) fn distance_profiles(g: &Graph, opts: &ClosenessOptions) -> Result<Vec<Profile>> {
    let mut clock = parallel::probe::Clock::start();
    // Unweighted searches read `g` in place through a view in the
    // requested mode; weighted ones need lengths aligned with the arcs of a
    // graph oriented for the mode, searched along its out-arcs.
    let (h, view_mode) = if opts.weighted {
        (oriented(g, opts.mode)?, DegreeMode::Out)
    } else {
        (Cow::Borrowed(g), opts.mode)
    };
    let h = h.as_ref();
    let view = View::new(h, view_mode);
    let n = h.node_count();
    let lengths = path_lengths(h, opts.weighted, opts.weights_as, opts.algorithm)?;
    clock.mark("setup");
    let partials = match &lengths {
        // Unweighted: breadth-first search from 64 sources at a time
        // (`crate::bitset::MultiBfs`); one batch is one work item, so the
        // result does not depend on the thread count.
        Lengths::Uniform => parallel::for_each_item(
            (n + BATCH - 1) / BATCH,
            opts.threads,
            || (MultiBfs::new(n), Vec::new()),
            |(ms, found), b| {
                let lo = b * BATCH;
                let hi = (lo + BATCH).min(n);
                let sources: Vec<u32> = (lo as u32..hi as u32).collect();
                let mut count = [0u32; BATCH];
                let mut total = [0u64; BATCH];
                let mut profiles = [Profile::default(); BATCH];
                ms.run(view, &sources, u32::MAX, |step| match step {
                    Step::Reached { sources, .. } => for_each_bit(sources, |i| count[i] += 1),
                    Step::LevelDone(level) => {
                        for i in 0..sources.len() {
                            let c = count[i];
                            if c == 0 {
                                continue;
                            }
                            count[i] = 0;
                            let p = &mut profiles[i];
                            p.reached += c as usize;
                            total[i] += c as u64 * level as u64;
                            p.harmonic += c as f64 / level as f64;
                            p.longest = level as f64;
                        }
                    }
                });
                for (i, &s) in sources.iter().enumerate() {
                    let mut p = profiles[i];
                    p.total = total[i] as f64;
                    found.push((s, p));
                }
            },
        )
        .into_iter()
        .map(|(_, found)| found)
        .collect::<Vec<_>>(),
        _ => parallel::for_each_item(
            n,
            opts.threads,
            || Search::new(n),
            |ws, s| {
                let p = match &lengths {
                    Lengths::Integer(lens, max) => ws.dial(h, lens, *max, s),
                    Lengths::Real(lens) => ws.dijkstra(h, lens, s),
                    Lengths::Uniform => unreachable!("handled above"),
                };
                ws.found.push((s as u32, p));
            },
        )
        .into_iter()
        .map(|ws| ws.found)
        .collect(),
    };
    // Uniform weights reduce to breadth-first search, which counts hops;
    // unlike betweenness, closeness depends on the common length itself.
    let unit = match (&lengths, opts.weighted) {
        (Lengths::Uniform, true) => {
            h.out()
                .all_weights()
                .and_then(|w| w.first())
                .map_or(1.0, |&w| match opts.weights_as {
                    WeightKind::Distance => w,
                    WeightKind::Strength => 1.0 / w,
                })
        }
        _ => 1.0,
    };
    clock.mark("search");
    let mut out = vec![Profile::default(); n];
    for found in partials {
        for (s, mut p) in found {
            if unit != 1.0 {
                p.total *= unit;
                p.harmonic /= unit;
                p.longest *= unit;
            }
            out[s as usize] = p;
        }
    }
    clock.mark("merge");
    Ok(out)
}

/// Per-thread scratch space, reused across sources.
struct Search {
    found: Vec<(u32, Profile)>,
    idist: Vec<u64>,
    dist: Vec<f64>,
    done: Vec<bool>,
    touched: Vec<u32>,
    buckets: Vec<Vec<u32>>,
    heap: IndexedHeap,
}

impl Search {
    fn new(n: usize) -> Search {
        Search {
            found: Vec::new(),
            idist: vec![u64::MAX; n],
            dist: vec![f64::INFINITY; n],
            done: vec![false; n],
            touched: Vec::new(),
            buckets: Vec::new(),
            heap: IndexedHeap::new(n),
        }
    }

    /// Dial's bucket queue with exact integer lengths (`dial1969algorithm`).
    fn dial(&mut self, g: &Graph, lens: &[u32], max: u32, s: usize) -> Profile {
        let adj = g.out();
        let width = max as usize + 1;
        if self.buckets.len() != width {
            self.buckets = vec![Vec::new(); width];
        }
        self.touched.clear();
        self.touched.push(s as u32);
        self.idist[s] = 0;
        self.buckets[0].push(s as u32);
        let mut pending = 1usize;
        let mut cur: u64 = 0;
        let (mut reached, mut total, mut harmonic, mut longest) = (0usize, 0u64, 0.0, 0u64);
        while pending > 0 {
            let slot = (cur % width as u64) as usize;
            let v = match self.buckets[slot].pop() {
                Some(v) => v as usize,
                None => {
                    cur += 1;
                    continue;
                }
            };
            pending -= 1;
            if self.done[v] || self.idist[v] != cur {
                continue;
            }
            self.done[v] = true;
            if v != s {
                reached += 1;
                total += cur;
                harmonic += 1.0 / cur as f64;
                longest = longest.max(cur);
            }
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                if self.done[w] {
                    continue;
                }
                let alt = cur + lens[k] as u64;
                if self.idist[w] == u64::MAX {
                    self.touched.push(w as u32);
                }
                if alt < self.idist[w] {
                    self.idist[w] = alt;
                    self.buckets[(alt % width as u64) as usize].push(w as u32);
                    pending += 1;
                }
            }
        }
        for &v in &self.touched {
            self.idist[v as usize] = u64::MAX;
            self.done[v as usize] = false;
        }
        Profile {
            reached,
            total: total as f64,
            harmonic,
            longest: longest as f64,
        }
    }

    /// Dijkstra's algorithm with an indexed 4-ary heap
    /// (`dijkstra1959note`).
    fn dijkstra(&mut self, g: &Graph, lens: &[f64], s: usize) -> Profile {
        let adj = g.out();
        self.touched.clear();
        self.touched.push(s as u32);
        self.dist[s] = 0.0;
        self.heap.push_or_decrease(s as u32, 0.0);
        let (mut reached, mut total, mut harmonic, mut longest) = (0usize, 0.0, 0.0, 0.0f64);
        while let Some((d, node)) = self.heap.pop() {
            let v = node as usize;
            self.done[v] = true;
            if v != s {
                reached += 1;
                total += d;
                harmonic += 1.0 / d;
                longest = longest.max(d);
            }
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                if self.done[w] {
                    continue;
                }
                let alt = d + lens[k];
                if self.dist[w].is_infinite() {
                    self.touched.push(w as u32);
                }
                if alt < self.dist[w] {
                    self.dist[w] = alt;
                    self.heap.push_or_decrease(w as u32, alt);
                }
            }
        }
        for &v in &self.touched {
            self.dist[v as usize] = f64::INFINITY;
            self.done[v as usize] = false;
        }
        Profile {
            reached,
            total,
            harmonic,
            longest,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, s: &[i64], t: &[i64], directed: bool) -> Graph {
        Graph::from_edges(n, s, t, None, directed).unwrap()
    }

    fn opts(unreachable: Unreachable, normalized: bool) -> ClosenessOptions {
        ClosenessOptions {
            unreachable,
            normalized,
            ..Default::default()
        }
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len()
            && a.iter().zip(b).all(|(x, y)| {
                (x.is_nan() && y.is_nan()) || (x - y).abs() <= 1e-12 * (1.0 + x.abs())
            })
    }

    #[test]
    fn path_graph() {
        // 0 - 1 - 2: farness 3, 2, 3.
        let g = graph(3, &[0, 1], &[1, 2], false);
        let raw = closeness(&g, opts(Unreachable::Ignore, false)).unwrap();
        assert!(close(&raw, &[1.0 / 3.0, 0.5, 1.0 / 3.0]));
        let norm = closeness(&g, opts(Unreachable::Ignore, true)).unwrap();
        assert!(close(&norm, &[2.0 / 3.0, 1.0, 2.0 / 3.0]));
        let h = harmonic_centrality(&g, opts(Unreachable::Ignore, false)).unwrap();
        assert!(close(&h, &[1.5, 2.0, 1.5]));
    }

    #[test]
    fn disconnected_rules() {
        // A path 0 - 1 - 2, an edge 3 - 4, and an isolate 5. n = 6.
        let g = graph(6, &[0, 1, 3], &[1, 2, 4], false);
        let ig = closeness(&g, opts(Unreachable::Ignore, true)).unwrap();
        assert!(close(&ig, &[2.0 / 3.0, 1.0, 2.0 / 3.0, 1.0, 1.0, f64::NAN]));
        let wf = closeness(&g, opts(Unreachable::WassermanFaust, true)).unwrap();
        let a = (2.0 / 3.0) * (2.0 / 5.0);
        assert!(close(&wf, &[a, 0.4, a, 0.2, 0.2, 0.0]));
        let zero = closeness(&g, opts(Unreachable::Zero, true)).unwrap();
        assert!(close(&zero, &[0.0; 6]));
        // Longest finite length is 2, so missing lengths become 3.
        let ucinet = closeness(&g, opts(Unreachable::MaxPlusOne, true)).unwrap();
        assert!(close(
            &ucinet,
            &[
                5.0 / 12.0,
                5.0 / 11.0,
                5.0 / 12.0,
                5.0 / 13.0,
                5.0 / 13.0,
                5.0 / 15.0
            ]
        ));
    }

    #[test]
    fn reproduces_published_ucinet_output() {
        // UCINET 6.849's Closeness log for its g9_disc example (a star on
        // nodes 0-4 with hub 0, a path 5-6-7-8, and a triangle 9-10-11), as
        // published with `borgatti2026xucinet` (inst/goldens/centrality/
        // log_menu.txt), printed to three decimals: FreeClo and RecipClo.
        let g = graph(
            12,
            &[0, 0, 0, 0, 5, 6, 7, 9, 10, 9],
            &[1, 2, 3, 4, 6, 7, 8, 10, 11, 11],
            false,
        );
        let free = [
            0.344, 0.314, 0.314, 0.314, 0.314, 0.289, 0.306, 0.306, 0.289, 0.289, 0.289, 0.289,
        ];
        let recip = [
            0.364, 0.227, 0.227, 0.227, 0.227, 0.167, 0.227, 0.227, 0.167, 0.182, 0.182, 0.182,
        ];
        let c = closeness(&g, opts(Unreachable::MaxPlusOne, true)).unwrap();
        let h = harmonic_centrality(&g, opts(Unreachable::Ignore, true)).unwrap();
        for i in 0..12 {
            assert!((c[i] - free[i]).abs() < 5e-4, "FreeClo {i}: {}", c[i]);
            assert!((h[i] - recip[i]).abs() < 5e-4, "RecipClo {i}: {}", h[i]);
        }
    }

    #[test]
    fn directed_modes() {
        // 0 -> 1 -> 2.
        let g = graph(3, &[0, 1], &[1, 2], true);
        let mut o = opts(Unreachable::Ignore, false);
        assert!(close(
            &closeness(&g, o).unwrap(),
            &[1.0 / 3.0, 1.0, f64::NAN]
        ));
        o.mode = DegreeMode::In;
        assert!(close(
            &closeness(&g, o).unwrap(),
            &[f64::NAN, 1.0, 1.0 / 3.0]
        ));
        o.mode = DegreeMode::All;
        assert!(close(
            &closeness(&g, o).unwrap(),
            &[1.0 / 3.0, 0.5, 1.0 / 3.0]
        ));
    }

    #[test]
    fn weighted_algorithms_agree() {
        let g = Graph::from_edges(
            5,
            &[0, 0, 1, 2, 3],
            &[1, 2, 3, 3, 4],
            Some(&[2.0, 1.0, 1.0, 3.0, 2.0]),
            true,
        )
        .unwrap();
        let run = |algorithm| {
            let o = ClosenessOptions {
                weighted: true,
                algorithm,
                ..Default::default()
            };
            (
                closeness(&g, o).unwrap(),
                harmonic_centrality(&g, o).unwrap(),
            )
        };
        let (c1, h1) = run(PathAlgorithm::Dial);
        let (c2, h2) = run(PathAlgorithm::Dijkstra);
        assert!(close(&c1, &c2) && close(&h1, &h2));
        // Node 0: lengths 2, 1, 3, 5.
        assert!((c1[0] - 1.0 / 11.0).abs() < 1e-15);
    }

    #[test]
    fn uniform_weights_scale_lengths() {
        // Every edge has weight 2.5: distances are 2.5 times the hop counts.
        let g = Graph::from_edges(3, &[0, 1], &[1, 2], Some(&[2.5, 2.5]), false).unwrap();
        let mut o = ClosenessOptions {
            weighted: true,
            ..Default::default()
        };
        assert!(close(
            &closeness(&g, o).unwrap(),
            &[1.0 / 7.5, 0.2, 1.0 / 7.5]
        ));
        assert!(close(
            &harmonic_centrality(&g, o).unwrap(),
            &[0.6, 0.8, 0.6]
        ));
        o.weights_as = WeightKind::Strength;
        assert!(close(
            &closeness(&g, o).unwrap(),
            &[1.0 / 1.2, 1.25, 1.0 / 1.2]
        ));
    }

    #[test]
    fn single_node_and_option_checks() {
        let g = graph(1, &[], &[], false);
        assert!(closeness(&g, opts(Unreachable::Ignore, true)).unwrap()[0].is_nan());
        assert_eq!(
            closeness(&g, opts(Unreachable::WassermanFaust, true)).unwrap(),
            vec![0.0]
        );
        assert_eq!(
            harmonic_centrality(&g, opts(Unreachable::Ignore, true)).unwrap(),
            vec![0.0]
        );
        assert!(closeness(&g, opts(Unreachable::WassermanFaust, false)).is_err());
        let w = Graph::from_edges(2, &[0], &[1], Some(&[1.0]), false).unwrap();
        let o = ClosenessOptions {
            weighted: true,
            unreachable: Unreachable::MaxPlusOne,
            ..Default::default()
        };
        assert!(closeness(&w, o).is_err());
    }
}
