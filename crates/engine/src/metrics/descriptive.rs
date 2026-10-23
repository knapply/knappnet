//! Descriptive and structural measures: reciprocity, transitivity and
//! clustering, shortest-path distances and the measures built on them
//! (mean distance, diameter, eccentricity, radius), reachability, cut
//! points and bridges, and the counts a network summary reports.
//!
//! Sources:
//! - Reciprocity. The share of arcs that are reciprocated (Newman, Forrest,
//!   and Balthrop 2002, `newman2002email`); the ratio of mutual to non-null
//!   dyads and the share of symmetric dyads (Wasserman and Faust 1994,
//!   `wasserman1994social`); and the correlation between the two arcs of a
//!   dyad, rho = (r - p) / (1 - p) with r the share of reciprocated arcs and
//!   p the density (Garlaschelli and Loffredo 2004,
//!   `garlaschelli2004patterns`).
//! - Transitivity. The global coefficient is three times the number of
//!   triangles divided by the number of connected triples (Newman 2003,
//!   `newman2003structure`; Wasserman and Faust 1994). The local clustering
//!   coefficient of a node is the share of pairs of its neighbors that are
//!   adjacent (Watts and Strogatz 1998, `watts1998collective`); how nodes
//!   with fewer than two neighbors enter an average is discussed by Kaiser
//!   (2008, `kaiser2008mean`). Weighted local coefficients follow Barrat et
//!   al. (2004, `barrat2004architecture`) and Onnela et al. (2005,
//!   `onnela2005intensity`). Triangles are listed once each by orienting
//!   every edge from the endpoint of lower degree, the ordering of Chiba
//!   and Nishizeki (1985, `chiba1985arboricity`), in O(m^1.5) time.
//! - Distances. One single-source search per node, chosen by
//!   [`crate::weights`]: breadth-first search, Dial's bucket queue
//!   (`dial1969algorithm`), or Dijkstra's algorithm (`dijkstra1959note`).
//!   The mean distance is the characteristic path length of Watts and
//!   Strogatz (1998). Eccentricity, radius, and diameter are as in Hage and
//!   Harary (1995, `hage1995eccentricity`) and Wasserman and Faust (1994).
//! - Cut points and bridges. One depth-first search computing low points
//!   (Hopcroft and Tarjan 1973, `hopcroft1973algorithm`), written with an
//!   explicit stack so deep networks cannot overflow the call stack.
//!
//! Conventions:
//! - Reciprocity reads the network as an adjacency matrix: parallel arcs
//!   count once and self-loops are ignored, which is igraph's
//!   `reciprocity()` with its default `ignore.loops = TRUE` on a simplified
//!   graph, sna's `grecip()`, and NetworkX's `overall_reciprocity()` on a
//!   `DiGraph` without self-loops. With M mutual, A asymmetric, and N null
//!   dyads among the n(n - 1)/2: `Default` is 2M / (2M + A) (igraph's
//!   default; sna's "edgewise"; NetworkX), `Ratio` is M / (M + A) (igraph's
//!   "ratio"; sna's "dyadic.nonnull"), `Dyadic` is (M + N) / (M + A + N)
//!   (sna's "dyadic"), and `Correlation` is Garlaschelli and Loffredo's rho,
//!   which equals the Pearson correlation between x_ij and x_ji over ordered
//!   pairs of distinct nodes (sna's "correlation"). Ratios with a zero
//!   denominator are NaN, including the correlation of a network with no
//!   ties or every tie, which sna reports as 1. An undirected network is
//!   read as a symmetric directed one, so every tie is reciprocated.
//! - Transitivity and clustering use the simple undirected reading of the
//!   network: direction is ignored, parallel edges count once, and
//!   self-loops are ignored, as igraph's `transitivity()` does; NetworkX's
//!   `transitivity()` and `clustering()` agree on undirected simple graphs.
//!   A node's degree k is its number of distinct neighbors other than
//!   itself. The local coefficient of a node with k < 2, and the global
//!   coefficient of a network without connected triples, are NaN
//!   ([`Undefined::Nan`], igraph) or 0 ([`Undefined::Zero`], NetworkX).
//!   The average coefficient either leaves those nodes out (igraph's
//!   `isolates = "NaN"`) or counts them as 0 (igraph's `isolates = "zero"`,
//!   NetworkX's `average_clustering()`).
//! - Weighted clustering uses the weighted adjacency matrix of the simple
//!   undirected reading: the weight between two nodes is the sum of the
//!   weights of every edge joining them, in either direction. Weights must
//!   be positive; they are tie strengths and are used as given. With s_i
//!   the sum of node i's weights, Barrat's coefficient is
//!   sum over triangles (i, j, h) of (w_ij + w_ih), divided by
//!   s_i (k_i - 1) (igraph's `type = "barrat"`). Onnela's is
//!   2 sum over triangles of (w_ij w_ih w_jh)^(1/3) / (k_i (k_i - 1)), with
//!   every weight divided by the largest weight in the network (NetworkX's
//!   weighted `clustering()`). On networks without parallel edges these are
//!   the published definitions.
//! - Distances: `Out` gives lengths of paths from each source, `In` lengths
//!   of paths to it, and `All` ignores direction; undirected networks
//!   ignore the mode. Unreachable nodes are at distance infinity and a node
//!   is at distance 0 from itself. Self-loops never shorten a path, and
//!   parallel edges count through the lightest of them. Weights are handled
//!   as in [`crate::metrics::closeness`]: distances or strengths (lengths
//!   1/w), strictly positive, with uniform weights giving hop counts times
//!   the common length. igraph's `distances()` and NetworkX's shortest path
//!   lengths agree.
//! - The mean distance averages over ordered pairs of distinct nodes joined
//!   by a path, and the diameter is the longest such distance, as igraph's
//!   `mean_distance()` and `diameter()` with `unconnected = TRUE`. With
//!   `unconnected = false`, a network in which some node does not reach
//!   every other has mean distance and diameter infinity, as igraph with
//!   `unconnected = FALSE`. A network without pairs has mean distance NaN;
//!   its diameter is 0 with one node and NaN with none.
//! - Eccentricity is a node's largest distance to the nodes it reaches, so
//!   isolated nodes have eccentricity 0, as in igraph's `eccentricity()`;
//!   the radius is the smallest eccentricity (NaN without nodes). NetworkX
//!   agrees on (strongly) connected networks and refuses others.
//! - Reachability counts the other nodes each node reaches (`Out`), is
//!   reached from (`In`), or is joined to ignoring direction (`All`).
//! - Cut points and bridges use the undirected reading in which every edge
//!   stays an edge: direction is ignored, two parallel edges (or two arcs
//!   in opposite directions) between the same nodes are two edges, so
//!   neither is a bridge, and self-loops are ignored. A cut point is a node
//!   whose removal increases the number of components; a bridge is an edge
//!   whose removal does. NetworkX's `articulation_points()` and `bridges()`
//!   on an undirected multigraph agree.
//! - Weights are ignored by everything except the distance-based measures
//!   and weighted clustering.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::bitset::{for_each_bit, MultiBfs, Step, BATCH};
use crate::error::{Error, Result};
use crate::graph::{Csr, Graph};
use crate::matrix::Dense;
use crate::metrics::closeness::oriented;
use crate::metrics::degree::DegreeMode;
use crate::metrics::topography::{components, ComponentKind};
use crate::parallel;
use crate::weights::{path_lengths, Lengths, PathAlgorithm, WeightKind};

// ---- Reciprocity ---------------------------------------------------------------

/// Which reciprocity index. See the module documentation for formulas.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReciprocityMode {
    /// Share of arcs that are reciprocated, 2M / (2M + A) (igraph's
    /// default, sna's "edgewise", NetworkX).
    #[default]
    Default,
    /// Mutual dyads over non-null dyads, M / (M + A) (igraph's "ratio",
    /// sna's "dyadic.nonnull").
    Ratio,
    /// Symmetric dyads (mutual or null) over all dyads (sna's "dyadic").
    Dyadic,
    /// Garlaschelli and Loffredo's rho: the correlation between the two
    /// arcs of a dyad.
    Correlation,
}

impl ReciprocityMode {
    /// Decodes the C interface codes: 0 = default, 1 = ratio, 2 = dyadic,
    /// 3 = correlation.
    pub fn from_code(code: i32) -> Result<ReciprocityMode> {
        match code {
            0 => Ok(ReciprocityMode::Default),
            1 => Ok(ReciprocityMode::Ratio),
            2 => Ok(ReciprocityMode::Dyadic),
            3 => Ok(ReciprocityMode::Correlation),
            _ => Err(Error::InvalidArgument(format!(
                "unknown reciprocity mode code {code}"
            ))),
        }
    }
}

/// Distinct arcs between distinct nodes (L) and how many of them are
/// reciprocated (2M).
fn dyad_counts(g: &Graph) -> (u128, u128) {
    let n = g.node_count();
    let (mut arcs, mut reciprocated) = (0u128, 0u128);
    if !g.is_directed() {
        let l = g.simple_undirected(false).arc_count() as u128;
        return (l, l);
    }
    let (out, inc) = (g.out(), g.inc());
    for v in 0..n {
        let (a, b) = (out.neighbors(v), inc.neighbors(v));
        let (mut i, mut j) = (0, 0);
        while i < a.len() {
            let u = a[i];
            while i < a.len() && a[i] == u {
                i += 1;
            }
            if u as usize == v {
                continue;
            }
            arcs += 1;
            while j < b.len() && b[j] < u {
                j += 1;
            }
            if j < b.len() && b[j] == u {
                reciprocated += 1;
            }
        }
    }
    (arcs, reciprocated)
}

/// Reciprocity of the network. See the module documentation.
pub fn reciprocity(g: &Graph, mode: ReciprocityMode) -> f64 {
    let n = g.node_count() as u128;
    let (l, r) = dyad_counts(g);
    let ordered = n * n.saturating_sub(1);
    match mode {
        ReciprocityMode::Default => {
            if l == 0 {
                f64::NAN
            } else {
                r as f64 / l as f64
            }
        }
        ReciprocityMode::Ratio => {
            // M = r / 2 mutual dyads; M + A = l - M non-null dyads.
            if l == 0 {
                f64::NAN
            } else {
                r as f64 / (2 * l - r) as f64
            }
        }
        ReciprocityMode::Dyadic => {
            // Symmetric dyads are all dyads but the A = l - r asymmetric ones.
            if ordered == 0 {
                f64::NAN
            } else {
                let dyads = ordered / 2;
                (dyads - (l - r)) as f64 / dyads as f64
            }
        }
        ReciprocityMode::Correlation => {
            // rho = (r/l - l/N) / (1 - l/N) = (r N - l^2) / (l (N - l)),
            // with N = n(n - 1) ordered pairs, evaluated in exact integers.
            if l == 0 || l == ordered {
                f64::NAN
            } else {
                let num = (r * ordered) as i128 - (l * l) as i128;
                let den = l * (ordered - l);
                num as f64 / den as f64
            }
        }
    }
}

// ---- Transitivity and clustering --------------------------------------------

/// The value of a clustering coefficient that is 0 / 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Undefined {
    /// NaN, and left out of averages (igraph).
    #[default]
    Nan,
    /// 0, and counted in averages (NetworkX).
    Zero,
}

impl Undefined {
    /// Decodes the C interface codes: 0 = NaN, 1 = zero.
    pub fn from_code(code: i32) -> Result<Undefined> {
        match code {
            0 => Ok(Undefined::Nan),
            1 => Ok(Undefined::Zero),
            _ => Err(Error::InvalidArgument(format!(
                "unknown code {code} for undefined clustering coefficients"
            ))),
        }
    }

    fn apply(self, x: f64) -> f64 {
        if x.is_nan() && self == Undefined::Zero {
            0.0
        } else {
            x
        }
    }
}

/// Which local clustering coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Weighting {
    /// The binary coefficient of Watts and Strogatz.
    #[default]
    None,
    /// Barrat et al. (2004), igraph's `type = "barrat"`.
    Barrat,
    /// Onnela et al. (2005), NetworkX's weighted `clustering()`.
    Onnela,
}

impl Weighting {
    /// Decodes the C interface codes: 0 = none, 1 = Barrat, 2 = Onnela.
    pub fn from_code(code: i32) -> Result<Weighting> {
        match code {
            0 => Ok(Weighting::None),
            1 => Ok(Weighting::Barrat),
            2 => Ok(Weighting::Onnela),
            _ => Err(Error::InvalidArgument(format!(
                "unknown clustering weighting code {code}"
            ))),
        }
    }
}

/// Options for [`transitivity`], [`local_clustering`], and
/// [`average_clustering`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ClusteringOptions {
    /// Binary or weighted local coefficient ([`transitivity`] is always
    /// binary).
    pub weighting: Weighting,
    /// Value of a coefficient that is 0 / 0.
    pub undefined: Undefined,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

/// Number of triangles at each node of the simple undirected reading of
/// the network.
pub fn triangles(g: &Graph, threads: usize) -> Vec<u64> {
    triangles_of(&g.simple_undirected(false), threads)
}

fn triangles_of(s: &Csr, threads: usize) -> Vec<u64> {
    let n = s.node_count();
    // Orient every edge toward the endpoint later in (degree, index) order;
    // each triangle is then found once, from its earliest node, and no node
    // has more than O(sqrt(m)) later neighbors.
    let later = |a: usize, b: usize| (s.degree(a), a) < (s.degree(b), b);
    let mut off = Vec::with_capacity(n + 1);
    off.push(0usize);
    let mut fwd: Vec<u32> = Vec::with_capacity(s.arc_count() / 2);
    for v in 0..n {
        fwd.extend(
            s.neighbors(v)
                .iter()
                .copied()
                .filter(|&u| later(v, u as usize)),
        );
        off.push(fwd.len());
    }
    let parts = parallel::for_each_item(
        n,
        threads,
        || (vec![u32::MAX; n], vec![0u64; n]),
        |(mark, count), u| {
            let fu = &fwd[off[u]..off[u + 1]];
            if fu.len() < 2 {
                return;
            }
            for &v in fu {
                mark[v as usize] = u as u32;
            }
            for &v in fu {
                let v = v as usize;
                for &w in &fwd[off[v]..off[v + 1]] {
                    if mark[w as usize] == u as u32 {
                        count[u] += 1;
                        count[v] += 1;
                        count[w as usize] += 1;
                    }
                }
            }
        },
    );
    let mut total = vec![0u64; n];
    for (_, count) in parts {
        for (t, c) in total.iter_mut().zip(count) {
            *t += c;
        }
    }
    total
}

/// Global transitivity: three times the number of triangles over the number
/// of connected triples. `weighting` is ignored. See the module
/// documentation.
pub fn transitivity(g: &Graph, opts: ClusteringOptions) -> f64 {
    let s = g.simple_undirected(false);
    let t = triangles_of(&s, opts.threads);
    let closed: u128 = t.iter().map(|&x| x as u128).sum();
    let triples: u128 = (0..s.node_count())
        .map(|v| {
            let k = s.degree(v) as u128;
            k * k.saturating_sub(1) / 2
        })
        .sum();
    let value = if triples == 0 {
        f64::NAN
    } else {
        closed as f64 / triples as f64
    };
    opts.undefined.apply(value)
}

/// Local clustering coefficient of every node, in node order. See the module
/// documentation.
pub fn local_clustering(g: &Graph, opts: ClusteringOptions) -> Result<Vec<f64>> {
    let raw = raw_clustering(g, opts)?;
    Ok(raw.into_iter().map(|x| opts.undefined.apply(x)).collect())
}

/// Mean local clustering coefficient. With [`Undefined::Nan`] nodes with
/// fewer than two neighbors are left out (NaN if no node has two); with
/// [`Undefined::Zero`] they count as 0 (NaN for a network without nodes).
pub fn average_clustering(g: &Graph, opts: ClusteringOptions) -> Result<f64> {
    let raw = raw_clustering(g, opts)?;
    let (mut total, mut count) = (0.0, 0usize);
    for x in raw {
        match (x.is_nan(), opts.undefined) {
            (true, Undefined::Nan) => {}
            (true, Undefined::Zero) => count += 1,
            (false, _) => {
                total += x;
                count += 1;
            }
        }
    }
    Ok(if count == 0 {
        f64::NAN
    } else {
        total / count as f64
    })
}

/// Local coefficients with NaN for nodes with fewer than two neighbors.
fn raw_clustering(g: &Graph, opts: ClusteringOptions) -> Result<Vec<f64>> {
    if opts.weighting == Weighting::None {
        let s = g.simple_undirected(false);
        let t = triangles_of(&s, opts.threads);
        return Ok((0..s.node_count())
            .map(|v| {
                let k = s.degree(v) as f64;
                if k < 2.0 {
                    f64::NAN
                } else {
                    t[v] as f64 / (k * (k - 1.0) / 2.0)
                }
            })
            .collect());
    }
    let raw = g.edge_weights().ok_or(Error::MissingWeights)?;
    if let Some(e) = raw.iter().position(|&w| w <= 0.0) {
        return Err(Error::InvalidArgument(format!(
            "weighted clustering needs positive weights, but edge {e} has weight {}",
            raw[e]
        )));
    }
    let s = g.simple_undirected(true);
    Ok(weighted_clustering(&s, opts.weighting, opts.threads))
}

/// Barrat or Onnela coefficients, NaN where k < 2. Each node's sum is
/// accumulated by one thread in a fixed order, so results do not depend on
/// the number of threads.
fn weighted_clustering(s: &Csr, weighting: Weighting, threads: usize) -> Vec<f64> {
    let n = s.node_count();
    let all = s.all_weights().unwrap_or(&[]);
    let top = all.iter().copied().fold(0.0f64, f64::max);
    let parts = parallel::for_each_item(
        n,
        threads,
        || (vec![u32::MAX; n], vec![0.0f64; n], Vec::new()),
        |(stamp, wi, found), i| {
            let ni = s.neighbors(i);
            let k = ni.len();
            if k < 2 {
                found.push((i, f64::NAN));
                return;
            }
            let wsi = s.weights(i).unwrap_or(&[]);
            for (&j, &w) in ni.iter().zip(wsi) {
                stamp[j as usize] = i as u32;
                wi[j as usize] = w;
            }
            let mut sum = 0.0;
            for (&j, &wij) in ni.iter().zip(wsi) {
                let (nj, wj) = (
                    s.neighbors(j as usize),
                    s.weights(j as usize).unwrap_or(&[]),
                );
                let start = nj.partition_point(|&h| h <= j);
                for (&h, &wjh) in nj[start..].iter().zip(&wj[start..]) {
                    if stamp[h as usize] != i as u32 {
                        continue;
                    }
                    let wih = wi[h as usize];
                    sum += match weighting {
                        Weighting::Barrat => wij + wih,
                        _ => 2.0 * ((wij / top) * (wih / top) * (wjh / top)).cbrt(),
                    };
                }
            }
            let k = k as f64;
            let value = match weighting {
                Weighting::Barrat => {
                    let strength: f64 = wsi.iter().sum();
                    sum / (strength * (k - 1.0))
                }
                _ => sum / (k * (k - 1.0)),
            };
            found.push((i, value));
        },
    );
    let mut out = vec![f64::NAN; n];
    for (_, _, found) in parts {
        for (i, v) in found {
            out[i] = v;
        }
    }
    out
}

// ---- Shortest-path distances ------------------------------------------------

/// Options for the distance-based measures.
#[derive(Debug, Clone, Copy)]
pub struct PathOptions {
    /// Direction of the paths in directed graphs.
    pub mode: DegreeMode,
    /// Use edge weights (they must be strictly positive).
    pub weighted: bool,
    /// Whether weights are distances (the default) or strengths.
    pub weights_as: WeightKind,
    /// Shortest-path algorithm; `Auto` picks the cheapest exact one.
    pub algorithm: PathAlgorithm,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for PathOptions {
    fn default() -> Self {
        PathOptions {
            mode: DegreeMode::Out,
            weighted: false,
            weights_as: WeightKind::Distance,
            algorithm: PathAlgorithm::Auto,
            threads: 0,
        }
    }
}

/// A graph oriented for the requested mode, its arc lengths, and the
/// common length when every weight is the same.
struct Prepared<'a> {
    h: Cow<'a, Graph>,
    lengths: Lengths,
    unit: f64,
}

fn prepare<'a>(g: &'a Graph, opts: &PathOptions) -> Result<Prepared<'a>> {
    let h = oriented(g, opts.mode)?;
    let lengths = path_lengths(h.as_ref(), opts.weighted, opts.weights_as, opts.algorithm)?;
    // Uniform weights reduce to breadth-first search, which counts hops;
    // distances are the hop counts times the common length.
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
    Ok(Prepared { h, lengths, unit })
}

fn check_sources(n: usize, sources: Option<&[u32]>) -> Result<Vec<u32>> {
    match sources {
        None => Ok((0..n as u32).collect()),
        Some(s) => {
            if let Some(&bad) = s.iter().find(|&&v| v as usize >= n) {
                return Err(Error::InvalidArgument(format!(
                    "source node index {bad} is out of range for a network with {n} nodes"
                )));
            }
            Ok(s.to_vec())
        }
    }
}

/// Calls `sink(k, row)` with the distances from `sources[k]` (all nodes
/// when `sources` is `None`) to every node, for k in increasing order.
/// Unreachable nodes are at distance infinity. Rows are computed in
/// parallel, in blocks, so memory beyond one block of rows is not needed.
pub fn distance_rows(
    g: &Graph,
    sources: Option<&[u32]>,
    opts: PathOptions,
    mut sink: impl FnMut(usize, &[f64]),
) -> Result<()> {
    let n = g.node_count();
    let sources = check_sources(n, sources)?;
    let p = prepare(g, &opts)?;
    let h = p.h.as_ref();
    let threads = parallel::resolve_threads(opts.threads, sources.len());
    // Keep a block of rows under about 64 MB, but give every thread work.
    let block = ((1usize << 23) / n.max(1)).max(threads).max(1);
    let uniform = matches!(p.lengths, Lengths::Uniform);
    let block = if uniform {
        // Whole batches of the multi-source search.
        block.max(BATCH * threads)
    } else {
        block
    };
    let mut start = 0;
    while start < sources.len() {
        let end = (start + block).min(sources.len());
        let batch = &sources[start..end];
        let parts = if uniform {
            // Breadth-first search from 64 sources at a time
            // (`crate::bitset::MultiBfs`), one batch per work item.
            parallel::for_each_item(
                (batch.len() + BATCH - 1) / BATCH,
                threads,
                || (MultiBfs::new(n), Vec::new()),
                |(ms, rows), b| {
                    let lo = b * BATCH;
                    let hi = (lo + BATCH).min(batch.len());
                    let group = &batch[lo..hi];
                    let mut out: Vec<Vec<f64>> = (0..group.len())
                        .map(|_| vec![f64::INFINITY; n])
                        .collect();
                    for (i, &s) in group.iter().enumerate() {
                        out[i][s as usize] = 0.0;
                    }
                    ms.run(h.out(), group, u32::MAX, |step| {
                        if let Step::Reached {
                            level,
                            node,
                            sources,
                        } = step
                        {
                            let d = level as f64 * p.unit;
                            for_each_bit(sources, |i| out[i][node as usize] = d);
                        }
                    });
                    for (i, row) in out.into_iter().enumerate() {
                        rows.push((start + lo + i, row));
                    }
                },
            )
            .into_iter()
            .map(|(_, rows)| rows)
            .collect::<Vec<_>>()
        } else {
            parallel::for_each_item(
                batch.len(),
                threads,
                || (Search::new(n), Vec::new()),
                |(ws, rows), k| {
                    ws.run(h, &p.lengths, batch[k] as usize);
                    let mut row = vec![f64::INFINITY; n];
                    for &v in &ws.order {
                        row[v as usize] = ws.distance(v as usize, &p.lengths) * p.unit;
                    }
                    ws.reset();
                    rows.push((start + k, row));
                },
            )
            .into_iter()
            .map(|(_, rows)| rows)
            .collect()
        };
        let mut rows: Vec<(usize, Vec<f64>)> = parts.into_iter().flatten().collect();
        rows.sort_unstable_by_key(|r| r.0);
        for (k, row) in rows {
            sink(k, &row);
        }
        start = end;
    }
    Ok(())
}

/// Distance matrix with one row per source (every node when `sources` is
/// `None`) and one column per node. See the module documentation.
pub fn distances(g: &Graph, sources: Option<&[u32]>, opts: PathOptions) -> Result<Dense> {
    let n = g.node_count();
    let rows = sources.map_or(n, |s| s.len());
    let mut m = Dense::zeros(rows, n);
    distance_rows(g, sources, opts, |k, row| {
        for (j, &d) in row.iter().enumerate() {
            m.set(k, j, d);
        }
    })?;
    Ok(m)
}

/// What one single-source search records.
#[derive(Debug, Clone, Copy, Default)]
struct PathProfile {
    /// Other nodes reached.
    reached: usize,
    /// Sum of the distances to them, exact when lengths are integers.
    total_int: u128,
    total: f64,
    /// Longest finite distance.
    longest: f64,
}

/// One profile per source, whether lengths are integers, and the common
/// length that scales them.
fn profiles(
    g: &Graph,
    sources: &[u32],
    opts: &PathOptions,
) -> Result<(Vec<PathProfile>, bool, f64)> {
    let n = g.node_count();
    let p = prepare(g, opts)?;
    let h = p.h.as_ref();
    let integral = !matches!(p.lengths, Lengths::Real(_));
    let parts: Vec<Vec<(usize, PathProfile)>> = if matches!(p.lengths, Lengths::Uniform) {
        // Breadth-first search from 64 sources at a time
        // (`crate::bitset::MultiBfs`), one batch per work item. Every
        // total is an exact integer count of hops.
        parallel::for_each_item(
            (sources.len() + BATCH - 1) / BATCH,
            opts.threads,
            || (MultiBfs::new(n), Vec::new()),
            |(ms, found), b| {
                let lo = b * BATCH;
                let hi = (lo + BATCH).min(sources.len());
                let group = &sources[lo..hi];
                let mut count = [0u32; BATCH];
                let mut profiles = [PathProfile::default(); BATCH];
                ms.run(h.out(), group, u32::MAX, |step| match step {
                    Step::Reached { sources, .. } => for_each_bit(sources, |i| count[i] += 1),
                    Step::LevelDone(level) => {
                        for (q, c) in profiles.iter_mut().zip(count.iter_mut()) {
                            if *c == 0 {
                                continue;
                            }
                            q.reached += *c as usize;
                            q.total_int += *c as u128 * level as u128;
                            q.longest = level as f64;
                            *c = 0;
                        }
                    }
                });
                for (i, q) in profiles.iter_mut().enumerate().take(group.len()) {
                    q.total = q.total_int as f64 * p.unit;
                    q.longest *= p.unit;
                    found.push((lo + i, *q));
                }
            },
        )
        .into_iter()
        .map(|(_, found)| found)
        .collect()
    } else {
        parallel::for_each_item(
            sources.len(),
            opts.threads,
            || (Search::new(n), Vec::new()),
            |(ws, found), k| {
                let s = sources[k] as usize;
                ws.run(h, &p.lengths, s);
                let mut q = PathProfile {
                    reached: ws.order.len() - 1,
                    ..Default::default()
                };
                for &v in &ws.order {
                    let v = v as usize;
                    if integral {
                        let d = ws.hops[v];
                        q.total_int += d as u128;
                        q.longest = q.longest.max(d as f64);
                    } else {
                        let d = ws.dist[v];
                        q.total += d;
                        q.longest = q.longest.max(d);
                    }
                }
                ws.reset();
                if integral {
                    q.total = q.total_int as f64 * p.unit;
                }
                q.longest *= p.unit;
                found.push((k, q));
            },
        )
        .into_iter()
        .map(|(_, found)| found)
        .collect()
    };
    let mut out = vec![PathProfile::default(); sources.len()];
    for found in parts {
        for (k, q) in found {
            out[k] = q;
        }
    }
    Ok((out, integral, p.unit))
}

/// Totals over the shortest-path searches from a set of sources.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathSummary {
    /// Number of sources searched.
    pub sources: usize,
    /// Ordered pairs (source, other node) joined by a path.
    pub pairs: u64,
    /// Sum of their distances.
    pub total: f64,
    /// Longest finite distance (0 if no pairs).
    pub longest: f64,
    /// Number of sources that do not reach every other node.
    pub incomplete: usize,
}

impl PathSummary {
    /// Mean distance over the pairs joined by a path (NaN without pairs),
    /// or infinity if `unconnected` is false and some source does not
    /// reach every other node.
    pub fn mean(&self, unconnected: bool) -> f64 {
        if !unconnected && self.incomplete > 0 {
            f64::INFINITY
        } else if self.pairs == 0 {
            f64::NAN
        } else {
            self.total / self.pairs as f64
        }
    }

    /// Longest finite distance (NaN without sources), or infinity if
    /// `unconnected` is false and some source does not reach every other
    /// node.
    pub fn diameter(&self, unconnected: bool) -> f64 {
        if self.sources == 0 {
            f64::NAN
        } else if !unconnected && self.incomplete > 0 {
            f64::INFINITY
        } else {
            self.longest
        }
    }
}

/// Shortest-path totals from `sources` (every node when `None`). Totals
/// are combined in source order, so results do not depend on the number of
/// threads; with integer lengths the sum is exact.
pub fn path_summary(g: &Graph, sources: Option<&[u32]>, opts: PathOptions) -> Result<PathSummary> {
    let n = g.node_count();
    let sources = check_sources(n, sources)?;
    let (profiles, integral, unit) = profiles(g, &sources, &opts)?;
    let mut s = PathSummary {
        sources: sources.len(),
        pairs: 0,
        total: 0.0,
        longest: 0.0,
        incomplete: 0,
    };
    let mut exact = 0u128;
    for q in &profiles {
        s.pairs += q.reached as u64;
        exact += q.total_int;
        s.total += q.total;
        s.longest = s.longest.max(q.longest);
        if q.reached + 1 < n {
            s.incomplete += 1;
        }
    }
    // Exact integer totals are scaled once by the common length.
    if integral {
        s.total = exact as f64 * unit;
    }
    Ok(s)
}

/// Mean distance over ordered pairs of distinct nodes joined by a path. See
/// the module documentation and [`PathSummary::mean`].
pub fn mean_distance(g: &Graph, opts: PathOptions, unconnected: bool) -> Result<f64> {
    Ok(path_summary(g, None, opts)?.mean(unconnected))
}

/// Longest shortest-path distance. See [`PathSummary::diameter`].
pub fn diameter(g: &Graph, opts: PathOptions, unconnected: bool) -> Result<f64> {
    Ok(path_summary(g, None, opts)?.diameter(unconnected))
}

/// Each node's largest distance to the nodes it reaches (0 if it reaches
/// none), in node order.
pub fn eccentricity(g: &Graph, opts: PathOptions) -> Result<Vec<f64>> {
    let sources = check_sources(g.node_count(), None)?;
    Ok(profiles(g, &sources, &opts)?
        .0
        .iter()
        .map(|q| q.longest)
        .collect())
}

/// The smallest eccentricity (NaN for a network without nodes).
pub fn radius(g: &Graph, opts: PathOptions) -> Result<f64> {
    let e = eccentricity(g, opts)?;
    Ok(e.iter().copied().reduce(f64::min).unwrap_or(f64::NAN))
}

/// Number of other nodes each node reaches (`Out`), is reached from
/// (`In`), or is joined to ignoring direction (`All`), in node order.
pub fn reachability(g: &Graph, mode: DegreeMode, threads: usize) -> Result<Vec<f64>> {
    if !g.is_directed() || mode == DegreeMode::All {
        let c = components(g, ComponentKind::Weak);
        return Ok(c
            .membership
            .iter()
            .map(|&m| (c.sizes[m as usize] - 1) as f64)
            .collect());
    }
    let h = oriented(g, mode)?;
    let h = h.as_ref();
    let counts = within_hops(h.out(), None, u32::MAX, threads);
    Ok(counts.into_iter().map(|c| c as f64).collect())
}

/// Number of other nodes within `max_level` steps of each source (every
/// node when `None`), by breadth-first search from 64 sources at a time.
fn within_hops(adj: &Csr, sources: Option<&[u32]>, max_level: u32, threads: usize) -> Vec<usize> {
    let n = adj.node_count();
    let all: Vec<u32>;
    let sources = match sources {
        Some(s) => s,
        None => {
            all = (0..n as u32).collect();
            &all
        }
    };
    let parts = parallel::for_each_item(
        (sources.len() + BATCH - 1) / BATCH,
        threads,
        || (MultiBfs::new(n), Vec::new()),
        |(ms, found), b| {
            let lo = b * BATCH;
            let hi = (lo + BATCH).min(sources.len());
            let mut count = [0usize; BATCH];
            ms.run(adj, &sources[lo..hi], max_level, |step| {
                if let Step::Reached { sources, .. } = step {
                    for_each_bit(sources, |i| count[i] += 1);
                }
            });
            found.push((lo, count));
        },
    );
    let mut out = vec![0usize; sources.len()];
    for (_, found) in parts {
        for (lo, count) in found {
            let hi = (lo + BATCH).min(sources.len());
            out[lo..hi].copy_from_slice(&count[..hi - lo]);
        }
    }
    out
}

/// Number of other nodes at distance at most `limit` from each node, with
/// distances as in [`distances`]. Searches stop once every node within the
/// limit is settled.
pub fn within_distance(g: &Graph, limit: f64, opts: PathOptions) -> Result<Vec<usize>> {
    let p = prepare(g, &opts)?;
    let h = p.h.as_ref();
    if let Lengths::Uniform = p.lengths {
        // Levels whose distance, level times the common length, is within
        // the limit.
        let mut max_level = 0u32;
        while max_level < u32::MAX && (max_level + 1) as f64 * p.unit <= limit {
            max_level += 1;
        }
        return Ok(within_hops(h.out(), None, max_level, opts.threads));
    }
    let n = h.node_count();
    let parts = parallel::for_each_item(
        n,
        opts.threads,
        || (Search::new(n), Vec::new()),
        |(ws, found), s| {
            ws.run_within(h, &p.lengths, s, limit);
            found.push((s, ws.order.len() - 1));
            ws.reset();
        },
    );
    let mut out = vec![0usize; n];
    for (_, found) in parts {
        for (s, c) in found {
            out[s] = c;
        }
    }
    Ok(out)
}

/// Per-thread scratch space for single-source searches, reused across
/// sources. After a search, `order` lists the nodes reached (the source
/// first) and `hops` (breadth-first search and Dial) or `dist` (Dijkstra)
/// holds their distances; `reset` clears them.
struct Search {
    hops: Vec<u64>,
    dist: Vec<f64>,
    done: Vec<bool>,
    order: Vec<u32>,
    buckets: Vec<Vec<u32>>,
    heap: BinaryHeap<HeapItem>,
}

impl Search {
    fn new(n: usize) -> Search {
        Search {
            hops: vec![u64::MAX; n],
            dist: vec![f64::INFINITY; n],
            done: vec![false; n],
            order: Vec::new(),
            buckets: Vec::new(),
            heap: BinaryHeap::new(),
        }
    }

    fn run(&mut self, g: &Graph, lengths: &Lengths, s: usize) {
        self.run_within(g, lengths, s, f64::INFINITY)
    }

    /// Like `run`, but settles only the nodes at distance at most `limit`
    /// (the tentative state of farther nodes is reset before returning).
    fn run_within(&mut self, g: &Graph, lengths: &Lengths, s: usize, limit: f64) {
        match lengths {
            Lengths::Uniform => self.bfs(g, s, limit),
            Lengths::Integer(lens, max) => self.dial(g, lens, *max, s, limit),
            Lengths::Real(lens) => self.dijkstra(g, lens, s, limit),
        }
    }

    fn distance(&self, v: usize, lengths: &Lengths) -> f64 {
        match lengths {
            Lengths::Real(_) => self.dist[v],
            _ => self.hops[v] as f64,
        }
    }

    fn reset(&mut self) {
        for &v in &self.order {
            let v = v as usize;
            self.hops[v] = u64::MAX;
            self.dist[v] = f64::INFINITY;
            self.done[v] = false;
        }
        self.order.clear();
    }

    fn bfs(&mut self, g: &Graph, s: usize, limit: f64) {
        let adj = g.out();
        self.hops[s] = 0;
        self.order.push(s as u32);
        let mut head = 0;
        while head < self.order.len() {
            let v = self.order[head] as usize;
            head += 1;
            let next = self.hops[v] + 1;
            if next as f64 > limit {
                break;
            }
            for &w in adj.neighbors(v) {
                let w = w as usize;
                if self.hops[w] == u64::MAX {
                    self.hops[w] = next;
                    self.order.push(w as u32);
                }
            }
        }
    }

    /// Dial's bucket queue with exact integer lengths (`dial1969algorithm`).
    /// Every node given a tentative distance is settled before the queue
    /// empties, so `order` lists every node touched.
    fn dial(&mut self, g: &Graph, lens: &[u32], max: u32, s: usize, limit: f64) {
        let adj = g.out();
        let width = max as usize + 1;
        if self.buckets.len() != width {
            self.buckets = vec![Vec::new(); width];
        }
        self.hops[s] = 0;
        self.buckets[0].push(s as u32);
        let mut pending = 1usize;
        let mut cur: u64 = 0;
        while pending > 0 {
            if cur as f64 > limit {
                // Every node within the limit is settled; forget the rest.
                for bucket in &mut self.buckets {
                    for v in bucket.drain(..) {
                        self.hops[v as usize] = u64::MAX;
                    }
                }
                break;
            }
            let slot = (cur % width as u64) as usize;
            let v = match self.buckets[slot].pop() {
                Some(v) => v as usize,
                None => {
                    cur += 1;
                    continue;
                }
            };
            pending -= 1;
            if self.done[v] || self.hops[v] != cur {
                continue;
            }
            self.done[v] = true;
            self.order.push(v as u32);
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                if self.done[w] {
                    continue;
                }
                let alt = cur + lens[k] as u64;
                if alt < self.hops[w] {
                    self.hops[w] = alt;
                    self.buckets[(alt % width as u64) as usize].push(w as u32);
                    pending += 1;
                }
            }
        }
    }

    /// Dijkstra's algorithm with a binary heap and lazy deletion
    /// (`dijkstra1959note`).
    fn dijkstra(&mut self, g: &Graph, lens: &[f64], s: usize, limit: f64) {
        let adj = g.out();
        self.heap.clear();
        self.dist[s] = 0.0;
        self.heap.push(HeapItem {
            dist: 0.0,
            node: s as u32,
        });
        while let Some(HeapItem { dist: d, node }) = self.heap.pop() {
            let v = node as usize;
            if d > limit {
                // Every node within the limit is settled; forget the rest.
                self.dist[v] = f64::INFINITY;
                for item in self.heap.drain() {
                    self.dist[item.node as usize] = f64::INFINITY;
                }
                break;
            }
            if self.done[v] || d > self.dist[v] {
                continue;
            }
            self.done[v] = true;
            self.order.push(v as u32);
            for (k, &w) in adj.arc_range(v).zip(adj.neighbors(v)) {
                let w = w as usize;
                if self.done[w] {
                    continue;
                }
                let alt = d + lens[k];
                if alt < self.dist[w] {
                    self.dist[w] = alt;
                    self.heap.push(HeapItem {
                        dist: alt,
                        node: w as u32,
                    });
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct HeapItem {
    dist: f64,
    node: u32,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for HeapItem {}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapItem {
    // Reversed so that BinaryHeap (a max-heap) pops the smallest distance.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .dist
            .total_cmp(&self.dist)
            .then_with(|| other.node.cmp(&self.node))
    }
}

// ---- Cut points and bridges -------------------------------------------------

/// Cut points and bridges of the undirected reading of a network.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutStructure {
    /// Whether each node is a cut point.
    pub cut_point: Vec<bool>,
    /// Whether each edge (in input order) is a bridge.
    pub bridge: Vec<bool>,
}

/// Cut points and bridges by one depth-first search with low points
/// (Hopcroft and Tarjan 1973). See the module documentation.
pub fn cut_structure(g: &Graph) -> CutStructure {
    let n = g.node_count();
    let m = g.edge_count();
    let (src, tgt) = (g.edge_sources(), g.edge_targets());

    // Incidence lists: for every node, (neighbor, edge) for each non-loop
    // edge, built with a counting sort.
    let mut off = vec![0usize; n + 1];
    for e in 0..m {
        if src[e] != tgt[e] {
            off[src[e] as usize + 1] += 1;
            off[tgt[e] as usize + 1] += 1;
        }
    }
    for v in 0..n {
        off[v + 1] += off[v];
    }
    let mut cursor = off.clone();
    let mut adj = vec![(0u32, 0usize); off[n]];
    for e in 0..m {
        let (a, b) = (src[e], tgt[e]);
        if a != b {
            adj[cursor[a as usize]] = (b, e);
            cursor[a as usize] += 1;
            adj[cursor[b as usize]] = (a, e);
            cursor[b as usize] += 1;
        }
    }

    const UNSEEN: usize = usize::MAX;
    const NO_EDGE: usize = usize::MAX;
    let mut disc = vec![UNSEEN; n];
    let mut low = vec![0usize; n];
    let mut cut_point = vec![false; n];
    let mut bridge = vec![false; m];
    // Each frame is a node, the edge it was entered by, and the position
    // of its next incidence to explore.
    let mut frames: Vec<(usize, usize, usize)> = Vec::new();
    let mut timer = 0usize;
    for root in 0..n {
        if disc[root] != UNSEEN {
            continue;
        }
        disc[root] = timer;
        low[root] = timer;
        timer += 1;
        let mut root_children = 0usize;
        frames.push((root, NO_EDGE, off[root]));
        while let Some(frame) = frames.last_mut() {
            let (v, entered_by, pos) = *frame;
            if pos < off[v + 1] {
                frame.2 += 1;
                let (w, e) = adj[pos];
                let w = w as usize;
                // Skip only the edge we came by, so a parallel edge back to
                // the parent counts as a cycle.
                if e == entered_by {
                    continue;
                }
                if disc[w] == UNSEEN {
                    disc[w] = timer;
                    low[w] = timer;
                    timer += 1;
                    if v == root {
                        root_children += 1;
                    }
                    frames.push((w, e, off[w]));
                } else {
                    low[v] = low[v].min(disc[w]);
                }
                continue;
            }
            frames.pop();
            if let Some(&(p, _, _)) = frames.last() {
                low[p] = low[p].min(low[v]);
                if low[v] > disc[p] {
                    bridge[entered_by] = true;
                }
                if p != root && low[v] >= disc[p] {
                    cut_point[p] = true;
                }
            }
        }
        if root_children >= 2 {
            cut_point[root] = true;
        }
    }
    CutStructure { cut_point, bridge }
}

/// Indices of the cut points, in increasing order.
pub fn cut_points(g: &Graph) -> Vec<u32> {
    let c = cut_structure(g);
    (0..g.node_count() as u32)
        .filter(|&v| c.cut_point[v as usize])
        .collect()
}

/// Indices of the bridges in the edge list, in increasing order.
pub fn bridges(g: &Graph) -> Vec<usize> {
    let c = cut_structure(g);
    (0..g.edge_count()).filter(|&e| c.bridge[e]).collect()
}

// ---- Counts for a summary ---------------------------------------------------

/// Counts of self-loops, multiple edges, and isolates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeCounts {
    /// Edges from a node to itself.
    pub loops: usize,
    /// Edges that repeat an earlier edge between the same (ordered, if
    /// directed) pair of nodes, loops included; igraph's
    /// `sum(which_multiple())`.
    pub multiple: usize,
    /// Nodes without an edge to another node (a node whose only edges are
    /// self-loops is an isolate).
    pub isolates: usize,
}

/// Counts of self-loops, multiple edges, and isolates.
pub fn edge_counts(g: &Graph) -> EdgeCounts {
    let n = g.node_count();
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut keys: Vec<u64> = src
        .iter()
        .zip(tgt)
        .map(|(&a, &b)| {
            let (a, b) = if g.is_directed() || a <= b {
                (a, b)
            } else {
                (b, a)
            };
            ((a as u64) << 32) | b as u64
        })
        .collect();
    keys.sort_unstable();
    let multiple = keys.windows(2).filter(|w| w[0] == w[1]).count();
    let loops = src.iter().zip(tgt).filter(|(a, b)| a == b).count();
    let mut joined = vec![false; n];
    for (&a, &b) in src.iter().zip(tgt) {
        if a != b {
            joined[a as usize] = true;
            joined[b as usize] = true;
        }
    }
    EdgeCounts {
        loops,
        multiple,
        isolates: joined.iter().filter(|&&j| !j).count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, s: &[i64], t: &[i64], directed: bool) -> Graph {
        Graph::from_edges(n, s, t, None, directed).unwrap()
    }

    fn close(a: f64, b: f64) -> bool {
        (a.is_nan() && b.is_nan()) || a == b || (a - b).abs() <= 1e-12 * (1.0 + a.abs())
    }

    #[test]
    fn reciprocity_indices() {
        // 0 <-> 1, 1 -> 2, 2 -> 2 (loop, ignored), 0 -> 1 again (parallel,
        // counts once). n = 4: M = 1, A = 1, N = 4.
        let g = graph(4, &[0, 1, 1, 2, 0], &[1, 0, 2, 2, 1], true);
        assert!(close(reciprocity(&g, ReciprocityMode::Default), 2.0 / 3.0));
        assert!(close(reciprocity(&g, ReciprocityMode::Ratio), 0.5));
        assert!(close(reciprocity(&g, ReciprocityMode::Dyadic), 5.0 / 6.0));
        // p = 3/12, r = 2/3: rho = (2/3 - 1/4) / (3/4) = 5/9.
        assert!(close(
            reciprocity(&g, ReciprocityMode::Correlation),
            5.0 / 9.0
        ));
        let empty = graph(3, &[], &[], true);
        assert!(reciprocity(&empty, ReciprocityMode::Default).is_nan());
        assert_eq!(reciprocity(&empty, ReciprocityMode::Dyadic), 1.0);
        let u = graph(3, &[0, 1], &[1, 2], false);
        assert_eq!(reciprocity(&u, ReciprocityMode::Default), 1.0);
        assert_eq!(reciprocity(&u, ReciprocityMode::Correlation), 1.0);
    }

    #[test]
    fn clustering_of_a_kite() {
        // Triangle 0-1-2 with a tail 2-3 and an isolate 4.
        let g = graph(5, &[0, 1, 2, 2], &[1, 2, 0, 3], false);
        let o = ClusteringOptions::default();
        assert_eq!(triangles(&g, 0), vec![1, 1, 1, 0, 0]);
        // 3 closed out of 1 + 1 + 3 = 5 triples.
        assert!(close(transitivity(&g, o), 3.0 / 5.0));
        let c = local_clustering(&g, o).unwrap();
        assert!(close(c[0], 1.0) && close(c[2], 1.0 / 3.0) && c[3].is_nan() && c[4].is_nan());
        assert!(close(
            average_clustering(&g, o).unwrap(),
            (2.0 + 1.0 / 3.0) / 3.0
        ));
        let z = ClusteringOptions {
            undefined: Undefined::Zero,
            ..o
        };
        assert_eq!(local_clustering(&g, z).unwrap()[4], 0.0);
        assert!(close(
            average_clustering(&g, z).unwrap(),
            (2.0 + 1.0 / 3.0) / 5.0
        ));
        let none = graph(2, &[0], &[1], false);
        assert!(transitivity(&none, o).is_nan());
        assert_eq!(transitivity(&none, z), 0.0);
    }

    #[test]
    fn weighted_clustering_by_hand() {
        // Triangle 0-1-2 with weights 1, 2, 3 and a pendant 0-3 of weight 4.
        let g = Graph::from_edges(
            4,
            &[0, 1, 2, 0],
            &[1, 2, 0, 3],
            Some(&[1.0, 2.0, 3.0, 4.0]),
            false,
        )
        .unwrap();
        let b = ClusteringOptions {
            weighting: Weighting::Barrat,
            ..Default::default()
        };
        let c = local_clustering(&g, b).unwrap();
        // Node 0: s = 8, k = 3, one triangle with w01 + w02 = 4.
        assert!(close(c[0], 4.0 / 16.0));
        assert!(close(c[1], 1.0));
        assert!(c[3].is_nan());
        let o = ClusteringOptions {
            weighting: Weighting::Onnela,
            ..Default::default()
        };
        let c = local_clustering(&g, o).unwrap();
        let t = (1.0f64 * 2.0 * 3.0 / 64.0).cbrt();
        assert!(close(c[0], 2.0 * t / 6.0));
        assert!(close(c[1], t));
    }

    #[test]
    fn distances_and_their_summaries() {
        // 0 -> 1 -> 2, 3 isolated.
        let g = graph(4, &[0, 1], &[1, 2], true);
        let o = PathOptions::default();
        let d = distances(&g, None, o).unwrap();
        assert_eq!(d.row(0), &[0.0, 1.0, 2.0, f64::INFINITY]);
        assert_eq!(
            d.row(2),
            &[f64::INFINITY, f64::INFINITY, 0.0, f64::INFINITY]
        );
        let din = distances(
            &g,
            Some(&[2]),
            PathOptions {
                mode: DegreeMode::In,
                ..o
            },
        )
        .unwrap();
        assert_eq!(din.row(0), &[2.0, 1.0, 0.0, f64::INFINITY]);
        // Pairs: (0,1) 1, (0,2) 2, (1,2) 1.
        assert!(close(mean_distance(&g, o, true).unwrap(), 4.0 / 3.0));
        assert_eq!(mean_distance(&g, o, false).unwrap(), f64::INFINITY);
        assert_eq!(diameter(&g, o, true).unwrap(), 2.0);
        assert_eq!(eccentricity(&g, o).unwrap(), vec![2.0, 1.0, 0.0, 0.0]);
        assert_eq!(radius(&g, o).unwrap(), 0.0);
        assert_eq!(
            reachability(&g, DegreeMode::Out, 0).unwrap(),
            vec![2.0, 1.0, 0.0, 0.0]
        );
        assert_eq!(
            reachability(&g, DegreeMode::In, 0).unwrap(),
            vec![0.0, 1.0, 2.0, 0.0]
        );
        assert_eq!(
            reachability(&g, DegreeMode::All, 0).unwrap(),
            vec![2.0, 2.0, 2.0, 0.0]
        );
        let one = graph(1, &[], &[], false);
        assert!(mean_distance(&one, o, true).unwrap().is_nan());
        assert_eq!(diameter(&one, o, true).unwrap(), 0.0);
        let none = graph(0, &[], &[], false);
        assert!(diameter(&none, o, true).unwrap().is_nan());
        assert!(radius(&none, o).unwrap().is_nan());
        assert!(distances(&g, Some(&[4]), o).is_err());
    }

    #[test]
    fn weighted_distances_use_every_algorithm() {
        let w = [2.0, 1.0, 1.0, 3.0, 2.0];
        let g = Graph::from_edges(5, &[0, 0, 1, 2, 3], &[1, 2, 3, 3, 4], Some(&w), true).unwrap();
        let run = |algorithm| {
            distances(
                &g,
                None,
                PathOptions {
                    weighted: true,
                    algorithm,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let a = run(PathAlgorithm::Dial);
        assert_eq!(a, run(PathAlgorithm::Dijkstra));
        assert_eq!(a.row(0), &[0.0, 2.0, 1.0, 3.0, 5.0]);
        // Uniform weights scale hop counts.
        let u = Graph::from_edges(3, &[0, 1], &[1, 2], Some(&[2.5, 2.5]), false).unwrap();
        let o = PathOptions {
            weighted: true,
            ..Default::default()
        };
        assert_eq!(distances(&u, None, o).unwrap().row(0), &[0.0, 2.5, 5.0]);
        assert!(close(mean_distance(&u, o, true).unwrap(), 2.5 * 8.0 / 6.0));
        assert_eq!(diameter(&u, o, true).unwrap(), 5.0);
        let s = PathOptions {
            weights_as: WeightKind::Strength,
            ..o
        };
        assert_eq!(diameter(&u, s, true).unwrap(), 0.8);
    }

    #[test]
    fn cut_points_and_bridges() {
        // Triangle 0-1-2, bridge 2-3, and 3-4 doubled (not a bridge), loop
        // at 4, isolate 5.
        let g = graph(6, &[0, 1, 2, 2, 3, 4, 4], &[1, 2, 0, 3, 4, 3, 4], false);
        assert_eq!(cut_points(&g), vec![2, 3]);
        assert_eq!(bridges(&g), vec![3]);
        // Mutual arcs are two edges in the undirected reading.
        let d = graph(3, &[0, 1, 1], &[1, 0, 2], true);
        assert_eq!(cut_points(&d), vec![1]);
        assert_eq!(bridges(&d), vec![2]);
        // A star: the hub is the only cut point, every edge a bridge.
        let star = graph(4, &[0, 0, 0], &[1, 2, 3], false);
        assert_eq!(cut_points(&star), vec![0]);
        assert_eq!(bridges(&star), vec![0, 1, 2]);
    }

    #[test]
    fn deep_path_does_not_overflow() {
        let n = 200_000;
        let s: Vec<i64> = (0..n as i64 - 1).collect();
        let t: Vec<i64> = (1..n as i64).collect();
        let g = graph(n, &s, &t, false);
        assert_eq!(cut_points(&g).len(), n - 2);
        assert_eq!(bridges(&g).len(), n - 1);
    }

    #[test]
    fn summary_counts() {
        let g = graph(5, &[0, 0, 1, 1, 2, 2], &[1, 1, 0, 1, 2, 2], true);
        let c = edge_counts(&g);
        assert_eq!(
            c,
            EdgeCounts {
                loops: 3,
                multiple: 2,
                isolates: 3
            }
        );
        let u = graph(5, &[0, 0, 1, 1, 2, 2], &[1, 1, 0, 1, 2, 2], false);
        assert_eq!(edge_counts(&u).multiple, 3);
    }
}
