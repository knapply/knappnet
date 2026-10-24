//! k-cores and coreness.
//!
//! The k-core of a graph (Seidman 1983, `seidman1983network`) is its largest
//! subgraph in which every node has degree at least k. The coreness (core
//! number) of a node is the largest k for which it belongs to the k-core.
//! Coreness is computed with the O(n + m) bucket algorithm of Batagelj and
//! Zaversnik (2003, `batagelj2003algorithm`): nodes are removed in order of
//! non-decreasing current degree, and each removal lowers the degree of the
//! node's remaining neighbors.
//!
//! Conventions (matching igraph's `coreness()`):
//! - Directed graphs take a mode: `All` counts in-arcs plus out-arcs (a
//!   mutual pair counts twice), `Out` builds cores on out-degree, and `In`
//!   on in-degree. Undirected graphs ignore the mode.
//! - Parallel edges each count, as in degree.
//! - A self-loop adds 2 to a node's degree in an undirected graph, and 1 to
//!   out-degree and 1 to in-degree in a directed graph. A loop stays inside
//!   every subgraph that contains its node, so it never stops counting.
//!   NetworkX's `core_number` rejects self-loops and multigraphs; on simple
//!   graphs without loops the two agree.
//! - Weights are ignored: cores are defined on binary ties.
//! - Isolated nodes have coreness 0; disconnected graphs need no special
//!   treatment.
//!
//! # Shells, crusts, coronas, trusses, and onion layers
//!
//! From the core numbers c (`DegreeMode::All`, NetworkX's `core_number`)
//! follow, as in NetworkX 3.7 (`hagberg2008exploring`): the k-core {v :
//! c_v ≥ k}, the k-shell {v : c_v = k} and the k-crust {v : c_v ≤ k}
//! (Carmi et al. 2007, `carmi2007model`), and the k-corona, the nodes of
//! the k-shell with exactly k neighbours in the k-core (Goltsev,
//! Dorogovtsev, and Mendes 2006, `goltsev2006k`). Without k, the core and
//! shell take the largest core number and the crust one less. NetworkX
//! counts corona neighbours over `G[v]`: the distinct neighbours, the
//! successors in a directed graph.
//!
//! The k-truss (Wang and Cheng 2012, `wang2012truss`) is the
//! largest subgraph in which every edge lies in at least k − 2 triangles
//! of the subgraph; the truss number of an edge is the largest k whose
//! k-truss contains it (2 for an edge in no triangle). [`truss_numbers`]
//! computes them all by peeling edges in order of support (triangle
//! count), in O(m^1.5), and [`k_truss`] keeps the edges with truss number
//! at least k and the nodes they touch, which is NetworkX's `k_truss`
//! (that removes isolated nodes, so for k ≤ 2 it keeps every edge and the
//! nodes with one).
//!
//! The onion decomposition (Hébert-Dufresne, Grochow, and Allard 2016,
//! `hebertdufresne2016multi`) refines the cores into layers: isolated
//! nodes form layer 1; then, round after round, the core value rises to
//! the smallest remaining degree if that is larger, and every node whose
//! degree is at most the core value forms the next layer and is removed
//! at once. [`onion_layers`] reproduces NetworkX's `onion_layers`.
//!
//! Conventions for trusses and onion layers: they are defined on the
//! simple undirected graph (direction ignored, parallel edges merged).
//! Self-loops are ignored (NetworkX raises on them); a loop has truss
//! number 0. Parallel copies of an edge share its truss number.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::degree::DegreeMode;

/// Coreness of every node, in node order. See the module documentation.
pub fn coreness(g: &Graph, mode: DegreeMode) -> Vec<usize> {
    let n = g.node_count();
    let mode = if g.is_directed() {
        mode
    } else {
        DegreeMode::All
    };
    let deg: Vec<usize> = (0..n)
        .map(|v| {
            if !g.is_directed() {
                // Undirected rows store a loop once; it contributes 2.
                let row = g.out().neighbors(v);
                row.len() + row.iter().filter(|&&u| u as usize == v).count()
            } else {
                match mode {
                    DegreeMode::Out => g.out().degree(v),
                    DegreeMode::In => g.inc().degree(v),
                    DegreeMode::All => g.out().degree(v) + g.inc().degree(v),
                }
            }
        })
        .collect();
    // Removing v lowers the degree of every node whose counted arcs reach v:
    // v's in-neighbors for out-cores, its out-neighbors for in-cores.
    peel(deg, |v, buf| {
        buf.clear();
        if !g.is_directed() {
            buf.extend_from_slice(g.out().neighbors(v));
            return;
        }
        if mode != DegreeMode::Out {
            buf.extend_from_slice(g.out().neighbors(v));
        }
        if mode != DegreeMode::In {
            buf.extend_from_slice(g.inc().neighbors(v));
        }
    })
    .0
}

/// The bucket algorithm of Batagelj and Zaversnik. `deg` holds initial
/// degrees; `affected(v, buf)` fills `buf` with the nodes whose degree drops
/// by one per entry when `v` is removed. Returns the core numbers and the
/// removal order, which is a smallest-last (degeneracy) ordering (Matula
/// and Beck 1983, `matula1983smallest`).
pub(crate) fn peel<F>(mut deg: Vec<usize>, affected: F) -> (Vec<usize>, Vec<u32>)
where
    F: Fn(usize, &mut Vec<u32>),
{
    let n = deg.len();
    let max_deg = deg.iter().copied().max().unwrap_or(0);
    // bin[d] = position in `vert` of the first node with current degree d.
    let mut bin = vec![0usize; max_deg + 1];
    for &d in &deg {
        bin[d] += 1;
    }
    let mut start = 0;
    for b in bin.iter_mut() {
        let count = *b;
        *b = start;
        start += count;
    }
    let mut pos = vec![0usize; n];
    let mut vert = vec![0u32; n];
    {
        let mut next = bin.clone();
        for v in 0..n {
            pos[v] = next[deg[v]];
            vert[pos[v]] = v as u32;
            next[deg[v]] += 1;
        }
    }
    let mut buf = Vec::new();
    for i in 0..n {
        let v = vert[i] as usize;
        affected(v, &mut buf);
        for &u in &buf {
            let u = u as usize;
            if deg[u] > deg[v] {
                // Move u to the front of its bin, then shrink the bin by one.
                let du = deg[u];
                let pu = pos[u];
                let pw = bin[du];
                let w = vert[pw] as usize;
                if u != w {
                    pos[u] = pw;
                    vert[pu] = w as u32;
                    pos[w] = pu;
                    vert[pw] = u as u32;
                }
                bin[du] += 1;
                deg[u] -= 1;
            }
        }
    }
    (deg, vert)
}

/// Which part of the core decomposition [`core_subset`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoreSubset {
    /// The k-core: core number at least k.
    Core,
    /// The k-shell: core number exactly k.
    Shell,
    /// The k-crust: core number at most k.
    Crust,
    /// The k-corona: the k-shell nodes with exactly k neighbours in the
    /// k-core.
    Corona,
}

impl CoreSubset {
    /// Decodes the C interface codes: 0 = core, 1 = shell, 2 = crust, 3 =
    /// corona.
    pub fn from_code(code: i32) -> Result<CoreSubset> {
        match code {
            0 => Ok(CoreSubset::Core),
            1 => Ok(CoreSubset::Shell),
            2 => Ok(CoreSubset::Crust),
            3 => Ok(CoreSubset::Corona),
            _ => Err(Error::InvalidArgument(format!(
                "unknown core subset code {code}"
            ))),
        }
    }
}

/// The nodes of the k-core, k-shell, k-crust, or k-corona, increasing. `k`
/// = `None` takes NetworkX's default (the largest core number, one less
/// for the crust); the corona needs a k. See the module documentation.
pub fn core_subset(g: &Graph, which: CoreSubset, k: Option<i64>) -> Result<Vec<u32>> {
    let c = coreness(g, DegreeMode::All);
    let max = c.iter().copied().max();
    let k = match (k, which) {
        (Some(k), _) => k,
        (None, CoreSubset::Corona) => {
            return Err(Error::InvalidArgument(
                "the k-corona needs a value of k".into(),
            ))
        }
        (None, CoreSubset::Crust) => match max {
            Some(m) => m as i64 - 1,
            None => return Ok(Vec::new()),
        },
        (None, _) => match max {
            Some(m) => m as i64,
            None => return Ok(Vec::new()),
        },
    };
    let n = g.node_count();
    let keep = |v: usize| {
        let cv = c[v] as i64;
        match which {
            CoreSubset::Core => cv >= k,
            CoreSubset::Shell => cv == k,
            CoreSubset::Crust => cv <= k,
            CoreSubset::Corona => {
                if cv != k {
                    return false;
                }
                let row = g.out().neighbors(v);
                let mut count = 0i64;
                for (i, &w) in row.iter().enumerate() {
                    if (i == 0 || row[i - 1] != w) && c[w as usize] as i64 >= k {
                        count += 1;
                    }
                }
                count == k
            }
        }
    };
    Ok((0..n).filter(|&v| keep(v)).map(|v| v as u32).collect())
}

/// Truss number of every edge of `g`, in edge order (0 for self-loops).
/// See the module documentation.
pub fn truss_numbers(g: &Graph) -> Vec<u32> {
    let n = g.node_count();
    // Distinct undirected pairs u < v and, per node, (neighbour, pair).
    let mut pairs: Vec<(u32, u32)> = g
        .edge_sources()
        .iter()
        .zip(g.edge_targets())
        .filter(|(s, t)| s != t)
        .map(|(&s, &t)| (s.min(t), s.max(t)))
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    let m = pairs.len();
    let mut rows: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    for (e, &(u, v)) in pairs.iter().enumerate() {
        rows[u as usize].push((v, e as u32));
        rows[v as usize].push((u, e as u32));
    }
    for r in rows.iter_mut() {
        r.sort_unstable();
    }
    let common = |e: usize, f: &mut dyn FnMut(u32, u32)| {
        let (u, v) = pairs[e];
        let (a, b) = (&rows[u as usize], &rows[v as usize]);
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].0.cmp(&b[j].0) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    f(a[i].1, b[j].1);
                    i += 1;
                    j += 1;
                }
            }
        }
    };
    let mut sup = vec![0usize; m];
    for (e, s) in sup.iter_mut().enumerate() {
        common(e, &mut |_, _| *s += 1);
    }
    // Bucket queue over supports, as in the core algorithm; an edge's
    // support counts only triangles whose other two edges remain.
    let max_sup = sup.iter().copied().max().unwrap_or(0);
    let mut bin = vec![0usize; max_sup + 1];
    for &d in &sup {
        bin[d] += 1;
    }
    let mut start = 0;
    for b in bin.iter_mut() {
        let count = *b;
        *b = start;
        start += count;
    }
    let mut pos = vec![0usize; m];
    let mut vert = vec![0u32; m];
    {
        let mut next = bin.clone();
        for e in 0..m {
            pos[e] = next[sup[e]];
            vert[pos[e]] = e as u32;
            next[sup[e]] += 1;
        }
    }
    let mut removed = vec![false; m];
    fn lower(
        f: usize,
        floor: usize,
        sup: &mut [usize],
        bin: &mut [usize],
        pos: &mut [usize],
        vert: &mut [u32],
    ) {
        if sup[f] > floor {
            let df = sup[f];
            let pf = pos[f];
            let pw = bin[df];
            let w = vert[pw] as usize;
            if f != w {
                pos[f] = pw;
                vert[pf] = w as u32;
                pos[w] = pf;
                vert[pw] = f as u32;
            }
            bin[df] += 1;
            sup[f] -= 1;
        }
    }
    for i in 0..m {
        let e = vert[i] as usize;
        let floor = sup[e];
        let mut hit: Vec<(u32, u32)> = Vec::new();
        common(e, &mut |a, b| hit.push((a, b)));
        for (a, b) in hit {
            let (a, b) = (a as usize, b as usize);
            if removed[a] || removed[b] {
                continue;
            }
            lower(a, floor, &mut sup, &mut bin, &mut pos, &mut vert);
            lower(b, floor, &mut sup, &mut bin, &mut pos, &mut vert);
        }
        removed[e] = true;
    }
    g.edge_sources()
        .iter()
        .zip(g.edge_targets())
        .map(|(&s, &t)| {
            if s == t {
                0
            } else {
                let e = pairs.binary_search(&(s.min(t), s.max(t))).expect("a pair");
                sup[e] as u32 + 2
            }
        })
        .collect()
}

/// The k-truss: the edges (indices into the graph's edges, increasing)
/// with truss number at least k, and the nodes they touch. See the module
/// documentation.
pub fn k_truss(g: &Graph, k: i64) -> (Vec<usize>, Vec<u32>) {
    let t = truss_numbers(g);
    let mut touched = vec![false; g.node_count()];
    let mut edges = Vec::new();
    for (e, &x) in t.iter().enumerate() {
        let s = g.edge_sources()[e];
        let d = g.edge_targets()[e];
        if s != d && x as i64 >= k.max(2) {
            edges.push(e);
            touched[s as usize] = true;
            touched[d as usize] = true;
        }
    }
    let nodes = (0..g.node_count())
        .filter(|&v| touched[v])
        .map(|v| v as u32)
        .collect();
    (edges, nodes)
}

/// The onion layer of every node (1 for the first). See the module
/// documentation.
pub fn onion_layers(g: &Graph) -> Vec<u32> {
    let s = g.simple_undirected(false);
    let n = s.node_count();
    let mut deg: Vec<usize> = (0..n).map(|v| s.degree(v)).collect();
    let mut layer = vec![0u32; n];
    let mut current = 1u32;
    let mut alive = n;
    for v in 0..n {
        if deg[v] == 0 {
            layer[v] = 1;
            alive -= 1;
        }
    }
    if alive < n {
        current = 2;
    }
    // Lazy buckets by degree for finding the smallest remaining degree.
    let max_deg = deg.iter().copied().max().unwrap_or(0);
    let mut bucket: Vec<Vec<u32>> = vec![Vec::new(); max_deg + 1];
    for v in 0..n {
        if layer[v] == 0 {
            bucket[deg[v]].push(v as u32);
        }
    }
    let mut core = 1usize;
    let mut queued = vec![false; n];
    let mut next: Vec<u32> = Vec::new();
    for v in 0..n {
        if layer[v] == 0 && deg[v] <= core {
            queued[v] = true;
            next.push(v as u32);
        }
    }
    let mut scan = core + 1;
    while alive > 0 {
        if next.is_empty() {
            // Every remaining degree exceeds the core value: raise it to
            // the smallest one.
            loop {
                let b = &mut bucket[scan];
                while let Some(&v) = b.last() {
                    let v = v as usize;
                    if layer[v] != 0 || deg[v] != scan {
                        b.pop();
                    } else {
                        break;
                    }
                }
                if !b.is_empty() {
                    break;
                }
                scan += 1;
            }
            core = scan;
            for &v in &bucket[scan] {
                let v = v as usize;
                if layer[v] == 0 && deg[v] == scan && !queued[v] {
                    queued[v] = true;
                    next.push(v as u32);
                }
            }
        }
        let this = std::mem::take(&mut next);
        for &v in &this {
            layer[v as usize] = current;
        }
        alive -= this.len();
        for &v in &this {
            for &w in s.neighbors(v as usize) {
                let w = w as usize;
                if layer[w] != 0 {
                    continue;
                }
                deg[w] -= 1;
                if deg[w] <= core {
                    if !queued[w] {
                        queued[w] = true;
                        next.push(w as u32);
                    }
                } else {
                    bucket[deg[w]].push(w as u32);
                }
            }
        }
        current += 1;
    }
    layer
}

#[cfg(test)]
mod tests {
    use super::*;

    fn und(n: usize, e: &[(i64, i64)]) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, false).unwrap()
    }

    #[test]
    fn hand_checked_undirected() {
        // A triangle 0-1-2 with a pendant 3 on node 2 and an isolate 4.
        let g = und(5, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 2, 2, 1, 0]);
        // The NetworkX documentation example (Havel-Hakimi graph of
        // degrees [0, 1, 2, 2, 2, 2, 3]) has core numbers
        // {0: 1, 1: 2, 2: 2, 3: 2, 4: 1, 5: 2, 6: 0}; a K4 is a 3-core.
        let k4 = und(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]);
        assert_eq!(coreness(&k4, DegreeMode::All), vec![3; 4]);
    }

    #[test]
    fn loops_and_parallel_edges_count() {
        // A loop alone gives degree 2 and keeps it.
        let g = und(2, &[(0, 0), (0, 1)]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 1]);
        // A doubled edge is a 2-core.
        let g = und(2, &[(0, 1), (0, 1)]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 2]);
    }

    #[test]
    fn directed_modes() {
        // 0 -> 1, 0 -> 2, 1 -> 2, 2 -> 0
        let g = Graph::from_edges(3, &[0, 0, 1, 2], &[1, 2, 2, 0], None, true).unwrap();
        // Out-degrees 2, 1, 1. Removing 1 drops nothing below 1; the whole
        // graph is a 1-out-core.
        assert_eq!(coreness(&g, DegreeMode::Out), vec![1, 1, 1]);
        // In-degrees 1, 1, 2.
        assert_eq!(coreness(&g, DegreeMode::In), vec![1, 1, 1]);
        // Total degrees 3, 2, 3: removing 1 leaves 0 and 2 joined by two
        // arcs, a 2-core.
        assert_eq!(coreness(&g, DegreeMode::All), vec![2, 2, 2]);
    }

    #[test]
    fn trusses_and_onion_layers() {
        // K4 on 0..4 with a pendant path 3-4-5: K4 edges have truss 4,
        // the path edges 2.
        let g = und(
            6,
            &[
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 2),
                (1, 3),
                (2, 3),
                (3, 4),
                (4, 5),
            ],
        );
        assert_eq!(truss_numbers(&g), vec![4, 4, 4, 4, 4, 4, 2, 2]);
        let (e, v) = k_truss(&g, 3);
        assert_eq!(e, vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(v, vec![0, 1, 2, 3]);
        // Onion: 5 (degree 1) first, then 4, then the K4.
        assert_eq!(onion_layers(&g), vec![3, 3, 3, 3, 2, 1]);
        // Shells: core numbers 3, 3, 3, 3, 1, 1.
        assert_eq!(
            core_subset(&g, CoreSubset::Shell, None).unwrap(),
            vec![0, 1, 2, 3]
        );
        assert_eq!(
            core_subset(&g, CoreSubset::Crust, None).unwrap(),
            vec![4, 5]
        );
        assert_eq!(
            core_subset(&g, CoreSubset::Corona, Some(1)).unwrap(),
            vec![5]
        );
    }

    #[test]
    fn empty_graph() {
        let g = und(0, &[]);
        assert!(coreness(&g, DegreeMode::All).is_empty());
        let g = und(3, &[]);
        assert_eq!(coreness(&g, DegreeMode::All), vec![0, 0, 0]);
    }
}
