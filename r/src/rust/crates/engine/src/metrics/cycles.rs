//! Cycles: cycle bases, minimum cycle bases, simple cycles, chordless
//! cycles, girth, and finding a cycle.
//!
//! Sources:
//! - A cycle basis of an undirected network is a set of cycles from which
//!   every cycle (every even subgraph) is a sum modulo 2; it has m - n + c
//!   cycles (Kirchhoff 1847, `kirchhoff1847ueber`). The basis returned by
//!   [`cycle_basis`] closes each non-tree edge of a depth-first spanning
//!   forest with tree paths (Paton 1969, `paton1969algorithm`), with the
//!   bookkeeping of NetworkX's `cycle_basis` (`hagberg2008exploring`), so
//!   it is NetworkX's basis, cycle for cycle, for the same neighbour order.
//! - A minimum cycle basis has least total weight. de Pina's (1995,
//!   `depina1995applications`; Kavitha et al. 2009, `kavitha2009cycle`)
//!   algorithm finds one cycle per non-tree edge: keeping a witness vector
//!   S_i for each, the i-th cycle is a shortest cycle meeting S_i in an odd
//!   number of edges, found as a shortest path between the two copies of a
//!   node in the graph doubled over S_i, and later witnesses are updated
//!   to stay orthogonal to the cycles found. Such a cycle contains an edge
//!   of S_i, so the searches start only from the lower endpoint of each
//!   edge of S_i, and run in parallel.
//! - Simple cycles by Johnson's (1975, `johnson1975finding`) algorithm,
//!   which spends O(n + m) per cycle: from each start node, in increasing
//!   order, cycles are searched in its strongly connected component among
//!   the nodes not yet used as starts, with blocking to skip dead ends.
//!   With a length bound, a bounded depth-first search from each start is
//!   used instead (correct for any bound; Gupta and Suzumura 2021,
//!   `gupta2021finding`, discuss the bounded problem).
//! - A chordless (induced) cycle is one whose nodes induce no edge besides
//!   the cycle's own. They are found by extending chordless paths from
//!   each start node with counters of path neighbours (Dias et al. 2013,
//!   `dias2013efficient`, as NetworkX does), so no non-chordless cycle is
//!   ever built.
//! - The girth is the length of a shortest cycle; breadth-first search
//!   from every node finds it (Itai and Rodeh 1978, `itai1978finding`).
//! - [`find_cycle`] follows NetworkX's `find_cycle`: a depth-first search
//!   over edges (NetworkX's `edge_dfs`) from each node in turn, reporting
//!   the first cycle it closes.
//!
//! Conventions:
//! - Cycles are node sequences; a self-loop is the cycle `[v]` and, in an
//!   undirected network, parallel edges between u and v give the cycle
//!   `[u, v]`. Simple and chordless cycles treat multigraphs as NetworkX
//!   does: one cycle per node sequence (and per direction for directed
//!   networks), a self-loop counted once per node; for chordless cycles a
//!   node with a self-loop lies on no longer cycle, two parallel edges are
//!   a chordless 2-cycle and more are not, and parallel arcs or a pair of
//!   opposite arcs are chords of any longer cycle. Where NetworkX's search
//!   would let parallel edges act as an unseen chord, this module follows
//!   the definition (the induced subgraph has exactly as many edges as the
//!   cycle has nodes).
//! - Cycle bases and minimum cycle bases are for undirected networks and
//!   use every edge (loops and parallel edges included), so the basis has
//!   m - n + c cycles; NetworkX refuses multigraphs for these, and on a
//!   simple graph the two coincide. [`cycle_basis`] lists the extra cycles
//!   of parallel edges and repeated loops after NetworkX's.
//! - Minimum cycle bases use weights (with `weighted`, non-negative) or 1.
//!   The minimum total weight, and the sorted list of cycle weights, are
//!   the same for every minimum cycle basis; the cycles themselves need not
//!   be unique.
//! - Girth counts edges, 1 for a self-loop and 2 for parallel edges, and is
//!   infinite for a forest. Directed networks are accepted (the shortest
//!   directed cycle); NetworkX accepts only undirected simple graphs.
//! - The order of simple and chordless cycles, and where each cycle
//!   starts, is the engine's (NetworkX's depends on Python set order);
//!   compare them as sets.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::order::{rows, Direction, NeighborOrder};
use crate::parallel;

fn undirected(g: &Graph, what: &str) -> Result<()> {
    if g.is_directed() {
        Err(Error::InvalidArgument(format!(
            "{what} is defined for undirected networks"
        )))
    } else {
        Ok(())
    }
}

/// Number of self-loops at each node and, for an undirected network, the
/// number of edges between each pair (as a sorted list of (u, v, count)
/// with u < v).
fn loops_and_multiplicity(g: &Graph) -> (Vec<u32>, Vec<(u32, u32, u32)>) {
    let n = g.node_count();
    let mut loops = vec![0u32; n];
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    for (&a, &b) in g.edge_sources().iter().zip(g.edge_targets()) {
        if a == b {
            loops[a as usize] += 1;
        } else if g.is_directed() {
            pairs.push((a, b));
        } else {
            pairs.push((a.min(b), a.max(b)));
        }
    }
    pairs.sort_unstable();
    let mut mult: Vec<(u32, u32, u32)> = Vec::new();
    for p in pairs {
        match mult.last_mut() {
            Some(last) if (last.0, last.1) == p => last.2 += 1,
            _ => mult.push((p.0, p.1, 1)),
        }
    }
    (loops, mult)
}

/// A cycle basis of an undirected network. With `root`, the first search
/// starts there; otherwise searches start, as in NetworkX, from the
/// highest-numbered node not yet reached. See the module documentation.
pub fn cycle_basis(g: &Graph, order: NeighborOrder, root: Option<usize>) -> Result<Vec<Vec<u32>>> {
    undirected(g, "a cycle basis")?;
    let n = g.node_count();
    if let Some(r) = root {
        if r >= n {
            return Err(Error::InvalidArgument(format!(
                "root {r} is not a node (n = {n})"
            )));
        }
    }
    let adj = rows(g, Direction::Out, order, false, true);
    let mut remaining = vec![true; n];
    let mut last = n; // search position for the highest remaining node
    let mut cycles: Vec<Vec<u32>> = Vec::new();
    const NONE: usize = usize::MAX;
    let mut pred = vec![NONE; n];
    // used[v]: the nodes from which v's cycles have already been closed
    // (NetworkX's `used`), kept as small sorted lists.
    let mut used: Vec<Option<Vec<u32>>> = vec![None; n];
    let mut next_root = root;
    loop {
        let r = match next_root.take() {
            Some(r) => r,
            None => {
                while last > 0 && !remaining[last - 1] {
                    last -= 1;
                }
                if last == 0 {
                    break;
                }
                last - 1
            }
        };
        let mut component = vec![r];
        pred[r] = r;
        used[r] = Some(Vec::new());
        let mut stack = vec![r];
        while let Some(z) = stack.pop() {
            for &nbr in adj.neighbors(z) {
                let nbr = nbr as usize;
                if used[nbr].is_none() {
                    pred[nbr] = z;
                    stack.push(nbr);
                    used[nbr] = Some(vec![z as u32]);
                    component.push(nbr);
                } else if nbr == z {
                    cycles.push(vec![z as u32]);
                } else if !used[z].as_ref().expect("visited").contains(&(nbr as u32)) {
                    let pn = used[nbr].as_ref().expect("visited");
                    let mut cycle = vec![nbr as u32, z as u32];
                    let mut p = pred[z];
                    while !pn.contains(&(p as u32)) {
                        cycle.push(p as u32);
                        p = pred[p];
                    }
                    cycle.push(p as u32);
                    cycles.push(cycle);
                    used[nbr].as_mut().expect("visited").push(z as u32);
                }
            }
        }
        for &v in &component {
            remaining[v] = false;
        }
    }
    // Cycles NetworkX's simple graph cannot hold: extra parallel edges and
    // repeated loops, in node order.
    let (loops, mult) = loops_and_multiplicity(g);
    let mut extra: Vec<(u32, Vec<u32>)> = Vec::new();
    for (v, &k) in loops.iter().enumerate() {
        for _ in 1..k {
            extra.push((v as u32, vec![v as u32]));
        }
    }
    for &(u, v, k) in &mult {
        for _ in 1..k {
            extra.push((u, vec![u, v]));
        }
    }
    extra.sort_by_key(|e| e.0);
    cycles.extend(extra.into_iter().map(|e| e.1));
    Ok(cycles)
}

/// A minimum cycle basis: its cycles (node sequences) and their weights.
#[derive(Debug, Clone, PartialEq)]
pub struct CycleBasis {
    /// The cycles.
    pub cycles: Vec<Vec<u32>>,
    /// The weight of each cycle.
    pub weights: Vec<f64>,
}

/// Turns a set of edges forming one cycle into its node sequence, starting
/// at its lowest node and continuing toward the lower of its two cycle
/// neighbours.
fn cycle_nodes(g: &Graph, edges: &[usize]) -> Vec<u32> {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    if edges.len() == 1 {
        return vec![src[edges[0]]];
    }
    // Incidences of each node on the cycle.
    let mut inc: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for &e in edges {
        inc.entry(src[e]).or_default().push(e);
        inc.entry(tgt[e]).or_default().push(e);
    }
    let start = *inc.keys().next().expect("non-empty cycle");
    let other = |e: usize, v: u32| if src[e] == v { tgt[e] } else { src[e] };
    let first = inc[&start]
        .iter()
        .copied()
        .min_by_key(|&e| (other(e, start), e))
        .expect("degree two");
    let mut out = vec![start];
    let (mut v, mut via) = (other(first, start), first);
    while v != start {
        out.push(v);
        let next = inc[&v]
            .iter()
            .copied()
            .find(|&e| e != via)
            .expect("degree two");
        via = next;
        v = other(next, v);
    }
    out
}

/// Splits an even subgraph (every node of even degree) into simple cycles,
/// returned as edge lists.
fn split_cycles(g: &Graph, edges: &[usize]) -> Vec<Vec<usize>> {
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    let mut inc: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for &e in edges {
        inc.entry(src[e]).or_default().push(e);
        if src[e] != tgt[e] {
            inc.entry(tgt[e]).or_default().push(e);
        }
    }
    let mut used: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut cycles = Vec::new();
    for &e0 in edges {
        if used.contains(&e0) {
            continue;
        }
        if src[e0] == tgt[e0] {
            used.insert(e0);
            cycles.push(vec![e0]);
            continue;
        }
        // Walk from src[e0] along unused edges until a node repeats.
        let mut walk_nodes = vec![src[e0]];
        let mut walk_edges: Vec<usize> = Vec::new();
        let mut v = src[e0];
        let mut e = e0;
        loop {
            used.insert(e);
            walk_edges.push(e);
            let w = if src[e] == v { tgt[e] } else { src[e] };
            if let Some(pos) = walk_nodes.iter().position(|&x| x == w) {
                cycles.push(walk_edges.split_off(pos));
                walk_nodes.truncate(pos + 1);
                if walk_edges.is_empty() {
                    break;
                }
                v = w;
            } else {
                walk_nodes.push(w);
                v = w;
            }
            match inc[&v].iter().copied().find(|x| !used.contains(x)) {
                Some(x) => e = x,
                None => break,
            }
        }
    }
    cycles
}

/// A minimum cycle basis of an undirected network (de Pina's algorithm).
/// See the module documentation.
pub fn minimum_cycle_basis(g: &Graph, weighted: bool, threads: usize) -> Result<CycleBasis> {
    undirected(g, "a minimum cycle basis")?;
    let n = g.node_count();
    let m = g.edge_count();
    let w: Vec<f64> = if weighted {
        let w = g.edge_weights().ok_or(Error::MissingWeights)?;
        if let Some(e) = w.iter().position(|&x| x < 0.0) {
            return Err(Error::InvalidArgument(format!(
                "edge {e} has negative weight {}; a minimum cycle basis needs non-negative \
                 weights",
                w[e]
            )));
        }
        w.to_vec()
    } else {
        vec![1.0; m]
    };
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    // Spanning forest in edge order; the other edges are chords.
    let mut forest = crate::union_find::UnionFind::new(n);
    let chords: Vec<usize> = (0..m)
        .filter(|&e| forest.union(src[e] as usize, tgt[e] as usize).is_none())
        .collect();
    let words = (m + 63) / 64;
    let bit = |set: &[u64], e: usize| set[e / 64] >> (e % 64) & 1 == 1;
    let mut witness: Vec<Vec<u64>> = chords
        .iter()
        .map(|&e| {
            let mut s = vec![0u64; words];
            s[e / 64] |= 1 << (e % 64);
            s
        })
        .collect();
    // Incidence lists: (neighbour, edge).
    let mut adj: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    for e in 0..m {
        adj[src[e] as usize].push((tgt[e], e as u32));
        if src[e] != tgt[e] {
            adj[tgt[e] as usize].push((src[e], e as u32));
        }
    }
    let mut cycles = Vec::with_capacity(chords.len());
    let mut weights = Vec::with_capacity(chords.len());
    for i in 0..chords.len() {
        let s = witness[i].clone();
        // A cycle meeting S_i oddly contains an edge of S_i, and so its
        // lower endpoint: shortest paths from (v, 0) to (v, 1) are needed
        // only from those nodes. They run in parallel; the best
        // (length, v) wins.
        let mut starts: Vec<usize> = (0..m)
            .filter(|&e| bit(&s, e))
            .map(|e| src[e].min(tgt[e]) as usize)
            .collect();
        starts.sort_unstable();
        starts.dedup();
        let parts = parallel::for_each_item(
            starts.len(),
            threads,
            || (f64::INFINITY, usize::MAX, Vec::<usize>::new()),
            |best, k| {
                let v = starts[k];
                if let Some((d, path)) = lifted_path(&adj, &w, &s, v, best.0) {
                    if d < best.0 || (d == best.0 && v < best.1) {
                        *best = (d, v, path);
                    }
                }
            },
        );
        let best = parts
            .into_iter()
            .filter(|p| p.1 != usize::MAX)
            .min_by(|a, b| a.0.partial_cmp(&b.0).expect("finite").then(a.1.cmp(&b.1)))
            .expect("every chord closes a cycle");
        // Edges used an odd number of times form an even subgraph; one of
        // its cycles meets S_i oddly and is no longer than the path.
        let mut count: std::collections::BTreeMap<usize, u32> = std::collections::BTreeMap::new();
        for e in best.2 {
            *count.entry(e).or_default() += 1;
        }
        let odd: Vec<usize> = count
            .into_iter()
            .filter(|&(_, c)| c % 2 == 1)
            .map(|(e, _)| e)
            .collect();
        let cycle = split_cycles(g, &odd)
            .into_iter()
            .find(|c| c.iter().filter(|&&e| bit(&s, e)).count() % 2 == 1)
            .expect("an odd cycle");
        let mut cvec = vec![0u64; words];
        for &e in &cycle {
            cvec[e / 64] |= 1 << (e % 64);
        }
        for wj in witness.iter_mut().skip(i + 1) {
            let parity = cvec
                .iter()
                .zip(wj.iter())
                .map(|(a, b)| (a & b).count_ones())
                .sum::<u32>();
            if parity % 2 == 1 {
                for (x, y) in wj.iter_mut().zip(&s) {
                    *x ^= y;
                }
            }
        }
        weights.push(cycle.iter().map(|&e| w[e]).sum());
        cycles.push(cycle_nodes(g, &cycle));
    }
    Ok(CycleBasis { cycles, weights })
}

/// Dijkstra in the graph doubled over the witness `s` (an edge in `s`
/// switches copies), from (v, 0) to (v, 1); the edges of a shortest path
/// and its length, if shorter than `bound` (ties at `bound` too).
fn lifted_path(
    adj: &[Vec<(u32, u32)>],
    w: &[f64],
    s: &[u64],
    v: usize,
    bound: f64,
) -> Option<(f64, Vec<usize>)> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let n = adj.len();
    let mut dist = vec![f64::INFINITY; 2 * n];
    let mut via = vec![(usize::MAX, usize::MAX); 2 * n];
    #[derive(PartialEq)]
    struct Item(f64, usize);
    impl Eq for Item {}
    impl PartialOrd for Item {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for Item {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            self.0
                .partial_cmp(&o.0)
                .expect("finite")
                .then(self.1.cmp(&o.1))
        }
    }
    let mut heap = BinaryHeap::new();
    dist[v] = 0.0;
    heap.push(Reverse(Item(0.0, v)));
    let target = n + v;
    while let Some(Reverse(Item(d, x))) = heap.pop() {
        if d > dist[x] || d > bound {
            continue;
        }
        if x == target {
            break;
        }
        let (node, layer) = (x % n, x / n);
        for &(y, e) in &adj[node] {
            let e = e as usize;
            let flip = s[e / 64] >> (e % 64) & 1 == 1;
            let ly = if flip { 1 - layer } else { layer };
            let yx = ly * n + y as usize;
            let nd = d + w[e];
            if nd < dist[yx] {
                dist[yx] = nd;
                via[yx] = (x, e);
                heap.push(Reverse(Item(nd, yx)));
            }
        }
    }
    if dist[target] > bound || dist[target].is_infinite() {
        return None;
    }
    let mut path = Vec::new();
    let mut x = target;
    while x != v {
        let (p, e) = via[x];
        path.push(e);
        x = p;
    }
    Some((dist[target], path))
}

/// Simple adjacency without loops, deduplicated, in sorted order: out- and
/// in-neighbours for a directed network, neighbours for an undirected one.
fn simple_adjacency(g: &Graph) -> (Vec<Vec<u32>>, Vec<Vec<u32>>) {
    let out = rows(g, Direction::Out, NeighborOrder::Sorted, true, true);
    let inc = rows(g, Direction::In, NeighborOrder::Sorted, true, true);
    let n = g.node_count();
    let o: Vec<Vec<u32>> = (0..n).map(|v| out.neighbors(v).to_vec()).collect();
    let i: Vec<Vec<u32>> = (0..n).map(|v| inc.neighbors(v).to_vec()).collect();
    (o, i)
}

/// Members of the strongly connected component of `s` among the nodes
/// `>= s` (forward and backward reachability).
fn component_from(out: &[Vec<u32>], inc: &[Vec<u32>], s: usize) -> Vec<bool> {
    let n = out.len();
    let reach = |adj: &[Vec<u32>]| {
        let mut seen = vec![false; n];
        seen[s] = true;
        let mut stack = vec![s];
        while let Some(v) = stack.pop() {
            for &w in &adj[v] {
                let w = w as usize;
                if w >= s && !seen[w] {
                    seen[w] = true;
                    stack.push(w);
                }
            }
        }
        seen
    };
    let f = reach(out);
    let b = reach(inc);
    f.iter().zip(&b).map(|(&x, &y)| x && y).collect()
}

/// Every simple cycle (at most `length_bound` nodes long, and at most
/// `limit` cycles, when given). See the module documentation.
pub fn simple_cycles(
    g: &Graph,
    length_bound: Option<usize>,
    limit: Option<usize>,
) -> Result<Vec<Vec<u32>>> {
    let n = g.node_count();
    let directed = g.is_directed();
    let mut cycles: Vec<Vec<u32>> = Vec::new();
    let full = |c: &Vec<Vec<u32>>| matches!(limit, Some(l) if c.len() >= l);
    if length_bound == Some(0) {
        return Ok(cycles);
    }
    let (loops, mult) = loops_and_multiplicity(g);
    for (v, &k) in loops.iter().enumerate() {
        if k > 0 && !full(&cycles) {
            cycles.push(vec![v as u32]);
        }
    }
    if length_bound == Some(1) {
        return Ok(cycles);
    }
    if !directed {
        for &(u, v, k) in &mult {
            if k > 1 && !full(&cycles) {
                cycles.push(vec![u, v]);
            }
        }
    }
    let (out, inc) = simple_adjacency(g);
    let min_len = if directed { 2 } else { 3 };
    for s in 0..n {
        if full(&cycles) {
            break;
        }
        let comp = component_from(&out, &inc, s);
        if comp.iter().filter(|&&x| x).count() < 2 {
            continue;
        }
        let keep = |path: &[u32]| -> bool {
            path.len() >= min_len && (directed || path[1] < path[path.len() - 1])
        };
        match length_bound {
            Some(bound) => {
                // Bounded depth-first search over the component.
                let mut on_path = vec![false; n];
                on_path[s] = true;
                let mut path = vec![s as u32];
                let mut frames: Vec<usize> = vec![0];
                while let Some(&i) = frames.last() {
                    let v = *path.last().expect("path") as usize;
                    if i < out[v].len() {
                        *frames.last_mut().expect("frame") += 1;
                        let w = out[v][i] as usize;
                        if w == s {
                            if keep(&path) && !full(&cycles) {
                                cycles.push(path.clone());
                            }
                        } else if comp[w] && !on_path[w] && path.len() < bound {
                            on_path[w] = true;
                            path.push(w as u32);
                            frames.push(0);
                        }
                    } else {
                        frames.pop();
                        let v = path.pop().expect("path") as usize;
                        on_path[v] = false;
                    }
                }
            }
            None => johnson_from(s, &out, &comp, &mut cycles, &keep, limit),
        }
    }
    Ok(cycles)
}

/// Johnson's circuit search from `s` within `comp`.
fn johnson_from(
    s: usize,
    out: &[Vec<u32>],
    comp: &[bool],
    cycles: &mut Vec<Vec<u32>>,
    keep: &dyn Fn(&[u32]) -> bool,
    limit: Option<usize>,
) {
    let n = out.len();
    let mut blocked = vec![false; n];
    let mut b: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut path: Vec<u32> = vec![s as u32];
    // Frames: (node, next arc index, found a cycle below).
    let mut frames: Vec<(usize, usize, bool)> = vec![(s, 0, false)];
    blocked[s] = true;
    let unblock = |u: usize, blocked: &mut Vec<bool>, b: &mut Vec<Vec<u32>>| {
        let mut stack = vec![u];
        while let Some(x) = stack.pop() {
            if !blocked[x] {
                continue;
            }
            blocked[x] = false;
            for w in std::mem::take(&mut b[x]) {
                if blocked[w as usize] {
                    stack.push(w as usize);
                }
            }
        }
    };
    while let Some(frame) = frames.last_mut() {
        let (v, i) = (frame.0, frame.1);
        if i < out[v].len() {
            frame.1 += 1;
            let w = out[v][i] as usize;
            if !comp[w] {
                continue;
            }
            if w == s {
                if keep(&path) {
                    if matches!(limit, Some(l) if cycles.len() >= l) {
                        return;
                    }
                    cycles.push(path.clone());
                }
                frame.2 = true;
            } else if !blocked[w] {
                blocked[w] = true;
                path.push(w as u32);
                frames.push((w, 0, false));
            }
            continue;
        }
        let (v, _, found) = frames.pop().expect("frame");
        path.pop();
        if found {
            unblock(v, &mut blocked, &mut b);
        } else {
            for &w in &out[v] {
                let w = w as usize;
                if comp[w] && !b[w].contains(&(v as u32)) {
                    b[w].push(v as u32);
                }
            }
        }
        if let Some(parent) = frames.last_mut() {
            parent.2 |= found;
        }
    }
}

/// Every chordless cycle (at most `length_bound` nodes long, and at most
/// `limit` cycles, when given). See the module documentation.
pub fn chordless_cycles(
    g: &Graph,
    length_bound: Option<usize>,
    limit: Option<usize>,
) -> Result<Vec<Vec<u32>>> {
    let n = g.node_count();
    let directed = g.is_directed();
    let mut cycles: Vec<Vec<u32>> = Vec::new();
    let full = |c: &Vec<Vec<u32>>| matches!(limit, Some(l) if c.len() >= l);
    if length_bound == Some(0) {
        return Ok(cycles);
    }
    let (loops, mult) = loops_and_multiplicity(g);
    for (v, &k) in loops.iter().enumerate() {
        if k == 1 && !full(&cycles) {
            cycles.push(vec![v as u32]);
        }
    }
    if length_bound == Some(1) {
        return Ok(cycles);
    }
    let alive: Vec<bool> = loops.iter().map(|&k| k == 0).collect();
    // B: the simple undirected adjacency among loop-free nodes; F: the
    // arcs (edges) that may lie on a chordless cycle of length 3 or more.
    let mut b: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut f: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut arc_count: std::collections::HashMap<(u32, u32), u32> =
        std::collections::HashMap::new();
    for &(u, v, k) in &mult {
        if alive[u as usize] && alive[v as usize] {
            arc_count.insert((u, v), k);
        }
    }
    let has = |a: u32, c: u32| arc_count.get(&(a, c)).copied().unwrap_or(0);
    for &(u, v) in arc_count.keys() {
        b[u as usize].push(v);
        b[v as usize].push(u);
    }
    for row in b.iter_mut() {
        row.sort_unstable();
        row.dedup();
    }
    let mut digons: Vec<(u32, u32)> = Vec::new();
    if directed {
        for (&(u, v), &k) in &arc_count {
            let back = has(v, u);
            if k == 1 && back == 1 {
                if u < v {
                    digons.push((u, v));
                }
            } else if k == 1 && back == 0 {
                f[u as usize].push(v);
            }
        }
    } else {
        for (&(u, v), &k) in &arc_count {
            if k == 2 {
                digons.push((u, v));
            } else if k == 1 {
                f[u as usize].push(v);
                f[v as usize].push(u);
            }
        }
    }
    digons.sort_unstable();
    for (u, v) in digons {
        if !full(&cycles) {
            cycles.push(vec![u, v]);
        }
    }
    if length_bound == Some(2) {
        return Ok(cycles);
    }
    for row in f.iter_mut() {
        row.sort_unstable();
    }
    let mut fin: Vec<Vec<u32>> = vec![Vec::new(); n];
    if directed {
        for (u, row) in f.iter().enumerate() {
            for &v in row {
                fin[v as usize].push(u as u32);
            }
        }
    } else {
        fin = f.clone();
    }
    let adjacent_b = |a: u32, c: u32| b[a as usize].binary_search(&c).is_ok();
    let in_f = |a: u32, c: u32| f[a as usize].binary_search(&c).is_ok();
    let mut blocked = vec![0u32; n];
    for s in 0..n {
        if !alive[s] {
            continue;
        }
        let s32 = s as u32;
        let preds: Vec<u32> = fin[s].iter().copied().filter(|&u| u > s32).collect();
        let succs: Vec<u32> = f[s].iter().copied().filter(|&w| w > s32).collect();
        for &u in &preds {
            for &w in &succs {
                if u == w || (!directed && u >= w) {
                    continue;
                }
                if full(&cycles) {
                    return Ok(cycles);
                }
                // A stem u - s - w. An arc u -> w is a chord of every cycle
                // through the stem; a lone arc (edge) w -> u closes a
                // triangle; any other tie between u and w is a chord.
                if directed && arc_present(&arc_count, u, w) {
                    continue;
                }
                if in_f(w, u) {
                    cycles.push(vec![u, s32, w]);
                    continue;
                }
                if adjacent_b(u, w) || matches!(length_bound, Some(l) if l < 4) {
                    continue;
                }
                chordless_search(
                    [u, s32, w],
                    s32,
                    &f,
                    &b,
                    &mut blocked,
                    length_bound,
                    limit,
                    &mut cycles,
                );
            }
        }
    }
    Ok(cycles)
}

fn arc_present(counts: &std::collections::HashMap<(u32, u32), u32>, a: u32, c: u32) -> bool {
    counts.contains_key(&(a, c))
}

/// Extends the chordless path `stem` = [u, s, w] to every chordless cycle
/// closing back at u, using nodes above `s`.
#[allow(clippy::too_many_arguments)]
fn chordless_search(
    stem: [u32; 3],
    s: u32,
    f: &[Vec<u32>],
    b: &[Vec<u32>],
    blocked: &mut [u32],
    length_bound: Option<usize>,
    limit: Option<usize>,
    cycles: &mut Vec<Vec<u32>>,
) {
    let target = stem[0];
    let mut path: Vec<u32> = stem.to_vec();
    // blocked[x] counts path nodes after the target adjacent to x (plus
    // one for s itself, as a marker).
    blocked[s as usize] += 1;
    for &p in &path[1..] {
        for &x in &b[p as usize] {
            blocked[x as usize] += 1;
        }
    }
    let in_f = |a: u32, c: u32| f[a as usize].binary_search(&c).is_ok();
    let adjacent_b = |a: u32, c: u32| b[a as usize].binary_search(&c).is_ok();
    let mut frames: Vec<usize> = vec![0];
    while let Some(&i) = frames.last() {
        let y = *path.last().expect("path");
        let row = &f[y as usize];
        if i < row.len() {
            *frames.last_mut().expect("frame") += 1;
            let x = row[i];
            if x <= s || blocked[x as usize] != 1 {
                continue;
            }
            if matches!(length_bound, Some(l) if path.len() >= l) {
                continue;
            }
            if in_f(x, target) {
                if matches!(limit, Some(l) if cycles.len() >= l) {
                    break;
                }
                let mut c = path.clone();
                c.push(x);
                cycles.push(c);
            } else if adjacent_b(x, target) {
                continue;
            } else {
                for &z in &b[x as usize] {
                    blocked[z as usize] += 1;
                }
                path.push(x);
                frames.push(0);
            }
            continue;
        }
        frames.pop();
        if frames.is_empty() {
            break;
        }
        let x = path.pop().expect("path");
        for &z in &b[x as usize] {
            blocked[z as usize] -= 1;
        }
    }
    // Undo the counts of any extension left after an early stop, then the
    // stem's.
    while path.len() > 3 {
        let x = path.pop().expect("path");
        for &z in &b[x as usize] {
            blocked[z as usize] -= 1;
        }
    }
    for &p in &stem[1..] {
        for &x in &b[p as usize] {
            blocked[x as usize] -= 1;
        }
    }
    blocked[s as usize] -= 1;
}

/// Length of a shortest cycle (infinite without one). See the module
/// documentation.
pub fn girth(g: &Graph) -> f64 {
    let n = g.node_count();
    let (loops, mult) = loops_and_multiplicity(g);
    if loops.iter().any(|&k| k > 0) {
        return 1.0;
    }
    if !g.is_directed() && mult.iter().any(|&(_, _, k)| k > 1) {
        return 2.0;
    }
    let (out, _) = simple_adjacency(g);
    let mut best = usize::MAX;
    let mut dist = vec![usize::MAX; n];
    let mut parent = vec![usize::MAX; n];
    let mut queue = Vec::with_capacity(n);
    for s in 0..n {
        queue.clear();
        dist[s] = 0;
        parent[s] = s;
        queue.push(s);
        let mut head = 0;
        'bfs: while head < queue.len() {
            let v = queue[head];
            head += 1;
            // No cycle found from here can be shorter than the best.
            let floor = if g.is_directed() {
                dist[v] + 1
            } else {
                2 * dist[v] + 1
            };
            if floor >= best {
                break;
            }
            for &w in &out[v] {
                let w = w as usize;
                if g.is_directed() {
                    if w == s {
                        best = best.min(dist[v] + 1);
                        break 'bfs;
                    }
                    if dist[w] == usize::MAX {
                        dist[w] = dist[v] + 1;
                        queue.push(w);
                    }
                } else if dist[w] == usize::MAX {
                    dist[w] = dist[v] + 1;
                    parent[w] = v;
                    queue.push(w);
                } else if parent[v] != w {
                    best = best.min(dist[v] + dist[w] + 1);
                }
            }
        }
        for &v in &queue {
            dist[v] = usize::MAX;
            parent[v] = usize::MAX;
        }
    }
    if best == usize::MAX {
        f64::INFINITY
    } else {
        best as f64
    }
}

/// How [`find_cycle`] follows arcs of a directed network (NetworkX's
/// `orientation`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orientation {
    /// Arcs in their own direction.
    #[default]
    Original,
    /// Arcs against their direction.
    Reverse,
    /// Arcs either way.
    Ignore,
}

impl Orientation {
    /// Decodes the C interface codes: 0 = original, 1 = reverse,
    /// 2 = ignore.
    pub fn from_code(code: i32) -> Result<Orientation> {
        match code {
            0 => Ok(Orientation::Original),
            1 => Ok(Orientation::Reverse),
            2 => Ok(Orientation::Ignore),
            _ => Err(Error::InvalidArgument(format!(
                "unknown orientation code {code}"
            ))),
        }
    }
}

/// One step of a cycle found by [`find_cycle`]: the edge (input position)
/// and whether it was followed against its direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CycleStep {
    /// Input position of the edge.
    pub edge: usize,
    /// Whether the edge was traversed from its target to its source.
    pub reversed: bool,
}

/// A cycle found by depth-first search over edges from `source` (or from
/// every node in turn), as NetworkX's `find_cycle` finds it; `None` when
/// there is no cycle. See the module documentation.
pub fn find_cycle(
    g: &Graph,
    source: Option<usize>,
    orientation: Orientation,
    order: NeighborOrder,
) -> Result<Option<Vec<CycleStep>>> {
    let n = g.node_count();
    if let Some(s) = source {
        if s >= n {
            return Err(Error::InvalidArgument(format!(
                "source {s} is not a node (n = {n})"
            )));
        }
    }
    let directed = g.is_directed();
    let orientation = if directed {
        orientation
    } else {
        Orientation::Original
    };
    // Edges from each node in NetworkX's order: neighbours by first
    // appearance (or sorted), parallel edges in input order.
    let grouped = |dir: Direction| -> Vec<Vec<(u32, u32)>> {
        let r = rows(g, dir, order, false, false);
        (0..n)
            .map(|v| {
                let mut list: Vec<(u32, u32)> = Vec::new();
                let mut seen: Vec<u32> = Vec::new();
                for &w in r.neighbors(v) {
                    if seen.contains(&w) {
                        continue;
                    }
                    seen.push(w);
                    for (&x, &e) in r.neighbors(v).iter().zip(r.edges(v)) {
                        if x == w {
                            list.push((w, e));
                        }
                    }
                }
                list
            })
            .collect()
    };
    let out = grouped(Direction::Out);
    let inc = if directed {
        grouped(Direction::In)
    } else {
        Vec::new()
    };
    let (src, tgt) = (g.edge_sources(), g.edge_targets());
    // (tail, head) of a step as find_cycle reads it.
    let tailhead = |st: CycleStep, from: u32| -> (u32, u32) {
        let (a, b) = (src[st.edge], tgt[st.edge]);
        if directed {
            if st.reversed {
                (b, a)
            } else {
                (a, b)
            }
        } else if a == from {
            (a, b)
        } else {
            (b, a)
        }
    };
    let starts: Vec<usize> = match source {
        Some(s) => vec![s],
        None => (0..n).collect(),
    };
    let mut explored = vec![false; n];
    for &start in &starts {
        if explored[start] {
            continue;
        }
        // edge_dfs from `start`, with its own visited sets.
        let mut visited_nodes = vec![false; n];
        let mut visited_edges = vec![false; g.edge_count()];
        let mut cursor: Vec<usize> = vec![0; n];
        let mut stack: Vec<usize> = vec![start];
        let mut steps: Vec<(CycleStep, u32)> = Vec::new();
        let mut seen = vec![false; n];
        let mut active = vec![false; n];
        seen[start] = true;
        active[start] = true;
        let mut previous_head: Option<u32> = None;
        // Edges out of `node` in edge_dfs order: forward ones, then (for
        // Ignore) reverse ones.
        let edge_at = |node: usize, k: usize| -> Option<(CycleStep, u32)> {
            let fwd = match orientation {
                Orientation::Reverse => 0,
                _ => out[node].len(),
            };
            if k < fwd {
                let (w, e) = out[node][k];
                return Some((
                    CycleStep {
                        edge: e as usize,
                        reversed: false,
                    },
                    w,
                ));
            }
            let back = match orientation {
                Orientation::Original => 0,
                _ => inc[node].len(),
            };
            if k - fwd < back {
                let (w, e) = inc[node][k - fwd];
                return Some((
                    CycleStep {
                        edge: e as usize,
                        reversed: true,
                    },
                    w,
                ));
            }
            None
        };
        let mut cycle: Option<(Vec<(CycleStep, u32)>, u32)> = None;
        while let Some(&current) = stack.last() {
            if !visited_nodes[current] {
                visited_nodes[current] = true;
            }
            let k = cursor[current];
            let Some((step, next)) = edge_at(current, k) else {
                stack.pop();
                continue;
            };
            cursor[current] += 1;
            if visited_edges[step.edge] {
                continue;
            }
            visited_edges[step.edge] = true;
            stack.push(next as usize);
            // find_cycle's processing of the edge.
            let (tail, head) = tailhead(step, current as u32);
            if explored[head as usize] {
                continue;
            }
            if let Some(ph) = previous_head {
                if tail != ph {
                    loop {
                        match steps.pop() {
                            None => {
                                for a in active.iter_mut() {
                                    *a = false;
                                }
                                active[tail as usize] = true;
                                break;
                            }
                            Some((popped, from)) => {
                                let (_, h) = tailhead(popped, from);
                                active[h as usize] = false;
                            }
                        }
                        if let Some(&(last, from)) = steps.last() {
                            if tailhead(last, from).1 == tail {
                                break;
                            }
                        }
                    }
                }
            }
            steps.push((step, current as u32));
            if active[head as usize] {
                cycle = Some((steps.clone(), head));
                break;
            }
            seen[head as usize] = true;
            active[head as usize] = true;
            previous_head = Some(head);
        }
        if let Some((steps, final_node)) = cycle {
            // NetworkX keeps the last edge when no edge leaves the final
            // node (which cannot happen for a closed cycle).
            let begin = steps
                .iter()
                .position(|&(st, from)| tailhead(st, from).0 == final_node)
                .unwrap_or(steps.len() - 1);
            return Ok(Some(steps[begin..].iter().map(|x| x.0).collect()));
        }
        for v in 0..n {
            if seen[v] {
                explored[v] = true;
            }
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, e: &[(i64, i64)], directed: bool) -> Graph {
        let (s, t): (Vec<i64>, Vec<i64>) = e.iter().copied().unzip();
        Graph::from_edges(n, &s, &t, None, directed).unwrap()
    }

    fn canon(mut c: Vec<u32>, directed: bool) -> Vec<u32> {
        let i = (0..c.len()).min_by_key(|&i| c[i]).unwrap();
        c.rotate_left(i);
        if !directed && c.len() > 2 && c[1] > c[c.len() - 1] {
            c[1..].reverse();
        }
        c
    }

    #[test]
    fn bases_and_girth() {
        // Two triangles sharing an edge (a diamond) plus a pendant.
        let g = graph(5, &[(0, 1), (1, 2), (2, 0), (1, 3), (3, 2), (3, 4)], false);
        let b = cycle_basis(&g, NeighborOrder::Input, None).unwrap();
        assert_eq!(b.len(), 2);
        let mb = minimum_cycle_basis(&g, false, 2).unwrap();
        let mut w = mb.weights.clone();
        w.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(w, vec![3.0, 3.0]);
        assert_eq!(girth(&g), 3.0);
        assert_eq!(girth(&graph(3, &[(0, 1), (1, 2)], false)), f64::INFINITY);
        assert_eq!(girth(&graph(3, &[(0, 1), (1, 2), (2, 0)], true)), 3.0);
        assert_eq!(girth(&graph(2, &[(0, 1), (0, 1)], false)), 2.0);
    }

    #[test]
    fn simple_and_chordless() {
        // A 4-cycle with one chord: three simple cycles, two chordless.
        let g = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)], false);
        let mut s: Vec<Vec<u32>> = simple_cycles(&g, None, None)
            .unwrap()
            .into_iter()
            .map(|c| canon(c, false))
            .collect();
        s.sort();
        assert_eq!(s, vec![vec![0, 1, 2], vec![0, 1, 2, 3], vec![0, 2, 3]]);
        let mut c: Vec<Vec<u32>> = chordless_cycles(&g, None, None)
            .unwrap()
            .into_iter()
            .map(|c| canon(c, false))
            .collect();
        c.sort();
        assert_eq!(c, vec![vec![0, 1, 2], vec![0, 2, 3]]);
        let d = graph(3, &[(0, 1), (1, 2), (2, 0), (1, 0)], true);
        assert_eq!(simple_cycles(&d, None, None).unwrap().len(), 2);
        assert_eq!(simple_cycles(&d, Some(2), None).unwrap().len(), 1);
    }

    #[test]
    fn finds_a_cycle() {
        let d = graph(4, &[(0, 1), (1, 2), (2, 0), (2, 3)], true);
        let c = find_cycle(&d, None, Orientation::Original, NeighborOrder::Input)
            .unwrap()
            .unwrap();
        assert_eq!(c.iter().map(|s| s.edge).collect::<Vec<_>>(), vec![0, 1, 2]);
        let t = graph(3, &[(0, 1), (1, 2)], true);
        assert!(
            find_cycle(&t, None, Orientation::Original, NeighborOrder::Input)
                .unwrap()
                .is_none()
        );
        assert!(
            find_cycle(&t, None, Orientation::Ignore, NeighborOrder::Input)
                .unwrap()
                .is_none()
        );
    }
}
