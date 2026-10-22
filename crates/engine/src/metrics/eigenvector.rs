//! Eigenvector centrality.
//!
//! A node's score is proportional to the sum of its neighbors' scores: the
//! entries of the principal eigenvector of the adjacency matrix (Bonacich
//! 1972, `bonacich1972factoring`; see also Bonacich 1987,
//! `bonacich1987power`). For a nonnegative matrix whose graph is strongly
//! connected, the Perron-Frobenius theorem guarantees that this eigenvector
//! is unique up to scale and positive.
//!
//! Algorithm: power iteration (von Mises and Pollaczek-Geiringer 1929,
//! `mises1929praktische`) on A' + I, where A' is the (transposed, for
//! directed graphs) adjacency matrix. The shift by the identity leaves the
//! eigenvectors unchanged and makes the dominant eigenvalue strictly largest
//! in modulus even for bipartite and periodic graphs, where plain power
//! iteration on A oscillates; NetworkX's `eigenvector_centrality` uses the
//! same shift. The iteration starts from the all-ones vector, rescales to a
//! largest entry of 1 after every step, and stops when no entry changes by
//! more than `tol` (default 1e-12). If that does not happen within
//! `max_iter` steps (default 10,000) an error is returned rather than an
//! inaccurate answer. Convergence is geometric with ratio
//! (lambda_2 + 1) / (lambda_1 + 1) for the two largest eigenvalues, so it
//! is slow when they are close. Each step costs O(n + m) and runs in
//! parallel on large graphs; the result does not depend on the thread
//! count.
//!
//! For symmetric matrices (undirected graphs, or mode `All`) each connected
//! component is iterated with its own normalization, and a component is
//! dropped (scores 0) as soon as Collatz's bounds, which place the largest
//! eigenvalue of an irreducible nonnegative matrix between the minimum and
//! maximum of (A'x)_i / x_i for any positive x (Collatz 1942,
//! `collatz1942einschliessungssatz`), show that its largest eigenvalue is
//! below another component's.
//!
//! Conventions:
//! - Directed graphs: with mode `In` (igraph with `directed = TRUE`,
//!   NetworkX) a node's score is the sum of the scores of nodes with arcs
//!   *to* it, the left eigenvector of A. With `Out` (sna's `evcent`) it sums
//!   the nodes it points to, the right eigenvector. `All` ignores direction,
//!   so two opposite arcs become two parallel edges (igraph's R default,
//!   `directed = FALSE`). Undirected graphs ignore the mode.
//! - Parallel edges add up: their count, or their weights.
//! - Self-loops ([`LoopWeight`]): igraph puts twice the number of loops on
//!   the diagonal of an undirected graph (each loop is two edge endpoints)
//!   and NetworkX once; in directed graphs both count a loop once. sna
//!   ignores the diagonal.
//! - Weights are used as given (larger means a stronger tie) and must be
//!   nonnegative.
//! - Scaling ([`EigenScale`]): `Max` makes the largest score 1 (igraph);
//!   `Unit` gives the vector Euclidean length 1 (NetworkX, sna, UCINET).
//!   Scores are nonnegative.
//! - Graphs without cycles (including graphs without edges): the only
//!   eigenvalue is 0. Following igraph, nodes whose score feeds no other
//!   node's score (nodes without outgoing arcs for mode `In`, without
//!   incoming arcs for `Out`, and every node in a graph without edges) score
//!   1 before scaling, every other node 0, and the eigenvalue is 0.
//! - Disconnected graphs (not strongly connected, for directed graphs) have
//!   no unique positive eigenvector, and [`Eigenvector::connected`] is false
//!   so callers can warn. igraph returns scores without checking; NetworkX's
//!   `eigenvector_centrality_numpy` refuses. Here:
//!   - Symmetric matrices: components whose largest eigenvalue is below the
//!     overall largest score 0. If several components attain it (to a
//!     relative 1e-9), each contributes its unit Perron vector u times
//!     sum(u), which is the limit of power iteration from equal scores.
//!   - Directed matrices: the result is the limit of the iteration, which
//!     is zero outside the strongly connected components that carry the
//!     principal eigenvalue and those they reach (`In`) or that reach them
//!     (`Out`). When several components share the principal eigenvalue and
//!     one reaches another, the matrix is defective, the iteration converges
//!     only sublinearly, and an error is returned.

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::metrics::closeness::oriented;
use crate::metrics::degree::DegreeMode;
use crate::metrics::topography::{components, ComponentKind};
use crate::parallel;

/// How eigenvector scores are scaled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EigenScale {
    /// Largest score 1 (igraph).
    #[default]
    Max,
    /// Euclidean norm 1 (NetworkX, sna, UCINET).
    Unit,
}

impl EigenScale {
    /// Decodes the C interface codes: 0 = max, 1 = unit.
    pub fn from_code(code: i32) -> Result<EigenScale> {
        match code {
            0 => Ok(EigenScale::Max),
            1 => Ok(EigenScale::Unit),
            _ => Err(Error::InvalidArgument(format!(
                "unknown eigenvector scale code {code}"
            ))),
        }
    }
}

/// What a self-loop contributes to the diagonal of an undirected graph's
/// adjacency matrix. Directed loops always count once unless ignored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LoopWeight {
    /// Twice its count or weight (igraph).
    #[default]
    Twice,
    /// Once (NetworkX).
    Once,
    /// Not at all (sna, UCINET).
    Ignore,
}

impl LoopWeight {
    /// Decodes the C interface codes: 0 = twice, 1 = once, 2 = ignore.
    pub fn from_code(code: i32) -> Result<LoopWeight> {
        match code {
            0 => Ok(LoopWeight::Twice),
            1 => Ok(LoopWeight::Once),
            2 => Ok(LoopWeight::Ignore),
            _ => Err(Error::InvalidArgument(format!(
                "unknown self-loop weight code {code}"
            ))),
        }
    }
}

/// Default convergence tolerance (largest change in any score, with the
/// largest score scaled to 1).
pub const DEFAULT_TOL: f64 = 1e-12;
/// Default iteration cap.
pub const DEFAULT_MAX_ITER: usize = 10_000;

/// Options for [`eigenvector_centrality`].
#[derive(Debug, Clone, Copy)]
pub struct EigenvectorOptions {
    /// Direction for directed graphs (`In` is the usual definition).
    pub mode: DegreeMode,
    /// Use edge weights as given (they must be nonnegative).
    pub weighted: bool,
    /// Scaling of the result.
    pub scale: EigenScale,
    /// Self-loops in undirected graphs.
    pub loops: LoopWeight,
    /// Convergence tolerance.
    pub tol: f64,
    /// Iteration cap.
    pub max_iter: usize,
    /// Worker threads; 0 uses all available cores.
    pub threads: usize,
}

impl Default for EigenvectorOptions {
    fn default() -> Self {
        EigenvectorOptions {
            mode: DegreeMode::In,
            weighted: false,
            scale: EigenScale::Max,
            loops: LoopWeight::Twice,
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
            threads: 0,
        }
    }
}

/// Result of [`eigenvector_centrality`].
#[derive(Debug, Clone, PartialEq)]
pub struct Eigenvector {
    /// Scores in node order.
    pub values: Vec<f64>,
    /// The principal eigenvalue (0 for graphs without cycles).
    pub eigenvalue: f64,
    /// Iterations used (0 when no iteration was needed).
    pub iterations: usize,
    /// Whether the graph is connected (strongly, if directed and the mode
    /// is not `All`). When false the scores depend on more than the
    /// network; see the module documentation.
    pub connected: bool,
}

/// Arcs whose number exceeds this use several threads per step.
const PARALLEL_ARCS: usize = 1 << 16;

/// Eigenvector centrality. See the module documentation for conventions.
pub fn eigenvector_centrality(g: &Graph, opts: EigenvectorOptions) -> Result<Eigenvector> {
    if opts.tol.is_nan() || opts.tol <= 0.0 || opts.max_iter == 0 {
        return Err(Error::InvalidArgument(
            "the tolerance and the iteration cap must be positive".into(),
        ));
    }
    let n = g.node_count();
    if n == 0 {
        return Ok(Eigenvector {
            values: Vec::new(),
            eigenvalue: f64::NAN,
            iterations: 0,
            connected: true,
        });
    }
    if opts.weighted {
        let w = g.edge_weights().ok_or(Error::MissingWeights)?;
        if let Some(edge) = w.iter().position(|&x| x < 0.0) {
            return Err(Error::InvalidArgument(format!(
                "edge {edge} has weight {}; eigenvector centrality requires \
                 nonnegative weights",
                w[edge]
            )));
        }
    }
    // Rows of the iteration matrix: x_i <- x_i + sum over row(i) of w x_j.
    // `In` sums over arcs arriving at i; `Out` over arcs leaving i.
    let mode = if g.is_directed() {
        opts.mode
    } else {
        DegreeMode::All
    };
    let h = oriented(
        g,
        if mode == DegreeMode::All {
            mode
        } else {
            DegreeMode::Out
        },
    )?;
    let h = h.as_ref();
    let rows = if h.is_directed() && mode == DegreeMode::In {
        h.inc()
    } else {
        h.out()
    };
    let loop_factor = match (h.is_directed(), opts.loops) {
        (_, LoopWeight::Ignore) => 0.0,
        (false, LoopWeight::Twice) => 2.0,
        _ => 1.0,
    };
    let weight = |v: usize, k: usize, j: usize| -> f64 {
        let w = if opts.weighted {
            rows.weights(v).map_or(1.0, |ws| ws[k])
        } else {
            1.0
        };
        if j == v {
            w * loop_factor
        } else {
            w
        }
    };

    let kind = if h.is_directed() {
        ComponentKind::Strong
    } else {
        ComponentKind::Weak
    };
    let connected = components(h, kind).count() == 1;

    // Acyclic case: every eigenvalue is 0. A graph has a cycle with
    // positive weight exactly when some strong component has a positive
    // arc inside it.
    let comp = components(h, ComponentKind::Strong).membership;
    let cyclic = (0..n).any(|v| {
        rows.neighbors(v)
            .iter()
            .enumerate()
            .any(|(k, &j)| comp[j as usize] == comp[v] && weight(v, k, j as usize) > 0.0)
    });
    if !cyclic {
        // Nodes with no positive entries in their column of the matrix
        // feed no other score: igraph's sinks, for mode `In`.
        let mut sink = vec![1.0; n];
        for v in 0..n {
            for (k, &j) in rows.neighbors(v).iter().enumerate() {
                if weight(v, k, j as usize) > 0.0 {
                    sink[j as usize] = 0.0;
                }
            }
        }
        scale(&mut sink, opts.scale);
        return Ok(Eigenvector {
            values: sink,
            eigenvalue: 0.0,
            iterations: 0,
            connected,
        });
    }

    let threads = if rows.arc_count() >= PARALLEL_ARCS {
        parallel::resolve_threads(opts.threads, n)
    } else {
        1
    };
    let step = |x: &[f64], y: &mut [f64]| {
        let row = |v: usize| -> f64 {
            rows.neighbors(v)
                .iter()
                .enumerate()
                .map(|(k, &j)| weight(v, k, j as usize) * x[j as usize])
                .sum()
        };
        if threads == 1 {
            for (v, yv) in y.iter_mut().enumerate() {
                *yv = row(v);
            }
            return;
        }
        // Contiguous blocks, one value per node, so the sums do not depend
        // on the thread count.
        let block = (n + threads - 1) / threads;
        let parts = parallel::for_each_item(
            threads,
            threads,
            Vec::new,
            |out: &mut Vec<(usize, Vec<f64>)>, b| {
                let lo = b * block;
                let hi = ((b + 1) * block).min(n);
                out.push((lo, (lo..hi).map(row).collect()));
            },
        );
        for part in parts {
            for (lo, vals) in part {
                y[lo..lo + vals.len()].copy_from_slice(&vals);
            }
        }
    };

    // Symmetric matrices iterate every connected component separately (each
    // normalized by its own largest score), so the result on a disconnected
    // graph is well defined and a slowly converging small component cannot
    // hold up the rest. Directed matrices iterate as one block.
    let (group, ngroups) = if h.is_directed() {
        (vec![0u32; n], 1)
    } else {
        let c = components(h, ComponentKind::Weak);
        let k = c.count();
        (c.membership, k)
    };
    let mut active = vec![true; ngroups];
    let mut x = vec![1.0; n];
    let mut y = vec![0.0; n];
    let mut top = vec![0.0f64; ngroups];
    let mut lower = vec![f64::INFINITY; ngroups];
    let mut upper = vec![0.0f64; ngroups];
    for iter in 1..=opts.max_iter {
        step(&x, &mut y);
        top.iter_mut().for_each(|t| *t = 0.0);
        lower.iter_mut().for_each(|t| *t = f64::INFINITY);
        upper.iter_mut().for_each(|t| *t = 0.0);
        for v in 0..n {
            let g = group[v] as usize;
            if !active[g] {
                continue;
            }
            // Collatz-Wielandt bounds: min and max of (A'x)_v / x_v bracket
            // the component's largest eigenvalue (x is positive on it).
            let ratio = y[v] / x[v];
            lower[g] = lower[g].min(ratio);
            upper[g] = upper[g].max(ratio);
            y[v] += x[v];
            top[g] = top[g].max(y[v]);
        }
        let mut change = 0.0f64;
        for v in 0..n {
            let g = group[v] as usize;
            if active[g] {
                y[v] /= top[g];
                change = change.max((y[v] - x[v]).abs());
            } else {
                y[v] = 0.0;
            }
        }
        std::mem::swap(&mut x, &mut y);
        if ngroups > 1 {
            // A component whose upper bound is below another's lower bound
            // cannot carry the principal eigenvalue; it scores 0.
            let best = (0..ngroups)
                .filter(|&g| active[g])
                .map(|g| lower[g])
                .fold(0.0f64, f64::max);
            for g in 0..ngroups {
                if active[g] && upper[g] < best * (1.0 - TIE) {
                    active[g] = false;
                }
            }
            for v in 0..n {
                if !active[group[v] as usize] {
                    x[v] = 0.0;
                }
            }
        }
        if change <= opts.tol {
            return Ok(finish(
                x, &group, ngroups, &active, step, opts, iter, connected,
            ));
        }
    }
    Err(Error::InvalidArgument(format!(
        "eigenvector centrality did not converge to tolerance {} within {} iterations. \
         Power iteration converges slowly when the largest eigenvalue is close to the \
         next one, and not at all to a single answer when several strongly connected \
         components share it; a larger iteration cap may help in the first case",
        opts.tol, opts.max_iter
    )))
}

/// Relative difference within which two components' largest eigenvalues
/// are treated as equal.
const TIE: f64 = 1e-9;

/// Combines converged components into the principal eigenvector.
#[allow(clippy::too_many_arguments)]
fn finish(
    mut x: Vec<f64>,
    group: &[u32],
    ngroups: usize,
    active: &[bool],
    step: impl Fn(&[f64], &mut [f64]),
    opts: EigenvectorOptions,
    iterations: usize,
    connected: bool,
) -> Eigenvector {
    let n = x.len();
    let mut y = vec![0.0; n];
    step(&x, &mut y);
    // Per component: sum of A'x, sum of x, and sum of squares of x.
    let mut sums = vec![(0.0f64, 0.0f64, 0.0f64); ngroups];
    for v in 0..n {
        let s = &mut sums[group[v] as usize];
        s.0 += y[v];
        s.1 += x[v];
        s.2 += x[v] * x[v];
    }
    let lambda: Vec<f64> = sums
        .iter()
        .zip(active)
        .map(|(s, &a)| if a && s.1 > 0.0 { s.0 / s.1 } else { 0.0 })
        .collect();
    let eigenvalue = lambda.iter().copied().fold(0.0f64, f64::max);
    // Each component that attains the largest eigenvalue contributes its
    // unit Perron vector u times sum(u), the limit of power iteration from
    // equal scores; the others contribute nothing.
    let weight: Vec<f64> = (0..ngroups)
        .map(|g| {
            let (_, s1, s2) = sums[g];
            if lambda[g] >= eigenvalue * (1.0 - TIE) && s2 > 0.0 {
                s1 / s2
            } else {
                0.0
            }
        })
        .collect();
    for v in 0..n {
        x[v] *= weight[group[v] as usize];
    }
    scale(&mut x, opts.scale);
    Eigenvector {
        values: x,
        eigenvalue,
        iterations,
        connected,
    }
}

fn scale(x: &mut [f64], how: EigenScale) {
    let d = match how {
        EigenScale::Max => x.iter().copied().fold(0.0f64, f64::max),
        EigenScale::Unit => x.iter().map(|v| v * v).sum::<f64>().sqrt(),
    };
    if d > 0.0 {
        for v in x {
            *v /= d;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(g: &Graph, opts: EigenvectorOptions) -> Eigenvector {
        eigenvector_centrality(g, opts).unwrap()
    }

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)
    }

    #[test]
    fn star_closed_form() {
        // Star with k leaves: eigenvalue sqrt(k), hub sqrt(k) times a leaf.
        let k = 4.0f64;
        let g = Graph::from_edges(5, &[0, 0, 0, 0], &[1, 2, 3, 4], None, false).unwrap();
        let e = run(&g, Default::default());
        assert!((e.eigenvalue - k.sqrt()).abs() < 1e-10);
        assert!(close(&e.values, &[1.0, 0.5, 0.5, 0.5, 0.5], 1e-10));
        assert!(e.connected);
        let u = run(
            &g,
            EigenvectorOptions {
                scale: EigenScale::Unit,
                ..Default::default()
            },
        );
        let norm: f64 = u.values.iter().map(|v| v * v).sum();
        assert!((norm - 1.0).abs() < 1e-12);
    }

    #[test]
    fn bipartite_graph_converges() {
        // A path is bipartite; plain power iteration on A would oscillate.
        let g = Graph::from_edges(4, &[0, 1, 2], &[1, 2, 3], None, false).unwrap();
        let e = run(&g, Default::default());
        let phi = (1.0 + 5f64.sqrt()) / 2.0;
        assert!((e.eigenvalue - phi).abs() < 1e-10);
        assert!(close(&e.values, &[1.0 / phi, 1.0, 1.0, 1.0 / phi], 1e-10));
    }

    #[test]
    fn directed_modes() {
        // 0 -> 1, 1 -> 2, 2 -> 0, 2 -> 1: in-scores favor 1.
        let g = Graph::from_edges(3, &[0, 1, 2, 2], &[1, 2, 0, 1], None, true).unwrap();
        let ins = run(&g, Default::default());
        let outs = run(
            &g,
            EigenvectorOptions {
                mode: DegreeMode::Out,
                ..Default::default()
            },
        );
        // Check A^T x = lambda x and A x = lambda x directly.
        let a = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0]];
        for (i, a_i) in a.iter().enumerate() {
            let col: f64 = (0..3).map(|j| a[j][i] * ins.values[j]).sum();
            let row: f64 = (0..3).map(|j| a_i[j] * outs.values[j]).sum();
            assert!((col - ins.eigenvalue * ins.values[i]).abs() < 1e-10);
            assert!((row - outs.eigenvalue * outs.values[i]).abs() < 1e-10);
        }
        assert!(ins.values[1] > ins.values[0]);
    }

    #[test]
    fn loops_and_weights() {
        // Undirected loop at 0 counts twice by default (A_00 = 2).
        let g = Graph::from_edges(2, &[0, 0], &[0, 1], Some(&[1.0, 3.0]), false).unwrap();
        let twice = run(&g, Default::default());
        // Binary A = [[2, 1], [1, 0]]: eigenvalue 1 + sqrt(2).
        assert!((twice.eigenvalue - (1.0 + 2f64.sqrt())).abs() < 1e-10);
        let once = run(
            &g,
            EigenvectorOptions {
                loops: LoopWeight::Once,
                weighted: true,
                ..Default::default()
            },
        );
        // Weighted A = [[1, 3], [3, 0]]: eigenvalue (1 + sqrt(37)) / 2.
        assert!((once.eigenvalue - (1.0 + 37f64.sqrt()) / 2.0).abs() < 1e-10);
        let neg = Graph::from_edges(2, &[0], &[1], Some(&[-1.0]), false).unwrap();
        let o = EigenvectorOptions {
            weighted: true,
            ..Default::default()
        };
        assert!(eigenvector_centrality(&neg, o).is_err());
    }

    #[test]
    fn acyclic_and_empty_graphs() {
        let empty = Graph::from_edges(3, &[], &[], None, false).unwrap();
        let e = run(&empty, Default::default());
        assert_eq!(e.values, vec![1.0, 1.0, 1.0]);
        assert_eq!(e.eigenvalue, 0.0);
        assert!(!e.connected);
        // 0 -> 1 -> 2: with in-ties the sink 2 scores 1.
        let dag = Graph::from_edges(3, &[0, 1], &[1, 2], None, true).unwrap();
        assert_eq!(run(&dag, Default::default()).values, vec![0.0, 0.0, 1.0]);
        let out = EigenvectorOptions {
            mode: DegreeMode::Out,
            ..Default::default()
        };
        assert_eq!(run(&dag, out).values, vec![1.0, 0.0, 0.0]);
    }

    #[test]
    fn disconnected_graph_is_flagged() {
        // Triangle plus an edge: the triangle (eigenvalue 2) dominates.
        let g = Graph::from_edges(5, &[0, 1, 2, 3], &[1, 2, 0, 4], None, false).unwrap();
        let e = run(&g, Default::default());
        assert!(!e.connected);
        assert!((e.eigenvalue - 2.0).abs() < 1e-10);
        assert!(close(&e.values, &[1.0, 1.0, 1.0, 0.0, 0.0], 1e-9));
    }

    #[test]
    fn thread_count_does_not_change_results() {
        let n = 3000usize;
        let mut s = Vec::new();
        let mut t = Vec::new();
        for i in 0..n {
            for d in [1, 7, 31, 101, 997, 1999] {
                s.push(i as i64);
                t.push(((i * d + 3) % n) as i64);
            }
        }
        let g = Graph::from_edges(n, &s, &t, None, true).unwrap();
        let one = run(
            &g,
            EigenvectorOptions {
                threads: 1,
                ..Default::default()
            },
        );
        let four = run(
            &g,
            EigenvectorOptions {
                threads: 4,
                ..Default::default()
            },
        );
        assert_eq!(one, four);
    }
}
