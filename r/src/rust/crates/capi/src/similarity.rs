//! C ABI for link prediction and node and graph similarity.
//!
//! Integer codes used by these functions:
//! - `method` (link prediction): 0 = common neighbours, 1 = Jaccard,
//!   2 = Adamic-Adar, 3 = resource allocation, 4 = preferential
//!   attachment, 5 = common-neighbour centrality, 6 = CN
//!   Soundarajan-Hopcroft, 7 = RA Soundarajan-Hopcroft, 8 =
//!   within-inter-cluster.
//! - `stop` (SimRank): 0 = NetworkX's `simrank_similarity` (allclose),
//!   1 = its pairwise iteration (relative change).
//! - Labels for Weisfeiler-Lehman hashes are NUL-terminated ASCII strings
//!   (an array of n node labels or m edge labels, or NULL); hashes come
//!   back as hexadecimal text without separators or terminators.

use std::ffi::{c_char, CStr};

use engine::metrics::edit_distance::{graph_edit_distance, EditCosts, EditOptions};
use engine::metrics::linkpred::{
    link_scores, link_scores_non_edges, non_edge_count, LinkMethod, LinkOptions,
};
use engine::metrics::similarity::{
    panther_similarity, simrank, PantherOptions, SimRankOptions, SimRankStop, Similar,
};
use engine::metrics::wl::{wl_graph_hash, wl_subgraph_hashes, WlOptions};
use engine::{Error, Graph};

use crate::{count, guard, handle, slice, slice_mut, thread_count, GraphHandle};

/// Node indices from a caller's 64-bit array, each checked against the
/// graph's node count.
///
/// # Safety
/// `ptr` must be null or point to `len` readable values.
pub(crate) unsafe fn node_indices(
    g: &Graph,
    ptr: *const i64,
    len: i64,
    what: &str,
) -> Result<Vec<u32>, Error> {
    let len = usize::try_from(len)
        .map_err(|_| Error::InvalidArgument(format!("the length of {what} is negative")))?;
    let n = g.node_count();
    slice(ptr, len)
        .iter()
        .map(|&v| {
            if v < 0 || v as u64 >= n as u64 {
                Err(Error::InvalidArgument(format!(
                    "{what} contains node index {v}, but valid indices are 0..{n}"
                )))
            } else {
                Ok(v as u32)
            }
        })
        .collect()
}

/// Community labels (n values) or none when `ptr` is null.
///
/// # Safety
/// `ptr` must be null or point to `n` readable values.
unsafe fn labels<'a>(g: &Graph, ptr: *const i64) -> Option<&'a [i64]> {
    (!ptr.is_null()).then(|| slice(ptr, g.node_count()))
}

fn link_options(method: i32, alpha: f64, delta: f64, threads: i32) -> Result<LinkOptions, Error> {
    Ok(LinkOptions {
        method: LinkMethod::from_code(method)?,
        alpha,
        delta,
        threads: thread_count(threads),
    })
}

/// Link prediction scores of `npairs` pairs `(us[i], vs[i])` into `out`.
/// `communities` holds one label per node (methods 6 to 8; may be null
/// otherwise).
///
/// # Safety
/// `h` must be a live handle; `us`, `vs`, and `out` must hold `npairs`
/// values; `communities` must be null or hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_link_prediction(
    h: *const GraphHandle,
    method: i32,
    alpha: f64,
    delta: f64,
    communities: *const i64,
    us: *const i64,
    vs: *const i64,
    npairs: i64,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let a = node_indices(g, us, npairs, "us")?;
        let b = node_indices(g, vs, npairs, "vs")?;
        let opts = link_options(method, alpha, delta, threads)?;
        let scores = link_scores(g, &a, &b, labels(g, communities), opts)?;
        slice_mut(out, scores.len()).copy_from_slice(&scores);
        Ok(())
    })
}

/// Number of non-edges (unordered pairs of distinct nodes adjacent in
/// neither direction) into `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold one value.
#[no_mangle]
pub unsafe extern "C" fn net_non_edge_count(
    h: *const GraphHandle,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        slice_mut(out, 1)[0] = non_edge_count(g) as i64;
        Ok(())
    })
}

/// Link prediction scores of every non-edge u < v in lexicographic order:
/// `us`, `vs`, and `scores` must each hold `len` values, where `len` is
/// the count from `net_non_edge_count`.
///
/// # Safety
/// `h` must be a live handle; `us`, `vs`, and `scores` must hold `len`
/// values; `communities` must be null or hold n values.
#[no_mangle]
pub unsafe extern "C" fn net_link_prediction_non_edges(
    h: *const GraphHandle,
    method: i32,
    alpha: f64,
    delta: f64,
    communities: *const i64,
    threads: i32,
    us: *mut i64,
    vs: *mut i64,
    scores: *mut f64,
    len: i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let len = usize::try_from(len)
            .map_err(|_| Error::InvalidArgument("len must be non-negative".into()))?;
        let opts = link_options(method, alpha, delta, threads)?;
        let (mut a, mut b) = (vec![0u32; len], vec![0u32; len]);
        let s = slice_mut(scores, len);
        link_scores_non_edges(g, labels(g, communities), opts, &mut a, &mut b, s)?;
        for (d, &x) in slice_mut(us, len).iter_mut().zip(&a) {
            *d = x as i64;
        }
        for (d, &x) in slice_mut(vs, len).iter_mut().zip(&b) {
            *d = x as i64;
        }
        Ok(())
    })
}

/// `len` NUL-terminated strings, or none when `ptr` is null.
///
/// # Safety
/// `ptr` must be null or point to `len` pointers to NUL-terminated strings.
unsafe fn strings(
    ptr: *const *const c_char,
    len: usize,
    what: &str,
) -> Result<Option<Vec<String>>, Error> {
    if ptr.is_null() {
        return Ok(None);
    }
    slice(ptr, len)
        .iter()
        .map(|&p| {
            if p.is_null() {
                return Err(Error::InvalidArgument(format!(
                    "{what} contains a null pointer"
                )));
            }
            CStr::from_ptr(p)
                .to_str()
                .map(str::to_string)
                .map_err(|_| Error::InvalidArgument(format!("{what} must be UTF-8 text")))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

/// Writes hexadecimal text into a caller's buffer of exactly `len` bytes.
///
/// # Safety
/// `out` must point to `len` writable bytes.
unsafe fn write_text(text: &[u8], out: *mut c_char, len: usize) -> Result<(), Error> {
    if text.len() != len {
        return Err(Error::LengthMismatch {
            what: "output buffer",
            expected: text.len(),
            found: len,
        });
    }
    let dst = slice_mut(out as *mut u8, len);
    dst.copy_from_slice(text);
    Ok(())
}

/// Weisfeiler-Lehman graph hash, as NetworkX's
/// `weisfeiler_lehman_graph_hash`: `out` receives 2 `digest_size`
/// hexadecimal characters (no terminator; `out_len` must equal that).
///
/// # Safety
/// `h` must be a live handle; `node_labels` must be null or hold n
/// strings, `edge_labels` null or m strings; `out` must hold `out_len`
/// bytes.
#[no_mangle]
pub unsafe extern "C" fn net_wl_graph_hash(
    h: *const GraphHandle,
    node_labels: *const *const c_char,
    edge_labels: *const *const c_char,
    iterations: i64,
    digest_size: i64,
    threads: i32,
    out: *mut c_char,
    out_len: usize,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let nl = strings(node_labels, g.node_count(), "node_labels")?;
        let el = strings(edge_labels, g.edge_count(), "edge_labels")?;
        let opts = WlOptions {
            iterations: count(iterations, "iterations")?,
            digest_size: count(digest_size, "digest_size")?,
            node_labels: nl.as_deref(),
            edge_labels: el.as_deref(),
            include_initial_labels: false,
            threads: thread_count(threads),
        };
        write_text(wl_graph_hash(g, opts)?.as_bytes(), out, out_len)
    })
}

/// Weisfeiler-Lehman subgraph hashes, as NetworkX's
/// `weisfeiler_lehman_subgraph_hashes`: every node has L = `iterations`
/// (+ 1 with `include_initial_labels`) hashes of 2 `digest_size`
/// characters; node v's hashes, oldest first, fill
/// out[v * L * 2 digest_size ..), and `out_len` must be n L 2 digest_size.
///
/// # Safety
/// As for `net_wl_graph_hash`.
#[no_mangle]
pub unsafe extern "C" fn net_wl_subgraph_hashes(
    h: *const GraphHandle,
    node_labels: *const *const c_char,
    edge_labels: *const *const c_char,
    iterations: i64,
    digest_size: i64,
    include_initial_labels: i32,
    threads: i32,
    out: *mut c_char,
    out_len: usize,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let nl = strings(node_labels, g.node_count(), "node_labels")?;
        let el = strings(edge_labels, g.edge_count(), "edge_labels")?;
        let opts = WlOptions {
            iterations: count(iterations, "iterations")?,
            digest_size: count(digest_size, "digest_size")?,
            node_labels: nl.as_deref(),
            edge_labels: el.as_deref(),
            include_initial_labels: include_initial_labels != 0,
            threads: thread_count(threads),
        };
        let hashes = wl_subgraph_hashes(g, opts)?;
        let text: String = hashes.iter().flatten().map(String::as_str).collect();
        write_text(text.as_bytes(), out, out_len)
    })
}

/// SimRank similarity of every pair of nodes: an n-by-n symmetric matrix
/// in `out`. `stop`: 0 = NetworkX's `simrank_similarity`, 1 = its pairwise
/// iteration (unweighted). Reaching `max_iterations` is an error.
///
/// # Safety
/// `h` must be a live handle; `out` must hold n * n doubles.
#[no_mangle]
pub unsafe extern "C" fn net_simrank(
    h: *const GraphHandle,
    importance: f64,
    max_iterations: i64,
    tolerance: f64,
    weighted: i32,
    stop: i32,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let opts = SimRankOptions {
            importance,
            max_iterations: count(max_iterations, "max_iterations")?,
            tolerance,
            weighted: weighted != 0,
            stop: SimRankStop::from_code(stop)?,
            threads: thread_count(threads),
        };
        let s = simrank(g, opts)?;
        let n = g.node_count();
        slice_mut(out, n * n).copy_from_slice(&s.to_col_major());
        Ok(())
    })
}

#[allow(clippy::too_many_arguments)]
fn panther_options(
    k: i64,
    path_length: i64,
    c: f64,
    delta: f64,
    eps: f64,
    weighted: i32,
    seed: u64,
    threads: i32,
) -> Result<PantherOptions, Error> {
    Ok(PantherOptions {
        k: count(k, "k")?,
        path_length: count(path_length, "path_length")?,
        c,
        delta,
        eps: (eps > 0.0).then_some(eps),
        weighted: weighted != 0,
        seed,
        threads: thread_count(threads),
    })
}

/// # Safety
/// `nodes` and `scores` must hold `capacity` values and `found` one.
unsafe fn write_similar(
    r: &Similar,
    nodes: *mut i64,
    scores: *mut f64,
    capacity: i64,
    found: *mut i64,
) -> Result<(), Error> {
    let cap = count(capacity, "capacity")?;
    if r.nodes.len() > cap {
        return Err(Error::InvalidArgument(format!(
            "the result has {} nodes but capacity is {cap}",
            r.nodes.len()
        )));
    }
    for (d, &v) in slice_mut(nodes, cap).iter_mut().zip(&r.nodes) {
        *d = v as i64;
    }
    slice_mut(scores, cap)[..r.scores.len()].copy_from_slice(&r.scores);
    slice_mut(found, 1)[0] = r.nodes.len() as i64;
    Ok(())
}

/// Panther similarity: the k nodes most similar to `source`, most similar
/// first, into `nodes` and `scores` (`*found` of them: k, or k + 1 when
/// more than k + 1 nodes tie with the source; capacity must be at least
/// k + 1). `eps <= 0` uses 1/sqrt(m). Paths come from stream (seed, path).
///
/// # Safety
/// `h` must be a live handle; `nodes` and `scores` must hold `capacity`
/// values and `found` one.
#[no_mangle]
pub unsafe extern "C" fn net_panther_similarity(
    h: *const GraphHandle,
    source: i64,
    k: i64,
    path_length: i64,
    c: f64,
    delta: f64,
    eps: f64,
    weighted: i32,
    seed: u64,
    threads: i32,
    nodes: *mut i64,
    scores: *mut f64,
    capacity: i64,
    found: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let s = node_indices(g, &source, 1, "source")?[0] as usize;
        let opts = panther_options(k, path_length, c, delta, eps, weighted, seed, threads)?;
        let r = panther_similarity(g, s, opts)?;
        write_similar(&r, nodes, scores, capacity, found)
    })
}

/// Exact graph edit distance between `h1` and `h2` into `out`, with
/// `*found` = 1, or `*found` = 0 when it exceeds `upper_bound` (pass
/// INFINITY for none). Node labels (n1 and n2 values) and edge labels (m1
/// and m2 values) are both NULL or both given; equal labels substitute at
/// no cost. `costs` holds node substitution, deletion, insertion, then
/// edge substitution, deletion, insertion costs (NULL: all 1).
/// `root1`, `root2` >= 0 fix that pair. The search fails after
/// `max_steps` mappings (0: 10^8).
///
/// # Safety
/// `h1` and `h2` must be live handles; label arrays null or of the sizes
/// above; `costs` null or 6 doubles; `out` and `found` one value each.
#[no_mangle]
pub unsafe extern "C" fn net_graph_edit_distance(
    h1: *const GraphHandle,
    h2: *const GraphHandle,
    node_labels1: *const i64,
    node_labels2: *const i64,
    edge_labels1: *const i64,
    edge_labels2: *const i64,
    costs: *const f64,
    root1: i64,
    root2: i64,
    upper_bound: f64,
    max_steps: i64,
    out: *mut f64,
    found: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g1 = handle(h1)?;
        let g2 = handle(h2)?;
        let pair = |a: *const i64, b: *const i64, n1: usize, n2: usize, what: &str| match (
            a.is_null(),
            b.is_null(),
        ) {
            (true, true) => Ok(None),
            (false, false) => Ok(Some((slice(a, n1), slice(b, n2)))),
            _ => Err(Error::InvalidArgument(format!(
                "{what} must be given for both networks or for neither"
            ))),
        };
        let nl = pair(
            node_labels1,
            node_labels2,
            g1.node_count(),
            g2.node_count(),
            "node labels",
        )?;
        let el = pair(
            edge_labels1,
            edge_labels2,
            g1.edge_count(),
            g2.edge_count(),
            "edge labels",
        )?;
        let costs = if costs.is_null() {
            EditCosts::default()
        } else {
            let c = slice(costs, 6);
            EditCosts {
                node_subst: c[0],
                node_del: c[1],
                node_ins: c[2],
                edge_subst: c[3],
                edge_del: c[4],
                edge_ins: c[5],
            }
        };
        let roots = if root1 >= 0 && root2 >= 0 {
            Some((root1 as usize, root2 as usize))
        } else {
            None
        };
        let opts = EditOptions {
            node_labels: nl,
            edge_labels: el,
            costs,
            roots,
            upper_bound: (!upper_bound.is_nan() && upper_bound.is_finite()).then_some(upper_bound),
            max_steps: count(max_steps, "max_steps")? as u64,
        };
        let d = graph_edit_distance(g1, g2, opts)?;
        slice_mut(out, 1)[0] = d.unwrap_or(f64::NAN);
        slice_mut(found, 1)[0] = d.is_some() as i32;
        Ok(())
    })
}
