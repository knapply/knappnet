//! C ABI for graph isomorphism, subgraph isomorphism, and monomorphism
//! (VF2++), the quick isomorphism checks, and tree isomorphism.
//!
//! Integer codes used by these functions:
//! - `problem`: 0 = isomorphism, 1 = induced subgraph isomorphism,
//!   2 = monomorphism.
//! - `properties` (quick checks): bit 1 = degrees, 2 = triangles, 4 =
//!   maximal cliques per node.

use std::ffi::c_char;

use engine::metrics::isomorphism::{
    could_be_isomorphic, count_mappings, rooted_tree_isomorphism, tree_isomorphism, Matcher,
    Problem,
};
use engine::Error;

use crate::{count, deref, guard, handle, slice, slice_mut, thread_count, GraphHandle};

/// A VF2++ search in progress; release with `net_vf2pp_free`.
pub struct MatcherHandle {
    matcher: Matcher,
    sg_nodes: usize,
}

/// Labels of the two graphs.
type Labels<'a> = (Option<&'a [i64]>, Option<&'a [i64]>);

/// Labels of both graphs (both null or both given).
///
/// # Safety
/// Each pointer must be null or hold the node count of its graph.
unsafe fn label_pair<'a>(
    a: *const i64,
    na: usize,
    b: *const i64,
    nb: usize,
) -> Result<Labels<'a>, Error> {
    match (a.is_null(), b.is_null()) {
        (true, true) => Ok((None, None)),
        (false, false) => Ok((Some(slice(a, na)), Some(slice(b, nb)))),
        _ => Err(Error::InvalidArgument(
            "labels must be given for both graphs or for neither".into(),
        )),
    }
}

/// # Safety
/// As for `net_vf2pp_new`.
unsafe fn matcher(
    fg: *const GraphHandle,
    sg: *const GraphHandle,
    fg_labels: *const i64,
    sg_labels: *const i64,
    problem: i32,
    empty_matches: i32,
) -> Result<(Matcher, usize), Error> {
    let f = handle(fg)?;
    let s = handle(sg)?;
    let (lf, ls) = label_pair(fg_labels, f.node_count(), sg_labels, s.node_count())?;
    let m = Matcher::new(
        f,
        s,
        lf,
        ls,
        Problem::from_code(problem)?,
        empty_matches != 0,
    )?;
    Ok((m, s.node_count()))
}

/// Starts a VF2++ search for the mappings of `sg` into `fg`. Labels are
/// both NULL or hold one value per node of each graph (equal labels
/// match). With `empty_matches == 0`, a graph without nodes gives no
/// mapping, as in NetworkX; otherwise an `sg` without nodes has one empty
/// mapping. Release `*out` with `net_vf2pp_free`.
///
/// # Safety
/// `fg` and `sg` must be live handles; labels null or of the node counts;
/// `out` writable.
#[no_mangle]
pub unsafe extern "C" fn net_vf2pp_new(
    fg: *const GraphHandle,
    sg: *const GraphHandle,
    fg_labels: *const i64,
    sg_labels: *const i64,
    problem: i32,
    empty_matches: i32,
    out: *mut *mut MatcherHandle,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let (m, sg_nodes) = matcher(fg, sg, fg_labels, sg_labels, problem, empty_matches)?;
        *out = Box::into_raw(Box::new(MatcherHandle {
            matcher: m,
            sg_nodes,
        }));
        Ok(())
    })
}

/// Continues a search: writes up to `max` further mappings into
/// `mappings` (mapping i occupies mappings[i * n_sg .. (i + 1) * n_sg),
/// the node of `fg` for every node of `sg`) and their number into
/// `*found`; fewer than `max` means the search is over.
///
/// # Safety
/// `h` must be a live handle; `mappings` must hold `max * n_sg` values and
/// `found` one.
#[no_mangle]
pub unsafe extern "C" fn net_vf2pp_next(
    h: *mut MatcherHandle,
    max: i64,
    mappings: *mut i64,
    found: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        deref(h as *const MatcherHandle, "matcher")?;
        let h = &mut *h;
        let max = count(max, "max")?;
        let ns = h.sg_nodes;
        let out = slice_mut(mappings, max * ns);
        let mut k = 0;
        while k < max {
            match h.matcher.next() {
                Some(m) => {
                    for (d, &v) in out[k * ns..(k + 1) * ns].iter_mut().zip(&m) {
                        *d = v as i64;
                    }
                    k += 1;
                }
                None => break,
            }
        }
        slice_mut(found, 1)[0] = k as i64;
        Ok(())
    })
}

/// Releases a search (NULL is allowed).
///
/// # Safety
/// `h` must be null or a live handle that is not used again.
#[no_mangle]
pub unsafe extern "C" fn net_vf2pp_free(h: *mut MatcherHandle) {
    if !h.is_null() {
        drop(Box::from_raw(h));
    }
}

/// Number of mappings of `sg` into `fg`, counting at most `limit` (1
/// answers whether one exists). Arguments as for `net_vf2pp_new`.
///
/// # Safety
/// As for `net_vf2pp_new`; `out` must hold one value.
#[no_mangle]
pub unsafe extern "C" fn net_vf2pp_count(
    fg: *const GraphHandle,
    sg: *const GraphHandle,
    fg_labels: *const i64,
    sg_labels: *const i64,
    problem: i32,
    empty_matches: i32,
    limit: i64,
    out: *mut i64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let (mut m, _) = matcher(fg, sg, fg_labels, sg_labels, problem, empty_matches)?;
        let limit = count(limit, "limit")? as u64;
        slice_mut(out, 1)[0] = count_mappings(&mut m, limit) as i64;
        Ok(())
    })
}

/// NetworkX's `could_be_isomorphic` checks: `*out` = 1 when the node
/// counts agree and, cumulatively, the sorted degree, (degree,
/// triangles), and (degree, triangles, cliques) sequences selected by
/// `properties` agree; 0 otherwise. Triangles and cliques need undirected
/// graphs.
///
/// # Safety
/// `g1` and `g2` must be live handles; `out` must hold one value.
#[no_mangle]
pub unsafe extern "C" fn net_could_be_isomorphic(
    g1: *const GraphHandle,
    g2: *const GraphHandle,
    properties: i32,
    threads: i32,
    out: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let (a, b) = (handle(g1)?, handle(g2)?);
        let ok = could_be_isomorphic(
            a,
            b,
            properties & 1 != 0,
            properties & 2 != 0,
            properties & 4 != 0,
            thread_count(threads),
        )?;
        slice_mut(out, 1)[0] = ok as i32;
        Ok(())
    })
}

/// Tree isomorphism (direction ignored): with `root1` and `root2` >= 0,
/// rooted trees with the roots mapped onto each other
/// (`rooted_tree_isomorphism`); with both < 0, free trees
/// (`tree_isomorphism`). `*found` = 1 and `mapping` (n1 values: the node
/// of `t2` for every node of `t1`) when isomorphic, else `*found` = 0. A
/// graph that is not a tree is an error.
///
/// # Safety
/// `t1` and `t2` must be live handles; `mapping` must hold n1 values and
/// `found` one.
#[no_mangle]
pub unsafe extern "C" fn net_tree_isomorphism(
    t1: *const GraphHandle,
    t2: *const GraphHandle,
    root1: i64,
    root2: i64,
    mapping: *mut i64,
    found: *mut i32,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let (a, b) = (handle(t1)?, handle(t2)?);
        let m = if root1 >= 0 && root2 >= 0 {
            rooted_tree_isomorphism(a, root1 as usize, b, root2 as usize)?
        } else if root1 < 0 && root2 < 0 {
            tree_isomorphism(a, b)?
        } else {
            return Err(Error::InvalidArgument("give both roots or neither".into()));
        };
        slice_mut(found, 1)[0] = m.is_some() as i32;
        if let Some(m) = m {
            for (d, &v) in slice_mut(mapping, a.node_count()).iter_mut().zip(&m) {
                *d = v as i64;
            }
        }
        Ok(())
    })
}
