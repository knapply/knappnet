//! C ABI for cut size, volume, normalized cut, conductance, and edge,
//! mixing, node, and boundary expansion.

use std::ffi::c_char;

use engine::metrics::cut_measures::cut_measures;

use crate::similarity::node_indices;
use crate::{guard, handle, slice_mut, GraphHandle};

/// Every measure of the cut between node sets `s` (`ns` nodes) and `t`
/// (`nt` nodes; NULL for the complement of `s`), into `out`: cut size,
/// volume of S, volume of T, normalized cut size, conductance, edge
/// expansion, mixing expansion, node expansion, boundary expansion (as
/// NetworkX's `algorithms.cuts`; IEEE inf or NaN where NetworkX divides by
/// zero). `weighted != 0` adds weights instead of counting edges.
///
/// # Safety
/// `h` must be a live handle; `s` must hold `ns` values, `t` be null or
/// hold `nt`; `out` must hold 9 values.
#[no_mangle]
pub unsafe extern "C" fn net_cut_measures(
    h: *const GraphHandle,
    s: *const i64,
    ns: i64,
    t: *const i64,
    nt: i64,
    weighted: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {
    guard(err, err_len, || {
        let g = handle(h)?;
        let sv = node_indices(g, s, ns, "S")?;
        let tv = if t.is_null() {
            None
        } else {
            Some(node_indices(g, t, nt, "T")?)
        };
        let c = cut_measures(g, &sv, tv.as_deref(), weighted != 0)?;
        slice_mut(out, 9).copy_from_slice(&[
            c.cut_size,
            c.volume_s,
            c.volume_t,
            c.normalized_cut,
            c.conductance,
            c.edge_expansion,
            c.mixing_expansion,
            c.node_expansion,
            c.boundary_expansion,
        ]);
        Ok(())
    })
}
