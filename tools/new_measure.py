#!/usr/bin/env python3
"""Scaffold a new node-level measure through every layer of the project.

    python3 tools/new_measure.py <area> <name> [--dry-run]

`area` is the module the measure belongs to (an existing one such as
`descriptive`, or a new one), and `name` the measure's snake_case name,
which becomes the engine function, the C symbol `net_<name>`, the C glue
entry point `C_<name>`, the R function `net_<name>`, and the manifest key. For each
layer the script either creates the file when the area is new, or prints
the snippet to paste when the file already exists. Files shared by every
area (`metrics/mod.rs`, `capi/src/lib.rs`, `netcapi.h`, `init.c`,
`NAMESPACE`, `_pkgdown.yml`, `references.bib`) are never edited; their
lines are printed. `--dry-run` prints everything and writes nothing.

The skeleton computes a stand-in value (each node's number of distinct
neighbors other than itself) at every layer, so that the whole chain
builds and the validation entry runs before the real measure is written.
Every place that must change is marked with TODO. CONTRIBUTING.md explains
the recipe the skeleton follows.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def engine_module(area: str, name: str) -> str:
    return f'''//! TODO: one line saying what this module computes.
//!
//! Sources:
//! - TODO: the paper that defines the measure (Author Year, `key`) and the
//!   paper that gives the algorithm, with their keys in docs/references.bib.
//!
//! Conventions:
//! - Self-loops: TODO (ignored | counted once | counted twice).
//! - Parallel edges: TODO (each counts | collapsed to one tie).
//! - Direction: TODO (out-ties | in-ties | ignored).
//! - Weights: TODO (ignored, defined on binary ties | distances and
//!   strengths through `crate::weights` | added as given).
//! - Normalization: TODO.
//! - Disconnected networks and undefined values: TODO (NaN where the
//!   definition has no value).
//! - Matches: TODO (igraph's `...`, NetworkX's `...`, sna's `...`).

use crate::error::Result;
use crate::graph::Graph;
use crate::parallel;

/// TODO: what the value of node `v` means, and what is returned for a
/// network without nodes. Returns one value per node, in node order.
pub fn {name}(g: &Graph, threads: usize) -> Result<Vec<f64>> {{
    let n = g.node_count();
    // TODO: replace the stand-in (distinct neighbors other than the node)
    // with the measure. `parallel::for_each_item` splits the nodes over
    // `threads` workers in a fixed pattern, so results do not depend on
    // the thread count.
    let parts =
        parallel::for_each_item(n, threads, Vec::new, |found: &mut Vec<(usize, f64)>, v| {{
            let mut count = 0usize;
            let mut last = None;
            for &w in g.out().neighbors(v) {{
                if w as usize != v && last != Some(w) {{
                    count += 1;
                }}
                last = Some(w);
            }}
            found.push((v, count as f64));
        }});
    let mut out = vec![0.0; n];
    for found in parts {{
        for (v, value) in found {{
            out[v] = value;
        }}
    }}
    Ok(out)
}}

#[cfg(test)]
mod tests {{
    use super::*;

    fn graph(n: usize, edges: &[(i64, i64)], directed: bool) -> Graph {{
        let src: Vec<i64> = edges.iter().map(|e| e.0).collect();
        let tgt: Vec<i64> = edges.iter().map(|e| e.1).collect();
        Graph::from_edges(n, &src, &tgt, None, directed).unwrap()
    }}

    #[test]
    fn hand_checked_values() {{
        // TODO: values worked out by hand from the definition.
        let g = graph(4, &[(0, 1), (1, 2), (1, 2), (2, 2)], false);
        assert_eq!({name}(&g, 1).unwrap(), vec![1.0, 2.0, 1.0, 0.0]);
    }}

    #[test]
    fn empty_network() {{
        let g = graph(0, &[], true);
        assert!({name}(&g, 1).unwrap().is_empty());
    }}
}}
'''


def engine_snippet(area: str, name: str) -> str:
    body = engine_module(area, name)
    start = body.index("/// TODO: what the value")
    end = body.index("#[cfg(test)]")
    return body[start:end]


def property_tests(area: str, name: str) -> str:
    return f'''//! Randomized property tests for {area}: invariants that must hold on
//! any network (relabeling nodes permutes results; the thread count does
//! not change them; TODO: the invariants of the measure itself).

extern crate knappnet as engine;

use engine::metrics::{area}::{name};
use engine::Graph;

/// SplitMix64 (Steele, Lea, and Flood 2014, `steele2014fast`).
struct Rng(u64);
impl Rng {{
    fn next(&mut self) -> u64 {{
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }}
    fn below(&mut self, n: u64) -> u64 {{
        self.next() % n
    }}
}}

struct Case {{
    n: usize,
    src: Vec<i64>,
    tgt: Vec<i64>,
    directed: bool,
}}

fn random_case(rng: &mut Rng) -> Case {{
    let n = 1 + rng.below(25) as usize;
    let m = rng.below((n * n / 3 + 1) as u64) as usize;
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    for _ in 0..m {{
        src.push(rng.below(n as u64) as i64);
        tgt.push(rng.below(n as u64) as i64);
    }}
    Case {{
        n,
        src,
        tgt,
        directed: rng.below(2) == 0,
    }}
}}

fn graph(c: &Case) -> Graph {{
    Graph::from_edges(c.n, &c.src, &c.tgt, None, c.directed).unwrap()
}}

fn same(a: f64, b: f64) -> bool {{
    (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}}

const CASES: u64 = 200;

#[test]
fn thread_count_does_not_change_results() {{
    let mut rng = Rng(1);
    for _ in 0..CASES {{
        let c = random_case(&mut rng);
        let g = graph(&c);
        let one = {name}(&g, 1).unwrap();
        let four = {name}(&g, 4).unwrap();
        assert!(one.iter().zip(&four).all(|(a, b)| same(*a, *b)));
    }}
}}

#[test]
fn relabeling_nodes_permutes_results() {{
    let mut rng = Rng(2);
    for _ in 0..CASES {{
        let c = random_case(&mut rng);
        // A permutation of the node labels (reversal keeps the test simple).
        let perm: Vec<i64> = (0..c.n as i64).rev().collect();
        let src: Vec<i64> = c.src.iter().map(|&v| perm[v as usize]).collect();
        let tgt: Vec<i64> = c.tgt.iter().map(|&v| perm[v as usize]).collect();
        let g = graph(&c);
        let h = Graph::from_edges(c.n, &src, &tgt, None, c.directed).unwrap();
        let a = {name}(&g, 1).unwrap();
        let b = {name}(&h, 1).unwrap();
        for v in 0..c.n {{
            assert!(same(a[v], b[perm[v] as usize]));
        }}
    }}
}}
'''


def capi_module(area: str, name: str) -> str:
    return f'''//! C ABI for {area}. TODO: list the integer codes these functions use,
//! as the other modules do.

use std::ffi::c_char;

use engine::metrics::{area}::{name};

use crate::{{guard, handle, slice_mut, GraphHandle}};

{capi_snippet(area, name)}'''


def capi_snippet(area: str, name: str) -> str:
    return f'''/// TODO: one line on the result. Writes one double per node to `out`.
///
/// # Safety
/// `h` must be a live handle; `out` must hold `n` doubles; `err` must
/// point to `err_len` writable bytes.
#[no_mangle]
pub unsafe extern "C" fn net_{name}(
    h: *const GraphHandle,
    threads: i32,
    out: *mut f64,
    err: *mut c_char,
    err_len: usize,
) -> i32 {{
    guard(err, err_len, || {{
        let g = handle(h)?;
        let values = {name}(g, threads.max(0) as usize)?;
        slice_mut(out, g.node_count()).copy_from_slice(&values);
        Ok(())
    }})
}}
'''


def header_snippet(area: str, name: str) -> str:
    return f'''/* TODO: one line on the result; out holds n doubles. */
int32_t net_{name}(const GraphHandle *h, int32_t threads, double *out,
                  char *err, size_t err_len);
'''


def c_glue(area: str, name: str) -> str:
    return f'''/*
 * R bindings for {area}. R errors are raised only after the engine call has
 * returned, so no engine memory is leaked by a longjmp (see glue.h).
 */
#include "glue.h"

{c_glue_snippet(area, name)}'''


def c_glue_snippet(area: str, name: str) -> str:
    return f'''SEXP C_{name}(SEXP ptr, SEXP threads) {{
  GraphHandle *h = graph_from_ptr(ptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)net_graph_node_count(h)));
  char err[ERR_LEN] = {{0}};
  int status = net_{name}(h, Rf_asInteger(threads), REAL(out), err, ERR_LEN);
  return vector_result(status, out, err);
}}
'''


def init_snippet(name: str) -> str:
    return f'''extern SEXP C_{name}(SEXP, SEXP);  /* in the group for its C file, sorted by name */
    {{"C_{name}", (DL_FUNC)&C_{name}, 2}},  /* in that file's group of call_methods[] */
'''


def r_function(area: str, name: str) -> str:
    title = name.replace("_", " ").capitalize()
    return f'''# {title} (man/net_{name}.Rd). Conventions and references are documented
# there; tools/check_r_docs.py checks that the usage matches this signature.
net_{name} <- function(x, threads = 0L) {{
  check_netframe(x)
  threads <- resolve_threads(threads)
  engine_call("net_{name}", .Call(C_{name}, engine_graph(x), threads))
}}
'''


def rd_file(area: str, name: str) -> str:
    title = name.replace("_", " ").capitalize()
    return f'''\\name{{net_{name}}}
\\alias{{net_{name}}}
\\title{{{title}}}
\\description{{
TODO: one sentence on what the measure is.
}}
\\usage{{
net_{name}(x, threads = 0L)
}}
\\arguments{{
\\item{{x}}{{A \\code{{\\link[=net_frame]{{netframe}}}}.}}
\\item{{threads}}{{Number of threads; \\code{{0}} uses all available cores.}}
}}
\\details{{
TODO: the conventions, in the same order as the engine module: loops,
parallel edges, direction, weights, normalization, undefined values
(\\code{{NaN}}), and which tools the result matches. If tools disagree, add a
\\code{{convention}} argument (see \\code{{\\link{{net_closeness}}}}) and one example
per convention below, with a comment naming the tool each reproduces.
}}
\\value{{A numeric vector aligned with the rows of \\code{{x$nodes}}.}}
\\references{{
TODO: Author, A. (Year). Title. \\emph{{Journal}}, volume(issue), pages.
\\doi{{10.xxxx/xxxxx}}
}}
\\examples{{
g <- net_frame(data.frame(source = c(1, 2, 2, 3), target = c(2, 3, 3, 3)),
               directed = FALSE)
net_{name}(g)
}}
'''


def r_tests(area: str, name: str) -> str:
    return f'''# {area}: hand-checked values first, then spot checks against the
# reference tools. The full comparison lives in validation/.

{r_tests_snippet(area, name)}'''


def r_tests_snippet(area: str, name: str) -> str:
    return f'''test_that("net_{name} matches hand-checked values", {{
  # TODO: values worked out by hand from the definition.
  g <- net_frame(data.frame(source = c(1, 2, 2, 3), target = c(2, 3, 3, 3)),
                 nodes = data.frame(id = 1:4), directed = FALSE)
  expect_equal(net_{name}(g, threads = 2L), c(1, 2, 1, 0))
}})

test_that("net_{name} rejects bad input", {{
  g <- net_frame(data.frame(source = 1, target = 2))
  expect_error(net_{name}(g, threads = -1L), "threads")
  expect_error(net_{name}(list()), "netframe")
}})

test_that("net_{name} agrees with igraph", {{
  skip_if_not_installed("igraph")
  # TODO: the igraph call that computes the same thing, if one exists.
  g <- net_from(igraph::make_graph("Zachary"))
  expect_length(net_{name}(g, threads = 2L), 34)
}})
'''


def manifest_file(area: str, name: str) -> str:
    return f'''# {area}: TODO one line on the measures here. Format: see ../README.md.
#
# Undefined values (NaN) are written as TODO (an empty line, or -1 when the
# measure is never negative) by the engine runner (../bindings/{area}.py)
# and by every oracle.

{manifest_snippet(area, name)}'''


def manifest_snippet(area: str, name: str) -> str:
    return f'''{name}:
  engine: {{fn: {name}}}
  applies_to: any
  tolerance: {{rel: 1.0e-12, abs: 1.0e-12}}
  oracles:
    definition:
      call: "TODO: the definition in one line ({area}_definition_oracle.py)"
    networkx:
      call: "TODO: nx.function(G, ...)"
    # igraph:
    #   call: "igraph::function(g, ...)"
    #   differs_when: [multi_edges]
    #   reason: TODO why igraph's convention differs on those fixtures.
'''


def bindings_file(area: str, name: str) -> str:
    return f'''"""Engine runners for {area} (validation/manifest.d/{area}.yaml).

Undefined values (NaN) are written as TODO, as the oracles do.
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

_declared = False


def _declare():
    """Declares the argument types of this area's C functions once."""
    global _declared
    lib = load()
    if _declared:
        return lib
    h, i32, f64p = ctypes.c_void_p, ctypes.c_int32, ctypes.POINTER(ctypes.c_double)
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_{name}.argtypes = [h, i32, f64p, cp, sz]
    lib.net_{name}.restype = i32
    _declared = True
    return lib


def _doubles(fn, handle, size, *args):
    out = (ctypes.c_double * max(size, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    if fn(handle, *args, out, err, ERR_LEN):
        raise EngineError(err.value.decode())
    return list(out)[:size]


{bindings_snippet(area, name)}'''


def bindings_snippet(area: str, name: str) -> str:
    return f'''@register("{name}")
def _{name}(g, spec, fx):
    lib = _declare()
    return _doubles(lib.net_{name}, g.handle, g.n, 0)
'''


def definition_oracle(area: str, name: str) -> str:
    return f'''#!/usr/bin/env python3
"""Golden values for {area} computed directly from the definitions.

This oracle shares no code or algorithmic idea with the engine: it is
deliberately slow and simple (exact arithmetic, brute force), so it runs
only on fixtures of at most MAX_NODES nodes.

    python3 validation/oracles/{area}_definition_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 80


def {name}(fx) -> list[float]:
    """TODO: the textbook definition, written as plainly as possible."""
    # Stand-in: distinct neighbors other than the node itself.
    adjacent = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            adjacent[a].add(b)
            if not fx.directed:
                adjacent[b].add(a)
    return [float(len(s)) for s in adjacent]


def compute(fx, engine: dict):
    if engine["fn"] == "{name}":
        return {name}(fx)
    return None


def main() -> None:
    spec = load_area("{area}")
    count = 0
    for fixture in fixture_names():
        fx = load_fixture(fixture)
        if fx.n > MAX_NODES:
            continue
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, fixture, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition ({area}): wrote {{count}} golden files")


if __name__ == "__main__":
    main()
'''


def networkx_oracle(area: str, name: str) -> str:
    return f'''#!/usr/bin/env python3
"""Golden values for {area} from NetworkX (Hagberg, Schult, and Swart
2008, `hagberg2008exploring`), as called in validation/manifest.d/{area}.yaml.

    python3 validation/oracles/{area}_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture,  # noqa: E402
                    write_values, write_version)

ORACLE = "networkx"


def build(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weights else 1.0)
    return g


def compute(fx, engine: dict):
    g = build(fx)
    if engine["fn"] == "{name}":
        # TODO: the NetworkX call named in the manifest.
        simple = nx.DiGraph(g) if fx.directed else nx.Graph(g)
        simple.remove_edges_from(nx.selfloop_edges(simple))
        return [float(simple.degree(v) if not fx.directed else simple.out_degree(v))
                for v in range(fx.n)]
    return None


def main() -> None:
    spec = load_area("{area}")
    count = 0
    for fixture in fixture_names():
        fx = load_fixture(fixture)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, fixture, variant, values)
            count += 1
    write_version(ORACLE, f"networkx {{nx.__version__}}")
    print(f"networkx {{nx.__version__}} ({area}): wrote {{count}} golden files")


if __name__ == "__main__":
    main()
'''


def bib_snippet(name: str) -> str:
    return '''@article{author2000key,
  author  = {Author, A. and Coauthor, B.},
  title   = {TODO: Title},
  journal = {TODO: Journal},
  year    = {2000},
  volume  = {1},
  number  = {1},
  pages   = {1--10},
  doi     = {10.xxxx/xxxxx}
}
'''


def plan(area: str, name: str) -> list[tuple[str, str, str]]:
    """(kind, path, content): kind is 'create' (file does not exist),
    'append' (print the snippet for an existing file), or 'shared' (a file
    the script never edits)."""
    items = [
        ("crates/engine/src/metrics/%s.rs" % area, engine_module, engine_snippet),
        ("crates/engine/tests/%s_properties.rs" % area, property_tests, None),
        ("crates/capi/src/%s.rs" % area, capi_module, capi_snippet),
        ("r/src/%s.c" % area, c_glue, c_glue_snippet),
        ("r/R/%s.R" % area, r_function, r_function),
        ("r/man/%s.Rd" % name, rd_file, None),
        ("r/tests/testthat/test-%s.R" % area, r_tests, r_tests_snippet),
        ("validation/manifest.d/%s.yaml" % area, manifest_file, manifest_snippet),
        ("validation/bindings/%s.py" % area, bindings_file, bindings_snippet),
        ("validation/oracles/%s_definition_oracle.py" % area, definition_oracle, None),
        ("validation/oracles/%s_networkx_oracle.py" % area, networkx_oracle, None),
    ]
    out = []
    for rel, whole, part in items:
        exists = (ROOT / rel).exists()
        if not exists:
            out.append(("create", rel, whole(area, name)))
        elif part is not None:
            out.append(("append", rel, part(area, name)))
        else:
            out.append(("append", rel, f"# {rel} exists; add the {name} parts by hand, "
                                       f"following the file's own pattern.\n"))
    shared = [
        ("crates/engine/src/metrics/mod.rs", f"pub mod {area};\n"),
        ("crates/capi/src/lib.rs", f"mod {area};\n"),
        ("crates/capi/include/netcapi.h", header_snippet(area, name)),
        ("r/src/init.c", init_snippet(name)),
        ("r/NAMESPACE", f"export(net_{name})\n"),
        ("r/_pkgdown.yml", f"  - net_{name}        # under the reference section the measure belongs to\n"),
        ("docs/references.bib", bib_snippet(name)),
    ]
    for rel, text in shared:
        if rel.endswith(("mod.rs", "lib.rs")) and re.search(
                rf"^(pub )?mod {re.escape(area)};$", (ROOT / rel).read_text(), re.M):
            continue
        out.append(("shared", rel, text))
    return out


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("-")]
    flags = {a for a in argv if a.startswith("-")}
    if {"-h", "--help"} & flags or len(args) != 2 or flags - {"--dry-run"}:
        print(__doc__.strip())
        return 0 if {"-h", "--help"} & flags else 2
    area, name = args
    for value, what in ((area, "area"), (name, "name")):
        if not re.match(r"^[a-z][a-z0-9_]*$", value):
            print(f"{what} must be snake_case: letters, digits, and underscores")
            return 2
    dry = "--dry-run" in flags
    for kind, rel, text in plan(area, name):
        path = ROOT / rel
        if kind == "create":
            print(f"==> {'would create' if dry else 'create'} {rel}")
            if not dry:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8", newline="\n")
        elif kind == "append":
            print(f"==> add to {rel}:\n")
            print(text)
        else:
            print(f"==> add to {rel} (shared file, not edited by this script):\n")
            print(text)
    print(f"""==> then
    cargo test --workspace --offline
    cargo clippy --workspace --all-targets --offline -- -D warnings
    ./tools/check_c_abi.sh
    gcc -std=c99 -Wall -Wextra -fsyntax-only -Ir/src r/src/{area}.c    # with R's headers
    python3 tools/check_r_docs.py
    python3 validation/check_manifest.py
    python3 validation/run_oracles.py --python-only
    python3 validation/compare.py
    python3 tools/check_bib.py --update
    python3 tools/sync_r_sources.py
    R CMD build r && R CMD check --as-cran knappnet_*.tar.gz
Every TODO marks a place that needs the real measure.""")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
