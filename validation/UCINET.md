# Producing the UCINET golden files

UCINET (Borgatti, Everett, and Freeman 2002) is a Windows program driven
through menus, so it cannot run in CI. Its results are checked in by hand
as golden files under `validation/golden/ucinet/`, and until they exist
`validation/compare.py` reports every UCINET comparison as `manual`. This
document is the complete procedure for producing them. It takes about two
hours the first time.

You need:

- a Windows machine with UCINET 6.849 or later (the engine's UCINET
  conventions follow 6.849; older builds use different maxima for directed
  degree centralization, see `crates/engine/src/conventions.rs`);
- a clone of this repository, with Python 3 and PyYAML (the same Python
  that runs `validation/compare.py`).

## 1. What has to be produced

One file per fixture and measure:

```
validation/golden/ucinet/<fixture>/<measure>.csv
validation/golden/ucinet/VERSION
```

`<fixture>` is a directory name under `validation/fixtures/` and
`<measure>` a manifest key with a `ucinet:` oracle entry. Each `.csv` file
is the line `value` followed by one number per line, exactly what
`common.write_values()` writes and `common.read_values()` reads:

- node-level measures: one line per node, in node order 0, 1, ..., n - 1;
- network-level measures: a single line;
- partitions (`compare: partition`): one line per node holding the node's
  group number; the numbers are arbitrary labels, only which nodes share
  one matters.

A missing value is written as `-1` for every measure in
`manifest.d/topography.yaml` (its header says so: none of those measures
can be negative, and the engine writes undefined results as -1 too), and
as `NA` elsewhere. Section 5 writes these files for you from plain text
files you fill in, so you never edit the CSVs by hand.

The measures with a `ucinet:` entry at the time of writing (all in
`validation/manifest.d/topography.yaml`) are listed below. Do not rely on
this list: step 2 prints the authoritative checklist from the manifests,
including any entry added later. If a new entry appears that section 4
does not describe, follow its `call:` text in the manifest the same way.

| Measure | Kind | Fixtures |
|---|---|---|
| `density_simple` | network | all |
| `components_weak` | partition | all |
| `components_strong` | partition | all |
| `fragmentation` | network | all |
| `closeness_ucinet` | node | all |
| `closeness_ucinet_in` | node | directed only |
| `harmonic_normalized` | node | all |
| `eigenvector_ucinet` | node | all |
| `centralization_degree_ucinet` | network | all |
| `centralization_betweenness` | network | all |
| `centralization_closeness_ucinet` | network | all |

All of these are defined on binary ties of the adjacency matrix (the
manifest entries use the engine's `simple`, `symmetrize`, or `ucinet`
convention, or are unaffected by weights), so UCINET is given each network
as a 0/1 matrix with an empty diagonal: parallel edges collapsed and
self-loops dropped, which is also what UCINET itself does with a
multigraph. Weights are not used.

## 2. Write the UCINET input files

From the repository root, save this as a file outside the repository (for
example `C:\ucinet-knappnet\make_inputs.py`) and run it with the folder to
use as its argument:

```
python C:\ucinet-knappnet\make_inputs.py C:\ucinet-knappnet
```

```python
"""Writes one UCINET DL file per fixture and the checklist of measures."""
import sys
from pathlib import Path

sys.path.insert(0, "validation")
import common

out = Path(sys.argv[1])
(out / "dl").mkdir(parents=True, exist_ok=True)
manifest = common.load_manifest()
lines = []
for name in common.fixture_names():
    fx = common.load_fixture(name)
    a = [[0] * fx.n for _ in range(fx.n)]
    for s, t in zip(fx.src, fx.tgt):
        if s != t:
            a[s][t] = 1
            if not fx.directed:
                a[t][s] = 1
    labels = ",".join(f"n{i}" for i in range(fx.n))
    rows = "\n".join(" ".join(map(str, r)) for r in a)
    (out / "dl" / f"{name}.txt").write_text(
        f"dl n={fx.n} format=fullmatrix\nlabels:\n{labels}\ndata:\n{rows}\n")
    for variant, spec in manifest.items():
        if "ucinet" in spec["oracles"] and common.applies(spec, fx):
            kind = "directed" if fx.directed else "undirected"
            lines.append(f"{name:<24} {kind:<10} n={fx.n:<3} {variant}")
(out / "checklist.txt").write_text("\n".join(lines) + "\n")
print(f"{len(lines)} results to produce; see {out / 'checklist.txt'}")
```

This writes `dl\<fixture>.txt` for every fixture and `checklist.txt`, one
line per result to produce (fixture, whether it is directed, its size,
measure). Nodes are labelled `n0`, `n1`, ...: row `nK` of every UCINET
output is node K, which is the order the golden files need.

## 3. Import the networks into UCINET

1. Start UCINET. Set the default folder (File > Change Default Folder, or
   the folder box on the toolbar) to `C:\ucinet-knappnet\dl`.
2. For each fixture: Data > Import text file > DL, choose
   `<fixture>.txt`, and accept the default output dataset name
   (`<fixture>`). UCINET writes `<fixture>.##h` and `<fixture>.##d`.
3. For each **directed** fixture (the checklist says `directed`), also
   make the symmetrized copy needed for fragmentation: Transform >
   Symmetrize, input `<fixture>`, method **Maximum**, output dataset
   `<fixture>-sym`.

Check the import once: open `karate` in the matrix editor (Data > Data
editors > Matrix editor) and confirm it is 34 x 34, symmetric, with 156
ones (78 ties) and labels `n0` to `n33`.

## 4. Run each measure

Record results in plain text files under `C:\ucinet-knappnet\results\`,
one per checklist line: `results\<fixture>\<measure>.txt`, holding the
numbers one per line, with nothing else. For node-level measures that is n
lines, in the order of the output rows (`n0` first), and for network-level
measures a single line. Copy numbers from the output log (select the
column, copy, paste) or from the routine's saved output dataset opened in
the matrix editor, which shows more digits; use every digit UCINET shows
(see section 7). Where this document says to write `undefined`, write that
word on the line instead of a number.

General rules:

- Always record what UCINET reports, even where it looks wrong. A
  difference from the engine is for the maintainers to judge (it either is
  a documented convention difference or a bug in one of the two programs);
  do not adjust values to match.
- Write `undefined` only when UCINET reports no value at all (blank,
  missing, `*`, or an error for that network).
- Some routines report proportions and others percentages. Every golden
  value is a proportion: where the output is a percentage, divide by 100.
  The checks at the end of each step tell you which one you have.
- Use the dataset `<fixture>` (not `-sym`) unless the step says otherwise.

### density_simple

Network > Cohesion > Density > Density Overall, input `<fixture>`. If the
dialog offers a diagonal option, choose to exclude (ignore) the diagonal.
Record the overall density (one number). Check: `complete6` and
`complete5_dir` give 1, `empty` gives 0. For `single_node` there is no
possible tie; record what UCINET shows, or `undefined` if it shows
nothing.

### components_weak and components_strong

Network > Regions > Components, input `<fixture>`. For directed data choose
**Weak** components for `components_weak` and **Strong** components for
`components_strong`; for undirected data the two are the same, so run it
once and record the same membership in both files. Record, for every node
in order, the number of the component it belongs to (from the membership
or partition the routine outputs, one number per node). Isolated nodes are
components of their own and must each get a number. Check: `pair` (a
single arc n0 -> n1) has one weak component and two strong components.

### fragmentation

Network > Regions > Components, input `<fixture>` for undirected data and
`<fixture>-sym` for directed data (fragmentation here is weak: pairs with
no path in either direction). Record the reported **fragmentation**; if
the routine labels it "normalized heterogeneity" instead, record that
value. Check: `complete6` gives 0, `empty` gives 1, and
`two_triangles_isolate` (two triangles and an isolate, 7 nodes) gives
1 - 2 * 6 / 42 = 0.7142857... If Components does not report it in your
build, use Network > Cohesion > Distance on the same dataset and record its
fragmentation value.

### closeness_ucinet

Network > Centrality > Closeness, input `<fixture>`, with:

- type: Freeman closeness (sum of geodesic distances);
- "Set undefined distances to": **Max observed distance plus 1**;
- "Divide totals into N-1" (normalized), if offered.

Record one value per node: column `FreeClo` (named `OutClose`, or the
out-closeness column, for directed data). Check: in `complete6` every node
gets 1 (or 100, a percentage: divide by 100); in `star8` the hub gets 1 and
each leaf 7/13 = 0.5384615... For `single_node` record what is shown or
`undefined`.

### closeness_ucinet_in (directed fixtures only)

The same run as `closeness_ucinet`; record the in-closeness column
(`InClose`). Check: `cycle7_dir` gives every node 6/21 = 0.2857142...

### harmonic_normalized

The same Closeness routine, type **reciprocal** closeness (sum of
reciprocal distances), with **Averages** (divide by n - 1). Record column
`RecipClo` (`OutRecipClo`, the out-direction column, for directed data).
Check: `complete6` gives 1 everywhere; `star8` gives the hub 1 and each
leaf (1 + 6/2) / 7 = 0.5714285... Nodes that reach no other node get 0.

### eigenvector_ucinet

Network > Centrality > Eigenvector, input `<fixture>`. UCINET symmetrizes
directed data by maximum itself; do not symmetrize first. Record the
**unscaled** eigenvector column (`Eigenvec`), not the normalized one
(`nEigenvec`). Check: the squares of the recorded values sum to 1 (to the
digits shown), and in `complete6` every node gets 1/sqrt(6) = 0.4082482...
Record the values as positive numbers (if UCINET shows them all negative,
drop the signs). For `empty`, `single_node`, and other networks without
ties, record what is shown, or `undefined` for each node if nothing is.
For networks with more than one component that has ties, UCINET reports
lambda-squared scores computed within components; record them as shown
(the manifest documents this difference).

### centralization_degree_ucinet

Network > Centrality > Degree, input `<fixture>`; for undirected data make
sure the data are treated as symmetric (the dialog's "Treat data as
symmetric" set to Yes, or Auto-detect). Record "Network Centralization"
for undirected data and "Network Centralization (Outdegree)" for directed
data. Check: `star8` gives 1 (or 100%: divide by 100), `complete6` 0.

### centralization_betweenness

Network > Centrality > Freeman Betweenness > Node Betweenness, input
`<fixture>`. Record "Network Centralization Index" divided by 100 (it is
a percentage). Check: `star8` gives 1, `cycle7_dir` 0.

### centralization_closeness_ucinet

The legacy Closeness routine (in 6.849, Network > Centrality > Closeness
(old version), or the routine labelled "legacy"), input `<fixture>`.
Record "Network Centralization" (for directed data the out-centralization)
divided by 100. UCINET reports it only when every node reaches every
other: if no centralization is shown, write `undefined`. Check: `star8`
gives 1, `complete6` 0, `two_triangles_isolate` `undefined`.

## 5. Write the golden files

When every checklist line has a results file, save this as
`C:\ucinet-knappnet\collect.py` and run it from the repository root with
the work folder and the exact UCINET version (Help > About, for example
`6.849 (2026-03-01)`):

```
python C:\ucinet-knappnet\collect.py C:\ucinet-knappnet "6.849 (2026-03-01)"
```

```python
"""Turns results\\<fixture>\\<measure>.txt into golden files."""
import sys
from pathlib import Path

sys.path.insert(0, "validation")
import common
import yaml

work, version = Path(sys.argv[1]), sys.argv[2]
# Manifests whose header says undefined results are written as -1.
minus_one = set()
for f in sorted(Path("validation/manifest.d").glob("*.yaml")):
    text = f.read_text()
    if "written as -1" in text:
        minus_one |= set((yaml.safe_load(text) or {}).keys())
manifest = common.load_manifest()
problems, written = [], 0
for name in common.fixture_names():
    fx = common.load_fixture(name)
    for variant, spec in manifest.items():
        if "ucinet" not in spec["oracles"] or not common.applies(spec, fx):
            continue
        path = work / "results" / name / f"{variant}.txt"
        if not path.exists():
            problems.append(f"missing {path}")
            continue
        words = path.read_text().split()
        # The definition oracle's golden file has the expected length
        # (n for node-level measures, 1 for network-level ones).
        expected = common.read_values("definition", name, variant)
        if expected is not None and len(words) != len(expected):
            problems.append(f"{path}: {len(words)} values, expected {len(expected)}")
            continue
        values = []
        for w in words:
            if w.lower() == "undefined":
                values.append(-1.0 if variant in minus_one else float("nan"))
            else:
                values.append(float(w.replace(",", "")))
        common.write_values("ucinet", name, variant, values)
        written += 1
(common.GOLDEN / "ucinet" / "VERSION").write_text(f"UCINET {version}\n")
print(f"wrote {written} golden files")
for p in problems:
    print("PROBLEM:", p)
```

Fix every `PROBLEM` line (a missing file, or a node-level result with the
wrong number of lines) and run it again. NaN values (for measures outside
`topography.yaml`) are written as `nan`, which `common.read_values()`
reads as undefined.

## 6. Check and commit

```
python3 validation/compare.py
```

The UCINET column of `validation/REPORT.md` should now show matches, and
no `manual` entries. For every `FAIL`:

1. confirm the value was copied from the right column and row, and that
   percentages were divided by 100;
2. if the value is right, keep it: it is a real difference. Report it
   (the fixture, measure, both values, and the UCINET output) rather than
   editing the golden file. The maintainers will either fix the engine or
   document the convention in the manifest (`differs_when` with a reason).

Then commit `validation/golden/ucinet/` (all fixture directories and
`VERSION`) and, for traceability, the UCINET output logs as
`validation/golden/ucinet/logs/<fixture>.txt` (File > Save of the output
log window, or copy and paste it). Do not commit the DL files or results
folder: `make_inputs.py` regenerates the inputs.

## 7. Precision

The manifest tolerances (as small as 1e-12 relative) assume full double
precision. Values copied from UCINET's output log carry only the digits
shown (often 3 or 4 decimals), so they will fail those comparisons even
when UCINET and the engine agree. Take values from the saved output
datasets (opened in the matrix editor, or exported with Data > Export to
a text or Excel file), which keep the stored precision, and record every
digit. If a comparison fails only by rounding at the last digit shown,
note the number of digits available in the commit message; the manifest
entry then needs a tolerance that reflects UCINET's precision.
