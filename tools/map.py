#!/usr/bin/env python3
"""Cross-layer index of the code base: for every exported R function, the
R file that defines it, its help page, the C glue routines it calls, the C
ABI functions those call, the Rust module that implements each, the engine
items the module uses, and the validation variants that exercise them.

    python3 tools/map.py                 # writes docs/map.json and prints a summary
    python3 tools/map.py --md FILE       # also writes a Markdown table
    python3 tools/map.py --find NAME     # every layer that mentions NAME

The index is rebuilt from the sources each time, so it never goes stale; it
is a navigation aid, not documentation, and is not committed.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")


def r_functions() -> dict[str, dict]:
    out = {}
    for f in sorted((ROOT / "r" / "R").glob("*.R")):
        text = read(f)
        for m in re.finditer(r"^([A-Za-z.][A-Za-z0-9._]*)\s*<-\s*function\s*\(", text, re.M):
            name = m.group(1)
            start = m.start()
            nxt = re.search(r"^[A-Za-z.][A-Za-z0-9._]*\s*<-\s*function\s*\(", text[m.end():], re.M)
            body = text[start: m.end() + (nxt.start() if nxt else len(text) - m.end())]
            out[name] = {
                "file": str(f.relative_to(ROOT)),
                "line": text.count("\n", 0, start) + 1,
                "calls_c": sorted(set(re.findall(r"\.Call\(\s*(C_[A-Za-z0-9_]+)", body))),
                "calls_r": sorted(set(re.findall(r"\b([A-Za-z][A-Za-z0-9._]*)\(", body))),
            }
    return out


def exports() -> set[str]:
    ns = read(ROOT / "r" / "NAMESPACE")
    return set(re.findall(r"^export\(([^)]+)\)", ns, re.M))


def rd_topics() -> dict[str, str]:
    out = {}
    for f in sorted((ROOT / "r" / "man").glob("*.Rd")):
        for a in re.findall(r"\\alias\{([^}]+)\}", read(f)):
            out[a] = str(f.relative_to(ROOT))
    return out


def c_glue() -> dict[str, dict]:
    out = {}
    for f in sorted((ROOT / "r" / "src").glob("*.c")):
        text = read(f)
        for m in re.finditer(r"^SEXP (C_[A-Za-z0-9_]+)\(", text, re.M):
            nxt = re.search(r"^SEXP C_", text[m.end():], re.M)
            body = text[m.start(): m.end() + (nxt.start() if nxt else len(text) - m.end())]
            out[m.group(1)] = {
                "file": str(f.relative_to(ROOT)),
                "abi": sorted(set(re.findall(r"\b(net_[a-z0-9_]+)\(", body))),
            }
    return out


def c_abi() -> dict[str, dict]:
    out = {}
    for f in sorted((ROOT / "crates" / "capi" / "src").glob("*.rs")):
        text = read(f)
        uses = sorted(set(re.findall(r"engine::([a-z_:]+)", text)))
        for m in re.finditer(r'pub unsafe extern "C" fn (net_[a-z0-9_]+)', text):
            nxt = re.search(r'pub unsafe extern "C" fn ', text[m.end():])
            body = text[m.start(): m.end() + (nxt.start() if nxt else len(text) - m.end())]
            calls = sorted(set(re.findall(r"\b([a-z_][a-z0-9_]*)::([a-z_][a-z0-9_]*)\(", body)))
            out[m.group(1)] = {
                "file": str(f.relative_to(ROOT)),
                "engine_calls": [f"{a}::{b}" for a, b in calls],
                "module_uses": uses,
            }
    return out


def engine_modules() -> dict[str, dict]:
    out = {}
    src = ROOT / "crates" / "engine" / "src"
    for f in sorted(src.rglob("*.rs")):
        text = read(f)
        mod = str(f.relative_to(src)).removesuffix(".rs").replace("/", "::")
        head = []
        for line in text.splitlines():
            if line.startswith("//!"):
                head.append(line[3:].strip())
            elif head:
                break
        out[mod] = {
            "file": str(f.relative_to(ROOT)),
            "summary": " ".join(head[:2])[:200],
            "pub_fns": re.findall(r"^pub fn ([a-z0-9_]+)", text, re.M),
            "cites": sorted(set(re.findall(r"`([a-z]+[0-9]{4}[a-z]+)`", text))),
        }
    return out


def manifest() -> dict[str, list[str]]:
    by_fn = defaultdict(list)
    for f in sorted((ROOT / "validation" / "manifest.d").glob("*.yaml")):
        cur = None
        for line in read(f).splitlines():
            m = re.match(r"^([A-Za-z0-9_]+):\s*$", line)
            if m:
                cur = m.group(1)
            m = re.match(r"^\s+engine:\s*\{fn:\s*([A-Za-z0-9_]+)", line)
            if m and cur:
                by_fn[m.group(1)].append(f"{f.name}:{cur}")
    return dict(by_fn)


def build() -> dict:
    rf, ex, rd, glue, abi, eng, man = (r_functions(), exports(), rd_topics(), c_glue(), c_abi(),
                                        engine_modules(), manifest())
    funcs = {}
    for name in sorted(ex):
        info = rf.get(name, {})
        cs = info.get("calls_c", [])
        # Helpers called by the function that themselves call C (one level).
        for callee in info.get("calls_r", []):
            if callee in rf and callee not in ex:
                cs = sorted(set(cs) | set(rf[callee]["calls_c"]))
        abis = sorted({a for c in cs for a in glue.get(c, {}).get("abi", [])})
        funcs[name] = {
            "r": f'{info.get("file")}:{info.get("line")}' if info else None,
            "rd": rd.get(name),
            "c_glue": {c: glue.get(c, {}).get("file") for c in cs},
            "c_abi": {a: abi.get(a, {}).get("file") for a in abis},
            "engine": sorted({e for a in abis for e in abi.get(a, {}).get("engine_calls", [])}),
            "validation": sorted({v for a in abis for v in man.get(a.removeprefix("net_"), [])}
                                 | set(man.get(name, []))),
        }
    return {"exports": funcs, "engine_modules": eng, "c_abi": abi, "c_glue": glue}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md", help="write a Markdown table of the exports here")
    ap.add_argument("--out", default=str(ROOT / "docs" / "map.json"))
    ap.add_argument("--find", help="print every layer that mentions this name")
    a = ap.parse_args()
    m = build()
    if a.find:
        hits = []
        for section, items in m.items():
            for k, v in items.items():
                if a.find in k or a.find in json.dumps(v):
                    hits.append((section, k, v))
        for section, k, v in hits:
            print(f"[{section}] {k}: {json.dumps(v)}")
        return
    Path(a.out).write_text(json.dumps(m, indent=1), encoding="utf-8")
    if a.md:
        rows = ["| R function | R | Rd | C ABI | engine | validation |", "|---|---|---|---|---|---|"]
        for k, v in m["exports"].items():
            rows.append(f"| `{k}` | {v['r']} | {v['rd']} | {', '.join(v['c_abi'])} | "
                        f"{', '.join(v['engine'][:4])} | {len(v['validation'])} |")
        Path(a.md).write_text("\n".join(rows) + "\n", encoding="utf-8")
    ex = m["exports"]
    print(f"{len(ex)} exports, {len(m['c_glue'])} C glue routines, {len(m['c_abi'])} C ABI functions, "
          f"{len(m['engine_modules'])} engine modules; "
          f"{sum(1 for v in ex.values() if not v['rd'])} exports without help, "
          f"{sum(1 for v in ex.values() if v['c_abi'] and not v['validation'])} engine-backed exports "
          f"without a validation variant -> {a.out}")


if __name__ == "__main__":
    sys.exit(main())
