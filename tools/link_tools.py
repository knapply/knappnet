#!/usr/bin/env python3
"""Links mentions of other network-analysis tools to their websites.

In each section of a Markdown or R Markdown file (text between headings),
the first plain-text mention of each tool becomes a link. Code blocks,
inline code, headings, and existing links are left alone, so running the
script again changes nothing.

    python3 tools/link_tools.py README.md r/vignettes/*.Rmd r/vignettes/articles/*.Rmd
"""

from __future__ import annotations

import re
import sys

# Name as written in prose -> website. Longer names first so that, for
# example, "python-igraph" is matched before "igraph".
TOOLS = [
    ("python-igraph", "https://python.igraph.org/"),
    ("tidygraph", "https://tidygraph.data-imaginist.com/"),
    ("NetworKit", "https://networkit.github.io/"),
    ("NetworkX", "https://networkx.org/"),
    ("igraph", "https://r.igraph.org/"),
    ("statnet", "https://statnet.org/"),
    ("sna", "https://cran.r-project.org/package=sna"),
    ("ergm", "https://cran.r-project.org/package=ergm"),
    ("RSiena", "https://www.stats.ox.ac.uk/~snijders/siena/"),
    ("UCINET", "https://sites.google.com/site/ucinetsoftware/home"),
    ("Pajek", "http://mrvar.fdv.uni-lj.si/pajek/"),
    ("Gephi", "https://gephi.org/"),
    ("ORA", "https://www.casos.cs.cmu.edu/projects/ora/"),
    ("Graphviz", "https://graphviz.org/"),
    ("Cytoscape", "https://cytoscape.org/"),
    ("yEd", "https://www.yworks.com/products/yed"),
    ("graph-tool", "https://graph-tool.skewed.de/"),
    ("keyplayer", "https://cran.r-project.org/package=keyplayer"),
    ("pkgdown", "https://pkgdown.r-lib.org/"),
]

# Protected spans: fenced code, inline code, links/images, raw URLs, HTML tags.
PROTECT = re.compile(
    r"(```.*?```|~~~.*?~~~|`[^`\n]*`|!?\[[^\]]*\]\([^)]*\)|<[^>\n]+>|https?://\S+)",
    re.S,
)
HEADING = re.compile(r"^(#{1,6} .*|<!--.*-->)$", re.M)


def link_section(text: str) -> str:
    done: set[str] = set()
    out = []
    for i, piece in enumerate(PROTECT.split(text)):
        if i % 2 == 1:  # protected span
            out.append(piece)
            for name, _ in TOOLS:  # a tool already linked counts as done
                if piece.startswith("[") and re.search(rf"\[{re.escape(name)}\b", piece):
                    done.add(name)
            continue
        for name, url in TOOLS:
            if name in done:
                continue
            # Word boundaries; "sna" and "ORA" must stand alone, and "igraph"
            # must not be part of "python-igraph" or an R call like igraph::.
            pat = re.compile(rf"(?<![\w/-]){re.escape(name)}(?![\w:-])")
            m = pat.search(piece)
            if m:
                piece = piece[:m.start()] + f"[{name}]({url})" + piece[m.end():]
                done.add(name)
                # Re-split so the new link is protected for later names.
                sub = PROTECT.split(piece)
                piece = "".join(
                    s if k % 2 else s for k, s in enumerate(sub))
        out.append(piece)
    return "".join(out)


def process(text: str) -> str:
    # Keep YAML front matter and R chunks untouched.
    front = ""
    if text.startswith("---\n"):
        end = text.index("\n---", 4) + 4
        front, text = text[:end], text[end:]
    parts = HEADING.split(text)
    out = []
    for i, part in enumerate(parts):
        if i % 2 == 1:  # heading line (separators captured by split)
            out.append(part)
        else:
            out.append(link_section(part))
    return front + "".join(out)


def main(paths: list[str]) -> None:
    for p in paths:
        with open(p, encoding="utf-8") as fh:
            old = fh.read()
        new = process(old)
        if new != old:
            with open(p, "w", encoding="utf-8", newline="") as fh:
                fh.write(new)
            print(f"linked {p}")


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    main(sys.argv[1:])
