#!/usr/bin/env python3
"""Builds the small example networks in r/inst/extdata/ (all but Noordin
Top, which data-raw/noordin.R builds).

    python3 data-raw/examples.py /path/to/knapply/homophily

- jemaah_islamiyah: Koschade (2006), Tables 1 and 2, as transcribed in the
  homophily package's data-raw/jemaah_islamiyah.R.
- karate, florentine, southern_women: the copies distributed with NetworkX.
"""

import csv
import os
import re
import sys

import networkx as nx

OUT = os.path.join("r", "inst", "extdata")


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def jemaah_islamiyah(homophily):
    src = open(os.path.join(homophily, "data-raw", "jemaah_islamiyah.R"), encoding="utf-8").read()

    def table(label):
        i = src.index(label)
        j = src.index('"', src.index("pasted_table_to_adj_mat(", i)) + 1
        k = src.index('"', j)
        lines = src[j:k].split("\n")
        hdr = lines[0].split("\t")
        m = {}
        for ln in lines[1:]:
            parts = ln.split("\t")
            for c, v in zip(hdr, parts[1:]):
                m[(parts[0], c)] = int(v) if v.strip() else None
        return hdr, m

    hdr, binary = table("# Table 1")
    _, strength = table("# Table 2")
    block = src.split("node_attrs <- tibble::tribble(")[1].split(")\n")[0]
    roles = dict(re.findall(r'^\s*"([A-Z]+)",\s+"([^"]+)"', block, re.M))
    write(os.path.join(OUT, "jemaah_islamiyah", "nodes.csv"), ["name", "role"],
          [[n.title(), roles[n]] for n in hdr])
    rows = []
    for i, a in enumerate(hdr):
        for c in hdr[i + 1:]:
            if binary[(a, c)]:
                # One cell of Table 2 is blank; use its symmetric counterpart.
                v = strength[(a, c)] if strength[(a, c)] is not None else strength[(c, a)]
                rows.append([a.title(), c.title(), v])
    write(os.path.join(OUT, "jemaah_islamiyah", "edges.csv"), ["from", "to", "strength"], rows)


def networkx_examples():
    g = nx.karate_club_graph()
    write(os.path.join(OUT, "karate", "nodes.csv"), ["name", "faction"],
          [[v + 1, g.nodes[v]["club"]] for v in g])
    write(os.path.join(OUT, "karate", "edges.csv"), ["from", "to", "contexts"],
          [[a + 1, b + 1, d.get("weight", 1)] for a, b, d in sorted(g.edges(data=True))])
    g = nx.florentine_families_graph()
    write(os.path.join(OUT, "florentine", "edges.csv"), ["from", "to"],
          [sorted([a, b]) for a, b in sorted(g.edges())])
    g = nx.davis_southern_women_graph()
    top = set(g.graph["top"])
    write(os.path.join(OUT, "southern_women", "edges.csv"), ["woman", "event"],
          [[a, b] if a in top else [b, a] for a, b in g.edges()])


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    jemaah_islamiyah(sys.argv[1])
    networkx_examples()
