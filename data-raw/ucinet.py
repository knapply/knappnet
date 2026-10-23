#!/usr/bin/env python3
"""Builds three of the UCINET standard data sets, read from the copies in
two MIT-licensed R packages by David Schoch (the UCINET data were hosted at
http://moreno.ss.uci.edu/data.html, which no longer exists; the UCINET
dataset pages at https://sites.google.com/site/ucinetsoftware/datasets
describe them):

    python3 data-raw/ucinet.py <networkdata checkout> <signnet checkout>

- krackhardt: Krackhardt's (1987) high-tech managers (networkdata
  `ht_advice`, `ht_friends`, `ht_reports`): 21 managers, three directed
  relations, four attributes.
- sampson: Sampson's (1968) monastery (networkdata `sampson`): 18 novices,
  ten ranked, signed relations, liking at three times.
- highland_tribes: Read's (1954) Gahuku-Gama alliance structure as drawn
  by Hage and Harary (1983) (signnet `tribes`): 16 tribes, signed ties.

The sources are git checkouts of https://github.com/schochastics/networkdata
and https://github.com/schochastics/signnet; the SHA-256 of each .rda file
read is checked against the value recorded below (the commits the data were
built from are named in each data set's README.md).
"""

import csv
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rda  # noqa: E402

OUT = os.path.join("r", "inst", "extdata")

CHECKSUMS = {
    "ht_advice.rda": "55e4d439b80e634c9ca531f57152bb4809cf9bda82b9c50be6dabd0ddadde549",
    "ht_friends.rda": "dcad3b35f79a2a4fd92019475b88ee181eb479d188321304ecd49fdc4ad71c16",
    "ht_reports.rda": "10624603f8863ec1069394d39473ece69e5f18b8ed31624a08dd42ccfa857f03",
    "sampson.rda": "f5644f75dff0ac5ef48ea2373ed63f261afdc666c8d76f30d922dca800b936b9",
    "tribes.rda": "d8f1537a7fa092ed699736f5d2eead38071e749adc4bf3558cbfe301908fa28d",
}


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def checked(path):
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    expected = CHECKSUMS.get(os.path.basename(path))
    print("%s  %s" % (digest, os.path.basename(path)))
    if expected is not None and digest != expected:
        raise SystemExit("checksum mismatch for %s (expected %s)" % (path, expected))
    return path


def load(path, name):
    return rda.load(checked(path))[name]


def krackhardt(networkdata):
    d = os.path.join(networkdata, "data")
    nodes = None
    edges = []
    for obj, relation in (("ht_advice", "advice"), ("ht_friends", "friendship"),
                          ("ht_reports", "reports_to")):
        n, e, directed = rda.igraph_tables(load(os.path.join(d, obj + ".rda"), obj))
        assert directed and len(n) == 21, (obj, len(n), directed)
        rows = [[r["id"], r["age"], r["tenure"], r["level"], r["dept"]] for r in n]
        if nodes is None:
            nodes = rows
        else:
            assert rows == nodes, "vertex attributes differ between relations"
        for r in e:
            assert r["from"] != r["to"]
            edges.append([r["from"], r["to"], relation])
    counts = {rel: sum(1 for e in edges if e[2] == rel) for rel in ("advice", "friendship", "reports_to")}
    assert counts == {"advice": 190, "friendship": 102, "reports_to": 20}, counts
    edges.sort(key=lambda r: (r[2], r[0], r[1]))
    write(os.path.join(OUT, "krackhardt", "nodes.csv"),
          ["name", "age", "tenure", "level", "department"], nodes)
    write(os.path.join(OUT, "krackhardt", "edges.csv"), ["from", "to", "relation"], edges)


# UCINET matrix name -> (relation, sign, time). The liking matrices were
# collected at three times; the other relations in the same wave as the
# third (UCINET's description of SAMPSON).
SAMPSON = {
    "SAMPLK1": ("liking", 1, 1),
    "SAMPLK2": ("liking", 1, 2),
    "SAMPLK3": ("liking", 1, 3),
    "SAMPDLK": ("disliking", -1, 3),
    "SAMPES": ("esteem", 1, 3),
    "SAMPDES": ("disesteem", -1, 3),
    "SAMPIN": ("positive_influence", 1, 3),
    "SAMPNIN": ("negative_influence", -1, 3),
    "SAMPPR": ("praise", 1, 3),
    "SAMNPR": ("blame", -1, 3),
}

# Full names for the abbreviated vertex names of the UCINET file
# ("Romul_10" is novice 10, Romuald), from the descriptions of the data in
# the UCINET and Pajek documentation.
SAMPSON_NAMES = {
    "Romul": "Romuald", "Bonaven": "Bonaventure", "Ambrose": "Ambrose",
    "Berth": "Berthold", "Peter": "Peter", "Louis": "Louis", "Victor": "Victor",
    "Winf": "Winfrid", "John": "John Bosco", "Greg": "Gregory", "Hugh": "Hugh",
    "Boni": "Boniface", "Mark": "Mark", "Albert": "Albert", "Amand": "Amand",
    "Basil": "Basil", "Elias": "Elias", "Simp": "Simplicius",
}


def sampson(networkdata):
    graphs = load(os.path.join(networkdata, "data", "sampson.rda"), "sampson")
    nodes = None
    edges = []
    for g in graphs.value:
        matrix = rda.igraph_graph_attrs(g)["name"]
        relation, sign, time = SAMPSON[matrix]
        n, e, directed = rda.igraph_tables(g)
        assert directed and len(n) == 18, (matrix, len(n), directed)
        rows = []
        for r in n:
            abbrev, number = r["name"].rsplit("_", 1)
            rows.append([SAMPSON_NAMES[abbrev], int(number)])
        if nodes is None:
            nodes = rows
        else:
            assert rows == nodes, "vertex names differ between relations"
        for r in e:
            assert r["from"] != r["to"] and r["rank"] in (1, 2, 3), r
            edges.append([nodes[r["from"] - 1][0], nodes[r["to"] - 1][0], relation, sign,
                          time, r["rank"]])
    assert {rda.igraph_graph_attrs(g)["name"] for g in graphs.value} == set(SAMPSON)
    order = [v[0] for v in SAMPSON.values()]
    number = {name: num for name, num in nodes}
    edges.sort(key=lambda r: (order.index(r[2]), r[4], number[r[0]], number[r[1]]))
    write(os.path.join(OUT, "sampson", "nodes.csv"), ["name", "number"], nodes)
    write(os.path.join(OUT, "sampson", "edges.csv"),
          ["from", "to", "relation", "sign", "time", "rank"], edges)


def highland_tribes(signnet):
    g = load(os.path.join(signnet, "data", "tribes.rda"), "tribes")
    n, e, directed = rda.igraph_tables(g)
    assert not directed and len(n) == 16 and len(e) == 58, (len(n), len(e), directed)
    names = [r["name"] for r in n]
    edges = []
    for r in e:
        a, b = sorted((r["from"], r["to"]))
        assert a != b and r["sign"] in (1, -1), r
        edges.append([names[a - 1], names[b - 1], r["sign"]])
    assert sum(1 for r in edges if r[2] == 1) == 29 and sum(1 for r in edges if r[2] == -1) == 29
    edges.sort(key=lambda r: (names.index(r[0]), names.index(r[1])))
    write(os.path.join(OUT, "highland_tribes", "nodes.csv"), ["name"], [[v] for v in names])
    write(os.path.join(OUT, "highland_tribes", "edges.csv"), ["from", "to", "sign"], edges)


if __name__ == "__main__":
    if len(sys.argv) != 3 or {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0 if len(sys.argv) > 1 else 2)
    networkdata, signnet = sys.argv[1:]
    krackhardt(networkdata)
    sampson(networkdata)
    highland_tribes(signnet)
