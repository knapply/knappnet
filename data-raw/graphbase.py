#!/usr/bin/env python3
"""Builds two example networks from the data files of Knuth's Stanford
GraphBase (Knuth 1993; the files are copyright Stanford University and
"may be freely copied but please do not change it in any way"):

    python3 data-raw/graphbase.py <sgb checkout>

The source is a checkout of https://github.com/ascherer/sgb, the maintained
copy of https://cs.stanford.edu/pub/sgb/sgb.tar.gz (the files read here are
identical in both, as their checksums below attest).

- les_miserables (jean.dat): the 77 characters of Victor Hugo's novel and
  the pairs that encounter each other. jean.dat lists, for each of the 356
  chapters, "cliques of encounters" (groups of characters who meet). Two
  characters are tied when they share a clique, and the tie's `encounters`
  value counts the cliques in which they do (a pair listed in two cliques
  of one chapter counts twice). This is the weighted form of the network
  distributed with D3 and NetworkX (`les_miserables_graph`), whose weights
  it reproduces exactly; SGB's own `book("jean", ...)` graph is the same
  set of ties with length 1.
- miles (miles.dat): 128 North American cities with 1980 populations and
  coordinates, and the highway mileage between every pair of them, adjusted
  by Knuth so that the mileages satisfy the triangle inequality. SGB's
  `miles(n, ..., max_distance, max_degree, seed)` keeps the pairs within
  `max_distance` miles or the `max_degree` nearest; the complete table is
  shipped and users choose their own threshold.
"""

import csv
import hashlib
import os
import re
import sys

OUT = os.path.join("r", "inst", "extdata")

CHECKSUMS = {
    "jean.dat": "2ebefa5a2eb6b461307d48e0445f664e657d85a4085efc89aaeb24e5e674746a",
    "miles.dat": "2888d2e90769c61abd63ba361204aaec7082b594288e0e4a41cf7151adc96542",
}


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read(sgb, name):
    path = os.path.join(sgb, name)
    data = open(path, "rb").read()
    digest = hashlib.sha256(data).hexdigest()
    print("%s  %s" % (digest, name))
    if digest != CHECKSUMS[name]:
        raise SystemExit("checksum mismatch for %s" % path)
    return data.decode("ascii")


# TeX accents used in the character descriptions of jean.dat.
TEX = [
    (r"\'e", "é"), (r"\`e", "è"), (r"\^o", "ô"), (r"\c{c}", "ç"),
    ("``", "“"), ("''", "”"), ("`", "‘"), ("'", "’"),
]


def untex(s):
    for a, b in TEX:
        s = s.replace(a, b)
    return s.replace("--", "—")


def les_miserables(sgb):
    lines = read(sgb, "jean.dat").split("\n")
    lines = [ln for ln in lines if not ln.startswith("*")]
    blank = lines.index("")
    chars = []
    for ln in lines[:blank]:
        code, rest = ln[:2], ln[3:]
        name, desc = rest.split(", ", 1)
        chars.append((code, untex(name), untex(desc)))
    codes = [c[0] for c in chars]
    assert len(chars) == 80 and len(set(codes)) == 80
    counts = {}
    chapters = 0
    for ln in lines[blank + 1:]:
        if not ln:
            continue
        chapters += 1
        if ":" not in ln:
            continue
        _, cliques = ln.split(":", 1)
        for clique in cliques.split(";"):
            members = clique.split(",")
            for i, a in enumerate(members):
                for b in members[i + 1:]:
                    assert a in codes and b in codes and a != b, ln
                    pair = tuple(sorted((a, b), key=codes.index))
                    counts[pair] = counts.get(pair, 0) + 1
    assert chapters == 356, chapters
    assert len(counts) == 254, len(counts)
    names = {c[0]: c[1] for c in chars}
    edges = sorted(counts.items(), key=lambda kv: (codes.index(kv[0][0]), codes.index(kv[0][1])))
    write(os.path.join(OUT, "les_miserables", "nodes.csv"), ["name", "code", "description"],
          [[name, code, desc] for code, name, desc in chars])
    write(os.path.join(OUT, "les_miserables", "edges.csv"), ["from", "to", "encounters"],
          [[names[a], names[b], n] for (a, b), n in edges])


def miles(sgb):
    lines = [ln for ln in read(sgb, "miles.dat").split("\n") if ln and not ln.startswith("*")]
    cities = []   # (name, state, lat, lon, population)
    dist = []     # dist[i] = distances from city i to cities 0..i-1
    pending = []
    for ln in lines:
        m = re.match(r"^(.*), ([A-Z]{2})\[(\d+),(\d+)\](\d+)$", ln)
        if m:
            if cities:
                assert len(pending) == len(cities) - 1, (cities[-1], len(pending))
                dist.append(pending)
            pending = []
            cities.append((m.group(1), m.group(2), int(m.group(3)) / 100.0,
                           -int(m.group(4)) / 100.0, int(m.group(5))))
        else:
            pending.extend(int(v) for v in ln.split())
    assert len(pending) == len(cities) - 1
    dist.append(pending)
    assert len(cities) == 128
    # miles.dat lists the cities in reverse alphabetical order; the distances
    # on the lines after city i refer to the cities listed before it, nearest
    # line first (so the first number is the distance to city i-1).
    edges = []
    for i in range(1, 128):
        for k, d in enumerate(dist[i]):
            j = i - 1 - k
            edges.append((j, i, d))
    assert len(edges) == 128 * 127 // 2
    # Knuth adjusted the data so that the triangle inequality holds.
    d = {}
    for i, j, v in edges:
        d[(i, j)] = d[(j, i)] = v
    for i in range(128):
        for j in range(128):
            if i == j:
                continue
            for k in range(128):
                if k in (i, j):
                    continue
                assert d[(i, j)] <= d[(i, k)] + d[(k, j)], (cities[i], cities[j], cities[k])
    order = sorted(range(128), key=lambda i: (cities[i][0], cities[i][1]))
    rank = {i: r for r, i in enumerate(order)}
    name = ["%s, %s" % (cities[i][0], cities[i][1]) for i in range(128)]
    assert len(set(name)) == 128
    write(os.path.join(OUT, "miles", "nodes.csv"),
          ["name", "city", "state", "latitude", "longitude", "population"],
          [[name[i], cities[i][0], cities[i][1], cities[i][2], cities[i][3], cities[i][4]]
           for i in order])
    rows = sorted((name[i], name[j], v) if rank[i] < rank[j] else (name[j], name[i], v)
                  for i, j, v in edges)
    write(os.path.join(OUT, "miles", "edges.csv"), ["from", "to", "miles"], rows)


if __name__ == "__main__":
    if len(sys.argv) != 2 or {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0 if len(sys.argv) > 1 else 2)
    les_miserables(sys.argv[1])
    miles(sys.argv[1])
