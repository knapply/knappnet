#!/usr/bin/env python3
"""Builds Knecht's classroom panel (r/inst/extdata/knecht/) from the
`klas12b` files distributed on the RSiena data page
(https://www.stats.ox.ac.uk/~snijders/siena/tutorial2010_data.htm, zip at
https://www.stats.ox.ac.uk/~snijders/siena/klas12b.zip):

    python3 data-raw/knecht.py [<download directory>]

The zip is downloaded into the download directory (default: a `knecht`
directory next to this script; an existing copy is reused) and its SHA-256
printed; once recorded in CHECKSUM, later runs verify it.

Redistribution: the RSiena page states no terms, and the btergm package
obtained Andrea Knecht's written permission to redistribute this excerpt
(its help page says so); the output of this script is therefore not
shipped in the package until the same permission is obtained (see
r/inst/extdata/DATASETS.md). The complete study data are archived at DANS
(doi:10.17026/dans-z9b-h2bp).

Source files and codes (quoted from the data page): friendship networks in
adjacency format, "klas12b-net-1.dat" to "klas12b-net-4.dat", coded "0 =
no, 1 = yes, 9 = missing, 10 = not a member of the classroom (structural
zero)"; "klas12b-delinquency.dat" (waves 1-4) and "klas12b-alcohol.dat"
(waves 2-4), "1 = never, 2 = once, 3 = 2-4 times, 4 = 5-10 times, 5 = more
than 10 times; 0 = missing"; "klas12b-demographics.dat" with sex (1 = girl,
2 = boy), age in years, ethnicity (1 = Dutch, 2 = other, 0 = missing) and
religion (1 = Christian, 2 = non-religious, 3 = non-Christian religion, 0 =
missing); "klas12b-advice.dat" (4 = low, 8 = high, 0 = missing);
"klas12b-primary.dat", whether two pupils attended the same primary school
(0 = no, 1 = yes); and "klas12b-present.dat", the changing composition of
the class (one column per wave, 1 = present).

Output, in the long-table layout of docs/DECISIONS.md (items 9, 14, 19):
nodes.csv (one row per pupil with the static attributes), edges.csv (one
row per friendship nomination with its wave, plus rows flagged `missing`
for pairs coded 9), node_spells.csv (presence spells in wave units: a pupil
present in waves 1-2 has start 1 and end 3), node_changes.csv (the
wave-specific attributes and whether the pupil answered the network
question in that wave) and same_primary_school.csv.
"""

import csv
import hashlib
import io
import os
import sys
import urllib.request
import zipfile

OUT = os.path.join("r", "inst", "extdata", "knecht")
URL = "https://www.stats.ox.ac.uk/~snijders/siena/klas12b.zip"
CHECKSUM = None  # SHA-256 of klas12b.zip, once recorded
N = 26


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fetch(directory):
    path = os.path.join(directory, "klas12b.zip")
    if not os.path.exists(path):
        os.makedirs(directory, exist_ok=True)
        print("downloading", URL)
        with urllib.request.urlopen(URL, timeout=120) as resp, open(path, "wb") as fh:
            fh.write(resp.read())
    data = open(path, "rb").read()
    digest = hashlib.sha256(data).hexdigest()
    print("%s  klas12b.zip" % digest)
    if CHECKSUM and digest != CHECKSUM:
        raise SystemExit("checksum mismatch for %s" % path)
    return zipfile.ZipFile(io.BytesIO(data))


def matrix(z, name):
    """Reads a whitespace-separated numeric table from the zip (the file may
    sit in a subdirectory)."""
    member = next((m for m in z.namelist() if os.path.basename(m) == name), None)
    if member is None:
        raise SystemExit("%s is not in the zip (members: %s)" % (name, ", ".join(z.namelist())))
    rows = []
    for ln in z.read(member).decode("utf-8", errors="replace").splitlines():
        parts = ln.replace(",", " ").split()
        if parts:
            rows.append([int(float(v)) for v in parts])
    widths = {len(r) for r in rows}
    assert len(widths) == 1, (name, widths)
    print("%s: %d x %d" % (name, len(rows), widths.pop()))
    return rows


def na(value, missing=(0,)):
    return "" if value in missing else value


def main(directory):
    z = fetch(directory)
    nets = [matrix(z, "klas12b-net-%d.dat" % w) for w in (1, 2, 3, 4)]
    demographics = matrix(z, "klas12b-demographics.dat")
    delinquency = matrix(z, "klas12b-delinquency.dat")
    alcohol = matrix(z, "klas12b-alcohol.dat")
    advice = matrix(z, "klas12b-advice.dat")
    primary = matrix(z, "klas12b-primary.dat")
    present = matrix(z, "klas12b-present.dat")
    for m in nets + [primary]:
        assert len(m) == N and len(m[0]) == N
    for m in (demographics, delinquency, alcohol, advice, present):
        assert len(m) == N, len(m)
    assert len(demographics[0]) == 4 and len(delinquency[0]) == 4
    assert len(alcohol[0]) == 3 and len(present[0]) == 4
    names = ["s%02d" % (i + 1) for i in range(N)]

    nodes = []
    for i in range(N):
        sex, age, ethnicity, religion = demographics[i]
        nodes.append([names[i], {1: "girl", 2: "boy"}[sex], age,
                      na({0: 0, 1: "Dutch", 2: "other"}[ethnicity]),
                      na({0: 0, 1: "Christian", 2: "non-religious",
                          3: "non-Christian religion"}[religion])])
    if len(advice[0]) == 1:
        for i in range(N):
            nodes[i].append(na(advice[i][0]))
        node_header = ["name", "sex", "age", "ethnicity", "religion", "advice"]
    else:
        node_header = ["name", "sex", "age", "ethnicity", "religion"]

    edges = []
    responded = {}
    for w, net in enumerate(nets, start=1):
        for i in range(N):
            row = [net[i][j] for j in range(N) if j != i]
            # A pupil who was in the class but did not answer has 9 in the
            # whole row (a non-respondent, decision 14); a pupil not in the
            # class has 10 in the whole row and column.
            responded[(i, w)] = any(v in (0, 1) for v in row)
            for j in range(N):
                if i == j:
                    continue
                v = net[i][j]
                assert v in (0, 1, 9, 10), (w, i, j, v)
                if v == 1:
                    edges.append([names[i], names[j], w, "FALSE"])
                elif v == 9 and responded[(i, w)]:
                    edges.append([names[i], names[j], w, "TRUE"])
    print("%d nominations, %d missing pairs among respondents" %
          (sum(1 for e in edges if e[3] == "FALSE"), sum(1 for e in edges if e[3] == "TRUE")))

    spells = []
    for i in range(N):
        w = 1
        while w <= 4:
            if present[i][w - 1] == 1:
                start = w
                while w <= 4 and present[i][w - 1] == 1:
                    w += 1
                spells.append([names[i], start, w])
            else:
                w += 1
    changes = []
    for i in range(N):
        for w in (1, 2, 3, 4):
            changes.append([names[i], w, na(delinquency[i][w - 1]),
                            na(alcohol[i][w - 2]) if w >= 2 else "",
                            "TRUE" if responded[(i, w)] else "FALSE"])
            if len(advice[0]) == 4:
                changes[-1].append(na(advice[i][w - 1]))
    change_header = ["name", "wave", "delinquency", "alcohol", "responded"]
    if len(advice[0]) == 4:
        change_header.append("advice")
    same_school = [[names[i], names[j]] for i in range(N) for j in range(i + 1, N)
                   if primary[i][j] == 1]
    print("%d presence spells, %d same-primary-school pairs" % (len(spells), len(same_school)))

    write(os.path.join(OUT, "nodes.csv"), node_header, nodes)
    write(os.path.join(OUT, "edges.csv"), ["from", "to", "wave", "missing"], edges)
    write(os.path.join(OUT, "node_spells.csv"), ["name", "start", "end"], spells)
    write(os.path.join(OUT, "node_changes.csv"), change_header, changes)
    write(os.path.join(OUT, "same_primary_school.csv"), ["from", "to"], same_school)


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.join(os.path.dirname(os.path.abspath(__file__)), "knecht"))
