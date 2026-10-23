#!/usr/bin/env python3
"""Builds the workplace contact network (r/inst/extdata/workplace_contacts/)
from the SocioPatterns "Contacts in a workplace" data set (Génois et al.
2015), which its page distributes under the Creative Commons CC0 Public
Domain Dedication (http://www.sociopatterns.org/datasets/contacts-in-a-workplace/):

    python3 data-raw/sociopatterns.py [<download directory>]

The two source files are downloaded from sociopatterns.org into the
download directory (default: a `sociopatterns` directory next to this
script; existing files are reused) and their SHA-256 printed. Recorded
checksums, once known, go in CHECKSUMS below so that later runs verify
them; a mismatch stops the build.

Source format (quoted from the data set's page): the contact list is a
tab-separated file whose lines have the form "t i j", where t is a
timestamp in seconds since 0:00 on 24 June 2013 and i and j are the
anonymous identifiers of two people in contact; "the contact was active
during the interval [t - 20s, t]"; simultaneous contacts give several lines
with the same t. The metadata file has lines "i Di" with the department of
each person.
"""

import csv
import gzip
import hashlib
import io
import os
import sys
import urllib.request
import zipfile

OUT = os.path.join("r", "inst", "extdata", "workplace_contacts")
BASE = "https://sociopatterns.org/assets/data/"
FILES = {
    "workplace_InVS_tij.dat.zip": None,
    "workplace_InVS_metadata.txt": None,
}
CHECKSUMS = {}  # file name -> SHA-256, once recorded


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fetch(directory, name):
    path = os.path.join(directory, name)
    if not os.path.exists(path):
        os.makedirs(directory, exist_ok=True)
        print("downloading", BASE + name)
        with urllib.request.urlopen(BASE + name, timeout=120) as resp, open(path, "wb") as fh:
            fh.write(resp.read())
    data = open(path, "rb").read()
    digest = hashlib.sha256(data).hexdigest()
    print("%s  %s" % (digest, name))
    if CHECKSUMS.get(name) and digest != CHECKSUMS[name]:
        raise SystemExit("checksum mismatch for %s" % path)
    return data


def text_members(data, name):
    """Yields the text of the data file(s) inside a download."""
    if name.endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for info in z.infolist():
                if not info.is_dir() and not info.filename.startswith("__MACOSX"):
                    yield info.filename, z.read(info).decode("utf-8")
    elif name.endswith(".gz"):
        yield name[:-3], gzip.decompress(data).decode("utf-8")
    else:
        yield name, data.decode("utf-8")


def main(directory):
    contacts = []
    for member, text in text_members(fetch(directory, "workplace_InVS_tij.dat.zip"),
                                     "workplace_InVS_tij.dat.zip"):
        print("reading", member)
        for ln in text.splitlines():
            parts = ln.split()
            if not parts:
                continue
            t, i, j = (int(v) for v in parts[:3])
            assert i != j
            contacts.append((t, i, j))
    departments = {}
    for member, text in text_members(fetch(directory, "workplace_InVS_metadata.txt"),
                                     "workplace_InVS_metadata.txt"):
        for ln in text.splitlines():
            parts = ln.split()
            if parts:
                departments[int(parts[0])] = parts[1]
    people = sorted({i for _, i, j in contacts} | {j for _, i, j in contacts} | set(departments))
    missing = [p for p in people if p not in departments]
    print("%d people (%d without a department), %d contacts, %d departments" %
          (len(people), len(missing), len(contacts), len(set(departments.values()))))
    contacts.sort()
    write(os.path.join(OUT, "nodes.csv"), ["name", "department"],
          [[p, departments.get(p, "")] for p in people])
    write(os.path.join(OUT, "contacts.csv"), ["from", "to", "time"],
          [[min(i, j), max(i, j), t] for t, i, j in contacts])


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.join(os.path.dirname(os.path.abspath(__file__)), "sociopatterns"))
