#!/usr/bin/env python3
"""Builds the country land-border network (r/inst/extdata/country_borders/)
from Natural Earth's 1:110m "Admin 0 - Countries" layer, which is in the
public domain (https://www.naturalearthdata.com/about/terms-of-use/):

    python3 data-raw/natural_earth.py <natural-earth-vector checkout>

The source is a checkout of https://github.com/nvkelso/natural-earth-vector
(only geojson/ne_110m_admin_0_countries.geojson is read; its SHA-256 is
checked against the value below, and the commit used is named in the data
set's README.md).

Nodes are the 177 country polygons of the layer (Natural Earth's de facto
view of the world, so it includes, for example, Kosovo, Somaliland,
Northern Cyprus, Taiwan and Western Sahara as separate units, and treats
Antarctica and Greenland as units). Two countries are adjacent when their
polygons share at least one boundary segment (identical consecutive vertex
pairs; the layer's polygons are topologically consistent, so a shared
border is a shared sequence of vertices). Countries that only touch at a
point (Zambia and Botswana at the Kazungula quadripoint, for instance) are
not adjacent. Each edge carries the shared border as a WKT MULTILINESTRING
in the layer's coordinate system (WGS 84 longitude and latitude, EPSG:4326)
and its length along the great-circle arcs between consecutive vertices
(`border_km`, mean Earth radius 6371.0088 km). Nodes carry Natural Earth's
label point (`longitude`, `latitude`) as their location.
"""

import csv
import hashlib
import json
import math
import os
import sys

OUT = os.path.join("r", "inst", "extdata", "country_borders")
SOURCE = os.path.join("geojson", "ne_110m_admin_0_countries.geojson")
CHECKSUM = "6866c877d39cba9c357620878839b336d569f8c662d3cfab4cb1dbe2d39c977f"
EARTH_RADIUS_KM = 6371.0088


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def great_circle_km(a, b):
    lon1, lat1 = map(math.radians, a)
    lon2, lat2 = map(math.radians, b)
    h = (math.sin((lat2 - lat1) / 2) ** 2
         + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2)
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(h))


def rings(geometry):
    if geometry["type"] == "Polygon":
        return geometry["coordinates"]
    assert geometry["type"] == "MultiPolygon", geometry["type"]
    return [ring for polygon in geometry["coordinates"] for ring in polygon]


def chains(segments):
    """Joins undirected segments (pairs of points) into polylines."""
    adj = {}
    for a, b in segments:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    unused = set(frozenset(s) for s in segments)
    lines = []
    # Start at the ends of open chains first, then pick up any closed loops.
    starts = [p for p, nb in adj.items() if len(nb) % 2 == 1] + list(adj)
    for start in starts:
        while any(frozenset((start, nb)) in unused for nb in adj[start]):
            line = [start]
            cur = start
            while True:
                nxt = next((nb for nb in adj[cur] if frozenset((cur, nb)) in unused), None)
                if nxt is None:
                    break
                unused.discard(frozenset((cur, nxt)))
                line.append(nxt)
                cur = nxt
            lines.append(line)
    assert not unused
    return lines


def fmt(value):
    return repr(value) if isinstance(value, float) else str(value)


def wkt(lines):
    return "MULTILINESTRING (%s)" % ", ".join(
        "(%s)" % ", ".join("%s %s" % (fmt(x), fmt(y)) for x, y in line) for line in lines)


def main(checkout):
    path = os.path.join(checkout, SOURCE)
    data = open(path, "rb").read()
    digest = hashlib.sha256(data).hexdigest()
    print("%s  %s" % (digest, SOURCE))
    if digest != CHECKSUM:
        raise SystemExit("checksum mismatch for %s" % path)
    features = json.loads(data.decode("utf-8"))["features"]
    nodes = []
    owner = {}  # undirected segment -> set of country codes
    for f in features:
        p = f["properties"]
        code = p["ADM0_A3"]
        nodes.append([p["ADMIN"], code, p["ISO_A3"] if p["ISO_A3"] != "-99" else "",
                      p["CONTINENT"], p["SUBREGION"], p["TYPE"], p["POP_EST"], p["POP_YEAR"],
                      p["LABEL_X"], p["LABEL_Y"]])
        for ring in rings(f["geometry"]):
            for a, b in zip(ring, ring[1:]):
                a, b = tuple(a), tuple(b)
                if a == b:
                    continue
                owner.setdefault(frozenset((a, b)), set()).add(code)
    codes = [n[1] for n in nodes]
    assert len(set(codes)) == len(codes) == 177, len(codes)
    shared = {}
    for seg, owners in owner.items():
        assert len(owners) <= 2, (seg, owners)
        if len(owners) == 2:
            shared.setdefault(frozenset(owners), []).append(tuple(seg))
    edges = []
    for pair, segments in shared.items():
        a, b = sorted(pair, key=codes.index)
        lines = chains(segments)
        km = sum(great_circle_km(p, q) for line in lines for p, q in zip(line, line[1:]))
        edges.append([a, b, round(km, 1), wkt(lines)])
    edges.sort(key=lambda e: (codes.index(e[0]), codes.index(e[1])))
    print("%d countries, %d borders" % (len(nodes), len(edges)))
    write(os.path.join(OUT, "nodes.csv"),
          ["name", "code", "iso_a3", "continent", "subregion", "type", "population",
           "population_year", "longitude", "latitude"], nodes)
    write(os.path.join(OUT, "edges.csv"), ["from", "to", "border_km", "geometry"], edges)


if __name__ == "__main__":
    if len(sys.argv) != 2 or {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0 if len(sys.argv) > 1 else 2)
    main(sys.argv[1])
