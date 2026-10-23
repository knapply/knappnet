# Highway mileage between North American cities

One-mode, undirected, weighted (distance), spatial (point locations): the
128 cities of the Stanford GraphBase file `miles.dat`, their coordinates
and 1980 populations, and the highway mileage between every pair of them.
Knuth (1993) took the cities (those whose names sort at or after "Ravenna,
Ohio") and mileages from the 1949 edition of Rand McNally's *Standard
Highway Mileage Guide*, corrected several dozen entries that were
obviously wrong, and adjusted about 5% of them slightly so that every
mileage satisfies the triangle inequality (the build script checks this).

## Files

`nodes.csv` (128 rows), one per city:

| column | meaning |
|---|---|
| `name` | `city, state` (unique), e.g. `San Francisco, CA` |
| `city`, `state` | the city and its two-letter state or province |
| `latitude`, `longitude` | degrees north and east (the file gives hundredths of degrees north and west; west longitudes are negative here), WGS 84 for practical purposes |
| `population` | 1980 census population |

`edges.csv` (8128 rows), one per pair of cities:

| column | meaning |
|---|---|
| `from`, `to` | the two cities |
| `miles` | highway mileage between them; a distance weight |

The table is complete (every pair). To get a sparse network, keep the
pairs within a chosen mileage or each city's nearest few, which is what
the GraphBase `miles(n, ..., max_distance, max_degree, seed)` generator
does ("if `max_distance != 0`, edges with more than `max_distance` miles
will not appear; if `max_degree != 0`, each vertex will be limited to at
most `max_degree` of its shortest edges").

## Source and build

Built by `data-raw/graphbase.py` from `miles.dat` in the Stanford
GraphBase (https://cs.stanford.edu/pub/sgb/sgb.tar.gz), read from the
maintained GitHub copy https://github.com/ascherer/sgb (commit
`bf93bd70fa1e0dbfb011214ca97c9c2cf8763edf`, 28 December 2025; SHA-256 of
`miles.dat` recorded in the script and identical to the original
distribution's). The module `gb_miles.w` documents the file.

## Licence

As for `les_miserables`: the Stanford GraphBase files are copyright 1993
Stanford University and "may be freely copied and distributed, provided
that no changes whatsoever are made"; derived files must carry a new name
and say they are not part of the GraphBase, which these tables do. Knuth's
page calls the sources and data public domain.

## Citation

Knuth, D. E. (1993). *The Stanford GraphBase: A Platform for
Combinatorial Computing*. Addison-Wesley.
