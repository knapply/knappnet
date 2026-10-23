# Les Misérables character encounters

One-mode, undirected, weighted (strength): the 80 characters of Victor
Hugo's *Les Misérables* and the pairs that encounter each other in the
novel, with the number of encounters as tie strength. Knuth (1993) coded,
for each of the 356 chapters, the groups of characters who meet
("cliques of encounters") in the file `jean.dat` of the Stanford
GraphBase; this is the co-appearance network drawn by D3's force-layout
example and distributed with NetworkX and igraph.

## Files

`nodes.csv` (80 rows), one per character:

| column | meaning |
|---|---|
| `name` | the character's name as given in `jean.dat` (TeX accents converted) |
| `code` | Knuth's two-letter code |
| `description` | Knuth's one-line description (which refers to other characters by code) |

Three characters (Bruneseau, Louis Philippe, Madame Victurnien) are
mentioned but never meet anyone, so they are isolates.

`edges.csv` (254 rows), one per pair that meets:

| column | meaning |
|---|---|
| `from`, `to` | the two characters |
| `encounters` | the number of cliques of encounters that list both (a pair listed in two cliques of one chapter counts twice); 1 to 31 |

No self-ties. The weights reproduce NetworkX's `les_miserables_graph()`
exactly (the networks are isomorphic with equal weights). Knuth's own
`book("jean", ...)` graph has the same 254 ties with length 1.

## Source and build

Built by `data-raw/graphbase.py` from `jean.dat` in the Stanford GraphBase
(https://cs.stanford.edu/pub/sgb/sgb.tar.gz), read from the maintained
GitHub copy https://github.com/ascherer/sgb (commit
`bf93bd70fa1e0dbfb011214ca97c9c2cf8763edf`, 28 December 2025; SHA-256 of
`jean.dat` recorded in the script and identical to the original
distribution's).

## Licence

The Stanford GraphBase files are copyright 1993 Stanford University. The
file header says "This file may be freely copied but please do not change
it in any way!", the README says "These files may be freely copied and
distributed, provided that no changes whatsoever are made" (changed
derivatives must carry a new name and be identified as not part of the
GraphBase), and Knuth's GraphBase page says "Public-domain sources for all
programs and data of The Stanford GraphBase are freely available". The
tables here are a derived form under a new name, identified as such.

## Citation

Knuth, D. E. (1993). *The Stanford GraphBase: A Platform for
Combinatorial Computing*. Addison-Wesley.
