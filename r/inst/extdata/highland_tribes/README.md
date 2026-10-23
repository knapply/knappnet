# Read's highland tribes (Gahuku-Gama)

One-mode, undirected, signed: the 16 sub-tribes of the Gahuku-Gama
alliance structure of the Eastern Central Highlands of New Guinea
described by Read (1954), with their alliance ("rova") and antagonistic
("hina") relations, as drawn by Hage and Harary (1983, pp. 56-60) to
illustrate a clusterable signed graph. This is the UCINET standard data
set GAMA (the two 16 x 16 matrices GAMAPOS and GAMANEG).

## Files

`nodes.csv` (16 rows): `name`, the tribe's name as abbreviated in the
UCINET file (Gavev, Kotun, Ove, Alika, Nagam, Gahuk, Masil, Ukudz, Notoh,
Kohik, Geham, Asaro, Uheto, Seuve, Nagad, Gama).

`edges.csv` (58 rows), one per pair with a relation:

| column | meaning |
|---|---|
| `from`, `to` | the two tribes |
| `sign` | +1 for an alliance (29 ties), -1 for an antagonistic relation (29 ties) |

No self-ties; no pair has both signs.

## Source and build

Built by `data-raw/ucinet.py` from the `tribes` object of the MIT-licensed
R package [signnet](https://github.com/schochastics/signnet) (commit
`ac3fbd028c9adf555d0b9953bd1825a431f74162`, 1 October 2026; SHA-256
recorded in the script), whose documentation names the Pajek copy of the
UCINET file (http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/gama.dat)
as its source. The UCINET dataset page
(https://sites.google.com/site/ucinetsoftware/datasets/read-highland-tribes)
describes the data.

## Licence

Published research data with no licence from their authors: the relations
are from Read (1954) and the signed graph is printed in Hage and Harary
(1983); the UCINET file is redistributed without stated terms by UCINET,
Pajek and signnet (MIT). Redistributed here with attribution on that
basis.

## Citation

Read, K. E. (1954). Cultures of the Central Highlands, New Guinea.
*Southwestern Journal of Anthropology*, 10(1), 1-43.
https://doi.org/10.1086/soutjanth.10.1.3629074

Hage, P., & Harary, F. (1983). *Structural Models in Anthropology*.
Cambridge University Press. https://doi.org/10.1017/CBO9780511659843
