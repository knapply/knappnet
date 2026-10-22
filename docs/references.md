# How the literature and other software are used

This project builds on decades of methodological work and on the reference
implementations that made network analysis widely available. This page says
where each source is used. Full citations are in
[references.bib](references.bib); keys are given in parentheses.

## Measures and algorithms

| Area | Source | Use |
|---|---|---|
| Betweenness definition | Freeman 1977 (`freeman1977set`), Freeman 1978 (`freeman1978centrality`) | Definition of shortest-path betweenness and its normalization. |
| Betweenness algorithm | Brandes 2001 (`brandes2001faster`) | O(nm) single-source accumulation used by `metrics/betweenness.rs`. |
| Betweenness variants | Brandes 2008 (`brandes2008variants`) | Successor-based accumulation without predecessor lists; future variants (endpoints, stress, load). |
| Triad types | Holland and Leinhardt 1976 (`holland1976local`) | The 16 triad types and MAN labels. |
| Triad census algorithm | Batagelj and Mrvar 2001 (`batagelj2001subquadratic`) | Subquadratic census and arc-code table used by `metrics/triad.rs`. |
| Shortest paths | Dijkstra 1959 (`dijkstra1959note`), Dial 1969 (`dial1969algorithm`) | Weighted searches; Dial's bucket queue for small integer lengths. |
| Weighted networks | Newman 2001 (`newman2001scientific`), Opsahl, Agneessens, and Skvoretz 2010 (`opsahl2010node`), Barrat et al. 2004 (`barrat2004architecture`) | Strength weights as lengths 1/w; node strength. |
| Random testing | Steele, Lea, and Flood 2014 (`steele2014fast`) | SplitMix64 generator in property tests. |
| Sparse linear algebra | Kepner et al. 2016 (`kepner2016mathematical`), Davis 2019 (`davis2019algorithm`) | Semiring formulation guiding the planned kernel layer. |
| Traversal | Beamer, Asanović, and Patterson 2012 (`beamer2012direction`) | Direction-optimizing search planned for BFS-based measures. |
| QAP (planned) | Krackhardt 1988 (`krackhardt1988predicting`), Dekker, Krackhardt, and Snijders 2007 (`dekker2007sensitivity`) | Permutation tests for dyadic regression. |
| ERGM (planned) | Snijders et al. 2006 (`snijders2006new`) | Model terms and estimation. |
| SAOM (planned) | Snijders 2001 (`snijders2001statistical`) | Actor-oriented models for network panels. |

## Reference software

These tools serve as validation oracles and design references. Their
conventions are documented in `validation/manifest.yaml`.

| Software | Citation | Role |
|---|---|---|
| igraph | `csardi2006igraph` | Oracle; default conventions for parallel edges and normalization. |
| sna | `butts2008sna` | Oracle; adjacency-matrix conventions (`diag = FALSE`). |
| network | `butts2008network` | Design reference for relational data with attributes. |
| NetworkX | `hagberg2008exploring` | Oracle; source of bundled classic datasets. |
| ergm, statnet | `hunter2008ergm`, `handcock2008statnet` | Planned oracles for ERGM. |
| UCINET | `borgatti2002ucinet` | Planned oracle (manual golden files). |
| Pajek | `batagelj1998pajek` | Planned oracle (manual golden files). |
| Gephi | `bastian2009gephi` | Planned oracle for community detection and layouts. |
| graph-tool | `peixoto2014graphtool` | Benchmark and design reference. |
| NetworKit | `staudt2016networkit` | Benchmark and design reference for parallel algorithms. |
| hellorust | `ooms_hellorust` | Pattern for linking Rust into an R package through C. |
| Arrow C Data Interface | `arrow_c_data_interface` | Planned zero-copy exchange of node and edge tables. |

## File formats

| Format | Source |
|---|---|
| GraphML | Brandes et al. 2002 (`brandes2002graphml`) |
| GML | Himsolt 1997 (`himsolt1997gml`) |
| GEXF | GEXF Working Group (`gexf2009`); Gephi (`bastian2009gephi`) |
| Pajek | Batagelj and Mrvar 1998 (`batagelj1998pajek`); de Nooy, Mrvar, and Batagelj 2018 (`denooy2018exploratory`) |
| UCINET DL | Borgatti, Everett, and Freeman 2002 (`borgatti2002ucinet`) |
| JSON node-link | NetworkX (`hagberg2008exploring`); D3 (`bostock2011d3`) |
| DOT | Gansner and North 2000 (`gansner2000open`) |

## Datasets

| Fixture | Original source |
|---|---|
| `karate` | Zachary 1977 (`zachary1977information`) |
| `florentine` | Breiger and Pattison 1986 (`breiger1986cumulated`); see also Padgett and Ansell 1993 (`padgett1993robust`) |
| `southern_women` | Davis, Gardner, and Gardner 1941 (`davis1941deep`) |
| `kite` | Krackhardt 1990 (`krackhardt1990assessing`) |
| `les_miserables` | Knuth 1993 (`knuth1993stanford`) |
