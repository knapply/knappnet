# Roadmap

What the package still needs in order to replace the workflows of the
teaching sources it is meant to supersede, consolidated from three surveys
of those sources (October 2026): the two online R books
(`rawlings2023network`, chapters cited as R3-R15; `shafie2026r4sna`,
S1-S23), the three tutorials (`wiederstein2021network`, W;
`almeidasilva2024network`, NAR; the qgraph cookbook, `epskamp2012qgraph`,
E), and the three textbooks (`borgatti2022analyzing`, A; `kadushin2012understanding`,
K; `cunningham2016understanding`, U). Each item says what it is, why
analysts need it (which source and chapter teaches it), the published
source to implement from, and a tier:

- **now**: the data-model changes the owner decided in October 2026
  ([DECISIONS.md](DECISIONS.md)), marked *decided*, and the gaps an
  analyst following any of the sources hits in the first sessions;
- **next**: standard descriptive methods the sources treat as routine;
- **later**: specialised methods and strategic additions;
- **out of scope**: things the package will document or hand off rather
  than implement.

Sources with a BibTeX key are in `docs/references.bib`; the others are
named by author and year and must be verified and added when the item is
implemented. Reference tools under copyleft licences ([igraph](https://r.igraph.org/), [sna](https://cran.r-project.org/package=sna),
network, [statnet](https://statnet.org/), [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/), egor, backbone, signnet, netrankr, blockmodeling)
are oracles only; everything is implemented from the published source
(CONTRIBUTING.md, "Licensing").

ERGM, TERGM, and SAOM are not on this roadmap: they live on the
`ergm-saom` branch until their API is redesigned (decision 13), and the
bridges below are the supported path meanwhile.

## Data model and validation

| What | Why | Implement from | Tier |
|---|---|---|---|
| `net_` prefix on every exported function (*decided*, 1) | Name clashes with [sna](https://cran.r-project.org/package=sna) and [igraph](https://r.igraph.org/); discoverability | stringr's naming | now |
| Printed card plus table previews (*decided*, 2) | W1 exists only to decode igraph's header; every survey asks for weight meaning, two-mode, temporal, spatial, signed status and loop and parallel-edge counts in the print | [tidygraph](https://tidygraph.data-imaginist.com/)'s preview for the layout | now |
| S3 subclasses for spatial, temporal, and spatio-temporal networks (*decided*, 3) | Dispatch for print, plot, and summary | | now |
| Engine validation in four groups with listed rows and advisory notes, `quiet = TRUE` (*decided*, 4) | W7, NAR appendix, R3: endpoints missing from the node table, misaligned attributes, 0-based ids cost readers whole sections | `butts2008network` (`as.network.data.frame()` checks) | now |
| `weight_type` required (*decided*, 5) | R11 and E: weights read as strengths by every source, as distances by igraph's path functions | `newman2001scientific`, `opsahl2010node` | now |
| `mode` node role of any discrete type, `type` accepted, separate two-mode and bipartite tests (*decided*, 6) | R11, S6, A13, U3: two-mode data; the DyNetML `type` collision | `wasserman1994social` ch. 8, `borgatti1997network`, `newman2018networks` ch. 6 | now |
| Role-aware parallel edges; `net_is_simple()`, `net_simplify(combine =)`, edge multiplicity (*decided*, 7) | NAR1, W8, S1 call igraph's `simplify()`; U3 relation tables repeat pairs | `newman2018networks` ch. 6 | now |
| `relation` edge role; `relation =` on measures; split helper (*decided*, 8) | U every chapter (pick or combine relations); A5.4.6; R7, R10 handle relations by hand | `kivela2014multilayer`, `battiston2014structural`, `hanneman2005introduction` ch. 16 | now |
| Explicit `source`/`target` endpoint arguments (*decided*, 10) | W1, W7, NAR1: position-based meaning | | now |
| Loop helpers (`net_has_loops()`, `net_remove_loops()`) | R3, NAR1 drop loops before measures | | next |
| `graph$design` metadata (boundary, sampling, instrument, fixed choice) | A3, K2, U3: research design recorded with the data | `marsden1990network`, `handcock2010modeling` | next |
| Node and edge editing with validation: `net_add_nodes()`, `net_add_edges()`, `net_filter_edges()`, `net_mutate_nodes()`, `net_mutate_edges()`, `net_edge_attr_from_nodes()` | NAR1.3, W3, S21-23: daily attribute editing and pruning; the `.N()[from]` trick in S22-23 and R8 | | next |
| Multilayer subclass with node-layers, `as_multilayer()`, supra-adjacency (*decided* as later, 8) | Multilayer measures when they are added | `kivela2014multilayer`, `dedomenico2013mathematical`, `magnani2021analysis` | later |
| Long `node_spells`/`edge_spells` tables for several node spells and time-varying attributes | R3 (TEAs in networkDynamic); [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/) behaviour variables | `butts2026networkdynamic`, `latapy2018stream` | later |
| Hypergraphs and simplicial complexes (*decided*, 12) | Documented two-mode mapping only | `battiston2020networks`, `torres2021why` | out of scope |

## Data management transforms

| What | Why | Implement from | Tier |
|---|---|---|---|
| `net_dyad_census()` and attribute-conditioned dyad and triad censuses with named MAN output | R4, R7, R13, S3, S16, NAR2: both R books start inference from the dyad census and hand-roll the attribute split | `holland1976local`; Holland and Leinhardt 1970; `wasserman1994social` ch. 13-14 | now |
| Shortest-path enumeration (node and edge sequences, all geodesics) | R3, S3, S21: a path is shown to the reader, not just its length | Brandes 2001 predecessor lists (`brandes2001faster`) | now |
| Vertex connectivity (pairwise and global) and cohesive blocking | R8: Moody and White cohesion | Even and Tarjan 1975; Moody and White 2003 | next |
| `net_normalize_ties()` (row, column, matrix; sum, mean, max, z, Euclidean, iterative CA) | A5.5.6, U3 | [UCINET](https://sites.google.com/site/ucinetsoftware/home) 6 reference manual (`borgatti2002ucinet`); `hanneman2005introduction` ch. 6 | next |
| `net_attribute_matrix()` (same, diff, absdiff, sum, product, sender, receiver) | A5.6, A14.5: inputs to QAP and MRQAP | `borgatti2002ucinet` | next |
| Full `symmetrize` method set (difference, product, division, comparisons) and a missing-value rule; `dichotomize(density =)` and a diagonal policy | A5.5.1, A5.5.4, U3 | `borgatti2002ucinet` | next |
| Walk counts, `net_neighbors()`, per-edge mutual flag | R3, S3, NAR2 | `wasserman1994social` ch. 4 | next |
| Deterministic generators (empty, complete, star, ring, lattice, tree) | S3, S4: reference graphs for density and centralization | [igraph](https://r.igraph.org/)'s `make_*` for conventions | next |
| Line graph | NAR2: edge-level measures | Harary and Norman 1960 | later |
| Multiplex coding of several relations into one coded matrix; sparsification (L-Spar); inverse-weighted degree | A5.9, A5.5.2 (2e) | Satuluri, Parthasarathy, and Ruan 2011 | later |
| Degree-preserving rewiring null (`net_rewire()`), connected simple configuration graphs | NAR3, S16: the degree-matched null as the informative one | Maslov and Sneppen 2002; Viger and Latapy 2005; Milo et al. 2003 | next |
| Isomorphism test, literal constructor, correspondence analysis | NAR1, NAR3, A6.3, R12 | Greenacre 2017 for CA | out of scope (CA: point to FactoMineR; isomorphism later) |

## Ego networks

| What | Why | Implement from | Tier |
|---|---|---|---|
| Egocentric survey data: an `ego` role (*decided*, article) and import from and export to the three-table shape (ego, alter, alter-alter tie) | R6, S8, A8: ego-network chapters are mostly data wrangling | `burt1984network`; `perry2018egocentric`; egor's data model as a reference for conventions only | now |
| `net_ego_composition()`: shares and means of alter attributes, Blau and Shannon heterogeneity, IQV, proportion matching ego, ego-alter mixing table with an odds ratio against random pairing, a `ties` selector (any, out, in, reciprocated) | A8.2-8.5, R6, U7, K2 | `perry2018egocentric` ch. 6-7; `krackhardt1988predicting` (node E-I); Blau 1977 | now |
| `net_ego_similarity()` for continuous attributes with additive, ratio, and interval normalisations | A8.5 | `perry2018egocentric` | next |
| Structural holes in [UCINET](https://sites.google.com/site/ucinetsoftware/home)'s ego-network model, pendant and isolate policies, dyadic redundancy and constraint tables | A8.6.1, U7: the book's and UCINET's default differs from [igraph](https://r.igraph.org/)'s whole-network model | `burt1992structural` ch. 2; Borgatti 1997 (Connections 20(1)) | next |
| Ego-network gaps in `ego_summary`: 2-step reach, reach efficiency, brokerage pairs, node-level reciprocity | A8.6.2, U7 | `borgatti2022analyzing` | next |
| Wide one-file survey to long tables | R6 | | later |

## Centrality

| What | Why | Implement from | Tier |
|---|---|---|---|
| One-call node-measure table (`net_centrality_table()`) with raw and normalised columns, and `net_top_nodes()` with tie handling | NAR2, E, A9, U6: the first thing every analyst builds, by hand | `freeman1978centrality` | now |
| Opsahl's alpha-tuned degree and closeness | R11: the tnet detour | `opsahl2010node` | next |
| Subgraph centrality, stress centrality | S4 | Estrada and Rodriguez-Velazquez 2005; Shimbel 1953 | next |
| Induced centrality (a node's contribution to any graph index) | A9.7 | Everett and Borgatti 2010 | next |
| PN centrality (needs signed ties) | A9.6 | Everett and Borgatti 2014 | later |
| Valente-Foreman integration and radiality; beta-reach closeness | A9.3.4-9.3.5 | Valente and Foreman 1998; Agneessens, Borgatti, and Everett 2017 | later |
| Neighbourhood inclusion, rank intervals, expected ranks (netrankr) | S4: the principled answer to "which index" | Schoch and Brandes 2016 | later |
| Bootstrap stability of node measures (node-drop, edge-drop, case) with a correlation-stability coefficient | E (bootnet); the question analysts ask constantly | Epskamp, Borsboom, and Fried 2018; Costenbader and Valente 2003 | later |
| [sna](https://cran.r-project.org/package=sna)'s weak and strong transitivity, directed and weighted clustering with `convention`, cyclicality | R7, NAR3, A10.2.3, U4 | Fagiolo 2007; `onnela2005intensity`; `barrat2004architecture`; Opsahl and Panzarasa 2009 | next |
| Krackhardt's graph-theoretic dimensions (connectedness, hierarchy, efficiency, LUBness) | U4; `everton2012disrupting` ch. 5 | Krackhardt 1994 | next |

## Cohesion and communities

| What | Why | Implement from | Tier |
|---|---|---|---|
| Leiden (modularity and CPM, with the objective documented), walktrap, label propagation, Infomap | R8, S5, NAR2, A11.5, W8 | `traag2019louvain`; Pons and Latapy 2005; Raghavan, Albert, and Kumara 2007; Rosvall and Bergstrom 2008 | now |
| `net_compare_partitions()` (ARI, NMI, VI, cross-table) with canonical labels | R8, S5, NAR2 | `hubert1987assignment`; Danon et al. 2005; Meila 2007 | now |
| Modularity-by-k scan and dendrogram plot for hierarchical methods | R8, R10, A6.4 | | next |
| Clique co-membership and overlap matrices; the bimodal method | A11.2, U5, K: social circles | Everett and Borgatti 1998 (Connections 21(1)) | next |
| Factions with modularity and entailment criteria; cliques with a maximum size and strong/weak directed variants | A11.4, A11.2 | `borgatti2002ucinet` | later |
| Community composition tables and within/between odds ratios | R8 | | next |
| Motif census for 3- and 4-node subgraphs with z-scores against a degree-preserving null; power-law fit with x_min selection | NAR2, NAR3 | Milo et al. 2002; Wernicke 2006; `clauset2009powerlaw` | later |

## Positions and roles

| What | Why | Implement from | Tier |
|---|---|---|---|
| Multi-relation profiles in structural equivalence; blockmodels with one image per relation; mixed image-graph plot | R10, A12: four standardised matrices bound by hand | `burt1976positions`; `wasserman1994social` ch. 9-10 | next |
| Role (ego-centred triad) census with 36 types | R10 | Hummell and Sodeur 1987 | later |
| Equivalence measures: matches, positive matches, overlaps, cross-products, coverage; diagonal rule; geodesic profiles | A12.2-3 | `borgatti2002ucinet` | later |
| Direct blockmodel search (structural and regular, generalised block types) | A12.5-6, S5, U8 | `doreian2005generalized` | later |
| Non-metric MDS with stress-1 and a Shepard diagram | A6.2 | `kruskal1964multidimensional`; `borg2005modern` | later |

## Two-mode

| What | Why | Implement from | Tier |
|---|---|---|---|
| Projection weights: covariance, correlation, Jaccard, Bonacich, cosine | A13.2-3, U3 | `borgatti1997network`; `everett2013dual` | now |
| Backbone extraction (SDSM with the BiCM null, FDSM by curveball sampling, fixed-row, fixed-column, fixed-fill) returning p-values as an edge column; disparity filter for one-mode weighted networks | S6, A13, U: dark-network data are mostly two-mode | Neal 2014; Saracco et al. 2015; Strona et al. 2014; Serrano, Boguna, and Vespignani 2009 | now |
| Two-mode clustering coefficients (six-cycle global and local) | S6, R11: [igraph](https://r.igraph.org/)'s transitivity returns 0 on bipartite graphs silently | Opsahl 2013 | next |
| Two-mode core-periphery, dual-projection Louvain, bi-cliques, two-mode equivalence | A13.5-13.7 | `everett2013dual`; `borgatti1997network` | later |
| Incidence-matrix constructor; two attribute tables (one per mode) in `twomode()` | W1, R11: the stacked attribute frame with NA blocks | | next |
| Share of projected ties matching on an attribute; per-group homogeneity | R11 | | later |

## Inference

| What | Why | Implement from | Tier |
|---|---|---|---|
| `cug_test()` with vector-valued statistics (per-element two-tailed p) and degree-sequence conditioning; a triad-census CUG table | R7, R9, R13, S16 | `holland1976local` (expected counts under U\|MAN); `wasserman1994social` ch. 14 | now |
| Tau statistic for ranked-clustering, clustering, and balance hypotheses | R9 | Holland and Leinhardt 1970, 1976; Johnsen 1985 | next |
| Node-level permutation tests (t-test, ANOVA, regression) and group-level permutation tests | A14.3-4, U9: "does centrality differ by attribute?" | `hubert1987assignment`; Hubert and Schultz 1976; `krackhardt1988predicting` | now |
| ANOVA density models (constant homophily, variable homophily, structural blockmodel) with permutation p | A14, A10.5, U9 | `hanneman2005introduction` ch. 18 | next |
| Network autocorrelation (lag and error) model by maximum likelihood | A14.4, K9 | Doreian 1980; Leenders 2002; `anselin1988spatial` | next |
| Observed-statistics helpers: edgewise and dyadwise shared partners, nodematch counts; a simulate-and-compare table for any model that simulates networks | R13, S18 | `hunter2008ergm` for the statistics | next |
| LR-QAP (logistic MRQAP) | A14.5.2 | `dekker2007sensitivity` | later |
| Diffusion simulation on a fixed network (simple and threshold contagion, exposure matrix) | R14, K9, U2 | Granovetter 1978; Valente 1995; Centola and Macy 2007 | later |
| ERGM, TERGM, SAOM, REM estimation | R13, R15, S18-20, A15, U10-11 | | out of scope on `main` (branch `ergm-saom`; bridges below) |

## Temporal

| What | Why | Implement from | Tier |
|---|---|---|---|
| Panel waves as spells `[k, k+1)`, with helpers from and to a list of networks (*decided*, 9) | R3, R13, R15, S19, U11: every longitudinal chapter starts from a list of waves | `snijders2010introduction`; `butts2026networkdynamic` for conventions | now |
| Node-by-wave results keyed by id with changing composition; Jaccard tie stability per transition | R3, R15: tSnaStats cannot align node measures; Jaccard computed by hand | `snijders2010introduction` (the rule of thumb) | now |
| Snapshot comparison table (density, centralization, fragmentation, key players per period) and node-set alignment with presence flags | U11, A7.6 | `everton2012disrupting` | next |
| Fixed or interpolated layout across periods | R5, R8, R13, S12, U11 | Gansner, Koren, and North 2005 (stress with fixed starts) | next |
| Half-open spells, censoring, and the node-presence rule stated in `?temporal_network` and in the definitions article (*decided*) | R3 | `holme2012temporal`; `latapy2018stream` | now |
| Dynamic network movie, relational-event and epidemic simulation | R5, R13, R14 | | out of scope (bridges) |

## Spatial

| What | Why | Implement from | Tier |
|---|---|---|---|
| Planarity test (graph) and geometry crossing count (spatial), under distinct names | Definitions article; road networks | `barthelemy2011spatial` | later |
| Great-circle curved arcs in spatial plots; map backgrounds | W10, S12 | | later |
| Opt-in `check_geometry()` stays opt-in (*decided*, article) | sfnetworks enforces at construction | `vandermeer2026sfnetworks` | now (documented) |

## Signed

| What | Why | Implement from | Tier |
|---|---|---|---|
| `sign` edge role (-1/+1) and `as_signed()` (*decided* provisionally, 11) | S7, A9.6, E (correlations), U3 ([Pajek](http://mrvar.fdv.uni-lj.si/pajek/) and [UCINET](https://sites.google.com/site/ucinetsoftware/home) values) | `harary1953notion`; `cartwright1956structural` | next |
| Structural balance (triangles, walk-based), frustration index, signed degree, signed blockmodels by Doreian-Mrvar relocation, PN index | S7, A9.6, K2 | `cartwright1956structural`; `doreian1996partitioning`; Aref and Wilson 2018; Estrada 2019; Everett and Borgatti 2014 | later |
| Signed backbones from two-mode projections | S6 | Neal 2014 | later |

## Missing data

| What | Why | Implement from | Tier |
|---|---|---|---|
| `na` node role (non-respondents) and `na` edge role (unobserved dyads); measures refuse unresolved missingness unless `na_action` is given (*decided* provisionally, 11) | R4, A4-5, U3: an entire chapter of matrix surgery that the constructor currently refuses to start | `kossinets2006effects`; `robins2004missing`; `butts2008network` (`na` attribute) | next |
| Imputation: asymmetric, symmetric, probabilistic (repeated draws summarised), reconstruction; bias-calculator inputs (size, percent missing, in-degree SD) | R4, A5.5.5 | Smith and Moody 2013; Smith, Moody, and Morgan 2017; Smith, Morgan, and Moody 2022; Huisman 2009; Znidarsic, Ferligoj, and Doreian 2012 | later |
| Structural zeros and ones as a rule or pair table; [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/) codes on export | R13, R15, A15 | `ripley2026manual` | later |

## Psychometric and correlation networks

| What | Why | Implement from | Tier |
|---|---|---|---|
| `net_cor_network()`: correlation and partial-correlation networks with significance thresholds (none, Bonferroni, Holm, FDR), weights as similarities with sign kept apart | E (the whole first half), NAR1 | Lauritzen 1996 (GGM definition); Epskamp and Fried 2018 (described only) | next |
| Standardised multi-measure, multi-network centrality comparison plot | E (`centralityPlot`); multiplex comparisons | `epskamp2012qgraph` for the figure | next |
| Two-network comparison by permutation (global strength, edge differences) | E (NCT) | van Borkulo et al. 2023 | later |
| EBICglasso, Ising models, time-series GGMs, polychoric correlations | E | Friedman, Hastie, and Tibshirani 2008; Chen and Chen 2008; Foygel and Drton 2010; Olsson 1979 | out of scope (opt-in area at most) |

## Plotting

| What | Why | Implement from | Tier |
|---|---|---|---|
| Node and edge coordinate tables from a layout (`as.data.frame()` of a layout joined to the node table, edges with x, y, xend, yend) so ggplot2 users draw in three `geom_*` calls | S10-14, NAR4, R5: the sources are entirely ggraph | | now |
| Arrows with end caps by default for directed networks; `groups =` hulls; halo labels; `label_top = k`; `edge_scale = c(min, cut, max)`; a layout argument that accepts saved coordinates or another network's layout | R5, S10, S13, W7, E | Epskamp et al. 2012 for the scale arguments | next |
| Layouts: Kamada-Kawai, MDS, Sugiyama, pivot MDS, backbone, focus and centrality (concentric), multilevel, dendrogram | S12, NAR4, W2 | `kamada1989algorithm`; `sugiyama1981methods`; Brandes and Pich 2007; Nocaj, Ortmann, and Brandes 2015 | next |
| Shared (averaged) layout and shared edge scale across several networks on one node set | E (`averageLayout`), W2 | Gansner, Koren, and North 2005 | next |
| Matrix plot orderings by core, degree, membership | S5 | | next |
| Interactive HTML widgets, edge bundling, arc and chord and hive plots, dynamic movies | W4, W6, W10, NAR4, R5 | | out of scope (JSON node-link export serves sigma, vis, and d3) |

## Interoperability and bridges

| What | Why | Implement from | Tier |
|---|---|---|---|
| `to_network()` carrying node, edge, and dyad covariates, documented as the [ergm](https://cran.r-project.org/package=ergm) and tergm hand-off | R13, S18, A15, U10 | `butts2008network`; `hunter2008ergm` | now |
| `to_siena()`: n x n x T arrays, constant and changing covariates, composition change, structural zeros | R15, S19, A15, U11 | `ripley2026manual` | next |
| `to_relevent()`: (time, sender, receiver) plus covariate arrays | R13, S20 | Butts 2008 (relational event framework) | next |
| `to_diffnet()` and from it | R14 | Valente 1995; netdiffuseR (MIT) as a reference | later |
| egor: three-table import and export (see Ego networks) | R6, S8 | | now |
| networkDynamic and tsna: spells in both directions, including TEAs | R3 | `butts2026networkdynamic` | next |
| Lossless `to_tbl_graph()` and back; documented dplyr idioms on `g$nodes` and `g$edges` | S21-23 | `pedersen2017tidygraph` | next |
| Extra file formats: ncol, lgl, dimacs, graphdb, adjacency-list files; two-mode and multi-matrix DL; [UCINET](https://sites.google.com/site/ucinetsoftware/home) native .##h | S1, R3, A5.4 | | later |
| Result provenance: an attribute on returned vectors recording convention and weighting, printed by `summary()` | A (xucinet's `$assumptions`) | | later |

## Teaching data sets, with licences

Every shipped data set needs a licence that allows redistribution under
MIT and a `data-raw/` script that rebuilds it from its published source.

| Data set | Taught in | Source and licence | Tier |
|---|---|---|---|
| Krackhardt high-tech managers (advice, friendship, reports-to) | R7, A, K | `krackhardt1987qap`; [UCINET](https://sites.google.com/site/ucinetsoftware/home) public Datafiles | now |
| Campnet | A (throughout) | UCINET public Datafiles (Borgatti, Everett, and Freeman); check terms | now |
| Padgett Florentine business ties (marriage is shipped) | S21-23, A | Padgett and Ansell 1993 (`padgett1993robust`); public | now |
| Krackhardt kite | NAR, W | `krackhardt1990assessing`; a published figure, no data licence issue | now |
| Lazega law firm (advice, friendship, co-work; 71 lawyers) | S9, S16, S18, A | Lazega 2001; distributed with permission in several packages; check terms | next |
| Sampson monastery, Davis southern women (shipped), Zachary karate (shipped) | A, S4, NAR | `davis1941deep`; `zachary1977information` | next |
| Les Misérables (Knuth), Game of Thrones season 1, Star Wars | S10, W4, W7 | Stanford GraphBase (`knuth1993stanford`), public; Beveridge and Shan 2016 and Gabasova 2016, redistributable with attribution | next |
| Glasgow (Teenage Friends and Lifestyle) and s50 panels | R15, S19 | [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/) data page; check the stated terms before shipping | later |
| igraphdata sets (UKfaculty, enron, yeast, macaque, immuno, USairports, Koenigsberg) | W, NAR | igraphdata is CC BY-SA 4.0: share-alike is incompatible with bundling in an MIT package; point to the package instead | out of scope |
| Classroom and GSS egocentric data of Rawlings et al. | R3-R15 | companion repository under CC BY-NC-ND: not redistributable; document how to load them | out of scope |
| UCNets, cowList, Senate cosponsorship | S7, S8 | restricted or large; document only | out of scope |

## Cross-cutting lessons

Not features, but rules the surveys converged on; most are already policy
and are listed so they are not lost.

- Node tables stay in their original order; results are vectors aligned
  with `nodes` or `edges`, census-like results are named vectors or data
  frames, never lists that need `$vector`.
- Every argument that exists to match a tool is a `convention` or a named
  option documented in `?conventions`, with one example per convention;
  the surveys add structural holes (ego versus whole model), directed
  clustering definitions, E-I expected values, reciprocity (dyad versus
  arc), QAP p-value counting, and sample versus population SD to the
  table.
- Reproducibility by argument: `seed` and `threads` on every stochastic
  or parallel function, with the seed recorded in the result.
- Print methods explain the number next to it (the Jaccard rule of
  thumb, modularity's range, "p < 1/reps, never 0").
- Help pages carry one sentence of interpretation and the pitfalls the
  tutorials hit (zero weights, ties in `which.max`, the Leiden objective,
  `decompose()[[1]]` not being the largest component).
- Validation against the books' own worked numbers on the shipped data
  sets, with `convention = "ucinet"`, documenting the cases where the
  book's figure is not reproducible and why.
