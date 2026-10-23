# Bundled example data

Every directory holds the CSV tables of one data set, built by a script in
`data-raw/` from the published source named in the directory's `README.md`
(source, licence, build commit and checksums, and the meaning of every
column). Nothing is fetched at install time. Sizes are the uncompressed
CSV bytes; everything below together is about 90 KB gzipped.

| name | configuration | nodes | edges | size | licence | citation |
|---|---|---|---|---|---|---|
| `karate` | one-mode, undirected, valued (`contexts`) | 34 | 78 | 0.9 KB | published data (Zachary 1977) | Zachary (1977), doi:10.1086/jar.33.4.3629752 |
| `florentine` | one-mode, undirected | 15 | 20 | 0.4 KB | published data | Padgett & Ansell (1993), doi:10.1086/230190 |
| `southern_women` | two-mode | 18 + 14 | 89 | 1.6 KB | published data | Davis, Gardner & Gardner (1941) |
| `jemaah_islamiyah` | one-mode, undirected, valued (`strength`) | 17 | 63 | 1.4 KB | published data (Koschade 2006, Tables 1-2) | Koschade (2006), doi:10.1080/10576100600798418 |
| `noordin` | one-mode, undirected, multiplex (10 relations); two-mode affiliations | 79; 79 + 45 | 1014; 275 | 44 KB | ARDA data via networkdata (MIT) and manet | Roberts & Everton (2011); Everton (2012), doi:10.1017/CBO9781139136877 |
| `krackhardt` | one-mode, directed, multiplex (`relation`: advice, friendship, reports_to), node attributes | 21 | 312 | 4.5 KB | published data, no stated licence (see README) | Krackhardt (1987), doi:10.1016/0378-8733(87)90009-8 |
| `sampson` | one-mode, directed, signed (`sign`), multiplex (8 relations), ranked (`rank`), waves (`time` 1-3 for liking) | 18 | 510 | 16 KB | published data, no stated licence (see README) | Sampson (1968); Breiger, Boorman & Arabie (1975), doi:10.1016/0022-2496(75)90028-0 |
| `highland_tribes` | one-mode, undirected, signed (`sign`) | 16 | 58 | 0.9 KB | published data, no stated licence (see README) | Read (1954), doi:10.1086/soutjanth.10.1.3629074; Hage & Harary (1983) |
| `les_miserables` | one-mode, undirected, weighted strength (`encounters`), isolates | 80 | 254 | 10 KB | Stanford GraphBase: freely copyable, derived tables renamed | Knuth (1993) |
| `miles` | one-mode, undirected, weighted distance (`miles`), spatial points (`latitude`, `longitude`), complete graph | 128 | 8128 | 297 KB | Stanford GraphBase: freely copyable, derived tables renamed | Knuth (1993) |
| `country_borders` | one-mode, undirected, spatial (node points; edge `geometry` as WKT, EPSG:4326), weighted strength (`border_km`), isolates | 177 | 313 | 87 KB | public domain (Natural Earth) | Natural Earth (2022), 1:110m Admin 0 Countries |
| `workplace_contacts` | one-mode, undirected, temporal contacts (`time`), node attribute (`department`) | up to 100 (tags worn; Génois et al. 2015, Table 1) | built on CI (see below) | built on CI | CC0 (SocioPatterns) | Génois et al. (2015), doi:10.1017/nws.2015.10 |

`workplace_contacts` is downloaded from sociopatterns.org by
`data-raw/sociopatterns.py`; the download is not possible from every
build environment, so its CSV files are produced by the "Rebuild example
data" workflow (`.github/workflows/data-raw.yml`) and committed from its
output branch.

## Coverage of the network configurations

| configuration | data sets |
|---|---|
| one-mode undirected | `karate`, `florentine`, `les_miserables`, `highland_tribes`, `miles`, `country_borders` |
| one-mode directed | `krackhardt`, `sampson` |
| two-mode | `southern_women`, `noordin` affiliations |
| weighted, strength | `karate` (`contexts`), `jemaah_islamiyah` (`strength`), `les_miserables` (`encounters`), `country_borders` (`border_km`) |
| weighted, distance | `miles` (`miles`) |
| signed | `highland_tribes`, `sampson` |
| multiplex | `noordin`, `krackhardt`, `sampson` |
| temporal: contacts | `workplace_contacts` |
| temporal: waves | `sampson` (liking at three times) |
| spatial | `miles` (points), `country_borders` (points and edge geometry) |
| disconnected, isolates | `les_miserables`, `country_borders` |
| ranked choices | `sampson` |

Not yet covered by a shipped data set (candidates examined, with the
reason each was not shipped):

- **k-mode** (people, organizations, events): the `manet` copy of the
  Noordin affiliations does not say which of its 45 columns are
  organizations, operations, trainings and meetings; the complete Noordin
  data at the ARDA carry the ARDA's own terms of use, not checked.
- **Multilayer with interlayer edges beyond couplings**: no small data set
  with a stated open licence was found; the AUCS (Aarhus CS department)
  multiplex of Rossi and Magnani (2015) has no licence statement on its
  page. `krackhardt`, `sampson` and `noordin` are multiplex (the coupling
  special case).
- **Panel waves with changing composition, time-varying attributes, node
  spells and missing data**: Knecht's (2008) classroom (RSiena's
  `klas12b`: four waves, delinquency and alcohol per wave, a leaver, 9 =
  missing codes) is the right data set and `data-raw/knecht.py` builds
  it, but the RSiena page states no terms and the btergm package needed
  Andrea Knecht's written permission (17 April 2014) to redistribute it;
  not shipped until the owner has the same permission. Lazega's lawyers
  (RSiena page, no stated terms) and the Teenage Friends and Lifestyle
  excerpt `s50` ("for illustration") were not shipped for the same
  reason.
- **Temporal edge spells**: can be derived from `workplace_contacts`
  (consecutive 20-second contacts); no separate data set.
- **Spatio-temporal**: no small data set with a clear licence found
  (SocioPatterns sets other than the workplace are CC BY-NC-SA; flight
  data from KONECT lack provenance).
- **Ego networks from a survey**: the 2004 General Social Survey
  extract in networkdata (`gss_egor`) is the natural candidate, but NORC's
  terms say its contents may not be reproduced or transmitted without
  written consent; egos are better extracted from the bundled whole
  networks.
