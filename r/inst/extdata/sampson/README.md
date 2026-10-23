# Sampson's monastery

One-mode, directed, signed, multiplex, with three waves of one relation:
18 novices in a New England monastery and the ranked sociometric choices
Sampson (1968) collected while living there as an experimenter. Four
relations were asked in positive and negative form, and each novice ranked
his top three choices on each (a few gave tied ranks for four). Most of
the data are retrospective, collected after the "crisis in the cloister"
in which four novices were expelled and most of the others then left;
liking was also gathered at three times. This is the UCINET standard data
set SAMPSON (ten 18 x 18 matrices).

## Files

`nodes.csv` (18 rows), one per novice:

| column | meaning |
|---|---|
| `name` | the novice's name (Sampson's pseudonyms, in full) |
| `number` | Sampson's number for him, 1-18 |

`edges.csv` (510 rows), one per ranked choice:

| column | meaning |
|---|---|
| `from`, `to` | the novice making and receiving the choice |
| `relation` | `liking`, `disliking`, `esteem`, `disesteem`, `positive_influence`, `negative_influence`, `praise`, `blame` |
| `sign` | +1 for the positive relations, -1 for `disliking`, `disesteem`, `negative_influence`, `blame` |
| `time` | 1, 2 or 3 for `liking` (the three times it was gathered); 3 for the other relations, collected in the same wave as the third liking measure |
| `rank` | 3 = first choice, 2 = second, 1 = third, as coded in the UCINET file |

The UCINET matrix behind each `relation` and `time`: SAMPLK1-3 (liking at
times 1-3), SAMPDLK, SAMPES, SAMPDES, SAMPIN, SAMPNIN, SAMPPR, SAMPNPR.
No self-ties; no parallel rows within a relation and time.

Sampson's groups (the Young Turks, the Loyal Opposition, the Outcasts and
an interstitial group) and the Cloisterville attribute are not in the
UCINET file and are not shipped; see White, Boorman and Breiger (1976)
and Breiger, Boorman and Arabie (1975), who analyse these data.

## Source and build

Built by `data-raw/ucinet.py` from the `sampson` object of the
MIT-licensed R package
[networkdata](https://github.com/schochastics/networkdata) (commit
`521594049f725c3962c11496d5b0f48386ee84b2`, 21 April 2026; SHA-256
recorded in the script), a list of the ten UCINET matrices copied from
Linton Freeman's UCINET data page. The UCINET dataset page
(https://sites.google.com/site/ucinetsoftware/datasets/sampson-monastery)
describes them. The abbreviated vertex names of that copy ("Romul_10") are
expanded to the full names used by the UCINET and Pajek descriptions.

## Licence

Published research data with no licence from their author: the data are
from Sampson's (1968) dissertation, were analysed in Breiger,
Boorman and Arabie (1975) and White, Boorman and Breiger (1976), and are
redistributed without stated terms by UCINET, Pajek and networkdata
(MIT). Redistributed here with attribution on that basis.

## Citation

Sampson, S. F. (1968). *A Novitiate in a Period of Change: An Experimental
and Case Study of Social Relationships*. PhD dissertation, Cornell
University.

Breiger, R. L., Boorman, S. A., & Arabie, P. (1975). An algorithm for
clustering relational data with applications to social network analysis
and comparison with multidimensional scaling. *Journal of Mathematical
Psychology*, 12(3), 328-383. https://doi.org/10.1016/0022-2496(75)90028-0

White, H. C., Boorman, S. A., & Breiger, R. L. (1976). Social structure
from multiple networks. I. Blockmodels of roles and positions. *American
Journal of Sociology*, 81(4), 730-780. https://doi.org/10.1086/226141
