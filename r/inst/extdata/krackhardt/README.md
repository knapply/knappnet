# Krackhardt's high-tech managers

One-mode, directed, multiplex: 21 managers of a small manufacturing firm
on the west coast of the United States and three relations among them.

Krackhardt (1987) asked every manager "To whom do you go to for advice?"
and "Who is your friend?", and took the formal reporting structure from
company records. The network shipped here is the *locally aggregated
structure* of the cognitive social structure study (the ties each manager
reported about himself), which is the form distributed as the UCINET
standard data set "Krackhardt High-Tech Managers" and printed in
Wasserman and Faust (1994, Appendix B).

## Files

`nodes.csv` (21 rows), one per manager:

| column | meaning |
|---|---|
| `name` | manager number 1-21, as in the source |
| `age` | age in years |
| `tenure` | years with the firm (whole years in this copy) |
| `level` | level in the hierarchy: 1 = CEO (manager 7), 2 = vice president (2, 14, 18, 21), 3 = manager |
| `department` | department 1-4; 0 for the CEO, who heads no department |

`edges.csv` (312 rows), one per directed tie, with the `relation` edge
role:

| `relation` | ties | meaning of a row `from -> to` |
|---|---|---|
| `advice` | 190 | `from` goes to `to` for advice |
| `friendship` | 102 | `from` names `to` as a friend |
| `reports_to` | 20 | `from` reports to `to` (the organization chart; every manager but the CEO has exactly one) |

No self-ties; no parallel ties within a relation.

## Source and build

Built by `data-raw/ucinet.py` from the `ht_advice`, `ht_friends` and
`ht_reports` objects of the MIT-licensed R package
[networkdata](https://github.com/schochastics/networkdata) (commit
`521594049f725c3962c11496d5b0f48386ee84b2`, 21 April 2026; SHA-256 of the
three `.rda` files recorded in the script), which copied them from Linton
Freeman's UCINET data page. The UCINET dataset page
(https://sites.google.com/site/ucinetsoftware/datasets/krackhardt-high-tech-managers)
describes the same three 21 x 21 matrices and the 21 x 4 attribute table.

## Licence

Published research data with no licence from their authors: the data
were published in Krackhardt (1987) and are printed in Wasserman and Faust
(1994, Appendix B, pages 738-755), and are redistributed without stated
terms by UCINET and by networkdata (MIT). (The UCI Network Data Repository
attaches a CC BY-NC-ND 2.5 notice to its own copy, taken from the Pajek
"Wasserman and Faust" collection; that copy is not used here.)
Redistributed here with attribution on that basis.

## Citation

Krackhardt, D. (1987). Cognitive social structures. *Social Networks*,
9(2), 109-134. https://doi.org/10.1016/0378-8733(87)90009-8

Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
Applications*. Cambridge University Press.
https://doi.org/10.1017/CBO9780511815478
