# Face-to-face contacts in a workplace (SocioPatterns)

One-mode, undirected, temporal (contacts): the face-to-face proximity
contacts recorded by wearable sensors among employees of an office
building of the French health observatory InVS over two weeks in June
2013 (Génois et al. 2015; 100 tags were worn, about two thirds of the
staff, in five departments on two floors), with each person's
department. A contact is a 20-second interval in which two sensors
detected each other at close
range; repeated contacts between the same pair are separate rows, so the
network is a contact sequence from which spells, aggregated weights or
time windows can be built.

## Files

`nodes.csv`, one per participant:

| column | meaning |
|---|---|
| `name` | the anonymous identifier from the source |
| `department` | the person's department (the source's `Di` code) |

`contacts.csv`, one per 20-second contact:

| column | meaning |
|---|---|
| `from`, `to` | the two people (`from` < `to`) |
| `time` | the end of the 20-second interval, in seconds since 0:00 on 24 June 2013 ("the contact was active during the interval [t - 20s, t]") |

Rows are sorted by time. Several rows can share a time (simultaneous
contacts), and a pair in contact for a minute gives three rows 20 seconds
apart.

## Source and build

Built by `data-raw/sociopatterns.py`, which downloads
`workplace_InVS_tij.dat.zip` and `workplace_InVS_metadata.txt` from
http://www.sociopatterns.org/datasets/contacts-in-a-workplace/ and
prints their SHA-256 (recorded in the script once known). The data set's
page describes the formats quoted above.

## Licence

CC0. The data set's page: "The data are distributed to the public under a
Creative Commons Public Domain Dedication license." The page asks that
publications using the data cite Génois et al. (2015).

## Citation

Génois, M., Vestergaard, C. L., Fournet, J., Panisson, A., Bonmarin, I.,
& Barrat, A. (2015). Data on face-to-face contacts in an office building
suggest a low-cost vaccination strategy based on community linkers.
*Network Science*, 3(3), 326-347. https://doi.org/10.1017/nws.2015.10
