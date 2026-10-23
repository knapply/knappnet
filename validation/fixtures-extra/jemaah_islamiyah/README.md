# Jemaah Islamiyah (Bali bombing cell, 2002)

The network of the 17 members of the Jemaah Islamiyah cell behind the 2002
Bali bombings, from Koschade (2006):

- `nodes.csv`: each member's `name` and `role` in the operation (command
  team, operation assistant, bomb maker, suicide bomber, or Team Lima).
- `edges.csv`: the 63 undirected ties of the article's Table 1 (binary
  relations: the two members were known to have interacted).
- `edges_weighted.csv`: the same ties with the article's Table 2 interaction
  strength (1 to 5) in `strength`. Table 2 leaves the SARIJO-DULMATIN cell
  blank; the matrix is otherwise symmetric and DULMATIN-SARIJO is 5, so 5
  is used.

The matrices and roles are those of the R package homophily (Knapp 2019,
`data-raw/jemaah_islamiyah.R`), which ships the binary network with roles
as `jemmah_islamiyah`. AZAHARI appears as HUSIN in the article's figure.
The files are kept here so that they can be packaged as an example dataset
later; the tests build the network from them inline.

Reference: Koschade, S. (2006). A social network analysis of Jemaah
Islamiyah: The applications to counterterrorism and intelligence. *Studies
in Conflict & Terrorism*, 29(6), 559-575.
https://doi.org/10.1080/10576100600798418
