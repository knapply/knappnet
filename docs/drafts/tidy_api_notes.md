# Notes: a tidy network API (not scheduled)

Motivation: problems users hit with tidygraph — activate() mode state;
centrality_*()/.G() only work in context; integer from/to re-index after
filtering; edge rules cannot reach node attributes without joins; igraph
round trips per verb; morph/unmorph complexity; no time, layers, egos,
missing data, or validation.

Principles
1. Explicit table verbs (mutate_nodes, filter_edges, summarise_nodes, ...);
   plain dplyr verbs on g$nodes / g$edges.
2. Measures are ordinary functions that know their graph inside a verb.
3. Edge expressions reach endpoint attributes: from(dept) == to(dept).
4. Node expressions aggregate over neighbours: neighbors(mean(age)),
   in_neighbors(), out_neighbors(), with .edge for incident edge columns.
5. Nodes identified by label; filtering nodes drops incident edges, no
   renumbering.
6. group_by() + summarise_*() collapses to a quotient network (blockmodel).
7. Time/layers/egos first class: at(), during(), slice_time(), across_time(),
   layer(), egos().
8. Immutable, cached engine graph across steps; validated; card printed.
9. measures(g, betweenness, closeness) -> one tibble in one engine pass;
   broom-style tidy() for models.
10. Works on igraph / tbl_graph inputs, returns the same class.

Open: naming (unprefixed verbs vs net_ prefix vs companion package).
