# Two-mode networks: construction, conversion, projections, and centrality.
# Hand-checked values first, then spot checks against igraph. The full
# comparison lives in validation/.

# People ann, bob, cat, dan, eve (first mode) and events gala, fair, show.
# Incidence matrix (gala, fair, show):
#   ann 1 1 0 / bob 1 0 1 / cat 0 0 1 / dan 0 1 0 / eve 1 0 1
party <- function() {
  net_twomode(data.frame(
    person = c("ann", "ann", "bob", "bob", "cat", "dan", "eve", "eve"),
    event  = c("gala", "fair", "gala", "show", "show", "fair", "gala", "show"),
    hours  = c(1, 2, 3, 4, 5, 6, 7, 8)
  ), source = "person", target = "event", weight = "hours", weight_type = "strength")
}

# 0 for nodes of the first mode, 1 for the second.
types <- function(g) knappnet:::twomode_types(g)

test_that("net_twomode() builds the node table from the two edge columns", {
  g <- party()
  expect_true(net_is_twomode(g))
  expect_false(g$directed)
  expect_equal(g$nodes$id, c("ann", "bob", "cat", "dan", "eve", "gala", "fair", "show"))
  expect_equal(g$nodes$mode, rep(c("person", "event"), c(5, 3)))
  expect_equal(g$roles$mode, "mode")
  expect_equal(g$roles$mode_levels, c("person", "event"))
  expect_true(g$declared$between_modes_only)
  expect_equal(types(g), rep(0:1, c(5, 3)))
  expect_false(net_is_twomode(net_frame(data.frame(source = 1, target = 2))))
  expect_error(net_twomode(data.frame(source = c("x", "y"), target = c("y", "z"))), "both")
  # Endpoint columns are never guessed.
  expect_error(net_twomode(data.frame(a = "p", b = "x")), "source")
})

test_that("modes from a node column", {
  nodes <- data.frame(id = c("p", "q", "x", "y"), kind = c("person", "person", "club", "club"))
  edges <- data.frame(source = c("p", "q", "q"), target = c("x", "x", "y"))
  g <- net_twomode(edges, nodes, mode = "kind")
  # The first mode is the first to appear in the node table, unless `first`
  # says otherwise.
  expect_equal(g$roles$mode_levels, c("person", "club"))
  expect_equal(types(g), c(0L, 0L, 1L, 1L))
  expect_equal(g$nodes$kind, nodes$kind)
  h <- net_twomode(edges, nodes, mode = "kind", first = "club")
  expect_equal(h$roles$mode_levels, c("club", "person"))
  expect_equal(types(h), c(1L, 1L, 0L, 0L))
  # igraph's logical `type` column is accepted as the mode column.
  lg <- net_twomode(edges, data.frame(id = nodes$id, type = c(FALSE, FALSE, TRUE, TRUE)))
  expect_equal(lg$roles$mode, "type")
  expect_equal(types(lg), c(0L, 0L, 1L, 1L))
  dyn <- data.frame(id = nodes$id, class = c("agent", "agent", "event", "event"))
  d <- net_as_twomode(net_frame(edges, dyn), mode = "class")
  expect_true(net_is_twomode(d))
  expect_equal(d$roles$mode_levels, c("agent", "event"))
  expect_error(net_twomode(edges, nodes, mode = "kind", first = "robot"), "one of")
  expect_error(net_twomode(edges, nodes, mode = "nope"), "not found")
  three <- transform(nodes, kind = c("a", "b", "c", "c"))
  expect_error(net_twomode(edges, three, mode = "kind"), "exactly two")
  within <- data.frame(source = c("p", "p"), target = c("x", "q"))
  expect_error(net_twomode(within, nodes, mode = "kind"), "same mode")
})

test_that("net_as_twomode() converts matrices, data frames, and net_frames", {
  m <- matrix(c(1, 0, 1,
                1, 1, 0), 2, byrow = TRUE, dimnames = list(c("p", "q"), c("x", "y", "z")))
  g <- net_as_twomode(m)
  expect_equal(g$nodes$id, c("p", "q", "x", "y", "z"))
  expect_equal(types(g), c(0L, 0L, 1L, 1L, 1L))
  expect_equal(nrow(g$edges), 4L)
  expect_null(g$roles$weight)
  expect_equal(unname(net_incidence(g)), unname(m))
  m[1, 1] <- 3
  w <- net_as_twomode(m)
  expect_equal(w$roles$weight, "weight")
  expect_equal(net_incidence(w, attr = "weight"), m)
  expect_equal(net_as_twomode(unname(m))$nodes$id, c("r1", "r2", "c1", "c2", "c3"))
  e <- net_as_twomode(data.frame(source = c("p", "q"), target = c("x", "x")))
  expect_equal(types(e), c(0L, 0L, 1L))
  expect_identical(net_as_twomode(e), e)
  expect_error(net_as_twomode(net_frame(data.frame(source = 1, target = 2))), "mode")
})

test_that("a DyNetML file with two node classes becomes two-mode in one call", {
  nodes <- data.frame(id = c("ann", "bob", "gala", "fair"),
                      class = c("agent", "agent", "event", "event"),
                      type = c("Agent", "Agent", "Event", "Event"))
  g <- net_frame(data.frame(source = c("ann", "bob", "ann"), target = c("gala", "gala", "fair")),
                 nodes, directed = FALSE)
  path <- tempfile(fileext = ".xml")
  net_write(g, path, format = "dynetml")
  h <- net_as_twomode(net_read(path), mode = "class")
  expect_true(net_is_twomode(h))
  expect_equal(types(h), c(0L, 0L, 1L, 1L))
  expect_equal(h$roles$mode_levels, c("agent", "event"))
  expect_equal(net_twomode_density(h), 3 / 4)
})

test_that("incidence matrices", {
  g <- party()
  b <- net_incidence(g)
  expect_equal(dim(b), c(5L, 3L))
  expect_equal(rownames(b), c("ann", "bob", "cat", "dan", "eve"))
  expect_equal(colnames(b), c("gala", "fair", "show"))
  expect_equal(unname(b["bob", ]), c(1, 0, 1))
  expect_equal(net_incidence(g, attr = "hours")["eve", "show"], 8)
  # A directed network may store ties either way.
  d <- net_twomode(data.frame(source = c("p", "x"), target = c("x", "q")),
               data.frame(id = c("p", "q", "x"), type = c(FALSE, FALSE, TRUE)),
               directed = TRUE)
  expect_equal(unname(net_incidence(d)[, 1]), c(1, 1))
})

test_that("projections count shared partners", {
  g <- party()
  p <- net_project(g)
  expect_s3_class(p, "net_frame")
  expect_false(p$directed)
  expect_equal(p$nodes$id, c("ann", "bob", "cat", "dan", "eve"))
  expect_equal(p$edges$source, c("ann", "ann", "ann", "bob", "bob", "cat"))
  expect_equal(p$edges$target, c("bob", "dan", "eve", "cat", "eve", "eve"))
  expect_equal(p$edges$weight, c(1, 1, 1, 1, 2, 1))
  expect_equal(p$roles$weight_type, "strength")
  # gala and show each have three participants, fair two.
  n <- net_project(g, method = "newman")
  expect_equal(n$edges$weight, c(1 / 2, 1, 1 / 2, 1 / 2, 1, 1 / 2))
  s <- net_project(g, "second", loops = TRUE)
  expect_equal(s$edges$source, c("gala", "gala", "gala", "fair", "show"))
  expect_equal(s$edges$target, c("gala", "fair", "show", "fair", "show"))
  expect_equal(s$edges$weight, c(3, 1, 2, 2, 3))
  # Weighted by hours: gala-show through bob (3 * 4) and eve (7 * 8).
  w <- net_project(g, "second", method = "product")
  expect_equal(w$edges$weight, c(1 * 2, 3 * 4 + 7 * 8))
  both <- net_project(g, "both", weighted = FALSE)
  expect_equal(names(both), c("first", "second"))
  expect_null(both$first$roles$weight)
  expect_equal(nrow(both$second$edges), 2L)
  expect_error(net_project(g, method = "newman", loops = TRUE), "diagonal")
  expect_error(net_project(net_frame(data.frame(source = 1, target = 2))), "two-mode")
  plain <- net_twomode(data.frame(source = "p", target = "x"))
  expect_error(net_project(plain, method = "product"), "weights")
  expect_equal(nrow(net_project(plain)$edges), 0L)
  b <- net_incidence(g)
  full <- b %*% t(b)
  pf <- net_project(g, loops = TRUE)
  ids <- p$nodes$id
  expect_equal(pf$edges$weight, full[cbind(match(pf$edges$source, ids),
                                           match(pf$edges$target, ids))])
})

test_that("two-mode centrality uses Borgatti and Everett's normalizations", {
  g <- party()
  expect_equal(net_twomode_density(g), 8 / 15)
  expect_equal(net_twomode_degree(g), c(c(2, 2, 1, 1, 2) / 3, c(3, 2, 3) / 5))
  expect_equal(net_twomode_degree(g, normalized = FALSE), c(2, 2, 1, 1, 2, 3, 2, 3))
  raw <- c(10, 4, 0, 0, 4, 12.5, 6, 6.5)
  expect_equal(net_twomode_betweenness(g, normalized = FALSE), raw)
  # Maxima: 16 for people (5 people, 3 events) and 19 for events.
  expect_equal(net_twomode_betweenness(g), raw / rep(c(16, 19), c(5, 3)))
  # cat: distances 1, 2, 2, 3, 4, 5, 6 sum to 23; (3 + 2 * 4) / 23.
  expect_equal(net_twomode_closeness(g)[3], 11 / 23)
  expect_equal(net_twomode_closeness(g, unreachable = "ignore")[3], 11 / 23)
  expect_equal(net_twomode_closeness(g, normalized = FALSE, unreachable = "ignore")[3], 1 / 23)
  expect_error(net_twomode_closeness(g, normalized = FALSE), "Wasserman-Faust")
  ev <- net_twomode_eigen(g)
  expect_equal(sum(ev[1:5]^2), 1)
  expect_equal(sum(ev[6:8]^2), 1)
  b <- net_incidence(g)
  sv <- svd(b)
  expect_equal(ev[1:5], abs(sv$u[, 1]), tolerance = 1e-8)
  expect_equal(ev[6:8], abs(sv$v[, 1]), tolerance = 1e-8)
})

test_that("disconnected two-mode networks", {
  g <- net_twomode(data.frame(source = c("p", "q"), target = c("x", "y")))
  # Each node reaches one other at distance 1: (2 + 2 * 1) / 1, times 1/3
  # under the Wasserman-Faust rule.
  expect_equal(net_twomode_closeness(g), rep(4 / 3, 4))
  expect_equal(net_twomode_closeness(g, unreachable = "ignore"), rep(4, 4))
  expect_warning(net_twomode_eigen(g), "not connected")
  lone <- net_twomode(data.frame(source = "p", target = "x"),
                      data.frame(id = c("p", "x", "z"), type = c(FALSE, TRUE, FALSE)))
  expect_equal(net_twomode_closeness(lone)[3], 0)
  expect_true(is.nan(net_twomode_closeness(lone, unreachable = "ignore")[3]))
})

test_that("results match igraph", {
  skip_if_not_installed("igraph")
  g <- party()
  ig <- net_to_igraph(g)
  expect_true(igraph::is_bipartite(ig))
  plain <- igraph::delete_edge_attr(igraph::delete_edge_attr(ig, "weight"), "hours")
  proj <- igraph::bipartite_projection(plain, multiplicity = TRUE)
  p <- net_project(g)
  ends <- igraph::as_edgelist(proj$proj1)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  expect_setequal(key(p$edges$source, p$edges$target), key(ends[, 1], ends[, 2]))
  w <- igraph::E(proj$proj1)$weight
  expect_equal(p$edges$weight[match(key(ends[, 1], ends[, 2]), key(p$edges$source, p$edges$target))],
               w)
  expect_equal(net_twomode_betweenness(g, normalized = FALSE),
               unname(igraph::betweenness(ig, weights = NA)))
  expect_equal(net_twomode_closeness(g, normalized = FALSE, unreachable = "ignore"),
               unname(igraph::closeness(ig, weights = NA)))
  expect_equal(igraph::V(ig)$type, types(g) == 1L)
  back <- net_as_twomode(ig)
  expect_true(net_is_twomode(back))
  expect_equal(back$nodes$type, types(g) == 1L)
  bip <- igraph::make_bipartite_graph(c(FALSE, FALSE, TRUE), c(1, 3, 2, 3))
  expect_equal(types(net_as_twomode(bip)), c(0L, 0L, 1L))
})

test_that("bipartite networks from the network package", {
  skip_if_not_installed("network")
  net <- network::network.initialize(4, directed = FALSE, bipartite = 2)
  network::add.edges(net, tail = c(1, 1, 2), head = c(3, 4, 3))
  network::set.vertex.attribute(net, "vertex.names", c("p", "q", "x", "y"))
  g <- net_as_twomode(net)
  expect_true(net_is_twomode(g))
  expect_equal(types(g), c(0L, 0L, 1L, 1L))
})
