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
  ), weight = "hours", weight_type = "strength")
}

test_that("net_twomode() builds the node table from the two edge columns", {
  g <- party()
  expect_true(net_is_twomode(g))
  expect_false(g$directed)
  expect_equal(g$nodes$name, c("ann", "bob", "cat", "dan", "eve", "gala", "fair", "show"))
  expect_equal(g$nodes$type, rep(c(FALSE, TRUE), c(5, 3)))
  expect_equal(g$roles$type, "type")
  expect_equal(g$roles$type_labels, c("FALSE", "TRUE"))
  expect_false(net_is_twomode(net_frame(data.frame(source = 1, target = 2))))
  expect_error(net_twomode(data.frame(a = c("x", "y"), b = c("y", "z"))), "both")
})

test_that("node types from two-valued columns", {
  nodes <- data.frame(id = c("p", "q", "x", "y"), kind = c("person", "person", "club", "club"))
  edges <- data.frame(source = c("p", "q", "q"), target = c("x", "x", "y"))
  g <- net_twomode(edges, nodes, type = "kind")
  # "club" sorts first, so clubs are the first mode.
  expect_equal(g$nodes$type, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(g$roles$type_labels, c("club", "person"))
  expect_equal(g$nodes$kind, nodes$kind)
  h <- net_twomode(edges, nodes, type = "kind", first = "person")
  expect_equal(h$nodes$type, c(FALSE, FALSE, TRUE, TRUE))
  f <- net_twomode(edges, transform(nodes, kind = factor(kind, levels = c("person", "club"))),
               type = "kind")
  expect_equal(f$nodes$type, c(FALSE, FALSE, TRUE, TRUE))
  # A `type` column that is not logical is kept; the new column is named
  # after the second mode.
  dyn <- data.frame(id = nodes$id, class = c("agent", "agent", "event", "event"),
                    type = c("Agent", "Agent", "Event", "Event"))
  d <- net_as_twomode(net_frame(edges, dyn))
  expect_equal(d$roles$type, "is_event")
  expect_equal(d$nodes$is_event, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(d$nodes$type, dyn$type)
  expect_error(net_twomode(edges, nodes, type = "kind", first = "robot"), "one of")
  expect_error(net_twomode(edges, nodes, type = "nope"), "not found")
  three <- transform(nodes, kind = c("a", "b", "c", "c"))
  expect_error(net_twomode(edges, three, type = "kind"), "exactly two")
  within <- data.frame(source = c("p", "p"), target = c("x", "q"))
  expect_error(net_twomode(within, nodes, type = "kind"), "same mode")
})

test_that("net_as_twomode() converts matrices, data frames, and netframes", {
  m <- matrix(c(1, 0, 1,
                1, 1, 0), 2, byrow = TRUE, dimnames = list(c("p", "q"), c("x", "y", "z")))
  g <- net_as_twomode(m)
  expect_equal(g$nodes$name, c("p", "q", "x", "y", "z"))
  expect_equal(g$nodes$type, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(nrow(g$edges), 4L)
  expect_null(g$roles$weight)
  expect_equal(unname(net_incidence(g)), unname(m))
  m[1, 1] <- 3
  w <- net_as_twomode(m)
  expect_equal(w$roles$weight, "weight")
  expect_equal(net_incidence(w, attr = "weight"), m)
  expect_equal(net_as_twomode(unname(m))$nodes$name, c("r1", "r2", "c1", "c2", "c3"))
  e <- net_as_twomode(data.frame(a = c("p", "q"), b = c("x", "x")))
  expect_equal(e$nodes$type, c(FALSE, FALSE, TRUE))
  expect_identical(net_as_twomode(e), e)
  expect_error(net_as_twomode(net_frame(data.frame(source = 1, target = 2))), "Cannot tell")
})

test_that("a DyNetML file with two node classes becomes two-mode in one call", {
  nodes <- data.frame(name = c("ann", "bob", "gala", "fair"),
                      class = c("agent", "agent", "event", "event"),
                      type = c("Agent", "Agent", "Event", "Event"))
  g <- net_frame(data.frame(source = c("ann", "bob", "ann"), target = c("gala", "gala", "fair")),
                nodes, directed = FALSE)
  path <- tempfile(fileext = ".xml")
  net_write(g, path, format = "dynetml")
  h <- net_as_twomode(net_read(path))
  expect_true(net_is_twomode(h))
  expect_equal(h$nodes[[h$roles$type]], c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(h$roles$type_labels, c("agent", "event"))
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
  expect_s3_class(p, "netframe")
  expect_false(p$directed)
  expect_equal(p$nodes$name, c("ann", "bob", "cat", "dan", "eve"))
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
  plain <- net_twomode(data.frame(a = "p", b = "x"))
  expect_error(net_project(plain, method = "product"), "weights")
  expect_equal(nrow(net_project(plain)$edges), 0L)
  b <- net_incidence(g)
  full <- b %*% t(b)
  pf <- net_project(g, loops = TRUE)
  ids <- p$nodes$name
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
  g <- net_twomode(data.frame(a = c("p", "q"), b = c("x", "y")))
  # Each node reaches one other at distance 1: (2 + 2 * 1) / 1, times 1/3
  # under the Wasserman-Faust rule.
  expect_equal(net_twomode_closeness(g), rep(4 / 3, 4))
  expect_equal(net_twomode_closeness(g, unreachable = "ignore"), rep(4, 4))
  expect_warning(net_twomode_eigen(g), "not connected")
  lone <- net_twomode(data.frame(a = "p", b = "x"),
                  data.frame(name = c("p", "x", "z"), type = c(FALSE, TRUE, FALSE)))
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
  back <- net_as_twomode(ig)
  expect_true(net_is_twomode(back))
  expect_equal(back$nodes$type, g$nodes$type)
  bip <- igraph::make_bipartite_graph(c(FALSE, FALSE, TRUE), c(1, 3, 2, 3))
  expect_equal(net_as_twomode(bip)$nodes$type, c(FALSE, FALSE, TRUE))
})

test_that("bipartite networks from the network package", {
  skip_if_not_installed("network")
  net <- network::network.initialize(4, directed = FALSE, bipartite = 2)
  network::add.edges(net, tail = c(1, 1, 2), head = c(3, 4, 3))
  network::set.vertex.attribute(net, "vertex.names", c("p", "q", "x", "y"))
  g <- net_as_twomode(net)
  expect_true(net_is_twomode(g))
  expect_equal(g$roles$type_labels, c("1", "2"))
  expect_equal(g$nodes$type, c(FALSE, FALSE, TRUE, TRUE))
})
