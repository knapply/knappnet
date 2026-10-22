# Cohesive subgroups. The full validation (every fixture, every oracle)
# lives in validation/; these tests check the R bindings end to end.

barbell <- function() {
  # Two triangles joined by the edge c-d; nodes appear in order a to f.
  netframe(data.frame(from = c("a", "a", "b", "c", "d", "d", "e"),
                      to   = c("b", "c", "c", "d", "e", "f", "f")),
           directed = FALSE)
}

kite_edges <- function() {
  data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
}

test_that("coreness follows the definition and its conventions", {
  g <- netframe(data.frame(from = c(1, 2, 3, 3), to = c(2, 3, 1, 4)),
                nodes = data.frame(id = 1:5), directed = FALSE)
  expect_equal(coreness(g), c(2, 2, 2, 1, 0))
  # A loop adds 2 in an undirected network; a doubled edge is a 2-core.
  loop <- netframe(data.frame(from = c(1, 1), to = c(1, 2)), directed = FALSE)
  expect_equal(coreness(loop), c(2, 1))
  d <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 3, 1)))
  expect_equal(coreness(d, "all"), c(2, 2, 2))
  expect_equal(coreness(d, "out"), c(1, 1, 1))
  expect_error(coreness(d, "both"))
})

test_that("maximal cliques return node identifiers", {
  g <- netframe(data.frame(from = c("a", "b", "a", "c"), to = c("b", "c", "c", "d")),
                nodes = data.frame(id = c("a", "b", "c", "d", "e")), directed = FALSE)
  expect_equal(max_cliques(g), list(c("a", "b", "c"), c("c", "d"), "e"))
  expect_equal(max_cliques(g, min_size = 3), list(c("a", "b", "c")))
  expect_identical(clique_number(g), 3L)
  expect_error(max_cliques(g, min_size = -1), "non-negative")
  # Direction, loops, and repeated edges do not matter.
  d <- netframe(data.frame(from = c(1, 2, 3, 3, 1), to = c(2, 3, 1, 3, 2)))
  expect_equal(max_cliques(d), list(c(1, 2, 3)))
})

test_that("modularity matches hand-computed values", {
  g <- barbell()
  # NetworkX documents 5/14 for this split of barbell_graph(3, 0).
  expect_equal(modularity(g, c(1, 1, 1, 2, 2, 2)), 5 / 14, tolerance = 1e-14)
  expect_equal(modularity(g, c("x", "x", "x", "y", "y", "y")), 5 / 14, tolerance = 1e-14)
  expect_equal(modularity(g, rep(1, 6), resolution = 0.25), 0.75, tolerance = 1e-14)
  expect_error(modularity(g, 1:5), "length")
  expect_error(modularity(g, c(1, 1, 1, 2, 2, NA)), "missing")
  expect_error(modularity(g, rep(1, 6), resolution = -1), "non-negative")
  empty <- netframe(data.frame(from = integer(), to = integer()),
                    nodes = data.frame(id = 1:3), directed = FALSE)
  expect_true(is.nan(modularity(empty, 1:3)))
})

test_that("Louvain is reproducible and reports its own modularity", {
  g <- barbell()
  cl <- cluster_louvain(g, seed = 7)
  expect_equal(as.vector(cl), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(attr(cl, "modularity"), modularity(g, cl))
  expect_identical(cluster_louvain(g, seed = 3), cluster_louvain(g, seed = 3))
  set.seed(1)
  a <- cluster_louvain(g)
  set.seed(1)
  expect_identical(cluster_louvain(g), a)
  expect_error(cluster_louvain(netframe(data.frame(from = 1, to = 2))), "undirected")
  expect_error(cluster_louvain(g, seed = -1), "whole number")
})

test_that("weights are used when present", {
  e <- data.frame(from = c(1, 1, 2, 3, 4, 4, 5), to = c(2, 3, 3, 4, 5, 6, 6),
                  w = c(1, 1, 1, 50, 1, 1, 1))
  g <- netframe(e, directed = FALSE, weight = "w")
  cl <- cluster_louvain(g, seed = 1)
  expect_equal(cl[[3]], cl[[4]])
  expect_equal(attr(cl, "modularity"), modularity(g, cl))
  expect_false(isTRUE(all.equal(modularity(g, c(1, 1, 1, 2, 2, 2)),
                                modularity(g, c(1, 1, 1, 2, 2, 2), weighted = FALSE))))
  expect_error(modularity(barbell(), rep(1, 6), weighted = TRUE), "no weights")
})

test_that("subgroups match igraph", {
  skip_if_not_installed("igraph")
  e <- kite_edges()
  g <- netframe(e, directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = g$nodes)
  expect_equal(coreness(g), unname(as.numeric(igraph::coreness(ig))))
  ours <- lapply(max_cliques(g), function(c) sort(match(c, g$nodes$name)))
  theirs <- lapply(igraph::max_cliques(ig), function(c) sort(as.integer(c)))
  key <- function(l) sort(vapply(l, paste, character(1), collapse = " "))
  expect_equal(key(ours), key(theirs))
  expect_equal(clique_number(g), igraph::clique_num(ig))
  memb <- rep(1:3, length.out = nrow(g$nodes))
  expect_equal(modularity(g, memb), igraph::modularity(ig, memb), tolerance = 1e-12)
  expect_equal(modularity(g, memb, resolution = 0.5),
               igraph::modularity(ig, memb, resolution = 0.5), tolerance = 1e-12)
  cl <- cluster_louvain(g, seed = 1)
  expect_equal(attr(cl, "modularity"), igraph::modularity(ig, as.vector(cl)),
               tolerance = 1e-12)

  # Directed coreness modes and directed modularity.
  set.seed(3)
  d <- data.frame(from = sample(25, 90, TRUE), to = sample(25, 90, TRUE))
  nodes <- data.frame(id = 1:25)
  gd <- netframe(d, nodes)
  igd <- igraph::graph_from_data_frame(d, vertices = nodes)
  for (mode in c("all", "out", "in")) {
    expect_equal(coreness(gd, mode), unname(as.numeric(igraph::coreness(igd, mode = mode))))
  }
  memb <- rep(1:4, length.out = 25)
  expect_equal(modularity(gd, memb), igraph::modularity(igd, memb, directed = TRUE),
               tolerance = 1e-12)
})
