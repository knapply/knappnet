# Spot checks against igraph and sna on classic data. The full validation
# suite (many fixtures, every oracle) lives in validation/ at the repository
# root; these tests make sure the R bindings agree with it end to end.

kite_edges <- function() {
  # Krackhardt's kite as an undirected edge list.
  data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
}

test_that("betweenness matches igraph", {
  skip_if_not_installed("igraph")
  e <- kite_edges()
  g <- netframe(e, directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = g$nodes)
  expect_equal(betweenness(g), unname(igraph::betweenness(ig)), tolerance = 1e-10)
  expect_equal(betweenness(g, normalized = TRUE),
               unname(igraph::betweenness(ig, normalized = TRUE)), tolerance = 1e-10)
})

test_that("triad census matches igraph and sna", {
  skip_if_not_installed("igraph")
  set.seed(1)
  e <- data.frame(from = sample(30, 120, TRUE), to = sample(30, 120, TRUE))
  e <- e[e$from != e$to, ]
  e <- e[!duplicated(e), ]
  g <- netframe(e, nodes = data.frame(id = 1:30))
  ig <- igraph::graph_from_data_frame(e, vertices = data.frame(id = 1:30))
  expect_equal(unname(triad_census(g)), igraph::triad_census(ig))

  skip_if_not_installed("sna")
  m <- matrix(0, 30, 30)
  m[cbind(e$from, e$to)] <- 1
  expect_equal(unname(triad_census(g)), unname(drop(sna::triad.census(m))))
})

test_that("degree matches sna when loops are ignored", {
  skip_if_not_installed("sna")
  e <- data.frame(from = c(1, 1, 2, 3), to = c(2, 1, 3, 1))
  g <- netframe(e, nodes = data.frame(id = 1:3))
  m <- matrix(0, 3, 3)
  m[cbind(e$from, e$to)] <- 1
  expect_equal(degree(g, "out", loops = FALSE), sna::degree(m, cmode = "outdegree"))
  expect_equal(degree(g, "in", loops = FALSE), sna::degree(m, cmode = "indegree"))
})
