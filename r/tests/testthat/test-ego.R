# Ego networks and egocentric measures. Hand-checked values first, then
# spot checks against igraph. The full comparison lives in validation/.

# Edges a-b, a-c, a-d, b-c, d-e. Nodes come in order of first appearance:
# a, b, d, c, e.
small <- function() {
  netframe(data.frame(from = c("a", "a", "a", "b", "d"),
                      to   = c("b", "c", "d", "c", "e"),
                      w    = c(1, 2, 3, 4, 5)),
           directed = FALSE, weight = "w", weight_type = "strength")
}

test_that("ego_summary gives hand-checked values", {
  g <- small()
  expect_equal(g$nodes$name, c("a", "b", "d", "c", "e"))
  s <- ego_summary(g)
  expect_equal(names(s), c("node", "size", "ties", "pairs", "density", "components",
                           "ego_betweenness"))
  expect_equal(s$node, g$nodes$name)
  expect_equal(s$size, c(3, 2, 2, 2, 1))
  expect_equal(s$ties, c(1, 1, 0, 1, 0))
  expect_equal(s$pairs, c(3, 1, 1, 1, 0))
  expect_equal(s$density, c(1 / 3, 1, 0, 1, NaN))
  expect_equal(s$components, c(2, 1, 2, 1, 1))
  # a brokers b-d and c-d; d brokers a-e.
  expect_equal(s$ego_betweenness, c(2, 0, 1, 0, 0))
  expect_equal(ego_density(g), s$density)
  expect_equal(ego_betweenness(g), s$ego_betweenness)
  expect_false("ego_betweenness" %in% names(ego_summary(g, betweenness = FALSE)))
  expect_equal(ego_summary(g, threads = 1L), ego_summary(g, threads = 3L))
})

test_that("ego_size counts ego unless asked not to", {
  g <- small()
  expect_equal(ego_size(g), c(4, 3, 3, 3, 2))
  expect_equal(ego_size(g, include_ego = FALSE), c(3, 2, 2, 2, 1))
  expect_equal(ego_size(g, order = 2), c(5, 4, 5, 4, 3))
  expect_equal(ego_size(g, order = 0), rep(1, 5))
  expect_error(ego_size(g, order = -1), "order")
  expect_error(ego_summary(g, order = 0), "at least 1")
})

test_that("directed alters follow mode and ties keep their direction", {
  d <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 3, 1)))
  out <- ego_summary(d, mode = "out")
  expect_equal(out$size, c(2, 1, 1))
  expect_equal(out$ties[1], 1)
  expect_equal(out$pairs[1], 2)
  expect_equal(out$density[1], 0.5)
  # Only 3 -> 1 -> 2 passes through node 1.
  expect_equal(out$ego_betweenness[1], 1)
  expect_equal(ego_summary(d, mode = "in")$size, c(1, 1, 2))
  expect_equal(ego_size(d, mode = "in", include_ego = FALSE), c(1, 1, 2))
  expect_equal(ego_summary(d)$size, c(2, 2, 2))
})

test_that("weights, parallel edges, and loops do not change the measures", {
  g <- small()
  plain <- netframe(rbind(g$edges[, 1:2], data.frame(from = c("a", "c"), to = c("b", "c"))),
                    directed = FALSE)
  expect_equal(plain$nodes$name, g$nodes$name)
  expect_equal(ego_summary(g)[, -1], ego_summary(plain)[, -1])
})

test_that("ego_network returns the induced subgraph with its attributes", {
  g <- small()
  g$nodes$age <- c(30, 40, 50, 60, 70)
  e <- ego_network(g, "a")
  expect_s3_class(e, "netframe")
  expect_equal(e$nodes$name, c("a", "b", "d", "c"))
  expect_equal(e$nodes$age, c(30, 40, 50, 60))
  expect_equal(e$edges$w, c(1, 2, 3, 4))
  expect_equal(e$roles$weight, "w")
  expect_equal(e$roles$weight_type, "strength")
  expect_false(e$directed)
  alters <- ego_network(g, "a", include_ego = FALSE)
  expect_equal(alters$nodes$name, c("b", "d", "c"))
  expect_equal(nrow(alters$edges), 1L)
  expect_equal(nrow(ego_network(g, "e", order = 2)$nodes), 3L)
  all <- ego_networks(g)
  expect_equal(names(all), g$nodes$name)
  expect_equal(unname(vapply(all, function(x) as.numeric(nrow(x$nodes)), 0)), ego_size(g))
  some <- ego_networks(g, c("e", "a"))
  expect_equal(names(some), c("e", "a"))
  expect_error(ego_network(g, "zz"), "not found")
  expect_error(ego_network(g, c("a", "b")), "single")
})

test_that("ego networks of a two-mode network stay two-mode", {
  g <- twomode(data.frame(p = c("ann", "ann", "bob"), e = c("x", "y", "y")))
  e <- ego_network(g, "y", order = 2)
  expect_true(is_twomode(e))
  expect_equal(e$nodes$type, c(FALSE, FALSE, TRUE, TRUE))
})

test_that("results match igraph", {
  skip_if_not_installed("igraph")
  ig <- igraph::make_graph("Zachary")
  g <- as_netframe(ig)
  expect_equal(ego_size(g), igraph::ego_size(ig, order = 1))
  expect_equal(ego_size(g, order = 2, include_ego = FALSE),
               igraph::ego_size(ig, order = 2, mindist = 1))
  s <- ego_summary(g)
  hoods <- igraph::ego(ig, order = 1, mindist = 1)
  ties <- vapply(hoods, function(h) {
    as.numeric(igraph::ecount(igraph::induced_subgraph(ig, h)))
  }, 0)
  expect_equal(s$ties, ties)
  comps <- vapply(hoods, function(h) {
    as.numeric(igraph::count_components(igraph::induced_subgraph(ig, h)))
  }, 0)
  expect_equal(s$components, comps)
  igraph::V(ig)$name <- as.character(seq_len(igraph::vcount(ig)))
  eb <- vapply(seq_len(igraph::vcount(ig)), function(v) {
    sub <- igraph::induced_subgraph(ig, igraph::ego(ig, order = 1, nodes = v)[[1]])
    as.numeric(igraph::betweenness(sub, v = as.character(v)))
  }, 0)
  expect_equal(s$ego_betweenness, eb)

  set.seed(5)
  dg <- igraph::sample_gnp(30, 0.1, directed = TRUE)
  igraph::V(dg)$name <- as.character(seq_len(30))
  d <- as_netframe(dg)
  expect_equal(ego_size(d, mode = "out"), igraph::ego_size(dg, mode = "out"))
  expect_equal(ego_size(d, mode = "in"), igraph::ego_size(dg, mode = "in"))
  e <- ego_network(d, "1", mode = "in")
  expect_equal(sort(e$nodes$name),
               sort(names(igraph::ego(dg, order = 1, nodes = "1", mode = "in")[[1]])))
})
