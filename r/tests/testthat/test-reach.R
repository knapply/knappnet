# Reach, proximity prestige, edge betweenness, information centrality,
# bridging centrality, hierarchy, eigenvector centralization, and
# centralize(). Hand-checked values first, then spot checks against igraph
# and sna. The full comparison over every fixture lives in validation/.

kite <- function() {
  e <- data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
  netframe(e, nodes = data.frame(id = 1:10), directed = FALSE)
}

path5 <- function() {
  netframe(data.frame(from = 1:4, to = 2:5), nodes = data.frame(id = 1:5), directed = FALSE)
}

test_that("reach counts nodes within k steps", {
  g <- path5()
  expect_equal(reach(g, k = 1), c(1, 2, 2, 2, 1))
  expect_equal(reach(g, k = 2), c(2, 3, 4, 3, 2))
  expect_equal(reach(g), reach(g, k = 2))
  expect_equal(reach(g, k = Inf), rep(4, 5))
  expect_equal(reach(g, k = Inf), reachability(g))
  expect_equal(reach(g, k = 2, normalized = TRUE), c(2, 3, 4, 3, 2) / 4)
  expect_equal(reach(g, k = 0), rep(0, 5))
  # Direction: 1 -> 2 -> 3.
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(reach(d, k = 2), c(2, 1, 0))
  expect_equal(reach(d, k = 2, mode = "in"), c(0, 1, 2))
  expect_equal(reach(d, k = 2, mode = "all"), c(2, 2, 2))
  # Weighted: k is a length.
  w <- netframe(data.frame(from = c(1, 2), to = c(2, 3), w = c(0.5, 2)), directed = FALSE,
                weight = "w")
  expect_equal(reach(w, k = 1), c(1, 1, 0))
  expect_equal(reach(w, k = 1, weighted = FALSE), c(1, 2, 1))
  expect_equal(reach(w, k = 1, weights_as = "strength"), c(0, 1, 1))
  expect_error(reach(g, k = -1), "`k`")
  expect_error(reach(g, k = "two"), "`k`")
  expect_error(reach(g, normalized = NA), "`normalized`")
})

test_that("distance-weighted reach is normalized harmonic centrality", {
  g <- path5()
  expect_equal(distance_weighted_reach(g), harmonic_centrality(g, normalized = TRUE))
  # Node 3 of the path: (1 + 1/2 + 1/2 + 1) / 4.
  expect_equal(distance_weighted_reach(g)[3], 0.75)
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(distance_weighted_reach(d, mode = "in"), c(0, 0.5, 0.75))
})

test_that("proximity prestige follows Wasserman and Faust", {
  # 1 -> 2 -> 3 and 4 -> 3: node 3 is reached by all three at distances
  # 2, 1, 1, so (3/3) / (4/3) = 3/4; node 2 by node 1 only: (1/3) / 1.
  g <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 3)))
  expect_equal(proximity_prestige(g), c(0, 1 / 3, 3 / 4, 0))
  # Undirected: equal to the Wasserman-Faust closeness index.
  u <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)), directed = FALSE)
  expect_equal(proximity_prestige(u), closeness(u, convention = "networkx"))
})

test_that("edge betweenness matches hand values", {
  # Path 1 - 2 - 3: each edge carries its own pair and the pair (1, 3).
  g <- netframe(data.frame(from = c(1, 2), to = c(2, 3)), directed = FALSE)
  expect_equal(edge_betweenness(g), c(2, 2))
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(edge_betweenness(d), c(2, 2))
  expect_equal(edge_betweenness(d, normalized = TRUE), c(2, 2) / 6)
  # A square: each edge carries its own pair and half of two diagonals.
  sq <- netframe(data.frame(from = c(1, 2, 3, 4), to = c(2, 3, 4, 1)), directed = FALSE)
  expect_equal(edge_betweenness(sq), rep(2, 4))
  expect_equal(edge_betweenness(sq, normalized = TRUE), rep(2 / 6, 4))
  # Parallel edges split the paths; a self-loop carries none.
  p <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 2, 3, 3)), directed = FALSE)
  expect_equal(edge_betweenness(p), c(1, 1, 2, 0))
  # Weights: the long edge of a triangle is bypassed.
  w <- netframe(data.frame(from = c(1, 2, 1), to = c(2, 3, 3), w = c(1, 1, 5)),
                directed = FALSE, weight = "w")
  expect_equal(edge_betweenness(w), c(2, 2, 0))
  expect_equal(edge_betweenness(w, weighted = FALSE), c(1, 1, 1))
  expect_equal(edge_betweenness(w, weights_as = "strength"), c(1, 1, 1))
  # Values sum to the total of all distances in a connected network.
  k <- kite()
  expect_equal(sum(edge_betweenness(k)), sum(distances(k)[upper.tri(distances(k))]))
  expect_error(edge_betweenness(g, normalized = "yes"), "`normalized`")
})

test_that("information centrality matches the formula", {
  # Star with centre 1: I = n / sum of effective resistances.
  star <- netframe(data.frame(from = 1, to = 2:4), directed = FALSE)
  expect_equal(information_centrality(star), c(4 / 3, 0.8, 0.8, 0.8))
  expect_equal(information_centrality(star, convention = "networkx"), c(4 / 3, 0.8, 0.8, 0.8) / 4)
  # Isolates score 0 and are left out of n.
  iso <- netframe(data.frame(from = 1, to = 2:4), nodes = data.frame(id = 1:5),
                  directed = FALSE)
  expect_equal(information_centrality(iso), c(4 / 3, 0.8, 0.8, 0.8, 0))
  expect_error(information_centrality(iso, convention = "networkx"), "not connected")
  # Two components with ties: undefined.
  two <- netframe(data.frame(from = c(1, 3), to = c(2, 4)), directed = FALSE)
  expect_true(all(is.nan(information_centrality(two))))
  expect_error(information_centrality(two, convention = "networkx"), "not connected")
  # Weights are conductances: doubling them doubles every score.
  e <- data.frame(from = c(1, 2, 1, 3), to = c(2, 3, 3, 4), w = c(1, 2, 3, 1.5))
  g1 <- netframe(e, directed = FALSE, weight = "w")
  g2 <- netframe(transform(e, w = 2 * w), directed = FALSE, weight = "w")
  expect_equal(information_centrality(g2), 2 * information_centrality(g1))
  expect_equal(information_centrality(g1, weighted = FALSE),
               information_centrality(netframe(e, directed = FALSE)))
  bad <- netframe(data.frame(from = 1, to = 2, w = -1), directed = FALSE, weight = "w")
  expect_error(information_centrality(bad), "strictly positive")
  # Direction is ignored and parallel edges collapse.
  d <- netframe(data.frame(from = c(1, 2, 1, 4, 4), to = c(2, 1, 3, 1, 1)))
  expect_equal(information_centrality(d), c(4 / 3, 0.8, 0.8, 0.8))
})

test_that("bridging centrality is betweenness times the bridging coefficient", {
  # Two triangles joined through node 4.
  g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 5, 5, 6), to = c(2, 3, 3, 4, 5, 6, 7, 7)),
                directed = FALSE)
  expect_equal(bridging_coefficient(g), c(0.6, 0.6, 2 / 9, 0.75, 2 / 9, 0.6, 0.6))
  expect_equal(bridging_centrality(g), betweenness(g) * bridging_coefficient(g))
  expect_equal(bridging_centrality(g)[3:5], c(16 / 9, 6.75, 16 / 9))
  # Isolates are NaN; parallel edges and loops do not count.
  h <- netframe(data.frame(from = c(1, 2, 2, 3), to = c(2, 3, 3, 3)),
                nodes = data.frame(id = 1:4), directed = FALSE)
  expect_equal(bridging_coefficient(h)[1:3], c(2, 0.25, 2))
  expect_true(is.nan(bridging_coefficient(h)[4]))
  expect_true(is.nan(bridging_centrality(h)[4]))
})

test_that("hierarchy follows Burt's formula", {
  star <- netframe(data.frame(from = 1, to = 2:4), directed = FALSE)
  h <- hierarchy(star)
  expect_equal(h[1], 0)
  expect_true(all(is.nan(h[2:4])))
  # Node 1 tied to 2, 3, 4 with 3 - 4 tied: c = (1/9, 1/4, 1/4).
  g <- netframe(data.frame(from = c(1, 1, 1, 3), to = c(2, 3, 4, 4)), directed = FALSE)
  cc <- c(1 / 9, 1 / 4, 1 / 4)
  r <- cc / mean(cc)
  expect_equal(hierarchy(g)[1], sum(r * log(r)) / (3 * log(3)))
  expect_equal(hierarchy(g)[1], 0.0522045454467525, tolerance = 1e-12)
  expect_equal(hierarchy(g, threads = 2), hierarchy(g))
  expect_true(all(is.nan(hierarchy(netframe(data.frame(from = 1, to = 2), directed = FALSE)))))
})

test_that("eigenvector centralization follows igraph's definition", {
  star10 <- netframe(data.frame(from = 1, to = 2:10), directed = FALSE)
  expect_equal(centralization(star10, "eigenvector"), 0.75, tolerance = 1e-9)
  expect_equal(centralization(star10, "eigenvector", normalized = FALSE), 6, tolerance = 1e-9)
  one <- netframe(data.frame(from = 1, to = 2), nodes = data.frame(id = 1:10),
                  directed = FALSE)
  expect_equal(centralization(one, "eigenvector"), 1, tolerance = 1e-9)
  cycle <- netframe(data.frame(from = 1:5, to = c(2:5, 1)), directed = FALSE)
  expect_equal(centralization(cycle, "eigenvector"), 0, tolerance = 1e-9)
  expect_true(is.nan(centralization(netframe(data.frame(from = 1, to = 2), directed = FALSE),
                                    "eigenvector")))
  expect_error(centralization(star10, "eigenvector", convention = "sna"), "igraph")
  expect_error(centralization(star10, "eigenvector", convention = "ucinet"), "igraph")
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_error(centralization(d, "eigenvector", mode = "in"), "direction ignored")
  two <- netframe(data.frame(from = c(1, 3), to = c(2, 4)), directed = FALSE)
  expect_warning(centralization(two, "eigenvector"), "component")
})

test_that("centralize() applies Freeman's index to any scores", {
  star <- netframe(data.frame(from = 1, to = 2:5), directed = FALSE)
  expect_equal(centralize(degree(star), 4 * 3), 1)
  expect_equal(centralize(degree(star), 4 * 3), centralization(star))
  expect_equal(centralize(degree(star)), 12)
  expect_equal(centralize(c(2, 2, 2)), 0)
  expect_true(is.nan(centralize(numeric(0))))
  expect_true(is.nan(centralize(c(1, NaN))))
  expect_error(centralize("a"), "numeric")
  expect_error(centralize(1:3, theoretical_max = 0), "positive")
})

test_that("results match igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(edge_betweenness(g), igraph::edge_betweenness(ig))
  expect_equal(centralization(g, "eigenvector"), igraph::centr_eigen(ig)$centralization,
               tolerance = 1e-8)
  expect_equal(reach(g, k = 2), igraph::ego_size(ig, order = 2) - 1)
  expect_equal(distance_weighted_reach(g),
               unname(igraph::harmonic_centrality(ig, normalized = TRUE)))

  set.seed(3)
  e <- data.frame(from = sample(30, 120, TRUE), to = sample(30, 120, TRUE))
  e <- e[e$from != e$to, ]
  e <- e[!duplicated(e), ]
  e$w <- round(runif(nrow(e), 0.5, 3), 2)
  nodes <- data.frame(id = 1:30)
  d <- netframe(e, nodes, weight = "w")
  id <- igraph::graph_from_data_frame(e, vertices = nodes)
  expect_equal(edge_betweenness(d), igraph::edge_betweenness(id, weights = e$w),
               tolerance = 1e-9)
  expect_equal(edge_betweenness(d, weighted = FALSE), igraph::edge_betweenness(id, weights = NA),
               tolerance = 1e-9)
  expect_equal(reach(d, k = 2, mode = "in"), igraph::ego_size(id, order = 2, mode = "in") - 1)
})

test_that("results match sna", {
  skip_if_not_installed("sna")
  g <- kite()
  m <- as.matrix(g)
  expect_equal(information_centrality(g), unname(sna::infocent(m, gmode = "graph")),
               tolerance = 1e-10)
  set.seed(4)
  e <- data.frame(from = sample(12, 40, TRUE), to = sample(12, 40, TRUE))
  e <- e[e$from != e$to, ]
  e <- e[!duplicated(e), ]
  d <- netframe(e, nodes = data.frame(id = 1:12))
  md <- as.matrix(d)
  if (all(components(d)$csize > 1L) && components(d)$no == 1L) {
    expect_equal(information_centrality(d), unname(sna::infocent(md, gmode = "digraph")),
                 tolerance = 1e-10)
  }
})
