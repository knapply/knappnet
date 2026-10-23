# Degree and strength, betweenness, closeness and harmonic centrality,
# eigenvector centrality, reach, proximity prestige, edge betweenness, and
# information and bridging centrality (R/centrality.R). Hand-checked values
# and the convention switch first, then spot checks against igraph, sna,
# and (through published numbers) UCINET. The full comparison over every
# fixture and oracle lives in validation/ at the repository root.

path5 <- function() {
  netframe(data.frame(from = 1:4, to = 2:5), nodes = data.frame(id = 1:5), directed = FALSE)
}

# A triangle where the direct edge a-c has weight 4 and the path through b
# has two edges of weight 1. Read as distances, the path through b is
# shorter; read as strengths (lengths 1/w), the direct tie is shorter.
triangle <- function(weight_type = "distance") {
  netframe(data.frame(from = c("a", "b", "a"), to = c("b", "c", "c"), w = c(1, 1, 4)),
           directed = FALSE, weight = "w", weight_type = weight_type)
}

test_that("conventions leave defaults of degree and betweenness unchanged", {
  e <- data.frame(from = c(1, 1, 2, 3, 3), to = c(2, 2, 3, 3, 1))
  g <- netframe(e, nodes = data.frame(id = 1:3))
  expect_equal(degree(g), c(3, 3, 4))
  expect_equal(degree(g, convention = "sna"), c(2, 2, 2))
  expect_equal(degree(g, convention = "networkx"), degree(g))
  # From 1 to 4 there are two shortest paths through 2 (parallel edges) and
  # one through 3 under igraph's convention, and one of each under the others.
  diamond <- netframe(data.frame(from = c(1, 1, 2, 1, 3), to = c(2, 2, 4, 3, 4)),
                      nodes = data.frame(id = 1:4))
  expect_equal(betweenness(diamond), c(0, 2 / 3, 1 / 3, 0))
  expect_equal(betweenness(diamond, convention = "networkx"), c(0, 1 / 2, 1 / 2, 0))
  expect_equal(betweenness(diamond, convention = "sna"), c(0, 1 / 2, 1 / 2, 0))
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

test_that("betweenness matches igraph", {
  skip_if_not_installed("igraph")
  e <- kite_edges()
  g <- netframe(e, directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = g$nodes)
  expect_equal(betweenness(g), unname(igraph::betweenness(ig)), tolerance = 1e-10)
  expect_equal(betweenness(g, normalized = TRUE),
               unname(igraph::betweenness(ig, normalized = TRUE)), tolerance = 1e-10)
})

test_that("binary networks are reported and analysed as unweighted", {
  g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")))
  expect_null(g$roles$weight)
  expect_output(print(g), "unweighted")
  expect_equal(betweenness(g), c(0, 1, 0))
  expect_error(strength(g), "no weights")
})

test_that("weighted networks use their weights by default", {
  g <- triangle()
  expect_output(print(g), "weighted by `w` \\(distances\\)")
  expect_equal(betweenness(g), c(0, 1, 0))
  expect_equal(betweenness(g, weighted = FALSE), c(0, 0, 0))
})

test_that("strength weights are inverted for paths but summed as given", {
  g <- triangle("strength")
  expect_equal(betweenness(g), c(0, 0, 0))
  expect_equal(betweenness(triangle(), weights_as = "strength"), c(0, 0, 0))
  expect_equal(strength(g), c(5, 2, 5))
  expect_equal(degree(g), c(2, 2, 2))
})

test_that("every weighting matches igraph", {
  skip_if_not_installed("igraph")
  set.seed(42)
  e <- data.frame(from = sample(40, 200, TRUE), to = sample(40, 200, TRUE))
  e <- e[e$from != e$to, ]
  e$int_w <- sample(1:5, nrow(e), TRUE)
  e$real_w <- round(runif(nrow(e), 0.1, 3), 3)
  nodes <- data.frame(id = 1:40)
  for (col in c("int_w", "real_w")) {
    g <- netframe(e, nodes, weight = col)
    ig <- igraph::graph_from_data_frame(e, vertices = nodes)
    w <- e[[col]]
    expect_equal(betweenness(g), unname(igraph::betweenness(ig, weights = w)))
    expect_equal(betweenness(g, weights_as = "strength"),
                 unname(igraph::betweenness(ig, weights = 1 / w)))
    expect_equal(strength(g, "out"), unname(igraph::strength(ig, mode = "out", weights = w)))
  }
})

test_that("closeness conventions on a disconnected network", {
  g <- pieces()
  expect_equal(closeness(g), c(1 / 3, 1 / 2, 1 / 3, 1, 1, NaN))
  expect_equal(closeness(g, normalized = TRUE), c(2 / 3, 1, 2 / 3, 1, 1, NaN))
  wf <- (2 / 3) * (2 / 5)
  expect_equal(closeness(g, convention = "networkx"), c(wf, 2 / 5, wf, 1 / 5, 1 / 5, 0))
  expect_equal(closeness(g, convention = "sna"), rep(0, 6))
  # The longest path has length 2, so missing lengths count as 3.
  expect_equal(closeness(g, convention = "ucinet"),
               c(5 / 12, 5 / 11, 5 / 12, 5 / 13, 5 / 13, 5 / 15))
  expect_error(closeness(g, normalized = FALSE, convention = "networkx"), "normalized")
  expect_equal(harmonic_centrality(g), c(1.5, 2, 1.5, 1, 1, 0))
  expect_equal(harmonic_centrality(g, convention = "sna"), c(1.5, 2, 1.5, 1, 1, 0) / 5)
})

test_that("closeness reproduces UCINET's published output", {
  # UCINET 6.849's Closeness report for its g9_disc example (a star, a path,
  # and a triangle), printed to three decimals, as published with the
  # authors' xucinet package (inst/goldens/centrality/log_menu.txt).
  e <- data.frame(from = c(1, 1, 1, 1, 6, 7, 8, 10, 11, 10),
                  to   = c(2, 3, 4, 5, 7, 8, 9, 11, 12, 12))
  g <- netframe(e, nodes = data.frame(id = 1:12), directed = FALSE)
  free <- c(0.344, 0.314, 0.314, 0.314, 0.314, 0.289, 0.306, 0.306, 0.289, 0.289,
            0.289, 0.289)
  recip <- c(0.364, 0.227, 0.227, 0.227, 0.227, 0.167, 0.227, 0.227, 0.167, 0.182,
             0.182, 0.182)
  # Within half a unit in the last printed place.
  expect_true(all(abs(closeness(g, convention = "ucinet") - free) <= 5e-4))
  expect_true(all(abs(harmonic_centrality(g, convention = "ucinet") - recip) <= 5e-4))
})

test_that("directed closeness modes and weights", {
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3), w = c(2, 3)), weight = "w")
  expect_equal(closeness(d), c(1 / 7, 1 / 3, NaN))
  expect_equal(closeness(d, weighted = FALSE), c(1 / 3, 1, NaN))
  expect_equal(closeness(d, mode = "in"), c(NaN, 1 / 2, 1 / 8))
  expect_equal(closeness(d, weights_as = "strength"), c(1 / (1 / 2 + 5 / 6), 3, NaN))
  expect_error(closeness(d, weighted = TRUE, convention = "ucinet"), "binary")
})

test_that("eigenvector centrality closed forms and conventions", {
  star <- netframe(data.frame(from = 1, to = 2:5), directed = FALSE)
  expect_equal(eigenvector_centrality(star), c(1, 0.5, 0.5, 0.5, 0.5), tolerance = 1e-10)
  u <- eigenvector_centrality(star, convention = "networkx")
  expect_equal(sum(u^2), 1, tolerance = 1e-12)
  expect_equal(u / max(u), c(1, 0.5, 0.5, 0.5, 0.5), tolerance = 1e-10)
  path <- netframe(data.frame(from = 1:3, to = 2:4), directed = FALSE)
  phi <- (1 + sqrt(5)) / 2
  expect_equal(eigenvector_centrality(path), c(1 / phi, 1, 1, 1 / phi), tolerance = 1e-10)
  expect_warning(eigenvector_centrality(pieces()), "not connected")
  expect_error(eigenvector_centrality(pieces(), convention = "networkx"), "connected")
  expect_error(eigenvector_centrality(pieces(), convention = "ucinet"), "lambda-squared")
  dag <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(suppressWarnings(eigenvector_centrality(dag)), c(0, 0, 1))
  expect_error(eigenvector_centrality(star, max_iter = 0), "max_iter")
})

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
  g <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 3)), nodes = data.frame(id = 1:4))
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

test_that("closeness, harmonic, and eigenvector centrality match igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(closeness(g), unname(igraph::closeness(ig)))
  expect_equal(closeness(g, normalized = TRUE), unname(igraph::closeness(ig, normalized = TRUE)))
  expect_equal(harmonic_centrality(g), unname(igraph::harmonic_centrality(ig)))
  expect_equal(eigenvector_centrality(g), unname(igraph::eigen_centrality(ig)$vector),
               tolerance = 1e-8)

  set.seed(3)
  e <- data.frame(from = sample(25, 80, TRUE), to = sample(25, 80, TRUE))
  e$w <- round(runif(nrow(e), 0.5, 4), 2)
  e <- e[e$from != e$to, ]
  d <- netframe(e, nodes = data.frame(id = 1:25), weight = "w")
  id <- igraph::graph_from_data_frame(e, vertices = d$nodes)
  for (m in c("out", "in", "all")) {
    expect_equal(closeness(d, mode = m, normalized = TRUE),
                 unname(igraph::closeness(id, mode = m, weights = e$w, normalized = TRUE)))
    expect_equal(harmonic_centrality(d, mode = m),
                 unname(igraph::harmonic_centrality(id, mode = m, weights = e$w)))
  }
})

test_that("reach and edge betweenness match igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(edge_betweenness(g), igraph::edge_betweenness(ig))
  expect_equal(reach(g, k = 2), igraph::ego_size(ig, order = 2) - 1)
  expect_equal(distance_weighted_reach(g),
               unname(igraph::harmonic_centrality(ig, normalized = TRUE)))

  e <- random_digraph(30, 120, 3)
  set.seed(3)
  e$w <- round(runif(nrow(e), 0.5, 3), 2)
  nodes <- data.frame(id = 1:30)
  d <- netframe(e, nodes, weight = "w")
  id <- igraph::graph_from_data_frame(e, vertices = nodes)
  expect_equal(edge_betweenness(d), igraph::edge_betweenness(id, weights = e$w),
               tolerance = 1e-9)
  expect_equal(edge_betweenness(d, weighted = FALSE), igraph::edge_betweenness(id, weights = NA),
               tolerance = 1e-9)
  expect_equal(reach(d, k = 2, mode = "in", weighted = FALSE),
               igraph::ego_size(id, order = 2, mode = "in") - 1)
})

test_that("closeness, eigenvector, and information centrality match sna", {
  skip_if_not_installed("sna")
  g <- kite()
  m <- as.matrix(g)
  expect_equal(closeness(g, convention = "sna"),
               unname(sna::closeness(m, gmode = "graph", cmode = "undirected")))
  expect_equal(harmonic_centrality(g, convention = "sna"),
               unname(sna::closeness(m, gmode = "graph", cmode = "suminvundir")))
  ev <- sna::evcent(m, gmode = "graph", use.eigen = TRUE)
  expect_equal(eigenvector_centrality(g, convention = "sna"), abs(unname(ev)), tolerance = 1e-8)
  expect_equal(information_centrality(g), unname(sna::infocent(m, gmode = "graph")),
               tolerance = 1e-10)
  d <- netframe(random_digraph(12, 40, 4), nodes = data.frame(id = 1:12))
  md <- as.matrix(d)
  if (all(components(d)$csize > 1L) && components(d)$no == 1L) {
    expect_equal(information_centrality(d), unname(sna::infocent(md, gmode = "digraph")),
                 tolerance = 1e-10)
  }
})
