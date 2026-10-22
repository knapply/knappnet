# Density, components, fragmentation, closeness, harmonic and eigenvector
# centrality, centralization, and the convention switch. Hand-checked values
# first, then spot checks against igraph, sna, and (through published
# numbers) UCINET. The full comparison lives in validation/.

kite <- function() {
  e <- data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
  netframe(e, nodes = data.frame(id = 1:10), directed = FALSE)
}

# A path 1 - 2 - 3, an edge 4 - 5, and an isolate 6.
pieces <- function() {
  netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)),
           nodes = data.frame(id = 1:6), directed = FALSE)
}

test_that("density counts edges over possible edges", {
  g <- netframe(data.frame(from = c(1, 2, 3), to = c(2, 3, 4)), directed = FALSE)
  expect_equal(graph_density(g), 0.5)
  expect_equal(graph_density(g, loops = TRUE), 0.3)
  multi <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 2, 3, 3)))
  expect_equal(graph_density(multi), 4 / 6)
  expect_equal(graph_density(multi, convention = "sna"), 2 / 6)
  expect_error(graph_density(multi, loops = TRUE, convention = "sna"), "self-loops")
  one <- netframe(data.frame(from = 1, to = 1))
  expect_true(is.nan(graph_density(one)))
  expect_equal(graph_density(one, convention = "networkx"), 0)
})

test_that("components and fragmentation", {
  g <- pieces()
  cc <- components(g)
  expect_equal(cc$membership, c(1L, 1L, 1L, 2L, 2L, 3L))
  expect_equal(cc$csize, c(3L, 2L, 1L))
  expect_equal(cc$no, 3L)
  # Ordered pairs in the same component: 6 + 2 = 8 of 30.
  expect_equal(fragmentation(g), 22 / 30)
  d <- netframe(data.frame(from = c(1, 2, 3, 3), to = c(2, 3, 1, 4)))
  expect_equal(components(d, "strong")$membership, c(1L, 1L, 1L, 2L))
  expect_equal(components(d, "weak")$no, 1L)
  expect_equal(fragmentation(d, "strong"), 6 / 12)
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

test_that("centralization of a star is 1 under every convention", {
  star <- netframe(data.frame(from = 1, to = 2:6), directed = FALSE)
  for (conv in c("igraph", "sna", "ucinet")) {
    for (m in c("degree", "betweenness", "closeness")) {
      expect_equal(centralization(star, m, loops = FALSE, convention = conv), 1,
                   tolerance = 1e-12, info = paste(conv, m))
    }
  }
  expect_error(centralization(star, convention = "networkx"), "NetworkX")
  ring <- netframe(data.frame(from = 1:5, to = c(2:5, 1)), directed = FALSE)
  expect_equal(centralization(ring), 0)
})

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

test_that("results match igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(graph_density(g), igraph::edge_density(ig))
  expect_equal(components(g)$membership, unname(igraph::components(ig)$membership))
  expect_equal(closeness(g), unname(igraph::closeness(ig)))
  expect_equal(closeness(g, normalized = TRUE), unname(igraph::closeness(ig, normalized = TRUE)))
  expect_equal(harmonic_centrality(g), unname(igraph::harmonic_centrality(ig)))
  expect_equal(eigenvector_centrality(g), unname(igraph::eigen_centrality(ig)$vector),
               tolerance = 1e-8)
  expect_equal(centralization(g), igraph::centr_degree(ig, loops = TRUE)$centralization)
  expect_equal(centralization(g, "betweenness"), igraph::centr_betw(ig)$centralization)
  expect_equal(centralization(g, "closeness"), igraph::centr_clo(ig)$centralization)

  set.seed(3)
  e <- data.frame(from = sample(25, 80, TRUE), to = sample(25, 80, TRUE))
  e$w <- round(runif(nrow(e), 0.5, 4), 2)
  e <- e[e$from != e$to, ]
  d <- netframe(e, nodes = data.frame(id = 1:25), weight = "w")
  id <- igraph::graph_from_data_frame(e, vertices = d$nodes)
  expect_equal(components(d, "strong")$no, igraph::components(id, "strong")$no)
  for (m in c("out", "in", "all")) {
    expect_equal(closeness(d, mode = m, normalized = TRUE),
                 unname(igraph::closeness(id, mode = m, weights = e$w, normalized = TRUE)))
    expect_equal(harmonic_centrality(d, mode = m),
                 unname(igraph::harmonic_centrality(id, mode = m, weights = e$w)))
  }
  expect_equal(centralization(d, mode = "in"),
               igraph::centr_degree(id, mode = "in", loops = TRUE)$centralization)
})

test_that("results match sna", {
  skip_if_not_installed("sna")
  g <- kite()
  m <- as.matrix(g)
  expect_equal(graph_density(g, convention = "sna"), unname(sna::gden(m, mode = "graph")))
  expect_equal(closeness(g, convention = "sna"),
               unname(sna::closeness(m, gmode = "graph", cmode = "undirected")))
  expect_equal(harmonic_centrality(g, convention = "sna"),
               unname(sna::closeness(m, gmode = "graph", cmode = "suminvundir")))
  ev <- sna::evcent(m, gmode = "graph", use.eigen = TRUE)
  expect_equal(eigenvector_centrality(g, convention = "sna"), abs(unname(ev)), tolerance = 1e-8)
  expect_equal(centralization(g, "degree", convention = "sna"),
               sna::centralization(m, sna::degree, mode = "graph"))
  expect_equal(centralization(g, "closeness", convention = "sna"),
               sna::centralization(m, sna::closeness, mode = "graph"))
})
