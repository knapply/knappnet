# Brokerage and positions: hand-checked values, then spot checks against
# igraph and sna. The full comparison over every fixture lives in
# validation/ at the repository root.

kite_edges <- function() {
  # Krackhardt's kite as an undirected edge list.
  data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
}

random_digraph <- function(n, m, seed) {
  set.seed(seed)
  e <- data.frame(from = sample(n, m, TRUE), to = sample(n, m, TRUE))
  e <- e[e$from != e$to, ]
  e[!duplicated(e), ]
}

test_that("structural holes match hand-computed values", {
  # Node 1 is tied to 2, 3, and 4; 2 and 3 are also tied.
  g <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 3)),
                nodes = data.frame(id = 1:4), directed = FALSE)
  # p = 1/3 to each contact; C = 2 (1/3 + 1/6)^2 + (1/3)^2 = 11/18.
  # For node 2: (1/2 + 1/4)^2 + (1/2 + 1/6)^2 = 145/144.
  expect_equal(constraint(g), c(11 / 18, 145 / 144, 145 / 144, 1))
  # Borgatti: 3 contacts, 1 tie among them: 3 - 2/3.
  expect_equal(effective_size(g), c(7 / 3, 1, 1, 1))
  expect_equal(efficiency(g), c(7 / 9, 0.5, 0.5, 1))
})

test_that("structural holes are NaN for isolates and ignore loops", {
  g <- netframe(data.frame(from = c("a", "c"), to = c("b", "c")),
                nodes = data.frame(id = c("a", "b", "c")), directed = FALSE)
  expect_equal(constraint(g)[1:2], c(1, 1))
  expect_true(is.nan(constraint(g)[3]))
  expect_true(is.nan(effective_size(g)[3]))
  expect_true(is.nan(efficiency(g)[3]))
})

test_that("weighted structural holes use weights as tie strengths", {
  e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3), w = c(1, 3, 2))
  g <- netframe(e, directed = FALSE, weight = "w")
  expect_false(isTRUE(all.equal(constraint(g), constraint(g, weighted = FALSE))))
  g2 <- netframe(transform(e, w = 10 * w), directed = FALSE, weight = "w")
  expect_equal(constraint(g), constraint(g2))
  # The declared weight type does not change how weights are read.
  g3 <- netframe(e, directed = FALSE, weight = "w", weight_type = "strength")
  expect_equal(effective_size(g), effective_size(g3))
})

test_that("constraint matches igraph", {
  skip_if_not_installed("igraph")
  e <- kite_edges()
  g <- netframe(e, directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = g$nodes)
  expect_equal(constraint(g), unname(igraph::constraint(ig)), tolerance = 1e-12)

  e <- random_digraph(30, 120, 7)
  e$w <- round(runif(nrow(e), 0.1, 3), 3)
  nodes <- data.frame(id = 1:30)
  g <- netframe(e, nodes, weight = "w")
  ig <- igraph::graph_from_data_frame(e, vertices = nodes)
  got <- constraint(g)
  want <- unname(igraph::constraint(ig, weights = e$w))
  expect_equal(got[!is.nan(got)], want[!is.nan(got)], tolerance = 1e-12)
})

test_that("brokerage roles match a hand-built example", {
  # a -> b -> c with no a -> c: b brokers (a, c).
  g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")),
                nodes = data.frame(id = c("a", "b", "c"),
                                   team = c("x", "y", "x")))
  b <- brokerage(g, "team")
  expect_named(b, c("node", "coordinator", "itinerant", "representative",
                    "gatekeeper", "liaison", "total"))
  expect_equal(b$node, c("a", "b", "c"))
  expect_equal(unlist(b[2, -1]), c(coordinator = 0, itinerant = 1, representative = 0,
                                   gatekeeper = 0, liaison = 0, total = 1))
  expect_equal(brokerage(g, c(1, 1, 2))$representative, c(0, 1, 0))
  expect_equal(brokerage(g, c(1, 2, 2))$gatekeeper, c(0, 1, 0))
  expect_equal(brokerage(g, c(1, 2, 3))$liaison, c(0, 1, 0))
  expect_error(brokerage(g, c(1, 2)), "one value per node")
  expect_error(brokerage(g, c(1, NA, 2)), "missing")
})

test_that("brokerage matches sna", {
  skip_if_not_installed("sna")
  e <- random_digraph(25, 90, 3)
  g <- netframe(e, nodes = data.frame(id = 1:25))
  cl <- rep(c("p", "q", "r"), length.out = 25)
  m <- as.matrix(g)
  ours <- as.matrix(brokerage(g, cl)[, -1])
  theirs <- sna::brokerage(m, cl)$raw.nli[, c("w_I", "w_O", "b_IO", "b_OI", "b_O", "t")]
  expect_equal(unname(ours), unname(theirs))
})

test_that("structural equivalence matches the definition and sna", {
  e <- random_digraph(12, 40, 11)
  g <- netframe(e, nodes = data.frame(id = 1:12))
  a <- as.matrix(g)
  d <- structural_equivalence(g)
  expect_equal(dim(d), c(12L, 12L))
  expect_equal(rownames(d), as.character(1:12))
  # Direct evaluation for one pair, leaving out k = i, j. The two busiest
  # nodes are used so both profiles certainly vary.
  busiest <- order(-(rowSums(a) + colSums(a)))
  i <- busiest[[1]]
  j <- busiest[[2]]
  k <- setdiff(1:12, c(i, j))
  expect_equal(d[i, j], sqrt(sum((a[i, k] - a[j, k])^2 + (a[k, i] - a[k, j])^2)))
  h <- structural_equivalence(g, "hamming")
  expect_equal(h[i, j], sum(abs(a[i, k] - a[j, k]) + abs(a[k, i] - a[k, j])))
  r <- structural_equivalence(g, "correlation")
  expect_equal(r[i, j], stats::cor(c(a[i, k], a[k, i]), c(a[j, k], a[k, j])))

  skip_if_not_installed("sna")
  for (method in c("euclidean", "correlation", "hamming")) {
    expect_equal(unname(structural_equivalence(g, method)),
                 unname(sna::sedist(a, method = method)), tolerance = 1e-12)
  }
})

test_that("concor separates the sides of a complete bipartite graph", {
  e <- expand.grid(from = 1:3, to = 4:6)
  g <- netframe(e, nodes = data.frame(id = 1:7), directed = FALSE)
  # Node 7 is an isolate and gets a block of its own.
  expect_equal(concor(g), c(1L, 1L, 1L, 2L, 2L, 2L, 3L))
  expect_error(concor(g, depth = 0), "depth")
})

test_that("concor warns when the partition is numerically fragile", {
  cycle <- netframe(data.frame(from = 1:7, to = c(2:7, 1)))
  expect_warning(concor(cycle), "did not converge")
  # The kite as labelled in validation/fixtures/kite: in the second round a
  # node's correlation decays to about 1e-9 before rounding decides its side.
  e <- data.frame(from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 6, 6, 7, 8, 9),
                  to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 7, 8, 8, 9, 10))
  g <- netframe(e, nodes = data.frame(id = 1:10), directed = FALSE)
  expect_silent(concor(g, depth = 1))
  expect_warning(concor(g, depth = 2), "rounding")
})
