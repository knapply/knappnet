# Structural holes, brokerage roles, structural, regular, and automorphic
# equivalence, CONCOR, and blockmodels (R/positions.R). Hand-checked values
# first, then spot checks against igraph and sna. The full comparison over
# every fixture lives in validation/.

# A two-level hierarchy: 1 commands 2 and 3, who each command two others.
tree7 <- function() {
  net_frame(data.frame(source = c(1, 1, 2, 2, 3, 3), target = c(2, 3, 4, 5, 6, 7)),
           nodes = data.frame(id = 1:7))
}

test_that("structural holes match hand-computed values", {
  # Node 1 is tied to 2, 3, and 4; 2 and 3 are also tied.
  g <- net_frame(data.frame(source = c(1, 1, 1, 2), target = c(2, 3, 4, 3)),
                nodes = data.frame(id = 1:4), directed = FALSE)
  # p = 1/3 to each contact; C = 2 (1/3 + 1/6)^2 + (1/3)^2 = 11/18.
  # For node 2: (1/2 + 1/4)^2 + (1/2 + 1/6)^2 = 145/144.
  expect_equal(net_constraint(g), c(11 / 18, 145 / 144, 145 / 144, 1))
  # Borgatti: 3 contacts, 1 tie among them: 3 - 2/3.
  expect_equal(net_effective_size(g), c(7 / 3, 1, 1, 1))
  expect_equal(net_efficiency(g), c(7 / 9, 0.5, 0.5, 1))
})

test_that("structural holes are NaN for isolates and ignore loops", {
  g <- net_frame(data.frame(source = c("a", "c"), target = c("b", "c")),
                nodes = data.frame(id = c("a", "b", "c")), directed = FALSE)
  expect_equal(net_constraint(g)[1:2], c(1, 1))
  expect_true(is.nan(net_constraint(g)[3]))
  expect_true(is.nan(net_effective_size(g)[3]))
  expect_true(is.nan(net_efficiency(g)[3]))
})

test_that("weighted structural holes use weights as tie strengths", {
  e <- data.frame(source = c(1, 1, 2), target = c(2, 3, 3), w = c(1, 3, 2))
  g <- net_frame(e, directed = FALSE, weight = "w", weight_type = "distance")
  expect_false(isTRUE(all.equal(net_constraint(g), net_constraint(g, weighted = FALSE))))
  g2 <- net_frame(transform(e, w = 10 * w), directed = FALSE, weight = "w", weight_type = "distance")
  expect_equal(net_constraint(g), net_constraint(g2))
  # The declared weight type does not change how weights are read.
  g3 <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  expect_equal(net_effective_size(g), net_effective_size(g3))
})

test_that("hierarchy follows Burt's formula", {
  star <- net_frame(data.frame(source = 1, target = 2:4), directed = FALSE)
  h <- net_hierarchy(star)
  expect_equal(h[1], 0)
  expect_true(all(is.nan(h[2:4])))
  # Node 1 tied to 2, 3, 4 with 3 - 4 tied: c = (1/9, 1/4, 1/4).
  g <- net_frame(data.frame(source = c(1, 1, 1, 3), target = c(2, 3, 4, 4)), directed = FALSE)
  cc <- c(1 / 9, 1 / 4, 1 / 4)
  r <- cc / mean(cc)
  expect_equal(net_hierarchy(g)[1], sum(r * log(r)) / (3 * log(3)))
  expect_equal(net_hierarchy(g)[1], 0.0522045454467525, tolerance = 1e-12)
  expect_equal(net_hierarchy(g, threads = 2), net_hierarchy(g))
  expect_true(all(is.nan(net_hierarchy(net_frame(data.frame(source = 1, target = 2), directed = FALSE)))))
})

test_that("constraint matches igraph", {
  skip_if_not_installed("igraph")
  e <- kite_edges()
  g <- net_frame(e, directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = g$nodes)
  expect_equal(net_constraint(g), unname(igraph::constraint(ig)), tolerance = 1e-12)

  e <- random_digraph(30, 120, 7)
  e$w <- round(runif(nrow(e), 0.1, 3), 3)
  nodes <- data.frame(id = 1:30)
  g <- net_frame(e, nodes, weight = "w", weight_type = "distance")
  ig <- igraph::graph_from_data_frame(e, vertices = nodes)
  got <- net_constraint(g)
  want <- unname(igraph::constraint(ig, weights = e$w))
  expect_equal(got[!is.nan(got)], want[!is.nan(got)], tolerance = 1e-12)
})

test_that("brokerage roles match a hand-built example", {
  # a -> b -> c with no a -> c: b brokers (a, c).
  g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")),
                nodes = data.frame(id = c("a", "b", "c"),
                                   team = c("x", "y", "x")))
  b <- net_brokerage(g, "team")
  expect_named(b, c("node", "coordinator", "itinerant", "representative",
                    "gatekeeper", "liaison", "total"))
  expect_equal(b$node, c("a", "b", "c"))
  expect_equal(unlist(b[2, -1]), c(coordinator = 0, itinerant = 1, representative = 0,
                                   gatekeeper = 0, liaison = 0, total = 1))
  expect_equal(net_brokerage(g, c(1, 1, 2))$representative, c(0, 1, 0))
  expect_equal(net_brokerage(g, c(1, 2, 2))$gatekeeper, c(0, 1, 0))
  expect_equal(net_brokerage(g, c(1, 2, 3))$liaison, c(0, 1, 0))
  expect_error(net_brokerage(g, c(1, 2)), "one value per node")
  expect_error(net_brokerage(g, c(1, NA, 2)), "missing")
})

test_that("brokerage matches sna", {
  skip_if_not_installed("sna")
  e <- random_digraph(25, 90, 3)
  g <- net_frame(e, nodes = data.frame(id = 1:25))
  cl <- rep(c("p", "q", "r"), length.out = 25)
  m <- as.matrix(g)
  ours <- as.matrix(net_brokerage(g, cl)[, -1])
  theirs <- sna::brokerage(m, cl)$raw.nli[, c("w_I", "w_O", "b_IO", "b_OI", "b_O", "t")]
  expect_equal(unname(ours), unname(theirs))
})

test_that("structural equivalence matches the definition and sna", {
  e <- random_digraph(12, 40, 11)
  g <- net_frame(e, nodes = data.frame(id = 1:12))
  a <- as.matrix(g)
  d <- net_struct_equiv(g)
  expect_equal(dim(d), c(12L, 12L))
  expect_equal(rownames(d), as.character(1:12))
  # Direct evaluation for one pair, leaving out k = i, j. The two busiest
  # nodes are used so both profiles certainly vary.
  busiest <- order(-(rowSums(a) + colSums(a)))
  i <- busiest[[1]]
  j <- busiest[[2]]
  k <- setdiff(1:12, c(i, j))
  expect_equal(d[i, j], sqrt(sum((a[i, k] - a[j, k])^2 + (a[k, i] - a[k, j])^2)))
  h <- net_struct_equiv(g, "hamming")
  expect_equal(h[i, j], sum(abs(a[i, k] - a[j, k]) + abs(a[k, i] - a[k, j])))
  r <- net_struct_equiv(g, "correlation")
  expect_equal(r[i, j], stats::cor(c(a[i, k], a[k, i]), c(a[j, k], a[k, j])))

  skip_if_not_installed("sna")
  for (method in c("euclidean", "correlation", "hamming")) {
    expect_equal(unname(net_struct_equiv(g, method)),
                 unname(sna::sedist(a, method = method)), tolerance = 1e-12)
  }
})

test_that("concor separates the sides of a complete bipartite graph", {
  e <- expand.grid(source = 1:3, target = 4:6)
  g <- net_frame(e, nodes = data.frame(id = 1:7), directed = FALSE)
  # Node 7 is an isolate and gets a block of its own.
  expect_equal(net_concor(g), c(1L, 1L, 1L, 2L, 2L, 2L, 3L))
  expect_error(net_concor(g, depth = 0), "depth")
})

test_that("concor warns when the partition is numerically fragile", {
  cycle <- net_frame(data.frame(source = 1:7, target = c(2:7, 1)))
  expect_warning(net_concor(cycle), "did not converge")
  # The kite as labelled in validation/fixtures/kite: in the second round a
  # node's correlation decays to about 1e-9 before rounding decides its side.
  e <- data.frame(source = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 6, 6, 7, 8, 9),
                  target   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 7, 8, 8, 9, 10))
  g <- net_frame(e, nodes = data.frame(id = 1:10), directed = FALSE)
  expect_silent(net_concor(g, depth = 1))
  expect_warning(net_concor(g, depth = 2), "rounding")
})

test_that("regular equivalence finds the maximal regular equivalence", {
  g <- tree7()
  r <- net_regular_equiv(g)
  expect_equal(as.vector(r), c(1, 2, 2, 3, 3, 3, 3))
  expect_true(attr(r, "rounds") >= 1L)
  # Undirected networks without a seed are degenerate: tied nodes in one
  # class, isolates in another.
  u <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)), nodes = data.frame(id = 1:4),
                directed = FALSE)
  expect_equal(as.vector(net_regular_equiv(u)), c(1, 1, 1, 2))
  # A seed partition is refined, never merged.
  expect_equal(as.vector(net_regular_equiv(u, seed = c("a", "b", "a", "a"))), c(1, 2, 1, 3))
  u$nodes$part <- c("a", "b", "a", "a")
  expect_equal(net_regular_equiv(u, seed = "part"), net_regular_equiv(u, seed = u$nodes$part))
  # A mutual tie differs from a one-way tie: 2 has a mutual tie, 3 does not.
  d <- net_frame(data.frame(source = c(1, 2, 1), target = c(2, 1, 3)))
  expect_equal(as.vector(net_regular_equiv(d)), c(1, 2, 3))
  # Weights are categories.
  w <- net_frame(data.frame(source = c(1, 1), target = c(2, 3), w = c(1, 2)), weight = "w", weight_type = "distance")
  expect_equal(as.vector(net_regular_equiv(w)), c(1, 2, 3))
  expect_equal(as.vector(net_regular_equiv(w, weighted = FALSE)), c(1, 2, 2))
  expect_error(net_regular_equiv(g, seed = 1:3), "one value per node")
  expect_error(net_regular_equiv(g, seed = c(1, NA, 1, 1, 1, 1, 1)), "missing")
})

test_that("automorphic equivalence gives the orbits", {
  g <- tree7()
  expect_equal(net_auto_equiv(g), c(1L, 2L, 2L, 3L, 3L, 3L, 3L))
  # Two triangles joined through node 4.
  h <- net_frame(data.frame(source = c(1, 1, 2, 3, 4, 5, 5, 6), target = c(2, 3, 3, 4, 5, 6, 7, 7)),
                directed = FALSE)
  expect_equal(net_auto_equiv(h), c(1L, 1L, 2L, 3L, 2L, 1L, 1L))
  # Direction matters: reversing one arc breaks the symmetry of a path.
  p <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)))
  expect_equal(net_auto_equiv(p), 1:3)
  pu <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)), directed = FALSE)
  expect_equal(net_auto_equiv(pu), c(1L, 2L, 1L))
  # A self-loop is a property of its node; a parallel edge is not.
  l <- net_frame(data.frame(source = c(1, 2, 2, 3), target = c(2, 3, 3, 3)), directed = FALSE)
  expect_equal(net_auto_equiv(l), 1:3)
  # Krackhardt's kite has no non-trivial automorphism.
  kite <- net_frame(data.frame(
    source = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    target   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)), directed = FALSE)
  expect_equal(net_auto_equiv(kite), 1:10)
  # Automorphic equivalence refines regular equivalence.
  expect_true(all(tapply(as.vector(net_regular_equiv(g)), net_auto_equiv(g),
                         function(v) length(unique(v))) == 1L))
  expect_error(net_auto_equiv(g, max_steps = 0), "`max_steps`")
})

test_that("blockmodel densities, image, and fit match hand values", {
  # Blocks {1, 2} and {3, 4}: 1-2 tied, 3-4 not, 1-3 and 2-3 tied.
  g <- net_frame(data.frame(source = c(1, 1, 2), target = c(2, 3, 3)), nodes = data.frame(id = 1:4),
                directed = FALSE)
  bm <- net_blockmodel(g, c("a", "a", "b", "b"))
  expect_s3_class(bm, "blockmodel")
  expect_equal(bm$blocks, c("a", "b"))
  expect_equal(bm$sizes, c(2L, 2L))
  expect_equal(unname(bm$density), matrix(c(1, 0.5, 0.5, 0), 2))
  expect_equal(dimnames(bm$density), list(c("a", "b"), c("a", "b")))
  expect_equal(unname(bm$ties), matrix(c(2, 2, 2, 0), 2))
  expect_equal(unname(bm$pairs), matrix(c(2, 4, 4, 2), 2))
  expect_equal(bm$overall_density, 0.5)
  expect_equal(bm$cutoff, 0.5)
  expect_equal(unname(bm$image), matrix(c(1, 1, 1, 0), 2))
  expect_equal(unname(net_blockmodel(g, c("a", "a", "b", "b"), cutoff = 0.75)$image),
               matrix(c(1, 0, 0, 0), 2))
  # R-squared: ties (1,1, 1,0,1,0, 1,0,1,0, 0,0) against densities.
  a <- as.matrix(g)
  b <- c(1, 1, 2, 2)
  ties <- a[row(a) != col(a)]
  fitted <- bm$density[cbind(b[row(a)[row(a) != col(a)]], b[col(a)[row(a) != col(a)]])]
  expect_equal(bm$r_squared, cor(ties, fitted)^2)
  s <- net_block_summary(bm)
  expect_equal(nrow(s), 4L)
  expect_equal(s$density, c(1, 0.5, 0.5, 0))
  expect_equal(s$from, c("a", "b", "a", "b"))
  expect_equal(s$size_to, c(2L, 2L, 2L, 2L))
  expect_output(print(bm), "2 blocks")
  # Diagonal with loops; a one-node block has no diagonal pairs.
  one <- net_blockmodel(g, c("a", "a", "b", "c"))
  expect_true(is.nan(one$density["b", "b"]))
  expect_true(is.nan(one$image["b", "b"]))
  expect_equal(net_blockmodel(g, c("a", "a", "b", "c"), loops = TRUE)$density["b", "b"], 0)
  # Column name and membership vectors with attributes.
  g$nodes$part <- c("a", "a", "b", "b")
  expect_equal(net_blockmodel(g, "part")$density, bm$density)
  expect_equal(net_blockmodel(g, net_regular_equiv(g))$blocks, c(1L, 2L))
  # Weights are mean tie values; direction gives an asymmetric matrix.
  w <- net_frame(data.frame(source = c(1, 1, 3), target = c(2, 3, 1), w = c(2, 4, 6)),
                nodes = data.frame(id = 1:3), weight = "w", weight_type = "distance")
  bw <- net_blockmodel(w, c(1, 1, 2))
  expect_equal(unname(bw$density), matrix(c(1, 3, 2, NaN), 2))
  expect_equal(unname(net_blockmodel(w, c(1, 1, 2), weighted = FALSE)$density),
               matrix(c(0.5, 0.5, 0.5, NaN), 2))
  expect_error(net_blockmodel(g, 1:3), "one value per node")
  expect_error(net_blockmodel(g, "part", cutoff = "x"), "`cutoff`")
  expect_error(net_block_summary(list()), "blockmodel")
})

test_that("results match sna", {
  skip_if_not_installed("sna")
  g <- tree7()
  m <- as.matrix(g)
  # Zero distance in redist(method = "catrege") marks equivalent nodes.
  d <- as.matrix(sna::redist(m, mode = "digraph", method = "catrege"))
  classes <- match(apply(d == 0, 1, paste, collapse = ""),
                   unique(apply(d == 0, 1, paste, collapse = "")))
  expect_equal(as.vector(net_regular_equiv(g)), classes)
  set.seed(5)
  e <- data.frame(source = sample(16, 60, TRUE), target = sample(16, 60, TRUE))
  e <- e[e$source != e$target, ]
  e <- e[!duplicated(e), ]
  h <- net_frame(e, nodes = data.frame(id = 1:16))
  mh <- as.matrix(h)
  seed <- rep(1:2, 8)
  dh <- as.matrix(sna::redist(mh, mode = "digraph", method = "catrege",
                              seed.partition = seed))
  keys <- apply(dh == 0, 1, paste, collapse = "")
  expect_equal(as.vector(net_regular_equiv(h, seed = seed)), match(keys, unique(keys)))
  bm <- net_blockmodel(h, seed)
  sb <- sna::blockmodel(mh, seed, block.content = "density")$block.model
  expect_equal(unname(bm$density), unname(sb))
})

test_that("plot.blockmodel() draws densities and the permuted matrix", {
  g <- net_example("florentine")
  bm <- net_blockmodel(g, net_concor(g, depth = 2))
  on_null_device({
    expect_silent(out <- plot(bm))
    expect_identical(out, bm)
    expect_silent(plot(bm, network = g))
    expect_silent(plot(bm, network = g, which = "matrix"))
    expect_silent(plot(bm, network = g, which = "density", values = FALSE, main = "Densities"))
    expect_restores_par(plot(bm, network = g))
    expect_error(plot(bm, which = "matrix"), "needs `network`")
    other <- net_frame(data.frame(source = 1, target = 2))
    expect_error(plot(bm, network = other), "one node per element")
    # A one-node block without loops has a NaN density and still draws.
    one <- net_blockmodel(g, c("a", rep("b", nrow(g$nodes) - 1L)))
    expect_true(is.nan(one$density["a", "a"]))
    expect_silent(plot(one, network = g))
  })
})
