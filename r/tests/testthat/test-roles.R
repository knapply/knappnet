# Regular and automorphic equivalence and blockmodels. Hand-checked values
# first, then spot checks against sna. The full comparison over every
# fixture lives in validation/.

# A two-level hierarchy: 1 commands 2 and 3, who each command two others.
tree7 <- function() {
  netframe(data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7)),
           nodes = data.frame(id = 1:7))
}

test_that("regular equivalence finds the maximal regular equivalence", {
  g <- tree7()
  r <- regular_equivalence(g)
  expect_equal(as.vector(r), c(1, 2, 2, 3, 3, 3, 3))
  expect_true(attr(r, "rounds") >= 1L)
  # Undirected networks without a seed are degenerate: tied nodes in one
  # class, isolates in another.
  u <- netframe(data.frame(from = c(1, 2), to = c(2, 3)), nodes = data.frame(id = 1:4),
                directed = FALSE)
  expect_equal(as.vector(regular_equivalence(u)), c(1, 1, 1, 2))
  # A seed partition is refined, never merged.
  expect_equal(as.vector(regular_equivalence(u, seed = c("a", "b", "a", "a"))), c(1, 2, 1, 3))
  u$nodes$part <- c("a", "b", "a", "a")
  expect_equal(regular_equivalence(u, seed = "part"), regular_equivalence(u, seed = u$nodes$part))
  # A mutual tie differs from a one-way tie: 2 has a mutual tie, 3 does not.
  d <- netframe(data.frame(from = c(1, 2, 1), to = c(2, 1, 3)))
  expect_equal(as.vector(regular_equivalence(d)), c(1, 2, 3))
  # Weights are categories.
  w <- netframe(data.frame(from = c(1, 1), to = c(2, 3), w = c(1, 2)), weight = "w")
  expect_equal(as.vector(regular_equivalence(w)), c(1, 2, 3))
  expect_equal(as.vector(regular_equivalence(w, weighted = FALSE)), c(1, 2, 2))
  expect_error(regular_equivalence(g, seed = 1:3), "one value per node")
  expect_error(regular_equivalence(g, seed = c(1, NA, 1, 1, 1, 1, 1)), "missing")
})

test_that("automorphic equivalence gives the orbits", {
  g <- tree7()
  expect_equal(automorphic_equivalence(g), c(1L, 2L, 2L, 3L, 3L, 3L, 3L))
  # Two triangles joined through node 4.
  h <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 5, 5, 6), to = c(2, 3, 3, 4, 5, 6, 7, 7)),
                directed = FALSE)
  expect_equal(automorphic_equivalence(h), c(1L, 1L, 2L, 3L, 2L, 1L, 1L))
  # Direction matters: reversing one arc breaks the symmetry of a path.
  p <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(automorphic_equivalence(p), 1:3)
  pu <- netframe(data.frame(from = c(1, 2), to = c(2, 3)), directed = FALSE)
  expect_equal(automorphic_equivalence(pu), c(1L, 2L, 1L))
  # A self-loop is a property of its node; a parallel edge is not.
  l <- netframe(data.frame(from = c(1, 2, 2, 3), to = c(2, 3, 3, 3)), directed = FALSE)
  expect_equal(automorphic_equivalence(l), 1:3)
  # Krackhardt's kite has no non-trivial automorphism.
  kite <- netframe(data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)), directed = FALSE)
  expect_equal(automorphic_equivalence(kite), 1:10)
  # Automorphic equivalence refines regular equivalence.
  expect_true(all(tapply(as.vector(regular_equivalence(g)), automorphic_equivalence(g),
                         function(v) length(unique(v))) == 1L))
  expect_error(automorphic_equivalence(g, max_steps = 0), "`max_steps`")
})

test_that("blockmodel densities, image, and fit match hand values", {
  # Blocks {1, 2} and {3, 4}: 1-2 tied, 3-4 not, 1-3 and 2-3 tied.
  g <- netframe(data.frame(from = c(1, 1, 2), to = c(2, 3, 3)), nodes = data.frame(id = 1:4),
                directed = FALSE)
  bm <- blockmodel(g, c("a", "a", "b", "b"))
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
  expect_equal(unname(blockmodel(g, c("a", "a", "b", "b"), cutoff = 0.75)$image),
               matrix(c(1, 0, 0, 0), 2))
  # R-squared: ties (1,1, 1,0,1,0, 1,0,1,0, 0,0) against densities.
  a <- as.matrix(g)
  b <- c(1, 1, 2, 2)
  ties <- a[row(a) != col(a)]
  fitted <- bm$density[cbind(b[row(a)[row(a) != col(a)]], b[col(a)[row(a) != col(a)]])]
  expect_equal(bm$r_squared, cor(ties, fitted)^2)
  s <- block_summary(bm)
  expect_equal(nrow(s), 4L)
  expect_equal(s$density, c(1, 0.5, 0.5, 0))
  expect_equal(s$from, c("a", "b", "a", "b"))
  expect_equal(s$size_to, c(2L, 2L, 2L, 2L))
  expect_output(print(bm), "2 blocks")
  # Diagonal with loops; a one-node block has no diagonal pairs.
  one <- blockmodel(g, c("a", "a", "b", "c"))
  expect_true(is.nan(one$density["b", "b"]))
  expect_true(is.nan(one$image["b", "b"]))
  expect_equal(blockmodel(g, c("a", "a", "b", "c"), loops = TRUE)$density["b", "b"], 0)
  # Column name and membership vectors with attributes.
  g$nodes$part <- c("a", "a", "b", "b")
  expect_equal(blockmodel(g, "part")$density, bm$density)
  expect_equal(blockmodel(g, regular_equivalence(g))$blocks, c(1L, 2L))
  # Weights are mean tie values; direction gives an asymmetric matrix.
  w <- netframe(data.frame(from = c(1, 1, 3), to = c(2, 3, 1), w = c(2, 4, 6)),
                nodes = data.frame(id = 1:3), weight = "w")
  bw <- blockmodel(w, c(1, 1, 2))
  expect_equal(unname(bw$density), matrix(c(1, 3, 2, NaN), 2))
  expect_equal(unname(blockmodel(w, c(1, 1, 2), weighted = FALSE)$density),
               matrix(c(0.5, 0.5, 0.5, NaN), 2))
  expect_error(blockmodel(g, 1:3), "one value per node")
  expect_error(blockmodel(g, "part", cutoff = "x"), "`cutoff`")
  expect_error(block_summary(list()), "blockmodel")
})

test_that("results match sna", {
  skip_if_not_installed("sna")
  g <- tree7()
  m <- as.matrix(g)
  # Zero distance in redist(method = "catrege") marks equivalent nodes.
  d <- as.matrix(sna::redist(m, mode = "digraph", method = "catrege"))
  classes <- match(apply(d == 0, 1, paste, collapse = ""),
                   unique(apply(d == 0, 1, paste, collapse = "")))
  expect_equal(as.vector(regular_equivalence(g)), classes)
  set.seed(5)
  e <- data.frame(from = sample(16, 60, TRUE), to = sample(16, 60, TRUE))
  e <- e[e$from != e$to, ]
  e <- e[!duplicated(e), ]
  h <- netframe(e, nodes = data.frame(id = 1:16))
  mh <- as.matrix(h)
  seed <- rep(1:2, 8)
  dh <- as.matrix(sna::redist(mh, mode = "digraph", method = "catrege",
                              seed.partition = seed))
  keys <- apply(dh == 0, 1, paste, collapse = "")
  expect_equal(as.vector(regular_equivalence(h, seed = seed)), match(keys, unique(keys)))
  bm <- blockmodel(h, seed)
  sb <- sna::blockmodel(mh, seed, block.content = "density")$block.model
  expect_equal(unname(bm$density), unname(sb))
})
