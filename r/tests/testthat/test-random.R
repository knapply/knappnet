# Random networks: exact constraints, reproducibility by seed, and the
# netframe shape of the results. The generators' distributions are checked
# in validation/cug_check.py.

simple_undirected <- function(g) {
  e <- g$edges
  expect_false(g$directed)
  expect_true(all(e$from < e$to))
  expect_false(anyDuplicated(paste(e$from, e$to)) > 0)
}

test_that("G(n, p) and G(n, m)", {
  g <- random_gnp(20, 0.3, seed = 1, threads = 2)
  expect_s3_class(g, "netframe")
  expect_equal(g$nodes$id, 1:20)
  simple_undirected(g)
  expect_equal(random_gnp(20, 0.3, seed = 1, threads = 1)$edges, g$edges)
  expect_equal(nrow(random_gnp(20, 1, seed = 1)$edges), 190)
  expect_equal(nrow(random_gnp(20, 0, seed = 1)$edges), 0)
  expect_equal(nrow(random_gnp(5, 1, directed = TRUE, loops = TRUE, seed = 1)$edges), 25)

  h <- random_gnm(20, 37, seed = 2, threads = 2)
  expect_equal(nrow(h$edges), 37)
  simple_undirected(h)
  d <- random_gnm(10, 40, directed = TRUE, seed = 2)
  expect_true(d$directed)
  expect_equal(nrow(d$edges), 40)
  expect_true(all(d$edges$from != d$edges$to))
  expect_error(random_gnm(5, 11, seed = 1), "edges")
  expect_error(random_gnp(5, 1.5, seed = 1), "probability")

  # Replicates: a list, the first of which is the single draw.
  reps <- random_gnm(10, 15, reps = 3, seed = 7, threads = 2)
  expect_length(reps, 3)
  expect_equal(reps[[1]]$edges, random_gnm(10, 15, seed = 7)$edges)
  expect_false(identical(reps[[1]]$edges, reps[[2]]$edges))
  expect_equal(random_gnm(10, 15, reps = 3, seed = 7, threads = 1)[[3]]$edges, reps[[3]]$edges)
  # NULL seeds come from R's generator.
  set.seed(1)
  a <- random_gnm(10, 15)
  set.seed(1)
  expect_equal(random_gnm(10, 15)$edges, a$edges)
})

test_that("configuration model keeps the degrees", {
  karate <- example_network("karate")
  d <- degree(karate)
  multi <- random_configuration(d, seed = 1, threads = 2)
  expect_equal(degree(multi), d)
  # Rejection rarely finds a simple graph for karate's two hubs, so the
  # simple variant is tested on a sequence without them.
  d2 <- degree(random_gnm(30, 45, seed = 2))
  simple <- random_configuration(d2, simple = TRUE, seed = 1, threads = 2)
  expect_equal(degree(simple), d2)
  simple_undirected(simple)
  dir <- random_configuration(c(2, 1, 0, 1), in_degrees = c(0, 1, 2, 1), seed = 1)
  expect_true(dir$directed)
  expect_equal(degree(dir, "out"), c(2, 1, 0, 1))
  expect_equal(degree(dir, "in"), c(0, 1, 2, 1))
  expect_error(random_configuration(c(1, 1, 1), seed = 1), "odd|sum")
  # Two nodes of degree 2 can only be joined by a double edge.
  expect_error(random_configuration(c(2, 2), simple = TRUE, max_tries = 5, seed = 1), "simple")
  expect_error(random_configuration(c(-1, 1), seed = 1), "non-negative")
})

test_that("U|MAN keeps the dyad census", {
  g <- random_uman(8, mutual = 5, asym = 10, seed = 1, threads = 2)
  expect_true(g$directed)
  expect_equal(nrow(g$edges), 2 * 5 + 10)
  m <- as.matrix(g)
  expect_equal(sum(m * t(m)) / 2, 5)
  expect_equal(sum(m) - 2 * 5, 10)
  expect_equal(random_uman(8, 5, 10, seed = 1, threads = 1)$edges, g$edges)
  expect_error(random_uman(4, 4, 3, seed = 1), "dyads")
})

test_that("small worlds and preferential attachment", {
  sw <- random_small_world(30, 4, 0.1, seed = 1, threads = 2)
  simple_undirected(sw)
  expect_equal(nrow(sw$edges), 60)
  lattice <- random_small_world(30, 4, 0, seed = 1)
  expect_equal(degree(lattice), rep(4, 30))
  expect_error(random_small_world(10, 3, 0.1, seed = 1), "even|k")
  ba <- random_preferential(30, 2, seed = 1, threads = 2)
  simple_undirected(ba)
  expect_equal(nrow(ba$edges), 2 * 28)
  expect_true(all(degree(ba) >= 1))
  expect_error(random_preferential(5, 5, seed = 1), "m")
})

test_that("stochastic block models", {
  blocks <- random_sbm(c(5, 5), matrix(c(1, 0, 0, 1), 2), seed = 1, threads = 2)
  expect_equal(nrow(blocks$edges), 2 * choose(5, 2))
  expect_true(all((blocks$edges$from <= 5) == (blocks$edges$to <= 5)))
  expect_equal(components(blocks)$no, 2)
  d <- random_sbm(c(3, 3), matrix(c(0, 1, 0, 0), 2, byrow = TRUE), directed = TRUE, seed = 1)
  # Block 1 -> block 2 only.
  expect_equal(nrow(d$edges), 9)
  expect_true(all(d$edges$from <= 3 & d$edges$to > 3))
  expect_error(random_sbm(c(2, 2), matrix(c(0.5, 0.2, 0.3, 0.5), 2), seed = 1), "symmetric")
  expect_error(random_sbm(c(2, 2), matrix(2, 2, 2), seed = 1), "probabilities")
})
