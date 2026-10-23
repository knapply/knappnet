# Random networks: exact constraints, reproducibility by seed, and the
# netframe shape of the results. The generators' distributions are checked
# in validation/cug_check.py.

simple_undirected <- function(g) {
  e <- g$edges
  expect_false(g$directed)
  expect_true(all(e$source < e$target))
  expect_false(anyDuplicated(paste(e$source, e$target)) > 0)
}

test_that("G(n, p) and G(n, m)", {
  g <- net_random_gnp(20, 0.3, seed = 1, threads = 2)
  expect_s3_class(g, "netframe")
  expect_equal(g$nodes$id, 1:20)
  simple_undirected(g)
  expect_equal(net_random_gnp(20, 0.3, seed = 1, threads = 1)$edges, g$edges)
  expect_equal(nrow(net_random_gnp(20, 1, seed = 1)$edges), 190)
  expect_equal(nrow(net_random_gnp(20, 0, seed = 1)$edges), 0)
  expect_equal(nrow(net_random_gnp(5, 1, directed = TRUE, loops = TRUE, seed = 1)$edges), 25)

  h <- net_random_gnm(20, 37, seed = 2, threads = 2)
  expect_equal(nrow(h$edges), 37)
  simple_undirected(h)
  d <- net_random_gnm(10, 40, directed = TRUE, seed = 2)
  expect_true(d$directed)
  expect_equal(nrow(d$edges), 40)
  expect_true(all(d$edges$source != d$edges$target))
  expect_error(net_random_gnm(5, 11, seed = 1), "edges")
  expect_error(net_random_gnp(5, 1.5, seed = 1), "probability")

  # Replicates: a list, the first of which is the single draw.
  reps <- net_random_gnm(10, 15, reps = 3, seed = 7, threads = 2)
  expect_length(reps, 3)
  expect_equal(reps[[1]]$edges, net_random_gnm(10, 15, seed = 7)$edges)
  expect_false(identical(reps[[1]]$edges, reps[[2]]$edges))
  expect_equal(net_random_gnm(10, 15, reps = 3, seed = 7, threads = 1)[[3]]$edges, reps[[3]]$edges)
  # NULL seeds come from R's generator.
  set.seed(1)
  a <- net_random_gnm(10, 15)
  set.seed(1)
  expect_equal(net_random_gnm(10, 15)$edges, a$edges)
})

test_that("configuration model keeps the degrees", {
  karate <- net_example("karate")
  d <- net_degree(karate)
  multi <- net_random_config(d, seed = 1, threads = 2)
  expect_equal(net_degree(multi), d)
  # Rejection rarely finds a simple graph for karate's two hubs, so the
  # simple variant is tested on a sequence without them.
  d2 <- net_degree(net_random_gnm(30, 45, seed = 2))
  simple <- net_random_config(d2, simple = TRUE, seed = 1, threads = 2)
  expect_equal(net_degree(simple), d2)
  simple_undirected(simple)
  dir <- net_random_config(c(2, 1, 0, 1), in_degrees = c(0, 1, 2, 1), seed = 1)
  expect_true(dir$directed)
  expect_equal(net_degree(dir, "out"), c(2, 1, 0, 1))
  expect_equal(net_degree(dir, "in"), c(0, 1, 2, 1))
  expect_error(net_random_config(c(1, 1, 1), seed = 1), "odd|sum")
  # Two nodes of degree 2 can only be joined by a double edge.
  expect_error(net_random_config(c(2, 2), simple = TRUE, max_tries = 5, seed = 1), "simple")
  expect_error(net_random_config(c(-1, 1), seed = 1), "non-negative")
})

test_that("U|MAN keeps the dyad census", {
  g <- net_random_uman(8, mutual = 5, asym = 10, seed = 1, threads = 2)
  expect_true(g$directed)
  expect_equal(nrow(g$edges), 2 * 5 + 10)
  m <- as.matrix(g)
  expect_equal(sum(m * t(m)) / 2, 5)
  expect_equal(sum(m) - 2 * 5, 10)
  expect_equal(net_random_uman(8, 5, 10, seed = 1, threads = 1)$edges, g$edges)
  expect_error(net_random_uman(4, 4, 3, seed = 1), "dyads")
})

test_that("small worlds and preferential attachment", {
  sw <- net_random_small_world(30, 4, 0.1, seed = 1, threads = 2)
  simple_undirected(sw)
  expect_equal(nrow(sw$edges), 60)
  lattice <- net_random_small_world(30, 4, 0, seed = 1)
  expect_equal(net_degree(lattice), rep(4, 30))
  expect_error(net_random_small_world(10, 3, 0.1, seed = 1), "even|k")
  ba <- net_random_pref(30, 2, seed = 1, threads = 2)
  simple_undirected(ba)
  expect_equal(nrow(ba$edges), 2 * 28)
  expect_true(all(net_degree(ba) >= 1))
  expect_error(net_random_pref(5, 5, seed = 1), "m")
})

test_that("stochastic block models", {
  blocks <- net_random_sbm(c(5, 5), matrix(c(1, 0, 0, 1), 2), seed = 1, threads = 2)
  expect_equal(nrow(blocks$edges), 2 * choose(5, 2))
  expect_true(all((blocks$edges$source <= 5) == (blocks$edges$target <= 5)))
  expect_equal(net_components(blocks)$no, 2)
  d <- net_random_sbm(c(3, 3), matrix(c(0, 1, 0, 0), 2, byrow = TRUE), directed = TRUE, seed = 1)
  # Block 1 -> block 2 only.
  expect_equal(nrow(d$edges), 9)
  expect_true(all(d$edges$source <= 3 & d$edges$target > 3))
  expect_error(net_random_sbm(c(2, 2), matrix(c(0.5, 0.2, 0.3, 0.5), 2), seed = 1), "symmetric")
  expect_error(net_random_sbm(c(2, 2), matrix(2, 2, 2), seed = 1), "probabilities")
})
