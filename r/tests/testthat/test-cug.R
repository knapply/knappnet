# Conditional uniform graph tests: the object, the p-values, the null
# models, and statistical agreement with sna::cug.test when it is
# installed.

test_that("cug_test returns the observed statistic and a null distribution", {
  g <- example_network("karate")
  ct <- cug_test(g, transitivity, conditioning = "edges", reps = 50, seed = 1, threads = 2)
  expect_s3_class(ct, "cug_test")
  expect_equal(ct$observed, transitivity(g))
  expect_length(ct$null, 50)
  expect_equal(ct$reps, 50)
  expect_equal(ct$dropped, 0)
  expect_equal(ct$p_greater, mean(ct$null >= ct$observed))
  expect_equal(ct$p_less, mean(ct$null <= ct$observed))
  expect_equal(ct$statistic, "transitivity")
  # The karate club is far more clustered than G(n, m).
  expect_equal(ct$p_greater, 0)
  expect_equal(ct$p_less, 1)
  # Reproducible by seed, whatever the thread count.
  again <- cug_test(g, transitivity, conditioning = "edges", reps = 50, seed = 1, threads = 1)
  expect_equal(again$null, ct$null)
  expect_output(print(ct), "Pr\\(X >= Obs\\)")
  s <- summary(ct)
  expect_s3_class(s, "summary.cug_test")
  expect_output(print(s), "Null distribution")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(ct))
  expect_invisible(plot(ct, main = "x", breaks = 10))
})

test_that("null models keep what they condition on", {
  g <- example_network("karate")
  m <- nrow(g$edges)
  edges <- cug_test(g, function(x) nrow(x$edges), "edges", reps = 20, seed = 1, threads = 2)
  expect_true(all(edges$null == m))
  size <- cug_test(g, function(x) nrow(x$nodes), "size", reps = 20, seed = 1, threads = 2)
  expect_true(all(size$null == 34))
  # Density 1/2 under the size null.
  dens <- cug_test(g, graph_density, "size", reps = 200, seed = 1, threads = 2)
  expect_equal(mean(dens$null), 0.5, tolerance = 0.03)
  # Undirected dyad census is the edge count.
  dc <- cug_test(g, function(x) nrow(x$edges), "dyad.census", reps = 20, seed = 1, threads = 2)
  expect_true(all(dc$null == m))
  expect_equal(dc$conditioning, "dyad.census")

  d <- random_uman(12, mutual = 10, asym = 15, seed = 3, threads = 2)
  census <- cug_test(d, function(x) {
    a <- as.matrix(x)
    sum(a * t(a)) / 2
  }, "dyad.census", reps = 20, seed = 1, threads = 2)
  expect_true(all(census$null == 10))
  expect_true(all(cug_test(d, function(x) nrow(x$edges), "dyad.census", reps = 20, seed = 1,
                           threads = 2)$null == 35))
  expect_true(all(cug_test(d, function(x) as.numeric(x$directed), "edges", reps = 5, seed = 1,
                           threads = 2)$null == 1))
})

test_that("statistics with extra arguments, undefined values, and bad input", {
  g <- example_network("karate")
  ct <- cug_test(g, centralization, "edges", measure = "degree", reps = 20, seed = 1,
                 threads = 2)
  expect_equal(ct$observed, centralization(g, "degree"))
  # A statistic that is undefined on some random networks is dropped there.
  dropped <- cug_test(g, function(x) if (nrow(x$edges) > 0) NaN else 1, "edges", reps = 10,
                      seed = 1, threads = 2)
  expect_equal(dropped$dropped, 10)
  expect_length(dropped$null, 0)
  expect_true(is.nan(dropped$p_greater))
  expect_error(cug_test(g, function(x) degree(x), reps = 2, seed = 1), "single number")
  expect_error(cug_test(g, "transitivity", reps = 2, seed = 1), "function")
  expect_error(cug_test(g, transitivity, reps = 0, seed = 1), "at least 1")
  loops <- netframe(data.frame(from = c(1, 1, 2), to = c(1, 2, 3)), directed = FALSE)
  expect_warning(cug_test(loops, transitivity, reps = 2, seed = 1), "self-loops")
})

test_that("p-values agree with sna::cug.test", {
  skip_if_not_installed("sna")
  g <- example_network("florentine")
  a <- as.matrix(g)
  set.seed(1)
  # Transitivity conditioned on the number of edges: Monte Carlo estimates
  # from both packages of the same probability.
  ours <- cug_test(g, transitivity, "edges", reps = 1000, seed = 1, threads = 2)
  theirs <- sna::cug.test(a, sna::gtrans, mode = "graph", cmode = "edges", reps = 1000,
                          FUN.args = list(mode = "graph"))
  expect_equal(ours$observed, theirs$obs.stat, tolerance = 1e-10)
  # Each p-value has a standard error of at most 0.016 at 1000 draws.
  expect_true(abs(ours$p_greater - theirs$pgteobs) < 0.08)
  expect_true(abs(ours$p_less - theirs$plteobs) < 0.08)
  expect_true(abs(mean(ours$null) - mean(theirs$rep.stat)) < 0.3 * sd(theirs$rep.stat))
})
