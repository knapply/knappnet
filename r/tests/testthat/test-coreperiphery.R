# Core-periphery fits, the small-world quotient, compactness, average
# degree, and the degree distribution. Hand-checked values; the full
# comparison with exhaustive search and NetworkX lives in validation/.

# Node 2 ties to 1, 4, 5, 6, 7; 5 ties to 2, 3, 4, 7; 7 ties to 1, 2, 5.
# The core {2, 5, 7} is a triangle and the periphery has no ties among
# itself, so the fit with core-periphery pairs left out is perfect.
cp_net <- function() {
  netframe(data.frame(from = c(1, 1, 2, 2, 2, 2, 3, 4, 5), to = c(2, 7, 4, 5, 6, 7, 5, 5, 7)),
           directed = FALSE)
}

test_that("discrete core-periphery fits by every link rule", {
  g <- cp_net()
  cp <- core_periphery(g, seed = 1, threads = 2)
  expect_s3_class(cp, "core_periphery")
  expect_equal(g$nodes$name[cp$core], c(2, 5, 7))
  expect_equal(cp$correlation, 1)
  expect_equal(core_periphery(g, seed = 1, threads = 1)$core, cp$core)
  # Exhaustive optima for the other rules (see validation/).
  expect_equal(g$nodes$name[core_periphery(g, links = "ones", seed = 1, threads = 2)$core],
               c(2, 5))
  expect_equal(g$nodes$name[core_periphery(g, links = "zeros", seed = 1, threads = 2)$core],
               c(1, 2, 4, 5, 7))
  expect_equal(core_periphery(g, links = "ones", seed = 1, threads = 2)$correlation,
               0.6330542, tolerance = 1e-6)
  expect_output(print(cp), "Correlation")
  # No structure to fit.
  complete <- netframe(data.frame(from = c(1, 1, 2), to = c(2, 3, 3)), directed = FALSE)
  expect_error(core_periphery(complete, seed = 1), "different values")
})

test_that("continuous core-periphery fits", {
  g <- cp_net()
  cc <- core_periphery(g, method = "continuous", threads = 2)
  expect_equal(cc$method, "continuous")
  expect_equal(cc$fit, "minres")
  expect_true(all(cc$coreness >= 0))
  expect_equal(sum(cc$coreness^2), 1, tolerance = 1e-9)
  expect_equal(cc$correlation, 0.6957548, tolerance = 1e-6)
  expect_equal(cc$residual, 2.669288, tolerance = 1e-6)
  expect_equal(which.max(cc$coreness), 2)
  # The reported correlation is that of the pattern c_i c_j, and the
  # residual that of the least-squares scaling of the pattern.
  a <- as.matrix(g)
  iu <- upper.tri(a)
  pattern <- outer(cc$coreness, cc$coreness)
  expect_equal(cor(a[iu], pattern[iu]), cc$correlation, tolerance = 1e-8)
  scale <- sum(a[iu] * pattern[iu]) / sum(pattern[iu]^2)
  expect_equal(sum((a[iu] - scale * pattern[iu])^2), cc$residual, tolerance = 1e-6)
  expect_output(print(cc), "MINRES")

  # The correlation fit: also unit norm, with its own correlation.
  cr <- core_periphery(g, method = "continuous", fit = "correlation", seed = 1, threads = 2)
  expect_equal(cr$fit, "correlation")
  expect_true(is.nan(cr$residual))
  expect_equal(sum(cr$coreness^2), 1, tolerance = 1e-9)
  pattern <- outer(cr$coreness, cr$coreness)
  expect_equal(cor(a[iu], pattern[iu]), cr$correlation, tolerance = 1e-6)
  expect_true(cr$correlation >= cc$correlation - 1e-9)
  expect_equal(core_periphery(g, method = "continuous", fit = "correlation", seed = 1,
                              threads = 1)$coreness, cr$coreness)

  # A star: the center is the core and the leaves are symmetric under both.
  star <- netframe(data.frame(from = 1, to = 2:6), directed = FALSE)
  for (fit in c("minres", "correlation")) {
    s <- core_periphery(star, method = "continuous", fit = fit, seed = 1, threads = 2)
    expect_equal(which.max(s$coreness), 1)
    expect_equal(diff(range(s$coreness[-1])), 0, tolerance = 1e-6)
    expect_equal(s$correlation, 1, tolerance = 1e-6)
  }
  expect_error(core_periphery(g, method = "continuous", tol = -1), "tol")
})

test_that("core-periphery plots draw", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  g <- cp_net()
  cp <- core_periphery(g, seed = 1, threads = 2)
  cc <- core_periphery(g, method = "continuous", threads = 2)
  expect_invisible(plot(cp))
  expect_invisible(plot(cp, what = "network", seed = 1))
  expect_invisible(plot(cc, seed = 1))
  expect_invisible(plot(cc, what = "matrix", main = "coreness"))
  expect_error(plot(cp, what = "x"), "what")
})

test_that("small-world quotient is a ratio of ratios", {
  g <- example_network("karate")
  sw <- small_world(g, reps = 20, seed = 1, threads = 2)
  expect_s3_class(sw, "small_world")
  expect_equal(sw$clustering, transitivity(g))
  expect_equal(sw$mean_distance, mean_distance(g))
  expect_equal(sw$quotient,
               (sw$clustering / sw$random_clustering) /
                 (sw$mean_distance / sw$random_mean_distance))
  expect_true(sw$quotient > 1)
  expect_equal(small_world(g, reps = 20, seed = 1, threads = 1), sw)
  local <- small_world(g, local = TRUE, reps = 20, seed = 1, threads = 2)
  expect_equal(local$clustering, average_clustering(g))
  expect_false(isTRUE(all.equal(local$quotient, sw$quotient)))
  expect_output(print(sw), "quotient")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(sw))
  expect_error(small_world(g, reps = 0), "at least 1")
})

test_that("compactness and distance-weighted fragmentation", {
  g <- netframe(data.frame(from = c(1, 2), to = c(2, 3)), nodes = data.frame(id = 1:4),
                directed = FALSE)
  expect_equal(compactness(g), (4 + 2 * 0.5) / 12)
  expect_equal(distance_weighted_fragmentation(g), 1 - 5 / 12)
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(compactness(d), (2 + 0.5) / 6)
  expect_equal(compactness(d, mode = "in"), (2 + 0.5) / 6)
  expect_equal(compactness(d, mode = "all"), 5 / 6)
  w <- netframe(data.frame(from = c(1, 2), to = c(2, 3), w = c(2, 2)), directed = FALSE,
                weight = "w")
  expect_equal(compactness(w), (4 / 2 + 2 / 4) / 6)
  expect_equal(compactness(w, weighted = FALSE), 5 / 6)
  expect_equal(compactness(w, weights_as = "strength"), (4 * 2 + 2 * 1) / 6)
  one <- netframe(data.frame(from = integer(), to = integer()), nodes = data.frame(id = 1))
  expect_true(is.nan(compactness(one)))
})

test_that("average degree and degree distribution", {
  g <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 3)), directed = FALSE)
  expect_equal(average_degree(g), 2)
  expect_equal(average_degree(g, mode = "out"), 2)
  dd <- degree_distribution(g, cumulative = TRUE)
  expect_equal(dd$degree, c(1, 2, 3))
  expect_equal(dd$count, c(1, 2, 1))
  expect_equal(dd$proportion, c(0.25, 0.5, 0.25))
  expect_equal(dd$at_least, c(1, 0.75, 0.25))
  expect_null(degree_distribution(g)$at_least)
  d <- netframe(data.frame(from = c(1, 1), to = c(2, 3)))
  expect_equal(average_degree(d, mode = "out"), 2 / 3)
  expect_equal(average_degree(d, mode = "in"), 2 / 3)
  expect_equal(average_degree(d), 4 / 3)
  empty <- netframe(data.frame(from = integer(), to = integer()),
                    nodes = data.frame(id = integer()))
  expect_true(is.nan(average_degree(empty)))
})
