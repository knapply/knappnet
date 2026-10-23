# Degree distributions, component sizes, the triad census against its
# expectation, and ranked scores (R/plot-distributions.R).

karate <- function() net_example("karate")

test_that("degree distributions are tabulated and drawn in both forms", {
  g <- karate()
  on_null_device({
    expect_silent(out <- net_plot_degree(g))
    d <- net_degree(g)
    expect_identical(out$degree, sort(unique(d)))
    expect_identical(sum(out$count), length(d))
    expect_equal(out$ccdf[1], 1)
    expect_equal(out$ccdf, rev(cumsum(rev(out$count))) / length(d))
    expect_silent(net_plot_degree(g, type = "ccdf"))
    expect_silent(net_plot_degree(g, type = "ccdf", log = TRUE))
    expect_silent(net_plot_degree(g, type = "histogram", log = TRUE, main = "Log"))
    expect_restores_par(net_plot_degree(g))
    directed <- net_frame(data.frame(source = c(1, 1, 2), target = c(2, 3, 3)),
                         nodes = data.frame(id = 1:4))
    expect_silent(out <- net_plot_degree(directed, mode = "in"))
    expect_identical(out$degree, c(0, 1, 2))
    # Isolates are left out of logarithmic axes with a note, not an error.
    expect_silent(net_plot_degree(directed, type = "ccdf", log = TRUE))
    only_isolates <- net_frame(data.frame(source = 1, target = 1)[0, ], nodes = data.frame(id = 1:3))
    expect_error(net_plot_degree(only_isolates, log = TRUE), "degree 0")
    expect_silent(net_plot_degree(only_isolates))
  })
})

test_that("component sizes are drawn", {
  g <- net_example("noordin", relation = "friendship")
  on_null_device({
    expect_silent(out <- net_plot_components(g))
    comp <- net_components(g)
    expect_identical(sum(out$count), comp$no)
    expect_identical(sum(out$size * out$count), length(comp$membership))
    expect_silent(net_plot_components(g, log = TRUE, main = "Sizes"))
    expect_silent(net_plot_components(karate(), mode = "strong"))
    expect_restores_par(net_plot_components(g))
  })
})

test_that("the expected triad census under U|MAN is right", {
  # Three nodes, one asymmetric arc: the only triad is 012 with certainty.
  g <- net_frame(data.frame(source = 1, target = 2), nodes = data.frame(id = 1:3))
  e <- expected_triad_census(g)
  expect_equal(e, c(0, 1, rep(0, 14)))
  # Four nodes and two arcs 1->2, 2->1 (one mutual dyad among six): each of
  # the four triads holds the mutual dyad with probability 3/6, giving 102
  # with probability 1/2 and 003 with 1/2.
  m <- net_frame(data.frame(source = c(1, 2), target = c(2, 1)), nodes = data.frame(id = 1:4))
  e <- expected_triad_census(m)
  expect_equal(e[c(1, 3)], c(2, 2))
  expect_equal(sum(e), choose(4, 3))
  # Expected counts add up to the number of triads, for any network, and
  # the mixture over orientations follows the 1/4, 1/4, 1/2 and 3/4, 1/4
  # shares.
  set.seed(3)
  r <- net_frame(data.frame(source = sample(12, 40, TRUE), target = sample(12, 40, TRUE)),
                nodes = data.frame(id = 1:12))
  e <- expected_triad_census(r)
  expect_equal(sum(e), choose(12, 3))
  expect_equal(e[4], e[5])
  expect_equal(e[6], 2 * e[4])
  expect_equal(e[9], 3 * e[10])
  expect_equal(e[12], e[13])
  expect_equal(e[14], 2 * e[12])
  # Undirected networks: every edge is mutual, so only 003, 102, 201, 300.
  u <- karate()
  e <- expected_triad_census(u)
  expect_equal(sum(e), choose(34, 3))
  expect_true(all(e[-c(1, 3, 11, 16)] == 0))
  # Fewer than three nodes: no triads.
  expect_equal(expected_triad_census(net_frame(data.frame(source = 1, target = 2))), rep(0, 16))
  # A complete digraph has only 300 triads, observed and expected.
  pairs <- expand.grid(source = 1:4, target = 1:4)
  pairs <- pairs[pairs$source != pairs$target, ]
  k <- net_frame(pairs)
  expect_equal(expected_triad_census(k), c(rep(0, 15), 4))
  expect_equal(unname(net_triad_census(k)), c(rep(0, 15), 4))
})

test_that("net_plot_triads() draws observed and expected counts", {
  g <- net_frame(data.frame(source = c(1, 2, 3, 3, 4, 5, 5, 6, 1),
                           target   = c(2, 3, 1, 4, 5, 3, 6, 4, 6)))
  on_null_device({
    expect_silent(out <- net_plot_triads(g, threads = 1))
    expect_identical(nrow(out), 16L)
    expect_equal(out$observed, unname(net_triad_census(g)))
    expect_equal(sum(out$expected), choose(6, 3))
    expect_silent(net_plot_triads(g, expected = FALSE, threads = 1))
    expect_silent(net_plot_triads(g, log = TRUE, threads = 1, main = "Triads"))
    expect_silent(out <- net_plot_triads(karate(), threads = 1))
    expect_identical(out$type, c("003", "102", "201", "300"))
    expect_restores_par(net_plot_triads(g, threads = 1))
  })
})

test_that("net_plot_scores() ranks scores and labels them", {
  g <- net_example("florentine")
  b <- net_betweenness(g)
  on_null_device({
    expect_silent(out <- net_plot_scores(b, g))
    expect_identical(out$label[1], "Medici")
    expect_identical(out$score, sort(b, decreasing = TRUE))
    expect_silent(out <- net_plot_scores(b, g, top = 5, main = "Top five", xlab = "betweenness"))
    expect_identical(nrow(out), 5L)
    expect_silent(out <- net_plot_scores(b, g, decreasing = FALSE, top = 3))
    expect_identical(out$score, sort(b)[1:3])
    expect_silent(net_plot_scores(c(a = 1, b = 3, c = 2)))
    expect_silent(net_plot_scores(c(1, NA, NaN, 4)))
    expect_silent(net_plot_scores(b, g$nodes$name))
    expect_restores_par(net_plot_scores(b, g))
    expect_error(net_plot_scores(b, c("x", "y")), "one label per score")
    expect_error(net_plot_scores("x"), "numeric")
    expect_error(net_plot_scores(c(NA_real_, NaN)), "no defined scores")
  })
})
