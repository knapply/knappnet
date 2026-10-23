# Plot methods and plotting functions are checked by behavior on a null
# device: they must draw without error or warning, return their input or
# the plotted data invisibly, restore graphical parameters, and compute the
# statistics they display correctly.

on_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

expect_restores_par <- function(code) {
  before <- graphics::par(c("mar", "mfrow", "xpd"))
  force(code)
  expect_identical(graphics::par(c("mar", "mfrow", "xpd")), before)
}

karate <- function() example_network("karate")

test_that("plot_matrix() draws adjacency and incidence matrices in any order", {
  g <- karate()
  on_null_device({
    expect_silent(m <- plot_matrix(g))
    expect_identical(dim(m), c(34L, 34L))
    expect_identical(m, as.matrix(g))
    expect_silent(m2 <- plot_matrix(g, order = "faction", main = "Sorted"))
    expect_setequal(rownames(m2), rownames(m))
    # Rows are sorted by group in order of first appearance, stably.
    fac <- g$nodes$faction
    expect_identical(rownames(m2), as.character(g$nodes$name[order(match(fac, unique(fac)))]))
    expect_silent(plot_matrix(g, order = cluster_louvain(g, seed = 1), labels = FALSE))
    expect_silent(plot_matrix(g, values = TRUE, key = TRUE, labels = TRUE))
    expect_restores_par(plot_matrix(g))
    # Plain matrices, with and without dimnames and with a rectangular shape.
    expect_silent(plot_matrix(matrix(c(0, 1, 2, 0), 2)))
    expect_silent(plot_matrix(matrix(1:6, 2), order = c("b", "a")))
    expect_error(plot_matrix(matrix(1:4, 2), order = 1:3), "one value per row")
    expect_error(plot_matrix("x"), "netframe or a numeric matrix")
    women <- example_network("southern_women")
    expect_silent(inc <- plot_matrix(women))
    expect_identical(inc, incidence_matrix(women))
    expect_silent(plot_matrix(women, order = rep(c("x", "y"), length.out = nrow(women$nodes))))
    # Weighted values fill the matrix.
    w <- example_network("karate", weighted = TRUE)
    expect_silent(wm <- plot_matrix(w, attr = "contexts"))
    expect_identical(wm, as.matrix(w, attr = "contexts"))
  })
})

test_that("plot_mixing_matrix() draws the mixing matrix of a network or a matrix", {
  g <- karate()
  on_null_device({
    expect_silent(m <- plot_mixing_matrix(g, "faction"))
    expect_identical(m, mixing_matrix(g, "faction"))
    expect_silent(plot_mixing_matrix(mixing_matrix(g, "faction"), main = "Factions"))
    expect_silent(plot_mixing_matrix(g, "faction", "name", values = FALSE))
    expect_error(plot_mixing_matrix(g), "dim1")
    expect_error(plot_mixing_matrix(list()), "netframe or a numeric matrix")
  })
})

test_that("plot.blockmodel() draws densities and the permuted matrix", {
  g <- example_network("florentine")
  bm <- blockmodel(g, concor(g, depth = 2))
  on_null_device({
    expect_silent(out <- plot(bm))
    expect_identical(out, bm)
    expect_silent(plot(bm, network = g))
    expect_silent(plot(bm, network = g, which = "matrix"))
    expect_silent(plot(bm, network = g, which = "density", values = FALSE, main = "Densities"))
    expect_restores_par(plot(bm, network = g))
    expect_error(plot(bm, which = "matrix"), "needs `network`")
    other <- netframe(data.frame(from = 1, to = 2))
    expect_error(plot(bm, network = other), "one node per element")
    # A one-node block without loops has a NaN density and still draws.
    one <- blockmodel(g, c("a", rep("b", nrow(g$nodes) - 1L)))
    expect_true(is.nan(one$density["a", "a"]))
    expect_silent(plot(one, network = g))
  })
})

test_that("QAP plots show the permutation distribution and refuse empty ones", {
  set.seed(1)
  a <- matrix(rbinom(100, 1, 0.3), 10)
  b <- (a + matrix(rbinom(100, 1, 0.1), 10)) > 0
  q <- qap_cor(a, b, reps = 200, seed = 1, threads = 1)
  on_null_device({
    expect_silent(out <- plot(q))
    expect_identical(out, q)
    expect_silent(plot(q, breaks = 10, main = "QAP", col = "gray50"))
    expect_restores_par(plot(q))
    none <- qap_cor(a, b, reps = 0)
    expect_error(plot(none), "no permutation draws")
  })
  set.seed(2)
  x1 <- matrix(rnorm(144), 12)
  x2 <- matrix(rnorm(144), 12)
  y <- 0.8 * x1 + matrix(rnorm(144), 12)
  fit <- mrqap(y, list(same = x1, other = x2), reps = 200, seed = 1, threads = 1)
  on_null_device({
    expect_silent(out <- plot(fit))
    expect_identical(out, fit)
    expect_silent(plot(fit, intercept = TRUE, level = 0.9))
    expect_silent(plot(fit, which = "distributions"))
    expect_silent(plot(fit, which = "distributions", intercept = TRUE, main = "Null"))
    expect_restores_par(plot(fit, which = "distributions"))
    beta <- mrqap(y, x1, reps = 50, statistic = "beta", method = "y", seed = 1, threads = 1)
    expect_silent(plot(beta))
    expect_error(plot(fit, level = 1), "`level`")
    expect_error(plot(mrqap(y, x1, reps = 0)), "no permutation draws")
  })
})

test_that("classical_mds() results print and plot", {
  g <- example_network("florentine")
  d <- mds_distances(g)
  fit <- classical_mds(d, k = 2)
  expect_s3_class(fit, "classical_mds")
  expect_named(fit, c("points", "eigenvalues", "gof", "stress"))
  expect_output(print(fit), "15 points in 2 dimensions")
  on_null_device({
    expect_silent(out <- plot(fit))
    expect_identical(out, fit)
    expect_silent(plot(fit, labels = FALSE, main = "MDS"))
    expect_silent(plot(fit, which = "shepard", d = d))
    expect_silent(plot(fit, which = "shepard", d = stats::as.dist(d)))
    expect_silent(plot(fit, which = "scree"))
    expect_restores_par(plot(fit, which = "scree"))
    expect_error(plot(fit, which = "shepard"), "pass `d`")
    expect_error(plot(fit, which = "shepard", d = d[1:5, 1:5]), "5 points")
    expect_error(plot(fit, dims = c(1, 3)), "at most 2")
    expect_error(plot(fit, dims = 1), "two dimension")
    one <- classical_mds(d, k = 1)
    expect_silent(plot(one))
    expect_output(print(one), "1 dimension\n")
  })
})

test_that("degree distributions are tabulated and drawn in both forms", {
  g <- karate()
  on_null_device({
    expect_silent(out <- plot_degree_distribution(g))
    d <- degree(g)
    expect_identical(out$degree, sort(unique(d)))
    expect_identical(sum(out$count), length(d))
    expect_equal(out$ccdf[1], 1)
    expect_equal(out$ccdf, rev(cumsum(rev(out$count))) / length(d))
    expect_silent(plot_degree_distribution(g, type = "ccdf"))
    expect_silent(plot_degree_distribution(g, type = "ccdf", log = TRUE))
    expect_silent(plot_degree_distribution(g, type = "histogram", log = TRUE, main = "Log"))
    expect_restores_par(plot_degree_distribution(g))
    directed <- netframe(data.frame(from = c(1, 1, 2), to = c(2, 3, 3)),
                         nodes = data.frame(id = 1:4))
    expect_silent(out <- plot_degree_distribution(directed, mode = "in"))
    expect_identical(out$degree, c(0, 1, 2))
    # Isolates are left out of logarithmic axes with a note, not an error.
    expect_silent(plot_degree_distribution(directed, type = "ccdf", log = TRUE))
    only_isolates <- netframe(data.frame(from = 1, to = 1)[0, ], nodes = data.frame(id = 1:3))
    expect_error(plot_degree_distribution(only_isolates, log = TRUE), "degree 0")
    expect_silent(plot_degree_distribution(only_isolates))
  })
})

test_that("component sizes are drawn", {
  g <- example_network("noordin", relation = "friendship")
  on_null_device({
    expect_silent(out <- plot_component_sizes(g))
    comp <- components(g)
    expect_identical(sum(out$count), comp$no)
    expect_identical(sum(out$size * out$count), length(comp$membership))
    expect_silent(plot_component_sizes(g, log = TRUE, main = "Sizes"))
    expect_silent(plot_component_sizes(karate(), mode = "strong"))
    expect_restores_par(plot_component_sizes(g))
  })
})

test_that("the expected triad census under U|MAN is right", {
  # Three nodes, one asymmetric arc: the only triad is 012 with certainty.
  g <- netframe(data.frame(from = 1, to = 2), nodes = data.frame(id = 1:3))
  e <- expected_triad_census(g)
  expect_equal(e, c(0, 1, rep(0, 14)))
  # Four nodes and two arcs 1->2, 2->1 (one mutual dyad among six): each of
  # the four triads holds the mutual dyad with probability 3/6, giving 102
  # with probability 1/2 and 003 with 1/2.
  m <- netframe(data.frame(from = c(1, 2), to = c(2, 1)), nodes = data.frame(id = 1:4))
  e <- expected_triad_census(m)
  expect_equal(e[c(1, 3)], c(2, 2))
  expect_equal(sum(e), choose(4, 3))
  # Expected counts add up to the number of triads, for any network, and
  # the mixture over orientations follows the 1/4, 1/4, 1/2 and 3/4, 1/4
  # shares.
  set.seed(3)
  r <- netframe(data.frame(from = sample(12, 40, TRUE), to = sample(12, 40, TRUE)),
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
  expect_equal(expected_triad_census(netframe(data.frame(from = 1, to = 2))), rep(0, 16))
  # A complete digraph has only 300 triads, observed and expected.
  pairs <- expand.grid(from = 1:4, to = 1:4)
  pairs <- pairs[pairs$from != pairs$to, ]
  k <- netframe(pairs)
  expect_equal(expected_triad_census(k), c(rep(0, 15), 4))
  expect_equal(unname(triad_census(k)), c(rep(0, 15), 4))
})

test_that("plot_triad_census() draws observed and expected counts", {
  g <- netframe(data.frame(from = c(1, 2, 3, 3, 4, 5, 5, 6, 1),
                           to   = c(2, 3, 1, 4, 5, 3, 6, 4, 6)))
  on_null_device({
    expect_silent(out <- plot_triad_census(g, threads = 1))
    expect_identical(nrow(out), 16L)
    expect_equal(out$observed, unname(triad_census(g)))
    expect_equal(sum(out$expected), choose(6, 3))
    expect_silent(plot_triad_census(g, expected = FALSE, threads = 1))
    expect_silent(plot_triad_census(g, log = TRUE, threads = 1, main = "Triads"))
    expect_silent(out <- plot_triad_census(karate(), threads = 1))
    expect_identical(out$type, c("003", "102", "201", "300"))
    expect_restores_par(plot_triad_census(g, threads = 1))
  })
})

test_that("plot_scores() ranks scores and labels them", {
  g <- example_network("florentine")
  b <- betweenness(g)
  on_null_device({
    expect_silent(out <- plot_scores(b, g))
    expect_identical(out$label[1], "Medici")
    expect_identical(out$score, sort(b, decreasing = TRUE))
    expect_silent(out <- plot_scores(b, g, top = 5, main = "Top five", xlab = "betweenness"))
    expect_identical(nrow(out), 5L)
    expect_silent(out <- plot_scores(b, g, decreasing = FALSE, top = 3))
    expect_identical(out$score, sort(b)[1:3])
    expect_silent(plot_scores(c(a = 1, b = 3, c = 2)))
    expect_silent(plot_scores(c(1, NA, NaN, 4)))
    expect_silent(plot_scores(b, g$nodes$name))
    expect_restores_par(plot_scores(b, g))
    expect_error(plot_scores(b, c("x", "y")), "one label per score")
    expect_error(plot_scores("x"), "numeric")
    expect_error(plot_scores(c(NA_real_, NaN)), "no defined scores")
  })
})

test_that("network overlays draw centrality, partitions, cores, and egos", {
  g <- karate()
  on_null_device({
    expect_silent(xy <- plot_centrality(g, betweenness(g), title = "betweenness", seed = 1))
    expect_identical(xy, layout_network(g, seed = 1))
    expect_silent(plot_centrality(g, degree(g), color = FALSE, layout = "circle"))
    g2 <- g
    g2$nodes$deg <- degree(g)
    expect_silent(plot_centrality(g2, "deg", seed = 1))
    expect_error(plot_centrality(g, 1:3), "one value per node")
    expect_error(plot_centrality(g, degree(g), title = 1), "single string")
    # The temporary legend column never touches the caller's network.
    expect_identical(names(g$nodes), c("name", "faction"))

    groups <- cluster_louvain(g, seed = 1)
    expect_silent(xy <- plot_partition(g, groups, seed = 1))
    expect_identical(xy, layout_network(g, seed = 1))
    expect_silent(plot_partition(g, groups, seed = 1, main = "Communities"))
    expect_silent(plot_partition(g, "faction", seed = 1))
    expect_silent(plot_partition(g, components(g)$membership, title = "component", layout = "circle"))
    expect_error(plot_partition(g, 1:3), "one value per node")

    expect_silent(xy <- plot_coreness(g, seed = 1, threads = 1))
    expect_identical(dim(xy), c(34L, 2L))
    k <- coreness(g)
    radius <- sqrt(rowSums(xy^2))
    # Higher cores sit on inner rings.
    expect_true(all(diff(tapply(radius, k, mean)) <= 0))
    expect_silent(plot_coreness(g, seed = 1, threads = 1, main = "Cores", node_label = FALSE))
    directed <- netframe(data.frame(from = c(1, 2, 3, 1), to = c(2, 3, 1, 4)))
    expect_silent(plot_coreness(directed, mode = "in", seed = 1, threads = 1))
    expect_restores_par(plot_coreness(g, seed = 1, threads = 1))

    f <- example_network("florentine")
    expect_silent(e <- plot_ego(f, "Medici", seed = 1))
    expect_s3_class(e, "netframe")
    expect_identical(sort(e$nodes$name), sort(ego_network(f, "Medici")$nodes$name))
    expect_silent(plot_ego(f, "Medici", order = 2, seed = 1, main = "Two steps"))
    expect_error(plot_ego(f, "Nobody"), "not in the network")
  })
})

test_that("plot_robustness() tracks the largest component as nodes are removed", {
  g <- example_network("jemaah_islamiyah")
  n <- nrow(g$nodes)
  on_null_device({
    expect_silent(out <- plot_robustness(g, steps = 5, seed = 1, threads = 1))
    expect_setequal(unique(out$curve), c("degree", "betweenness", "random"))
    first <- out[out$removed == 0, ]
    expect_true(all(first$value == max(components(g)$csize) / n))
    expect_true(all(out$value >= 0 & out$value <= 1))
    # Removing the top-degree node by hand gives the same value as the curve.
    top <- g$nodes$name[which.max(degree(g))]
    expect_silent(by_hand <- plot_robustness(g, nodes = top, strategy = NULL, threads = 1))
    expect_identical(by_hand$curve, rep("selected nodes", 2))
    rest <- netframe(g$edges[g$edges$from != top & g$edges$to != top, ],
                     g$nodes[g$nodes$name != top, ], directed = FALSE)
    expect_equal(by_hand$value[2], max(components(rest)$csize) / n)
    kp <- key_players(g, 2, "fragmentation", seed = 1, threads = 1)
    expect_silent(out <- plot_robustness(g, nodes = kp$set, strategy = "degree",
                                         measure = "fragmentation", label = "key players",
                                         steps = 4, threads = 1, main = "Attack"))
    expect_setequal(unique(out$curve), c("degree", "key players"))
    expect_equal(out$value[out$curve == "key players"][3], kp$score)
    expect_restores_par(plot_robustness(g, strategy = "degree", steps = 3, threads = 1))
    expect_error(plot_robustness(g, strategy = "closeness"), "`strategy`")
    expect_error(plot_robustness(g, nodes = NULL, strategy = NULL), "Give `nodes`")
    expect_error(plot_robustness(g, nodes = "x", label = "degree"), "not the name of a strategy")
    expect_error(plot_robustness(g, steps = 1), "`steps`")
    expect_error(plot_robustness(netframe(data.frame(from = 1, to = 1))), "at least two nodes")
  })
})

test_that("temporal plots draw timelines and measures over time", {
  calls <- data.frame(from = c("a", "a", "b", "c", "a", "b"),
                      to   = c("b", "c", "c", "d", "b", "c"),
                      time = c(1, 2, 2, 5, 7, 9), end = c(3, 2, 4, 6, 7, 10))
  g <- temporal_network(calls, end = "end", directed = FALSE)
  on_null_device({
    expect_silent(out <- plot_timeline(g))
    expect_identical(out$row[1:2], c("a -- b", "b -- c"))
    expect_identical(out$contacts[1:2], c(2L, 2L))
    expect_identical(sum(out$contacts), 6L)
    expect_silent(out <- plot_timeline(g, by = "node", max_rows = 2, main = "Nodes"))
    expect_identical(nrow(out), 2L)
    expect_identical(out$row[1], "b")
    expect_restores_par(plot_timeline(g))
    directed <- temporal_network(calls, end = "end")
    expect_silent(out <- plot_timeline(directed))
    expect_identical(out$row[1], "a -> b")
    instant <- temporal_network(data.frame(from = "a", to = "a", time = 3))
    expect_silent(plot_timeline(instant))
    dated <- temporal_network(data.frame(from = c("a", "b"), to = c("b", "c"),
                                         time = as.Date(c("2024-01-01", "2024-03-01"))))
    expect_silent(plot_timeline(dated))
    stamped <- temporal_network(data.frame(from = c("a", "b"), to = c("b", "c"),
                                           time = as.POSIXct(c("2024-01-01 10:00", "2024-01-01 12:30"),
                                                             tz = "UTC")))
    expect_silent(plot_timeline(stamped))
    empty <- temporal_network(calls[0, ], nodes = data.frame(id = c("a", "b")), end = "end")
    expect_error(plot_timeline(empty), "no contacts")
    expect_error(plot_timeline(netframe(calls)), "temporal network")
  })
  logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3", "u1"),
                       host = c("h1", "h1", "h2", "h2", "h3", "h1", "h1"),
                       time = c(0, 1, 5, 6, 12, 13, 14))
  h <- temporal_network(logins)
  r <- rolling(h, degree, width = 5, step = 2)
  on_null_device({
    expect_silent(out <- plot_rolling(r))
    expect_identical(out, r)
    expect_silent(plot_rolling(r, nodes = c("u1", "h1"), main = "Two nodes"))
    expect_silent(plot_rolling(r, top = 2))
    expect_silent(plot_rolling(rolling(h, function(w) nrow(w$edges), width = 5), ylab = "contacts"))
    expect_restores_par(plot_rolling(r))
    expect_error(plot_rolling(r, nodes = "zz"), "not in the results")
    expect_error(plot_rolling(data.frame(a = 1)), "rolling()")
    expect_error(plot_rolling(r[0, ]), "no windows")
    dated <- temporal_network(data.frame(from = c("a", "b", "a"), to = c("b", "c", "c"),
                                         time = as.Date(c("2024-01-01", "2024-01-05", "2024-01-09"))))
    expect_silent(plot_rolling(rolling(dated, degree, width = 4)))
  })
})
