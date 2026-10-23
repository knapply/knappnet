# Reciprocity, transitivity and clustering, distances, mean distance,
# diameter, eccentricity, radius, reachability, cut points, bridges, and
# summary() (R/descriptive.R). Hand-checked values first, then spot checks
# against igraph and sna. The full comparison lives in validation/.

random_edges <- function(n, m, seed) {
  set.seed(seed)
  e <- data.frame(source = sample(n, m, TRUE), target = sample(n, m, TRUE))
  e <- e[e$source != e$target, ]
  e[!duplicated(data.frame(a = pmin(e$source, e$target), b = pmax(e$source, e$target))), ]
}

test_that("reciprocity indices count dyads", {
  # a <-> b, b -> c, a loop at c and a repeated a -> b (both ignored), and
  # an isolate d: M = 1, A = 1, N = 4.
  g <- net_frame(data.frame(source = c("a", "b", "b", "c", "a"),
                           target   = c("b", "a", "c", "c", "b")),
                nodes = data.frame(id = c("a", "b", "c", "d")))
  expect_equal(net_reciprocity(g), 2 / 3)
  expect_equal(net_reciprocity(g, "ratio"), 1 / 2)
  expect_equal(net_reciprocity(g, "dyadic"), 5 / 6)
  # p = 3 / 12 and r = 2 / 3, so rho = (2/3 - 1/4) / (3/4).
  expect_equal(net_reciprocity(g, "correlation"), 5 / 9)
  empty <- net_frame(data.frame(source = character(0), target = character(0)),
                    nodes = data.frame(id = c("a", "b")))
  expect_true(is.nan(net_reciprocity(empty)))
  expect_equal(net_reciprocity(empty, "dyadic"), 1)
  expect_equal(net_reciprocity(pieces()), 1)
  expect_error(net_reciprocity(g, "edgewise"))
})

test_that("reciprocity matches igraph and sna", {
  skip_if_not_installed("igraph")
  set.seed(3)
  e <- data.frame(source = sample(25, 120, TRUE), target = sample(25, 120, TRUE))
  g <- net_frame(e, nodes = data.frame(id = 1:25))
  ig <- igraph::simplify(igraph::graph_from_data_frame(e, vertices = data.frame(id = 1:25)))
  expect_equal(net_reciprocity(g), igraph::reciprocity(ig))
  expect_equal(net_reciprocity(g, "ratio"), igraph::reciprocity(ig, mode = "ratio"))
  skip_if_not_installed("sna")
  m <- matrix(0, 25, 25)
  m[cbind(e$source, e$target)] <- 1
  diag(m) <- 0
  expect_equal(net_reciprocity(g), unname(sna::grecip(m, measure = "edgewise")))
  expect_equal(net_reciprocity(g, "ratio"), unname(sna::grecip(m, measure = "dyadic.nonnull")))
  expect_equal(net_reciprocity(g, "dyadic"), unname(sna::grecip(m, measure = "dyadic")))
  pairs <- which(row(m) != col(m))
  expect_equal(net_reciprocity(g, "correlation"), stats::cor(m[pairs], t(m)[pairs]))
})

test_that("transitivity and clustering by hand", {
  g <- pieces()
  # One connected triple (centred on 2), not closed.
  expect_equal(net_transitivity(g), 0)
  expect_equal(net_clustering(g), c(NaN, 0, NaN, NaN, NaN, NaN))
  expect_equal(net_clustering(g, convention = "networkx"), rep(0, 6))
  expect_equal(net_avg_clustering(g), 0)
  # A triangle with a pendant node and an isolate.
  t <- net_frame(data.frame(source = c(1, 2, 3, 3), target = c(2, 3, 1, 4)),
                nodes = data.frame(id = 1:5), directed = FALSE)
  expect_equal(net_transitivity(t), 3 / 5)
  expect_equal(net_clustering(t), c(1, 1, 1 / 3, NaN, NaN))
  expect_equal(net_avg_clustering(t), (2 + 1 / 3) / 3)
  expect_equal(net_avg_clustering(t, convention = "networkx"), (2 + 1 / 3) / 5)
  # Direction, parallel edges, and loops do not matter.
  d <- net_frame(data.frame(source = c(1, 2, 2, 3, 3, 3, 1), target = c(2, 1, 3, 1, 3, 4, 2)),
                nodes = data.frame(id = 1:4))
  expect_equal(net_clustering(d), c(1, 1, 1 / 3, NaN))
  no_triples <- net_frame(data.frame(source = 1, target = 2), directed = FALSE)
  expect_true(is.nan(net_transitivity(no_triples)))
  expect_equal(net_transitivity(no_triples, convention = "networkx"), 0)
})

test_that("weighted clustering by hand", {
  # Triangle 1-2-3 with weights 1, 2, 3 and a pendant 1-4 of weight 4.
  g <- net_frame(data.frame(source = c(1, 2, 3, 1), target = c(2, 3, 1, 4), w = c(1, 2, 3, 4)),
                directed = FALSE, weight = "w", weight_type = "distance")
  # Node 1: strength 8, three neighbours, one triangle with w12 + w13 = 4.
  expect_equal(net_clustering(g, weighted = TRUE), c(4 / 16, 1, 1, NaN))
  tri <- (1 * 2 * 3 / 4^3)^(1 / 3)
  expect_equal(net_clustering(g, weighted = TRUE, convention = "networkx"),
               c(2 * tri / 6, tri, tri, 0))
  expect_equal(net_clustering(g, weighted = TRUE, method = "onnela"),
               c(2 * tri / 6, tri, tri, NaN))
  expect_equal(net_avg_clustering(g, weighted = TRUE, convention = "networkx"),
               (2 * tri / 6 + 2 * tri) / 4)
  # Binary by default, whatever the weights.
  expect_equal(net_clustering(g), c(1 / 3, 1, 1, NaN))
  expect_error(net_clustering(pieces(), weighted = TRUE), "no weights")
  # Zero distances are refused when the network is built.
  expect_error(net_frame(data.frame(source = 1, target = 2, w = 0), weight = "w",
                         weight_type = "distance", directed = FALSE), "positive")
})

test_that("clustering matches igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(net_transitivity(g), igraph::transitivity(ig, type = "global"))
  expect_equal(net_clustering(g), unname(igraph::transitivity(ig, type = "local", isolates = "NaN")))
  expect_equal(net_avg_clustering(g), igraph::transitivity(ig, type = "average"))
  e <- random_edges(30, 90, 5)
  e$w <- round(stats::runif(nrow(e), 0.5, 4), 2)
  wg <- net_frame(e, nodes = data.frame(id = 1:30), directed = FALSE, weight = "w", weight_type = "distance")
  wig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(id = 1:30))
  expect_equal(net_clustering(wg, weighted = TRUE),
               unname(igraph::transitivity(wig, type = "barrat", weights = igraph::E(wig)$w,
                                           isolates = "NaN")))
})

test_that("distances and their summaries by hand", {
  g <- pieces()
  d <- net_distances(g)
  expect_equal(dimnames(d), list(as.character(1:6), as.character(1:6)))
  expect_equal(unname(d[1, ]), c(0, 1, 2, Inf, Inf, Inf))
  expect_equal(unname(d[6, ]), c(Inf, Inf, Inf, Inf, Inf, 0))
  expect_equal(net_distances(g, from = 3, to = c(1, 5)),
               matrix(c(2, Inf), 1, dimnames = list("3", c("1", "5"))))
  expect_error(net_distances(g, from = 7), "not in the network")
  # Ordered pairs joined by a path: 6 in the path (total 8) and 2 in the edge.
  expect_equal(net_mean_distance(g), 10 / 8)
  expect_equal(net_mean_distance(g, unconnected = FALSE), Inf)
  expect_equal(net_diameter(g), 2)
  expect_equal(net_diameter(g, unconnected = FALSE), Inf)
  expect_equal(net_eccentricity(g), c(2, 1, 2, 1, 1, 0))
  expect_equal(net_radius(g), 0)
  expect_error(net_mean_distance(g, unconnected = NA), "unconnected")

  # Directed, weighted: a -> b -> c -> d with minutes 5, 5, 20.
  w <- net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"),
                           minutes = c(5, 5, 20)), weight = "minutes", weight_type = "distance")
  expect_equal(unname(net_distances(w)[1, ]), c(0, 5, 10, 30))
  expect_equal(unname(net_distances(w, mode = "in")[4, ]), c(30, 25, 20, 0))
  expect_equal(unname(net_distances(w, weighted = FALSE, mode = "all")[4, ]), c(3, 2, 1, 0))
  expect_equal(net_mean_distance(w), (5 + 10 + 30 + 5 + 25 + 20) / 6)
  expect_equal(net_diameter(w), 30)
  expect_equal(net_diameter(w, weighted = FALSE), 3)
  expect_equal(net_eccentricity(w, mode = "out", weighted = FALSE), c(3, 2, 1, 0))
  # As strengths, lengths are 1/5, 1/5, and 1/20.
  expect_equal(net_diameter(w, weights_as = "strength"), 0.2 + 0.2 + 0.05)
})

test_that("distances match igraph and sna", {
  skip_if_not_installed("igraph")
  set.seed(9)
  e <- data.frame(source = sample(30, 70, TRUE), target = sample(30, 70, TRUE))
  e$w <- sample(1:6, nrow(e), TRUE)
  e$r <- round(stats::runif(nrow(e), 0.1, 3), 3)
  nodes <- data.frame(id = 1:30)
  ig <- igraph::graph_from_data_frame(e, vertices = nodes)
  g <- net_frame(e, nodes = nodes)
  for (mode in c("out", "in", "all")) {
    expect_equal(unname(net_distances(g, mode = mode)),
                 unname(igraph::distances(ig, mode = mode, weights = NA)))
  }
  gw <- net_frame(e, nodes = nodes, weight = "w", weight_type = "distance")
  expect_equal(unname(net_distances(gw)),
               unname(igraph::distances(ig, mode = "out", weights = e$w)))
  gr <- net_frame(e, nodes = nodes, weight = "r", weight_type = "distance")
  expect_equal(unname(net_distances(gr)),
               unname(igraph::distances(ig, mode = "out", weights = e$r)))
  expect_equal(net_mean_distance(g), igraph::mean_distance(ig))
  expect_equal(net_mean_distance(gr), igraph::mean_distance(ig, weights = e$r))
  expect_equal(net_diameter(g), igraph::diameter(ig))
  expect_equal(net_diameter(gr), igraph::diameter(ig, weights = e$r))
  expect_equal(net_diameter(g, unconnected = FALSE),
               igraph::diameter(ig, unconnected = FALSE))
  iu <- igraph::delete_edge_attr(igraph::delete_edge_attr(ig, "w"), "r")
  expect_equal(net_eccentricity(g), unname(igraph::eccentricity(iu, mode = "all")))
  expect_equal(net_eccentricity(g, mode = "out"), unname(igraph::eccentricity(iu, mode = "out")))
  expect_equal(net_radius(g), igraph::radius(iu, mode = "all"))
  skip_if_not_installed("sna")
  m <- matrix(0, 30, 30)
  m[cbind(e$source, e$target)] <- 1
  diag(m) <- 0
  expect_equal(unname(net_distances(g)), sna::geodist(m)$gdist)
})

test_that("reachability counts", {
  d <- net_frame(data.frame(source = c(1, 2, 4), target = c(2, 3, 5)),
                nodes = data.frame(id = 1:6))
  expect_equal(net_reachability(d), c(2, 1, 0, 1, 0, 0))
  expect_equal(net_reachability(d, "in"), c(0, 1, 2, 0, 1, 0))
  expect_equal(net_reachability(d, "all"), c(2, 2, 2, 1, 1, 0))
  expect_equal(net_reachability(pieces()), c(2, 2, 2, 1, 1, 0))
  expect_equal(unname(rowSums(is.finite(net_distances(d)))) - 1, net_reachability(d))
})

test_that("cut points and bridges", {
  g <- pieces()
  expect_equal(net_cut_points(g), 2L)
  expect_equal(net_bridges(g), 1:3)
  # Triangle a-b-c, bridge c-d, and d-e doubled (not a bridge).
  h <- net_frame(data.frame(source = c("a", "b", "c", "c", "d", "e"),
                           target   = c("b", "c", "a", "d", "e", "d")),
                directed = FALSE)
  expect_equal(net_cut_points(h), c("c", "d"))
  expect_equal(net_bridges(h), 4L)
  # Arcs in both directions are two edges.
  d <- net_frame(data.frame(source = c(1, 2, 2), target = c(2, 1, 3)))
  expect_equal(net_cut_points(d), 2)
  expect_equal(net_bridges(d), 3L)
  empty <- net_frame(data.frame(source = 1, target = 1), directed = FALSE)
  expect_equal(net_bridges(empty), integer(0))
  expect_length(net_cut_points(empty), 0L)
})

test_that("cut points match igraph", {
  skip_if_not_installed("igraph")
  e <- random_edges(60, 70, 11)
  g <- net_frame(e, nodes = data.frame(id = 1:60), directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(id = 1:60))
  expect_equal(net_cut_points(g), sort(as.integer(igraph::articulation_points(ig))))
  ours <- e[net_bridges(g), ]
  theirs <- igraph::ends(ig, igraph::bridges(ig), names = FALSE)
  expect_equal(sort(paste(pmin(ours$source, ours$target), pmax(ours$source, ours$target))),
               sort(paste(pmin(theirs[, 1], theirs[, 2]), pmax(theirs[, 1], theirs[, 2]))))
})

test_that("summary reports the descriptive table", {
  g <- pieces()
  s <- summary(g)
  expect_s3_class(s, "summary.net_frame")
  expect_equal(s$nodes, 6L)
  expect_equal(s$edges, 3L)
  expect_equal(s$isolates, 1)
  expect_equal(s$loops, 0)
  expect_equal(s$components, 3L)
  expect_equal(s$largest_component, 3L)
  expect_equal(s$density, net_density(g))
  expect_equal(s$transitivity, 0)
  # On the largest component, the path 1 - 2 - 3.
  expect_equal(s$mean_distance, 4 / 3)
  expect_equal(s$diameter, 2)
  expect_equal(s$degree, c(mean = 1, min = 0, max = 2))
  expect_equal(s$centralization, net_centralization(g, "degree"))
  expect_true(is.na(s$reciprocity))
  out <- capture.output(print(s))
  expect_true(any(grepl("^Transitivity", out)))
  expect_true(any(grepl("largest component \\(3 of 6 nodes\\)", out)))
  expect_false(any(grepl("Reciprocity", out)))

  skipped <- capture.output(print(summary(g, paths = FALSE)))
  expect_true(any(grepl("skipped", skipped)))
  expect_true(any(grepl("not computed", skipped)))
  expect_error(summary(g, paths = NA), "paths")

  d <- net_frame(data.frame(source = c("a", "b", "b", "a", "c"), target = c("b", "a", "c", "b", "c"),
                           w = c(1, 2, 3, 4, 5)), weight = "w", weight_type = "distance")
  s <- summary(d)
  expect_equal(s$loops, 1)
  expect_equal(s$multiple, 1)
  expect_equal(s$reciprocity, net_reciprocity(d))
  expect_equal(s$mean_distance, net_mean_distance(d))
  out <- capture.output(print(s))
  expect_true(any(grepl("^Reciprocity", out)))
  expect_true(any(grepl("yes, `w` as distances (1 to 5)", out, fixed = TRUE)))
  expect_true(any(grepl("weights as distances", out)))
  expect_equal(summary(d, weighted = FALSE)$mean_distance, net_mean_distance(d, weighted = FALSE))
})
