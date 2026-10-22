# Reciprocity, transitivity and clustering, distances, mean distance,
# diameter, eccentricity, radius, reachability, cut points, bridges, and
# summary(). Hand-checked values first, then spot checks against igraph and
# sna. The full comparison lives in validation/.

# A path 1 - 2 - 3, an edge 4 - 5, and an isolate 6.
pieces <- function() {
  netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)),
           nodes = data.frame(id = 1:6), directed = FALSE)
}

kite <- function() {
  e <- data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
  netframe(e, nodes = data.frame(id = 1:10), directed = FALSE)
}

random_edges <- function(n, m, seed) {
  set.seed(seed)
  e <- data.frame(from = sample(n, m, TRUE), to = sample(n, m, TRUE))
  e <- e[e$from != e$to, ]
  e[!duplicated(data.frame(a = pmin(e$from, e$to), b = pmax(e$from, e$to))), ]
}

test_that("reciprocity indices count dyads", {
  # a <-> b, b -> c, a loop at c and a repeated a -> b (both ignored), and
  # an isolate d: M = 1, A = 1, N = 4.
  g <- netframe(data.frame(from = c("a", "b", "b", "c", "a"),
                           to   = c("b", "a", "c", "c", "b")),
                nodes = data.frame(name = c("a", "b", "c", "d")))
  expect_equal(reciprocity(g), 2 / 3)
  expect_equal(reciprocity(g, "ratio"), 1 / 2)
  expect_equal(reciprocity(g, "dyadic"), 5 / 6)
  # p = 3 / 12 and r = 2 / 3, so rho = (2/3 - 1/4) / (3/4).
  expect_equal(reciprocity(g, "correlation"), 5 / 9)
  empty <- netframe(data.frame(from = character(0), to = character(0)),
                    nodes = data.frame(name = c("a", "b")))
  expect_true(is.nan(reciprocity(empty)))
  expect_equal(reciprocity(empty, "dyadic"), 1)
  expect_equal(reciprocity(pieces()), 1)
  expect_error(reciprocity(g, "edgewise"))
})

test_that("reciprocity matches igraph and sna", {
  skip_if_not_installed("igraph")
  set.seed(3)
  e <- data.frame(from = sample(25, 120, TRUE), to = sample(25, 120, TRUE))
  g <- netframe(e, nodes = data.frame(id = 1:25))
  ig <- igraph::simplify(igraph::graph_from_data_frame(e, vertices = data.frame(id = 1:25)))
  expect_equal(reciprocity(g), igraph::reciprocity(ig))
  expect_equal(reciprocity(g, "ratio"), igraph::reciprocity(ig, mode = "ratio"))
  skip_if_not_installed("sna")
  m <- matrix(0, 25, 25)
  m[cbind(e$from, e$to)] <- 1
  diag(m) <- 0
  expect_equal(reciprocity(g), unname(sna::grecip(m, measure = "edgewise")))
  expect_equal(reciprocity(g, "ratio"), unname(sna::grecip(m, measure = "dyadic.nonnull")))
  expect_equal(reciprocity(g, "dyadic"), unname(sna::grecip(m, measure = "dyadic")))
  pairs <- which(row(m) != col(m))
  expect_equal(reciprocity(g, "correlation"), stats::cor(m[pairs], t(m)[pairs]))
})

test_that("transitivity and clustering by hand", {
  g <- pieces()
  # One connected triple (centred on 2), not closed.
  expect_equal(transitivity(g), 0)
  expect_equal(clustering(g), c(NaN, 0, NaN, NaN, NaN, NaN))
  expect_equal(clustering(g, convention = "networkx"), rep(0, 6))
  expect_equal(average_clustering(g), 0)
  # A triangle with a pendant node and an isolate.
  t <- netframe(data.frame(from = c(1, 2, 3, 3), to = c(2, 3, 1, 4)),
                nodes = data.frame(id = 1:5), directed = FALSE)
  expect_equal(transitivity(t), 3 / 5)
  expect_equal(clustering(t), c(1, 1, 1 / 3, NaN, NaN))
  expect_equal(average_clustering(t), (2 + 1 / 3) / 3)
  expect_equal(average_clustering(t, convention = "networkx"), (2 + 1 / 3) / 5)
  # Direction, parallel edges, and loops do not matter.
  d <- netframe(data.frame(from = c(1, 2, 2, 3, 3, 3, 1), to = c(2, 1, 3, 1, 3, 4, 2)),
                nodes = data.frame(id = 1:4))
  expect_equal(clustering(d), c(1, 1, 1 / 3, NaN))
  no_triples <- netframe(data.frame(from = 1, to = 2), directed = FALSE)
  expect_true(is.nan(transitivity(no_triples)))
  expect_equal(transitivity(no_triples, convention = "networkx"), 0)
})

test_that("weighted clustering by hand", {
  # Triangle 1-2-3 with weights 1, 2, 3 and a pendant 1-4 of weight 4.
  g <- netframe(data.frame(from = c(1, 2, 3, 1), to = c(2, 3, 1, 4), w = c(1, 2, 3, 4)),
                directed = FALSE, weight = "w")
  # Node 1: strength 8, three neighbours, one triangle with w12 + w13 = 4.
  expect_equal(clustering(g, weighted = TRUE), c(4 / 16, 1, 1, NaN))
  tri <- (1 * 2 * 3 / 4^3)^(1 / 3)
  expect_equal(clustering(g, weighted = TRUE, convention = "networkx"),
               c(2 * tri / 6, tri, tri, 0))
  expect_equal(clustering(g, weighted = TRUE, method = "onnela"),
               c(2 * tri / 6, tri, tri, NaN))
  expect_equal(average_clustering(g, weighted = TRUE, convention = "networkx"),
               (2 * tri / 6 + 2 * tri) / 4)
  # Binary by default, whatever the weights.
  expect_equal(clustering(g), c(1 / 3, 1, 1, NaN))
  expect_error(clustering(pieces(), weighted = TRUE), "no weights")
  bad <- netframe(data.frame(from = 1, to = 2, w = 0), weight = "w", directed = FALSE)
  expect_error(clustering(bad, weighted = TRUE), "positive")
})

test_that("clustering matches igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(transitivity(g), igraph::transitivity(ig, type = "global"))
  expect_equal(clustering(g), igraph::transitivity(ig, type = "local", isolates = "NaN"))
  expect_equal(average_clustering(g), igraph::transitivity(ig, type = "average"))
  e <- random_edges(30, 90, 5)
  e$w <- round(stats::runif(nrow(e), 0.5, 4), 2)
  wg <- netframe(e, nodes = data.frame(id = 1:30), directed = FALSE, weight = "w")
  wig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(id = 1:30))
  expect_equal(clustering(wg, weighted = TRUE),
               igraph::transitivity(wig, type = "barrat", weights = igraph::E(wig)$w,
                                    isolates = "NaN"))
})

test_that("distances and their summaries by hand", {
  g <- pieces()
  d <- distances(g)
  expect_equal(dimnames(d), list(as.character(1:6), as.character(1:6)))
  expect_equal(unname(d[1, ]), c(0, 1, 2, Inf, Inf, Inf))
  expect_equal(unname(d[6, ]), c(Inf, Inf, Inf, Inf, Inf, 0))
  expect_equal(distances(g, from = 3, to = c(1, 5)),
               matrix(c(2, Inf), 1, dimnames = list("3", c("1", "5"))))
  expect_error(distances(g, from = 7), "not in the network")
  # Ordered pairs joined by a path: 6 in the path (total 8) and 2 in the edge.
  expect_equal(mean_distance(g), 10 / 8)
  expect_equal(mean_distance(g, unconnected = FALSE), Inf)
  expect_equal(diameter(g), 2)
  expect_equal(diameter(g, unconnected = FALSE), Inf)
  expect_equal(eccentricity(g), c(2, 1, 2, 1, 1, 0))
  expect_equal(radius(g), 0)
  expect_error(mean_distance(g, unconnected = NA), "unconnected")

  # Directed, weighted: a -> b -> c -> d with minutes 5, 5, 20.
  w <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"),
                           minutes = c(5, 5, 20)), weight = "minutes")
  expect_equal(unname(distances(w)[1, ]), c(0, 5, 10, 30))
  expect_equal(unname(distances(w, mode = "in")[4, ]), c(30, 25, 20, 0))
  expect_equal(unname(distances(w, weighted = FALSE, mode = "all")[4, ]), c(3, 2, 1, 0))
  expect_equal(mean_distance(w), (5 + 10 + 30 + 5 + 25 + 20) / 6)
  expect_equal(diameter(w), 30)
  expect_equal(diameter(w, weighted = FALSE), 3)
  expect_equal(eccentricity(w, mode = "out", weighted = FALSE), c(3, 2, 1, 0))
  # As strengths, lengths are 1/5, 1/5, and 1/20.
  expect_equal(diameter(w, weights_as = "strength"), 0.2 + 0.2 + 0.05)
})

test_that("distances match igraph and sna", {
  skip_if_not_installed("igraph")
  set.seed(9)
  e <- data.frame(from = sample(30, 70, TRUE), to = sample(30, 70, TRUE))
  e$w <- sample(1:6, nrow(e), TRUE)
  e$r <- round(stats::runif(nrow(e), 0.1, 3), 3)
  nodes <- data.frame(id = 1:30)
  ig <- igraph::graph_from_data_frame(e, vertices = nodes)
  g <- netframe(e, nodes = nodes)
  for (mode in c("out", "in", "all")) {
    expect_equal(unname(distances(g, mode = mode)),
                 unname(igraph::distances(ig, mode = mode, weights = NA)))
  }
  gw <- netframe(e, nodes = nodes, weight = "w")
  expect_equal(unname(distances(gw)),
               unname(igraph::distances(ig, mode = "out", weights = e$w)))
  gr <- netframe(e, nodes = nodes, weight = "r")
  expect_equal(unname(distances(gr)),
               unname(igraph::distances(ig, mode = "out", weights = e$r)))
  expect_equal(mean_distance(g), igraph::mean_distance(ig))
  expect_equal(mean_distance(gr), igraph::mean_distance(ig, weights = e$r))
  expect_equal(diameter(g), igraph::diameter(ig))
  expect_equal(diameter(gr), igraph::diameter(ig, weights = e$r))
  expect_equal(diameter(g, unconnected = FALSE),
               igraph::diameter(ig, unconnected = FALSE))
  iu <- igraph::delete_edge_attr(igraph::delete_edge_attr(ig, "w"), "r")
  expect_equal(eccentricity(g), unname(igraph::eccentricity(iu, mode = "all")))
  expect_equal(eccentricity(g, mode = "out"), unname(igraph::eccentricity(iu, mode = "out")))
  expect_equal(radius(g), igraph::radius(iu, mode = "all"))
  skip_if_not_installed("sna")
  m <- matrix(0, 30, 30)
  m[cbind(e$from, e$to)] <- 1
  diag(m) <- 0
  expect_equal(unname(distances(g)), sna::geodist(m)$gdist)
})

test_that("reachability counts", {
  d <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)),
                nodes = data.frame(id = 1:6))
  expect_equal(reachability(d), c(2, 1, 0, 1, 0, 0))
  expect_equal(reachability(d, "in"), c(0, 1, 2, 0, 1, 0))
  expect_equal(reachability(d, "all"), c(2, 2, 2, 1, 1, 0))
  expect_equal(reachability(pieces()), c(2, 2, 2, 1, 1, 0))
  expect_equal(unname(rowSums(is.finite(distances(d)))) - 1, reachability(d))
})

test_that("cut points and bridges", {
  g <- pieces()
  expect_equal(cut_points(g), 2L)
  expect_equal(bridges(g), 1:3)
  # Triangle a-b-c, bridge c-d, and d-e doubled (not a bridge).
  h <- netframe(data.frame(from = c("a", "b", "c", "c", "d", "e"),
                           to   = c("b", "c", "a", "d", "e", "d")),
                directed = FALSE)
  expect_equal(cut_points(h), c("c", "d"))
  expect_equal(bridges(h), 4L)
  # Arcs in both directions are two edges.
  d <- netframe(data.frame(from = c(1, 2, 2), to = c(2, 1, 3)))
  expect_equal(cut_points(d), 2)
  expect_equal(bridges(d), 3L)
  empty <- netframe(data.frame(from = 1, to = 1), directed = FALSE)
  expect_equal(bridges(empty), integer(0))
  expect_length(cut_points(empty), 0L)
})

test_that("cut points match igraph", {
  skip_if_not_installed("igraph")
  e <- random_edges(60, 70, 11)
  g <- netframe(e, nodes = data.frame(id = 1:60), directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(id = 1:60))
  expect_equal(cut_points(g), sort(as.integer(igraph::articulation_points(ig))))
  ours <- e[bridges(g), ]
  theirs <- igraph::ends(ig, igraph::bridges(ig), names = FALSE)
  expect_equal(sort(paste(pmin(ours$from, ours$to), pmax(ours$from, ours$to))),
               sort(paste(pmin(theirs[, 1], theirs[, 2]), pmax(theirs[, 1], theirs[, 2]))))
})

test_that("summary reports the descriptive table", {
  g <- pieces()
  s <- summary(g)
  expect_s3_class(s, "summary.netframe")
  expect_equal(s$nodes, 6L)
  expect_equal(s$edges, 3L)
  expect_equal(s$isolates, 1)
  expect_equal(s$loops, 0)
  expect_equal(s$components, 3L)
  expect_equal(s$largest_component, 3L)
  expect_equal(s$density, graph_density(g))
  expect_equal(s$transitivity, 0)
  # On the largest component, the path 1 - 2 - 3.
  expect_equal(s$mean_distance, 4 / 3)
  expect_equal(s$diameter, 2)
  expect_equal(s$degree, c(mean = 1, min = 0, max = 2))
  expect_equal(s$centralization, centralization(g, "degree"))
  expect_true(is.na(s$reciprocity))
  out <- capture.output(print(s))
  expect_true(any(grepl("^Transitivity", out)))
  expect_true(any(grepl("largest component \\(3 of 6 nodes\\)", out)))
  expect_false(any(grepl("Reciprocity", out)))

  skipped <- capture.output(print(summary(g, paths = FALSE)))
  expect_true(any(grepl("skipped", skipped)))
  expect_true(any(grepl("not computed", skipped)))
  expect_error(summary(g, paths = NA), "paths")

  d <- netframe(data.frame(from = c("a", "b", "b", "a", "c"), to = c("b", "a", "c", "b", "c"),
                           w = c(1, 2, 3, 4, 5)), weight = "w")
  s <- summary(d)
  expect_equal(s$loops, 1)
  expect_equal(s$multiple, 1)
  expect_equal(s$reciprocity, reciprocity(d))
  expect_equal(s$mean_distance, mean_distance(d))
  out <- capture.output(print(s))
  expect_true(any(grepl("^Reciprocity", out)))
  expect_true(any(grepl("yes, by `w` \\(distances\\)", out)))
  expect_true(any(grepl("weights as distances", out)))
  expect_equal(summary(d, weighted = FALSE)$mean_distance, mean_distance(d, weighted = FALSE))
})
