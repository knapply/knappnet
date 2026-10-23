# Layouts are not unique, so these tests check properties: shape of the
# result, reproducibility, separation of components, and geometry that a
# correct layout must have. validation/layout_check.py compares stress
# against NetworkX.

# A triangle 1-2-3, a pendant 4 on node 3, a separate edge 5-6, and an
# isolate 7.
pieces <- function(directed = FALSE) {
  netframe(data.frame(from = c(1, 2, 3, 3, 5), to = c(2, 3, 1, 4, 6)),
           nodes = data.frame(id = 1:7), directed = directed)
}

dist_xy <- function(xy, a, b) sqrt(sum((xy[a, ] - xy[b, ])^2))

test_that("every method returns finite coordinates aligned with the nodes", {
  g <- pieces()
  g$nodes$type <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
  for (method in c("auto", "stress", "sparse_stress", "fr", "circle", "star", "grid",
                   "tree", "twomode")) {
    xy <- layout_network(g, method = method, seed = 1)
    expect_true(is.matrix(xy) && is.double(xy), info = method)
    expect_equal(dim(xy), c(7L, 2L), info = method)
    expect_equal(colnames(xy), c("x", "y"), info = method)
    expect_true(all(is.finite(xy)), info = method)
  }
})

test_that("layouts are reproducible from the seed and independent of threads", {
  g <- pieces(directed = TRUE)
  for (method in c("stress", "sparse_stress", "fr")) {
    a <- layout_network(g, method, seed = 42, threads = 1)
    expect_identical(a, layout_network(g, method, seed = 42, threads = 2), info = method)
    expect_identical(a, layout_network(g, method, seed = 42, threads = 1), info = method)
  }
  set.seed(3)
  a <- layout_network(g, "fr")
  set.seed(3)
  expect_identical(a, layout_network(g, "fr"))
  expect_false(identical(layout_network(g, "fr", seed = 1), layout_network(g, "fr", seed = 2)))
  expect_error(layout_network(g, seed = -1), "seed")
  expect_error(layout_network(g, seed = 1.5), "seed")
})

test_that("stress layouts reproduce distances that can be drawn exactly", {
  path <- netframe(data.frame(from = 1:4, to = 2:5), directed = FALSE)
  xy <- layout_network(path, seed = 1)
  for (i in 1:4) expect_equal(dist_xy(xy, i, i + 1), 1, tolerance = 1e-3)
  expect_equal(dist_xy(xy, 1, 5), 4, tolerance = 1e-3)
  expect_lt(layout_stress(path, xy), 1e-5)
  # The longest extent is horizontal, and node 1 is on the left.
  expect_lt(xy[1, "x"], xy[5, "x"])
  expect_true(all(abs(xy[, "y"]) < 1e-2))

  tri <- netframe(data.frame(from = c(1, 2, 3), to = c(2, 3, 1)), directed = FALSE)
  xy <- layout_network(tri, seed = 1)
  expect_equal(c(dist_xy(xy, 1, 2), dist_xy(xy, 2, 3), dist_xy(xy, 1, 3)), c(1, 1, 1),
               tolerance = 1e-6)
})

test_that("weights are lengths: distances as given, strengths as 1/weight", {
  e <- data.frame(from = c(1, 2, 3), to = c(2, 3, 1), w = c(2, 2, 2))
  far <- netframe(e, directed = FALSE, weight = "w")
  xy <- layout_network(far, seed = 1)
  expect_equal(dist_xy(xy, 1, 2), 2, tolerance = 1e-6)
  near <- netframe(e, directed = FALSE, weight = "w", weight_type = "strength")
  expect_equal(dist_xy(layout_network(near, seed = 1), 1, 2), 0.5, tolerance = 1e-6)
  expect_equal(dist_xy(layout_network(far, seed = 1, weighted = FALSE), 1, 2), 1,
               tolerance = 1e-6)
  expect_equal(dist_xy(layout_network(far, seed = 1, weights_as = "strength"), 1, 2), 0.5,
               tolerance = 1e-6)
  zero <- netframe(data.frame(from = 1, to = 2, w = 0), weight = "w")
  expect_error(layout_network(zero), "positive")
  expect_error(layout_network(pieces(), weighted = TRUE), "no weights")
})

test_that("components are packed apart and isolates are placed", {
  g <- pieces()
  for (method in c("stress", "fr", "tree")) {
    xy <- layout_network(g, method, seed = 7)
    comp <- components(g)$membership
    for (a in 1:7) for (b in 1:7) {
      if (comp[a] != comp[b]) {
        expect_gte(dist_xy(xy, a, b), 1 - 1e-9, label = paste(method, a, b))
      }
    }
  }
})

test_that("simple layouts", {
  g <- pieces()
  xy <- layout_network(g, "circle")
  expect_equal(sqrt(rowSums(xy^2)), rep(1, 7))
  expect_equal(unname(xy[1, ]), c(1, 0))
  star <- layout_network(g, "star", center = 3)
  expect_equal(unname(star[3, ]), c(0, 0))
  expect_equal(sqrt(rowSums(star[-3, ]^2)), rep(1, 6))
  # The default center is the node with the most neighbors (node 3).
  expect_identical(layout_network(g, "star"), star)
  expect_error(layout_network(g, "star", center = 99), "not nodes")
  grid <- layout_network(g, "grid")
  expect_equal(unname(grid[, "x"]), c(0, 1, 2, 0, 1, 2, 0))
  expect_equal(unname(grid[, "y"]), c(0, 0, 0, -1, -1, -1, -2))
})

test_that("tree layouts are layered", {
  tree <- netframe(data.frame(from = c("r", "r", "a", "a"), to = c("a", "b", "c", "d")),
                   directed = TRUE)
  xy <- layout_network(tree, "tree", root = "r")
  expect_equal(unname(xy[, "y"] - xy[1, "y"]), c(0, -1, -1, -2, -2))
  # Parents are centered over their first and last child.
  expect_equal(unname(xy[2, "x"]), mean(xy[4:5, "x"]))
  expect_equal(unname(xy[1, "x"]), mean(xy[2:3, "x"]))
  expect_error(layout_network(tree, "tree", root = c("r", "a")), "more than one root")
  expect_error(layout_network(tree, "tree", root = "zz"), "not nodes")
})

test_that("two-mode layouts put the modes in two rows", {
  g <- netframe(data.frame(from = c("p1", "p1", "p2", "p3"), to = c("e2", "e1", "e1", "e2")),
                directed = FALSE)
  g$nodes$kind <- ifelse(startsWith(g$nodes$name, "p"), "person", "event")
  xy <- layout_network(g, "twomode", types = "kind")
  # "event" sorts first, so events form the top row.
  events <- g$nodes$kind == "event"
  expect_true(min(xy[events, "y"]) > max(xy[!events, "y"]))
  expect_identical(layout_network(g, "twomode", types = g$nodes$kind), xy)
  expect_error(layout_network(g, "twomode"), "types")
  expect_error(layout_network(g, "twomode", types = "nope"), "not found")
  expect_error(layout_network(g, "twomode", types = seq_len(nrow(g$nodes))), "two distinct")
  g$nodes$type <- g$nodes$kind == "person"
  top <- layout_network(g, "twomode")
  expect_true(min(top[!g$nodes$type, "y"]) > max(top[g$nodes$type, "y"]))
})

test_that("arguments are checked", {
  g <- pieces()
  expect_error(layout_network(g, "nope"))
  expect_error(layout_network(g, max_iter = 0), "max_iter")
  expect_error(layout_network(g, pivots = -1), "pivots")
  expect_error(layout_network(g, tol = 0), "tol")
  expect_error(layout_network(data.frame()), "netframe")
  expect_error(layout_stress(g, matrix(0, 2, 2)), "7 rows")
})

test_that("a network with one node or no edges", {
  one <- netframe(data.frame(from = character(), to = character()),
                  nodes = data.frame(id = "a"))
  expect_equal(unname(layout_network(one)), matrix(0, 1, 2))
  empty <- netframe(data.frame(from = 1:0, to = 1:0)[0, ], nodes = data.frame(id = 1:4))
  xy <- layout_network(empty, seed = 1)
  expect_true(all(is.finite(xy)))
  expect_gte(min(dist(xy)), 1 - 1e-9)
})

test_that("stress is no worse than igraph's stress layout on a classic network", {
  skip_if_not_installed("igraph")
  karate <- as_netframe(igraph::make_graph("Zachary"))
  ours <- layout_network(karate, seed = 1)
  theirs <- igraph::layout_with_stress(igraph::make_graph("Zachary"))
  # Compare at the scale that minimizes each layout's stress.
  best <- function(xy) {
    stats::optimize(function(s) layout_stress(karate, xy * s), c(1e-3, 1e3))$objective
  }
  expect_lte(best(ours), best(theirs) * 1.02)
})

test_that("two-mode layouts use the node type of twomode() networks", {
  attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
                         event  = c("gala", "fair", "gala", "fair"))
  g <- twomode(attended)
  xy <- layout_network(g, "twomode")
  second <- g$nodes[[g$roles$type]]
  expect_true(min(xy[!second, "y"]) > max(xy[second, "y"]))
  # The type role is used even when the column is not called `type`.
  expect_identical(g$roles$type, "type")
  names(g$nodes)[names(g$nodes) == "type"] <- "is_event"
  g$roles$type <- "is_event"
  expect_identical(layout_network(g, "twomode"), xy)
})
