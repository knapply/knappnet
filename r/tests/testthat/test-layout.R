# Layouts are not unique, so these tests check properties: shape of the
# result, reproducibility, separation of components, and geometry that a
# correct layout must have. validation/layout_check.py compares stress
# against NetworkX.

# A triangle 1-2-3, a pendant 4 on node 3, a separate edge 5-6, and an
# isolate 7.
three_pieces <- function(directed = FALSE) {
  net_frame(data.frame(source = c(1, 2, 3, 3, 5), target = c(2, 3, 1, 4, 6)),
           nodes = data.frame(id = 1:7), directed = directed)
}

dist_xy <- function(xy, a, b) sqrt(sum((xy[a, ] - xy[b, ])^2))

test_that("every method returns finite coordinates aligned with the nodes", {
  g <- three_pieces()
  g$nodes$type <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
  for (method in c("auto", "stress", "sparse_stress", "fr", "circle", "star", "grid",
                   "tree", "twomode")) {
    xy <- net_layout(g, method = method, seed = 1)
    expect_true(is.matrix(xy) && is.double(xy), info = method)
    expect_equal(dim(xy), c(7L, 2L), info = method)
    expect_equal(colnames(xy), c("x", "y"), info = method)
    expect_true(all(is.finite(xy)), info = method)
  }
})

test_that("layouts are reproducible from the seed and independent of threads", {
  g <- three_pieces(directed = TRUE)
  for (method in c("stress", "sparse_stress", "fr")) {
    a <- net_layout(g, method, seed = 42, threads = 1)
    expect_identical(a, net_layout(g, method, seed = 42, threads = 2), info = method)
    expect_identical(a, net_layout(g, method, seed = 42, threads = 1), info = method)
  }
  set.seed(3)
  a <- net_layout(g, "fr")
  set.seed(3)
  expect_identical(a, net_layout(g, "fr"))
  expect_false(identical(net_layout(g, "fr", seed = 1), net_layout(g, "fr", seed = 2)))
  expect_error(net_layout(g, seed = -1), "seed")
  expect_error(net_layout(g, seed = 1.5), "seed")
})

test_that("stress layouts reproduce distances that can be drawn exactly", {
  path <- net_frame(data.frame(source = 1:4, target = 2:5), directed = FALSE)
  xy <- net_layout(path, seed = 1)
  for (i in 1:4) expect_equal(dist_xy(xy, i, i + 1), 1, tolerance = 1e-3)
  expect_equal(dist_xy(xy, 1, 5), 4, tolerance = 1e-3)
  expect_lt(net_layout_stress(path, xy), 1e-5)
  # The longest extent is horizontal, and node 1 is on the left.
  expect_lt(xy[1, "x"], xy[5, "x"])
  expect_true(all(abs(xy[, "y"]) < 1e-2))

  tri <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 1)), directed = FALSE)
  xy <- net_layout(tri, seed = 1)
  expect_equal(c(dist_xy(xy, 1, 2), dist_xy(xy, 2, 3), dist_xy(xy, 1, 3)), c(1, 1, 1),
               tolerance = 1e-6)
})

test_that("weights are lengths: distances as given, strengths as 1/weight", {
  e <- data.frame(source = c(1, 2, 3), target = c(2, 3, 1), w = c(2, 2, 2))
  far <- net_frame(e, directed = FALSE, weight = "w", weight_type = "distance")
  xy <- net_layout(far, seed = 1)
  expect_equal(dist_xy(xy, 1, 2), 2, tolerance = 1e-6)
  near <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  expect_equal(dist_xy(net_layout(near, seed = 1), 1, 2), 0.5, tolerance = 1e-6)
  expect_equal(dist_xy(net_layout(far, seed = 1, weighted = FALSE), 1, 2), 1,
               tolerance = 1e-6)
  expect_equal(dist_xy(net_layout(far, seed = 1, weights_as = "strength"), 1, 2), 0.5,
               tolerance = 1e-6)
  zero <- net_frame(data.frame(source = 1, target = 2, w = 0), weight = "w", weight_type = "distance")
  expect_error(net_layout(zero), "positive")
  expect_error(net_layout(three_pieces(), weighted = TRUE), "no weights")
})

test_that("components are packed apart and isolates are placed", {
  g <- three_pieces()
  for (method in c("stress", "fr", "tree")) {
    xy <- net_layout(g, method, seed = 7)
    comp <- net_components(g)$membership
    for (a in 1:7) for (b in 1:7) {
      if (comp[a] != comp[b]) {
        expect_gte(dist_xy(xy, a, b), 1 - 1e-9, label = paste(method, a, b))
      }
    }
  }
})

test_that("simple layouts", {
  g <- three_pieces()
  xy <- net_layout(g, "circle")
  expect_equal(sqrt(rowSums(xy^2)), rep(1, 7))
  expect_equal(unname(xy[1, ]), c(1, 0))
  star <- net_layout(g, "star", center = 3)
  expect_equal(unname(star[3, ]), c(0, 0))
  expect_equal(sqrt(rowSums(star[-3, ]^2)), rep(1, 6))
  # The default center is the node with the most neighbors (node 3).
  expect_identical(net_layout(g, "star"), star)
  expect_error(net_layout(g, "star", center = 99), "not in the network")
  grid <- net_layout(g, "grid")
  expect_equal(unname(grid[, "x"]), c(0, 1, 2, 0, 1, 2, 0))
  expect_equal(unname(grid[, "y"]), c(0, 0, 0, -1, -1, -1, -2))
})

test_that("tree layouts are layered", {
  tree <- net_frame(data.frame(source = c("r", "r", "a", "a"), target = c("a", "b", "c", "d")),
                   directed = TRUE)
  xy <- net_layout(tree, "tree", root = "r")
  expect_equal(unname(xy[, "y"] - xy[1, "y"]), c(0, -1, -1, -2, -2))
  # Parents are centered over their first and last child.
  expect_equal(unname(xy[2, "x"]), mean(xy[4:5, "x"]))
  expect_equal(unname(xy[1, "x"]), mean(xy[2:3, "x"]))
  expect_error(net_layout(tree, "tree", root = c("r", "a")), "more than one root")
  expect_error(net_layout(tree, "tree", root = "zz"), "not in the network")
})

test_that("two-mode layouts put the modes in two rows", {
  g <- net_frame(data.frame(source = c("p1", "p1", "p2", "p3"), target = c("e2", "e1", "e1", "e2")),
                directed = FALSE)
  g$nodes$kind <- ifelse(startsWith(g$nodes$name, "p"), "person", "event")
  xy <- net_layout(g, "twomode", types = "kind")
  # "event" sorts first, so events form the top row.
  events <- g$nodes$kind == "event"
  expect_true(min(xy[events, "y"]) > max(xy[!events, "y"]))
  expect_identical(net_layout(g, "twomode", types = g$nodes$kind), xy)
  expect_error(net_layout(g, "twomode"), "types")
  expect_error(net_layout(g, "twomode", types = "nope"), "not found")
  expect_error(net_layout(g, "twomode", types = seq_len(nrow(g$nodes))), "two distinct")
  g$nodes$type <- g$nodes$kind == "person"
  top <- net_layout(g, "twomode")
  expect_true(min(top[!g$nodes$type, "y"]) > max(top[g$nodes$type, "y"]))
})

test_that("arguments are checked", {
  g <- three_pieces()
  expect_error(net_layout(g, "nope"))
  expect_error(net_layout(g, max_iter = 0), "max_iter")
  expect_error(net_layout(g, pivots = -1), "pivots")
  expect_error(net_layout(g, tol = 0), "tol")
  expect_error(net_layout(data.frame()), "net_frame")
  expect_error(net_layout_stress(g, matrix(0, 2, 2)), "7 rows")
})

test_that("a network with one node or no edges", {
  one <- net_frame(data.frame(source = character(), target = character()),
                  nodes = data.frame(id = "a"))
  expect_equal(unname(net_layout(one)), matrix(0, 1, 2))
  empty <- net_frame(data.frame(source = 1:0, target = 1:0)[0, ], nodes = data.frame(id = 1:4))
  xy <- net_layout(empty, seed = 1)
  expect_true(all(is.finite(xy)))
  expect_gte(min(dist(xy)), 1 - 1e-9)
})

test_that("stress is no worse than igraph's Kamada-Kawai layout on a classic network", {
  skip_if_not_installed("igraph")
  karate <- as_net_frame(igraph::make_graph("Zachary"))
  ours <- net_layout(karate, seed = 1)
  # Kamada-Kawai minimizes a weighted stress of the same kind.
  set.seed(1)
  theirs <- igraph::layout_with_kk(igraph::make_graph("Zachary"))
  # Compare at the scale that minimizes each layout's stress.
  best <- function(xy) {
    stats::optimize(function(s) net_layout_stress(karate, xy * s), c(1e-3, 1e3))$objective
  }
  expect_lte(best(ours), best(theirs) * 1.02)
})

test_that("two-mode layouts use the node type of net_twomode() networks", {
  attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
                         event  = c("gala", "fair", "gala", "fair"))
  g <- net_twomode(attended)
  xy <- net_layout(g, "twomode")
  second <- g$nodes[[g$roles$type]]
  expect_true(min(xy[!second, "y"]) > max(xy[second, "y"]))
  # The type role is used even when the column is not called `type`.
  expect_identical(g$roles$type, "type")
  names(g$nodes)[names(g$nodes) == "type"] <- "is_event"
  g$roles$type <- "is_event"
  expect_identical(net_layout(g, "twomode"), xy)
})
