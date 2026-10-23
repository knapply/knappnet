# plot.net_frame() (R/plot-net-frame.R): it must draw without error or
# warning for every kind of input, return the layout it used, restore
# graphical parameters, and map columns as documented.

office <- function() {
  g <- net_frame(
    data.frame(source = c("ann", "bob", "cat", "dan", "eve", "bob", "fay", "ann"),
               target   = c("bob", "cat", "dan", "eve", "cat", "fay", "dan", "cat"),
               calls = c(1, 5, 2, 8, 3, 1, 4, 2)),
    nodes = data.frame(id = c("ann", "bob", "cat", "dan", "eve", "fay", "gus"),
                       dept = c("sales", "sales", "eng", "eng", "eng", "legal", NA),
                       tenure = c(1, 4, 2, 9, 3, NA, 5))
  )
  g$edges$kind <- rep(c("advice", "friend"), 4)
  g
}

test_that("the default plot draws and returns the layout it used", {
  g <- office()
  on_null_device({
    expect_silent(xy <- plot(g, seed = 1))
    expect_identical(xy, net_layout(g, seed = 1))
    undirected <- net_frame(g$edges, nodes = g$nodes, directed = FALSE)
    expect_silent(plot(undirected, seed = 1, main = "Office"))
  })
})

test_that("plots are reproducible from the seed", {
  g <- office()
  on_null_device({
    a <- plot(g, layout = "fr", seed = 5)
    b <- plot(g, layout = "fr", seed = 5)
    expect_identical(a, b)
    set.seed(9)
    c1 <- plot(g, layout = "fr")
    set.seed(9)
    c2 <- plot(g, layout = "fr")
    expect_identical(c1, c2)
  })
})

test_that("graphical parameters are restored", {
  g <- office()
  on_null_device({
    before <- graphics::par(c("mar", "xpd"))
    plot(g, node_color = "dept", seed = 1)
    expect_identical(graphics::par(c("mar", "xpd")), before)
  })
})

test_that("columns, vectors, constants, and I() values are accepted", {
  g <- office()
  n <- nrow(g$nodes)
  m <- nrow(g$edges)
  on_null_device({
    expect_silent(plot(g, node_color = "dept", node_size = "tenure", edge_color = "kind",
                       edge_width = "calls", seed = 1))
    expect_silent(plot(g, node_color = "tenure", seed = 1))
    expect_silent(plot(g, node_color = factor(g$nodes$dept), seed = 1))
    expect_silent(plot(g, node_color = g$nodes$tenure > 2, seed = 1))
    expect_silent(plot(g, node_color = "tomato", edge_color = "gray20", seed = 1))
    expect_silent(plot(g, node_color = I(rainbow(n)), node_size = I(rep(2, n)),
                       edge_width = I(seq_len(m)), seed = 1))
    expect_silent(plot(g, node_size = 3, edge_width = 2, seed = 1))
    expect_silent(plot(g, node_label = "dept", seed = 1))
    expect_silent(plot(g, node_label = FALSE, legend = FALSE, node_color = "dept", seed = 1))
    expect_silent(plot(g, node_label = TRUE, seed = 1))
  })
})

test_that("layouts can be given in every documented form", {
  g <- office()
  n <- nrow(g$nodes)
  g$nodes$lx <- seq_len(n)
  g$nodes$ly <- rev(seq_len(n))
  fixed <- cbind(x = seq_len(n), y = seq_len(n)^2)
  on_null_device({
    expect_equal(unname(plot(g, layout = c("lx", "ly"))), unname(cbind(g$nodes$lx, g$nodes$ly)))
    expect_equal(unname(plot(g, layout = fixed)), unname(fixed))
    expect_equal(unname(plot(g, layout = as.data.frame(fixed))), unname(fixed))
    expect_equal(unname(plot(g, layout = function(x) fixed)), unname(fixed))
    expect_identical(plot(g, layout = "circle"), net_layout(g, "circle"))
    expect_error(plot(g, layout = "spiral"), "Unknown layout")
    expect_error(plot(g, layout = c("lx", "zz")), "not found")
    expect_error(plot(g, layout = fixed[-1, ]), "rows")
    expect_error(plot(g, layout = fixed * NA), "finite")
  })
})

test_that("bad appearance arguments give clear errors", {
  g <- office()
  on_null_device({
    expect_error(plot(g, node_color = "nope", seed = 1), "neither a column")
    expect_error(plot(g, node_size = "dept", seed = 1), "numeric")
    expect_error(plot(g, edge_width = "kind", seed = 1), "numeric")
    expect_error(plot(g, node_color = c("red", "blue"), seed = 1), "one value per node")
    expect_error(plot(g, node_size = I(c(1, 2)), seed = 1), "one value per node")
    expect_error(plot(g, node_size = -1, seed = 1), "not negative")
  })
})

test_that("weighted, empty, tiny, and larger networks plot", {
  on_null_device({
    strong <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 1), w = c(1, 5, 10)),
                       directed = FALSE, weight = "w", weight_type = "strength")
    expect_silent(plot(strong, seed = 1))
    far <- net_frame(data.frame(source = c(1, 2), target = c(2, 3), w = c(1, 5)), weight = "w", weight_type = "distance")
    expect_silent(plot(far, seed = 1))
    empty <- net_frame(data.frame(source = 1:0, target = 1:0)[0, ], nodes = data.frame(id = 1:3))
    expect_silent(plot(empty, seed = 1))
    none <- net_frame(data.frame(source = character(), target = character()),
                     nodes = data.frame(id = character()))
    expect_identical(dim(plot(none)), c(0L, 2L))
    one <- net_frame(data.frame(source = "a", target = "a"))
    expect_silent(plot(one, seed = 1))
    set.seed(1)
    big <- net_frame(data.frame(source = sample(500, 1500, TRUE), target = sample(500, 1500, TRUE)),
                    nodes = data.frame(id = 1:500))
    expect_silent(xy <- plot(big, seed = 1))
    expect_true(all(is.finite(xy)))
  })
})

test_that("column mappings follow the documented palettes and order", {
  cols <- map_colors(c("b", "a", NA, "b"))
  expect_identical(cols$labels, c("a", "b", "NA"))
  pal <- plot_qualitative(2)
  expect_identical(cols$col, c(pal[2], pal[1], "gray80", pal[2]))
  f <- map_colors(factor(c("lo", "hi"), levels = c("lo", "hi", "unused")))
  expect_identical(f$labels, c("lo", "hi"))
  num <- map_colors(c(0, 5, 10))
  viridis <- plot_sequential(100)
  expect_identical(num$col, viridis[c(1, 51, 100)])
  sizes <- map_numeric(c(0, 1, 4, NA), 1, 2, area = TRUE)$values
  expect_equal(sizes, c(1, 1.5, 2, 1))
})

test_that("two-mode networks are colored by mode by default", {
  attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
                         event  = c("gala", "fair", "gala", "fair"))
  g <- net_twomode(attended)
  fill <- twomode_fill(g)
  pal <- plot_qualitative(2)
  expect_identical(fill$col, pal[g$nodes[[g$roles$type]] + 1L])
  expect_identical(fill$legend$labels, c("first", "second"))
  nodes <- data.frame(id = c("ann", "bob", "gala"), kind = c("person", "person", "event"))
  labelled <- net_twomode(data.frame(source = c("ann", "bob"), target = c("gala", "gala")), nodes,
                      type = "kind", first = "person")
  expect_identical(twomode_fill(labelled)$legend$labels, c("person", "event"))
  on_null_device({
    expect_silent(plot(g, seed = 1))
    expect_silent(plot(g, layout = "twomode"))
    expect_silent(plot(g, node_color = "red", seed = 1))
  })
})

test_that("spatial networks are drawn at their locations by default", {
  nodes <- data.frame(id = c("London", "Paris", "Berlin", "Madrid"),
                      lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                      lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                      region = c("north", "north", "north", "south"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(source = c("London", "Paris", "Paris"), target = c("Paris", "Berlin", "Madrid"),
                      stringsAsFactors = FALSE)
  plain <- net_frame(edges, nodes, directed = FALSE)
  g <- net_spatial(plain, coords = c("lon", "lat"), radius = 6371.0088)
  on_null_device({
    # The default and "geo" hand over to net_plot_spatial(), which returns the
    # network; a named layout method still draws an abstract layout.
    expect_silent(out <- plot(g))
    expect_identical(out, g)
    expect_silent(out <- plot(g, layout = "geo", node_color = "region", main = "Cities"))
    expect_identical(out, g)
    expect_silent(xy <- plot(g, layout = "circle"))
    expect_identical(xy, net_layout(g, "circle"))
    expect_error(plot(plain, layout = "geo"), "node locations")
    expect_silent(xy <- plot(plain, seed = 1))
    expect_identical(dim(xy), c(4L, 2L))
  })
})
