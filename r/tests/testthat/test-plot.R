# Plots are checked by behavior on a null device, not by pixels: they must
# draw without error or warning for every kind of input, return the layout
# they used, restore graphical parameters, and map columns as documented.

on_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

office <- function() {
  g <- netframe(
    data.frame(from = c("ann", "bob", "cat", "dan", "eve", "bob", "fay", "ann"),
               to   = c("bob", "cat", "dan", "eve", "cat", "fay", "dan", "cat"),
               calls = c(1, 5, 2, 8, 3, 1, 4, 2)),
    nodes = data.frame(name = c("ann", "bob", "cat", "dan", "eve", "fay", "gus"),
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
    expect_identical(xy, layout_network(g, seed = 1))
    undirected <- netframe(g$edges, nodes = g$nodes, directed = FALSE)
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
    expect_identical(plot(g, layout = "circle"), layout_network(g, "circle"))
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
    strong <- netframe(data.frame(from = c(1, 2, 3), to = c(2, 3, 1), w = c(1, 5, 10)),
                       directed = FALSE, weight = "w", weight_type = "strength")
    expect_silent(plot(strong, seed = 1))
    far <- netframe(data.frame(from = c(1, 2), to = c(2, 3), w = c(1, 5)), weight = "w")
    expect_silent(plot(far, seed = 1))
    empty <- netframe(data.frame(from = 1:0, to = 1:0)[0, ], nodes = data.frame(id = 1:3))
    expect_silent(plot(empty, seed = 1))
    none <- netframe(data.frame(from = character(), to = character()),
                     nodes = data.frame(id = character()))
    expect_identical(dim(plot(none)), c(0L, 2L))
    one <- netframe(data.frame(from = "a", to = "a"))
    expect_silent(plot(one, seed = 1))
    set.seed(1)
    big <- netframe(data.frame(from = sample(500, 1500, TRUE), to = sample(500, 1500, TRUE)),
                    nodes = data.frame(id = 1:500))
    expect_silent(xy <- plot(big, seed = 1))
    expect_true(all(is.finite(xy)))
  })
})

test_that("column mappings follow the documented palettes and order", {
  cols <- map_colors(c("b", "a", NA, "b"))
  expect_identical(cols$labels, c("a", "b", "NA"))
  pal <- grDevices::hcl.colors(2, "Dark 3")
  expect_identical(cols$col, c(pal[2], pal[1], "gray80", pal[2]))
  f <- map_colors(factor(c("lo", "hi"), levels = c("lo", "hi", "unused")))
  expect_identical(f$labels, c("lo", "hi"))
  num <- map_colors(c(0, 5, 10))
  viridis <- grDevices::hcl.colors(100, "viridis")
  expect_identical(num$col, viridis[c(1, 51, 100)])
  sizes <- map_numeric(c(0, 1, 4, NA), 1, 2, area = TRUE)$values
  expect_equal(sizes, c(1, 1.5, 2, 1))
})
