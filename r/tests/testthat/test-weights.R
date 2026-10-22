# A triangle where the direct edge a-c has weight 4 and the path through b
# has two edges of weight 1. Read as distances, the path through b is
# shorter; read as strengths (lengths 1/w), the direct tie is shorter.
triangle <- function(weight_type = "distance") {
  netframe(data.frame(from = c("a", "b", "a"), to = c("b", "c", "c"), w = c(1, 1, 4)),
           directed = FALSE, weight = "w", weight_type = weight_type)
}

test_that("binary networks are reported and analysed as unweighted", {
  g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")))
  expect_null(g$roles$weight)
  expect_output(print(g), "unweighted")
  expect_equal(betweenness(g), c(0, 1, 0))
  expect_error(strength(g), "no weights")
})

test_that("weighted networks use their weights by default", {
  g <- triangle()
  expect_output(print(g), "weighted by `w` \\(distances\\)")
  expect_equal(betweenness(g), c(0, 1, 0))
  expect_equal(betweenness(g, weighted = FALSE), c(0, 0, 0))
})

test_that("strength weights are inverted for paths but summed as given", {
  g <- triangle("strength")
  expect_equal(betweenness(g), c(0, 0, 0))
  expect_equal(betweenness(triangle(), weights_as = "strength"), c(0, 0, 0))
  expect_equal(strength(g), c(5, 2, 5))
  expect_equal(degree(g), c(2, 2, 2))
})

test_that("every weighting matches igraph", {
  skip_if_not_installed("igraph")
  set.seed(42)
  e <- data.frame(from = sample(40, 200, TRUE), to = sample(40, 200, TRUE))
  e <- e[e$from != e$to, ]
  e$int_w <- sample(1:5, nrow(e), TRUE)
  e$real_w <- round(runif(nrow(e), 0.1, 3), 3)
  nodes <- data.frame(id = 1:40)
  for (col in c("int_w", "real_w")) {
    g <- netframe(e, nodes, weight = col)
    ig <- igraph::graph_from_data_frame(e, vertices = nodes)
    w <- e[[col]]
    expect_equal(betweenness(g), unname(igraph::betweenness(ig, weights = w)))
    expect_equal(betweenness(g, weights_as = "strength"),
                 unname(igraph::betweenness(ig, weights = 1 / w)))
    expect_equal(strength(g, "out"), unname(igraph::strength(ig, mode = "out", weights = w)))
  }
})
