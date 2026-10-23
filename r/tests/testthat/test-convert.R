sample_net <- function(directed = TRUE) {
  nodes <- data.frame(name = c("ann", "bob", "cy", "dee"),
                      age = c(31L, 45L, NA, 22L),
                      member = c(TRUE, FALSE, TRUE, NA),
                      stringsAsFactors = FALSE)
  edges <- data.frame(source = c("ann", "bob", "cy", "ann"),
                      target = c("bob", "cy", "cy", "bob"),
                      weight = c(1.5, 2, 0.5, 3),
                      kind = c("a", "b", "c", "d"),
                      stringsAsFactors = FALSE)
  net_frame(edges, nodes, directed = directed, weight = "weight", weight_type = "distance",
           graph = list(title = "sample"))
}

expect_same_network <- function(a, b) {
  expect_equal(as.character(a$nodes[[a$roles$id]]), as.character(b$nodes[[b$roles$id]]))
  expect_equal(a$directed, b$directed)
  expect_equal(as.character(a$edges[[a$roles$source]]), as.character(b$edges[[b$roles$source]]))
  expect_equal(as.character(a$edges[[a$roles$target]]), as.character(b$edges[[b$roles$target]]))
}

test_that("igraph round trip keeps structure and attributes", {
  skip_if_not_installed("igraph")
  for (directed in c(TRUE, FALSE)) {
    g <- sample_net(directed)
    ig <- net_to_igraph(g)
    expect_equal(igraph::vcount(ig), 4)
    expect_equal(igraph::ecount(ig), 4)
    expect_equal(igraph::V(ig)$age, c(31L, 45L, NA, 22L))
    expect_equal(igraph::E(ig)$kind, c("a", "b", "c", "d"))
    expect_equal(igraph::graph_attr(ig, "title"), "sample")
    back <- net_from(ig)
    expect_same_network(g, back)
    expect_equal(back$roles$weight, "weight")
    expect_equal(back$nodes$member, g$nodes$member)
    expect_equal(back$graph$title, "sample")
    # Measures agree across the conversion, weighted (igraph's default when a
    # weight attribute exists) and unweighted.
    expect_equal(net_betweenness(back), unname(igraph::betweenness(ig)))
    expect_equal(net_betweenness(back, weighted = FALSE),
                 unname(igraph::betweenness(ig, weights = NA)))
  }
  expect_s3_class(igraph::as.igraph(sample_net()), "igraph")
})

test_that("unnamed igraph graphs use vertex positions", {
  skip_if_not_installed("igraph")
  g <- net_from(igraph::make_ring(5))
  expect_equal(g$nodes$name, 1:5)
  expect_false(g$directed)
})

test_that("tidygraph round trip", {
  skip_if_not_installed("tidygraph")
  tg <- net_to_tbl_graph(sample_net())
  expect_s3_class(tg, "tbl_graph")
  expect_same_network(sample_net(), net_from(tg))
})

test_that("network (statnet) round trip keeps attributes", {
  skip_if_not_installed("network")
  g <- sample_net()
  net <- net_to_network(g)
  expect_equal(network::network.size(net), 4)
  expect_equal(network::network.edgecount(net), 4)
  expect_true(network::is.multiplex(net))
  expect_true(network::has.loops(net))
  expect_equal(network::get.vertex.attribute(net, "age"), c(31L, 45L, NA, 22L))
  back <- net_from(net)
  expect_same_network(g, back)
  expect_equal(back$edges$kind, g$edges$kind)
  expect_equal(back$edges$weight, g$edges$weight)
  expect_equal(back$graph$title, "sample")
})

test_that("bipartite network objects keep their modes", {
  skip_if_not_installed("network")
  m <- matrix(c(1, 0, 1, 1, 0, 1), 2, 3)
  net <- network::network(m, bipartite = 2, directed = FALSE)
  g <- net_from(net)
  expect_equal(g$nodes$mode, c(1L, 1L, 2L, 2L, 2L))
  back <- net_to_network(g)
  expect_equal(network::get.network.attribute(back, "bipartite"), 2)
})

test_that("adjacency matrices convert both ways", {
  m <- matrix(c(0, 2, 0,
                2, 0, 1,
                0, 1, 0), 3, byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  g <- net_from(m)
  expect_false(g$directed)
  expect_equal(g$edges$weight, c(2, 1))
  expect_equal(as.matrix(g, attr = "weight"), m)

  d <- matrix(c(0, 1, 0, 0), 2, dimnames = list(c("x", "y"), c("x", "y")))
  gd <- net_from(d)
  expect_true(gd$directed)
  expect_null(gd$roles$weight)
  expect_equal(as.matrix(gd), d)

  # Parallel edges are summed.
  multi <- net_frame(data.frame(source = c(1, 1), target = c(2, 2)))
  expect_equal(unname(as.matrix(multi))[1, 2], 2)
})

test_that("sparse matrices convert", {
  skip_if_not_installed("Matrix")
  m <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 3), x = c(1, 1), dims = c(3, 3))
  g <- net_from(m)
  expect_true(g$directed)
  expect_equal(nrow(g$edges), 2)
  s <- Matrix::forceSymmetric(m)
  expect_false(net_from(s)$directed)
  expect_equal(as.matrix(as.matrix(g, sparse = TRUE)), unname(as.matrix(g)), ignore_attr = TRUE)
})

test_that("unsupported objects give a helpful error", {
  expect_error(net_from(1:3), "Supported:")
})
