rich_net <- function(directed = TRUE) {
  nodes <- data.frame(id = c("ana", "bo \"b\"", "céline", "d&d", "isolate"),
                      age = c(31L, 45L, NA, 22L, 60L),
                      score = c(0.5, NA, 1.25, -3, 1e-7),
                      member = c(TRUE, FALSE, NA, TRUE, FALSE),
                      team = c("red", "blue", "<blue>", NA, "green & gold"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(from = c("ana", "bo \"b\"", "céline", "ana"),
                      to = c("bo \"b\"", "céline", "céline", "d&d"),
                      weight = c(1, 2.5, 0.75, 3),
                      stringsAsFactors = FALSE)
  netframe(edges, nodes, directed = directed, weight = "weight",
           graph = list(title = "rich"))
}

round_trip <- function(g, ext, ...) {
  path <- tempfile(fileext = ext)
  on.exit(unlink(path))
  write_network(g, path, ...)
  read_network(path)
}

test_that("attribute-rich formats round-trip exactly", {
  for (ext in c(".graphml", ".gexf", ".json", ".gml")) {
    for (directed in c(TRUE, FALSE)) {
      g <- rich_net(directed)
      back <- round_trip(g, ext)
      expect_equal(back$directed, directed, info = ext)
      expect_equal(back$nodes$id, g$nodes$id, info = ext)
      expect_equal(back$edges$from, g$edges$from, info = ext)
      expect_equal(back$edges$weight, g$edges$weight, info = ext)
      expect_equal(back$nodes$score, g$nodes$score, info = ext)
      expect_equal(back$nodes$team, g$nodes$team, info = ext)
      if (ext != ".gml") {
        # GML has no boolean type, so logical columns come back as 0/1.
        expect_equal(back$nodes$member, g$nodes$member, info = ext)
        expect_identical(back$nodes$age, g$nodes$age, info = ext)
      }
    }
  }
})

test_that("structure-only formats keep structure and weights", {
  g <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "c"), w = c(2, 0.5, 1)),
                nodes = data.frame(id = c("a", "b", "c", "lonely")), weight = "w")
  for (ext in c(".net", ".dl")) {
    back <- round_trip(g, ext)
    expect_equal(back$nodes$id, g$nodes$id, info = ext)
    expect_equal(back$edges$to, g$edges$to, info = ext)
    expect_equal(back$edges$weight, g$edges$w, info = ext)
  }
})

test_that("edge lists read and write with optional node files", {
  g <- rich_net()
  path <- tempfile(fileext = ".csv")
  nodes_path <- tempfile(fileext = ".csv")
  write_network(g, path, nodes_file = nodes_path)
  back <- read_network(path, nodes = nodes_path)
  expect_equal(nrow(back$nodes), 5)
  expect_equal(back$edges$from, g$edges$from)
  expect_equal(back$roles$weight, "weight")
})

test_that("integer identifiers come back as integers", {
  g <- netframe(data.frame(from = c(1L, 2L), to = c(2L, 3L)))
  back <- round_trip(g, ".graphml")
  expect_identical(back$nodes$id, 1:3)
  expect_identical(back$edges$to, c(2L, 3L))
})

test_that("DOT is write-only and formats are validated", {
  g <- rich_net()
  path <- tempfile(fileext = ".dot")
  write_network(g, path)
  expect_match(readLines(path)[1], "^digraph")
  expect_error(read_network(path), "written but not read")
  expect_error(write_network(g, tempfile(fileext = ".xyz")), "Cannot tell the format")
  expect_error(read_network("does-not-exist.graphml"), "does not exist")
})

test_that("files written here are read correctly by igraph", {
  skip_if_not_installed("igraph")
  g <- rich_net()
  for (fmt in c("graphml", "gml", "pajek")) {
    path <- tempfile()
    write_network(g, path, format = fmt)
    ig <- igraph::read_graph(path, format = fmt)
    expect_equal(igraph::vcount(ig), nrow(g$nodes), info = fmt)
    expect_equal(igraph::ecount(ig), nrow(g$edges), info = fmt)
    expect_equal(sort(igraph::E(ig)$weight), sort(g$edges$weight), info = fmt)
  }
  path <- tempfile(fileext = ".graphml")
  write_network(g, path)
  ig <- igraph::read_graph(path, format = "graphml")
  expect_equal(igraph::V(ig)$id, g$nodes$id)
  expect_equal(igraph::V(ig)$age[-3], g$nodes$age[-3])
})

test_that("files written by igraph are read correctly here", {
  skip_if_not_installed("igraph")
  ig <- igraph::make_graph("Zachary")
  igraph::E(ig)$weight <- seq_len(igraph::ecount(ig)) / 10
  igraph::V(ig)$club <- rep(c("A", "B"), length.out = igraph::vcount(ig))
  for (fmt in c("graphml", "gml", "pajek")) {
    path <- tempfile()
    igraph::write_graph(ig, path, format = fmt)
    g <- read_network(path, format = fmt)
    expect_equal(nrow(g$nodes), igraph::vcount(ig), info = fmt)
    expect_equal(nrow(g$edges), igraph::ecount(ig), info = fmt)
    expect_equal(g$edges$weight, igraph::E(ig)$weight, info = fmt)
    expect_equal(betweenness(g), unname(igraph::betweenness(ig, weights = NA)),
                 info = fmt)
  }
})

test_that("UCINET DL files are read the same way igraph reads them", {
  skip_if_not_installed("igraph")
  path <- tempfile(fileext = ".dl")
  writeLines(c("dl n=4 format=fullmatrix", "labels:", "a,b,c,d", "data:",
               "0 1 0 0", "0 0 2 0", "1 0 0 0", "0 0 0 0"), path)
  g <- read_network(path)
  ig <- igraph::read_graph(path, format = "dl")
  expect_equal(nrow(g$edges), igraph::ecount(ig))
  expect_equal(g$nodes$id, igraph::V(ig)$name)
  expect_equal(g$edges$weight, igraph::E(ig)$weight)
})
