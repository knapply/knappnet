test_that("node set is derived from edges in order of appearance", {
  g <- netframe(data.frame(from = c("b", "a"), to = c("c", "b")))
  expect_equal(g$nodes$name, c("b", "a", "c"))
  expect_equal(g$index$src, c(0L, 1L))
  expect_equal(g$index$tgt, c(2L, 0L))
})

test_that("explicit nodes keep their order and attributes", {
  nodes <- data.frame(id = c("x", "y", "z"), group = c(1, 1, 2))
  edges <- data.frame(src = "z", dst = "x")
  g <- netframe(edges, nodes, from = "src", to = "dst", id = "id")
  expect_equal(g$nodes$group, c(1, 1, 2))
  expect_equal(degree(g), c(1, 0, 1))
})

test_that("bad input gives clear errors", {
  nodes <- data.frame(id = c("a", "b"))
  expect_error(netframe(data.frame(f = "a", t = "q"), nodes), "not in `nodes`")
  expect_error(netframe(data.frame(f = "a", t = "b"), data.frame(id = c("a", "a"))),
               "unique")
  expect_error(netframe(data.frame(f = NA, t = "b")), "missing")
  expect_error(netframe(data.frame(f = "a", t = "b"), weight = "w"), "not found")
  expect_error(betweenness(netframe(data.frame(f = "a", t = "b")), weighted = TRUE),
               "no weights")
})

test_that("engine handles survive serialization", {
  g <- netframe(data.frame(from = 1:3, to = 2:4), directed = FALSE)
  first <- betweenness(g)
  path <- tempfile(fileext = ".rds")
  saveRDS(g, path)
  h <- readRDS(path)
  expect_equal(betweenness(h), first)
})

test_that("engine errors surface as R errors", {
  g <- netframe(data.frame(from = 1, to = 2, w = 0), weight = "w")
  expect_error(betweenness(g, weighted = TRUE), "strictly positive")
})
