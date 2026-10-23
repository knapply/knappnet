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

test_that("missing and empty identifiers are errors that name the rows", {
  expect_error(netframe(data.frame(from = c("a", ""), to = c("b", "c"))),
               "empty strings; row 2 of `edges`")
  expect_error(netframe(data.frame(from = c("a", "b", NA), to = c("b", NA, "c"))),
               "missing values \\(NA\\); row 3 of `edges`")
  expect_error(netframe(data.frame(from = "a", to = "b"),
                        data.frame(id = c("a", "b", NA))),
               "row 3 of `nodes`")
  expect_error(netframe(data.frame(from = "a", to = "b"),
                        data.frame(id = c("a", "b", ""))),
               "empty strings")
})

test_that("factor, numeric, and date identifiers work", {
  f <- netframe(data.frame(from = factor(c("x", "y")), to = factor(c("y", "z"))))
  expect_identical(f$nodes$name, c("x", "y", "z"))
  expect_equal(degree(f), c(1, 2, 1))

  n <- netframe(data.frame(from = c(1.5, 2.5), to = c(2.5, 3.5)))
  expect_identical(n$nodes$name, c(1.5, 2.5, 3.5))

  d <- as.Date("2024-01-01") + 0:2
  g <- netframe(data.frame(from = d[1:2], to = d[2:3]))
  expect_s3_class(g$nodes$name, "Date")
  expect_equal(degree(g), c(1, 2, 1))
  h <- netframe(data.frame(from = d[1:2], to = d[2:3]), nodes = data.frame(day = d))
  expect_equal(degree(h), c(1, 2, 1))
  expect_error(netframe(data.frame(from = d[1:2], to = c("a", "b"))), "different types")
})

test_that("networks without edges are allowed", {
  e <- data.frame(from = character(), to = character())
  g <- netframe(e)
  expect_equal(nrow(g$nodes), 0)
  expect_equal(degree(g), numeric(0))
  expect_equal(betweenness(g), numeric(0))
  iso <- netframe(e, nodes = data.frame(id = c("a", "b", "c")), directed = FALSE)
  expect_output(print(iso), "3 nodes and 0 edges, unweighted")
  expect_equal(degree(iso), c(0, 0, 0))
  expect_equal(betweenness(iso), c(0, 0, 0))
  expect_equal(components(iso)$no, 3)
})

test_that("other arguments of netframe() are checked", {
  e <- data.frame(from = "a", to = "b")
  expect_error(netframe(e, directed = "yes"), "`directed` must be TRUE or FALSE")
  expect_error(netframe(e, directed = NA), "`directed` must be TRUE or FALSE")
  expect_error(netframe(list(from = "a", to = "b")), "`edges` must be a data frame")
  expect_error(netframe(e, nodes = c("a", "b")), "`nodes` must be a data frame or NULL")
  expect_error(netframe(e, from = "src"), "Column for `from` not found in `edges`")
  expect_error(netframe(e, graph = list("untitled")), "must be named")
  expect_error(netframe(data.frame(from = "a", to = "b", w = "1"), weight = "w"),
               "must be numeric")
  expect_error(netframe(e, data.frame(id = c("a", "b", "a"))),
               "`a` appears more than once")
})

test_that("self-loops, parallel edges, and missing weights are reported", {
  plain <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")))
  expect_identical(capture.output(print(plain)),
                   "<netframe> directed network with 3 nodes and 2 edges, unweighted")
  g <- netframe(data.frame(from = c("a", "a", "a", "b"), to = c("a", "b", "b", "c"),
                           w = c(1, 2, NA, 1)), weight = "w")
  expect_output(print(g), "including 1 self-loop and 1 parallel edge")
  expect_output(print(g), "1 value in weight column `w` is missing or not finite")
  # In an undirected network a-b and b-a are parallel.
  u <- netframe(data.frame(from = c("a", "b"), to = c("b", "a")), directed = FALSE)
  expect_output(print(u), "including 1 parallel edge")
})

test_that("missing and non-finite weights stop every measure with the rows", {
  g <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
                           w = c(1, NaN, Inf)), weight = "w")
  expect_error(degree(g), "rows 2 and 3 of `x\\$edges`")
  expect_error(betweenness(g, weighted = FALSE), "non-finite")
  # Functions that only read the tables carry them through.
  expect_false(is.finite(as.matrix(g, attr = "w")["b", "c"]))
})

test_that("zero and negative weights are checked by the measures that use them", {
  g <- netframe(data.frame(from = c(1, 2), to = c(2, 3), w = c(1, -2)), weight = "w")
  expect_equal(strength(g), c(1, -1, -2))
  expect_error(betweenness(g),
               "betweenness\\(\\): the edge in row 2 of `x\\$edges` has weight -2")
  expect_error(betweenness(g), "weighted = FALSE")
  expect_equal(betweenness(g, weighted = FALSE), c(0, 1, 0))
  expect_error(eigenvector_centrality(g), "eigenvector_centrality\\(\\): the edge in row 2")
  expect_error(closeness(g), "closeness\\(\\): the edge in row 2")
})

test_that("engine messages are rewritten for R users", {
  expect_identical(
    engine_message("degree", "edge 0 has weight -1; this requires nonnegative weights"),
    paste0("degree(): the edge in row 1 of `x$edges` has weight -1; this requires ",
           "nonnegative weights. Use `weighted = FALSE` to ignore weights, or remove ",
           "or transform the offending weights."))
  expect_identical(engine_message("f", "edge 4 has weight 0; strictly positive", FALSE),
                   paste0("f(): an edge has weight 0; strictly positive. Use ",
                          "`weighted = FALSE` to ignore weights, or remove or ",
                          "transform the offending weights."))
  expect_identical(engine_message("read_network", "gml file, line 3: bad token."),
                   "read_network(): gml file, line 3: bad token.")
})

test_that("cached engine graphs follow changes to the weights", {
  g <- netframe(data.frame(from = c("a", "b", "a"), to = c("b", "c", "c"), w = c(1, 1, 4)),
                directed = FALSE, weight = "w")
  expect_equal(betweenness(g), c(0, 1, 0))
  h <- g
  h$edges$w <- c(1, 4, 1)
  expect_equal(betweenness(h), c(1, 0, 0))
  expect_equal(betweenness(g), c(0, 1, 0))
  expect_equal(betweenness(g, convention = "networkx"), c(0, 1, 0))
  expect_equal(betweenness(h, convention = "networkx"), c(1, 0, 0))
})

test_that("changed rows are detected", {
  g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")))
  h <- g
  h$edges <- h$edges[1, , drop = FALSE]
  expect_error(degree(h), "Rebuild it with netframe")
  k <- g
  k$nodes <- k$nodes[1:2, , drop = FALSE]
  expect_error(degree(k), "Rebuild it with netframe")
})

test_that("measures check their arguments", {
  g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")))
  expect_error(degree(data.frame(from = "a", to = "b")), "`x` must be a netframe")
  expect_error(betweenness(g, threads = -1), "`threads` must be")
  expect_error(betweenness(g, threads = 1.5), "`threads` must be")
  expect_error(betweenness(g, threads = "2"), "`threads` must be")
  expect_error(betweenness(g, threads = NA), "`threads` must be")
  expect_error(betweenness(g, normalized = "yes"), "`normalized` must be TRUE or FALSE")
  expect_error(betweenness(g, weights_as = "cost"), "`weights_as` must be one of")
  expect_error(degree(g, loops = NA), "`loops` must be TRUE or FALSE")
  expect_error(closeness(g, mode = "sideways"), "`mode` must be one of")
  expect_error(eigenvector_centrality(g, scale = "l2"), "`scale` must be one of")
  expect_error(constraint(g, weighted = 1), "`weighted` must be TRUE, FALSE, or NULL")
  expect_equal(betweenness(g, threads = 1), betweenness(g, threads = 2))
})

test_that("the default thread count honours CRAN's limit", {
  vars <- c("_R_CHECK_LIMIT_CORES_", "OMP_THREAD_LIMIT")
  old <- Sys.getenv(vars, unset = NA)
  on.exit({
    for (v in vars) {
      if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
    }
  })
  Sys.unsetenv(vars)
  expect_identical(resolve_threads(0), 0L)
  expect_identical(resolve_threads(3), 3L)
  Sys.setenv(OMP_THREAD_LIMIT = "3")
  expect_identical(resolve_threads(0), 3L)
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_identical(resolve_threads(0), 2L)
  expect_identical(resolve_threads(1), 1L)
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "false")
  Sys.unsetenv("OMP_THREAD_LIMIT")
  expect_identical(resolve_threads(0), 0L)
})
