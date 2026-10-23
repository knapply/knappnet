test_that("node set is derived from edges in order of appearance", {
  g <- net_frame(data.frame(source = c("b", "a"), target = c("c", "b")))
  expect_equal(g$nodes$id, c("b", "a", "c"))
  expect_equal(g$index$src, c(0L, 1L))
  expect_equal(g$index$tgt, c(2L, 0L))
  expect_s3_class(g, "net_frame")
  expect_identical(class(g), "net_frame")
})

test_that("explicit nodes keep their order and attributes", {
  nodes <- data.frame(id = c("x", "y", "z"), group = c(1, 1, 2))
  edges <- data.frame(src = "z", dst = "x")
  g <- net_frame(edges, nodes, source = "src", target = "dst", id = "id", quiet = TRUE)
  expect_equal(g$nodes$group, c(1, 1, 2))
  expect_equal(net_degree(g), c(1, 0, 1))
  expect_identical(g$roles[c("source", "target", "id")],
                   list(source = "src", target = "dst", id = "id"))
})

test_that("endpoint columns are named explicitly, never guessed", {
  e <- data.frame(f = "a", t = "b")
  expect_error(net_frame(e), "Column for `source` not found in `edges`")
  expect_error(net_frame(e, source = "f"), "Column for `target` not found in `edges`")
  expect_error(net_frame(e, source = 1L), "Column for `source` not found")
  g <- net_frame(e, source = "f", target = "t")
  expect_equal(g$nodes$id, c("a", "b"))
  nodes <- data.frame(who = c("a", "b"))
  expect_error(net_frame(e, nodes, source = "f", target = "t"),
               "Column for `id` not found in `nodes`")
  expect_equal(nrow(net_frame(e, nodes, source = "f", target = "t", id = "who")$nodes), 2L)
})

test_that("every error is reported at once, with the rows", {
  nodes <- data.frame(id = c("a", "b"))
  expect_error(net_frame(data.frame(source = "a", target = "q"), nodes),
               "must be identifiers in column `id` of `nodes`")
  expect_error(net_frame(data.frame(source = "a", target = "q"), nodes), "not in `nodes`: q")
  expect_error(net_frame(data.frame(source = "a", target = "b"), data.frame(id = c("a", "a"))),
               "must be unique")
  expect_error(net_frame(data.frame(source = NA, target = "b")), "must not be missing")
  expect_error(net_frame(data.frame(source = "a", target = "b"), weight = "w"),
               "Column for `weight` not found")
  err <- tryCatch(net_frame(data.frame(source = c("a", NA), target = c("q", "b")),
                            data.frame(id = c("a", "b", "b"))),
                  error = function(e) conditionMessage(e))
  expect_match(err, "problems")
  expect_match(err, "must not be missing; row 2 of `edges`")
  expect_match(err, "repeat an earlier identifier; row 3 of `nodes`")
  expect_match(err, "repeated: b")
  expect_error(net_betweenness(net_frame(data.frame(source = "a", target = "b")),
                               weighted = TRUE), "no weights")
})

test_that("weights must declare their meaning", {
  e <- data.frame(source = c("a", "b"), target = c("b", "c"), w = c(1, 2))
  expect_error(net_frame(e, weight = "w"), "`weight_type` must say what the weights mean")
  expect_error(net_frame(e, weight = "w", weight_type = "cost"), "`weight_type` must be one of")
  g <- net_frame(e, weight = "w", weight_type = "strength")
  expect_identical(g$roles$weight_type, "strength")
  expect_identical(net_frame(e, weight = "w", weight_type = "distance")$roles$weight_type,
                   "distance")
  expect_error(net_frame(data.frame(source = "a", target = "b", w = "1"), weight = "w",
                         weight_type = "strength"), "must be numeric")
  expect_error(net_frame(data.frame(source = "a", target = "b", w = 0), weight = "w",
                         weight_type = "distance"), "distances and must be positive")
  expect_error(net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "a"),
                                    w = c(1, NaN, Inf)), weight = "w", weight_type = "strength"),
               "rows 2 and 3 of `edges`")
})

test_that("nodes added from edges are appended and noted", {
  nodes <- data.frame(id = c("a", "b"), age = c(30, 40))
  e <- data.frame(source = c("a", "c"), target = c("b", "d"))
  expect_message(g <- net_frame(e, nodes, add_nodes = TRUE), "added from the edge endpoints")
  expect_equal(g$nodes$id, c("a", "b", "c", "d"))
  expect_equal(g$nodes$age, c(30, 40, NA, NA))
  expect_silent(net_frame(e, nodes, add_nodes = TRUE, quiet = TRUE))
  expect_error(net_frame(e, nodes), "not in `nodes`: c, d")
})

test_that("notes report isolates, loops, parallel edges, and 0-based ids", {
  expect_message(net_frame(data.frame(source = 1L, target = 2L), data.frame(id = 1:3)),
                 "have no edges")
  expect_message(net_frame(data.frame(source = c(1, 1), target = c(1, 2))), "self-loop")
  expect_message(net_frame(data.frame(source = c(1, 1), target = c(2, 2))), "parallel edge")
  expect_message(net_frame(data.frame(source = 0L, target = 1L)), "run from 0 to 1")
  d <- net_diagnostics(net_frame(data.frame(source = c(1, 1), target = c(2, 2)), quiet = TRUE))
  expect_s3_class(d, "data.frame")
  expect_true("N_PARALLEL" %in% d$code)
  expect_true(all(d$severity == "note"))
  expect_equal(d$rows[[match("N_PARALLEL", d$code)]], 2L)
})

test_that("declared properties are enforced", {
  expect_error(net_frame(data.frame(source = c(1, 1), target = c(1, 2)), loops = FALSE),
               "declared absent \\(loops = FALSE\\)")
  expect_error(net_frame(data.frame(source = c(1, 1), target = c(2, 2)), multiple = FALSE),
               "declared absent \\(multiple = FALSE\\)")
  g <- net_frame(data.frame(source = c(1, 1), target = c(2, 2)), multiple = TRUE)
  expect_true(g$declared$multiple)
  expect_output(print(net_frame(data.frame(source = 1, target = 2), loops = FALSE,
                                multiple = FALSE)),
                "no loops (declared), no parallel edges (declared)", fixed = TRUE)
})

test_that("engine handles survive serialization", {
  g <- net_frame(data.frame(source = 1:3, target = 2:4), directed = FALSE)
  first <- net_betweenness(g)
  path <- tempfile(fileext = ".rds")
  saveRDS(g, path)
  h <- readRDS(path)
  expect_equal(net_betweenness(h), first)
  expect_output(print(h), "4 nodes, 3 edges")
})

test_that("engine errors surface as R errors", {
  g <- net_frame(data.frame(source = 1, target = 2, w = 0.5), weight = "w",
                 weight_type = "strength")
  g$edges$w <- 0
  expect_error(net_betweenness(g, weighted = TRUE), "strictly positive")
})

test_that("missing and empty identifiers are errors that name the rows", {
  expect_error(net_frame(data.frame(source = c("a", ""), target = c("b", "c"))),
               "empty strings; row 2 of `edges`")
  expect_error(net_frame(data.frame(source = c("a", "b", NA), target = c("b", NA, "c"))),
               "must not be missing; rows 2 and 3 of `edges`")
  expect_error(net_frame(data.frame(source = "a", target = "b"),
                         data.frame(id = c("a", "b", NA))),
               "row 3 of `nodes`")
  expect_error(net_frame(data.frame(source = "a", target = "b"),
                         data.frame(id = c("a", "b", ""))),
               "empty strings")
})

test_that("factor, numeric, logical, and date identifiers work", {
  f <- net_frame(data.frame(source = factor(c("x", "y")), target = factor(c("y", "z"))))
  expect_identical(f$nodes$id, c("x", "y", "z"))
  expect_equal(net_degree(f), c(1, 2, 1))

  n <- net_frame(data.frame(source = c(1.5, 2.5), target = c(2.5, 3.5)))
  expect_identical(n$nodes$id, c(1.5, 2.5, 3.5))

  l <- net_frame(data.frame(source = c(TRUE, FALSE), target = c(FALSE, TRUE)))
  expect_equal(nrow(l$nodes), 2L)

  d <- as.Date("2024-01-01") + 0:2
  g <- net_frame(data.frame(source = d[1:2], target = d[2:3]))
  expect_s3_class(g$nodes$id, "Date")
  expect_equal(net_degree(g), c(1, 2, 1))
  h <- net_frame(data.frame(source = d[1:2], target = d[2:3]), nodes = data.frame(day = d),
                 id = "day")
  expect_equal(net_degree(h), c(1, 2, 1))
  expect_error(net_frame(data.frame(source = d[1:2], target = c("a", "b"))), "different types")
})

test_that("networks without edges are allowed", {
  e <- data.frame(source = character(), target = character())
  g <- net_frame(e)
  expect_equal(nrow(g$nodes), 0)
  expect_identical(names(g$nodes), "id")
  expect_equal(net_degree(g), numeric(0))
  expect_equal(net_betweenness(g), numeric(0))
  iso <- net_frame(e, nodes = data.frame(id = c("a", "b", "c")), directed = FALSE,
                   quiet = TRUE)
  expect_output(print(iso), "3 nodes, 0 edges")
  expect_output(print(iso), "isolates: 3")
  expect_equal(net_degree(iso), c(0, 0, 0))
  expect_equal(net_betweenness(iso), c(0, 0, 0))
  expect_equal(net_components(iso)$no, 3)
})

test_that("other arguments of net_frame() are checked", {
  e <- data.frame(source = "a", target = "b")
  expect_error(net_frame(e, directed = "yes"), "`directed` must be TRUE or FALSE")
  expect_error(net_frame(e, directed = NA), "`directed` must be TRUE or FALSE")
  expect_error(net_frame(list(source = "a", target = "b")), "`edges` must be a data frame")
  expect_error(net_frame(e, nodes = c("a", "b")), "`nodes` must be a data frame or NULL")
  expect_error(net_frame(e, graph = list("untitled")), "must be named")
  expect_error(net_frame(e, data.frame(id = c("a", "b", "a"))), "repeated: a")
  expect_error(net_frame(e, loops = "no"), "`loops` must be TRUE, FALSE, or NULL")
  expect_error(net_frame(e, mode = "m"), "`mode` names a column of `nodes`")
  expect_error(net_frame(e, relation = "r"), "Column for `relation` not found")
})

test_that("modes, relations, and signs are validated and stored", {
  nodes <- data.frame(id = c("a", "b", "c"), kind = c("p", "e", NA))
  expect_error(net_frame(data.frame(source = "a", target = "b"), nodes, mode = "kind"),
               "every node needs a mode")
  nodes$kind[3] <- "e"
  g <- net_frame(data.frame(source = "a", target = "b"), nodes, mode = "kind", quiet = TRUE)
  expect_identical(g$roles$mode, "kind")
  expect_false(net_is_twomode(g))
  expect_output(print(g), "mode: 2 modes (1 p, 2 e)", fixed = TRUE)

  e <- data.frame(source = c("a", "a"), target = c("b", "b"), rel = c("x", "y"))
  m <- net_frame(e, relation = "rel")
  expect_true(net_is_simple(m))
  expect_output(print(m), "relations: x (1), y (1) (multiplex)", fixed = TRUE)
  expect_error(net_frame(e, relation = "rel", sign = "rel"), "sign column")
  s <- net_frame(data.frame(source = c("a", "a"), target = c("b", "c"), s = c(1L, -1L)),
                 sign = "s")
  expect_output(print(s), "signs:     s (1 positive, 1 negative)", fixed = TRUE)
  expect_error(net_frame(data.frame(source = "a", target = "b", s = 2L), sign = "s"),
               "must be -1 or 1")
})

test_that("the print card shows a one-mode network", {
  g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")))
  out <- sub(" +$", "", capture.output(print(g)))
  expect_identical(out[1:7], c(
    "# net_frame: 3 nodes, 2 edges",
    "#   directed:  yes       mode: one-mode (3 actors)",
    "#   weights:   unweighted",
    "#   time:      static    space: none",
    "#   simple:    no loops, no parallel edges",
    "#   components: 1    isolates: 0    density: 0.333",
    "#"))
  expect_identical(out[8:13], c(
    "# Nodes (3 x 1)   Edges (2 x 2)",
    "#   id              source target",
    "#   <chr>           <chr>  <chr>",
    "# 1 a             1 a      b",
    "# 2 b             2 b      c",
    "# 3 c"))
})

test_that("the print card shows weights, loops, parallel edges, and empty networks", {
  g <- net_frame(data.frame(source = c("a", "a", "a", "b"), target = c("a", "b", "b", "c"),
                            w = c(1, 2, 4, 1)), weight = "w", weight_type = "strength",
                 quiet = TRUE)
  out <- capture.output(print(g))
  expect_identical(out[3], "#   weights:   w (strength, 1 to 4)")
  expect_identical(out[5], "#   simple:    1 loop, 1 parallel edge")
  d <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c"), w = c(2, 5)),
                 weight = "w", weight_type = "distance", loops = FALSE, multiple = FALSE)
  out <- capture.output(print(d))
  expect_identical(out[3], "#   weights:   w (distance, 2 to 5)")
  expect_identical(out[5], "#   simple:    no loops (declared), no parallel edges (declared)")
  u <- net_frame(data.frame(source = c("a", "b"), target = c("b", "a")), directed = FALSE,
                 quiet = TRUE)
  expect_output(print(u), "no loops, 1 parallel edge")
  empty <- net_frame(data.frame(source = character(), target = character()))
  out <- sub(" +$", "", capture.output(print(empty)))
  expect_identical(out[1], "# net_frame: 0 nodes, 0 edges")
  expect_identical(out[6], "#   components: 0    isolates: 0    density: NaN")
  expect_true(any(out == "# Nodes (0 x 1)   Edges (0 x 2)"))
})

test_that("the print card stacks the tables on a narrow console", {
  g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")))
  out <- capture.output(print(g, width = 30))
  expect_true("# Nodes (3 x 1)" %in% out)
  expect_true("# Edges (2 x 2)" %in% out)
  out <- capture.output(print(g, n = 1))
  expect_true(any(grepl("2 more rows", out)))
  wide <- net_frame(data.frame(source = "a", target = "b", c1 = 1, c2 = 2, c3 = 3, c4 = 4,
                               c5 = 5, c6 = 6, c7 = 7, c8 = 8))
  out <- capture.output(print(wide, width = 50))
  expect_true(any(grepl("more columns", out)))
})

test_that("the print card uses only ASCII", {
  g <- net_example("karate")
  out <- capture.output(print(g))
  expect_true(all(vapply(out, function(l) all(utf8ToInt(l) < 128), TRUE)))
  expect_identical(out[1], "# net_frame: 34 nodes, 78 edges")
  expect_identical(out[2], "#   directed:  no        mode: one-mode (34 actors)")
})

test_that("engine messages are rewritten for R users", {
  expect_identical(
    engine_message("net_degree", "edge 0 has weight -1; this requires nonnegative weights"),
    paste0("net_degree(): the edge in row 1 of `x$edges` has weight -1; this requires ",
           "nonnegative weights. Use `weighted = FALSE` to ignore weights, or remove ",
           "or transform the offending weights."))
  expect_identical(engine_message("f", "edge 4 has weight 0; strictly positive", FALSE),
                   paste0("f(): an edge has weight 0; strictly positive. Use ",
                          "`weighted = FALSE` to ignore weights, or remove or ",
                          "transform the offending weights."))
  expect_identical(engine_message("net_read", "gml file, line 3: bad token."),
                   "net_read(): gml file, line 3: bad token.")
})

test_that("cached engine graphs follow changes to the weights", {
  g <- net_frame(data.frame(source = c("a", "b", "a"), target = c("b", "c", "c"),
                            w = c(1, 1, 4)), directed = FALSE, weight = "w",
                 weight_type = "distance")
  expect_equal(net_betweenness(g), c(0, 1, 0))
  h <- g
  h$edges$w <- c(1, 4, 1)
  expect_equal(net_betweenness(h), c(1, 0, 0))
  expect_equal(net_betweenness(g), c(0, 1, 0))
  expect_equal(net_betweenness(g, convention = "networkx"), c(0, 1, 0))
  expect_equal(net_betweenness(h, convention = "networkx"), c(1, 0, 0))
  # Edited weights are validated again.
  k <- g
  k$edges$w <- c(1, NA, 1)
  expect_error(net_degree(k), "row 2 of `edges`")
})

test_that("changed rows are detected", {
  g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")))
  h <- g
  h$edges <- h$edges[1, , drop = FALSE]
  expect_error(net_degree(h), "Rebuild it with net_frame")
  k <- g
  k$nodes <- k$nodes[1:2, , drop = FALSE]
  expect_error(net_degree(k), "Rebuild it with net_frame")
})

test_that("net_simplify() follows the role-aware parallel rule", {
  e <- data.frame(source = c("a", "a", "a", "b", "b"), target = c("b", "b", "a", "c", "c"),
                  w = c(1, 2, 5, 1, 1), relation = c("x", "x", "x", "x", "y"))
  g <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength",
                 relation = "relation", quiet = TRUE)
  expect_equal(net_multiplicity(g), c(2L, 2L, 1L, 1L, 1L))
  expect_false(net_is_simple(g))
  s <- net_simplify(g)
  expect_equal(nrow(s$edges), 4L)
  expect_equal(s$edges$w, c(3, 5, 1, 1))
  expect_true(net_is_simple(net_simplify(g, remove_loops = TRUE)))
  expect_equal(net_simplify(g, combine = "max")$edges$w[1], 2)
  expect_equal(net_simplify(g, combine = "mean")$edges$w[1], 1.5)
  expect_equal(net_simplify(g, combine = "first")$edges$w[1], 1)
  expect_false(s$declared$multiple)
  u <- net_frame(data.frame(source = c("a", "b", "a"), target = c("b", "a", "c")),
                 directed = FALSE, quiet = TRUE)
  c_ <- net_simplify(u, combine = "count")
  expect_equal(c_$edges$count, c(2, 1))
  expect_identical(c_$roles$weight, "count")
  expect_identical(c_$roles$weight_type, "strength")
  # In a temporal network two spells of one tie are not parallel.
  t <- net_temporal(data.frame(source = c("a", "a"), target = c("b", "b"), time = c(1, 2)))
  expect_true(net_is_simple(t))
})

test_that("measures check their arguments", {
  g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")))
  expect_error(net_degree(data.frame(source = "a", target = "b")), "`x` must be a net_frame")
  expect_error(net_betweenness(g, threads = -1), "`threads` must be")
  expect_error(net_betweenness(g, threads = 1.5), "`threads` must be")
  expect_error(net_betweenness(g, threads = "2"), "`threads` must be")
  expect_error(net_betweenness(g, threads = NA), "`threads` must be")
  expect_error(net_betweenness(g, normalized = "yes"), "`normalized` must be TRUE or FALSE")
  expect_error(net_betweenness(g, weights_as = "cost"), "`weights_as` must be one of")
  expect_error(net_degree(g, loops = NA), "`loops` must be TRUE or FALSE")
  expect_error(net_closeness(g, mode = "sideways"), "`mode` must be one of")
  expect_error(net_eigen(g, scale = "l2"), "`scale` must be one of")
  expect_error(net_constraint(g, weighted = 1), "`weighted` must be TRUE, FALSE, or NULL")
  expect_equal(net_betweenness(g, threads = 1), net_betweenness(g, threads = 2))
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
