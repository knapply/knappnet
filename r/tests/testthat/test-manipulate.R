# Data manipulation: combining relations, dichotomizing, symmetrizing,
# transposing, subnetworks, isolates and pendants, components, recoding and
# joining attributes.

multiplex <- function() {
  e <- data.frame(from = c("a", "a", "b", "a", "c", "a"),
                  to = c("b", "c", "c", "b", "d", "b"),
                  relation = c("friend", "friend", "friend", "work", "work", "friend"),
                  stringsAsFactors = FALSE)
  netframe(e, nodes = data.frame(name = c("a", "b", "c", "d")), directed = FALSE)
}

test_that("combine_relations by union, intersection, sum, difference, and select", {
  g <- multiplex()
  u <- combine_relations(g, method = "union")
  expect_equal(nrow(u$edges), 4)
  expect_equal(paste(u$edges$from, u$edges$to), c("a b", "a c", "b c", "c d"))
  expect_equal(u$edges$relations, c(2, 1, 1, 1))
  expect_null(u$roles$weight)
  i <- combine_relations(g, method = "intersection")
  expect_equal(paste(i$edges$from, i$edges$to), "a b")
  s <- combine_relations(g, method = "sum")
  expect_equal(s$roles$weight, "weight")
  expect_equal(s$roles$weight_type, "strength")
  expect_equal(s$edges$weight, c(2, 1, 1, 1))
  d <- combine_relations(g, method = "difference", which = c("friend", "work"))
  expect_equal(paste(d$edges$from, d$edges$to), c("a c", "b c"))
  w <- combine_relations(g, method = "select", which = "work")
  expect_equal(nrow(w$edges), 2)
  expect_equal(names(w$edges), names(g$edges))
  expect_equal(u$nodes, g$nodes)
  expect_error(combine_relations(g, method = "difference"), "exactly two")
  expect_error(combine_relations(g, which = "foe"), "Unknown relation")
  expect_error(combine_relations(netframe(data.frame(from = 1, to = 2))), "relation")
  # Weighted sums add weights; directed pairs stay ordered.
  e <- data.frame(from = c(1, 2, 1), to = c(2, 1, 2), rel = c("x", "x", "y"), w = c(1, 2, 3))
  dw <- netframe(e, weight = "w")
  sw <- combine_relations(dw, relation = "rel", method = "sum")
  expect_equal(sw$edges$weight, c(4, 2))
  expect_equal(sw$roles$weight_type, "distance")
})

test_that("dichotomize keeps edges passing the test and drops the weight role", {
  e <- data.frame(from = 1:4, to = c(2, 3, 4, 1), w = c(1, 2, 3, 4))
  g <- netframe(e, weight = "w", directed = FALSE)
  d <- dichotomize(g, threshold = 2)
  expect_null(d$roles$weight)
  expect_equal(d$edges$w, c(3, 4))
  expect_equal(dichotomize(g, 2, ">=")$edges$w, c(2, 3, 4))
  expect_equal(dichotomize(g, 2, "<")$edges$w, 1)
  expect_equal(dichotomize(g, 2, "==")$edges$w, 2)
  expect_equal(nrow(dichotomize(g, 2, "!=", column = "w")$edges), 3)
  expect_equal(nrow(dichotomize(g, 10, "<=")$edges), 4)
  expect_error(dichotomize(netframe(e)), "no weight column")
  expect_error(dichotomize(g, "a"), "single number")
})

test_that("symmetrize by every method and transpose", {
  d <- netframe(data.frame(from = c(1, 2, 1, 3), to = c(2, 1, 3, 2)))
  key <- function(g) paste(g$edges$from, g$edges$to)
  mx <- symmetrize(d, "max")
  expect_false(mx$directed)
  expect_equal(key(mx), c("1 2", "1 3", "2 3"))
  expect_null(mx$roles$weight)
  expect_equal(key(symmetrize(d, "min")), "1 2")
  sm <- symmetrize(d, "sum")
  expect_equal(key(sm), c("1 2", "1 3", "2 3"))
  expect_equal(sm$edges$weight, c(2, 1, 1))
  expect_equal(sm$roles$weight_type, "strength")
  expect_equal(symmetrize(d, "average")$edges$weight, c(1, 0.5, 0.5))
  expect_equal(key(symmetrize(d, "upper")), c("1 2", "1 3"))
  expect_equal(key(symmetrize(d, "lower")), c("1 2", "2 3"))
  expect_error(symmetrize(mx), "already undirected")
  # Weights: max and min of the two directions, parallel arcs collapse.
  w <- netframe(data.frame(from = c(1, 2, 1, 1), to = c(2, 1, 3, 3), w = c(5, 2, 1, 4)),
                weight = "w")
  expect_equal(symmetrize(w, "max")$edges$weight, c(5, 4))
  expect_equal(symmetrize(w, "min")$edges$weight, 2)
  expect_equal(symmetrize(w, "sum")$edges$weight, c(7, 5))
  expect_equal(symmetrize(w, "max")$roles$weight_type, "distance")
  # A self-loop survives once.
  loop <- netframe(data.frame(from = c(1, 1), to = c(1, 2)))
  expect_equal(key(symmetrize(loop)), c("1 1", "1 2"))

  t1 <- transpose_network(d)
  expect_equal(t1$edges$from, d$edges$to)
  expect_equal(t1$edges$to, d$edges$from)
  expect_true(t1$directed)
  expect_identical(transpose_network(mx)$edges, mx$edges)
})

test_that("subnetworks, isolates, pendants, and the largest component", {
  # A star 1 - {2, 3, 4} with a path 4 - 5 - 6 and an isolate 7.
  g <- netframe(data.frame(from = c(1, 1, 1, 4, 5), to = c(2, 3, 4, 5, 6), w = 1:5),
                nodes = data.frame(id = 1:7, grp = c("a", "a", "a", "b", "b", "b", "c")),
                directed = FALSE, weight = "w")
  s <- subnetwork(g, c(1, 2, 3))
  expect_equal(s$nodes$id, 1:3)
  expect_equal(s$edges$w, 1:2)
  expect_equal(s$roles$weight, "w")
  expect_equal(subnetwork(g, g$nodes$grp == "b")$edges$w, 4:5)
  expect_equal(subnetwork(g, c(3, 1))$nodes$id, c(1, 3))
  expect_error(subnetwork(g, c(1, 9)), "not in")
  expect_equal(filter_nodes(g, grp != "a")$nodes$id, 4:7)
  expect_equal(filter_nodes(g, id > 3 & grp == "b")$nodes$id, 4:6)
  expect_equal(nrow(filter_nodes(g, FALSE)$nodes), 0)
  expect_equal(remove_isolates(g)$nodes$id, 1:6)
  expect_equal(remove_pendants(g)$nodes$id, c(1, 4, 5, 7))
  expect_equal(remove_pendants(g, recursive = TRUE)$nodes$id, c(4, 7))
  expect_equal(largest_component(g)$nodes$id, 1:6)
  # A self-loop does not make a node non-isolated.
  loop <- netframe(data.frame(from = c(1, 2), to = c(1, 3)), directed = FALSE)
  expect_equal(remove_isolates(loop)$nodes$name, c(2, 3))
  d <- netframe(data.frame(from = c(1, 2, 3, 4), to = c(2, 3, 1, 1)))
  expect_equal(largest_component(d, "strong")$nodes$name, 1:3)
  expect_equal(nrow(largest_component(d, "weak")$nodes), 4)
})

test_that("recode_attribute and join_attributes", {
  g <- example_network("jemaah_islamiyah")
  r <- recode_attribute(g, "role", c(`command team` = "command"), default = "operative")
  expect_equal(sort(unique(r$nodes$role)), c("command", "operative"))
  expect_equal(sum(r$nodes$role == "command"), 3)
  keep <- recode_attribute(g, "role", c(`command team` = "command"))
  expect_equal(sum(keep$nodes$role == "bomb maker"), 5)
  f <- recode_attribute(g, "role", toupper)
  expect_equal(f$nodes$role, toupper(g$nodes$role))
  # Recoding identifiers rewrites the endpoints.
  ids <- recode_attribute(g, "name", function(v) paste0("p", seq_along(v)))
  expect_equal(ids$nodes$name[1], "p1")
  expect_equal(nrow(ids$edges), nrow(g$edges))
  expect_equal(degree(ids), degree(g))
  expect_error(recode_attribute(g, "role", c(1, 2)), "named")
  e <- recode_attribute(g, "strength", function(v) v * 2, table = "edges")
  expect_equal(e$edges$strength, g$edges$strength * 2)

  extra <- data.frame(name = g$nodes$name[c(3, 1)], arrested = c(TRUE, FALSE), role = "x")
  j <- join_attributes(g, extra)
  expect_equal(j$nodes$arrested, c(FALSE, NA, TRUE, rep(NA, nrow(g$nodes) - 3)))
  expect_equal(j$nodes$role[1:3], c("x", NA, "x"))
  expect_equal(j$nodes$name, g$nodes$name)
  j2 <- join_attributes(g, data.frame(who = g$nodes$name[2], score = 7), by = "who")
  expect_equal(j2$nodes$score[2], 7)
  expect_error(join_attributes(g, data.frame(name = rep("Muklas", 2), v = 1:2)), "duplicated")
  expect_error(join_attributes(g, 1), "data frame")
})
