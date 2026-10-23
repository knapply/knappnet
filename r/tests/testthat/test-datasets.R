test_that("every example network loads with the documented size", {
  info <- example_networks()
  for (i in seq_len(nrow(info))) {
    g <- example_network(info$name[i])
    expect_s3_class(g, "netframe")
    expect_equal(nrow(g$nodes), info$nodes[i], info = info$name[i])
    expect_equal(nrow(g$edges), info$edges[i], info = info$name[i])
  }
})

test_that("Noordin Top relations can be selected", {
  g <- example_network("noordin")
  expect_setequal(unique(g$edges$relation),
                  c("business_finance", "classmates", "communications", "friendship",
                    "kinship", "logistics", "meetings", "operations", "soulmates",
                    "training"))
  expect_false(g$directed)
  ops <- example_network("noordin", relation = "operations")
  expect_true(all(ops$edges$relation == "operations"))
  expect_equal(nrow(ops$nodes), 79L)
  expect_error(example_network("noordin", relation = "nope"), "Unknown `relation`")
})

test_that("two-mode examples are two-mode", {
  aff <- example_network("noordin_affiliations")
  expect_true(is_twomode(aff))
  expect_equal(sum(aff$nodes$event), 45L)
  sw <- example_network("southern_women")
  expect_true(is_twomode(sw))
})

test_that("weights are opt-in", {
  expect_null(example_network("jemaah_islamiyah")$roles$weight)
  ji <- example_network("jemaah_islamiyah", weighted = TRUE)
  expect_identical(ji$roles$weight, "strength")
  expect_true(all(ji$edges$strength %in% 1:5))
  expect_null(example_network("karate")$roles$weight)
})

test_that("the karate club matches igraph's copy", {
  skip_if_not_installed("igraph")
  k <- example_network("karate")
  z <- igraph::make_graph("Zachary")
  expect_equal(betweenness(k), igraph::betweenness(z))
})
