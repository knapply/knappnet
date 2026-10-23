test_that("every example network loads with the documented size", {
  info <- net_examples()
  for (i in seq_len(nrow(info))) {
    g <- net_example(info$name[i])
    expect_s3_class(g, "net_frame")
    expect_equal(nrow(g$nodes), info$nodes[i], info = info$name[i])
    expect_equal(nrow(g$edges), info$edges[i], info = info$name[i])
  }
})

test_that("Noordin Top relations can be selected", {
  g <- net_example("noordin")
  expect_setequal(unique(g$edges$relation),
                  c("business_finance", "classmates", "communications", "friendship",
                    "kinship", "logistics", "meetings", "operations", "soulmates",
                    "training"))
  expect_false(g$directed)
  ops <- net_example("noordin", relation = "operations")
  expect_true(all(ops$edges$relation == "operations"))
  expect_equal(nrow(ops$nodes), 79L)
  expect_error(net_example("noordin", relation = "nope"), "Unknown `relation`")
})

test_that("two-mode examples are two-mode", {
  aff <- net_example("noordin_affiliations")
  expect_true(net_is_twomode(aff))
  expect_equal(sum(aff$nodes$event), 45L)
  sw <- net_example("southern_women")
  expect_true(net_is_twomode(sw))
})

test_that("weights are opt-in", {
  expect_null(net_example("jemaah_islamiyah")$roles$weight)
  ji <- net_example("jemaah_islamiyah", weighted = TRUE)
  expect_identical(ji$roles$weight, "strength")
  expect_true(all(ji$edges$strength %in% 1:5))
  expect_null(net_example("karate")$roles$weight)
})

test_that("the karate club matches igraph's copy", {
  skip_if_not_installed("igraph")
  k <- net_example("karate")
  z <- igraph::make_graph("Zachary")
  expect_equal(net_betweenness(k), igraph::betweenness(z))
})
