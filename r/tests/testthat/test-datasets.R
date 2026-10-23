test_that("every installed example network loads with the documented size", {
  info <- net_examples()
  expect_true(all(info$available[info$name != "workplace_contacts"]))
  for (i in which(info$available)) {
    g <- net_example(info$name[i])
    expect_s3_class(g, "net_frame")
    if (!is.na(info$nodes[i])) expect_equal(nrow(g$nodes), info$nodes[i], info = info$name[i])
    if (!is.na(info$edges[i])) expect_equal(nrow(g$edges), info$edges[i], info = info$name[i])
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
  expect_error(net_example("karate", relation = "x"), "multiplex examples")
})

test_that("two-mode examples are two-mode", {
  aff <- net_example("noordin_affiliations")
  expect_true(net_is_twomode(aff))
  expect_equal(sum(aff$nodes$mode == "event"), 45L)
  sw <- net_example("southern_women")
  expect_true(net_is_twomode(sw))
})

test_that("weights are opt-in and carry their meaning", {
  expect_null(net_example("jemaah_islamiyah")$roles$weight)
  ji <- net_example("jemaah_islamiyah", weighted = TRUE)
  expect_identical(ji$roles$weight, "strength")
  expect_true(all(ji$edges$strength %in% 1:5))
  expect_null(net_example("karate")$roles$weight)
  expect_identical(net_example("les_miserables", weighted = TRUE)$roles$weight_type, "strength")
  expect_identical(net_example("miles", weighted = TRUE)$roles$weight_type, "distance")
})

test_that("the new examples have the configurations DATASETS.md lists", {
  k <- net_example("krackhardt")
  expect_true(k$directed)
  expect_setequal(unique(k$edges$relation), c("advice", "friendship", "reports_to"))
  s <- net_example("sampson")
  expect_s3_class(s, "temporal_net_frame")
  expect_identical(s$roles$sign, "sign")
  expect_identical(s$roles$wave, "time")
  h <- net_example("highland_tribes")
  expect_false(h$directed)
  expect_equal(sum(h$edges$sign == -1), 29L)
  m <- net_example("miles")
  expect_s3_class(m, "spatial_net_frame")
  cb <- net_example("country_borders")
  expect_s3_class(cb, "spatial_net_frame")
  expect_gt(net_components(cb)$no, 1L)
  if (net_examples()$available[net_examples()$name == "workplace_contacts"]) {
    w <- net_example("workplace_contacts")
    expect_s3_class(w, "temporal_net_frame")
  } else {
    expect_error(net_example("workplace_contacts"), "not in this installation")
  }
})

test_that("the karate club matches igraph's copy", {
  skip_if_not_installed("igraph")
  k <- net_example("karate")
  z <- igraph::make_graph("Zachary")
  expect_equal(net_betweenness(k), igraph::betweenness(z))
})
