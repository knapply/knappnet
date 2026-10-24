
test_that("the cached engine graph notices edited endpoints", {
  g <- net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "d")))
  expect_equal(unname(net_degree(g)), c(1, 2, 2, 1))
  g2 <- g
  g2$edges$target[3] <- "a"
  expect_error(net_degree(g2), "rebuild|edited")
  # The original, which shares the cache, is unaffected.
  expect_equal(unname(net_degree(g)), c(1, 2, 2, 1))
})
