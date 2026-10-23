# Heat maps of adjacency, incidence, and mixing matrices (R/plot-matrix.R).

karate <- function() example_network("karate")

test_that("plot_matrix() draws adjacency and incidence matrices in any order", {
  g <- karate()
  on_null_device({
    expect_silent(m <- plot_matrix(g))
    expect_identical(dim(m), c(34L, 34L))
    expect_identical(m, as.matrix(g))
    expect_silent(m2 <- plot_matrix(g, order = "faction", main = "Sorted"))
    expect_setequal(rownames(m2), rownames(m))
    # Rows are sorted by group in order of first appearance, stably.
    fac <- g$nodes$faction
    expect_identical(rownames(m2), as.character(g$nodes$name[order(match(fac, unique(fac)))]))
    expect_silent(plot_matrix(g, order = cluster_louvain(g, seed = 1), labels = FALSE))
    expect_silent(plot_matrix(g, values = TRUE, key = TRUE, labels = TRUE))
    expect_restores_par(plot_matrix(g))
    # Plain matrices, with and without dimnames and with a rectangular shape.
    expect_silent(plot_matrix(matrix(c(0, 1, 2, 0), 2)))
    expect_silent(plot_matrix(matrix(1:6, 2), order = c("b", "a")))
    expect_error(plot_matrix(matrix(1:4, 2), order = 1:3), "one value per row")
    expect_error(plot_matrix("x"), "netframe or a numeric matrix")
    women <- example_network("southern_women")
    expect_silent(inc <- plot_matrix(women))
    expect_identical(inc, incidence_matrix(women))
    expect_silent(plot_matrix(women, order = rep(c("x", "y"), length.out = nrow(women$nodes))))
    # Weighted values fill the matrix.
    w <- example_network("karate", weighted = TRUE)
    expect_silent(wm <- plot_matrix(w, attr = "contexts"))
    expect_identical(wm, as.matrix(w, attr = "contexts"))
  })
})

test_that("plot_mixing_matrix() draws the mixing matrix of a network or a matrix", {
  g <- karate()
  on_null_device({
    expect_silent(m <- plot_mixing_matrix(g, "faction"))
    expect_identical(m, mixing_matrix(g, "faction"))
    expect_silent(plot_mixing_matrix(mixing_matrix(g, "faction"), main = "Factions"))
    expect_silent(plot_mixing_matrix(g, "faction", "name", values = FALSE))
    expect_error(plot_mixing_matrix(g), "dim1")
    expect_error(plot_mixing_matrix(list()), "netframe or a numeric matrix")
  })
})
