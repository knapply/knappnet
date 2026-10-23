# Heat maps of adjacency, incidence, and mixing matrices (R/plot-matrix.R).

karate <- function() net_example("karate")

test_that("net_plot_matrix() draws adjacency and incidence matrices in any order", {
  g <- karate()
  on_null_device({
    expect_silent(m <- net_plot_matrix(g))
    expect_identical(dim(m), c(34L, 34L))
    expect_identical(m, as.matrix(g))
    expect_silent(m2 <- net_plot_matrix(g, order = "faction", main = "Sorted"))
    expect_setequal(rownames(m2), rownames(m))
    # Rows are sorted by group in order of first appearance, stably.
    fac <- g$nodes$faction
    expect_identical(rownames(m2), as.character(g$nodes$id[order(match(fac, unique(fac)))]))
    expect_silent(net_plot_matrix(g, order = net_louvain(g, seed = 1), labels = FALSE))
    expect_silent(net_plot_matrix(g, values = TRUE, key = TRUE, labels = TRUE))
    expect_restores_par(net_plot_matrix(g))
    # Plain matrices, with and without dimnames and with a rectangular shape.
    expect_silent(net_plot_matrix(matrix(c(0, 1, 2, 0), 2)))
    expect_silent(net_plot_matrix(matrix(1:6, 2), order = c("b", "a")))
    expect_error(net_plot_matrix(matrix(1:4, 2), order = 1:3), "one value per row")
    expect_error(net_plot_matrix("x"), "net_frame or a numeric matrix")
    women <- net_example("southern_women")
    expect_silent(inc <- net_plot_matrix(women))
    expect_identical(inc, net_incidence(women))
    expect_silent(net_plot_matrix(women, order = rep(c("x", "y"), length.out = nrow(women$nodes))))
    # Weighted values fill the matrix.
    w <- net_example("karate", weighted = TRUE)
    expect_silent(wm <- net_plot_matrix(w, attr = "contexts"))
    expect_identical(wm, as.matrix(w, attr = "contexts"))
  })
})

test_that("net_plot_mixing() draws the mixing matrix of a network or a matrix", {
  g <- karate()
  on_null_device({
    expect_silent(m <- net_plot_mixing(g, "faction"))
    expect_identical(m, net_mixing(g, "faction"))
    expect_silent(net_plot_mixing(net_mixing(g, "faction"), main = "Factions"))
    expect_silent(net_plot_mixing(g, "faction", "id", values = FALSE))
    expect_error(net_plot_mixing(g), "dim1")
    expect_error(net_plot_mixing(list()), "net_frame or a numeric matrix")
  })
})
