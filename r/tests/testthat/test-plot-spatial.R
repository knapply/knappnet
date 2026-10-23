# net_plot_spatial() and net_plot_edge_lengths() on networks with coordinate
# columns (R/plot-spatial.R); sf geometries are drawn in test-sf.R.

test_that("net_plot_spatial() draws point-based networks without sf", {
  nodes <- data.frame(id = c("London", "Paris", "Berlin", "Madrid"),
                      lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                      lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                      region = c("north", "north", "north", "south"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(source = c("London", "Paris", "Paris"), target = c("Paris", "Berlin", "Madrid"),
                      stringsAsFactors = FALSE)
  g <- net_spatial(net_frame(edges, nodes, directed = TRUE), radius = 6371.0088)
  on_null_device({
    expect_silent(out <- net_plot_spatial(g))
    expect_identical(out, g)
    expect_silent(net_plot_spatial(g, node_color = "region", node_size = "lat", main = "Cities"))
    expect_silent(net_plot_spatial(g, node_color = "red", node_label = FALSE, edge_width = 2,
                               edge_color = c("red", "blue", "green"), legend = FALSE))
    expect_silent(net_plot_spatial(net_spatial(net_frame(edges[0, ], nodes, directed = FALSE, quiet = TRUE), quiet = TRUE)))
    expect_error(net_plot_spatial(net_frame(edges, directed = FALSE)), "no node locations")
    expect_error(net_plot_spatial(g, node_color = "nope"), "color")
    expect_error(net_plot_spatial(g, node_size = 1:2), "node_size")
    expect_silent(net_plot_spatial(g, node_color = c("a", "b", "c", "d")))
    nine <- net_spatial(net_frame(data.frame(source = 1:9, target = c(2:9, 1)),
                                     data.frame(id = 1:9, x = cos(1:9), y = sin(1:9),
                                                cat = letters[1:9]),
                                     directed = FALSE), coords = c("x", "y"), crs = "planar")
    # Beyond eight categories the qualitative palette falls back to HCL
    # "Dark 3" without a warning.
    expect_silent(net_plot_spatial(nine, node_color = "cat"))
    expect_silent(out <- net_plot_edge_lengths(g))
    expect_identical(out, g)
    expect_silent(net_plot_edge_lengths(g, scale = 500, form = "gaussian", main = "Lengths"))
    expect_silent(net_plot_edge_lengths(net_spatial(net_frame(edges[0, ], nodes, directed = FALSE, quiet = TRUE), quiet = TRUE)))
    expect_equal(graphics::par("mar"), c(5.1, 4.1, 4.1, 2.1))
  })
})
