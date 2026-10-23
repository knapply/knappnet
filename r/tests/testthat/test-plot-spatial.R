# plot_spatial() and plot_edge_lengths() on networks with coordinate
# columns (R/plot-spatial.R); sf geometries are drawn in test-sf.R.

test_that("plot_spatial() draws point-based networks without sf", {
  nodes <- data.frame(name = c("London", "Paris", "Berlin", "Madrid"),
                      lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                      lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                      region = c("north", "north", "north", "south"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(from = c("London", "Paris", "Paris"), to = c("Paris", "Berlin", "Madrid"),
                      stringsAsFactors = FALSE)
  g <- spatial_network(netframe(edges, nodes, directed = TRUE), radius = 6371.0088)
  on_null_device({
    expect_silent(out <- plot_spatial(g))
    expect_identical(out, g)
    expect_silent(plot_spatial(g, node_color = "region", node_size = "lat", main = "Cities"))
    expect_silent(plot_spatial(g, node_color = "red", node_label = FALSE, edge_width = 2,
                               edge_color = c("red", "blue", "green"), legend = FALSE))
    expect_silent(plot_spatial(spatial_network(netframe(edges[0, ], nodes, directed = FALSE))))
    expect_error(plot_spatial(netframe(edges, directed = FALSE)), "no node locations")
    expect_error(plot_spatial(g, node_color = "nope"), "color")
    expect_error(plot_spatial(g, node_size = 1:2), "node_size")
    expect_silent(plot_spatial(g, node_color = c("a", "b", "c", "d")))
    nine <- spatial_network(netframe(data.frame(from = 1:9, to = c(2:9, 1)),
                                     data.frame(id = 1:9, x = cos(1:9), y = sin(1:9),
                                                cat = letters[1:9]),
                                     directed = FALSE), coords = c("x", "y"), crs = "planar")
    # Beyond eight categories the qualitative palette falls back to HCL
    # "Dark 3" without a warning.
    expect_silent(plot_spatial(nine, node_color = "cat"))
    expect_silent(out <- plot_edge_lengths(g))
    expect_identical(out, g)
    expect_silent(plot_edge_lengths(g, scale = 500, form = "gaussian", main = "Lengths"))
    expect_silent(plot_edge_lengths(spatial_network(netframe(edges[0, ], nodes, directed = FALSE))))
    expect_equal(graphics::par("mar"), c(5.1, 4.1, 4.1, 2.1))
  })
})
