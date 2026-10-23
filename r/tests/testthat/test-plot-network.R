# Overlays on the network plot (R/plot-network.R): centrality, partitions,
# k-core shells, ego networks, and robustness curves.

karate <- function() net_example("karate")

test_that("network overlays draw centrality, partitions, cores, and egos", {
  g <- karate()
  on_null_device({
    expect_silent(xy <- net_plot_centrality(g, net_betweenness(g), title = "betweenness", seed = 1))
    expect_identical(xy, net_layout(g, seed = 1))
    expect_silent(net_plot_centrality(g, net_degree(g), color = FALSE, layout = "circle"))
    g2 <- g
    g2$nodes$deg <- net_degree(g)
    expect_silent(net_plot_centrality(g2, "deg", seed = 1))
    expect_error(net_plot_centrality(g, 1:3), "one value per node")
    expect_error(net_plot_centrality(g, net_degree(g), title = 1), "single string")
    # The temporary legend column never touches the caller's network.
    expect_identical(names(g$nodes), c("id", "faction"))

    groups <- net_louvain(g, seed = 1)
    expect_silent(xy <- net_plot_partition(g, groups, seed = 1))
    expect_identical(xy, net_layout(g, seed = 1))
    expect_silent(net_plot_partition(g, groups, seed = 1, main = "Communities"))
    expect_silent(net_plot_partition(g, "faction", seed = 1))
    expect_silent(net_plot_partition(g, net_components(g)$membership, title = "component", layout = "circle"))
    expect_error(net_plot_partition(g, 1:3), "one value per node")

    expect_silent(xy <- net_plot_coreness(g, seed = 1, threads = 1))
    expect_identical(dim(xy), c(34L, 2L))
    k <- net_coreness(g)
    radius <- sqrt(rowSums(xy^2))
    # Higher cores sit on inner rings.
    expect_true(all(diff(tapply(radius, k, mean)) <= 0))
    expect_silent(net_plot_coreness(g, seed = 1, threads = 1, main = "Cores", node_label = FALSE))
    directed <- net_frame(data.frame(source = c(1, 2, 3, 1), target = c(2, 3, 1, 4)))
    expect_silent(net_plot_coreness(directed, mode = "in", seed = 1, threads = 1))
    expect_restores_par(net_plot_coreness(g, seed = 1, threads = 1))

    f <- net_example("florentine")
    expect_silent(e <- net_plot_ego(f, "Medici", seed = 1))
    expect_s3_class(e, "net_frame")
    expect_identical(sort(e$nodes$id), sort(net_ego(f, "Medici")$nodes$id))
    expect_silent(net_plot_ego(f, "Medici", order = 2, seed = 1, main = "Two steps"))
    expect_error(net_plot_ego(f, "Nobody"), "not in the network")
  })
})

test_that("net_plot_partition() takes the result of a grouping routine", {
  # Two triangles {1, 2, 3} and {4, 5, 6} joined by the bridge 3 - 4, with
  # a pendant 7 on 6.
  g <- net_frame(data.frame(source = c(1, 1, 2, 3, 4, 4, 5, 6), target = c(2, 3, 3, 4, 5, 6, 6, 7)),
                directed = FALSE)
  f <- net_factions(g, k = 2, seed = 1, threads = 2)
  on_null_device({
    expect_invisible(net_plot_partition(g, f, seed = 1))
    expect_invisible(net_plot_partition(g, letters[c(1, 1, 2, 2, 3, 3, 3)], seed = 1, main = "x"))
    expect_error(net_plot_partition(g, 1:3), "one value per node")
  })
})

test_that("net_plot_robustness() tracks the largest component as nodes are removed", {
  g <- net_example("jemaah_islamiyah")
  n <- nrow(g$nodes)
  on_null_device({
    expect_silent(out <- net_plot_robustness(g, steps = 5, seed = 1, threads = 1))
    expect_setequal(unique(out$curve), c("degree", "betweenness", "random"))
    first <- out[out$removed == 0, ]
    expect_true(all(first$value == max(net_components(g)$csize) / n))
    expect_true(all(out$value >= 0 & out$value <= 1))
    # Removing the top-degree node by hand gives the same value as the curve.
    top <- g$nodes$id[which.max(net_degree(g))]
    expect_silent(by_hand <- net_plot_robustness(g, nodes = top, strategy = NULL, threads = 1))
    expect_identical(by_hand$curve, rep("selected nodes", 2))
    rest <- net_frame(g$edges[g$edges$source != top & g$edges$target != top, ],
                     g$nodes[g$nodes$id != top, ], directed = FALSE)
    expect_equal(by_hand$value[2], max(net_components(rest)$csize) / n)
    kp <- net_key_players(g, 2, "fragmentation", seed = 1, threads = 1)
    expect_silent(out <- net_plot_robustness(g, nodes = kp$set, strategy = "degree",
                                         measure = "fragmentation", label = "key players",
                                         steps = 4, threads = 1, main = "Attack"))
    expect_setequal(unique(out$curve), c("degree", "key players"))
    expect_equal(out$value[out$curve == "key players"][3], kp$score)
    expect_restores_par(net_plot_robustness(g, strategy = "degree", steps = 3, threads = 1))
    expect_error(net_plot_robustness(g, strategy = "closeness"), "`strategy`")
    expect_error(net_plot_robustness(g, nodes = NULL, strategy = NULL), "Give `nodes`")
    expect_error(net_plot_robustness(g, nodes = "x", label = "degree"), "not the name of a strategy")
    expect_error(net_plot_robustness(g, steps = 1), "`steps`")
    expect_error(net_plot_robustness(net_frame(data.frame(source = 1, target = 1))), "at least two nodes")
  })
})
