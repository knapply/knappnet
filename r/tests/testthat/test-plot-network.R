# Overlays on the network plot (R/plot-network.R): centrality, partitions,
# k-core shells, ego networks, and robustness curves.

karate <- function() example_network("karate")

test_that("network overlays draw centrality, partitions, cores, and egos", {
  g <- karate()
  on_null_device({
    expect_silent(xy <- plot_centrality(g, betweenness(g), title = "betweenness", seed = 1))
    expect_identical(xy, layout_network(g, seed = 1))
    expect_silent(plot_centrality(g, degree(g), color = FALSE, layout = "circle"))
    g2 <- g
    g2$nodes$deg <- degree(g)
    expect_silent(plot_centrality(g2, "deg", seed = 1))
    expect_error(plot_centrality(g, 1:3), "one value per node")
    expect_error(plot_centrality(g, degree(g), title = 1), "single string")
    # The temporary legend column never touches the caller's network.
    expect_identical(names(g$nodes), c("name", "faction"))

    groups <- cluster_louvain(g, seed = 1)
    expect_silent(xy <- plot_partition(g, groups, seed = 1))
    expect_identical(xy, layout_network(g, seed = 1))
    expect_silent(plot_partition(g, groups, seed = 1, main = "Communities"))
    expect_silent(plot_partition(g, "faction", seed = 1))
    expect_silent(plot_partition(g, components(g)$membership, title = "component", layout = "circle"))
    expect_error(plot_partition(g, 1:3), "one value per node")

    expect_silent(xy <- plot_coreness(g, seed = 1, threads = 1))
    expect_identical(dim(xy), c(34L, 2L))
    k <- coreness(g)
    radius <- sqrt(rowSums(xy^2))
    # Higher cores sit on inner rings.
    expect_true(all(diff(tapply(radius, k, mean)) <= 0))
    expect_silent(plot_coreness(g, seed = 1, threads = 1, main = "Cores", node_label = FALSE))
    directed <- netframe(data.frame(from = c(1, 2, 3, 1), to = c(2, 3, 1, 4)))
    expect_silent(plot_coreness(directed, mode = "in", seed = 1, threads = 1))
    expect_restores_par(plot_coreness(g, seed = 1, threads = 1))

    f <- example_network("florentine")
    expect_silent(e <- plot_ego(f, "Medici", seed = 1))
    expect_s3_class(e, "netframe")
    expect_identical(sort(e$nodes$name), sort(ego_network(f, "Medici")$nodes$name))
    expect_silent(plot_ego(f, "Medici", order = 2, seed = 1, main = "Two steps"))
    expect_error(plot_ego(f, "Nobody"), "not in the network")
  })
})

test_that("plot_partition() takes the result of a grouping routine", {
  # Two triangles {1, 2, 3} and {4, 5, 6} joined by the bridge 3 - 4, with
  # a pendant 7 on 6.
  g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5, 6), to = c(2, 3, 3, 4, 5, 6, 6, 7)),
                directed = FALSE)
  f <- factions(g, k = 2, seed = 1, threads = 2)
  on_null_device({
    expect_invisible(plot_partition(g, f, seed = 1))
    expect_invisible(plot_partition(g, letters[c(1, 1, 2, 2, 3, 3, 3)], seed = 1, main = "x"))
    expect_error(plot_partition(g, 1:3), "one value per node")
  })
})

test_that("plot_robustness() tracks the largest component as nodes are removed", {
  g <- example_network("jemaah_islamiyah")
  n <- nrow(g$nodes)
  on_null_device({
    expect_silent(out <- plot_robustness(g, steps = 5, seed = 1, threads = 1))
    expect_setequal(unique(out$curve), c("degree", "betweenness", "random"))
    first <- out[out$removed == 0, ]
    expect_true(all(first$value == max(components(g)$csize) / n))
    expect_true(all(out$value >= 0 & out$value <= 1))
    # Removing the top-degree node by hand gives the same value as the curve.
    top <- g$nodes$name[which.max(degree(g))]
    expect_silent(by_hand <- plot_robustness(g, nodes = top, strategy = NULL, threads = 1))
    expect_identical(by_hand$curve, rep("selected nodes", 2))
    rest <- netframe(g$edges[g$edges$from != top & g$edges$to != top, ],
                     g$nodes[g$nodes$name != top, ], directed = FALSE)
    expect_equal(by_hand$value[2], max(components(rest)$csize) / n)
    kp <- key_players(g, 2, "fragmentation", seed = 1, threads = 1)
    expect_silent(out <- plot_robustness(g, nodes = kp$set, strategy = "degree",
                                         measure = "fragmentation", label = "key players",
                                         steps = 4, threads = 1, main = "Attack"))
    expect_setequal(unique(out$curve), c("degree", "key players"))
    expect_equal(out$value[out$curve == "key players"][3], kp$score)
    expect_restores_par(plot_robustness(g, strategy = "degree", steps = 3, threads = 1))
    expect_error(plot_robustness(g, strategy = "closeness"), "`strategy`")
    expect_error(plot_robustness(g, nodes = NULL, strategy = NULL), "Give `nodes`")
    expect_error(plot_robustness(g, nodes = "x", label = "degree"), "not the name of a strategy")
    expect_error(plot_robustness(g, steps = 1), "`steps`")
    expect_error(plot_robustness(netframe(data.frame(from = 1, to = 1))), "at least two nodes")
  })
})
