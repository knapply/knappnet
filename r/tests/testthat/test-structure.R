# Density, average degree and the degree distribution, components and
# fragmentation, the triad census, centralization, compactness, the
# small-world quotient, and core-periphery structure (R/structure.R).
# Hand-checked values first, then spot checks against igraph and sna. The
# full comparison lives in validation/.

# Node 2 ties to 1, 4, 5, 6, 7; 5 ties to 2, 3, 4, 7; 7 ties to 1, 2, 5.
# The core {2, 5, 7} is a triangle and the periphery has no ties among
# itself, so the fit with core-periphery pairs left out is perfect.
cp_net <- function() {
  net_frame(data.frame(source = c(1, 1, 2, 2, 2, 2, 3, 4, 5), target = c(2, 7, 4, 5, 6, 7, 5, 5, 7)),
           directed = FALSE)
}

test_that("density counts edges over possible edges", {
  g <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 4)), directed = FALSE)
  expect_equal(net_density(g), 0.5)
  expect_equal(net_density(g, loops = TRUE), 0.3)
  multi <- net_frame(data.frame(source = c(1, 1, 2, 3), target = c(2, 2, 3, 3)))
  expect_equal(net_density(multi), 4 / 6)
  expect_equal(net_density(multi, convention = "sna"), 2 / 6)
  expect_error(net_density(multi, loops = TRUE, convention = "sna"), "self-loops")
  one <- net_frame(data.frame(source = 1, target = 1))
  expect_true(is.nan(net_density(one)))
  expect_equal(net_density(one, convention = "networkx"), 0)
})

test_that("average degree and degree distribution", {
  g <- net_frame(data.frame(source = c(1, 1, 1, 2), target = c(2, 3, 4, 3)), directed = FALSE)
  expect_equal(net_avg_degree(g), 2)
  expect_equal(net_avg_degree(g, mode = "out"), 2)
  dd <- net_degree_dist(g, cumulative = TRUE)
  expect_equal(dd$degree, c(1, 2, 3))
  expect_equal(dd$count, c(1, 2, 1))
  expect_equal(dd$proportion, c(0.25, 0.5, 0.25))
  expect_equal(dd$at_least, c(1, 0.75, 0.25))
  expect_null(net_degree_dist(g)$at_least)
  d <- net_frame(data.frame(source = c(1, 1), target = c(2, 3)))
  expect_equal(net_avg_degree(d, mode = "out"), 2 / 3)
  expect_equal(net_avg_degree(d, mode = "in"), 2 / 3)
  expect_equal(net_avg_degree(d), 4 / 3)
  empty <- net_frame(data.frame(source = integer(), target = integer()),
                    nodes = data.frame(id = integer()))
  expect_true(is.nan(net_avg_degree(empty)))
})

test_that("components and fragmentation", {
  g <- pieces()
  cc <- net_components(g)
  expect_equal(cc$membership, c(1L, 1L, 1L, 2L, 2L, 3L))
  expect_equal(cc$csize, c(3L, 2L, 1L))
  expect_equal(cc$no, 3L)
  # Ordered pairs in the same component: 6 + 2 = 8 of 30.
  expect_equal(net_fragmentation(g), 22 / 30)
  d <- net_frame(data.frame(source = c(1, 2, 3, 3), target = c(2, 3, 1, 4)))
  expect_equal(net_components(d, "strong")$membership, c(1L, 1L, 1L, 2L))
  expect_equal(net_components(d, "weak")$no, 1L)
  expect_equal(net_fragmentation(d, "strong"), 6 / 12)
})

test_that("triad census matches igraph and sna", {
  skip_if_not_installed("igraph")
  set.seed(1)
  e <- data.frame(source = sample(30, 120, TRUE), target = sample(30, 120, TRUE))
  e <- e[e$source != e$target, ]
  e <- e[!duplicated(e), ]
  g <- net_frame(e, nodes = data.frame(id = 1:30))
  ig <- igraph::graph_from_data_frame(e, vertices = data.frame(id = 1:30))
  expect_equal(unname(net_triad_census(g)), igraph::triad_census(ig))

  skip_if_not_installed("sna")
  m <- matrix(0, 30, 30)
  m[cbind(e$source, e$target)] <- 1
  expect_equal(unname(net_triad_census(g)), unname(drop(sna::triad.census(m))))
})

test_that("centralization of a star is 1 under every convention", {
  star <- net_frame(data.frame(source = 1, target = 2:6), directed = FALSE)
  for (conv in c("igraph", "sna", "ucinet")) {
    for (m in c("degree", "betweenness", "closeness")) {
      expect_equal(net_centralization(star, m, loops = FALSE, convention = conv), 1,
                   tolerance = 1e-12, info = paste(conv, m))
    }
  }
  expect_error(net_centralization(star, convention = "networkx"), "NetworkX")
  ring <- net_frame(data.frame(source = 1:5, target = c(2:5, 1)), directed = FALSE)
  expect_equal(net_centralization(ring), 0)
})

test_that("eigenvector centralization follows igraph's definition", {
  star10 <- net_frame(data.frame(source = 1, target = 2:10), directed = FALSE)
  expect_equal(net_centralization(star10, "eigenvector"), 0.75, tolerance = 1e-9)
  expect_equal(net_centralization(star10, "eigenvector", normalized = FALSE), 6, tolerance = 1e-9)
  one <- net_frame(data.frame(source = 1, target = 2), nodes = data.frame(id = 1:10),
                  directed = FALSE)
  expect_equal(net_centralization(one, "eigenvector"), 1, tolerance = 1e-9)
  cycle <- net_frame(data.frame(source = 1:5, target = c(2:5, 1)), directed = FALSE)
  expect_equal(net_centralization(cycle, "eigenvector"), 0, tolerance = 1e-9)
  expect_true(is.nan(net_centralization(net_frame(data.frame(source = 1, target = 2), directed = FALSE),
                                    "eigenvector")))
  expect_error(net_centralization(star10, "eigenvector", convention = "sna"), "igraph")
  expect_error(net_centralization(star10, "eigenvector", convention = "ucinet"), "igraph")
  d <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)))
  expect_error(net_centralization(d, "eigenvector", mode = "in"), "direction ignored")
  two <- net_frame(data.frame(source = c(1, 3), target = c(2, 4)), directed = FALSE)
  expect_warning(net_centralization(two, "eigenvector"), "component")
})

test_that("net_centralize() applies Freeman's index to any scores", {
  star <- net_frame(data.frame(source = 1, target = 2:5), directed = FALSE)
  expect_equal(net_centralize(net_degree(star), 4 * 3), 1)
  # Freeman's maximum (n - 1)(n - 2) is igraph's without self-loops.
  expect_equal(net_centralize(net_degree(star), 4 * 3), net_centralization(star, loops = FALSE))
  expect_equal(net_centralize(net_degree(star)), 12)
  expect_equal(net_centralize(c(2, 2, 2)), 0)
  expect_true(is.nan(net_centralize(numeric(0))))
  expect_true(is.nan(net_centralize(c(1, NaN))))
  expect_error(net_centralize("a"), "numeric")
  expect_error(net_centralize(1:3, theoretical_max = 0), "positive")
})

test_that("compactness and distance-weighted fragmentation", {
  g <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)), nodes = data.frame(id = 1:4),
                directed = FALSE)
  expect_equal(net_compactness(g), (4 + 2 * 0.5) / 12)
  expect_equal(net_dw_fragmentation(g), 1 - 5 / 12)
  d <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)))
  expect_equal(net_compactness(d), (2 + 0.5) / 6)
  expect_equal(net_compactness(d, mode = "in"), (2 + 0.5) / 6)
  expect_equal(net_compactness(d, mode = "all"), 5 / 6)
  w <- net_frame(data.frame(source = c(1, 2), target = c(2, 3), w = c(2, 2)), directed = FALSE,
                weight = "w", weight_type = "distance")
  expect_equal(net_compactness(w), (4 / 2 + 2 / 4) / 6)
  expect_equal(net_compactness(w, weighted = FALSE), 5 / 6)
  expect_equal(net_compactness(w, weights_as = "strength"), (4 * 2 + 2 * 1) / 6)
  one <- net_frame(data.frame(source = integer(), target = integer()), nodes = data.frame(id = 1))
  expect_true(is.nan(net_compactness(one)))
})

test_that("small-world quotient is a ratio of ratios", {
  g <- net_example("karate")
  sw <- net_small_world(g, reps = 20, seed = 1, threads = 2)
  expect_s3_class(sw, "small_world")
  expect_equal(sw$clustering, net_transitivity(g))
  expect_equal(sw$mean_distance, net_mean_distance(g))
  expect_equal(sw$quotient,
               (sw$clustering / sw$random_clustering) /
                 (sw$mean_distance / sw$random_mean_distance))
  expect_true(sw$quotient > 1)
  expect_equal(net_small_world(g, reps = 20, seed = 1, threads = 1), sw)
  local <- net_small_world(g, local = TRUE, reps = 20, seed = 1, threads = 2)
  # Nodes with fewer than two neighbours count as 0, as in NetworkX.
  expect_equal(local$clustering, net_avg_clustering(g, convention = "networkx"))
  expect_false(isTRUE(all.equal(local$quotient, sw$quotient)))
  expect_output(print(sw), "quotient")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(sw))
  expect_error(net_small_world(g, reps = 0), "at least 1")
})

test_that("discrete core-periphery fits by every link rule", {
  g <- cp_net()
  cp <- net_core_periphery(g, seed = 1, threads = 2)
  expect_s3_class(cp, "core_periphery")
  expect_equal(g$nodes$name[cp$core], c(2, 5, 7))
  expect_equal(cp$correlation, 1)
  expect_equal(net_core_periphery(g, seed = 1, threads = 1)$core, cp$core)
  # Exhaustive optima for the other rules (see validation/).
  expect_equal(g$nodes$name[net_core_periphery(g, links = "ones", seed = 1, threads = 2)$core],
               c(2, 5))
  expect_equal(g$nodes$name[net_core_periphery(g, links = "zeros", seed = 1, threads = 2)$core],
               c(1, 2, 4, 5, 7))
  expect_equal(net_core_periphery(g, links = "ones", seed = 1, threads = 2)$correlation,
               0.6330542, tolerance = 1e-6)
  expect_output(print(cp), "Correlation")
  # No structure to fit.
  complete <- net_frame(data.frame(source = c(1, 1, 2), target = c(2, 3, 3)), directed = FALSE)
  expect_error(net_core_periphery(complete, seed = 1), "different values")
})

test_that("continuous core-periphery fits", {
  g <- cp_net()
  cc <- net_core_periphery(g, method = "continuous", threads = 2)
  expect_equal(cc$method, "continuous")
  expect_equal(cc$fit, "minres")
  expect_true(all(cc$coreness >= 0))
  expect_equal(sum(cc$coreness^2), 1, tolerance = 1e-9)
  expect_equal(cc$correlation, 0.6957548, tolerance = 1e-6)
  expect_equal(cc$residual, 2.669288, tolerance = 1e-6)
  expect_equal(which.max(cc$coreness), 2)
  # The reported correlation is that of the pattern c_i c_j, and the
  # residual that of the least-squares scaling of the pattern.
  a <- as.matrix(g)
  iu <- upper.tri(a)
  pattern <- outer(cc$coreness, cc$coreness)
  expect_equal(cor(a[iu], pattern[iu]), cc$correlation, tolerance = 1e-8)
  scale <- sum(a[iu] * pattern[iu]) / sum(pattern[iu]^2)
  expect_equal(sum((a[iu] - scale * pattern[iu])^2), cc$residual, tolerance = 1e-6)
  expect_output(print(cc), "MINRES")

  # The correlation fit: also unit norm, with its own correlation.
  cr <- net_core_periphery(g, method = "continuous", fit = "correlation", seed = 1, threads = 2)
  expect_equal(cr$fit, "correlation")
  expect_true(is.nan(cr$residual))
  expect_equal(sum(cr$coreness^2), 1, tolerance = 1e-9)
  pattern <- outer(cr$coreness, cr$coreness)
  expect_equal(cor(a[iu], pattern[iu]), cr$correlation, tolerance = 1e-6)
  expect_true(cr$correlation >= cc$correlation - 1e-9)
  expect_equal(net_core_periphery(g, method = "continuous", fit = "correlation", seed = 1,
                              threads = 1)$coreness, cr$coreness)

  # A star: the center is the core and the leaves are symmetric under both.
  star <- net_frame(data.frame(source = 1, target = 2:6), directed = FALSE)
  for (fit in c("minres", "correlation")) {
    s <- net_core_periphery(star, method = "continuous", fit = fit, seed = 1, threads = 2)
    expect_equal(which.max(s$coreness), 1)
    expect_equal(diff(range(s$coreness[-1])), 0, tolerance = 1e-6)
    expect_equal(s$correlation, 1, tolerance = 1e-6)
  }
  expect_error(net_core_periphery(g, method = "continuous", tol = -1), "tol")
})

test_that("core-periphery plots draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- cp_net()
  cp <- net_core_periphery(g, seed = 1, threads = 2)
  cc <- net_core_periphery(g, method = "continuous", threads = 2)
  expect_invisible(plot(cp))
  expect_invisible(plot(cp, what = "network", seed = 1))
  expect_invisible(plot(cc, seed = 1))
  expect_invisible(plot(cc, what = "matrix", main = "coreness"))
  expect_error(plot(cp, what = "x"), "what")
})

test_that("density, components, and centralization match igraph", {
  skip_if_not_installed("igraph")
  g <- kite()
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  expect_equal(net_density(g), igraph::edge_density(ig))
  expect_equal(net_components(g)$membership, unname(igraph::components(ig)$membership))
  expect_equal(net_centralization(g), igraph::centr_degree(ig, loops = TRUE)$centralization)
  expect_equal(net_centralization(g, "betweenness"), igraph::centr_betw(ig)$centralization)
  expect_equal(net_centralization(g, "closeness"), igraph::centr_clo(ig)$centralization)
  expect_equal(net_centralization(g, "eigenvector"), igraph::centr_eigen(ig)$centralization,
               tolerance = 1e-8)

  e <- random_digraph(25, 80, 3)
  nodes <- data.frame(id = 1:25)
  d <- net_frame(e, nodes)
  id <- igraph::graph_from_data_frame(e, vertices = nodes)
  expect_equal(net_components(d, "strong")$no, igraph::components(id, "strong")$no)
  expect_equal(net_centralization(d, mode = "in"),
               igraph::centr_degree(id, mode = "in", loops = TRUE)$centralization)
})

test_that("density and centralization match sna", {
  skip_if_not_installed("sna")
  g <- kite()
  m <- as.matrix(g)
  expect_equal(net_density(g, convention = "sna"), unname(sna::gden(m, mode = "graph")))
  expect_equal(net_centralization(g, "degree", convention = "sna"),
               sna::centralization(m, sna::degree, mode = "graph"))
  expect_equal(net_centralization(g, "closeness", convention = "sna"),
               sna::centralization(m, sna::closeness, mode = "graph"))
})
