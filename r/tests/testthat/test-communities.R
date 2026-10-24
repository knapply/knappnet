# Community detection beyond Louvain (R/communities.R): Walktrap, label
# propagation, Infomap and the map equation, Leiden, spinglass, exact
# maximum modularity, fluid and Voronoi communities, partition comparison,
# and dendrograms. Hand-checked values and invariants first, then igraph
# where it has a deterministic counterpart. The full comparison with the
# reference tools lives in validation/.

# Two triangles {a, b, c} and {d, e, f} joined by the edge c - d.
comm_barbell <- function() {
  net_frame(data.frame(source = c("a", "a", "b", "c", "d", "d", "e"),
                       target = c("b", "c", "c", "d", "e", "f", "f")),
            directed = FALSE)
}

# Two disjoint triangles {1, 2, 3} and {4, 5, 6}.
comm_triangles <- function() {
  net_frame(data.frame(source = c(1, 1, 2, 4, 4, 5), target = c(2, 3, 3, 5, 6, 6)),
            nodes = data.frame(id = 1:6), directed = FALSE)
}

comm_groups <- function(memb) as.vector(memb)

test_that("every method returns a membership aligned with the nodes", {
  g <- comm_barbell()
  results <- list(
    net_walktrap(g, threads = 1),
    net_label_prop(g, seed = 1),
    net_infomap(g, seed = 1),
    net_leiden(g, seed = 1),
    net_spinglass(g, seed = 1),
    net_optimal(g),
    net_fluid(g, 2, seed = 1),
    net_voronoi(g, threads = 1)
  )
  for (r in results) {
    expect_type(r, "integer")
    expect_length(r, 6L)
    expect_identical(sort(unique(as.vector(r))), seq_len(max(r)))
    # Numbered in order of first appearance.
    expect_identical(as.vector(r), match(as.vector(r), unique(as.vector(r))))
  }
})

test_that("results follow the order of the node table", {
  e <- data.frame(source = c(1, 1, 2, 4, 4, 5), target = c(2, 3, 3, 5, 6, 6))
  g <- net_frame(e, nodes = data.frame(id = c(3, 1, 5, 2, 6, 4)), directed = FALSE)
  for (r in list(net_infomap(g, seed = 1), net_label_prop(g, seed = 1),
                 net_leiden(g, seed = 1), net_optimal(g))) {
    ids <- g$nodes$id
    expect_equal(length(unique(r[ids %in% 1:3])), 1L)
    expect_equal(length(unique(r[ids %in% 4:6])), 1L)
    expect_false(r[ids == 1] == r[ids == 4])
  }
})

# ---- Walktrap ---------------------------------------------------------------

test_that("Walktrap reports its dendrogram and the modularity of every level", {
  g <- comm_barbell()
  w <- net_walktrap(g, threads = 1)
  merges <- attr(w, "merges")
  expect_identical(dim(merges), c(5L, 2L))
  expect_identical(attr(merges, "leaves"), 6L)
  # Row i creates cluster 6 + i, and each cluster is joined at most once.
  expect_true(all(merges[, 1L] < 6L + seq_len(5L)) && all(merges[, 2L] < 6L + seq_len(5L)))
  expect_false(anyDuplicated(as.vector(merges)) > 0L)
  path <- attr(w, "modularity_path")
  expect_length(path, 6L)
  expect_equal(path[[1L]], net_modularity(g, 1:6))
  expect_equal(path[[6L]], 0)
  expect_equal(attr(w, "modularity"), max(path))
  expect_equal(attr(w, "modularity"), net_modularity(g, w))
  # Cutting the dendrogram at the best level reproduces the partition.
  best <- net_cut_dendrogram(merges, x = g)
  expect_equal(attr(best, "modularity"), attr(w, "modularity"))
  expect_equal(attr(best, "modularity_path"), path)
})

test_that("Walktrap thread count, weights, and components", {
  g <- net_example("karate")
  expect_equal(net_walktrap(g, threads = 1), net_walktrap(g, threads = 2))
  # Unit weights give the unweighted result.
  e <- g$edges
  e$w <- 1
  gw <- net_frame(e, g$nodes, directed = FALSE, weight = "w", weight_type = "strength")
  a <- net_walktrap(g, threads = 1)
  b <- net_walktrap(gw, threads = 1)
  expect_equal(comm_groups(b), comm_groups(a))
  expect_equal(attr(b, "merges"), attr(a, "merges"))
  expect_equal(attr(b, "modularity"), attr(a, "modularity"))
  # Two components: n - 2 merges, an incomplete dendrogram.
  w <- net_walktrap(comm_triangles(), threads = 1)
  expect_identical(nrow(attr(w, "merges")), 4L)
  expect_equal(comm_groups(w), c(1, 1, 1, 2, 2, 2))
  expect_error(net_cut_dendrogram(attr(w, "merges"), k = 1), "net_cut_dendrogram")
  expect_equal(comm_groups(net_cut_dendrogram(attr(w, "merges"), k = 2)), c(1, 1, 1, 2, 2, 2))
})

test_that("Walktrap checks its arguments", {
  g <- comm_barbell()
  expect_error(net_walktrap(net_frame(data.frame(source = 1, target = 2))), "undirected")
  expect_error(net_walktrap(g, steps = 0), "`steps`")
  expect_error(net_walktrap(g, steps = 2.5), "`steps`")
  expect_error(net_walktrap(g, weighted = TRUE), "no weights")
  expect_error(net_walktrap(g, threads = -1), "`threads`")
  expect_error(net_walktrap(1), "net_frame")
})

# ---- Label propagation --------------------------------------------------------

test_that("label propagation ends in a fixed point on cliques, for every variant", {
  g <- comm_triangles()
  for (v in c("dominance", "retention", "fast", "semi_synchronous")) {
    r <- net_label_prop(g, v, seed = 3)
    expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2), info = v)
  }
  # The semi-synchronous variant uses no random numbers.
  b <- comm_barbell()
  expect_identical(net_label_prop(b, "semi_synchronous", seed = 1),
                   net_label_prop(b, "semi_synchronous", seed = 99))
  expect_identical(net_label_prop(b, seed = 5), net_label_prop(b, seed = 5))
  set.seed(1)
  a <- net_label_prop(b)
  set.seed(1)
  expect_identical(net_label_prop(b), a)
})

test_that("label propagation takes starting labels and fixed nodes", {
  g <- comm_triangles()
  start <- c("x", NA, NA, "y", NA, NA)
  r <- net_label_prop(g, initial = start, fixed = !is.na(start), seed = 1)
  expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2))
  g$nodes$start <- start
  expect_equal(comm_groups(net_label_prop(g, initial = "start", seed = 1)), c(1, 1, 1, 2, 2, 2))
  expect_error(net_label_prop(g, initial = 1:3), "one value per node")
  expect_error(net_label_prop(g, fixed = c(TRUE, NA, rep(FALSE, 4))), "`fixed`")
  expect_error(net_label_prop(g, fixed = 1:6), "`fixed`")
})

test_that("label propagation reads directions according to `mode`", {
  d <- net_frame(data.frame(source = c(1, 2, 3, 4), target = c(2, 3, 1, 3)))
  for (m in c("all", "out", "in")) {
    expect_length(net_label_prop(d, mode = m, seed = 1), 4L)
  }
  expect_error(net_label_prop(d, mode = "both"), "should be one of")
  expect_error(net_label_prop(d, "semi_synchronous"), "undirected")
  expect_length(net_label_prop(d, "fast", seed = 1), 4L)
})

# ---- Infomap and the map equation -----------------------------------------------

test_that("the map equation of small networks follows the definition", {
  # A cycle of four nodes as one module costs the entropy of the visit
  # rates, log2(4) = 2 bits per step.
  cycle <- net_frame(data.frame(source = 1:4, target = c(2:4, 1)), directed = FALSE)
  expect_equal(net_map_equation(cycle, rep(1, 4)), 2)
  # Two disjoint edges as two modules: the walker never leaves a module, so
  # only the position within a module (one bit) is coded.
  pairs <- net_frame(data.frame(source = c(1, 3), target = c(2, 4)), directed = FALSE)
  expect_equal(net_map_equation(pairs, c(1, 1, 2, 2)), 1)
  expect_equal(net_map_equation(pairs, c("a", "a", "b", "b")), 1)
  expect_gt(net_map_equation(pairs, c(1, 2, 3, 4)), 1)
  expect_error(net_map_equation(pairs, 1:3), "one value per node")
  expect_error(net_map_equation(pairs, c(1, 1, NA, 2)), "missing")
})

test_that("Infomap minimizes the map equation it reports", {
  g <- comm_triangles()
  r <- net_infomap(g, seed = 1)
  expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2))
  expect_equal(attr(r, "codelength"), log2(3))
  expect_equal(attr(r, "codelength"), net_map_equation(g, r))
  expect_identical(net_infomap(g, seed = 4), net_infomap(g, seed = 4))
  b <- comm_barbell()
  one <- net_map_equation(b, rep(1, 6))
  expect_lte(attr(net_infomap(b, seed = 1), "codelength"), one + 1e-12)
  # Directed networks use the directions and the node weights.
  d <- net_frame(data.frame(source = c(1, 2, 3, 3, 4, 5, 6, 4), target = c(2, 3, 1, 4, 5, 6, 4, 1)))
  r <- net_infomap(d, seed = 1)
  expect_equal(attr(r, "codelength"), net_map_equation(d, r))
  r2 <- net_infomap(d, node_weights = c(5, 1, 1, 1, 1, 1), seed = 1)
  expect_equal(attr(r2, "codelength"), net_map_equation(d, r2, node_weights = c(5, 1, 1, 1, 1, 1)))
  expect_error(net_infomap(d, node_weights = 1:3), "`node_weights`")
  expect_error(net_infomap(d, trials = 0), "`trials`")
})

# ---- Leiden -----------------------------------------------------------------------

test_that("Leiden reports the quality of its own partition", {
  g <- comm_barbell()
  r <- net_leiden(g, seed = 1)
  expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2))
  expect_equal(attr(r, "quality"), net_modularity(g, r))
  expect_equal(attr(r, "quality"), 5 / 14)
  expect_gte(attr(r, "iterations"), 1L)
  expect_identical(net_leiden(g, seed = 2), net_leiden(g, seed = 2))
  # Modularity resolution matches net_modularity().
  r2 <- net_leiden(g, resolution = 0.5, seed = 1)
  expect_equal(attr(r2, "quality"), net_modularity(g, r2, resolution = 0.5))
  # The Constant Potts Model with a high resolution splits into singletons
  # or edges; with a low one it joins the connected network.
  expect_equal(length(unique(net_leiden(g, "cpm", resolution = 0.01, seed = 1))), 1L)
  expect_equal(length(unique(net_leiden(g, "cpm", resolution = 5, seed = 1))), 6L)
  # A fixed number of iterations, and a starting partition.
  expect_identical(attr(net_leiden(g, n_iterations = 2, seed = 1), "iterations"), 2L)
  expect_equal(comm_groups(net_leiden(g, initial = c(1, 1, 1, 2, 2, 2), n_iterations = 1,
                                      seed = 1)), c(1, 1, 1, 2, 2, 2))
})

test_that("Leiden groups are connected and weights count", {
  g <- net_example("karate")
  r <- net_leiden(g, seed = 3)
  for (k in unique(r)) {
    members <- which(r == k)
    sub <- net_subnetwork(g, g$nodes$id[members])
    expect_identical(max(net_components(sub)$membership), 1L)
  }
  e <- data.frame(source = c(1, 1, 2, 3, 4, 4, 5), target = c(2, 3, 3, 4, 5, 6, 6),
                  w = c(1, 1, 1, 50, 1, 1, 1))
  gw <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  rw <- net_leiden(gw, seed = 1)
  expect_equal(rw[[3]], rw[[4]])
  expect_equal(attr(rw, "quality"), net_modularity(gw, rw))
})

test_that("Leiden checks its arguments", {
  g <- comm_barbell()
  expect_error(net_leiden(net_frame(data.frame(source = 1, target = 2))), "undirected")
  expect_error(net_leiden(g, objective = "potts"), "should be one of")
  expect_error(net_leiden(g, resolution = -1), "`resolution`")
  expect_error(net_leiden(g, beta = 0), "`beta`")
  expect_error(net_leiden(g, n_iterations = 1.5), "`n_iterations`")
  expect_error(net_leiden(g, seed = -1), "whole number")
  expect_error(net_leiden(g, node_weights = c(1, 2)), "`node_weights`")
  expect_error(net_leiden(g, initial = 1:2), "one value per node")
})

# ---- Spinglass ----------------------------------------------------------------------

test_that("spinglass finds the two triangles and reports its modularity", {
  g <- comm_barbell()
  for (rule in c("config", "simple")) {
    r <- net_spinglass(g, update_rule = rule, seed = 1)
    expect_equal(attr(r, "modularity"), net_modularity(g, r))
    expect_true(is.finite(attr(r, "energy")))
  }
  r <- net_spinglass(g, seed = 1)
  expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2))
  expect_identical(net_spinglass(g, seed = 7), net_spinglass(g, seed = 7))
  expect_equal(comm_groups(net_spinglass(g, spins = 1, seed = 1)), rep(1, 6))
  expect_length(net_spinglass(g, parupdate = TRUE, seed = 1), 6L)
  # Directions are ignored.
  d <- net_frame(data.frame(source = c("a", "a", "b", "c", "d", "d", "e"),
                            target = c("b", "c", "c", "d", "e", "f", "f")))
  expect_equal(comm_groups(net_spinglass(d, seed = 1)), c(1, 1, 1, 2, 2, 2))
})

test_that("spinglass with negative weights", {
  e <- data.frame(source = c(1, 1, 2, 3, 4, 4, 5), target = c(2, 3, 3, 4, 5, 6, 6),
                  w = c(1, 1, 1, -1, 1, 1, 1))
  g <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  r <- net_spinglass(g, implementation = "neg", seed = 1)
  expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2))
  expect_true(is.nan(attr(r, "modularity")))
  expect_error(net_spinglass(g, implementation = "orig", seed = 1), "net_spinglass")
})

test_that("spinglass checks its arguments and needs a connected network", {
  g <- comm_barbell()
  expect_error(net_spinglass(pieces(), seed = 1), "net_spinglass")
  expect_error(net_spinglass(g, spins = 0), "`spins`")
  expect_error(net_spinglass(g, start_temp = 0.001, stop_temp = 1), "`stop_temp`")
  expect_error(net_spinglass(g, cool_fact = 1), "`cool_fact`")
  expect_error(net_spinglass(g, gamma = -1), "`gamma`")
  expect_error(net_spinglass(g, parupdate = NA), "`parupdate`")
  expect_error(net_spinglass(g, update_rule = "both"), "should be one of")
})

# ---- Exact maximum modularity --------------------------------------------------------

test_that("the exact optimum is found and beats the heuristics", {
  g <- comm_barbell()
  r <- net_optimal(g)
  expect_equal(comm_groups(r), c(1, 1, 1, 2, 2, 2))
  expect_equal(attr(r, "modularity"), 5 / 14)
  expect_equal(attr(r, "modularity"), net_modularity(g, r))
  k <- net_example("karate")
  best <- net_optimal(k)
  expect_equal(attr(best, "modularity"), 0.4197896, tolerance = 1e-6)
  expect_gte(attr(best, "modularity") + 1e-12, attr(net_louvain(k, seed = 1), "modularity"))
  # Resolution, and a strong bridge that joins its endpoints.
  expect_equal(attr(net_optimal(g, resolution = 0.5), "modularity"),
               net_modularity(g, net_optimal(g, resolution = 0.5), resolution = 0.5))
  e <- data.frame(source = c(1, 1, 2, 3, 4, 4, 5), target = c(2, 3, 3, 4, 5, 6, 6),
                  w = c(1, 1, 1, 50, 1, 1, 1))
  gw <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  rw <- net_optimal(gw)
  expect_equal(rw[[3]], rw[[4]])
})

test_that("the exact search refuses what it cannot solve", {
  path <- net_frame(data.frame(source = 1:69, target = 2:70), directed = FALSE)
  expect_error(net_optimal(path), "net_optimal")
  d <- net_frame(data.frame(source = c(1, 2, 3, 3), target = c(2, 3, 1, 4)))
  r <- net_optimal(d)
  expect_equal(attr(r, "modularity"), net_modularity(d, r))
  empty <- net_frame(data.frame(source = integer(), target = integer()),
                     nodes = data.frame(id = 1:3), directed = FALSE)
  expect_identical(comm_groups(net_optimal(empty)), 1:3)
  expect_error(net_optimal(comm_barbell(), resolution = -1), "`resolution`")
})

# ---- Fluid communities -----------------------------------------------------------------

test_that("fluid communities return exactly k groups", {
  g <- net_example("karate")
  for (k in c(1, 2, 4, 34)) {
    r <- net_fluid(g, k, seed = 1)
    expect_identical(max(r), as.integer(k))
    expect_length(unique(r), k)
  }
  expect_identical(net_fluid(g, 3, seed = 5), net_fluid(g, 3, seed = 5))
  expect_error(net_fluid(g, 0), "`k`")
  expect_error(net_fluid(g, 35), "net_fluid")
  expect_error(net_fluid(g), "\"k\"")
  expect_error(net_fluid(comm_triangles(), 2), "net_fluid")
  expect_error(net_fluid(net_frame(data.frame(source = 1, target = 2)), 1), "undirected")
})

# ---- Voronoi communities -------------------------------------------------------------------

test_that("Voronoi communities report generators, modularity, and radius", {
  g <- comm_barbell()
  r <- net_voronoi(g, threads = 1)
  expect_equal(attr(r, "modularity"), net_modularity(g, r))
  expect_true(all(attr(r, "generators") %in% 1:6))
  expect_length(attr(r, "generators"), max(r))
  expect_true(attr(r, "radius") >= 0)
  # A radius of 0 makes every node a generator; an infinite one leaves a
  # single generator in a connected network.
  expect_identical(comm_groups(net_voronoi(g, radius = 0, threads = 1)), 1:6)
  one <- net_voronoi(g, radius = Inf, threads = 1)
  expect_identical(comm_groups(one), rep(1L, 6L))
  expect_length(attr(one, "generators"), 1L)
  expect_equal(attr(net_voronoi(g, radius = 0.7, threads = 1), "radius"), 0.7)
})

test_that("Voronoi communities take lengths, weights, and directions", {
  g <- comm_barbell()
  g$edges$len <- c(1, 1, 1, 5, 1, 1, 1)
  by_name <- net_voronoi(g, lengths = "len", threads = 1)
  by_value <- net_voronoi(g, lengths = g$edges$len, threads = 1)
  expect_identical(by_name, by_value)
  expect_error(net_voronoi(g, lengths = "nope"), "names no column")
  expect_error(net_voronoi(g, lengths = 1:3), "one per edge")
  expect_error(net_voronoi(g, lengths = c(-1, rep(1, 6))), "non-negative")
  expect_error(net_voronoi(g, radius = -1), "`radius`")
  d <- net_example("krackhardt", relation = "advice")
  for (m in c("out", "in", "all")) {
    r <- net_voronoi(d, mode = m, threads = 1)
    expect_equal(attr(r, "modularity"), net_modularity(d, r))
  }
  expect_error(net_voronoi(d, mode = "both"), "should be one of")
})

# ---- Comparing partitions ---------------------------------------------------------------------

test_that("partition comparisons follow their definitions", {
  a <- c(1, 1, 2, 2)
  b <- c(1, 2, 1, 2)
  expect_equal(net_compare_partitions(a, b, "vi"), 2 * log(2))
  expect_equal(net_compare_partitions(a, b, "nmi"), 0)
  expect_equal(net_compare_partitions(a, b, "split_join"), 4)
  expect_equal(net_compare_partitions(a, b, "rand"), 1 / 3)
  expect_equal(net_compare_partitions(a, b, "adjusted_rand"), -0.5)
  # Identical partitions, whatever the labels or their type.
  same <- c("x", "x", "y", "y")
  expect_equal(net_compare_partitions(a, same, "vi"), 0)
  expect_equal(net_compare_partitions(a, same, "nmi"), 1)
  expect_equal(net_compare_partitions(a, same, "split_join"), 0)
  expect_equal(net_compare_partitions(a, same, "rand"), 1)
  expect_equal(net_compare_partitions(a, factor(same), "adjusted_rand"), 1)
  expect_equal(net_compare_partitions(a, same, "ami"), 1)
  expect_equal(net_compare_partitions(a, 7 * a), 0)
  expect_true(is.nan(net_compare_partitions(integer(), integer())))
  expect_error(net_compare_partitions(a, 1:3), "same length")
  expect_error(net_compare_partitions(c(a[1:3], NA), b), "missing")
  expect_error(net_compare_partitions(a, b, "kappa"), "should be one of")
  expect_error(net_compare_partitions(list(1), b), "vector of group labels")
})

test_that("partition comparisons agree with igraph", {
  skip_if_not_installed("igraph")
  set.seed(11)
  a <- sample(4, 40, TRUE)
  b <- sample(3, 40, TRUE)
  b[1:15] <- a[1:15]
  expect_equal(net_compare_partitions(a, b, "vi"), igraph::compare(a, b, method = "vi"))
  expect_equal(net_compare_partitions(a, b, "nmi"), igraph::compare(a, b, method = "nmi"))
  expect_equal(net_compare_partitions(a, b, "split_join"),
               igraph::compare(a, b, method = "split.join"))
  expect_equal(net_compare_partitions(a, b, "rand"), igraph::compare(a, b, method = "rand"))
  expect_equal(net_compare_partitions(a, b, "adjusted_rand"),
               igraph::compare(a, b, method = "adjusted.rand"))
})

# ---- Dendrograms -----------------------------------------------------------------------------------

test_that("dendrograms are cut by groups, steps, and modularity", {
  m <- rbind(c(1, 2), c(3, 4), c(5, 6))
  expect_equal(comm_groups(net_cut_dendrogram(m, k = 2)), c(1, 1, 2, 2))
  expect_equal(comm_groups(net_cut_dendrogram(m, k = 1)), c(1, 1, 1, 1))
  expect_equal(comm_groups(net_cut_dendrogram(m, k = 4)), 1:4)
  expect_equal(comm_groups(net_cut_dendrogram(m, steps = 1)), c(1, 1, 2, 3))
  expect_equal(comm_groups(net_cut_dendrogram(m, steps = 0)), 1:4)
  expect_equal(comm_groups(net_cut_dendrogram(m, steps = 3)), c(1, 1, 1, 1))
  # Cut where modularity is highest.
  g <- net_frame(data.frame(source = c(1, 3), target = c(2, 4)), directed = FALSE)
  best <- net_cut_dendrogram(m, x = g)
  expect_equal(comm_groups(best), c(1, 1, 2, 2))
  expect_equal(attr(best, "modularity"), 0.5)
  expect_equal(attr(best, "modularity_path"), c(-0.25, 0.125, 0.5, 0))
  expect_error(net_cut_dendrogram(m), "Give `k`")
  expect_error(net_cut_dendrogram(m, k = 2, steps = 1), "at most one")
  expect_error(net_cut_dendrogram(m, k = 0), "`k`")
  expect_error(net_cut_dendrogram(m, k = 5), "net_cut_dendrogram")
  expect_error(net_cut_dendrogram(m, steps = 4), "net_cut_dendrogram")
  expect_error(net_cut_dendrogram(rbind(c(1, 9), c(3, 4)), k = 2), "net_cut_dendrogram")
  expect_error(net_cut_dendrogram(c(1, 2), k = 1), "two-column")
  expect_error(net_cut_dendrogram(rbind(c(1, 0)), k = 1), "two-column")
  other <- net_frame(data.frame(source = 1, target = 2), nodes = data.frame(id = 1:3),
                     directed = FALSE)
  expect_error(net_cut_dendrogram(structure(m, leaves = 4L), x = other), "describes 4 leaves")
})

test_that("dendrograms convert to hclust objects", {
  m <- rbind(c(1, 2), c(3, 4), c(5, 6))
  hc <- net_to_hclust(m)
  expect_s3_class(hc, "hclust")
  expect_equal(hc$merge, rbind(c(-1L, -2L), c(-3L, -4L), c(1L, 2L)))
  expect_equal(hc$height, c(1, 2, 3))
  expect_equal(hc$order, 1:4)
  expect_identical(hc$labels, as.character(1:4))
  for (k in 1:4) {
    expect_equal(net_compare_partitions(stats::cutree(hc, k = k), net_cut_dendrogram(m, k = k)), 0)
  }
  # A singleton is listed before the cluster it joins.
  expect_equal(net_to_hclust(rbind(c(1, 2), c(4, 3)))$merge, rbind(c(-1L, -2L), c(-3L, 1L)))
  # An incomplete dendrogram is completed above its last merge.
  part <- structure(rbind(c(1, 2)), leaves = 3L)
  expect_equal(net_to_hclust(part)$merge, rbind(c(-1L, -2L), c(-3L, 1L)))
  g <- comm_triangles()
  w <- net_walktrap(g, threads = 1)
  hc <- net_to_hclust(attr(w, "merges"), g)
  expect_identical(hc$labels, as.character(1:6))
  expect_identical(nrow(hc$merge), 5L)
  expect_equal(net_compare_partitions(stats::cutree(hc, k = 2), w), 0)
  expect_error(net_to_hclust(rbind(c(1, 2), c(2, 3))), "not a dendrogram")
  expect_error(net_to_hclust(rbind(c(1, 5), c(3, 4))), "not a dendrogram")
  expect_error(net_to_hclust(rbind(c(1, 1))), "not a dendrogram")
})

test_that("communities can be drawn with net_plot_partition()", {
  g <- net_example("karate")
  on_null_device({
    expect_silent(net_plot_partition(g, net_walktrap(g, threads = 1), seed = 1))
    expect_silent(net_plot_partition(g, net_leiden(g, seed = 1), seed = 1))
    expect_silent(plot(net_to_hclust(attr(net_walktrap(g, threads = 1), "merges"), g)))
  })
})
