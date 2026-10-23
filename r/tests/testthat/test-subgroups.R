# Cohesive subgroups and communities: k-cores, cliques and their
# relaxations, bicomponents, lambda sets, modularity, the community
# detection methods, factions, and hierarchical clustering
# (R/subgroups.R). Hand-checked values first, then igraph where it has the
# routine. The full comparison lives in validation/.

barbell <- function() {
  # Two triangles joined by the edge c-d; nodes appear in order a to f.
  netframe(data.frame(from = c("a", "a", "b", "c", "d", "d", "e"),
                      to   = c("b", "c", "c", "d", "e", "f", "f")),
           directed = FALSE)
}

# Two triangles {1, 2, 3} and {4, 5, 6} joined by the bridge 3 - 4, with a
# pendant 7 on 6.
two_triangles <- function() {
  netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5, 6), to = c(2, 3, 3, 4, 5, 6, 6, 7)),
           directed = FALSE)
}

hexagon <- function() netframe(data.frame(from = 1:6, to = c(2:6, 1)), directed = FALSE)

test_that("coreness follows the definition and its conventions", {
  g <- netframe(data.frame(from = c(1, 2, 3, 3), to = c(2, 3, 1, 4)),
                nodes = data.frame(id = 1:5), directed = FALSE)
  expect_equal(coreness(g), c(2, 2, 2, 1, 0))
  # A loop adds 2 in an undirected network; a doubled edge is a 2-core.
  loop <- netframe(data.frame(from = c(1, 1), to = c(1, 2)), directed = FALSE)
  expect_equal(coreness(loop), c(2, 1))
  d <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 3, 1)))
  expect_equal(coreness(d, "all"), c(2, 2, 2))
  expect_equal(coreness(d, "out"), c(1, 1, 1))
  expect_error(coreness(d, "both"))
})

test_that("maximal cliques return node identifiers", {
  g <- netframe(data.frame(from = c("a", "b", "a", "c"), to = c("b", "c", "c", "d")),
                nodes = data.frame(id = c("a", "b", "c", "d", "e")), directed = FALSE)
  expect_equal(max_cliques(g), list(c("a", "b", "c"), c("c", "d"), "e"))
  expect_equal(max_cliques(g, min_size = 3), list(c("a", "b", "c")))
  expect_identical(clique_number(g), 3L)
  expect_error(max_cliques(g, min_size = -1), "non-negative")
  # Direction, loops, and repeated edges do not matter.
  d <- netframe(data.frame(from = c(1, 2, 3, 3, 1), to = c(2, 3, 1, 3, 2)))
  expect_equal(max_cliques(d), list(c(1, 2, 3)))
})

test_that("n-cliques, n-clans, and k-plexes on a hexagon", {
  g <- hexagon()
  expect_equal(n_cliques(g, n = 1, threads = 2), max_cliques(g))
  cl2 <- n_cliques(g, n = 2, threads = 2)
  expect_equal(cl2, list(c(1, 2, 3), c(1, 2, 6), c(1, 3, 5), c(1, 5, 6), c(2, 3, 4),
                         c(2, 4, 6), c(3, 4, 5), c(4, 5, 6)))
  expect_equal(n_clans(g, n = 2, threads = 2),
               list(c(1, 2, 3), c(1, 2, 6), c(1, 5, 6), c(2, 3, 4), c(3, 4, 5), c(4, 5, 6)))
  expect_equal(k_plexes(g, k = 2, threads = 2), n_clans(g, n = 2, threads = 2))
  expect_equal(k_plexes(g, k = 1, min_size = 2, threads = 2), max_cliques(g))
  # min_size only filters, and the whole hexagon is a 3-clique.
  expect_equal(n_cliques(g, n = 2, min_size = 4, threads = 2), list())
  expect_equal(n_cliques(g, n = 3, threads = 2), list(1:6))
  expect_error(n_cliques(g, n = 0), "at least 1")
  expect_error(k_plexes(g, k = 2, max_steps = 0), "max_steps")
})

test_that("bicomponents cover every edge once", {
  g <- two_triangles()
  b <- bicomponents(g)
  expect_equal(b, list(c(1, 2, 3), c(3, 4), c(4, 5, 6), c(6, 7)))
  expect_equal(sort(cut_points(g)), c(3, 4, 6))
  # Directions, parallel edges, and loops do not matter; isolates belong
  # to no block.
  d <- netframe(data.frame(from = c(1, 2, 2, 3, 3), to = c(2, 1, 3, 1, 3)),
                nodes = data.frame(id = 1:4))
  expect_equal(bicomponents(d), list(c(1, 2, 3)))
  expect_equal(bicomponents(netframe(data.frame(from = integer(), to = integer()),
                                     nodes = data.frame(id = 1:3))), list())
})

test_that("lambda sets nest by connectivity", {
  g <- two_triangles()
  ls <- lambda_sets(g)
  expect_s3_class(ls, "lambda_sets")
  m <- ls$connectivity
  expect_equal(dim(m), c(7, 7))
  expect_equal(m["1", "2"], 2)
  expect_equal(m["1", "5"], 1)
  expect_equal(m["6", "7"], 1)
  expect_equal(diag(m), rep(0, 7), ignore_attr = TRUE)
  expect_equal(lapply(ls$sets, `[[`, "members"), list(1:7, c(1, 2, 3), c(4, 5, 6)))
  expect_equal(vapply(ls$sets, `[[`, numeric(1), "lambda"), c(1, 2, 2))
  expect_output(print(ls), "lambda = 2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(ls))
})

test_that("modularity matches hand-computed values", {
  g <- barbell()
  # NetworkX documents 5/14 for this split of barbell_graph(3, 0).
  expect_equal(modularity(g, c(1, 1, 1, 2, 2, 2)), 5 / 14, tolerance = 1e-14)
  expect_equal(modularity(g, c("x", "x", "x", "y", "y", "y")), 5 / 14, tolerance = 1e-14)
  expect_equal(modularity(g, rep(1, 6), resolution = 0.25), 0.75, tolerance = 1e-14)
  expect_error(modularity(g, 1:5), "length")
  expect_error(modularity(g, c(1, 1, 1, 2, 2, NA)), "missing")
  expect_error(modularity(g, rep(1, 6), resolution = -1), "non-negative")
  empty <- netframe(data.frame(from = integer(), to = integer()),
                    nodes = data.frame(id = 1:3), directed = FALSE)
  expect_true(is.nan(modularity(empty, 1:3)))
})

test_that("Louvain is reproducible and reports its own modularity", {
  g <- barbell()
  cl <- cluster_louvain(g, seed = 7)
  expect_equal(as.vector(cl), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(attr(cl, "modularity"), modularity(g, cl))
  expect_identical(cluster_louvain(g, seed = 3), cluster_louvain(g, seed = 3))
  set.seed(1)
  a <- cluster_louvain(g)
  set.seed(1)
  expect_identical(cluster_louvain(g), a)
  expect_error(cluster_louvain(netframe(data.frame(from = 1, to = 2))), "undirected")
  expect_error(cluster_louvain(g, seed = -1), "whole number")
})

test_that("weights are used when present", {
  e <- data.frame(from = c(1, 1, 2, 3, 4, 4, 5), to = c(2, 3, 3, 4, 5, 6, 6),
                  w = c(1, 1, 1, 50, 1, 1, 1))
  g <- netframe(e, directed = FALSE, weight = "w")
  cl <- cluster_louvain(g, seed = 1)
  expect_equal(cl[[3]], cl[[4]])
  expect_equal(attr(cl, "modularity"), modularity(g, cl))
  expect_false(isTRUE(all.equal(modularity(g, c(1, 1, 1, 2, 2, 2)),
                                modularity(g, c(1, 1, 1, 2, 2, 2), weighted = FALSE))))
  expect_error(modularity(barbell(), rep(1, 6), weighted = TRUE), "no weights")
})

test_that("Girvan-Newman, leading eigenvector, and greedy modularity", {
  g <- two_triangles()
  gn <- cluster_girvan_newman(g, threads = 2)
  expect_equal(as.vector(gn), c(1, 1, 1, 2, 2, 2, 2))
  expect_equal(attr(gn, "modularity"), modularity(g, gn))
  expect_equal(attr(gn, "removed")[1, ], data.frame(from = 3, to = 4), ignore_attr = TRUE)
  expect_equal(length(attr(gn, "path")), 9)
  expect_equal(max(attr(gn, "path")), attr(gn, "modularity"))
  tree <- attr(gn, "tree")
  expect_s3_class(tree, "hclust")
  expect_equal(sort(tree$order), 1:7)
  expect_equal(unname(cutree(tree, 2)), as.vector(gn))

  le <- cluster_leading_eigen(g)
  expect_equal(as.vector(le), c(1, 1, 1, 2, 2, 2, 2))
  expect_equal(attr(le, "modularity"), modularity(g, le))
  fg <- cluster_fast_greedy(g)
  expect_equal(as.vector(fg), c(1, 1, 1, 2, 2, 2, 2))
  expect_equal(attr(fg, "modularity"), modularity(g, fg))
  expect_equal(attr(cluster_fast_greedy(g, resolution = 0.5), "modularity"),
               modularity(g, cluster_fast_greedy(g, resolution = 0.5), resolution = 0.5))

  d <- netframe(data.frame(from = 1:3, to = c(2, 3, 1)))
  expect_error(cluster_girvan_newman(d), "undirected")
  expect_error(cluster_leading_eigen(d), "undirected")
  expect_error(cluster_fast_greedy(d), "undirected")
  empty <- netframe(data.frame(from = integer(), to = integer()), nodes = data.frame(id = 1:3),
                    directed = FALSE)
  expect_equal(as.vector(cluster_fast_greedy(empty)), 1:3)
  expect_true(is.nan(attr(cluster_fast_greedy(empty), "modularity")))
})

test_that("factions find the obvious split and report their fit", {
  g <- two_triangles()
  f <- factions(g, k = 2, seed = 1, threads = 2)
  expect_equal(as.vector(f), c(1, 1, 1, 2, 2, 2, 2))
  # Missing ties within {4, 5, 6, 7}: 4-7, 5-7; the bridge 3-4 between.
  expect_equal(attr(f, "errors"), 3)
  expect_equal(attr(f, "correlation"), 0.7076731, tolerance = 1e-6)
  expect_equal(factions(g, k = 2, fit = "correlation", seed = 1, threads = 2), f)
  expect_equal(factions(g, k = 2, seed = 1, threads = 1), f)
  # k = 1 has no between-group ties and an undefined correlation.
  one <- factions(g, k = 1, seed = 1, threads = 2)
  expect_equal(attr(one, "errors"), choose(7, 2) - 8)
  expect_true(is.nan(attr(one, "correlation")))
  expect_error(factions(g, k = 1, fit = "correlation", seed = 1), "more than one faction")
  expect_error(factions(g, k = 8, seed = 1), "between 1")
  # Node 2 is an isolate: the error count pairs it with the pendant 3, the
  # correlation criterion leaves it alone.
  h <- netframe(data.frame(from = c(1, 1, 3, 4, 4), to = c(5, 6, 4, 5, 6)),
                nodes = data.frame(id = 1:6), directed = FALSE)
  expect_equal(as.vector(factions(h, k = 2, seed = 1, threads = 2)), c(1, 2, 2, 1, 1, 1))
  expect_equal(as.vector(factions(h, k = 2, fit = "correlation", seed = 1, threads = 2)),
               c(1, 2, 1, 1, 1, 1))
})

test_that("hierarchical clustering wraps hclust", {
  # Four points on a line at 0, 1, 3, and 7: every distance is distinct.
  d <- as.matrix(dist(c(a = 0, b = 1, c = 3, d = 7)))
  tree <- cluster_hierarchical(d)
  expect_s3_class(tree, "hclust")
  expect_equal(tree$labels, c("a", "b", "c", "d"))
  # Average linkage: {a, b} at 1, then c at 2.5, then d at 17/3.
  expect_equal(tree$height, c(1, 2.5, 17 / 3))
  expect_equal(cut_clusters(tree, k = 2), c(1, 1, 1, 2))
  expect_equal(cut_clusters(tree, height = 2), c(1, 1, 2, 3))
  expect_equal(cluster_hierarchical(10 - d, similarity = TRUE)$height, tree$height)
  expect_equal(cluster_hierarchical(d, method = "single")$height, c(1, 2, 4))
  expect_equal(cluster_hierarchical(d, method = "complete")$height, c(1, 3, 7))
  # Unreachable pairs join last.
  pieces <- netframe(data.frame(from = c(1, 2), to = c(2, 3)), nodes = data.frame(id = 1:4),
                     directed = FALSE)
  tree2 <- cluster_hierarchical(distances(pieces), method = "single")
  expect_equal(max(tree2$height), 3)
  expect_equal(cut_clusters(tree2, k = 2), c(1, 1, 1, 2))
  expect_error(cluster_hierarchical(matrix(1:6, 2)), "square")
  expect_error(cluster_hierarchical(d, similarity = TRUE, method = "ward"), NA)
  expect_error(cut_clusters(tree), "exactly one")
  expect_error(cut_clusters(tree, k = 9), "between 1")
})

test_that("subgroups match igraph", {
  skip_if_not_installed("igraph")
  e <- kite_edges()
  g <- netframe(e, directed = FALSE)
  ig <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = g$nodes)
  expect_equal(coreness(g), unname(as.numeric(igraph::coreness(ig))))
  ours <- lapply(max_cliques(g), function(c) sort(match(c, g$nodes$name)))
  theirs <- lapply(igraph::max_cliques(ig), function(c) sort(as.integer(c)))
  key <- function(l) sort(vapply(l, paste, character(1), collapse = " "))
  expect_equal(key(ours), key(theirs))
  expect_equal(clique_number(g), igraph::clique_num(ig))
  memb <- rep(1:3, length.out = nrow(g$nodes))
  expect_equal(modularity(g, memb), igraph::modularity(ig, memb), tolerance = 1e-12)
  expect_equal(modularity(g, memb, resolution = 0.5),
               igraph::modularity(ig, memb, resolution = 0.5), tolerance = 1e-12)
  cl <- cluster_louvain(g, seed = 1)
  expect_equal(attr(cl, "modularity"), igraph::modularity(ig, as.vector(cl)),
               tolerance = 1e-12)

  # Directed coreness modes and directed modularity.
  set.seed(3)
  d <- data.frame(from = sample(25, 90, TRUE), to = sample(25, 90, TRUE))
  nodes <- data.frame(id = 1:25)
  gd <- netframe(d, nodes)
  igd <- igraph::graph_from_data_frame(d, vertices = nodes)
  for (mode in c("all", "out", "in")) {
    expect_equal(coreness(gd, mode), unname(as.numeric(igraph::coreness(igd, mode = mode))))
  }
  memb <- rep(1:4, length.out = 25)
  expect_equal(modularity(gd, memb), igraph::modularity(igd, memb, directed = TRUE),
               tolerance = 1e-12)
})

test_that("results match igraph", {
  skip_if_not_installed("igraph")
  g <- example_network("karate")
  ig <- igraph::graph_from_data_frame(g$edges, directed = FALSE, vertices = g$nodes)
  theirs <- vapply(igraph::biconnected_components(ig)$components,
                   function(c) paste(sort(g$nodes$name[as.integer(c)]), collapse = "|"),
                   character(1))
  ours <- vapply(bicomponents(g), function(b) paste(sort(b), collapse = "|"), character(1))
  expect_setequal(ours, theirs)
  lam <- lambda_sets(g)$connectivity
  for (pair in list(c(1, 34), c(1, 12), c(3, 10))) {
    expect_equal(lam[pair[1], pair[2]],
                 igraph::edge_connectivity(ig, source = pair[1], target = pair[2]))
  }
  fg <- cluster_fast_greedy(g)
  expect_equal(attr(fg, "modularity"),
               igraph::modularity(igraph::cluster_fast_greedy(ig)), tolerance = 1e-10)
  le <- cluster_leading_eigen(g)
  expect_equal(attr(le, "modularity"),
               igraph::modularity(igraph::cluster_leading_eigen(ig)), tolerance = 1e-10)
  gn <- cluster_girvan_newman(g, threads = 2)
  eb <- igraph::cluster_edge_betweenness(ig)
  expect_equal(attr(gn, "modularity"), max(eb$modularity), tolerance = 1e-10)
})
