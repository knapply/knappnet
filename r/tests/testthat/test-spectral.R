# PageRank, Katz and Bonacich power centrality, and hubs and authorities.
# Hand-checked values come from direct dense solutions of each definition
# (validation/oracles/spectral_definition_oracle.py); then spot checks
# against igraph and sna. The full comparison lives in validation/.

# 1 -> 2, 1 -> 3, 2 -> 3, 3 -> 1, 4 -> 3, 3 -> 5; node 5 has no out-ties.
small <- function() {
  netframe(data.frame(from = c(1, 1, 2, 3, 4, 3), to = c(2, 3, 3, 1, 3, 5)),
           nodes = data.frame(id = 1:5))
}

# A star around node 1 with an extra edge 2 - 5.
kite_tail <- function() {
  netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 5)),
           nodes = data.frame(id = 1:5), directed = FALSE)
}

test_that("PageRank solves its definition", {
  g <- small()
  pr <- pagerank(g)
  expect_equal(sum(pr), 1)
  expect_equal(pr, c(0.2142011097, 0.1574496602, 0.3477339318, 0.0664141886,
                     0.2142011097), tolerance = 1e-9)
  # Node 5 is dangling; it jumps by the personalization vector by default.
  p <- c(1, 0, 0, 0, 0)
  expect_equal(pagerank(g, personalized = p),
               c(0.3928645968, 0.1669674536, 0.3088897892, 0, 0.1312781604),
               tolerance = 1e-9)
  expect_equal(pagerank(g, personalized = p, dangling = "uniform"),
               c(0.3166510816, 0.1629073936, 0.3254597594, 0.0283306839, 0.1666510816),
               tolerance = 1e-9)
  g$nodes$p <- p
  expect_equal(pagerank(g, personalized = "p"), pagerank(g, personalized = p))
  # A directed 3-cycle is uniform.
  cyc <- netframe(data.frame(from = 1:3, to = c(2, 3, 1)))
  expect_equal(pagerank(cyc), rep(1 / 3, 3))
  expect_error(pagerank(g, damping = 2), "damping")
  expect_error(pagerank(g, personalized = c(1, 2)), "one value per node")
})

test_that("PageRank conventions differ only on undirected self-loops", {
  # A loop at node 1, then 1 - 2 - 3.
  g <- netframe(data.frame(from = c(1, 1, 2), to = c(1, 2, 3)), directed = FALSE)
  expect_equal(pagerank(g), c(0.4556354916, 0.3469224620, 0.1974420464),
               tolerance = 1e-9)
  expect_equal(pagerank(g, convention = "networkx"),
               c(0.3817177298, 0.3987945756, 0.2194876946), tolerance = 1e-9)
  h <- kite_tail()
  expect_equal(pagerank(h), pagerank(h, convention = "networkx"))
})

test_that("Katz centrality counts attenuated walks", {
  g <- small()
  expect_equal(katz_centrality(g, alpha = 0.2),
               c(1.3445378151, 1.2689075630, 1.7226890756, 1, 1.3445378151),
               tolerance = 1e-9)
  expect_equal(katz_centrality(g, alpha = 0.2, convention = "networkx"),
               c(0.4434413349, 0.4184977598, 0.5681592104, 0.3298094929, 0.4434413349),
               tolerance = 1e-9)
  # lambda_max is about 1.3247, so alpha must be below about 0.755.
  expect_error(katz_centrality(g, alpha = 0.8), "lambda_max")
  # Without cycles every alpha is allowed (the walk series is finite).
  dag <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(katz_centrality(dag, alpha = 10), c(1, 11, 111))
  expect_equal(katz_centrality(dag, alpha = 1, beta = c(1, 2, 3)), c(1, 3, 6))
})

test_that("power centrality follows Bonacich's definition", {
  g <- kite_tail()
  pos <- power_centrality(g, exponent = 0.5)
  expect_equal(sum(pos^2), 5)
  expect_equal(pos, c(1.4696735043, 1.1238679739, 0.7780624435, 0.7780624435,
                      0.6051596782), tolerance = 1e-9)
  expect_equal(power_centrality(g, exponent = -0.5),
               c(1.8605210188, -0.6201736729, -0.6201736729, -0.6201736729, 0.6201736729),
               tolerance = 1e-9)
  expect_equal(power_centrality(g, exponent = 0.5, rescale = TRUE),
               c(0.3090909091, 0.2363636364, 0.1636363636, 0.1636363636, 0.1272727273),
               tolerance = 1e-9)
  # Directed: walks along out-ties, so the sink scores 0.
  expect_equal(power_centrality(small(), exponent = 0.5),
               c(1.4766784404, 0.8306316227, 1.1998012328, 0.8306316227, 0),
               tolerance = 1e-9)
  # exponent 1 on a directed cycle: 1 is an eigenvalue, I - A is singular.
  cyc <- netframe(data.frame(from = 1:3, to = c(2, 3, 1)))
  expect_error(power_centrality(cyc), "singular")
  empty <- netframe(data.frame(from = integer(0), to = integer(0)),
                    nodes = data.frame(id = 1:3))
  expect_true(all(is.nan(power_centrality(empty))))
})

test_that("hub and authority scores", {
  g <- small()
  expect_equal(hub_score(g), c(1, sqrt(0.5), 0, sqrt(0.5), 0), tolerance = 1e-9)
  expect_equal(authority_score(g), c(0, sqrt(2) - 1, 1, 0, 0), tolerance = 1e-9)
  expect_equal(sum(authority_score(g, convention = "networkx")), 1)
  cyc <- netframe(data.frame(from = 1:4, to = c(2, 3, 4, 1)))
  expect_warning(h <- hub_score(cyc), "not unique")
  expect_equal(h, rep(1, 4))
})

test_that("undirected eigenvector centrality is unchanged by the Lanczos solver", {
  # Closed forms: star with k leaves, and a path (bipartite).
  star <- netframe(data.frame(from = rep(1, 4), to = 2:5), directed = FALSE)
  expect_equal(eigenvector_centrality(star), c(1, 0.5, 0.5, 0.5, 0.5), tolerance = 1e-10)
  path <- netframe(data.frame(from = 1:3, to = 2:4), directed = FALSE)
  phi <- (1 + sqrt(5)) / 2
  expect_equal(eigenvector_centrality(path), c(1 / phi, 1, 1, 1 / phi), tolerance = 1e-10)
  set.seed(7)
  e <- data.frame(from = sample(200, 800, TRUE), to = sample(200, 800, TRUE))
  g <- netframe(e, nodes = data.frame(id = 1:200), directed = FALSE)
  expect_identical(eigenvector_centrality(g, threads = 1L),
                   eigenvector_centrality(g, threads = 2L))
})

test_that("spectral measures match igraph", {
  skip_if_not_installed("igraph")
  set.seed(3)
  e <- data.frame(from = sample(40, 160, TRUE), to = sample(40, 160, TRUE))
  e <- e[e$from != e$to, ]
  g <- netframe(e, nodes = data.frame(id = 1:40))
  ig <- igraph::graph_from_data_frame(e, vertices = data.frame(id = 1:40))
  expect_equal(pagerank(g), igraph::page_rank(ig)$vector, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(pagerank(g, directed = FALSE),
               igraph::page_rank(ig, directed = FALSE)$vector, tolerance = 1e-8,
               ignore_attr = TRUE)
  a <- as.matrix(igraph::as_adjacency_matrix(ig, sparse = FALSE))
  lambda <- max(Mod(eigen(a, only.values = TRUE)$values))
  expect_equal(katz_centrality(g, alpha = 0.5 / lambda),
               igraph::alpha_centrality(ig, alpha = 0.5 / lambda), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(power_centrality(g, exponent = 0.5 / lambda),
               igraph::power_centrality(ig, exponent = 0.5 / lambda), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(power_centrality(g, exponent = -0.3),
               igraph::power_centrality(ig, exponent = -0.3), tolerance = 1e-8,
               ignore_attr = TRUE)
  # igraph 2.1 replaced hub_score() and authority_score() with hits_scores().
  if (exists("hits_scores", envir = asNamespace("igraph"))) {
    hits <- igraph::hits_scores(ig)
    hub <- hits$hub
    auth <- hits$authority
  } else {
    hub <- igraph::hub_score(ig)$vector
    auth <- igraph::authority_score(ig)$vector
  }
  expect_equal(hub_score(g), hub, tolerance = 1e-7, ignore_attr = TRUE)
  expect_equal(authority_score(g), auth, tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("power centrality matches sna", {
  skip_if_not_installed("sna")
  g <- kite_tail()
  m <- matrix(0, 5, 5)
  m[cbind(c(1, 1, 1, 2), c(2, 3, 4, 5))] <- 1
  m <- m + t(m)
  expect_equal(power_centrality(g, exponent = 0.5, convention = "sna"),
               sna::bonpow(m, gmode = "graph", exponent = 0.5), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(power_centrality(g, exponent = -0.5, convention = "sna"),
               sna::bonpow(m, gmode = "graph", exponent = -0.5), tolerance = 1e-8,
               ignore_attr = TRUE)
})
