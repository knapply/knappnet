# Mixing matrices, assortativity, the E-I index, and within-group tie
# shares. The Jemaah Islamiyah network (Koschade 2006) reproduces the
# examples of the R package homophily (Knapp 2019), with the numbers its
# README shows; the same data are in validation/fixtures-extra/.

jemaah <- function() {
  names <- c("MUKLAS", "AMROZI", "IMRON", "SAMUDRA", "DULMATIN", "IDRIS", "MUBAROK",
             "AZAHARI", "GHONI", "ARNASAN", "RAUF", "OCTAVIA", "HIDAYAT", "JUNAEDI",
             "PATEK", "FERI", "SARIJO")
  roles <- c("command team", "operation assistant", "operation assistant", "command team",
             "bomb maker", "command team", "operation assistant", "bomb maker", "bomb maker",
             "suicide bomber", "Team Lima", "Team Lima", "Team Lima", "Team Lima",
             "bomb maker", "suicide bomber", "bomb maker")
  from <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 8, 8, 8, 8, 9, 9, 9, 10, 10,
            10, 10, 11, 11, 11, 12, 12, 13, 15, 15, 16)
  to <- c(2, 3, 4, 5, 6, 8, 9, 15, 17, 4, 6, 7, 4, 5, 6, 8, 9, 15, 16, 17, 5, 6, 7, 8, 9,
          10, 11, 12, 13, 14, 15, 17, 6, 8, 9, 15, 16, 17, 7, 8, 9, 15, 17, 9, 15, 16, 17,
          15, 16, 17, 11, 12, 13, 14, 12, 13, 14, 13, 14, 14, 16, 17, 17)
  netframe(data.frame(from = names[from], to = names[to]),
           nodes = data.frame(name = names, role = roles), directed = FALSE)
}

# Interaction strengths of Koschade's Table 2, in the edge order of jemaah().
# The table leaves SARIJO-DULMATIN blank; DULMATIN-SARIJO is 5, and the
# table is otherwise symmetric, so 5 is used.
jemaah_weighted <- function() {
  g <- jemaah()
  strength <- c(2, 2, 1, 1, 5, 1, 1, 1, 1, 2, 4, 5, 3, 5, 3, 5, 5, 5, 1, 5, 2, 5, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 1, 5, 2, 2, 2, 2, 2, 5, 2, 1, 2, 5, 1, 5,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 5, 1)
  netframe(cbind(g$edges, strength = strength), nodes = g$nodes, directed = FALSE,
           weight = "strength", weight_type = "strength")
}

roles <- c("command team", "operation assistant", "bomb maker", "suicide bomber", "Team Lima")

test_that("mixing matrix of the Jemaah Islamiyah network matches homophily", {
  g <- jemaah()
  expect_equal(nrow(g$edges), 63L)
  m <- mixing_matrix(g, "role")
  expected <- matrix(c(6, 16, 30, 2, 8,
                       0, 2, 10, 2, 0,
                       0, 0, 20, 10, 0,
                       0, 0, 0, 0, 8,
                       0, 0, 0, 0, 12), nrow = 5, byrow = TRUE,
                     dimnames = list(roles, roles))
  expect_identical(m, expected)
  expect_equal(sum(m), 2 * 63)
  by_node <- mixing_matrix(g, "name", "role")
  expect_identical(dim(by_node), c(17L, 5L))
  expect_identical(rowSums(by_node), setNames(degree(g), g$nodes$name))
})

test_that("E-I index of the Jemaah Islamiyah network matches homophily", {
  g <- jemaah()
  expect_equal(ei_index(g, "role"), 0.3650794, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(ei_index(g, "role", scope = "group", convention = "homophily"),
               setNames(c(0.8064516, 0.7142857, -0.3333333, 1, -1), roles),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_identical(names(ei_index(g, "role", scope = "group")), roles)
  node <- c(0.5555556, 0.5, 1, 0.7333333, 0.1111111, 0.6, 0.3333333, 0.1111111, 0.1111111,
            1, -0.2, -0.2, -0.2, -0.2, 0.1111111, 1, 0.1111111)
  expect_equal(ei_index(g, "role", scope = "node"), node, tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(ei_index(g, "role", scope = "node", convention = "homophily"), node,
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("group E-I counts each tie for both groups (Krackhardt and Stern)", {
  g <- jemaah()
  # Operation assistants: 1 tie among themselves, 14 with other roles; the
  # homophily convention counts only the 6 ties to roles listed later.
  ei <- ei_index(g, "role", scope = "group")
  expect_equal(unname(ei), c(25 / 31, 13 / 15, 3 / 7, 1, 1 / 7))
  expect_equal(attr(ei, "internal"), c(3, 1, 10, 0, 6))
  expect_equal(attr(ei, "external"), c(28, 14, 25, 11, 8))
  # Every tie is counted once per side, and the network total once.
  expect_equal(sum(attr(ei, "internal")) + sum(attr(ei, "external")) / 2, 63)
})

test_that("assortativity of the Jemaah Islamiyah network matches homophily", {
  g <- jemaah()
  expect_equal(assortativity_nominal(g, "role", convention = "homophily"), 0.09078704,
               tolerance = 1e-6)
  expect_equal(assortativity_nominal(g, "role"), 0.09078704, tolerance = 1e-6)
})

test_that("strength assortativity of the weighted Jemaah Islamiyah network", {
  # Values from NetworkX's degree_assortativity_coefficient(), with and
  # without weight = "strength".
  g <- jemaah_weighted()
  expect_equal(sum(g$edges$strength), 166)
  expect_equal(assortativity_degree(g, weighted = TRUE), 0.175836797, tolerance = 1e-8)
  expect_equal(assortativity_degree(g), -0.022650693, tolerance = 1e-8)
  # Weights do not enter the E-I index or the mixing matrix.
  expect_identical(mixing_matrix(g, "role"), mixing_matrix(jemaah(), "role"))
  expect_equal(ei_index(g, "role"), ei_index(jemaah(), "role"))
})

test_that("directed mixing matrices follow homophily's layout", {
  # a -> b, a -> c, b -> c, c -> a; groups x, x, y.
  g <- netframe(data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "a")),
                nodes = data.frame(name = c("a", "b", "c"), grp = c("x", "x", "y")))
  m <- mixing_matrix(g, "grp")
  expect_identical(m, matrix(c(1, 2, 1, 0), 2, byrow = TRUE,
                             dimnames = list(c("x", "y"), c("x", "y"))))
  expect_identical(mixing_matrix(g, "grp", direction = "in"), t(m))
  out <- mixing_matrix(g, "name", "grp", direction = "out")
  expect_identical(out, matrix(c(1, 1, 0, 1, 1, 0), 3, byrow = TRUE,
                               dimnames = list(c("a", "b", "c"), c("x", "y"))))
  inn <- mixing_matrix(g, "name", "grp", direction = "in")
  # Rows: the sender's group; columns: the receiving node.
  expect_identical(inn, matrix(c(0, 1, 2, 1, 0, 0), 2, byrow = TRUE,
                               dimnames = list(c("x", "y"), c("a", "b", "c"))))
  all <- mixing_matrix(g, "name", "grp")
  expect_identical(all, matrix(c(1, 2, 1, 1, 3, 0), 3, byrow = TRUE,
                               dimnames = list(c("a", "b", "c"), c("x", "y"))))
})

test_that("assortativity coefficients", {
  # Newman's degree assortativity of a path with four nodes is -1/2.
  path <- netframe(data.frame(from = 1:3, to = 2:4), directed = FALSE)
  expect_equal(assortativity_degree(path), -0.5)
  expect_equal(assortativity(path, c(1, 2, 3, 4)), 5 / 11)
  expect_equal(assortativity(path, c(1, 2, 3, 4), normalized = FALSE), 5 / 12)
  # Perfect sorting by category.
  two <- netframe(data.frame(from = c(1, 3), to = c(2, 4)), directed = FALSE)
  expect_equal(assortativity_nominal(two, c("a", "a", "b", "b")), 1)
  expect_equal(assortativity_nominal(two, c("a", "a", "b", "b"), normalized = FALSE), 0.5)
  expect_true(is.nan(assortativity_nominal(two, c("a", "a", "a", "a"))))
  expect_error(assortativity_nominal(two, c("a", NA, "b", "b")), "missing")
  # Undirected self-loops: two arcs (igraph, homophily) or one (NetworkX).
  loop <- netframe(data.frame(from = c(1, 1, 2), to = c(1, 2, 3)), directed = FALSE)
  expect_false(isTRUE(all.equal(assortativity_degree(loop),
                                assortativity_degree(loop, convention = "networkx"))))
})

test_that("within-group tie shares", {
  g <- netframe(data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "d", "d")),
                nodes = data.frame(name = c("a", "b", "c", "d"), grp = c("x", "x", "y", "y")),
                directed = FALSE)
  h <- homophily(g, "grp")
  expect_identical(h$group, c("x", "y"))
  expect_identical(h$size, c(2L, 2L))
  expect_equal(h$internal, c(2, 2))
  expect_equal(h$external, c(2, 2))
  expect_equal(h$observed, c(0.5, 0.5))
  expect_equal(h$expected, c(1 / 3, 1 / 3))
})

test_that("parallel edges are reported", {
  g <- netframe(data.frame(from = c(1, 1, 2), to = c(2, 2, 3)), directed = FALSE)
  expect_warning(ei_index(g, c("a", "a", "b")), "parallel edges")
})

test_that("mixing measures match igraph", {
  skip_if_not_installed("igraph")
  set.seed(11)
  e <- data.frame(from = sample(30, 90, TRUE), to = sample(30, 90, TRUE))
  g <- netframe(e, nodes = data.frame(id = 1:30, grp = rep(1:3, 10),
                                      val = (7 * (0:29)) %% 11))
  ig <- igraph::graph_from_data_frame(e, vertices = data.frame(id = 1:30))
  expect_equal(assortativity_degree(g), igraph::assortativity_degree(ig), tolerance = 1e-10)
  expect_equal(assortativity_degree(g, directed = FALSE),
               igraph::assortativity_degree(ig, directed = FALSE), tolerance = 1e-10)
  expect_equal(assortativity_nominal(g, "grp"),
               igraph::assortativity_nominal(ig, types = rep(1:3, 10)), tolerance = 1e-10)
  expect_equal(assortativity(g, "val"),
               igraph::assortativity(ig, values = (7 * (0:29)) %% 11), tolerance = 1e-10)
})

test_that("homophily's samplike examples are reproduced", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")
  env <- new.env()
  utils::data("sampson", package = "ergm", envir = env)
  s <- as_netframe(env$samplike)
  groups <- c("Turks", "Outcasts", "Loyal")
  expect_identical(mixing_matrix(s, "group"),
                   matrix(c(30, 1, 5, 7, 10, 1, 9, 2, 23), 3, byrow = TRUE,
                          dimnames = list(groups, groups)))
  monks <- c("John Bosco", "Gregory", "Basil", "Peter", "Bonaventure", "Berthold", "Mark",
             "Victor", "Ambrose", "Romauld", "Louis", "Winfrid", "Amand", "Hugh", "Boniface",
             "Albert", "Elias", "Simplicius")
  by_monk <- matrix(c(9, 3, 5, 11, 3, 1, 3, 5, 0, 0, 0, 9, 3, 2, 8, 1, 0, 5, 8, 2, 1, 4, 0, 7,
                      2, 0, 6, 1, 1, 6, 3, 0, 5, 10, 0, 1, 1, 4, 3, 8, 0, 3, 8, 0, 1, 6, 0, 2,
                      1, 5, 0, 3, 6, 0), 18, byrow = TRUE, dimnames = list(monks, groups))
  expect_equal(mixing_matrix(s, s$roles$id, "group"), by_monk)
  expect_equal(mixing_matrix(s, "cloisterville", "group"),
               matrix(c(34, 15, 24, 48, 16, 39), 2, byrow = TRUE,
                      dimnames = list(c("TRUE", "FALSE"), groups)))
  expect_equal(mixing_matrix(s, "cloisterville", "group", direction = "out"),
               matrix(c(15, 5, 10, 31, 8, 19), 2, byrow = TRUE,
                      dimnames = list(c("TRUE", "FALSE"), groups)))
  expect_equal(mixing_matrix(s, "cloisterville", "group", direction = "in"),
               matrix(c(19, 17, 10, 8, 14, 20), 3, byrow = TRUE,
                      dimnames = list(groups, c("TRUE", "FALSE"))))
  expect_equal(ei_index(s, "group"), -0.4318182, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(ei_index(s, "group", scope = "group", convention = "homophily"),
               c(-0.6666667, -0.8181818, -1), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(ei_index(s, "group", scope = "node"),
               c(-0.05882353, -0.46666667, -0.25, -1, -0.23076923, -0.66666667,
                 -0.45454545, -0.27272727, -0.5, -0.5, -0.25, -0.81818182, 0,
                 -0.45454545, -0.77777778, -0.5, -0.66666667, -0.33333333),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(assortativity_nominal(s, "group", convention = "homophily"), 0.5445606,
               tolerance = 1e-6)
  expect_equal(assortativity_degree(s), 0.05569702, tolerance = 1e-6)
})
