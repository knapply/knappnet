# Key players, group centrality, and flow betweenness. Hand-checked values;
# the comparisons with NetworkX, sna, keyplayer, and brute-force evaluation
# of the definitions live in validation/ (manifest.d/resilience.yaml).

path5 <- function() {
  net_frame(data.frame(source = c("a", "b", "c", "d"), target = c("b", "c", "d", "e")),
           directed = FALSE)
}

# Two triangles joined through node 4: 1-2-3 and 5-6-7.
barbell <- function() {
  net_frame(data.frame(source = c(1, 2, 3, 3, 4, 5, 6, 7), target = c(2, 3, 1, 4, 5, 6, 7, 5)),
           nodes = data.frame(id = 1:7), directed = FALSE)
}

test_that("key player criteria match hand calculations", {
  g <- path5()
  # Removing c leaves {a, b} and {d, e}: 4 of 12 ordered pairs connected.
  expect_equal(net_key_player_score(g, "c"), 8 / 12)
  expect_equal(net_key_player_score(g, "c", "distance_fragmentation"), 8 / 12)
  # Removing a leaves a path of four: sum of 1/d over ordered pairs 26/3.
  expect_equal(net_key_player_score(g, "a", "distance_fragmentation"), 1 - 26 / 36)
  # From c: itself 1, two nodes at distance 1, two at distance 2.
  expect_equal(net_key_player_score(g, "c", "distance_reach"), 4 / 5)
  expect_equal(net_key_player_score(g, c("a", "e"), "m_reach", m = 1), 2)
  expect_equal(net_key_player_score(g, "a", "m_reach", m = 2), 2)
  expect_equal(net_key_player_score(g, character(), "fragmentation"), 0)
  # Removing all but one node leaves no pairs.
  expect_true(is.nan(net_key_player_score(g, c("a", "b", "c", "d"))))
})

test_that("key player criteria follow arc directions", {
  d <- net_frame(data.frame(source = c(1, 2), target = c(2, 3)))
  expect_equal(net_key_player_score(d, integer()), 0.5)
  expect_equal(net_key_player_score(d, integer(), mode = "all"), 0)
  expect_equal(net_key_player_score(d, 1, "distance_reach"), (1 + 1 + 1 / 2) / 3)
  expect_equal(net_key_player_score(d, 3, "distance_reach"), 1 / 3)
  expect_equal(net_key_player_score(d, 3, "distance_reach", mode = "in"), (1 + 1 + 1 / 2) / 3)
})

test_that("weights are path lengths", {
  w <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c"), len = c(1, 3)),
                directed = FALSE, weight = "len", weight_type = "distance")
  expect_equal(net_key_player_score(w, "a", "distance_reach"), (1 + 1 + 1 / 4) / 3)
  expect_equal(net_key_player_score(w, "a", "distance_reach", weighted = FALSE),
               (1 + 1 + 1 / 2) / 3)
  expect_equal(net_key_player_score(w, "a", "m_reach", m = 3), 1)
  expect_equal(net_key_player_score(w, "a", "distance_reach", weights_as = "strength"),
               (1 + 1 + 1 / (1 + 1 / 3)) / 3)
})

test_that("greedy and exact searches find the barbell's centre", {
  g <- barbell()
  for (method in c("greedy", "exact")) {
    kp <- net_key_players(g, 1, "fragmentation", method = method, seed = 3, threads = 1)
    expect_equal(kp$set, 4L)
    expect_equal(kp$index, 4L)
    expect_equal(kp$score, 18 / 30)
    expect_equal(kp$measure, "fragmentation")
  }
  # Removing 3 or 4 leaves DF = 0.6 either way; ties go to the earlier node.
  expect_equal(net_key_players(g, 1, "distance_fragmentation", method = "exact")$set, 3L)
  # 3 reaches 1, 2, and 4 in one step.
  expect_equal(net_key_players(g, 1, "distance_reach", method = "exact")$score, 31 / 42)
  kp <- net_key_players(g, 2, "m_reach", m = 1, method = "exact")
  expect_equal(kp$score, 5)
  expect_equal(kp$evaluations, choose(7, 2))
})

test_that("greedy search is reproducible", {
  set.seed(1)
  e <- data.frame(source = sample(1:30, 60, TRUE), target = sample(1:30, 60, TRUE))
  g <- net_frame(e, nodes = data.frame(id = 1:30), directed = FALSE)
  a <- net_key_players(g, 3, "distance_fragmentation", restarts = 4, seed = 11, threads = 1)
  b <- net_key_players(g, 3, "distance_fragmentation", restarts = 4, seed = 11, threads = 2)
  expect_identical(a, b)
  ex <- net_key_players(g, 2, "distance_reach", method = "exact", threads = 2)
  gr <- net_key_players(g, 2, "distance_reach", restarts = 20, seed = 2, threads = 2)
  expect_lte(gr$score, ex$score + 1e-12)
  expect_equal(net_key_player_score(g, gr$set, "distance_reach"), gr$score)
  set.seed(5)
  r1 <- net_key_players(g, 3, restarts = 2, threads = 1)
  set.seed(5)
  r2 <- net_key_players(g, 3, restarts = 2, threads = 1)
  expect_identical(r1, r2)
})

test_that("key player arguments are checked", {
  g <- path5()
  expect_error(net_key_players(g, 0), "`k`")
  expect_error(net_key_players(g, 6), "`k`")
  expect_error(net_key_players(g, 1, restarts = 0), "`restarts`")
  expect_error(net_key_players(g, 1, m = -1), "`m`")
  expect_error(net_key_player_score(g, "z"), "not in the network")
  expect_error(net_key_player_score(g, c("a", "a")), "more than once")
  expect_error(net_key_player_score(g, NA), "missing")
  expect_error(net_key_player_score(g, "a", weighted = TRUE), "no weights")
  expect_error(net_key_players(g, 1, seed = -1), "`seed`")
})

test_that("group centrality matches hand calculations", {
  g <- path5()
  expect_equal(net_group_centrality(g, c("b", "d"), "degree"), 1)
  expect_equal(net_group_centrality(g, c("b", "d"), "degree", normalized = FALSE), 3)
  expect_equal(net_group_centrality(g, c("b", "d"), "closeness"), 1)
  expect_equal(net_group_centrality(g, c("b", "d"), "closeness", normalized = FALSE), 1 / 3)
  # Pairs {a, c}, {a, e}, {c, e}: every shortest path crosses the group.
  expect_equal(net_group_centrality(g, c("b", "d"), "betweenness"), 1)
  expect_equal(net_group_centrality(g, c("b", "d"), "betweenness", normalized = FALSE), 3)
  # A 4-cycle: one of the two shortest paths from 1 to 3 passes through 2.
  c4 <- net_frame(data.frame(source = 1:4, target = c(2:4, 1)), directed = FALSE)
  expect_equal(net_group_centrality(c4, 2, "betweenness", normalized = FALSE), 0.5)
  # Too few non-members to form a pair.
  expect_true(is.nan(net_group_centrality(g, c("a", "b", "c", "d"), "betweenness")))
})

test_that("group centrality in directed networks", {
  d <- net_frame(data.frame(source = c(1, 1, 2, 4), target = c(2, 3, 4, 1)),
                nodes = data.frame(id = 1:4))
  expect_equal(net_group_centrality(d, 1, "degree"), 2 / 3)
  expect_equal(net_group_centrality(d, 1, "degree", mode = "in"), 1 / 3)
  expect_equal(net_group_centrality(d, 1, "degree", mode = "all"), 1)
  # From 1: 2 and 3 at distance 1, 4 at distance 2.
  expect_equal(net_group_centrality(d, 1, "closeness"), 3 / 4)
  # To 1: 4 at 1, 2 at 2; 3 cannot reach 1 and is left out.
  expect_equal(net_group_centrality(d, 1, "closeness", mode = "in"), 3 / 3)
  expect_error(net_group_centrality(d, 1, "degree", weighted = TRUE), "binary")
  expect_error(net_group_centrality(d, integer(), "degree"), "at least one node")
})

test_that("flow betweenness matches hand calculations", {
  p <- net_frame(data.frame(source = 1:2, target = 2:3), directed = FALSE)
  expect_equal(net_flow_betweenness(p), c(0, 1, 0))
  expect_equal(net_flow_betweenness(p, "normalized"), c(0, 1, 0))
  # A 4-cycle: each node loses one of two units between the opposite pair
  # and between each of its neighbours and the far node.
  c4 <- net_frame(data.frame(source = 1:4, target = c(2:4, 1)), directed = FALSE)
  expect_equal(net_flow_betweenness(c4, threads = 2), rep(3, 4))
  expect_equal(net_flow_betweenness(c4, "fractional"), rep(1.5, 4))
  # Total flow between the other three pairs is 2 + 2 + 2 = 6.
  expect_equal(net_flow_betweenness(c4, "normalized"), rep(0.5, 4))
  # Directed: 1 -> 2 -> 3, pairs ordered.
  d <- net_frame(data.frame(source = 1:2, target = 2:3))
  expect_equal(net_flow_betweenness(d), c(0, 1, 0))
})

test_that("flow betweenness uses capacities when asked", {
  # Two routes from a to d: through b (capacity 2) and through c (1).
  e <- data.frame(source = c("a", "b", "a", "c"), target = c("b", "d", "c", "d"),
                  w = c(2, 2, 1, 1))
  s <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  raw <- net_flow_betweenness(s)
  unit <- net_flow_betweenness(s, capacity = FALSE)
  expect_equal(unit, rep(3, 4))
  expect_false(isTRUE(all.equal(raw, unit)))
  dist <- net_frame(e, directed = FALSE, weight = "w", weight_type = "distance")
  expect_equal(net_flow_betweenness(dist), unit)
  expect_equal(net_flow_betweenness(dist, capacity = TRUE), raw)
  expect_error(net_flow_betweenness(path5(), capacity = TRUE), "no weights")
  expect_error(net_frame(transform(e, w = -w), directed = FALSE, weight = "w",
                         weight_type = "distance"), "positive")
})

test_that("flow betweenness matches sna", {
  skip_if_not_installed("sna")
  g <- net_frame(kite_edges(), nodes = data.frame(id = 1:10), directed = FALSE)
  m <- as.matrix(g)
  for (scale in c("raw", "normalized", "fractional")) {
    cmode <- switch(scale, raw = "rawflow", normalized = "normflow", fractional = "fracflow")
    expect_equal(net_flow_betweenness(g, scale, threads = 2),
                 unname(sna::flowbet(m, gmode = "graph", cmode = cmode, ignore.eval = TRUE)),
                 tolerance = 1e-10, info = scale)
  }
  # Directed, with capacities.
  set.seed(2)
  e <- data.frame(source = sample(8, 24, TRUE), target = sample(8, 24, TRUE), w = sample(1:3, 24, TRUE))
  e <- e[e$source != e$target & !duplicated(e[c("source", "target")]), ]
  d <- net_frame(e, nodes = data.frame(id = 1:8), weight = "w", weight_type = "strength")
  md <- matrix(0, 8, 8)
  md[cbind(e$source, e$target)] <- e$w
  expect_equal(net_flow_betweenness(d, threads = 2),
               unname(sna::flowbet(md, gmode = "digraph", cmode = "rawflow")), tolerance = 1e-10)
  expect_equal(net_flow_betweenness(d, capacity = FALSE, threads = 2),
               unname(sna::flowbet(md, gmode = "digraph", cmode = "rawflow", ignore.eval = TRUE)),
               tolerance = 1e-10)
})

test_that("distance-weighted fragmentation matches keyplayer's fragmentation centrality", {
  skip_if_not_installed("keyplayer")
  g <- net_frame(kite_edges(), nodes = data.frame(id = 1:10), directed = FALSE)
  m <- as.matrix(g)
  ours <- vapply(1:10, function(v) net_key_player_score(g, v, "distance_fragmentation"), numeric(1))
  theirs <- vapply(1:10, function(v) keyplayer::fragment(m, v, binary = TRUE), numeric(1))
  expect_equal(ours, unname(theirs), tolerance = 1e-10)
})

test_that("the example networks give the published key players", {
  ji <- net_example("jemaah_islamiyah")
  kp <- net_key_players(ji, 2, "fragmentation", seed = 1, threads = 2)
  expect_equal(sort(kp$set), c("Muklas", "Samudra"))
  expect_equal(kp$score, net_key_player_score(ji, kp$set))
  ex <- net_key_players(ji, 2, "fragmentation", method = "exact", threads = 2)
  expect_equal(ex$set, kp$set)
  expect_equal(ex$evaluations, choose(17, 2))
  expect_equal(net_key_players(ji, 1, "distance_reach", method = "exact", threads = 2)$set, "Samudra")
  command <- ji$nodes$id[ji$nodes$role == "command team"]
  expect_equal(net_group_centrality(ji, command, "degree"), 13 / 14)
  fb <- net_flow_betweenness(ji, threads = 2)
  expect_equal(ji$nodes$id[which.max(fb)], "Samudra")
  expect_equal(fb[ji$nodes$id == "Samudra"], 317)
})
