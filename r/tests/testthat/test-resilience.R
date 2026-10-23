# Key players, group centrality, and flow betweenness. Hand-checked values;
# the comparisons with NetworkX, sna, keyplayer, and brute-force evaluation
# of the definitions live in validation/ (manifest.d/resilience.yaml).

path5 <- function() {
  netframe(data.frame(from = c("a", "b", "c", "d"), to = c("b", "c", "d", "e")),
           directed = FALSE)
}

# Two triangles joined through node 4: 1-2-3 and 5-6-7.
barbell <- function() {
  netframe(data.frame(from = c(1, 2, 3, 3, 4, 5, 6, 7), to = c(2, 3, 1, 4, 5, 6, 7, 5)),
           nodes = data.frame(id = 1:7), directed = FALSE)
}

test_that("key player criteria match hand calculations", {
  g <- path5()
  # Removing c leaves {a, b} and {d, e}: 4 of 12 ordered pairs connected.
  expect_equal(key_player_score(g, "c"), 8 / 12)
  expect_equal(key_player_score(g, "c", "distance_fragmentation"), 8 / 12)
  # Removing a leaves a path of four: sum of 1/d over ordered pairs 26/3.
  expect_equal(key_player_score(g, "a", "distance_fragmentation"), 1 - 26 / 36)
  # From c: itself 1, two nodes at distance 1, two at distance 2.
  expect_equal(key_player_score(g, "c", "distance_reach"), 4 / 5)
  expect_equal(key_player_score(g, c("a", "e"), "m_reach", m = 1), 2)
  expect_equal(key_player_score(g, "a", "m_reach", m = 2), 2)
  expect_equal(key_player_score(g, character(), "fragmentation"), 0)
  # Removing all but one node leaves no pairs.
  expect_true(is.nan(key_player_score(g, c("a", "b", "c", "d"))))
})

test_that("key player criteria follow arc directions", {
  d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(key_player_score(d, integer()), 0.5)
  expect_equal(key_player_score(d, integer(), mode = "all"), 0)
  expect_equal(key_player_score(d, 1, "distance_reach"), (1 + 1 + 1 / 2) / 3)
  expect_equal(key_player_score(d, 3, "distance_reach"), 1 / 3)
  expect_equal(key_player_score(d, 3, "distance_reach", mode = "in"), (1 + 1 + 1 / 2) / 3)
})

test_that("weights are path lengths", {
  w <- netframe(data.frame(from = c("a", "b"), to = c("b", "c"), len = c(1, 3)),
                directed = FALSE, weight = "len")
  expect_equal(key_player_score(w, "a", "distance_reach"), (1 + 1 + 1 / 4) / 3)
  expect_equal(key_player_score(w, "a", "distance_reach", weighted = FALSE),
               (1 + 1 + 1 / 2) / 3)
  expect_equal(key_player_score(w, "a", "m_reach", m = 3), 1)
  expect_equal(key_player_score(w, "a", "distance_reach", weights_as = "strength"),
               (1 + 1 + 1 / (1 + 1 / 3)) / 3)
})

test_that("greedy and exact searches find the barbell's centre", {
  g <- barbell()
  for (method in c("greedy", "exact")) {
    kp <- key_players(g, 1, "fragmentation", method = method, seed = 3, threads = 1)
    expect_equal(kp$set, 4L)
    expect_equal(kp$index, 4L)
    expect_equal(kp$score, 18 / 30)
    expect_equal(kp$measure, "fragmentation")
  }
  # Removing 3 or 4 leaves DF = 0.6 either way; ties go to the earlier node.
  expect_equal(key_players(g, 1, "distance_fragmentation", method = "exact")$set, 3L)
  # 3 reaches 1, 2, and 4 in one step.
  expect_equal(key_players(g, 1, "distance_reach", method = "exact")$score, 31 / 42)
  kp <- key_players(g, 2, "m_reach", m = 1, method = "exact")
  expect_equal(kp$score, 5)
  expect_equal(kp$evaluations, choose(7, 2))
})

test_that("greedy search is reproducible", {
  set.seed(1)
  e <- data.frame(from = sample(1:30, 60, TRUE), to = sample(1:30, 60, TRUE))
  g <- netframe(e, nodes = data.frame(id = 1:30), directed = FALSE)
  a <- key_players(g, 3, "distance_fragmentation", restarts = 4, seed = 11, threads = 1)
  b <- key_players(g, 3, "distance_fragmentation", restarts = 4, seed = 11, threads = 2)
  expect_identical(a, b)
  ex <- key_players(g, 2, "distance_reach", method = "exact", threads = 2)
  gr <- key_players(g, 2, "distance_reach", restarts = 20, seed = 2, threads = 2)
  expect_lte(gr$score, ex$score + 1e-12)
  expect_equal(key_player_score(g, gr$set, "distance_reach"), gr$score)
  set.seed(5)
  r1 <- key_players(g, 3, restarts = 2, threads = 1)
  set.seed(5)
  r2 <- key_players(g, 3, restarts = 2, threads = 1)
  expect_identical(r1, r2)
})

test_that("key player arguments are checked", {
  g <- path5()
  expect_error(key_players(g, 0), "`k`")
  expect_error(key_players(g, 6), "`k`")
  expect_error(key_players(g, 1, restarts = 0), "`restarts`")
  expect_error(key_players(g, 1, m = -1), "`m`")
  expect_error(key_player_score(g, "z"), "not in the network")
  expect_error(key_player_score(g, c("a", "a")), "more than once")
  expect_error(key_player_score(g, NA), "missing")
  expect_error(key_player_score(g, "a", weighted = TRUE), "no weights")
  expect_error(key_players(g, 1, seed = -1), "`seed`")
})

test_that("group centrality matches hand calculations", {
  g <- path5()
  expect_equal(group_centrality(g, c("b", "d"), "degree"), 1)
  expect_equal(group_centrality(g, c("b", "d"), "degree", normalized = FALSE), 3)
  expect_equal(group_centrality(g, c("b", "d"), "closeness"), 1)
  expect_equal(group_centrality(g, c("b", "d"), "closeness", normalized = FALSE), 1 / 3)
  # Pairs {a, c}, {a, e}, {c, e}: every shortest path crosses the group.
  expect_equal(group_centrality(g, c("b", "d"), "betweenness"), 1)
  expect_equal(group_centrality(g, c("b", "d"), "betweenness", normalized = FALSE), 3)
  # A 4-cycle: one of the two shortest paths from 1 to 3 passes through 2.
  c4 <- netframe(data.frame(from = 1:4, to = c(2:4, 1)), directed = FALSE)
  expect_equal(group_centrality(c4, 2, "betweenness", normalized = FALSE), 0.5)
  # Too few non-members to form a pair.
  expect_true(is.nan(group_centrality(g, c("a", "b", "c", "d"), "betweenness")))
})

test_that("group centrality in directed networks", {
  d <- netframe(data.frame(from = c(1, 1, 2, 4), to = c(2, 3, 4, 1)),
                nodes = data.frame(id = 1:4))
  expect_equal(group_centrality(d, 1, "degree"), 2 / 3)
  expect_equal(group_centrality(d, 1, "degree", mode = "in"), 1 / 3)
  expect_equal(group_centrality(d, 1, "degree", mode = "all"), 1)
  # From 1: 2 and 3 at distance 1, 4 at distance 2.
  expect_equal(group_centrality(d, 1, "closeness"), 3 / 4)
  # To 1: 4 at 1, 2 at 2; 3 cannot reach 1 and is left out.
  expect_equal(group_centrality(d, 1, "closeness", mode = "in"), 3 / 3)
  expect_error(group_centrality(d, 1, "degree", weighted = TRUE), "binary")
  expect_error(group_centrality(d, integer(), "degree"), "node identifiers")
})

test_that("flow betweenness matches hand calculations", {
  p <- netframe(data.frame(from = 1:2, to = 2:3), directed = FALSE)
  expect_equal(flow_betweenness(p), c(0, 1, 0))
  expect_equal(flow_betweenness(p, "normalized"), c(0, 1, 0))
  # A 4-cycle: each node loses one of two units between the opposite pair
  # and between each of its neighbours and the far node.
  c4 <- netframe(data.frame(from = 1:4, to = c(2:4, 1)), directed = FALSE)
  expect_equal(flow_betweenness(c4, threads = 2), rep(3, 4))
  expect_equal(flow_betweenness(c4, "fractional"), rep(1.5, 4))
  # Total flow between the other three pairs is 2 + 2 + 2 = 6.
  expect_equal(flow_betweenness(c4, "normalized"), rep(0.5, 4))
  # Directed: 1 -> 2 -> 3, pairs ordered.
  d <- netframe(data.frame(from = 1:2, to = 2:3))
  expect_equal(flow_betweenness(d), c(0, 1, 0))
})

test_that("flow betweenness uses capacities when asked", {
  # Two routes from a to d: through b (capacity 2) and through c (1).
  e <- data.frame(from = c("a", "b", "a", "c"), to = c("b", "d", "c", "d"),
                  w = c(2, 2, 1, 1))
  s <- netframe(e, directed = FALSE, weight = "w", weight_type = "strength")
  raw <- flow_betweenness(s)
  unit <- flow_betweenness(s, capacity = FALSE)
  expect_equal(unit, rep(3, 4))
  expect_false(isTRUE(all.equal(raw, unit)))
  dist <- netframe(e, directed = FALSE, weight = "w")
  expect_equal(flow_betweenness(dist), unit)
  expect_equal(flow_betweenness(dist, capacity = TRUE), raw)
  expect_error(flow_betweenness(path5(), capacity = TRUE), "no weights")
  neg <- netframe(transform(e, w = -w), directed = FALSE, weight = "w")
  expect_error(flow_betweenness(neg, capacity = TRUE), "non-negative")
})

test_that("flow betweenness matches sna", {
  skip_if_not_installed("sna")
  g <- netframe(kite_edges(), nodes = data.frame(id = 1:10), directed = FALSE)
  m <- as.matrix(g)
  for (scale in c("raw", "normalized", "fractional")) {
    cmode <- switch(scale, raw = "rawflow", normalized = "normflow", fractional = "fracflow")
    expect_equal(flow_betweenness(g, scale, threads = 2),
                 unname(sna::flowbet(m, gmode = "graph", cmode = cmode, ignore.eval = TRUE)),
                 tolerance = 1e-10, info = scale)
  }
  # Directed, with capacities.
  set.seed(2)
  e <- data.frame(from = sample(8, 24, TRUE), to = sample(8, 24, TRUE), w = sample(1:3, 24, TRUE))
  e <- e[e$from != e$to & !duplicated(e[c("from", "to")]), ]
  d <- netframe(e, nodes = data.frame(id = 1:8), weight = "w", weight_type = "strength")
  md <- matrix(0, 8, 8)
  md[cbind(e$from, e$to)] <- e$w
  expect_equal(flow_betweenness(d, threads = 2),
               unname(sna::flowbet(md, gmode = "digraph", cmode = "rawflow")), tolerance = 1e-10)
  expect_equal(flow_betweenness(d, capacity = FALSE, threads = 2),
               unname(sna::flowbet(md, gmode = "digraph", cmode = "rawflow", ignore.eval = TRUE)),
               tolerance = 1e-10)
})

test_that("distance-weighted fragmentation matches keyplayer's fragmentation centrality", {
  skip_if_not_installed("keyplayer")
  g <- netframe(kite_edges(), nodes = data.frame(id = 1:10), directed = FALSE)
  m <- as.matrix(g)
  ours <- vapply(1:10, function(v) key_player_score(g, v, "distance_fragmentation"), numeric(1))
  theirs <- vapply(1:10, function(v) keyplayer::fragment(m, v, binary = TRUE), numeric(1))
  expect_equal(ours, unname(theirs), tolerance = 1e-10)
})

test_that("the example networks give the published key players", {
  ji <- example_network("jemaah_islamiyah")
  kp <- key_players(ji, 2, "fragmentation", seed = 1, threads = 2)
  expect_equal(sort(kp$set), c("Muklas", "Samudra"))
  expect_equal(kp$score, key_player_score(ji, kp$set))
  ex <- key_players(ji, 2, "fragmentation", method = "exact", threads = 2)
  expect_equal(ex$set, kp$set)
  expect_equal(ex$evaluations, choose(17, 2))
  expect_equal(key_players(ji, 1, "distance_reach", method = "exact", threads = 2)$set, "Samudra")
  command <- ji$nodes$name[ji$nodes$role == "command team"]
  expect_equal(group_centrality(ji, command, "degree"), 13 / 14)
  fb <- flow_betweenness(ji, threads = 2)
  expect_equal(ji$nodes$name[which.max(fb)], "Samudra")
  expect_equal(fb[ji$nodes$name == "Samudra"], 317)
})
