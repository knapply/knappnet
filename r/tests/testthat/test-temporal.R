# Temporal networks: construction, slicing, aggregation, time-respecting
# paths, reachability, temporal closeness, rolling measures, and event
# statistics. Values are checked by hand; the full comparison against a
# brute-force search over walks lives in validation/.

# a -> b and b -> c at time 1 (they chain unless times must increase),
# b -> c again at 3, and a -> c at 6.
msgs <- function() {
  net_temporal(data.frame(source = c("a", "b", "b", "a"), target = c("b", "c", "c", "c"),
                              time = c(1, 1, 3, 6)))
}

test_that("temporal networks are net_frames with time roles", {
  g <- msgs()
  expect_s3_class(g, "net_frame")
  expect_true(net_is_temporal(g))
  expect_false(net_is_temporal(net_frame(data.frame(source = "a", target = "b"))))
  expect_false(net_is_temporal(data.frame(time = 1)))
  expect_equal(g$roles$start, "time")
  expect_equal(g$roles$time_class, "numeric")
  expect_equal(g$nodes$id, c("a", "b", "c"))
  expect_output(print(g), "time:      contacts (4 events, 1 to 6)", fixed = TRUE)
  expect_output(print(g), "clock:     numbers", fixed = TRUE)
  # Static measures see every contact as an edge.
  expect_equal(net_degree(g), c(2, 3, 3))
})

test_that("start, end, and duration columns are found and checked", {
  d <- data.frame(source = c("a", "b"), target = c("b", "c"), start = c(0, 5), stop = c(3, 5),
                  len = c(3, 0))
  g <- net_temporal(d, start = "start", end = "stop")
  expect_equal(g$roles$start, "start")
  expect_output(print(g), "edge spells, 0 to 5")
  expect_equal(nrow(net_at(g, 2)$edges), 1L)
  h <- net_temporal(d, start = "start", duration = "len")
  expect_equal(nrow(net_at(h, 2)$edges), 1L)
  expect_equal(net_aggregate(h)$edges$duration, c(3, 0))
  expect_error(net_temporal(d, start = "start", end = "stop", duration = "len"), "not both")
  bad <- data.frame(source = c("a", "b"), target = c("b", "c"), time = c(2, 1), end = c(1, 1))
  expect_error(net_temporal(bad, end = "end"), "must not end before they start.*row 1")
  expect_error(net_temporal(data.frame(source = "a", target = "b", time = NA_real_)),
               "must not be missing")
  expect_error(net_temporal(data.frame(source = "a", target = "b", when = 1)), "`start`")
  expect_error(net_temporal(data.frame(source = "a", target = "b", time = "noon")),
               "numbers, dates, or date-times")
  mixed <- data.frame(source = "a", target = "b", time = 1, end = as.Date("2024-01-01"))
  expect_error(net_temporal(mixed, end = "end"), "must be numbers")
  expect_error(net_temporal_reach(net_frame(data.frame(source = "a", target = "b"))),
               "temporal network")
})

test_that("a net_frame whose edges carry a time column can be made temporal", {
  x <- net_frame(data.frame(src = c("a", "b"), dst = c("b", "c"), start = c(1, 2),
                            bytes = c(10, 20)), source = "src", target = "dst", weight = "bytes",
                 weight_type = "distance")
  g <- net_temporal(x, start = "start")
  expect_true(net_is_temporal(g))
  expect_equal(g$roles$weight, "bytes")
  expect_equal(net_temporal_reach(g), c(2, 1, 0))
  expect_identical(net_temporal(g)$roles, g$roles)
  expect_error(net_temporal(x, nodes = data.frame(id = "a")), "must be NULL")
})

test_that("slicing follows half-open intervals", {
  calls <- data.frame(source = c("a", "a", "b", "c"), target = c("b", "c", "c", "d"),
                      time = c(1, 2, 2, 5), end = c(3, 2, 4, 6))
  g <- net_temporal(calls, end = "end", directed = FALSE)
  at2 <- net_at(g, 2)
  expect_true(net_is_temporal(at2))
  expect_equal(nrow(at2$edges), 3L)
  expect_equal(nrow(at2$nodes), 4L)
  expect_equal(nrow(net_at(g, 3)$edges), 1L)
  w <- net_window(g, 3, 5, nodes = "active")
  expect_equal(w$edges$source, "b")
  expect_equal(w$nodes$id, c("b", "c"))
  expect_equal(nrow(net_window(g, 0, 10)$edges), 4L)
  expect_error(net_window(g, 5, 5), "earlier")
  expect_error(net_at(g, as.Date("2024-01-01")), "single time")
})

test_that("node presence intervals select nodes", {
  nodes <- data.frame(id = c("a", "b", "c"), joined = c(0, 5, NA))
  g <- net_temporal(data.frame(source = c("a", "b", "b", "a"), target = c("c", "c", "c", "b"),
                               time = c(1, 6, 7, 8)),
                    nodes = nodes, node_start = "joined")
  expect_equal(net_at(g, 2, nodes = "present")$nodes$id, c("a", "c"))
  expect_equal(nrow(net_at(g, 2, nodes = "active")$nodes), 0L)
  expect_equal(net_at(g, 2)$nodes$id, c("a", "b", "c"))
  expect_equal(net_at(g, 6, nodes = "present")$nodes$id, c("a", "b", "c"))
  expect_equal(net_at(g, 6, nodes = "active")$nodes$id, c("b", "c"))
  # A contact before a node joins is an error that names the row.
  early <- data.frame(source = c("a", "b"), target = c("c", "c"), time = c(1, 1))
  expect_error(net_temporal(early, nodes, node_start = "joined"), "row 2")
})

test_that("snapshots and rolling measures", {
  logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2"),
                       host = c("h1", "h1", "h2", "h2", "h3"), time = c(0, 1, 5, 6, 12))
  g <- net_temporal(logins, source = "user", target = "host")
  s <- net_snapshots(g, width = 5)
  expect_equal(attr(s, "windows")$start, c(0, 5, 10))
  expect_equal(attr(s, "windows")$end, c(5, 10, 15))
  expect_equal(vapply(s, function(w) nrow(w$edges), numeric(1)), c(2, 2, 1))
  slide <- net_snapshots(g, width = 6, step = 3)
  expect_equal(vapply(slide, function(w) nrow(w$edges), numeric(1)), c(3, 2, 1, 1, 1))
  b <- net_snapshots(g, breaks = c(0, 6, 13))
  expect_equal(vapply(b, function(w) nrow(w$edges), numeric(1)), c(3, 2))
  expect_error(net_snapshots(g), "width")
  expect_error(net_snapshots(g, width = 5, breaks = c(0, 1)), "not both")
  expect_error(net_snapshots(g, breaks = c(3, 1)), "increasing")
  r <- net_rolling(g, net_degree, width = 5)
  expect_equal(names(r), c("start", "end", "node", "value"))
  expect_equal(nrow(r), 18L)
  expect_equal(r$value[r$start == 0], c(1, 2, 1, 0, 0, 0))
  expect_equal(r$node[1:6], c("u1", "h1", "u2", "h2", "u3", "h3"))
  dens <- net_rolling(g, function(w) nrow(w$edges), width = 5)
  expect_equal(dens$value, c(2, 2, 1))
  expect_true(all(is.na(dens$node)))
  expect_error(net_rolling(g, function(w) c(1, 2), width = 5), "one value per node")
})

test_that("aggregation counts and clips contacts", {
  calls <- data.frame(source = c("a", "b", "a", "c"), target = c("b", "a", "c", "a"),
                      time = c(0, 10, 20, 25), end = c(5, 12, 20, 40))
  g <- net_temporal(calls, end = "end", directed = FALSE)
  a <- net_aggregate(g)
  expect_false(net_is_temporal(a))
  expect_equal(a$roles$weight, "count")
  expect_equal(a$roles$weight_type, "strength")
  expect_equal(a$edges$source, c("a", "a"))
  expect_equal(a$edges$target, c("b", "c"))
  expect_equal(a$edges$count, c(2, 2))
  expect_equal(a$edges$duration, c(7, 15))
  expect_equal(a$edges$first, c(0, 20))
  expect_equal(a$edges$last, c(12, 40))
  w <- net_aggregate(g, weight = "duration", from = 0, to = 30)
  expect_equal(w$roles$weight, "duration")
  expect_equal(w$edges$duration, c(7, 5))
  expect_error(net_aggregate(g, from = 0), "both")
  d <- net_aggregate(net_temporal(calls, end = "end"))
  expect_equal(nrow(d$edges), 4L)
})

test_that("time-respecting paths", {
  g <- msgs()
  ids <- c("a", "b", "c")
  ea <- matrix(c(1, 1, 1, Inf, 1, 1, Inf, Inf, 1), 3, byrow = TRUE,
               dimnames = list(ids, ids))
  expect_equal(net_temporal_distances(g), ea)
  # Strictly increasing times: b -> c at 1 no longer follows a -> b at 1.
  expect_equal(unname(net_temporal_distances(g, from = "a", strict = TRUE)[1, ]), c(1, 1, 3))
  ld <- matrix(c(6, 1, 6, -Inf, 6, 3, -Inf, -Inf, 6), 3, byrow = TRUE,
               dimnames = list(ids, ids))
  expect_equal(net_temporal_distances(g, "latest_departure"), ld)
  expect_equal(net_temporal_distances(g, "latest_departure", from = "b", to = c("c", "b")),
               matrix(c(3, 6), 1, dimnames = list("b", c("c", "b"))))
  fast <- net_temporal_distances(g, "fastest")
  expect_equal(unname(fast[1, ]), c(0, 0, 0))
  expect_equal(unname(fast[3, ]), c(Inf, Inf, 0))
  expect_equal(unname(net_temporal_distances(g, "shortest")["a", ]), c(0, 1, 1))
  # Waiting at most 1 at b: the walk a -> b (1), b -> c (3) is not allowed.
  expect_equal(net_temporal_distances(g, from = "a", strict = TRUE, max_wait = 1)["a", "c"], 6)
  # One time unit to traverse each contact.
  expect_equal(unname(net_temporal_distances(g, from = "a", traversal = 1)[1, ]), c(1, 2, 4))
  expect_equal(net_temporal_distances(g, from = "a", start = 2)["a", "c"], 6)
  expect_equal(net_temporal_distances(g, from = "a", end = 2)["a", "c"], 1)
  avg <- net_temporal_distances(g, "average", start = 0, end = 10)
  expect_equal(unname(avg[1, ]), c(0, 5, 2.5))
  expect_equal(avg["b", "c"], 3.4)
  expect_error(net_temporal_distances(g, from = "z"), "not in the network")
  expect_error(net_temporal_distances(g, start = 5, end = 1), "later")
  expect_error(net_temporal_distances(g, min_wait = 2, max_wait = 1), "smaller")
  expect_error(net_temporal_distances(g, traversal = "flight"), "duration")
  expect_error(net_temporal_distances(g, min_wait = -1), "non-negative")
  expect_equal(net_temporal_distances(g, "fastest", threads = 1),
               net_temporal_distances(g, "fastest", threads = 2))
})

test_that("contact durations as traversal times", {
  # A 2-unit call from a to b, then b calls c at 1 (during the first call)
  # and at 3 (after it).
  d <- data.frame(source = c("a", "b", "b"), target = c("b", "c", "c"), time = c(0, 1, 3),
                  end = c(2, 2, 4))
  g <- net_temporal(d, end = "end")
  expect_equal(net_temporal_distances(g, from = "a")["a", "c"], 1)
  expect_equal(net_temporal_distances(g, from = "a", traversal = "duration")["a", "c"], 4)
  expect_equal(net_temporal_distances(g, "fastest", traversal = "duration")["a", "c"], 4)
})

test_that("reachability and temporal closeness", {
  g <- net_temporal(data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"),
                                   time = c(1, 2, 1)))
  expect_equal(net_temporal_reach(g), c(2, 1, 1, 0))
  expect_equal(net_temporal_reach(g, mode = "in"), c(0, 1, 2, 1))
  expect_equal(net_reach_ratio(g), 1 / 3)
  expect_equal(net_temporal_reach(g, max_wait = 0.5), c(1, 1, 1, 0))
  one <- net_temporal(data.frame(source = "a", target = "a", time = 1))
  expect_true(is.nan(net_reach_ratio(one)))
  m <- msgs()
  # a reaches b on average in 5 and c in 2.5 over starting times in [0, 10].
  expect_equal(net_temporal_closeness(m, start = 0, end = 10), c(0.3, 5 / 34, 0))
  expect_equal(net_temporal_closeness(m, start = 0, end = 10, strict = TRUE)[1], 17 / 70)
  expect_true(all(is.nan(net_temporal_closeness(m, start = 1, end = 1))))
})

test_that("dates and date-times", {
  d <- data.frame(source = c("a", "b"), target = c("b", "c"),
                  day = as.Date(c("2024-01-01", "2024-01-03")))
  g <- net_temporal(d, start = "day")
  expect_equal(g$roles$time_class, "Date")
  expect_output(print(g), "2024-01-01 to 2024-01-03")
  expect_equal(net_temporal_distances(g, "fastest")["a", "c"], 2)
  expect_equal(net_temporal_distances(g, from = "a")["a", "c"], 19725)
  expect_equal(net_temporal_distances(g, "fastest", max_wait = 1)["a", "c"], Inf)
  expect_equal(net_temporal_distances(g, "fastest", max_wait = as.difftime(2, units = "days"))["a", "c"], 2)
  s <- net_snapshots(g, width = 1)
  expect_s3_class(attr(s, "windows")$start, "Date")
  expect_equal(length(s), 3L)
  expect_equal(nrow(net_at(g, as.Date("2024-01-03"))$edges), 1L)
  expect_error(net_at(g, 19725), "a date")
  expect_equal(net_inter_event(g)$gap, 2)

  p <- data.frame(src = c("a", "b"), dst = c("b", "c"),
                  t = as.POSIXct(c("2024-05-01 10:00:00", "2024-05-01 11:00:00"), tz = "UTC"))
  h <- net_temporal(p, start = "t", source = "src", target = "dst")
  expect_equal(h$roles$time_class, "POSIXct")
  expect_equal(net_temporal_distances(h, "fastest")["a", "c"], 3600)
  expect_equal(net_temporal_distances(h, "fastest",
                                  max_wait = as.difftime(30, units = "mins"))["a", "c"], Inf)
  w <- net_window(h, as.POSIXct("2024-05-01 10:30:00", tz = "UTC"),
                      as.POSIXct("2024-05-01 12:00:00", tz = "UTC"))
  expect_equal(nrow(w$edges), 1L)
  a <- net_aggregate(h)
  expect_s3_class(a$edges$first, "POSIXct")
  expect_equal(format(a$edges$first[2], tz = "UTC"), "2024-05-01 11:00:00")
})

test_that("inter-event times, burstiness, and memory", {
  ev <- data.frame(source = c("a", "a", "a", "a", "b", "b", "b"),
                   target = c("b", "b", "c", "c", "c", "c", "c"),
                   time = c(0, 1, 2, 10, 0, 5, 10))
  g <- net_temporal(ev, directed = FALSE)
  ie <- net_inter_event(g)
  expect_equal(ie$node, rep(c("a", "b", "c"), c(3, 4, 4)))
  expect_equal(ie$gap, c(1, 1, 8, 0, 1, 4, 5, 2, 3, 5, 0))
  b <- net_burstiness(g)
  expect_equal(b$node, c("a", "b", "c"))
  expect_equal(b$events, c(4, 5, 5))
  gaps <- c(1, 1, 8)
  s <- sqrt(mean((gaps - mean(gaps))^2))
  expect_equal(b$burstiness[1], (s - mean(gaps)) / (s + mean(gaps)))
  expect_true(is.nan(b$memory[1]))
  # The memory coefficient is the correlation of consecutive gaps.
  expect_equal(b$memory[2], cor(c(0, 1, 4), c(1, 4, 5)))
  p <- net_burstiness(g, by = "pair")
  expect_equal(p$source, c("a", "a", "b"))
  expect_equal(p$target, c("b", "c", "c"))
  expect_equal(p$events, c(2, 2, 3))
  expect_equal(p$burstiness, c(NaN, NaN, -1))
  pi <- net_inter_event(g, by = "pair")
  expect_equal(pi$gap, c(1, 8, 5, 5))
  d <- net_temporal(ev)
  expect_equal(net_burstiness(d, mode = "out")$events, c(4, 3, 0))
  expect_equal(net_burstiness(d, mode = "in")$events, c(0, 2, 5))
})

test_that("earliest arrival agrees with tsna on instantaneous contacts", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("network")
  # Simultaneous contacts chain, repeated contacts, a loop, and an isolate.
  ev <- data.frame(source = c(1, 2, 3, 1, 4, 2), target = c(2, 3, 4, 2, 4, 1),
                   time = c(1, 1, 3, 5, 5, 6))
  g <- net_temporal(ev, nodes = data.frame(id = 1:5))
  base <- network::network.initialize(5, directed = TRUE, loops = TRUE, multiple = FALSE)
  nd <- networkDynamic::networkDynamic(
    base.net = base,
    edge.spells = data.frame(onset = ev$time, terminus = ev$time, tail = ev$source,
                             head = ev$target),
    verbose = FALSE)
  ours <- net_temporal_distances(g, threads = 1)
  for (v in 1:5) {
    p <- tsna::tPath(nd, v = v, direction = "fwd", type = "earliest.arrive", start = 1,
                     end = 7, graph.step.time = 0)
    expect_equal(unname(ours[v, ]), ifelse(is.finite(p$tdist), 1 + p$tdist, Inf),
                 info = sprintf("origin %d", v))
  }
})
