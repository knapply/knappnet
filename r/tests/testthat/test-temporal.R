# Temporal networks: construction, slicing, aggregation, time-respecting
# paths, reachability, temporal closeness, rolling measures, and event
# statistics. Values are checked by hand; the full comparison against a
# brute-force search over walks lives in validation/.

# a -> b and b -> c at time 1 (they chain unless times must increase),
# b -> c again at 3, and a -> c at 6.
msgs <- function() {
  temporal_network(data.frame(from = c("a", "b", "b", "a"), to = c("b", "c", "c", "c"),
                              time = c(1, 1, 3, 6)))
}

test_that("temporal networks are netframes with time roles", {
  g <- msgs()
  expect_s3_class(g, "netframe")
  expect_true(is_temporal(g))
  expect_false(is_temporal(netframe(data.frame(from = "a", to = "b"))))
  expect_false(is_temporal(data.frame(time = 1)))
  expect_equal(g$roles$start, "time")
  expect_equal(g$roles$time_class, "numeric")
  expect_equal(g$nodes$name, c("a", "b", "c"))
  expect_output(print(g), "temporal: contacts from 1 to 6")
  expect_output(print(g), "all instantaneous")
  # Static measures see every contact as an edge.
  expect_equal(degree(g), c(2, 3, 3))
})

test_that("start, end, and duration columns are found and checked", {
  d <- data.frame(from = c("a", "b"), to = c("b", "c"), start = c(0, 5), stop = c(3, 5),
                  len = c(3, 0))
  g <- temporal_network(d, end = "stop")
  expect_equal(g$roles$start, "start")
  expect_output(print(g), "1 instantaneous")
  expect_equal(nrow(network_at(g, 2)$edges), 1L)
  h <- temporal_network(d, duration = "len")
  expect_equal(nrow(network_at(h, 2)$edges), 1L)
  expect_equal(aggregate_network(h)$edges$duration, c(3, 0))
  expect_error(temporal_network(d, end = "stop", duration = "len"), "not both")
  bad <- data.frame(from = c("a", "b"), to = c("b", "c"), time = c(2, 1), end = c(1, 1))
  expect_error(temporal_network(bad, end = "end"), "end before they start; see row 1")
  expect_error(temporal_network(data.frame(from = "a", to = "b", time = NA_real_)),
               "must not be missing")
  expect_error(temporal_network(data.frame(from = "a", to = "b", when = 1)), "`start`")
  expect_error(temporal_network(data.frame(from = "a", to = "b", time = "noon")),
               "numbers, dates, or date-times")
  mixed <- data.frame(from = "a", to = "b", time = 1, end = as.Date("2024-01-01"))
  expect_error(temporal_network(mixed, end = "end"), "must be numbers")
  expect_error(temporal_reach(netframe(data.frame(from = "a", to = "b"))),
               "temporal network")
})

test_that("a netframe whose edges carry a time column can be made temporal", {
  x <- netframe(data.frame(src = c("a", "b"), dst = c("b", "c"), start = c(1, 2),
                           bytes = c(10, 20)), weight = "bytes")
  g <- temporal_network(x)
  expect_true(is_temporal(g))
  expect_equal(g$roles$weight, "bytes")
  expect_equal(temporal_reach(g), c(2, 1, 0))
  expect_identical(temporal_network(g)$roles, g$roles)
  expect_error(temporal_network(x, nodes = data.frame(id = "a")), "must be NULL")
})

test_that("slicing follows half-open intervals", {
  calls <- data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d"),
                      time = c(1, 2, 2, 5), end = c(3, 2, 4, 6))
  g <- temporal_network(calls, end = "end", directed = FALSE)
  at2 <- network_at(g, 2)
  expect_true(is_temporal(at2))
  expect_equal(nrow(at2$edges), 3L)
  expect_equal(nrow(at2$nodes), 4L)
  expect_equal(nrow(network_at(g, 3)$edges), 1L)
  w <- network_window(g, 3, 5, nodes = "active")
  expect_equal(w$edges$from, "b")
  expect_equal(w$nodes$name, c("b", "c"))
  expect_equal(nrow(network_window(g, 0, 10)$edges), 4L)
  expect_error(network_window(g, 5, 5), "earlier")
  expect_error(network_at(g, as.Date("2024-01-01")), "single time")
})

test_that("node presence intervals select nodes", {
  nodes <- data.frame(name = c("a", "b", "c"), joined = c(0, 5, NA))
  g <- temporal_network(data.frame(from = c("a", "b", "b", "a"), to = c("b", "c", "c", "c"),
                                   time = c(1, 1, 3, 6)),
                        nodes = nodes, node_start = "joined")
  expect_equal(network_at(g, 2, nodes = "present")$nodes$name, c("a", "c"))
  expect_equal(nrow(network_at(g, 2, nodes = "active")$nodes), 0L)
  expect_equal(network_at(g, 2)$nodes$name, c("a", "b", "c"))
  # An endpoint of an active contact is always kept.
  expect_equal(network_at(g, 1, nodes = "present")$nodes$name, c("a", "b", "c"))
})

test_that("snapshots and rolling measures", {
  logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2"),
                       host = c("h1", "h1", "h2", "h2", "h3"), time = c(0, 1, 5, 6, 12))
  g <- temporal_network(logins)
  s <- snapshots(g, width = 5)
  expect_equal(attr(s, "windows")$start, c(0, 5, 10))
  expect_equal(attr(s, "windows")$end, c(5, 10, 15))
  expect_equal(vapply(s, function(w) nrow(w$edges), numeric(1)), c(2, 2, 1))
  slide <- snapshots(g, width = 6, step = 3)
  expect_equal(vapply(slide, function(w) nrow(w$edges), numeric(1)), c(3, 2, 1, 1, 1))
  b <- snapshots(g, breaks = c(0, 6, 13))
  expect_equal(vapply(b, function(w) nrow(w$edges), numeric(1)), c(3, 2))
  expect_error(snapshots(g), "width")
  expect_error(snapshots(g, width = 5, breaks = c(0, 1)), "not both")
  expect_error(snapshots(g, breaks = c(3, 1)), "increasing")
  r <- rolling(g, degree, width = 5)
  expect_equal(names(r), c("start", "end", "node", "value"))
  expect_equal(nrow(r), 18L)
  expect_equal(r$value[r$start == 0], c(1, 1, 0, 2, 0, 0))
  expect_equal(r$node[1:6], c("u1", "u2", "u3", "h1", "h2", "h3"))
  dens <- rolling(g, function(w) nrow(w$edges), width = 5)
  expect_equal(dens$value, c(2, 2, 1))
  expect_true(all(is.na(dens$node)))
  expect_error(rolling(g, function(w) c(1, 2), width = 5), "one value per node")
})

test_that("aggregation counts and clips contacts", {
  calls <- data.frame(from = c("a", "b", "a", "c"), to = c("b", "a", "c", "a"),
                      time = c(0, 10, 20, 25), end = c(5, 12, 20, 40))
  g <- temporal_network(calls, end = "end", directed = FALSE)
  a <- aggregate_network(g)
  expect_false(is_temporal(a))
  expect_equal(a$roles$weight, "count")
  expect_equal(a$roles$weight_type, "strength")
  expect_equal(a$edges$from, c("a", "a"))
  expect_equal(a$edges$to, c("b", "c"))
  expect_equal(a$edges$count, c(2, 2))
  expect_equal(a$edges$duration, c(7, 15))
  expect_equal(a$edges$first, c(0, 20))
  expect_equal(a$edges$last, c(12, 40))
  w <- aggregate_network(g, weight = "duration", from = 0, to = 30)
  expect_equal(w$roles$weight, "duration")
  expect_equal(w$edges$duration, c(7, 5))
  expect_error(aggregate_network(g, from = 0), "both")
  d <- aggregate_network(temporal_network(calls, end = "end"))
  expect_equal(nrow(d$edges), 4L)
})

test_that("time-respecting paths", {
  g <- msgs()
  ids <- c("a", "b", "c")
  ea <- matrix(c(1, 1, 1, Inf, 1, 1, Inf, Inf, 1), 3, byrow = TRUE,
               dimnames = list(ids, ids))
  expect_equal(temporal_distances(g), ea)
  # Strictly increasing times: b -> c at 1 no longer follows a -> b at 1.
  expect_equal(unname(temporal_distances(g, from = "a", strict = TRUE)[1, ]), c(1, 1, 3))
  ld <- matrix(c(6, 1, 6, -Inf, 6, 3, -Inf, -Inf, 6), 3, byrow = TRUE,
               dimnames = list(ids, ids))
  expect_equal(temporal_distances(g, "latest_departure"), ld)
  expect_equal(temporal_distances(g, "latest_departure", from = "b", to = c("c", "b")),
               matrix(c(3, 6), 1, dimnames = list("b", c("c", "b"))))
  fast <- temporal_distances(g, "fastest")
  expect_equal(unname(fast[1, ]), c(0, 0, 0))
  expect_equal(unname(fast[3, ]), c(Inf, Inf, 0))
  expect_equal(unname(temporal_distances(g, "shortest")["a", ]), c(0, 1, 1))
  # Waiting at most 1 at b: the walk a -> b (1), b -> c (3) is not allowed.
  expect_equal(temporal_distances(g, from = "a", strict = TRUE, max_wait = 1)["a", "c"], 6)
  # One time unit to traverse each contact.
  expect_equal(unname(temporal_distances(g, from = "a", traversal = 1)[1, ]), c(1, 2, 4))
  expect_equal(temporal_distances(g, from = "a", start = 2)["a", "c"], 6)
  expect_equal(temporal_distances(g, from = "a", end = 2)["a", "c"], 1)
  avg <- temporal_distances(g, "average", start = 0, end = 10)
  expect_equal(unname(avg[1, ]), c(0, 5, 2.5))
  expect_equal(avg["b", "c"], 3.4)
  expect_error(temporal_distances(g, from = "z"), "not in the network")
  expect_error(temporal_distances(g, start = 5, end = 1), "later")
  expect_error(temporal_distances(g, min_wait = 2, max_wait = 1), "smaller")
  expect_error(temporal_distances(g, traversal = "flight"), "duration")
  expect_error(temporal_distances(g, min_wait = -1), "non-negative")
  expect_equal(temporal_distances(g, "fastest", threads = 1),
               temporal_distances(g, "fastest", threads = 2))
})

test_that("contact durations as traversal times", {
  # A 2-unit call from a to b, then b calls c at 1 (during the first call)
  # and at 3 (after it).
  d <- data.frame(from = c("a", "b", "b"), to = c("b", "c", "c"), time = c(0, 1, 3),
                  end = c(2, 2, 4))
  g <- temporal_network(d, end = "end")
  expect_equal(temporal_distances(g, from = "a")["a", "c"], 1)
  expect_equal(temporal_distances(g, from = "a", traversal = "duration")["a", "c"], 4)
  expect_equal(temporal_distances(g, "fastest", traversal = "duration")["a", "c"], 4)
})

test_that("reachability and temporal closeness", {
  g <- temporal_network(data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"),
                                   time = c(1, 2, 1)))
  expect_equal(temporal_reach(g), c(2, 1, 1, 0))
  expect_equal(temporal_reach(g, mode = "in"), c(0, 1, 2, 1))
  expect_equal(reachability_ratio(g), 1 / 3)
  expect_equal(temporal_reach(g, max_wait = 0.5), c(1, 1, 1, 0))
  one <- temporal_network(data.frame(from = "a", to = "a", time = 1))
  expect_true(is.nan(reachability_ratio(one)))
  m <- msgs()
  # a reaches b on average in 5 and c in 2.5 over starting times in [0, 10].
  expect_equal(temporal_closeness(m, start = 0, end = 10), c(0.3, 5 / 34, 0))
  expect_equal(temporal_closeness(m, start = 0, end = 10, strict = TRUE)[1], 17 / 70)
  expect_true(all(is.nan(temporal_closeness(m, start = 1, end = 1))))
})

test_that("dates and date-times", {
  d <- data.frame(from = c("a", "b"), to = c("b", "c"),
                  day = as.Date(c("2024-01-01", "2024-01-03")))
  g <- temporal_network(d, start = "day")
  expect_equal(g$roles$time_class, "Date")
  expect_output(print(g), "2024-01-01 to 2024-01-03")
  expect_equal(temporal_distances(g, "fastest")["a", "c"], 2)
  expect_equal(temporal_distances(g, from = "a")["a", "c"], 19725)
  expect_equal(temporal_distances(g, "fastest", max_wait = 1)["a", "c"], Inf)
  expect_equal(temporal_distances(g, "fastest", max_wait = as.difftime(2, units = "days"))["a", "c"], 2)
  s <- snapshots(g, width = 1)
  expect_s3_class(attr(s, "windows")$start, "Date")
  expect_equal(length(s), 3L)
  expect_equal(nrow(network_at(g, as.Date("2024-01-03"))$edges), 1L)
  expect_error(network_at(g, 19725), "a date")
  expect_equal(inter_event_times(g)$gap, 2)

  p <- data.frame(src = c("a", "b"), dst = c("b", "c"),
                  t = as.POSIXct(c("2024-05-01 10:00:00", "2024-05-01 11:00:00"), tz = "UTC"))
  h <- temporal_network(p, start = "t")
  expect_equal(h$roles$time_class, "POSIXct")
  expect_equal(temporal_distances(h, "fastest")["a", "c"], 3600)
  expect_equal(temporal_distances(h, "fastest",
                                  max_wait = as.difftime(30, units = "mins"))["a", "c"], Inf)
  w <- network_window(h, as.POSIXct("2024-05-01 10:30:00", tz = "UTC"),
                      as.POSIXct("2024-05-01 12:00:00", tz = "UTC"))
  expect_equal(nrow(w$edges), 1L)
  a <- aggregate_network(h)
  expect_s3_class(a$edges$first, "POSIXct")
  expect_equal(format(a$edges$first[2], tz = "UTC"), "2024-05-01 11:00:00")
})

test_that("inter-event times, burstiness, and memory", {
  ev <- data.frame(from = c("a", "a", "a", "a", "b", "b", "b"),
                   to = c("b", "b", "c", "c", "c", "c", "c"),
                   time = c(0, 1, 2, 10, 0, 5, 10))
  g <- temporal_network(ev, directed = FALSE)
  ie <- inter_event_times(g)
  expect_equal(ie$node, rep(c("a", "b", "c"), c(3, 4, 4)))
  expect_equal(ie$gap, c(1, 1, 8, 0, 1, 4, 5, 2, 3, 5, 0))
  b <- burstiness(g)
  expect_equal(b$node, c("a", "b", "c"))
  expect_equal(b$events, c(4, 5, 5))
  gaps <- c(1, 1, 8)
  s <- sqrt(mean((gaps - mean(gaps))^2))
  expect_equal(b$burstiness[1], (s - mean(gaps)) / (s + mean(gaps)))
  expect_true(is.nan(b$memory[1]))
  # The memory coefficient is the correlation of consecutive gaps.
  expect_equal(b$memory[2], cor(c(0, 1, 4), c(1, 4, 5)))
  p <- burstiness(g, by = "pair")
  expect_equal(p$from, c("a", "a", "b"))
  expect_equal(p$to, c("b", "c", "c"))
  expect_equal(p$events, c(2, 2, 3))
  expect_equal(p$burstiness, c(NaN, NaN, -1))
  pi <- inter_event_times(g, by = "pair")
  expect_equal(pi$gap, c(1, 8, 5, 5))
  d <- temporal_network(ev)
  expect_equal(burstiness(d, mode = "out")$events, c(4, 3, 0))
  expect_equal(burstiness(d, mode = "in")$events, c(0, 2, 5))
})
