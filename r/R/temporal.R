# Temporal networks: the constructor, slices and windows, aggregation,
# time-respecting paths and reachability, temporal closeness, rolling
# measures, and burstiness. Plots are in plot-temporal.R.

# Temporal networks (man/temporal_network.Rd).
temporal_network <- function(edges, nodes = NULL, start = "time", end = NULL,
                             duration = NULL, directed = TRUE, from = 1L, to = 2L, id = 1L,
                             node_start = NULL, node_end = NULL, weight = NULL,
                             weight_type = c("distance", "strength"), graph = list()) {
  weight_type <- match.arg(weight_type)
  default_start <- missing(start)
  if (inherits(edges, "netframe")) {
    check_netframe(edges)
    if (!is.null(nodes)) {
      stop("`nodes` must be NULL when `edges` is a netframe; it already has nodes.",
           call. = FALSE)
    }
    x <- edges
    if (default_start && !is.null(x$roles$start)) {
      start <- x$roles$start
      if (missing(end) && missing(duration)) {
        end <- x$roles$end
        duration <- x$roles$duration
      }
      default_start <- FALSE
    }
    if (is.null(node_start) && is.null(node_end)) {
      node_start <- x$roles$node_start
      node_end <- x$roles$node_end
    }
  } else {
    x <- netframe(edges, nodes, directed = directed, from = from, to = to, id = id,
                  weight = weight, weight_type = weight_type, graph = graph)
  }
  if (default_start && !"time" %in% names(x$edges) && "start" %in% names(x$edges)) {
    start <- "start"
  }
  if (!is.null(end) && !is.null(duration)) {
    stop("Give `end` or `duration`, not both.", call. = FALSE)
  }
  start <- resolve_column(x$edges, start, "start", "edges")
  kind <- time_kind(x$edges[[start]])
  if (is.na(kind)) {
    stop(sprintf(paste0("Start times (column `%s` of `edges`) must be numbers, dates, ",
                        "or date-times, not %s."), start, describe_class(x$edges[[start]])),
         call. = FALSE)
  }
  if (!is.null(end)) {
    end <- resolve_column(x$edges, end, "end", "edges")
    check_same_kind(x$edges[[end]], kind, sprintf("End times (column `%s` of `edges`)", end))
  }
  if (!is.null(duration)) {
    duration <- resolve_column(x$edges, duration, "duration", "edges")
    d <- x$edges[[duration]]
    if (!(inherits(d, "difftime") || (is.numeric(d) && !is.object(d)))) {
      stop(sprintf("Durations (column `%s` of `edges`) must be numbers or a difftime, not %s.",
                   duration, describe_class(d)), call. = FALSE)
    }
  }
  if (!is.null(node_start)) {
    node_start <- resolve_column(x$nodes, node_start, "node_start", "nodes")
    check_same_kind(x$nodes[[node_start]], kind,
                    sprintf("Node start times (column `%s` of `nodes`)", node_start))
  }
  if (!is.null(node_end)) {
    node_end <- resolve_column(x$nodes, node_end, "node_end", "nodes")
    check_same_kind(x$nodes[[node_end]], kind,
                    sprintf("Node end times (column `%s` of `nodes`)", node_end))
  }
  v <- x$edges[[start]]
  tz <- if (kind == "POSIXct") {
    z <- attr(as.POSIXct(v), "tzone")
    if (is.null(z)) "" else z[[1L]]
  }
  x$roles$start <- start
  x$roles$end <- end
  x$roles$duration <- duration
  x$roles$node_start <- node_start
  x$roles$node_end <- node_end
  x$roles$time_class <- kind
  x$roles$tz <- tz
  temporal_times(x)
  if (!is.null(node_start) || !is.null(node_end)) node_presence(x)
  x
}

# See man/temporal_network.Rd.
is_temporal <- function(x) inherits(x, "netframe") && !is.null(x$roles$start)

# A temporal network at a time or in a window (man/network_at.Rd).
network_at <- function(x, time, nodes = c("all", "present", "active")) {
  check_temporal(x)
  nodes <- match.arg(nodes)
  t <- time_arg(time, x, "time")
  ptr <- engine_temporal(x)
  idx <- engine_call("network_at", .Call(C_temporal_active_at, ptr, t)) + 1
  present <- if (nodes == "present") node_present(x, t, t) else NULL
  slice_netframe(x, idx, nodes, present)
}

# See man/network_at.Rd.
network_window <- function(x, from, to, nodes = c("all", "present", "active")) {
  check_temporal(x)
  nodes <- match.arg(nodes)
  a <- time_arg(from, x, "from")
  b <- time_arg(to, x, "to")
  if (!(a < b)) stop("`from` must be earlier than `to`.", call. = FALSE)
  ptr <- engine_temporal(x)
  w <- engine_call("network_window", .Call(C_temporal_windows, ptr, a, b))
  present <- if (nodes == "present") node_present(x, a, b) else NULL
  slice_netframe(x, w[[2L]] + 1, nodes, present)
}

# Snapshots of a temporal network (man/snapshots.Rd).
snapshots <- function(x, width = NULL, step = width, breaks = NULL, from = NULL, to = NULL,
                      nodes = c("all", "present", "active")) {
  check_temporal(x)
  nodes <- match.arg(nodes)
  w <- window_series(x, width, step, breaks, from, to)
  ptr <- engine_temporal(x)
  res <- engine_call("snapshots", .Call(C_temporal_windows, ptr, w$from, w$to))
  offsets <- res[[1L]]
  idx <- res[[2L]] + 1
  out <- lapply(seq_along(w$from), function(k) {
    keep <- if (offsets[k + 1L] > offsets[k]) idx[(offsets[k] + 1):offsets[k + 1L]] else numeric()
    present <- if (nodes == "present") node_present(x, w$from[k], w$to[k]) else NULL
    slice_netframe(x, keep, nodes, present)
  })
  attr(out, "windows") <- data.frame(start = from_seconds(w$from, x),
                                     end = from_seconds(w$to, x))
  out
}

# Aggregate a temporal network into a weighted network (man/aggregate_network.Rd).
aggregate_network <- function(x, weight = c("count", "duration"), from = NULL, to = NULL) {
  check_temporal(x)
  weight <- match.arg(weight)
  if (is.null(from) != is.null(to)) {
    stop("Give both `from` and `to`, or neither.", call. = FALSE)
  }
  window <- NULL
  if (!is.null(from)) {
    window <- c(time_arg(from, x, "from"), time_arg(to, x, "to"))
    if (!(window[1L] < window[2L])) stop("`from` must be earlier than `to`.", call. = FALSE)
  }
  ptr <- engine_temporal(x)
  a <- engine_call("aggregate_network", .Call(C_temporal_aggregate, ptr, window))
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(ids[a[[1L]] + 1], ids[a[[2L]] + 1], stringsAsFactors = FALSE)
  ends <- c(x$roles$from, x$roles$to)
  if (any(ends %in% c("count", "duration", "first", "last"))) ends <- c("from", "to")
  names(edges) <- ends
  edges$count <- a[[3L]]
  edges$duration <- a[[4L]] / time_scale(x)
  edges$first <- from_seconds(a[[5L]], x)
  edges$last <- from_seconds(a[[6L]], x)
  netframe(edges, x$nodes, directed = x$directed, from = 1L, to = 2L, id = x$roles$id,
           weight = weight, weight_type = "strength", graph = x$graph)
}

# Time-respecting paths (man/temporal_distances.Rd).
temporal_distances <- function(x, type = c("earliest_arrival", "latest_departure", "fastest",
                                           "shortest", "average"),
                               from = NULL, to = NULL, start = NULL, end = NULL,
                               min_wait = 0, max_wait = Inf, strict = FALSE, traversal = 0,
                               threads = 0L) {
  check_temporal(x)
  type <- match.arg(type)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  n <- nrow(x$nodes)
  rows <- if (is.null(from)) seq_len(n) else node_positions(x, from, "from", empty_ok = FALSE)
  cols <- if (is.null(to)) seq_len(n) else node_positions(x, to, "to", empty_ok = FALSE)
  ptr <- engine_temporal(x)
  scale <- time_scale(x)
  if (type == "average") {
    r <- engine_call("temporal_distances", .Call(C_temporal_closeness, ptr, opts, TRUE))
    m <- matrix(r[[2L]], n, n, byrow = TRUE)[rows, cols, drop = FALSE] / scale
  } else if (type == "latest_departure") {
    # One pass per destination, walks into it: row k holds the latest
    # departure from every node to destination cols[k].
    v <- engine_call("temporal_distances",
                     .Call(C_temporal_paths, ptr, 1L, 2L, as.integer(cols - 1L), opts))
    m <- t(matrix(v, length(cols), n, byrow = TRUE))[rows, , drop = FALSE] / scale
  } else {
    code <- switch(type, earliest_arrival = 0L, fastest = 2L, shortest = 3L)
    v <- engine_call("temporal_distances",
                     .Call(C_temporal_paths, ptr, code, 1L, as.integer(rows - 1L), opts))
    m <- matrix(v, length(rows), n, byrow = TRUE)[, cols, drop = FALSE]
    if (type != "shortest") m <- m / scale
  }
  ids <- as.character(x$nodes[[x$roles$id]])
  dimnames(m) <- list(ids[rows], ids[cols])
  m
}

# Temporal reachability (man/temporal_reach.Rd).
temporal_reach <- function(x, mode = c("out", "in"), start = NULL, end = NULL, min_wait = 0,
                           max_wait = Inf, strict = FALSE, traversal = 0, threads = 0L) {
  check_temporal(x)
  mode <- match.arg(mode)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  engine_call("temporal_reach", .Call(C_temporal_reach, engine_temporal(x),
                                      if (mode == "out") 1L else 2L, opts))
}

# See man/temporal_reach.Rd.
reachability_ratio <- function(x, start = NULL, end = NULL, min_wait = 0, max_wait = Inf,
                               strict = FALSE, traversal = 0, threads = 0L) {
  check_temporal(x)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  n <- nrow(x$nodes)
  if (n < 2L) return(NaN)
  r <- engine_call("reachability_ratio", .Call(C_temporal_reach, engine_temporal(x), 1L, opts))
  sum(r) / (n * (n - 1))
}

# Temporal closeness (man/temporal_closeness.Rd).
temporal_closeness <- function(x, start = NULL, end = NULL, min_wait = 0, max_wait = Inf,
                               strict = FALSE, traversal = 0, threads = 0L) {
  check_temporal(x)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  r <- engine_call("temporal_closeness",
                   .Call(C_temporal_closeness, engine_temporal(x), opts, FALSE))
  r[[1L]] * time_scale(x)
}

# Apply a measure over a series of windows (man/rolling.Rd).
rolling <- function(x, fun, width = NULL, step = width, ..., breaks = NULL, from = NULL,
                    to = NULL, nodes = c("all", "present", "active")) {
  check_temporal(x)
  fun <- match.fun(fun)
  nodes <- match.arg(nodes)
  snaps <- snapshots(x, width = width, step = step, breaks = breaks, from = from, to = to,
                     nodes = nodes)
  w <- attr(snaps, "windows")
  parts <- lapply(seq_along(snaps), function(k) {
    s <- snaps[[k]]
    v <- fun(s, ...)
    if (!is.numeric(v) && !is.logical(v)) {
      stop(sprintf("`fun` must return numbers; it returned %s for window %d.",
                   describe_class(v), k), call. = FALSE)
    }
    v <- unname(as.vector(v))
    ids <- s$nodes[[s$roles$id]]
    if (length(v) == 1L && nrow(s$nodes) != 1L) {
      list(k = k, node = NA, value = as.double(v))
    } else if (length(v) == nrow(s$nodes)) {
      list(k = rep(k, length(v)), node = ids, value = as.double(v))
    } else {
      stop(sprintf(paste0("`fun` must return one value per node or a single value; it ",
                          "returned %d values for window %d, which has %d nodes."),
                   length(v), k, nrow(s$nodes)), call. = FALSE)
    }
  })
  k <- unlist(lapply(parts, `[[`, "k"))
  node <- do.call(c, lapply(parts, function(p) {
    if (is.factor(p$node)) as.character(p$node) else p$node
  }))
  value <- unlist(lapply(parts, `[[`, "value"))
  if (is.null(k)) k <- integer()
  if (is.null(node)) node <- character()
  if (is.null(value)) value <- numeric()
  data.frame(start = w$start[k], end = w$end[k], node = node, value = value,
             stringsAsFactors = FALSE)
}

# Inter-event times and burstiness (man/burstiness.Rd).
burstiness <- function(x, by = c("node", "pair"), mode = c("all", "out", "in")) {
  check_temporal(x)
  by <- match.arg(by)
  mode <- match.arg(mode)
  ptr <- engine_temporal(x)
  r <- engine_call("burstiness", .Call(C_temporal_burstiness, ptr,
                                       if (by == "node") 0L else 1L, mode_code(mode)))
  out <- event_groups(x, ptr, by, seq_along(r[[1L]]))
  out$events <- r[[1L]]
  out$burstiness <- r[[2L]]
  out$memory <- r[[3L]]
  out
}

# See man/burstiness.Rd.
inter_event_times <- function(x, by = c("node", "pair"), mode = c("all", "out", "in")) {
  check_temporal(x)
  by <- match.arg(by)
  mode <- match.arg(mode)
  ptr <- engine_temporal(x)
  r <- engine_call("inter_event_times", .Call(C_temporal_inter_event, ptr,
                                              if (by == "node") 0L else 1L, mode_code(mode)))
  out <- event_groups(x, ptr, by, r[[1L]] + 1)
  out$gap <- r[[2L]] / time_scale(x)
  out
}

# ---- Internal helpers ----------------------------------------------------------

check_temporal <- function(x) {
  check_netframe(x)
  if (is.null(x$roles$start)) {
    stop(paste0("`x` must be a temporal network, with times on its edges; create one ",
                "with temporal_network()."), call. = FALSE)
  }
}

# "numeric", "Date", "POSIXct", or NA for anything else.
time_kind <- function(v) {
  if (inherits(v, "POSIXct") || inherits(v, "POSIXlt")) return("POSIXct")
  if (inherits(v, "Date")) return("Date")
  if (is.numeric(v) && !is.object(v)) return("numeric")
  NA_character_
}

check_same_kind <- function(v, kind, what) {
  k <- time_kind(v)
  if (is.na(k) || k != kind) {
    want <- switch(kind, numeric = "numbers", Date = "dates", POSIXct = "date-times")
    stop(sprintf("%s must be %s, like the start times, not %s.", what, want,
                 describe_class(v)), call. = FALSE)
  }
}

# Seconds per unit of the data's times: 86400 for dates, 1 otherwise.
time_scale <- function(x) if (identical(x$roles$time_class, "Date")) 86400 else 1

# Times of the data's type as seconds (or raw numbers for numeric times).
to_seconds <- function(v, kind) {
  switch(kind,
    numeric = as.double(v),
    Date = as.double(unclass(v)) * 86400,
    POSIXct = as.double(unclass(as.POSIXct(v)))
  )
}

# Seconds back to the data's type.
from_seconds <- function(s, x) {
  switch(x$roles$time_class,
    numeric = s,
    Date = structure(s / 86400, class = "Date"),
    POSIXct = .POSIXct(s, tz = x$roles$tz)
  )
}

# Start and end times (seconds), checked. Instantaneous contacts end when
# they start.
temporal_times <- function(x) {
  kind <- x$roles$time_class
  s <- to_seconds(x$edges[[x$roles$start]], kind)
  bad <- which(!is.finite(s))
  if (length(bad)) {
    stop(sprintf(paste0("Start times (column `%s` of `x$edges`) must not be missing or ",
                        "infinite; %s %s one."), x$roles$start, rows_text(bad),
                 if (length(bad) == 1L) "has" else "have"), call. = FALSE)
  }
  e <- NULL
  if (!is.null(x$roles$end)) {
    e <- to_seconds(x$edges[[x$roles$end]], kind)
    what <- sprintf("End times (column `%s` of `x$edges`)", x$roles$end)
  } else if (!is.null(x$roles$duration)) {
    d <- x$edges[[x$roles$duration]]
    d <- if (inherits(d, "difftime")) {
      as.double(d, units = "secs")
    } else {
      as.double(d) * time_scale(x)
    }
    e <- s + d
    what <- sprintf("Durations (column `%s` of `x$edges`)", x$roles$duration)
  }
  if (!is.null(e)) {
    bad <- which(!is.finite(e))
    if (length(bad)) {
      stop(sprintf("%s must not be missing or infinite; %s %s one.", what, rows_text(bad),
                   if (length(bad) == 1L) "has" else "have"), call. = FALSE)
    }
    bad <- which(e < s)
    if (length(bad)) {
      stop(sprintf("Contacts must not end before they start; see %s of `x$edges`.",
                   rows_text(bad)), call. = FALSE)
    }
  }
  list(start = s, end = e)
}

# Node presence intervals (seconds), with -Inf and Inf for missing values.
node_presence <- function(x) {
  n <- nrow(x$nodes)
  kind <- x$roles$time_class
  a <- if (is.null(x$roles$node_start)) rep(-Inf, n) else to_seconds(x$nodes[[x$roles$node_start]], kind)
  b <- if (is.null(x$roles$node_end)) rep(Inf, n) else to_seconds(x$nodes[[x$roles$node_end]], kind)
  a[is.na(a)] <- -Inf
  b[is.na(b)] <- Inf
  bad <- which(b < a)
  if (length(bad)) {
    stop(sprintf("Nodes must not leave before they arrive; see %s of `x$nodes`.",
                 rows_text(bad)), call. = FALSE)
  }
  list(start = a, end = b)
}

# Which nodes are present at time a (a == b) or in the window [a, b).
node_present <- function(x, a, b) {
  if (is.null(x$roles$node_start) && is.null(x$roles$node_end)) {
    return(rep(TRUE, nrow(x$nodes)))
  }
  p <- node_presence(x)
  if (a == b) {
    (p$start <= a & a < p$end) | (p$start == a & p$end == a)
  } else {
    ifelse(p$start == p$end, p$start >= a & p$start < b, p$start < b & p$end > a)
  }
}

# The engine's temporal network for `x`, built on first use and cached
# alongside the graph handle (see engine_graph()).
engine_temporal <- function(x) {
  check_index(x)
  tt <- temporal_times(x)
  cache <- x$cache
  ptr <- cache$temporal_ptr
  fresh <- !is.null(ptr) && .Call(C_graph_alive, ptr) &&
    identical(cache$temporal_n, nrow(x$nodes)) && identical(cache$temporal_times, tt) &&
    identical(cache$temporal_directed, x$directed)
  if (!fresh) {
    m <- length(x$index$src)
    if (m && max(x$index$src, x$index$tgt) >= nrow(x$nodes)) index_error()
    ptr <- engine_call("temporal_network",
      .Call(C_temporal_new, nrow(x$nodes), x$index$src, x$index$tgt, tt$start, tt$end,
            x$directed))
    cache$temporal_ptr <- ptr
    cache$temporal_n <- nrow(x$nodes)
    cache$temporal_times <- tt
    cache$temporal_directed <- x$directed
  }
  ptr
}

# A time argument as seconds: numbers for numeric times; dates or
# date-times for dates and date-times.
time_arg <- function(value, x, name) {
  kind <- x$roles$time_class
  k <- time_kind(value)
  ok <- length(value) == 1L && !is.na(k) &&
    (k == kind || (kind != "numeric" && k != "numeric"))
  if (!ok) {
    want <- switch(kind, numeric = "a number", Date = "a date", POSIXct = "a date-time")
    stop(sprintf("`%s` must be a single time: %s, like the network's times.", name, want),
         call. = FALSE)
  }
  s <- to_seconds(value, k)
  if (is.na(s)) stop(sprintf("`%s` must not be missing.", name), call. = FALSE)
  s
}

# A length of time as seconds: a difftime, or a number in the data's unit.
duration_arg <- function(value, x, name, allow_inf = FALSE) {
  if (inherits(value, "difftime")) {
    s <- as.double(value, units = "secs")
  } else if (is.numeric(value) && !is.object(value)) {
    s <- as.double(value) * time_scale(x)
  } else {
    s <- NA_real_
  }
  if (length(s) != 1L || is.na(s) || s < 0 || (!allow_inf && !is.finite(s))) {
    stop(sprintf("`%s` must be a single non-negative %snumber or difftime.", name,
                 if (allow_inf) "" else "finite "), call. = FALSE)
  }
  s
}

# Engine path options: c(start, end, min_wait, max_wait, strict, traversal,
# traversal_time, threads), seconds for times.
path_options <- function(x, start, end, min_wait, max_wait, strict, traversal, threads) {
  s <- if (is.null(start)) NA_real_ else time_arg(start, x, "start")
  e <- if (is.null(end)) NA_real_ else time_arg(end, x, "end")
  if (!is.na(s) && !is.na(e) && s > e) {
    stop("`start` must not be later than `end`.", call. = FALSE)
  }
  lo <- duration_arg(min_wait, x, "min_wait")
  hi <- duration_arg(max_wait, x, "max_wait", allow_inf = TRUE)
  if (hi < lo) stop("`max_wait` must not be smaller than `min_wait`.", call. = FALSE)
  strict <- check_flag(strict, "strict")
  if (is.character(traversal)) {
    if (length(traversal) != 1L || is.na(traversal) || traversal != "duration") {
      stop("`traversal` must be a number, a difftime, or \"duration\".", call. = FALSE)
    }
    tcode <- 1
    tval <- 0
  } else {
    tcode <- 0
    tval <- duration_arg(traversal, x, "traversal")
  }
  c(s, e, lo, hi, as.double(strict), tcode, tval, as.double(resolve_threads(threads)))
}

# Window starts and ends (seconds) for snapshots() and rolling().
window_series <- function(x, width, step, breaks, from, to) {
  if (!is.null(breaks)) {
    if (!is.null(width) || !is.null(from) || !is.null(to)) {
      stop("Give `breaks`, or `width` (with `step`, `from`, and `to`), not both.",
           call. = FALSE)
    }
    b <- vapply(seq_along(breaks), function(k) time_arg(breaks[k], x, "breaks"), numeric(1))
    if (length(b) < 2L || any(diff(b) <= 0)) {
      stop("`breaks` must hold at least two strictly increasing times.", call. = FALSE)
    }
    return(list(from = b[-length(b)], to = b[-1L]))
  }
  if (is.null(width)) stop("Give `width` or `breaks`.", call. = FALSE)
  w <- duration_arg(width, x, "width")
  st <- duration_arg(step, x, "step")
  if (w <= 0 || st <= 0) stop("`width` and `step` must be positive.", call. = FALSE)
  tt <- temporal_times(x)
  last <- if (is.null(tt$end)) tt$start else tt$end
  a <- if (is.null(from)) {
    if (length(tt$start)) min(tt$start) else 0
  } else {
    time_arg(from, x, "from")
  }
  b <- if (is.null(to)) {
    if (length(last)) max(last) else a
  } else {
    time_arg(to, x, "to")
  }
  if (b < a) stop("`to` must not be earlier than `from`.", call. = FALSE)
  count <- floor((b - a) / st) + 1
  if (count > 1e7) {
    stop(sprintf("That gives %.0f windows; use a larger `step`.", count), call. = FALSE)
  }
  starts <- a + (seq_len(count) - 1) * st
  starts <- starts[starts <= b]
  list(from = starts, to = starts + w)
}

# The contacts at 1-based positions `keep` (ascending), with the nodes
# chosen by `nodes`, as a temporal netframe with every column and role.
slice_netframe <- function(x, keep, nodes, present) {
  keep <- as.integer(keep)
  node_keep <- switch(nodes,
    all = rep(TRUE, nrow(x$nodes)),
    present = present,
    active = rep(FALSE, nrow(x$nodes))
  )
  node_keep[x$index$src[keep] + 1L] <- TRUE
  node_keep[x$index$tgt[keep] + 1L] <- TRUE
  nodes_df <- x$nodes[node_keep, , drop = FALSE]
  edges_df <- x$edges[keep, , drop = FALSE]
  rownames(nodes_df) <- NULL
  rownames(edges_df) <- NULL
  weight_type <- if (is.null(x$roles$weight_type)) "distance" else x$roles$weight_type
  out <- netframe(edges_df, nodes_df, directed = x$directed, from = x$roles$from,
                  to = x$roles$to, id = x$roles$id, weight = x$roles$weight,
                  weight_type = weight_type, graph = x$graph)
  extra <- setdiff(names(x$roles), names(out$roles))
  if (length(extra)) out$roles[extra] <- x$roles[extra]
  out
}

# Node or pair identifiers for rows of event statistics (1-based groups).
event_groups <- function(x, ptr, by, groups) {
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (by == "node") {
    return(data.frame(node = ids[groups], stringsAsFactors = FALSE))
  }
  a <- engine_call("burstiness", .Call(C_temporal_aggregate, ptr, NULL))
  data.frame(from = ids[a[[1L]][groups] + 1], to = ids[a[[2L]][groups] + 1],
             stringsAsFactors = FALSE)
}

# One line for print.netframe().
temporal_summary <- function(x) {
  tt <- tryCatch(temporal_times(x), error = function(e) NULL)
  if (is.null(tt)) {
    return(sprintf("temporal: times in `%s` (some are invalid)\n", x$roles$start))
  }
  m <- length(tt$start)
  cols <- sprintf("`%s`", x$roles$start)
  if (!is.null(x$roles$end)) cols <- sprintf("%s to `%s`", cols, x$roles$end)
  if (!is.null(x$roles$duration)) cols <- sprintf("%s plus `%s`", cols, x$roles$duration)
  if (m == 0L) return(sprintf("temporal: no contacts (times in %s)\n", cols))
  last <- if (is.null(tt$end)) max(tt$start) else max(tt$end)
  inst <- if (is.null(tt$end)) m else sum(tt$end == tt$start)
  sprintf("temporal: contacts from %s to %s (times in %s), %s instantaneous\n",
          format(from_seconds(min(tt$start), x)), format(from_seconds(last, x)), cols,
          if (inst == m) "all" else format(inst))
}
