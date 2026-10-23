#' Temporal networks
#'
#' `temporal_network()` creates a network whose edges are timed contacts:
#' each edge has a start time and, optionally, an end time or a duration,
#' and nodes may have presence intervals. The result is an ordinary
#' [netframe()] with extra column roles, so every static measure still
#' works on it (treating each contact as an edge), and the temporal
#' functions listed below use the times. `is_temporal()` tells whether a
#' netframe has times.
#'
#' @param edges A data frame of contacts, one row per contact, or a
#'   [netframe()] whose edge table holds the times; its nodes, endpoints,
#'   weights, and other roles are kept.
#' @param nodes A data frame of nodes, or `NULL`; see [netframe()]. Must be
#'   `NULL` when `edges` is a netframe.
#' @param start Column of `edges` holding each contact's start time, by
#'   name or position. Times may be numbers (in any unit), `Date`, or
#'   `POSIXct` date-times. By default the column named `"time"`, or, if
#'   there is none, the column named `"start"`; a netframe that already
#'   has times keeps them.
#' @param end Optional column of end times, of the same type as `start`.
#' @param duration Optional column of durations, an alternative to `end`:
#'   numbers in the unit of the times (seconds for date-times, days for
#'   dates) or a `difftime`. Without `end` or `duration`, every contact is
#'   instantaneous.
#' @param directed Whether contacts are directed. Ignored when `edges` is a
#'   netframe, which keeps its own direction.
#' @param from,to,id,weight,weight_type,graph As in [netframe()].
#' @param node_start,node_end Optional columns of `nodes` holding when each
#'   node is present, of the same type as `start`. A missing value means
#'   present from the beginning (start) or until the end (end).
#' @param x Any object.
#'
#' @details
#' **Times.** Times are stored as given in `x$edges` and `x$nodes`; the
#' column roles are recorded in `x$roles` (`start`, `end` or `duration`,
#' `node_start`, `node_end`, and the time type). Internally dates and
#' date-times are converted to seconds since 1970-01-01 UTC, so a date is
#' its midnight in UTC. Arguments that are times (such as `from` and `to`
#' of [network_window()]) take the type of the data: numbers for numeric
#' times, and dates or date-times for dates and date-times. Arguments that
#' are lengths of time (window widths, waiting times) take a `difftime`, or
#' a number in the unit of the data: seconds for date-times, days for
#' dates. Results that are lengths of time use the same unit.
#'
#' **When a contact is active.** Intervals are half-open, as the spells of
#' networkDynamic (Butts et al. 2026): a contact from \eqn{s} to \eqn{e > s}
#' is active at times \eqn{s \le t < e}{s <= t < e}, and an instantaneous
#' contact (\eqn{e = s}) only at \eqn{t = s}. A contact is active in a
#' window \eqn{[a, b)} when it is active at some time in it: \eqn{s < b} and
#' \eqn{e > a} for an interval, \eqn{a \le s < b}{a <= s < b} for an
#' instantaneous contact. So back-to-back contacts \eqn{[1, 2)} and
#' \eqn{[2, 3)} never overlap, consecutive windows of a tumbling series never
#' share an instantaneous contact, and an interval belongs to every window it
#' overlaps. Node presence intervals follow the same rule.
#'
#' **Contacts and edges.** Self-loops and repeated contacts between the same
#' nodes are allowed and kept; each contact is one edge of the netframe.
#' Static measures applied directly to a temporal network see every contact
#' as an edge (a multigraph); use [aggregate_network()] for one weighted edge
#' per pair of nodes, or [network_window()] and [snapshots()] for the
#' contacts of a period.
#'
#' **Log data.** Tables of events such as authentication logs or network
#' flows (one row per event, with source, destination, and time) are
#' contact sequences: pass them, or a netframe built from them, straight to
#' `temporal_network()`. Millions of contacts are fine; slicing
#' and paths run in compiled code after a single sort by time.
#'
#' @section Missing values and invalid input:
#' Start times must not be missing and, with `end` or `duration`, neither
#' may end times or durations; ends must not precede starts. Errors name the
#' rows of `edges`. Start and end columns must have the same type. Everything
#' else is checked as in [netframe()].
#'
#' @references
#' Holme, P., & Saramäki, J. (2012). Temporal networks. *Physics Reports*,
#' 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
#'
#' Butts, C. T., Leslie-Cook, A., Krivitsky, P. N., & Bender-deMoll, S.
#' (2026). networkDynamic: Dynamic extensions for network objects. R
#' package version 0.12.0. \doi{10.32614/CRAN.package.networkDynamic}
#'
#' @seealso [network_at()], [network_window()], [snapshots()],
#'   [aggregate_network()], [temporal_distances()], [temporal_reach()],
#'   [temporal_closeness()], [rolling()], [burstiness()].
#'
#' @return `temporal_network()` returns a [netframe()] with time roles;
#'   `is_temporal()` returns `TRUE` or `FALSE`.
#' @export
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

#' @rdname temporal_network
#' @export
is_temporal <- function(x) inherits(x, "netframe") && !is.null(x$roles$start)

#' A temporal network at a time or in a window
#'
#' `network_at()` keeps the contacts active at one time and
#' `network_window()` those active at some time in the window
#' \eqn{[from, to)}, following the half-open rule of [temporal_network()].
#' Both return a temporal network with all columns and roles, contacts in
#' their original order, and contacts unchanged (not clipped to the
#' window).
#'
#' @param x A temporal network ([temporal_network()]).
#' @param time A time, of the type of the network's times.
#' @param from,to Start (included) and end (excluded) of the window.
#' @param nodes Which nodes to keep: `"all"` (the default) keeps every node,
#'   so results of measures stay aligned across times; `"present"` keeps the
#'   nodes present according to the network's node presence intervals (every
#'   node if it has none); `"active"` keeps only nodes with a contact in the
#'   result. Endpoints of the contacts kept are always kept.
#'
#' @return A temporal network ([netframe()]).
#' @examples
#' calls <- data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d"),
#'                     time = c(1, 2, 2, 5), end = c(3, 2, 4, 6))
#' g <- temporal_network(calls, end = "end", directed = FALSE)
#' network_at(g, 2)$edges
#' network_window(g, 3, 5, nodes = "active")$edges
#' @export
network_at <- function(x, time, nodes = c("all", "present", "active")) {
  check_temporal(x)
  nodes <- match.arg(nodes)
  t <- time_arg(time, x, "time")
  ptr <- engine_temporal(x)
  idx <- engine_call("network_at", .Call(C_temporal_active_at, ptr, t)) + 1
  present <- if (nodes == "present") node_present(x, t, t) else NULL
  slice_netframe(x, idx, nodes, present)
}

#' @rdname network_at
#' @export
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

#' Snapshots of a temporal network
#'
#' Splits a temporal network into a series of windows, sliding (`step`
#' smaller than `width`), tumbling (`step` equal to `width`, the default),
#' or with gaps (`step` larger), or between given `breaks`. Each window
#' \eqn{[a, b)} holds the contacts active in it, as in [network_window()].
#'
#' Contacts are sorted by time once and the whole series is computed in one
#' pass, with a heap of active contacts, in \eqn{O(m \log m)}{O(m log m)}
#' time plus the size of the output, so series over millions of contacts
#' are fast. An interval contact appears in every window it overlaps.
#'
#' @inheritParams network_at
#' @param width Window length: a number in the unit of the times (seconds
#'   for date-times, days for dates) or a `difftime`.
#' @param step Distance between window starts, as `width`.
#' @param breaks Alternatively to `width`, increasing times: windows are
#'   \eqn{[breaks_k, breaks_{k+1})}.
#' @param from,to With `width`: the first window starts at `from` (default:
#'   the earliest start time), and windows start up to and including `to`
#'   (default: the latest end time), so with `step <= width` every contact
#'   falls in some window.
#'
#' @return A list of temporal networks, one per window, with attribute
#'   `"windows"`, a data frame of window `start` and `end` times.
#' @examples
#' logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2"),
#'                      host = c("h1", "h1", "h2", "h2", "h3"),
#'                      time = c(0, 1, 5, 6, 12))
#' g <- temporal_network(logins)
#' s <- snapshots(g, width = 5)
#' attr(s, "windows")
#' vapply(s, function(w) nrow(w$edges), numeric(1))
#' # Sliding windows of width 6 every 3 time units.
#' vapply(snapshots(g, width = 6, step = 3), function(w) nrow(w$edges), numeric(1))
#' @export
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

#' Aggregate a temporal network into a weighted network
#'
#' Collapses the contacts between each pair of nodes (ordered pairs in
#' directed networks, unordered otherwise) into one edge carrying the
#' number of contacts, their total duration, and the first and last times.
#' With `from` and `to`, only contacts active in \eqn{[from, to)} count and
#' durations are clipped to the window, so an instantaneous contact adds to
#' the count but not to the duration.
#'
#' Self-loops are kept as edges from a node to itself. Edges appear in the
#' order of each pair's first contact; in undirected networks the endpoints
#' are given in the order of that first contact. The original weights, if
#' any, are not used. All nodes are kept, with their columns.
#'
#' @inheritParams network_at
#' @param weight Which column becomes the weight: `"count"` or
#'   `"duration"`. Either is a tie strength (`weight_type = "strength"`), so
#'   path-based measures use 1/weight as the length.
#' @param from,to Optional window; both or neither.
#'
#' @return A static [netframe()] whose edges have columns for the
#'   endpoints, `count`, `duration` (in the unit of the times), `first`
#'   (earliest start), and `last` (latest end).
#' @examples
#' calls <- data.frame(from = c("a", "b", "a", "c"), to = c("b", "a", "c", "a"),
#'                     time = c(0, 10, 20, 25), end = c(5, 12, 20, 40))
#' g <- temporal_network(calls, end = "end", directed = FALSE)
#' aggregate_network(g)$edges
#' aggregate_network(g, weight = "duration", from = 0, to = 30)$edges
#' @export
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

#' Time-respecting paths
#'
#' Optimal time-respecting walks between nodes of a temporal network:
#' earliest arrival, latest departure, fastest, and shortest (fewest
#' contacts), or Pan and Saramäki's (2011) average temporal distance.
#'
#' @inheritParams network_at
#' @param type `"earliest_arrival"`, `"latest_departure"`, `"fastest"`,
#'   `"shortest"`, or `"average"`; see Details.
#' @param from,to Identifiers of the origins (rows) and destinations
#'   (columns); `NULL` uses every node.
#' @param start,end The window: walks leave at or after `start` and arrive
#'   at or before `end`. `NULL` uses the earliest start time and the latest
#'   arrival of any contact.
#' @param min_wait,max_wait Waiting times allowed at each intermediate node,
#'   from arriving to leaving again (a number in the unit of the times, or a
#'   `difftime`). `max_wait = Inf` sets no limit.
#' @param strict If `TRUE`, waiting times must be strictly greater than
#'   `min_wait`; with `min_wait = 0`, consecutive contacts must have strictly
#'   increasing times, so contacts at the same time never chain.
#' @param traversal Time to traverse a contact: a number (or `difftime`),
#'   `0` by default, or `"duration"` for each contact's own duration (its end
#'   minus its start), as for a flight or a phone call that must finish
#'   before the information moves on.
#' @param threads Number of threads; `0` uses all available cores (see
#'   "Threads" in [netframe()]).
#'
#' @details
#' **Walks.** Following Wu et al. (2014), each contact from \eqn{u} to
#' \eqn{v} starting at \eqn{t} is a *temporal arc* that leaves \eqn{u} at
#' \eqn{t} and arrives at \eqn{v} at \eqn{t + \lambda}, where \eqn{\lambda}
#' is the traversal time. A time-respecting walk is a sequence of arcs, each
#' leaving the node where the previous one arrived, such that the waiting
#' time \eqn{w} between arriving and leaving satisfies
#' \eqn{w \ge} `min_wait` (or \eqn{w >} `min_wait` with `strict = TRUE`) and
#' \eqn{w \le} `max_wait`. The default, non-decreasing times
#' (`min_wait = 0`, `strict = FALSE`), is the definition of Kempe,
#' Kleinberg, and Kumar (2002) and of Holme and Saramäki (2012), under which
#' contacts at the same time can chain (\eqn{a \to b}{a -> b} and
#' \eqn{b \to c}{b -> c} at the same time let \eqn{a} reach \eqn{c}). With
#' `strict = TRUE` times strictly increase, as in Pan and Saramäki (2011),
#' who also require each contact to start after the previous one ends
#' (`strict = TRUE, traversal = "duration"`). The walk may leave its origin
#' at any time in the window; waiting rules apply only at intermediate
#' nodes. Undirected contacts can be traversed either way; self-loops never
#' help. Interval contacts are traversed from their start time; to let walks
#' join a long contact part-way, split it into shorter contacts. Node
#' presence intervals and weights are not used.
#'
#' **Values** for walks from \eqn{i} to \eqn{j} (row \eqn{i}, column
#' \eqn{j}):
#' \describe{
#'   \item{`"earliest_arrival"`}{the earliest arrival time at \eqn{j}
#'     (Wu et al. 2014); `start` (the window start) on the diagonal.}
#'   \item{`"latest_departure"`}{the latest time to leave \eqn{i} and still
#'     reach \eqn{j} by `end`; `end` on the diagonal, `-Inf` without a
#'     walk.}
#'   \item{`"fastest"`}{the smallest elapsed time from leaving \eqn{i} to
#'     arriving at \eqn{j}.}
#'   \item{`"shortest"`}{the fewest contacts.}
#'   \item{`"average"`}{Pan and Saramäki's (2011, Eq. 1) average temporal
#'     distance: the time from a starting time \eqn{t} to the earliest
#'     arrival at \eqn{j}, averaged over all \eqn{t} in the window, where a
#'     start after the last walk waits for the first walk of the window
#'     repeated after its end (their periodic boundary condition). It needs a
#'     window of positive length (otherwise `NaN`).}
#' }
#' Unreachable pairs get `Inf` (`-Inf` for latest departure). Times are
#' numbers on the scale of the data: as given for numeric times, seconds
#' since 1970-01-01 UTC for date-times, and days since then for dates.
#' Lengths of time are in the unit of the times. Without `max_wait`, an
#' optimal walk never needs to repeat a node, so the values are those of
#' time-respecting paths; with `max_wait`, waiting can only be cut short by
#' moving on, and the optimum is taken over walks.
#'
#' **Algorithms.** Contacts are sorted once by time into an edge stream and
#' each origin costs one pass over it (Wu et al. 2014): latest departure is
#' earliest arrival in reversed time, fastest walks keep the latest
#' departure that reaches each arrival, and fewest contacts the smallest
#' count. Waiting limits use a heap of pending arrivals and a monotone queue
#' of usable ones. Earliest arrival without `max_wait` runs 64 origins per
#' pass. Origins run in parallel. Latest departure costs one pass per
#' destination, and `"average"` one pass per node of the network. Results
#' are checked against a brute-force search over all walks.
#'
#' @references
#' Kempe, D., Kleinberg, J., & Kumar, A. (2002). Connectivity and inference
#' problems for temporal networks. *Journal of Computer and System
#' Sciences*, 64(4), 820-842. \doi{10.1006/jcss.2002.1829}
#'
#' Holme, P., & Saramäki, J. (2012). Temporal networks. *Physics Reports*,
#' 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
#'
#' Pan, R. K., & Saramäki, J. (2011). Path lengths, correlations, and
#' centrality in temporal networks. *Physical Review E*, 84(1), 016105.
#' \doi{10.1103/PhysRevE.84.016105}
#'
#' Wu, H., Cheng, J., Huang, S., Ke, Y., Lu, Y., & Xu, Y. (2014). Path
#' problems in temporal graphs. *Proceedings of the VLDB Endowment*, 7(9),
#' 721-732. \doi{10.14778/2732939.2732945}
#'
#' @return A numeric matrix with one row per origin and one column per
#'   destination, named by node identifiers.
#' @examples
#' # b passes a message from a to c only if it hears it first.
#' msgs <- data.frame(from = c("a", "b", "b", "a"), to = c("b", "c", "c", "c"),
#'                    time = c(1, 1, 3, 6))
#' g <- temporal_network(msgs)
#' temporal_distances(g)
#' # Non-decreasing times (Kempe et al. 2002; Holme and Saramaki 2012): the
#' # contacts a -> b and b -> c at time 1 chain, so a reaches c at time 1.
#' temporal_distances(g, from = "a", strict = FALSE)
#' # Strictly increasing times (Pan and Saramaki 2011): a reaches c at 3.
#' temporal_distances(g, from = "a", strict = TRUE)
#' temporal_distances(g, "fastest")
#' temporal_distances(g, "shortest", max_wait = 1)
#' @export
temporal_distances <- function(x, type = c("earliest_arrival", "latest_departure", "fastest",
                                           "shortest", "average"),
                               from = NULL, to = NULL, start = NULL, end = NULL,
                               min_wait = 0, max_wait = Inf, strict = FALSE, traversal = 0,
                               threads = 0L) {
  check_temporal(x)
  type <- match.arg(type)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  n <- nrow(x$nodes)
  rows <- if (is.null(from)) seq_len(n) else temporal_positions(x, from, "from")
  cols <- if (is.null(to)) seq_len(n) else temporal_positions(x, to, "to")
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

#' Temporal reachability
#'
#' `temporal_reach()` counts, for every node, the other nodes it can reach
#' by time-respecting walks (`mode = "out"`: the size of its *set of
#' influence*) or that can reach it (`mode = "in"`: its *source set*), after
#' Holme (2005). `reachability_ratio()` is the share of ordered pairs of
#' distinct nodes joined by a time-respecting walk (Holme 2005), which is 1
#' in a network where everyone can eventually reach everyone in time order.
#'
#' Walks follow the rules described in [temporal_distances()], including
#' the window and waiting times. Each node costs one pass over the
#' contacts (64 nodes per pass without `max_wait`), in parallel.
#'
#' @inheritParams temporal_distances
#' @param mode `"out"` or `"in"`.
#'
#' @references
#' Holme, P. (2005). Network reachability of real-world contact sequences.
#' *Physical Review E*, 71(4), 046119. \doi{10.1103/PhysRevE.71.046119}
#'
#' Holme, P., & Saramäki, J. (2012). Temporal networks. *Physics Reports*,
#' 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
#'
#' @return `temporal_reach()` returns a numeric vector aligned with the rows
#'   of `x$nodes`; `reachability_ratio()` a number (`NaN` with fewer than two
#'   nodes).
#' @examples
#' msgs <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"), time = c(1, 2, 1))
#' g <- temporal_network(msgs)
#' temporal_reach(g)
#' temporal_reach(g, mode = "in")
#' reachability_ratio(g)
#' # Simultaneous contacts: a -> b and b -> c both at time 1.
#' same <- temporal_network(data.frame(from = c("a", "b"), to = c("b", "c"), time = c(1, 1)))
#' # Non-decreasing times (Kempe et al. 2002; Holme and Saramaki 2012): a reaches c.
#' temporal_reach(same, strict = FALSE)
#' # Strictly increasing times (Pan and Saramaki 2011): a reaches only b.
#' temporal_reach(same, strict = TRUE)
#' @export
temporal_reach <- function(x, mode = c("out", "in"), start = NULL, end = NULL, min_wait = 0,
                           max_wait = Inf, strict = FALSE, traversal = 0, threads = 0L) {
  check_temporal(x)
  mode <- match.arg(mode)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  engine_call("temporal_reach", .Call(C_temporal_reach, engine_temporal(x),
                                      if (mode == "out") 1L else 2L, opts))
}

#' @rdname temporal_reach
#' @export
reachability_ratio <- function(x, start = NULL, end = NULL, min_wait = 0, max_wait = Inf,
                               strict = FALSE, traversal = 0, threads = 0L) {
  check_temporal(x)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  n <- nrow(x$nodes)
  if (n < 2L) return(NaN)
  r <- engine_call("reachability_ratio", .Call(C_temporal_reach, engine_temporal(x), 1L, opts))
  sum(r) / (n * (n - 1))
}

#' Temporal closeness
#'
#' Pan and Saramäki's (2011) temporal closeness centrality
#' \deqn{C_i = \frac{1}{n - 1} \sum_{j \ne i} \frac{1}{\tau_{ij}}}{C_i = (1 / (n - 1)) sum_{j != i} 1 / tau_ij}
#' where \eqn{\tau_{ij}} is the average temporal distance from \eqn{i} to
#' \eqn{j} over all starting times in the window (their Eq. 1; see
#' `type = "average"` in [temporal_distances()]). Nodes that cannot reach
#' \eqn{j} at all contribute 0, as in harmonic centrality. High values mean
#' that information leaving the node at a random time reaches the others
#' quickly.
#'
#' The average temporal distance is computed exactly from each pair's
#' fastest-path profile (the walks that leave latest for each arrival), with
#' one pass over the contacts per node, in parallel. Values are in reciprocal
#' units of time (per second for date-times, per day for dates).
#'
#' @inheritParams temporal_distances
#' @param strict,traversal Path rules, as in [temporal_distances()]. Pan and
#'   Saramäki's own rule is `strict = TRUE, traversal = "duration"`: each
#'   contact starts strictly after the previous one ends. The defaults allow
#'   simultaneous contacts to chain, as in Holme and Saramäki (2012).
#'
#' @references
#' Pan, R. K., & Saramäki, J. (2011). Path lengths, correlations, and
#' centrality in temporal networks. *Physical Review E*, 84(1), 016105.
#' \doi{10.1103/PhysRevE.84.016105}
#'
#' Holme, P., & Saramäki, J. (2012). Temporal networks. *Physics Reports*,
#' 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`; `NaN` for
#'   every node when there are fewer than two nodes or the window has no
#'   length.
#' @examples
#' msgs <- data.frame(from = c("a", "b", "b", "a"), to = c("b", "c", "c", "c"),
#'                    time = c(1, 1, 3, 6))
#' g <- temporal_network(msgs)
#' # Simultaneous contacts chain (Holme and Saramaki 2012).
#' temporal_closeness(g, start = 0, end = 10)
#' # Strictly increasing times (Pan and Saramaki 2011).
#' temporal_closeness(g, start = 0, end = 10, strict = TRUE)
#' @export
temporal_closeness <- function(x, start = NULL, end = NULL, min_wait = 0, max_wait = Inf,
                               strict = FALSE, traversal = 0, threads = 0L) {
  check_temporal(x)
  opts <- path_options(x, start, end, min_wait, max_wait, strict, traversal, threads)
  r <- engine_call("temporal_closeness",
                   .Call(C_temporal_closeness, engine_temporal(x), opts, FALSE))
  r[[1L]] * time_scale(x)
}

#' Apply a measure over a series of windows
#'
#' Computes a node-level (or network-level) measure on every snapshot of a
#' temporal network ([snapshots()]) and returns the results as one tidy
#' data frame, ready for plotting or modelling change over time.
#'
#' Snapshots keep every node by default, so node-level results line up
#' across windows. `fun` receives each snapshot (a temporal network holding
#' the contacts active in the window) and must return either one value per
#' node of the snapshot or a single value; any knappnet measure works, for
#' example [degree()], [betweenness()], or [graph_density()]. Static
#' measures see each contact as an edge; aggregate inside `fun` (as in the
#' examples) to count each pair of nodes once.
#'
#' @inheritParams snapshots
#' @param fun A function, or the name of one, taking a netframe as its
#'   first argument.
#' @param ... Further arguments to `fun`.
#'
#' @return A data frame with columns `start` and `end` (of each window),
#'   `node` (the node identifier, or `NA` for single values), and `value`.
#' @examples
#' logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3"),
#'                      host = c("h1", "h1", "h2", "h2", "h3", "h1"),
#'                      time = c(0, 1, 5, 6, 12, 13))
#' g <- temporal_network(logins)
#' rolling(g, degree, width = 5)
#' rolling(g, function(s) graph_density(aggregate_network(s)), width = 10, step = 5)
#' @export
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

#' Inter-event times and burstiness
#'
#' `inter_event_times()` lists the times between consecutive events of each
#' node or pair of nodes; `burstiness()` summarizes them with Goh and
#' Barabási's (2008) burstiness and memory coefficients.
#'
#' The events of a node are the start times of its contacts: all of them
#' (`mode = "all"`), or in directed networks those it sends (`"out"`) or
#' receives (`"in"`); a self-loop is one event. The events of a pair are the
#' start times of the contacts between the two nodes (ordered pairs in
#' directed networks). Simultaneous events give inter-event times of 0.
#'
#' With inter-event times \eqn{\tau_1, \dots, \tau_k} of mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}, the burstiness is
#' \deqn{B = \frac{\sigma - \mu}{\sigma + \mu},}{B = (sigma - mu) / (sigma + mu),}
#' which is -1 for perfectly regular events, about 0 for a Poisson process,
#' and approaches 1 for very bursty activity (Goh and Barabási 2008, Eq. 1).
#' The memory coefficient is the correlation of consecutive inter-event
#' times (their Eq. 4),
#' \deqn{M = \frac{1}{k - 1} \sum_{i=1}^{k-1} \frac{(\tau_i - \mu_1)(\tau_{i+1} - \mu_2)}{\sigma_1 \sigma_2},}{M = (1 / (k - 1)) sum_i (tau_i - mu_1)(tau_{i+1} - mu_2) / (sigma_1 sigma_2),}
#' where \eqn{\mu_1, \sigma_1} (\eqn{\mu_2, \sigma_2}) are the mean and
#' standard deviation of \eqn{\tau_1, \dots, \tau_{k-1}}
#' (\eqn{\tau_2, \dots, \tau_k}). Standard deviations divide by the number
#' of values (the standard deviation of the observed distribution), so
#' \eqn{M} is exactly the Pearson correlation; `sd()` divides by one less,
#' which gives a slightly different \eqn{B}. \eqn{B} is `NaN` with fewer than
#' two inter-event times or when all are 0, and \eqn{M} with fewer than
#' three or when either sequence is constant. Neither coefficient corrects
#' for the number of events (Kim and Jo 2016 discuss finite-size effects),
#' so compare nodes with similar activity.
#'
#' @inheritParams network_at
#' @param by `"node"` or `"pair"`.
#' @param mode For nodes of directed networks: `"all"`, `"out"`, or `"in"`.
#'
#' @references
#' Goh, K.-I., & Barabási, A.-L. (2008). Burstiness and memory in complex
#' systems. *EPL (Europhysics Letters)*, 81(4), 48002.
#' \doi{10.1209/0295-5075/81/48002}
#'
#' Kim, E.-K., & Jo, H.-H. (2016). Measuring burstiness for finite event
#' sequences. *Physical Review E*, 94(3), 032311.
#' \doi{10.1103/PhysRevE.94.032311}
#'
#' @return `inter_event_times()` returns a data frame with the node (or the
#'   pair, as `from` and `to`) and `gap`, the inter-event time in the unit
#'   of the times, ordered by node (or pair) and then time.
#'   `burstiness()` returns a data frame with one row per node (aligned with
#'   `x$nodes`) or per pair (in order of first contact), with `events`,
#'   `burstiness`, and `memory`.
#' @examples
#' ev <- data.frame(from = c("a", "a", "a", "a", "b", "b", "b"),
#'                  to = c("b", "b", "c", "c", "c", "c", "c"),
#'                  time = c(0, 1, 2, 10, 0, 5, 10))
#' g <- temporal_network(ev, directed = FALSE)
#' inter_event_times(g)
#' burstiness(g)
#' burstiness(g, by = "pair")
#' @export
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

#' @rdname burstiness
#' @export
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

# Positions (1-based) of node identifiers in x$nodes.
temporal_positions <- function(x, nodes, what) {
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (is.factor(nodes)) nodes <- as.character(nodes)
  pos <- match(nodes, ids)
  if (anyNA(pos) || !length(pos)) {
    if (!length(pos)) stop(sprintf("`%s` must name at least one node.", what), call. = FALSE)
    shown <- utils::head(unique(nodes[is.na(pos)]), 5L)
    stop(sprintf("%d node(s) in `%s` are not in the network, for example: %s.",
                 length(unique(nodes[is.na(pos)])), what,
                 paste(format(shown), collapse = ", ")), call. = FALSE)
  }
  pos
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
