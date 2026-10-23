#' Timeline of contacts in a temporal network
#'
#' Draws when each pair of nodes (or each node) is in contact: one row per
#' pair, with a tick for every instantaneous contact and a bar for every
#' interval contact, the raw-data view of a contact sequence (Holme and
#' Saramäki 2012, Fig. 1).
#'
#' @param x A temporal network ([temporal_network()]).
#' @param by `"pair"` for one row per pair of nodes (ordered pairs in a
#'   directed network) or `"node"` for one row per node, on which every
#'   contact of the node appears.
#' @param max_rows How many rows to draw at most: the pairs or nodes with
#'   the most contacts are kept and the rest are counted in the subtitle.
#' @param ... Further arguments passed to [graphics::plot()], such as
#'   `main` or `xlab`.
#'
#' @details
#' Rows are sorted by number of contacts, the most active at the top. The
#' time axis is labelled in the type of the network's times (numbers,
#' dates, or date-times). Interval contacts are drawn from their start to
#' their end; an instantaneous contact is a tick at its time.
#'
#' @references
#' Holme, P., & Saramäki, J. (2012). Temporal networks. *Physics Reports*,
#' 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
#'
#' @return A data frame of the rows drawn, from the top, invisibly: `row`
#'   (the pair or node label) and `contacts`.
#' @seealso [temporal_network()], [plot_rolling()], [burstiness()].
#' @examples
#' calls <- data.frame(from = c("a", "a", "b", "c", "a", "b"),
#'                     to   = c("b", "c", "c", "d", "b", "c"),
#'                     time = c(1, 2, 2, 5, 7, 9), end = c(3, 2, 4, 6, 7, 10))
#' g <- temporal_network(calls, end = "end", directed = FALSE)
#' plot_timeline(g)
#' plot_timeline(g, by = "node", main = "Contacts per node")
#' @export
plot_timeline <- function(x, by = c("pair", "node"), max_rows = 50, ...) {
  check_temporal(x)
  by <- match.arg(by)
  if (!is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) || max_rows < 1) {
    stop("`max_rows` must be a positive number.", call. = FALSE)
  }
  tt <- temporal_times(x)
  m <- length(tt$start)
  if (m == 0L) stop("The network has no contacts to draw.", call. = FALSE)
  s <- tt$start
  e <- if (is.null(tt$end)) s else tt$end
  ids <- node_ids(x)
  src <- x$index$src + 1L
  tgt <- x$index$tgt + 1L
  if (by == "pair") {
    a <- if (x$directed) src else pmin(src, tgt)
    b <- if (x$directed) tgt else pmax(src, tgt)
    label <- paste(ids[a], ids[b], sep = if (x$directed) " -> " else " -- ")
    contact <- seq_len(m)
  } else {
    loop <- src == tgt
    label <- c(ids[src], ids[tgt[!loop]])
    contact <- c(seq_len(m), which(!loop))
  }
  counts <- table(label)
  rows <- names(counts)[order(-as.integer(counts), match(names(counts), label))]
  total <- length(rows)
  rows <- rows[seq_len(min(total, max_rows))]
  keep <- label %in% rows
  y <- match(label[keep], rows)
  y <- length(rows) + 1L - y
  cs <- s[contact[keep]]
  ce <- e[contact[keep]]
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- "Contact timeline"
  if (is.null(dots$xlab)) dots$xlab <- "time"
  left <- label_lines(rows, 0.8) + 1.5
  old <- graphics::par(mar = c(4.6, max(4.1, left), 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  rng <- range(c(cs, ce))
  if (!(diff(rng) > 0)) rng <- rng + c(-1, 1)
  do.call(graphics::plot, c(list(x = rng, y = c(0.5, length(rows) + 0.5), type = "n",
                                 axes = FALSE, ylab = ""), dots))
  graphics::abline(h = seq_along(rows), col = "gray92")
  time_axis(x, rng)
  graphics::axis(2, at = seq_along(rows), labels = rev(rows), las = 1, tick = FALSE,
                 cex.axis = 0.8, line = -0.5)
  graphics::box(col = "gray40")
  inst <- ce == cs
  if (any(!inst)) {
    graphics::segments(cs[!inst], y[!inst], ce[!inst], y[!inst], col = plot_blue, lwd = 4,
                       lend = 1)
  }
  if (any(inst)) {
    graphics::points(cs[inst], y[inst], pch = "|", col = plot_accent, cex = 1.2)
  }
  if (total > length(rows)) {
    graphics::mtext(sprintf("%d of %d %ss shown", length(rows), total, by), side = 3,
                    line = 0.3, cex = 0.8, col = "gray30")
  }
  invisible(data.frame(row = rows, contacts = as.integer(counts[rows]),
                       stringsAsFactors = FALSE))
}

# A time axis in the type of the network's times for a plot whose x axis
# holds seconds (or raw numbers).
time_axis <- function(x, rng) {
  at <- pretty(rng, 6)
  at <- at[at >= rng[1] & at <= rng[2]]
  if (length(at) < 2L) at <- rng
  labels <- format(from_seconds(at, x))
  graphics::axis(1, at = at, labels = labels, cex.axis = 0.8)
}

#' Plot a measure over time
#'
#' Draws the result of [rolling()], a measure computed on every window of a
#' temporal network, as lines over time: one line for a network-level
#' measure, or one per node for a node-level measure.
#'
#' @param x A data frame from [rolling()], with columns `start`, `end`,
#'   `node`, and `value`.
#' @param nodes For node-level results, the identifiers of the nodes to
#'   draw. `NULL` draws the `top` nodes by mean value.
#' @param top How many nodes to draw when `nodes` is `NULL`.
#' @param ... Further arguments passed to [graphics::plot()], such as
#'   `main` or `ylab`.
#'
#' @details
#' Each value is drawn at the start of its window. Windows where a node's
#' value is undefined (`NaN`) leave a gap in its line. The time axis uses
#' the type of the window times (numbers, dates, or date-times).
#'
#' @return `x`, invisibly.
#' @seealso [rolling()], [snapshots()], [plot_timeline()].
#' @examples
#' logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3", "u1"),
#'                      host = c("h1", "h1", "h2", "h2", "h3", "h1", "h1"),
#'                      time = c(0, 1, 5, 6, 12, 13, 14))
#' g <- temporal_network(logins)
#' plot_rolling(rolling(g, degree, width = 5, step = 2), main = "Degree over time")
#' plot_rolling(rolling(g, function(w) nrow(w$edges), width = 5, step = 2),
#'              ylab = "contacts in window")
#' @export
plot_rolling <- function(x, nodes = NULL, top = 8, ...) {
  if (!is.data.frame(x) || !all(c("start", "end", "node", "value") %in% names(x))) {
    stop("`x` must be the data frame returned by rolling().", call. = FALSE)
  }
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top < 1) {
    stop("`top` must be a positive number.", call. = FALSE)
  }
  if (nrow(x) == 0L) stop("There are no windows to draw.", call. = FALSE)
  dots <- list(...)
  starts <- sort(unique(x$start))
  network_level <- all(is.na(x$node))
  if (network_level) {
    series <- list(value = x$value[match(starts, x$start)])
  } else {
    ids <- as.character(x$node)
    if (is.null(nodes)) {
      means <- tapply(x$value, ids, function(v) if (all(is.na(v))) -Inf else mean(v, na.rm = TRUE))
      chosen <- names(sort(means, decreasing = TRUE))[seq_len(min(top, length(means)))]
    } else {
      chosen <- as.character(nodes)
      missing <- setdiff(chosen, ids)
      if (length(missing)) {
        stop(sprintf("`nodes` not in the results: %s.",
                     paste(utils::head(missing, 5L), collapse = ", ")), call. = FALSE)
      }
    }
    series <- lapply(chosen, function(id) {
      rows <- x[ids == id, , drop = FALSE]
      rows$value[match(starts, rows$start)]
    })
    names(series) <- chosen
  }
  values <- unlist(series, use.names = FALSE)
  finite <- values[is.finite(values)]
  ylim <- if (length(finite)) range(finite) else c(0, 1)
  if (!(diff(ylim) > 0)) ylim <- ylim + c(-0.5, 0.5)
  if (is.null(dots$main)) dots$main <- if (network_level) "Measure over time" else "Node measure over time"
  if (is.null(dots$xlab)) dots$xlab <- "window start"
  if (is.null(dots$ylab)) dots$ylab <- "value"
  if (is.null(dots$ylim)) dots$ylim <- ylim
  cols <- if (network_level) plot_blue else plot_qualitative(length(series))
  old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  do.call(graphics::plot, c(list(x = starts, y = rep(ylim[1], length(starts)), type = "n"), dots))
  graphics::grid(col = "gray90", lty = 1)
  for (i in seq_along(series)) {
    graphics::lines(starts, series[[i]], col = cols[i], lwd = 2)
    graphics::points(starts, series[[i]], col = cols[i], pch = 16, cex = 0.8)
  }
  if (!network_level) {
    graphics::legend("topright", legend = names(series), col = cols, lwd = 2, pch = 16,
                     bty = "n", cex = 0.8, title = "node", title.adj = 0)
  }
  invisible(x)
}
