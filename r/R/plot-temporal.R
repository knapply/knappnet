# Plots of temporal networks: contact timelines and measures over
# windows.

# Timeline of contacts in a temporal network (man/plot_timeline.Rd).
plot_timeline <- function(x, by = c("pair", "node"), max_rows = 50, ...) {
  check_temporal(x)
  by <- match.arg(by)
  max_rows <- check_number(max_rows, "max_rows", lower = 1)
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

# Plot a measure over time (man/plot_rolling.Rd).
plot_rolling <- function(x, nodes = NULL, top = 8, ...) {
  if (!is.data.frame(x) || !all(c("start", "end", "node", "value") %in% names(x))) {
    stop("`x` must be the data frame returned by rolling().", call. = FALSE)
  }
  top <- check_number(top, "top", lower = 1)
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
