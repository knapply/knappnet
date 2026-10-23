# Plot a QAP correlation test (man/plot.qap_cor.Rd).
plot.qap_cor <- function(x, breaks = "Sturges", ...) {
  d <- x$dist[is.finite(x$dist)]
  if (!length(d)) {
    stop("There are no permutation draws to plot; run qap_cor() with `reps > 0`.",
         call. = FALSE)
  }
  dots <- list(...)
  args <- c(list(x = d, breaks = breaks, plot = FALSE), dots[intersect(names(dots), "include.lowest")])
  h <- do.call(graphics::hist, args)
  observed <- x$statistic
  xlim <- range(c(h$breaks, observed), na.rm = TRUE)
  if (is.null(dots$main)) dots$main <- "QAP correlation test"
  if (is.null(dots$xlab)) dots$xlab <- "Correlation under permutation"
  if (is.null(dots$col)) dots$col <- "gray80"
  if (is.null(dots$border)) dots$border <- "white"
  if (is.null(dots$xlim)) dots$xlim <- xlim
  old <- graphics::par(mar = c(5.1, 4.1, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  do.call(graphics::hist, c(list(x = d, breaks = breaks), dots))
  if (is.finite(observed)) {
    graphics::abline(v = observed, col = plot_accent, lwd = 2.5)
    usr <- graphics::par("usr")
    right_side <- observed > mean(usr[1:2])
    graphics::text(observed, usr[4] - 0.04 * diff(usr[3:4]),
                   sprintf("observed r = %s", format(observed, digits = 3)),
                   adj = c(if (right_side) 1.05 else -0.05, 1), col = plot_accent, cex = 0.85)
  }
  undefined <- length(x$dist) - length(d)
  sub <- sprintf("Pr(>=|r|) = %s over %s; %s dyads", format(x$p.value[[3L]], digits = 3),
                 qap_draws_text(x), format(as.integer(x$dyads), big.mark = ","))
  if (undefined > 0L) sub <- sprintf("%s (%d draws undefined)", sub, undefined)
  graphics::mtext(sub, side = 3, line = 0.3, cex = 0.8, col = "gray30")
  invisible(x)
}

# Plot a multiple regression QAP (man/plot.mrqap.Rd).
plot.mrqap <- function(x, which = c("coefficients", "distributions"), level = 0.95,
                       intercept = FALSE, breaks = "Sturges", ...) {
  which <- match.arg(which)
  intercept <- check_flag(intercept, "intercept")
  level <- check_number(level, "level", lower = 0, upper = 1, open = TRUE)
  observed <- if (x$statistic == "t") x$t else x$coefficients
  what <- if (x$statistic == "t") "t statistic" else "coefficient"
  keep <- rep(TRUE, length(observed))
  if (!intercept && isTRUE(x$intercept) && length(observed) > 1L) keep[1L] <- FALSE
  if (!any(keep)) stop("There are no coefficients to plot.", call. = FALSE)
  names_ <- names(observed)[keep]
  observed <- observed[keep]
  dist <- x$dist[, keep, drop = FALSE]
  p <- x$p.value[keep, 3L]
  if (nrow(dist) == 0L) {
    stop("There are no permutation draws to plot; run mrqap() with `reps > 0`.", call. = FALSE)
  }
  dots <- list(...)
  k <- length(observed)
  if (which == "coefficients") {
    alpha <- (1 - level) / 2
    lo <- apply(dist, 2L, stats::quantile, probs = alpha, na.rm = TRUE, names = FALSE)
    hi <- apply(dist, 2L, stats::quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)
    xlim <- range(c(lo, hi, observed, 0), na.rm = TRUE)
    xlim <- xlim + c(-0.08, 0.3) * diff(xlim)
    left <- label_lines(names_, 0.9) + 1.5
    old <- graphics::par(mar = c(5.1, max(left, 4.1), 4.1, 1.1))
    on.exit(graphics::par(old), add = TRUE)
    if (is.null(dots$main)) dots$main <- "Multiple regression QAP"
    if (is.null(dots$xlab)) dots$xlab <- sprintf("%s (observed point, permutation interval)", what)
    if (is.null(dots$xlim)) dots$xlim <- xlim
    y <- rev(seq_len(k))
    do.call(graphics::plot, c(list(x = observed, y = y, type = "n", ylim = c(0.5, k + 0.5),
                                   yaxt = "n", ylab = ""), dots))
    graphics::abline(v = 0, col = "gray75", lty = 3)
    graphics::segments(lo, y, hi, y, col = plot_gray, lwd = 3, lend = 1)
    graphics::points(observed, y, pch = 21, bg = plot_accent, col = "white", cex = 1.6)
    graphics::axis(2, at = y, labels = names_, las = 1, tick = FALSE, cex.axis = 0.9)
    right <- pmax(hi, observed, na.rm = TRUE)
    graphics::text(right, y, sprintf("p = %s", format(p, digits = 2)), pos = 4, cex = 0.75,
                   col = "gray30")
    graphics::mtext(sprintf("%s; %d%% permutation intervals", qap_method_text(x),
                            round(100 * level)), side = 3, line = 0.3, cex = 0.75, col = "gray30")
  } else {
    rows <- ceiling(sqrt(k))
    cols <- ceiling(k / rows)
    old <- graphics::par(mfrow = c(rows, cols), mar = c(4.1, 4.1, 3.1, 1.1))
    on.exit(graphics::par(old), add = TRUE)
    for (j in seq_len(k)) {
      d <- dist[, j]
      d <- d[is.finite(d)]
      if (!length(d)) {
        graphics::plot.new()
        graphics::title(main = names_[j])
        next
      }
      h <- graphics::hist(d, breaks = breaks, plot = FALSE)
      args <- dots
      args$main <- names_[j]
      if (is.null(args$xlab)) args$xlab <- what
      if (is.null(args$col)) args$col <- "gray80"
      if (is.null(args$border)) args$border <- "white"
      if (is.null(args$xlim)) args$xlim <- range(c(h$breaks, observed[j]), na.rm = TRUE)
      do.call(graphics::hist, c(list(x = d, breaks = breaks), args))
      graphics::abline(v = observed[j], col = plot_accent, lwd = 2.5)
      graphics::mtext(sprintf("Pr(>=|b|) = %s", format(p[j], digits = 2)), side = 3,
                      line = 0.2, cex = 0.75, col = "gray30")
    }
  }
  invisible(x)
}
