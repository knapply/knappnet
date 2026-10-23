# Quadratic assignment procedure tests: QAP correlation and multiple
# regression QAP, with their print, summary, and plot methods.

# QAP correlation test (man/net_qap_cor.Rd).
net_qap_cor <- function(y, x, reps = 1000, seed = NULL, diag = FALSE,
                    mode = c("digraph", "graph"), exact = FALSE, attr = NULL,
                    threads = 0L) {
  mode <- match.arg(mode)
  diag <- check_flag(diag, "diag")
  exact <- check_flag(exact, "exact")
  threads <- resolve_threads(threads)
  m <- qap_matrices(list(y, x), attr)
  n <- m$n
  mats <- m$values
  draws <- qap_draws(n, reps, exact)
  seed <- check_seed(seed)
  res <- engine_call("net_qap_cor",
    .Call(C_qap_cor, mats[[1]], mats[[2]], n, qap_mode_code(mode), diag, exact,
          as.double(reps), seed, threads, as.double(draws)))
  s <- res[[1]]
  counts <- c(le = s[4], ge = s[5], abs = s[6])
  valid <- s[7]
  structure(list(
    statistic = s[1],
    p.value = c(`Pr(<=r)` = s[4] / valid, `Pr(>=r)` = s[5] / valid,
                `Pr(>=|r|)` = s[6] / valid),
    counts = counts,
    valid = valid,
    draws = s[3],
    dyads = s[2],
    dist = res[[2]],
    mode = mode,
    diag = isTRUE(diag),
    exact = isTRUE(exact),
    seed = seed,
    call = match.call()
  ), class = "qap_cor")
}

# See man/net_qap_cor.Rd.
print.qap_cor <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nQAP correlation test\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(sprintf("Correlation: %s over %d dyads (%s)\n",
              format(x$statistic, digits = digits), as.integer(x$dyads), x$mode))
  cat(qap_draws_text(x), "\n", sep = "")
  p <- x$p.value
  cat(sprintf("%s: %s\n", names(p), format(p, digits = digits)), sep = "")
  invisible(x)
}

# See man/net_qap_cor.Rd.
summary.qap_cor <- function(object, ...) {
  q <- stats::quantile(object$dist, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),
                       na.rm = TRUE, names = TRUE)
  structure(c(object, list(quantiles = q)), class = c("summary.qap_cor", "qap_cor"))
}

# See man/net_qap_cor.Rd.
print.summary.qap_cor <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.qap_cor(x, digits = digits)
  cat("\nReference distribution of r:\n")
  print(x$quantiles, digits = digits)
  invisible(x)
}

# Plot a QAP correlation test (man/plot.qap_cor.Rd).
plot.qap_cor <- function(x, breaks = "Sturges", ...) {
  d <- x$dist[is.finite(x$dist)]
  if (!length(d)) {
    stop("There are no permutation draws to plot; run net_qap_cor() with `reps > 0`.",
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

# Multiple regression QAP (man/net_mrqap.Rd).
net_mrqap <- function(y, x, reps = 1000, method = c("dsp", "y"), seed = NULL,
                  intercept = TRUE, mode = c("digraph", "graph"), diag = FALSE,
                  statistic = c("t", "beta"), exact = FALSE, attr = NULL,
                  tol = 1e-7, threads = 0L) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  statistic <- match.arg(statistic)
  intercept <- check_flag(intercept, "intercept")
  diag <- check_flag(diag, "diag")
  exact <- check_flag(exact, "exact")
  tol <- check_number(tol, "tol", lower = 0, upper = 1)
  threads <- resolve_threads(threads)
  if (is.matrix(x) || inherits(x, "netframe") || is.data.frame(x)) x <- list(x)
  if (!is.list(x)) stop("`x` must be a matrix, a netframe, or a list of them.", call. = FALSE)
  k <- length(x)
  xnames <- names(x)
  if (is.null(xnames)) xnames <- rep("", k)
  xnames[xnames == ""] <- paste0("x", seq_len(k))[xnames == ""]
  m <- qap_matrices(c(list(y), unname(x)), attr)
  n <- m$n
  mats <- m$values
  draws <- qap_draws(n, reps, exact)
  seed <- check_seed(seed)
  xs <- if (k) unlist(mats[-1], use.names = FALSE) else double(0)
  res <- engine_call("net_mrqap",
    .Call(C_mrqap, mats[[1]], as.double(xs), n, k, intercept, qap_mode_code(mode), diag,
          if (method == "dsp") 0L else 1L, if (statistic == "t") 0L else 1L,
          as.double(tol), exact, as.double(reps), seed, threads, as.double(draws)))
  coef_names <- c(if (isTRUE(intercept)) "(Intercept)", xnames)
  tab <- res[[1]]
  fit <- res[[2]]
  dist <- res[[3]]
  colnames(dist) <- coef_names
  valid <- tab[, 7]
  p <- cbind(`Pr(<=b)` = tab[, 4] / valid, `Pr(>=b)` = tab[, 5] / valid,
             `Pr(>=|b|)` = tab[, 6] / valid)
  rownames(p) <- coef_names
  counts <- cbind(le = tab[, 4], ge = tab[, 5], abs = tab[, 6], valid = valid)
  rownames(counts) <- coef_names
  structure(list(
    coefficients = stats::setNames(tab[, 1], coef_names),
    std.error = stats::setNames(tab[, 2], coef_names),
    t = stats::setNames(tab[, 3], coef_names),
    p.value = p,
    counts = counts,
    dyads = fit[1],
    df.residual = fit[2],
    rss = fit[3],
    mss = fit[4],
    r.squared = fit[5],
    adj.r.squared = fit[6],
    sigma = fit[7],
    fstatistic = c(value = fit[8], numdf = length(coef_names) - isTRUE(intercept),
                   dendf = fit[2]),
    method = if (fit[9] == 0) "dsp" else "y",
    statistic = statistic,
    draws = fit[10],
    dist = dist,
    mode = mode,
    diag = isTRUE(diag),
    intercept = isTRUE(intercept),
    exact = isTRUE(exact),
    seed = seed,
    call = match.call()
  ), class = "mrqap")
}

# See man/net_mrqap.Rd.
print.mrqap <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nMultiple regression QAP\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  tab <- cbind(Estimate = x$coefficients, `Pr(>=|b|)` = x$p.value[, 3])
  print(tab, digits = digits)
  cat("\n", qap_method_text(x), "\n", sep = "")
  invisible(x)
}

# See man/net_mrqap.Rd.
summary.mrqap <- function(object, ...) {
  tab <- cbind(Estimate = object$coefficients, `Std. Error` = object$std.error,
               `t value` = object$t, object$p.value)
  structure(c(object, list(table = tab)), class = c("summary.mrqap", "mrqap"))
}

# See man/net_mrqap.Rd.
print.summary.mrqap <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nMultiple regression QAP\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print(x$table, digits = digits)
  cat("\n", qap_method_text(x), "\n", sep = "")
  cat(sprintf("Residual standard error: %s on %d degrees of freedom (%d dyads)\n",
              format(signif(x$sigma, digits)), as.integer(x$df.residual),
              as.integer(x$dyads)))
  cat(sprintf("Multiple R-squared: %s,\tAdjusted R-squared: %s\n",
              formatC(x$r.squared, digits = digits),
              formatC(x$adj.r.squared, digits = digits)))
  if (is.finite(x$fstatistic[["value"]])) {
    cat(sprintf("F-statistic: %s on %d and %d DF (no classical p-value: dyads are not independent)\n",
                formatC(x$fstatistic[["value"]], digits = digits),
                as.integer(x$fstatistic[["numdf"]]), as.integer(x$fstatistic[["dendf"]])))
  }
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
    stop("There are no permutation draws to plot; run net_mrqap() with `reps > 0`.", call. = FALSE)
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

qap_mode_code <- function(mode) if (mode == "graph") 1L else 0L

qap_method_text <- function(x) {
  how <- if (x$method == "dsp") "double semi-partialling" else "Y-permutation"
  stat <- if (x$statistic == "t") "t-statistics" else "coefficients"
  paste0("p-values: ", how, " on ", stat, "; ", qap_draws_text(x))
}

qap_draws_text <- function(x) {
  if (isTRUE(x$exact)) {
    sprintf("all %s permutations (exact)", format(x$draws, big.mark = ","))
  } else {
    sprintf("%s random permutations (seed %s)", format(x$draws, big.mark = ","),
            format(x$seed, scientific = FALSE))
  }
}

# Number of draws whose statistics are returned.
qap_draws <- function(n, reps, exact) {
  if (isTRUE(exact)) {
    if (n > 10) stop("`exact = TRUE` needs at most 10 nodes.", call. = FALSE)
    return(factorial(n))
  }
  if (!is.numeric(reps) || length(reps) != 1L || is.na(reps) || reps < 0 ||
      reps != floor(reps) || reps > .Machine$integer.max) {
    stop("`reps` must be a whole number between 0 and .Machine$integer.max.", call. = FALSE)
  }
  reps
}

# Converts inputs to square matrices of one size, aligns node order by
# names when every input has them, and returns them as double vectors.
qap_matrices <- function(inputs, attr) {
  mats <- lapply(seq_along(inputs), function(i) {
    m <- inputs[[i]]
    if (inherits(m, "netframe")) {
      a <- if (is.null(attr)) m$roles$weight else attr
      m <- as.matrix(m, attr = a)
    } else if (is.data.frame(m)) {
      m <- as.matrix(m)
    }
    if (!is.matrix(m) || nrow(m) != ncol(m)) {
      stop(sprintf("Input %d is not a square matrix or a netframe.", i), call. = FALSE)
    }
    if (!(is.numeric(m) || is.logical(m))) {
      stop(sprintf("Input %d is not numeric.", i), call. = FALSE)
    }
    m
  })
  n <- nrow(mats[[1]])
  if (any(vapply(mats, nrow, 1L) != n)) {
    stop("All matrices must have the same number of nodes.", call. = FALSE)
  }
  ids <- rownames(mats[[1]])
  named <- vapply(mats, function(m) !is.null(rownames(m)) && !is.null(colnames(m)), TRUE)
  if (all(named) && !anyDuplicated(ids)) {
    mats <- lapply(mats, function(m) {
      if (!setequal(rownames(m), ids) || !setequal(colnames(m), ids)) {
        stop("Matrices have different node names.", call. = FALSE)
      }
      m[ids, ids, drop = FALSE]
    })
  }
  vecs <- lapply(mats, as.double)
  if (any(vapply(vecs, function(v) any(is.infinite(v)), TRUE))) {
    stop("Matrices must not contain infinite values.", call. = FALSE)
  }
  list(n = n, values = vecs)
}
