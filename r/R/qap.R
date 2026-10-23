# QAP correlation test (man/qap_cor.Rd).
qap_cor <- function(y, x, reps = 1000, seed = NULL, diag = FALSE,
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
  seed <- qap_seed(seed)
  res <- engine_call("qap_cor",
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

# See man/qap_cor.Rd.
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

# See man/qap_cor.Rd.
summary.qap_cor <- function(object, ...) {
  q <- stats::quantile(object$dist, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),
                       na.rm = TRUE, names = TRUE)
  structure(c(object, list(quantiles = q)), class = c("summary.qap_cor", "qap_cor"))
}

# See man/qap_cor.Rd.
print.summary.qap_cor <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.qap_cor(x, digits = digits)
  cat("\nReference distribution of r:\n")
  print(x$quantiles, digits = digits)
  invisible(x)
}

# Multiple regression QAP (man/mrqap.Rd).
mrqap <- function(y, x, reps = 1000, method = c("dsp", "y"), seed = NULL,
                  intercept = TRUE, mode = c("digraph", "graph"), diag = FALSE,
                  statistic = c("t", "beta"), exact = FALSE, attr = NULL,
                  tol = 1e-7, threads = 0L) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  statistic <- match.arg(statistic)
  intercept <- check_flag(intercept, "intercept")
  diag <- check_flag(diag, "diag")
  exact <- check_flag(exact, "exact")
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol)) {
    stop("`tol` must be a single number in [0, 1).", call. = FALSE)
  }
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
  seed <- qap_seed(seed)
  xs <- if (k) unlist(mats[-1], use.names = FALSE) else double(0)
  res <- engine_call("mrqap",
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

# See man/mrqap.Rd.
print.mrqap <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nMultiple regression QAP\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  tab <- cbind(Estimate = x$coefficients, `Pr(>=|b|)` = x$p.value[, 3])
  print(tab, digits = digits)
  cat("\n", qap_method_text(x), "\n", sep = "")
  invisible(x)
}

# See man/mrqap.Rd.
summary.mrqap <- function(object, ...) {
  tab <- cbind(Estimate = object$coefficients, `Std. Error` = object$std.error,
               `t value` = object$t, object$p.value)
  structure(c(object, list(table = tab)), class = c("summary.mrqap", "mrqap"))
}

# See man/mrqap.Rd.
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

# ---- helpers ----------------------------------------------------------------

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

qap_seed <- function(seed) {
  if (is.null(seed)) {
    # Two draws give 52 random bits, so set.seed() controls the result.
    return(floor(stats::runif(1L) * 2^26) * 2^26 + floor(stats::runif(1L) * 2^26))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0 ||
      seed >= 2^53 || seed != floor(seed)) {
    stop("`seed` must be a whole number in [0, 2^53).", call. = FALSE)
  }
  as.double(seed)
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
