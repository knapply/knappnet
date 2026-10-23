# Temporal exponential random graph models by pooled pseudolikelihood (man/tergm_fit.Rd).
tergm_fit <- function(formula, boot = 200, seed = NULL, level = 0.95, threads = 0L,
                      control = ergm_control()) {
  parts <- ergm_formula_parts(formula)
  nets <- parts$lhs
  if (!is.list(nets) || inherits(nets, "netframe") || length(nets) < 2L ||
      !all(vapply(nets, inherits, TRUE, "netframe"))) {
    stop("The left-hand side of `formula` must be a list of at least two netframes.", call. = FALSE)
  }
  n <- nrow(nets[[1]]$nodes)
  ids <- nets[[1]]$nodes[[nets[[1]]$roles$id]]
  directed <- nets[[1]]$directed
  for (t in seq_along(nets)) {
    check_netframe(nets[[t]])
    if (nrow(nets[[t]]$nodes) != n || !identical(as.character(nets[[t]]$nodes[[nets[[t]]$roles$id]]),
                                                 as.character(ids))) {
      stop(sprintf("Network %d has different nodes (or a different order) than the first.", t),
           call. = FALSE)
    }
    if (nets[[t]]$directed != directed) {
      stop("All networks must be directed, or all undirected.", call. = FALSE)
    }
  }
  if (!is.numeric(boot) || length(boot) != 1L || is.na(boot) || boot < 0 || boot != floor(boot)) {
    stop("`boot` must be a non-negative whole number.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number in (0, 1).", call. = FALSE)
  }
  check_ergm_control(control)
  threads <- resolve_threads(threads)
  seed <- qap_seed(seed)

  xs <- list()
  ys <- list()
  period <- list()
  dyads <- integer(0)
  names <- NULL
  for (t in 2:length(nets)) {
    model <- ergm_model(formula, NULL, x = nets[[t]], lag = nets[[t - 1L]], env = parts$env)
    cs <- ergm_changestats(model)
    if (is.null(names)) names <- model$names
    xs[[t - 1L]] <- cs$x
    ys[[t - 1L]] <- cs$y
    period[[t - 1L]] <- rep(t - 2L, length(cs$y))
    dyads <- c(dyads, length(cs$y))
  }
  x <- do.call(rbind, xs)
  y <- unlist(ys, use.names = FALSE)
  per <- as.integer(unlist(period, use.names = FALSE))
  res <- engine_call("tergm_fit",
    .Call(C_tergm_mple, x, as.double(y), per, length(nets) - 1L, as.double(boot), seed, threads,
          as.integer(control$mple_max_iter), as.double(control$mple_tol)), edge_rows = FALSE)
  coef <- stats::setNames(res[[1]], names)
  fit <- res[[3]]
  if (fit[5] != 1) {
    warning(paste0("The pooled pseudolikelihood did not converge; a statistic may separate ties ",
                   "from non-ties, in which case no finite estimate exists."), call. = FALSE)
  }
  bootm <- res[[4]]
  colnames(bootm) <- names
  ci <- tergm_ci(bootm, level)
  structure(list(
    coefficients = coef,
    boot = bootm,
    ci = ci,
    level = level,
    transitions = length(nets) - 1L,
    dyads = dyads,
    deviance = fit[1],
    null_deviance = fit[2],
    seed = seed,
    names = names,
    call = match.call()
  ), class = "tergm_fit")
}

tergm_ci <- function(bootm, level) {
  a <- (1 - level) / 2
  ok <- stats::complete.cases(bootm)
  ci <- t(apply(bootm[ok, , drop = FALSE], 2, function(v) {
    if (!length(v)) return(c(NA_real_, NA_real_))
    stats::quantile(v, c(a, 1 - a), names = FALSE, type = 7)
  }))
  if (!ncol(bootm)) ci <- matrix(NA_real_, 0, 2)
  dimnames(ci) <- list(colnames(bootm), paste0(format(100 * c(a, 1 - a), trim = TRUE), " %"))
  ci
}

# See man/tergm_fit.Rd.
print.tergm_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nTemporal exponential random graph model (pooled pseudolikelihood)\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(sprintf("%d transition%s, %s dyads; %d bootstrap replicates (seed %s)\n",
              x$transitions, if (x$transitions == 1L) "" else "s",
              format(sum(x$dyads), big.mark = ","), nrow(x$boot),
              format(x$seed, scientific = FALSE)))
  cat("\nCoefficients:\n")
  print(x$coefficients, digits = digits)
  invisible(x)
}

# See man/tergm_fit.Rd.
summary.tergm_fit <- function(object, ...) {
  tab <- cbind(Estimate = object$coefficients, object$ci)
  valid <- sum(stats::complete.cases(object$boot))
  structure(list(call = object$call, coefficients = tab, transitions = object$transitions,
                 dyads = object$dyads, boot = nrow(object$boot), valid = valid,
                 level = object$level, seed = object$seed), class = "summary.tergm_fit")
}

# See man/tergm_fit.Rd.
print.summary.tergm_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nTemporal exponential random graph model (pooled pseudolikelihood)\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(sprintf("%d transition%s, %s dyads\n", x$transitions, if (x$transitions == 1L) "" else "s",
              format(sum(x$dyads), big.mark = ",")))
  cat(sprintf("Bootstrap: %d replicates resampling transitions (%d converged; seed %s)\n\n",
              x$boot, x$valid, format(x$seed, scientific = FALSE)))
  cat(sprintf("Coefficients with %s%% percentile confidence intervals:\n", format(100 * x$level)))
  print(x$coefficients, digits = digits)
  invisible(x)
}

# See man/tergm_fit.Rd.
coef.tergm_fit <- function(object, ...) object$coefficients

# See man/tergm_fit.Rd.
confint.tergm_fit <- function(object, parm, level = NULL, ...) {
  ci <- if (is.null(level) || level == object$level) object$ci else tergm_ci(object$boot, level)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
