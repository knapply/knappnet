#' Temporal exponential random graph models by pooled pseudolikelihood
#'
#' Fits a discrete-time temporal ERGM (TERGM; Hanneke, Fu, and Xing 2010)
#' to a sequence of networks on the same nodes by maximum pseudolikelihood
#' pooled over the transitions, with bootstrap confidence intervals from
#' resampling the transitions (Desmarais and Cranmer 2012), the method of
#' the `btergm` package (Leifeld, Cranmer, and Desmarais 2018).
#'
#' @param formula A formula whose left-hand side is a list of
#'   [netframe()]s in time order (for example from [snapshots()]), all with
#'   the same nodes in the same order, and whose right-hand side sums the
#'   model terms of [ergm_fit()] plus the temporal terms below.
#' @param boot Number of bootstrap replicates.
#' @param seed Seed for the bootstrap, a whole number in \[0, 2^53); `NULL`
#'   draws one from R's random number generator.
#' @param level Confidence level of the percentile intervals.
#' @param threads Number of threads; `0` uses all available cores.
#' @param control Pseudolikelihood settings (`mple_max_iter`, `mple_tol`)
#'   from [ergm_control()].
#' @param object,x A `tergm_fit` object.
#' @param parm Coefficients to report (names or positions); default all.
#' @param ... Not used.
#'
#' @details
#' Each network \eqn{y_t} (\eqn{t = 2, \ldots, T}) is modeled conditional
#' on the previous one,
#' \deqn{P(Y_t = y_t | y_{t-1}) = \exp(\theta' g(y_t, y_{t-1})) / c(\theta, y_{t-1}),}{P(Y_t = y_t | y_{t-1}) = exp(theta' g(y_t, y_{t-1})) / c(theta, y_{t-1}),}
#' with the usual ERGM terms evaluated on \eqn{y_t} and temporal terms that
#' read \eqn{y_{t-1}} as a dyadic covariate:
#' \describe{
#'   \item{`memory(type = "stability")`}{The number of dyads in the same
#'     state as at the previous time (ties that persist and non-ties that
#'     persist); as a change statistic this is `edgecov(2 * y_prev - 1)`.
#'     `type = "autoregression"` counts persisting ties (`edgecov(y_prev)`),
#'     `"innovation"` new ties (`edgecov(1 - y_prev)`), and `"loss"`
#'     dissolved ties (`edgecov(-y_prev)`).}
#'   \item{`delrecip`}{Delayed reciprocity: ties \eqn{i \to j}{i -> j} whose
#'     reverse \eqn{j \to i}{j -> i} existed at the previous time
#'     (`edgecov(t(y_prev))`; directed networks).}
#' }
#' The pooled pseudolikelihood stacks the change statistics of every free
#' dyad of every transition and fits one logistic regression; its point
#' estimates equal `btergm`'s. Because the dyads are dependent, the
#' logistic standard errors are not used; instead the transitions are
#' resampled with replacement `boot` times, the pooled fit is repeated on
#' each resample, and confidence intervals are the percentiles of the
#' replicates (with few transitions the intervals are coarse, as the method
#' implies). Replicate \eqn{r} draws from stream \eqn{r} of `seed`, so
#' results do not depend on `threads`.
#'
#' @return An object of class `"tergm_fit"`: a list with `coefficients`,
#'   `boot` (a boot-by-p matrix of replicate estimates, `NA` where a
#'   replicate did not converge), `ci` (percentile intervals), `level`,
#'   `transitions`, `dyads` (per transition), `deviance`, `null_deviance`,
#'   `seed`, `names`, and `call`. Methods: `print()`, `summary()`,
#'   `coef()`, and `confint()`.
#'
#' @references
#' Hanneke, S., Fu, W., & Xing, E. P. (2010). Discrete temporal models of
#' social networks. *Electronic Journal of Statistics*, 4, 585-605.
#' \doi{10.1214/09-EJS548}
#'
#' Desmarais, B. A., & Cranmer, S. J. (2012). Statistical mechanics of
#' networks: Estimation and uncertainty. *Physica A*, 391(4), 1865-1876.
#' \doi{10.1016/j.physa.2011.10.018}
#'
#' Leifeld, P., Cranmer, S. J., & Desmarais, B. A. (2018). Temporal
#' exponential random graph models with btergm: Estimation and bootstrap
#' confidence intervals. *Journal of Statistical Software*, 83(6), 1-36.
#' \doi{10.18637/jss.v083.i06}
#'
#' @seealso [ergm_fit()], [snapshots()]
#' @export
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

#' @rdname tergm_fit
#' @export
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

#' @rdname tergm_fit
#' @export
summary.tergm_fit <- function(object, ...) {
  tab <- cbind(Estimate = object$coefficients, object$ci)
  valid <- sum(stats::complete.cases(object$boot))
  structure(list(call = object$call, coefficients = tab, transitions = object$transitions,
                 dyads = object$dyads, boot = nrow(object$boot), valid = valid,
                 level = object$level, seed = object$seed), class = "summary.tergm_fit")
}

#' @rdname tergm_fit
#' @export
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

#' @export
coef.tergm_fit <- function(object, ...) object$coefficients

#' @rdname tergm_fit
#' @export
confint.tergm_fit <- function(object, parm, level = NULL, ...) {
  ci <- if (is.null(level) || level == object$level) object$ci else tergm_ci(object$boot, level)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
