# Goodness of fit of an exponential random graph model (man/ergm_gof.Rd).
ergm_gof <- function(object, nsim = 100, seed = NULL,
                     stats = c("degree", "esp", "dsp", "distance", "triad", "model"),
                     burnin = object$control$burnin, interval = object$control$interval,
                     threads = object$control$threads) {
  if (!inherits(object, "ergm_fit")) stop("`object` must come from ergm_fit().", call. = FALSE)
  stats <- match.arg(stats, several.ok = TRUE)
  model <- object$model
  n <- model$n
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1 || nsim != floor(nsim)) {
    stop("`nsim` must be a whole number of at least 1.", call. = FALSE)
  }
  threads <- resolve_threads(threads)
  seed <- check_seed(seed)
  ptr <- ergm_handle(model)
  res <- engine_call("ergm_gof",
    .Call(C_ergm_simulate, ptr, as.double(object$coefficients), as.double(burnin),
          as.double(max(interval, 1)), as.double(nsim), 1L,
          proposal_code(object$control$proposal), seed, threads, TRUE, FALSE),
    edge_rows = FALSE)
  sample <- res[[3]]
  out <- list()
  kinds <- list()
  if ("degree" %in% stats) {
    if (model$directed) {
      kinds$idegree <- list(code = 1L, len = n, labels = as.character(seq_len(n) - 1L))
      kinds$odegree <- list(code = 2L, len = n, labels = as.character(seq_len(n) - 1L))
    } else {
      kinds$degree <- list(code = 0L, len = n, labels = as.character(seq_len(n) - 1L))
    }
  }
  if ("esp" %in% stats) kinds$esp <- list(code = 3L, len = n - 1L, labels = as.character(seq_len(n - 1L) - 1L))
  if ("dsp" %in% stats) kinds$dsp <- list(code = 4L, len = n - 1L, labels = as.character(seq_len(n - 1L) - 1L))
  if ("distance" %in% stats) {
    kinds$distance <- list(code = 5L, len = n, labels = c(as.character(seq_len(n - 1L)), "Inf"))
  }
  if ("triad" %in% stats && model$directed) {
    kinds$triadcensus <- list(code = 6L, len = 16L,
                              labels = c("003", "012", "102", "021D", "021U", "021C", "111D", "111U",
                                         "030T", "030C", "201", "120D", "120U", "120C", "210", "300"))
  }
  for (name in names(kinds)) {
    k <- kinds[[name]]
    obs <- engine_call("ergm_gof", .Call(C_ergm_gof_observed, ptr, k$code, k$len), edge_rows = FALSE)
    sim <- engine_call("ergm_gof", .Call(C_ergm_gof_sample, sample, k$code, k$len), edge_rows = FALSE)
    colnames(sim) <- k$labels
    names(obs) <- k$labels
    out[[name]] <- gof_entry(obs, sim)
  }
  if ("model" %in% stats) {
    sim <- res[[1]]
    colnames(sim) <- model$names
    out$model <- gof_entry(object$stats, sim)
  }
  structure(c(out, list(nsim = nsim, seed = seed)), class = "ergm_gof")
}

gof_entry <- function(obs, sim) {
  p <- vapply(seq_along(obs), function(j) {
    lo <- mean(sim[, j] <= obs[j])
    hi <- mean(sim[, j] >= obs[j])
    min(1, 2 * min(lo, hi))
  }, 1)
  list(observed = obs, simulated = sim,
       summary = data.frame(obs = obs, min = apply(sim, 2, min), mean = colMeans(sim),
                            max = apply(sim, 2, max), p = p, row.names = names(obs)))
}

gof_stat_names <- function(x) setdiff(names(x), c("nsim", "seed"))

gof_title <- function(name) {
  switch(name,
    degree = "Degree", idegree = "In-degree", odegree = "Out-degree",
    esp = "Edgewise shared partners", dsp = "Dyadwise shared partners",
    distance = "Minimum geodesic distance", triadcensus = "Triad census",
    model = "Model statistics", name)
}

# See man/ergm_gof.Rd.
print.ergm_gof <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("\nGoodness of fit: %d simulated networks (seed %s)\n", as.integer(x$nsim),
              format(x$seed, scientific = FALSE)))
  for (name in gof_stat_names(x)) {
    cat("\n", gof_title(name), ":\n", sep = "")
    tab <- x[[name]]$summary
    names(tab)[5] <- "MC p-value"
    print(tab, digits = digits)
  }
  invisible(x)
}

# See man/ergm_gof.Rd.
plot.ergm_gof <- function(x, which = NULL, ask = NULL, ...) {
  names <- gof_stat_names(x)
  if (!is.null(which)) {
    bad <- setdiff(which, names)
    if (length(bad)) {
      stop(sprintf("Unknown statistic(s) %s; the fit has %s.",
                   paste0("\"", bad, "\"", collapse = ", "),
                   paste0("\"", names, "\"", collapse = ", ")), call. = FALSE)
    }
    names <- which
  }
  if (is.null(ask)) ask <- length(names) > 1L && grDevices::dev.interactive()
  old <- graphics::par(ask = ask)
  on.exit(graphics::par(old), add = TRUE)
  for (name in names) {
    e <- x[[name]]
    sim <- e$simulated
    obs <- e$observed
    # Model statistics have different scales; the others are counts on a
    # common scale, so trailing all-zero cells are dropped for legibility.
    if (name != "model") {
      keep <- seq_len(max(base::which(obs > 0 | colSums(sim) > 0), 1L))
      sim <- sim[, keep, drop = FALSE]
      obs <- obs[keep]
    }
    graphics::boxplot(sim, main = gof_title(name), xlab = if (name == "model") "" else name,
                      ylab = if (name == "model") "value" else "count", las = 2, cex.axis = 0.8, ...)
    graphics::lines(seq_along(obs), obs, lwd = 2)
    graphics::points(seq_along(obs), obs, pch = 16)
  }
  invisible(x)
}

# Markov chain diagnostics for a Monte Carlo fit (man/mcmc_diagnostics.Rd).
mcmc_diagnostics <- function(object, lags = c(1, 2, 5, 10, 20)) {
  if (!inherits(object, "ergm_fit")) stop("`object` must come from ergm_fit().", call. = FALSE)
  if (is.null(object$sample)) {
    stop(paste0("The fit has no Monte Carlo sample; diagnostics need estimate = \"MCMLE\"."),
         call. = FALSE)
  }
  if (!is.numeric(lags) || !length(lags) || anyNA(lags) || any(lags < 1) || any(lags != floor(lags))) {
    stop("`lags` must be whole numbers of at least 1.", call. = FALSE)
  }
  chains <- object$control$chains
  nsim <- as.integer(object$control$nsim)
  sample <- sweep(object$sample, 2, object$stats)
  p <- ncol(sample)
  z <- numeric(p)
  ess <- numeric(p)
  for (j in seq_len(p)) {
    zs <- numeric(chains)
    es <- numeric(chains)
    for (c in seq_len(chains)) {
      v <- sample[(c - 1L) * nsim + seq_len(nsim), j]
      zs[c] <- geweke_z(v)
      es[c] <- mcmc_effective_size(v)
    }
    # Chains are independent: their Geweke statistics combine as a sum
    # of independent normals, and effective sizes add.
    z[j] <- sum(zs) / sqrt(chains)
    ess[j] <- sum(es)
  }
  ac <- vapply(seq_len(p), function(j) {
    vapply(lags, function(l) {
      mean(vapply(seq_len(chains), function(c) {
        v <- sample[(c - 1L) * nsim + seq_len(nsim), j]
        if (l >= length(v) || stats::sd(v) == 0) return(NA_real_)
        stats::cor(v[seq_len(length(v) - l)], v[-seq_len(l)])
      }, 1))
    }, 1)
  }, numeric(length(lags)))
  ac <- matrix(ac, nrow = length(lags), dimnames = list(paste0("lag ", lags), colnames(sample)))
  structure(list(
    summary = data.frame(mean = colMeans(sample), sd = apply(sample, 2, stats::sd),
                         geweke_z = z, geweke_p = 2 * stats::pnorm(-abs(z)), ess = ess,
                         row.names = colnames(sample)),
    autocorrelation = ac,
    sample = sample,
    chains = chains,
    nsim = nsim
  ), class = "ergm_mcmc")
}

# Variance of the mean of a chain segment by batch means.
batch_variance <- function(v) {
  n <- length(v)
  b <- max(2L, floor(sqrt(n)))
  k <- n %/% b
  if (k < 1L) return(stats::var(v) / n)
  means <- vapply(seq_len(b), function(i) mean(v[(i - 1L) * k + seq_len(k)]), 1)
  stats::var(means) / b
}

geweke_z <- function(v) {
  n <- length(v)
  a <- v[seq_len(max(1L, floor(0.1 * n)))]
  b <- v[(n - floor(0.5 * n) + 1L):n]
  se <- sqrt(batch_variance(a) + batch_variance(b))
  if (!is.finite(se) || se == 0) return(0)
  (mean(a) - mean(b)) / se
}

mcmc_effective_size <- function(v) {
  n <- length(v)
  s2 <- stats::var(v)
  vm <- batch_variance(v)
  if (!is.finite(s2) || s2 == 0 || !is.finite(vm) || vm == 0) return(n)
  min(n, s2 / vm)
}

# See man/mcmc_diagnostics.Rd.
print.ergm_mcmc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("\nMarkov chain diagnostics: %d chain%s of %d samples\n\n", x$chains,
              if (x$chains == 1L) "" else "s", x$nsim))
  cat("Sample statistics minus observed:\n")
  print(x$summary, digits = digits)
  cat("\nAutocorrelation:\n")
  print(x$autocorrelation, digits = digits)
  invisible(x)
}

# See man/mcmc_diagnostics.Rd.
plot.ergm_mcmc <- function(x, which = NULL, ask = NULL, ...) {
  stats <- colnames(x$sample)
  if (!is.null(which)) stats <- if (is.numeric(which)) stats[which] else intersect(which, stats)
  if (!length(stats)) stop("No statistics to plot.", call. = FALSE)
  if (is.null(ask)) ask <- length(stats) > 1L && grDevices::dev.interactive()
  old <- graphics::par(mfrow = c(1, 3), ask = ask)
  on.exit(graphics::par(old), add = TRUE)
  for (s in stats) {
    v <- x$sample[, s]
    graphics::plot(seq_along(v), v, type = "l", xlab = "sample", ylab = s, main = "Trace")
    for (c in seq_len(x$chains)[-1]) graphics::abline(v = (c - 1L) * x$nsim + 0.5, lty = 3)
    graphics::abline(h = 0, col = "grey50")
    if (stats::sd(v) > 0) {
      graphics::plot(stats::density(v), main = "Density", xlab = s)
    } else {
      graphics::plot(c(v[1], v[1]), c(0, 1), type = "h", main = "Density", xlab = s, ylab = "")
    }
    graphics::abline(v = 0, col = "grey50")
    first <- v[seq_len(x$nsim)]
    if (stats::sd(first) > 0) {
      stats::acf(first, main = "Autocorrelation (first chain)")
    } else {
      graphics::plot.new()
      graphics::title(main = "Autocorrelation (constant)")
    }
  }
  invisible(x)
}
