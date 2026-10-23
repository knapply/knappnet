# Exponential random graph models (man/ergm_fit.Rd).
ergm_fit <- function(formula, estimate = c("MCMLE", "MPLE", "exact"), constraints = NULL,
                     control = ergm_control(), start = NULL) {
  estimate <- match.arg(estimate)
  check_ergm_control(control)
  model <- ergm_model(formula, constraints)
  p <- model$width
  threads <- resolve_threads(control$threads)
  mple <- ergm_mple_fit(model, control)
  call <- match.call()

  if (model$dyad_independent && estimate == "MCMLE") estimate <- "MPLE"

  if (!is.null(start)) {
    if (!is.numeric(start) || length(start) != p || anyNA(start) || !all(is.finite(start))) {
      stop(sprintf("`start` must be a finite numeric vector of length %d.", p), call. = FALSE)
    }
    start <- as.double(start)
  } else {
    start <- mple$coefficients
    if (!mple$converged || any(!is.finite(start)) || any(abs(start) > 10)) {
      start <- rep(0, p)
      if (estimate != "MPLE") {
        warning(paste0("The pseudolikelihood has no finite maximizer (a statistic separates ",
                       "ties from non-ties), so estimation starts from zero coefficients."),
                call. = FALSE)
      }
    }
  }

  fit <- switch(estimate,
    MPLE = {
      if (!mple$converged) {
        warning(paste0("The pseudolikelihood did not converge in ", control$mple_max_iter,
                       " iterations; a statistic may separate ties from non-ties, in which ",
                       "case no finite estimate exists."), call. = FALSE)
      }
      list(coefficients = mple$coefficients, vcov = mple$vcov, mcmc_vcov = NULL,
           loglik = if (model$dyad_independent) mple$loglik else NULL,
           loglik_se = if (model$dyad_independent) 0 else NULL,
           pseudo_loglik = mple$loglik, converged = mple$converged,
           iterations = mple$iterations, sample = NULL)
    },
    exact = {
      ptr <- ergm_handle(model)
      res <- engine_call("ergm_fit", .Call(C_ergm_exact_mle, ptr, start, 200L), edge_rows = FALSE)
      list(coefficients = res[[1]], vcov = res[[2]], mcmc_vcov = NULL, loglik = res[[3]][1],
           loglik_se = 0, pseudo_loglik = mple$loglik, converged = TRUE,
           iterations = as.integer(res[[3]][2]), sample = NULL)
    },
    MCMLE = {
      ptr <- ergm_handle(model)
      seed <- qap_seed(control$seed)
      critical <- stats::qchisq(control$confidence, df = p)
      res <- engine_call("ergm_fit",
        .Call(C_ergm_mcmle, ptr, start, as.double(control$burnin), as.double(control$interval),
              as.double(control$nsim), as.integer(control$chains),
              proposal_code(control$proposal), seed, threads, as.integer(control$max_iter),
              as.double(control$margin), critical), edge_rows = FALSE)
      info <- res[[4]]
      if (info[2] != 1) {
        warning(sprintf(paste0("The Monte Carlo MLE did not converge in %d iterations (last ",
                               "step length %.2f); increase `max_iter`, `nsim`, or `burnin`, ",
                               "or reconsider the model."), as.integer(info[1]), info[3]),
                call. = FALSE)
      }
      ll <- NULL
      ll_se <- NULL
      if (control$loglik) {
        lr <- engine_call("ergm_fit",
          .Call(C_ergm_loglik, ptr, res[[1]], as.integer(control$loglik_steps),
                as.double(control$burnin), as.double(control$interval), as.double(control$nsim),
                as.integer(control$chains), proposal_code(control$proposal),
                (seed + 1) %% 2^53, threads), edge_rows = FALSE)
        ll <- lr[1] - model$dyads * log(2)
        ll_se <- lr[2]
      }
      sample <- res[[5]]
      colnames(sample) <- model$names
      list(coefficients = res[[1]], vcov = res[[2]] + res[[3]], mcmc_vcov = res[[3]],
           loglik = ll, loglik_se = ll_se, pseudo_loglik = mple$loglik,
           converged = info[2] == 1, iterations = as.integer(info[1]), sample = sample,
           seed = seed, step_length = info[3], statistic = info[4])
    }
  )
  names(fit$coefficients) <- model$names
  dimnames(fit$vcov) <- list(model$names, model$names)
  if (!is.null(fit$mcmc_vcov)) dimnames(fit$mcmc_vcov) <- list(model$names, model$names)
  structure(c(fit, list(
    estimate = estimate,
    stats = stats::setNames(model$stats, model$names),
    null_loglik = -model$dyads * log(2),
    dyads = model$dyads,
    mple = mple,
    model = model,
    control = control,
    call = call
  )), class = "ergm_fit")
}

# See man/ergm_fit.Rd.
ergm_stats <- function(formula) {
  model <- ergm_model(formula, NULL)
  stats::setNames(model$stats, model$names)
}

# See man/ergm_fit.Rd.
ergm_control <- function(seed = NULL, threads = 0L, burnin = 10000, interval = 100,
                         nsim = 1000, chains = 2L, proposal = c("tnt", "random"),
                         max_iter = 60L, margin = 0.05, confidence = 0.95, loglik = TRUE,
                         loglik_steps = 16L, mple_max_iter = 50L, mple_tol = 1e-12) {
  ctl <- list(seed = seed, threads = threads, burnin = burnin, interval = interval, nsim = nsim,
              chains = chains, proposal = match.arg(proposal), max_iter = max_iter,
              margin = margin, confidence = confidence, loglik = loglik,
              loglik_steps = loglik_steps, mple_max_iter = mple_max_iter, mple_tol = mple_tol)
  check_ergm_control(ctl)
  ctl
}

check_ergm_control <- function(ctl) {
  if (!is.list(ctl)) stop("`control` must come from ergm_control().", call. = FALSE)
  whole <- function(v, name, min = 0) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < min || v != floor(v)) {
      stop(sprintf("`%s` must be a whole number of at least %d.", name, min), call. = FALSE)
    }
  }
  whole(ctl$burnin, "burnin")
  whole(ctl$interval, "interval", 1)
  whole(ctl$nsim, "nsim", 1)
  whole(ctl$chains, "chains", 1)
  whole(ctl$max_iter, "max_iter", 1)
  whole(ctl$loglik_steps, "loglik_steps", 1)
  whole(ctl$mple_max_iter, "mple_max_iter", 1)
  check_choice(ctl$proposal, c("tnt", "random"), "proposal")
  check_flag(ctl$loglik, "loglik")
  if (!is.numeric(ctl$margin) || length(ctl$margin) != 1L || is.na(ctl$margin) ||
      ctl$margin < 0 || ctl$margin >= 1) {
    stop("`margin` must be a number in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(ctl$confidence) || length(ctl$confidence) != 1L || is.na(ctl$confidence) ||
      ctl$confidence <= 0 || ctl$confidence >= 1) {
    stop("`confidence` must be a number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(ctl$mple_tol) || length(ctl$mple_tol) != 1L || is.na(ctl$mple_tol) ||
      ctl$mple_tol <= 0) {
    stop("`mple_tol` must be a positive number.", call. = FALSE)
  }
  if (!is.null(ctl$seed)) qap_seed(ctl$seed)
  resolve_threads(ctl$threads)
  invisible(ctl)
}

proposal_code <- function(proposal) if (proposal == "random") 1L else 0L

# ---- Methods -----------------------------------------------------------------

# See man/ergm_fit.Rd.
print.ergm_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nExponential random graph model\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(estimate_text(x), "\n", sep = "")
  cat("\nCoefficients:\n")
  print(x$coefficients, digits = digits)
  invisible(x)
}

estimate_text <- function(x) {
  switch(x$estimate,
    MPLE = if (x$model$dyad_independent) {
      "Maximum likelihood (dyad-independent model, fitted by logistic regression)"
    } else {
      "Maximum pseudolikelihood (standard errors are not valid for dependent terms)"
    },
    exact = sprintf("Exact maximum likelihood by enumeration of %s networks",
                    format(2^x$dyads, big.mark = ",", scientific = FALSE)),
    MCMLE = sprintf(paste0("Monte Carlo maximum likelihood: %d iteration%s, %d chain%s of %d ",
                           "samples (interval %d, burn-in %d, %s proposals, seed %s)%s"),
                    x$iterations, if (x$iterations == 1L) "" else "s", x$control$chains,
                    if (x$control$chains == 1L) "" else "s", as.integer(x$control$nsim),
                    as.integer(x$control$interval), as.integer(x$control$burnin),
                    if (x$control$proposal == "tnt") "tie/no-tie" else "random-dyad",
                    format(x$seed, scientific = FALSE),
                    if (x$converged) "" else "; NOT CONVERGED")
  )
}

# See man/ergm_fit.Rd.
coef.ergm_fit <- function(object, ...) object$coefficients

# See man/ergm_fit.Rd.
vcov.ergm_fit <- function(object, ...) object$vcov

# See man/ergm_fit.Rd.
logLik.ergm_fit <- function(object, ...) {
  if (is.null(object$loglik)) {
    stop(paste0(
      if (object$estimate == "MPLE") {
        paste0("The log-likelihood of a pseudolikelihood fit with dependent terms is not ",
               "available (`object$pseudo_loglik` holds the maximized pseudolikelihood); ")
      } else {
        "The log-likelihood was not estimated; "
      },
      "fit with estimate = \"MCMLE\" and control = ergm_control(loglik = TRUE), or ",
      "estimate = \"exact\" on a small network."), call. = FALSE)
  }
  structure(object$loglik, df = length(object$coefficients), nobs = object$dyads,
            se = object$loglik_se, class = "logLik")
}

# See man/ergm_fit.Rd.
summary.ergm_fit <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se)
  if (!is.null(object$mcmc_vcov)) {
    tab <- cbind(tab, `MCMC %` = round(100 * diag(object$mcmc_vcov) / diag(object$vcov)))
  }
  tab <- cbind(tab, `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(tab) <- names(est)
  p <- length(est)
  ll <- object$loglik
  structure(list(
    call = object$call,
    estimate = object$estimate,
    coefficients = tab,
    dyads = object$dyads,
    null_deviance = -2 * object$null_loglik,
    deviance = if (!is.null(ll)) -2 * ll else NULL,
    aic = if (!is.null(ll)) -2 * ll + 2 * p else NULL,
    bic = if (!is.null(ll)) -2 * ll + p * log(object$dyads) else NULL,
    loglik_se = object$loglik_se,
    pseudo_loglik = object$pseudo_loglik,
    dyad_independent = object$model$dyad_independent,
    fit = object
  ), class = "summary.ergm_fit")
}

# See man/ergm_fit.Rd.
print.summary.ergm_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nExponential random graph model\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(estimate_text(x$fit), "\n\n", sep = "")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, P.values = TRUE, has.Pvalue = TRUE,
                      signif.legend = TRUE)
  if (x$estimate == "MPLE" && !x$dyad_independent) {
    cat("\nStandard errors, z values, and p-values come from the logistic regression of the\n",
        "pseudolikelihood and are not valid for dependent terms.\n", sep = "")
  }
  cat(sprintf("\n%s dyads; null deviance %s on %d degrees of freedom\n",
              format(x$dyads, big.mark = ","), format(x$null_deviance, digits = digits),
              as.integer(x$dyads)))
  if (!is.null(x$deviance)) {
    cat(sprintf("Residual deviance %s on %d degrees of freedom%s\n",
                format(x$deviance, digits = digits), as.integer(x$dyads - nrow(x$coefficients)),
                if (!is.null(x$loglik_se) && x$loglik_se > 0) {
                  sprintf(" (log-likelihood Monte Carlo s.e. %s)", format(x$loglik_se, digits = 2))
                } else ""))
    cat(sprintf("AIC: %s   BIC: %s\n", format(x$aic, digits = digits), format(x$bic, digits = digits)))
  } else if (!is.null(x$pseudo_loglik)) {
    cat(sprintf("Maximized pseudolikelihood: %s (not a likelihood; no AIC or BIC)\n",
                format(x$pseudo_loglik, digits = digits)))
  }
  invisible(x)
}

# Simulate networks from an exponential random graph model (man/simulate.ergm_fit.Rd).
simulate.ergm_fit <- function(object, nsim = 1, seed = NULL, ..., coef = NULL,
                              burnin = object$control$burnin,
                              interval = object$control$interval,
                              proposal = object$control$proposal, threads = object$control$threads,
                              output = c("netframe", "stats")) {
  output <- match.arg(output)
  theta <- if (is.null(coef)) object$coefficients else coef
  ergm_draw(object$model, theta, nsim, seed, burnin, interval, proposal, threads, output)
}

# See man/simulate.ergm_fit.Rd.
ergm_simulate <- function(formula, coef, nsim = 1, seed = NULL, burnin = 10000, interval = 100,
                          proposal = c("tnt", "random"), constraints = NULL, threads = 0L,
                          output = c("netframe", "stats")) {
  output <- match.arg(output)
  proposal <- match.arg(proposal)
  model <- ergm_model(formula, constraints)
  ergm_draw(model, coef, nsim, seed, burnin, interval, proposal, threads, output)
}

ergm_draw <- function(model, theta, nsim, seed, burnin, interval, proposal, threads, output) {
  p <- model$width
  if (!is.numeric(theta) || length(theta) != p || anyNA(theta) || !all(is.finite(theta))) {
    stop(sprintf("`coef` must be a finite numeric vector with one value per statistic (%d).", p),
         call. = FALSE)
  }
  check_choice(proposal, c("tnt", "random"), "proposal")
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1 || nsim != floor(nsim)) {
    stop("`nsim` must be a whole number of at least 1.", call. = FALSE)
  }
  for (v in c("burnin", "interval")) {
    val <- get(v)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val < 0 || val != floor(val)) {
      stop(sprintf("`%s` must be a non-negative whole number.", v), call. = FALSE)
    }
  }
  threads <- resolve_threads(threads)
  seed <- qap_seed(seed)
  ptr <- ergm_handle(model)
  res <- engine_call("simulate",
    .Call(C_ergm_simulate, ptr, as.double(theta), as.double(burnin), as.double(max(interval, 1)),
          as.double(nsim), 1L, proposal_code(proposal), seed, threads, output == "netframe", FALSE),
    edge_rows = FALSE)
  stats <- res[[1]]
  colnames(stats) <- model$names
  if (output == "stats") return(stats)
  nets <- .Call(C_ergm_sample_networks, res[[3]])
  out <- lapply(nets, function(m) edge_netframe(model$x, m))
  attr(out, "stats") <- stats
  attr(out, "seed") <- seed
  out
}

# A netframe with the nodes of `x` and the ties in the two-column matrix
# `m` (1-based node positions).
edge_netframe <- function(x, m) {
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(from = ids[m[, 1]], to = ids[m[, 2]], stringsAsFactors = FALSE)
  names(edges) <- c(x$roles$from, x$roles$to)
  netframe(edges, x$nodes, directed = x$directed, from = 1L, to = 2L, id = x$roles$id)
}

# Conditional tie probabilities from an exponential random graph model (man/predict.ergm_fit.Rd).
predict.ergm_fit <- function(object, type = c("response", "link"), coef = NULL, ...) {
  type <- match.arg(type)
  model <- object$model
  theta <- if (is.null(coef)) object$coefficients else coef
  if (!is.numeric(theta) || length(theta) != model$width || anyNA(theta)) {
    stop(sprintf("`coef` must be a numeric vector with one value per statistic (%d).", model$width),
         call. = FALSE)
  }
  cs <- ergm_changestats(model)
  eta <- drop(cs$x %*% as.double(theta))
  ids <- model$x$nodes[[model$x$roles$id]]
  out <- data.frame(from = ids[cs$tail], to = ids[cs$head], tie = cs$y == 1,
                    stringsAsFactors = FALSE)
  out[[if (type == "link") "link" else "p"]] <- if (type == "link") eta else 1 / (1 + exp(-eta))
  out
}

ergm_changestats <- function(model) {
  ptr <- ergm_handle(model)
  res <- engine_call("ergm_fit", .Call(C_ergm_changestats, ptr), edge_rows = FALSE)
  x <- res[[1]]
  colnames(x) <- model$names
  list(x = x, y = res[[2]], tail = res[[3]], head = res[[4]])
}

# ---- Pseudolikelihood ---------------------------------------------------------

ergm_mple_fit <- function(model, control) {
  ptr <- ergm_handle(model)
  res <- engine_call("ergm_fit",
    .Call(C_ergm_mple, ptr, as.integer(control$mple_max_iter), as.double(control$mple_tol)),
    edge_rows = FALSE)
  fit <- res[[3]]
  list(coefficients = stats::setNames(res[[1]], model$names),
       vcov = structure(res[[2]], dimnames = list(model$names, model$names)),
       deviance = fit[1], null_deviance = fit[2], loglik = -fit[1] / 2,
       iterations = as.integer(fit[4]), converged = fit[5] == 1)
}

# ---- Model construction -------------------------------------------------------

TERM_CODES <- c(
  edges = 0L, mutual = 1L, nodecov = 2L, nodeicov = 3L, nodeocov = 4L, nodefactor = 5L,
  nodeifactor = 6L, nodeofactor = 7L, nodematch = 8L, nodematch_diff = 9L, absdiff = 10L,
  edgecov = 11L, kstar = 12L, istar = 13L, ostar = 14L, triangle = 15L, ttriple = 16L,
  ctriple = 17L, gwesp = 18L, gwdsp = 19L, gwnsp = 20L, gwdegree = 21L, gwidegree = 22L,
  gwodegree = 23L, degree = 24L, idegree = 25L, odegree = 26L, isolates = 27L
)

# Splits the right-hand side of a formula into its summands.
rhs_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    if (length(expr) == 2L) return(rhs_terms(expr[[2]]))
    return(c(rhs_terms(expr[[2]]), rhs_terms(expr[[3]])))
  }
  if (is.call(expr) && identical(expr[[1]], as.name("("))) return(rhs_terms(expr[[2]]))
  list(expr)
}

ergm_formula_parts <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula such as `g ~ edges + triangle`.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame(2)
  list(lhs = eval(formula[[2]], env), terms = rhs_terms(formula[[3]]), env = env)
}

# Builds the term table and engine handle for a formula on netframe `x`.
# `lag` is the previous network of a temporal model, for memory terms.
ergm_model <- function(formula, constraints = NULL, x = NULL, lag = NULL, env = NULL) {
  parts <- ergm_formula_parts(formula)
  if (is.null(x)) x <- parts$lhs
  if (is.null(env)) env <- parts$env
  if (!inherits(x, "netframe")) {
    stop("The left-hand side of `formula` must be a netframe.", call. = FALSE)
  }
  check_netframe(x)
  n <- nrow(x$nodes)
  if (any(x$index$src == x$index$tgt)) {
    stop(paste0("The network has self-loops (", rows_text(which(x$index$src == x$index$tgt)),
                " of `x$edges`); exponential random graph models are defined on networks ",
                "without them."), call. = FALSE)
  }
  spec <- ergm_term_table(parts$terms, x, env, lag)
  fixed <- NULL
  if (!is.null(constraints)) {
    if (!is.matrix(constraints) || nrow(constraints) != n || ncol(constraints) != n ||
        !(is.logical(constraints) || is.numeric(constraints)) || anyNA(constraints)) {
      stop(sprintf("`constraints` must be an %d-by-%d logical or 0/1 matrix without NA.", n, n),
           call. = FALSE)
    }
    fixed <- as.integer(t(constraints) != 0) # row-major for the engine
  }
  model <- list(x = x, n = n, directed = x$directed, terms = spec$terms, names = spec$names,
                numeric = spec$numeric, categorical = spec$categorical, dyadic = spec$dyadic,
                fixed = fixed, cache = new.env(parent = emptyenv()))
  ptr <- ergm_handle(model)
  info <- .Call(C_ergm_info, ptr)
  model$width <- as.integer(info[1])
  model$dyads <- info[2]
  model$dyad_independent <- info[3] == 1
  model$stats <- engine_call("ergm_stats", .Call(C_ergm_stats, ptr), edge_rows = FALSE)
  if (length(model$names) != model$width) stop("statistic names do not match the engine (internal error)")
  model
}

# The engine handle of a model, rebuilt if it did not survive serialization.
ergm_handle <- function(model) {
  ptr <- model$cache$ptr
  if (!is.null(ptr) && .Call(C_ergm_alive, ptr)) return(ptr)
  x <- model$x
  t <- model$terms
  n <- model$n
  flat <- function(blocks, mode) {
    if (!length(blocks)) return(vector(mode, 0L))
    as.vector(unlist(lapply(blocks, function(b) as.vector(b)), use.names = FALSE), mode = mode)
  }
  levels_flat <- unlist(lapply(t, function(s) s$levels), use.names = FALSE)
  if (is.null(levels_flat)) levels_flat <- integer(0)
  lens <- vapply(t, function(s) length(s$levels), 1L)
  offs <- cumsum(c(0L, lens))[seq_along(t)]
  ptr <- engine_call("ergm_fit",
    .Call(C_ergm_new, n, x$index$src, x$index$tgt, x$directed,
          vapply(t, function(s) s$code, 1L), vapply(t, function(s) as.double(s$param), 1),
          vapply(t, function(s) s$attr, 1L), as.integer(offs), as.integer(lens),
          as.integer(levels_flat),
          flat(model$numeric, "double"), length(model$numeric),
          flat(model$categorical, "integer"), length(model$categorical),
          flat(model$dyadic, "double"), length(model$dyadic),
          model$fixed))
  model$cache$ptr <- ptr
  ptr
}

# Parses the terms of a formula into engine term specifications.
ergm_term_table <- function(terms, x, env, lag) {
  numeric <- list()
  categorical <- list()
  dyadic <- list()
  out <- list()
  names <- character(0)
  directed <- x$directed
  add <- function(code, param = 0, attr = -1L, levels = integer(0), name) {
    out[[length(out) + 1L]] <<- list(code = as.integer(code), param = param,
                                    attr = as.integer(attr), levels = as.integer(levels))
    names <<- c(names, name)
  }
  add_numeric <- function(attr, what) {
    v <- node_attribute(x, attr, what)
    if (!(is.numeric(v) || is.logical(v)) || anyNA(v)) {
      stop(sprintf("Term `%s` needs a numeric node attribute without missing values; `%s` is %s.",
                   what, attr, if (anyNA(v)) "missing for some nodes" else describe_class(v)),
           call. = FALSE)
    }
    numeric[[length(numeric) + 1L]] <<- as.double(v)
    length(numeric) - 1L
  }
  add_categorical <- function(attr, what) {
    v <- node_attribute(x, attr, what)
    if (anyNA(v)) {
      stop(sprintf("Term `%s` needs a node attribute without missing values; `%s` has NA.",
                   what, attr), call. = FALSE)
    }
    f <- if (is.factor(v)) factor(v) else factor(v, levels = sort(unique(v)))
    categorical[[length(categorical) + 1L]] <<- as.integer(f) - 1L
    list(index = length(categorical) - 1L, levels = levels(f))
  }
  add_dyadic <- function(m) {
    dyadic[[length(dyadic) + 1L]] <<- m
    length(dyadic) - 1L
  }
  need_directed <- function(what) {
    if (!directed) {
      stop(sprintf("Term `%s` is defined only for directed networks.", what), call. = FALSE)
    }
  }
  need_undirected <- function(what, alt) {
    if (directed) {
      stop(sprintf("Term `%s` is defined only for undirected networks; use %s.", what, alt),
           call. = FALSE)
    }
  }
  # The k-th unnamed argument, or NA.
  positional <- function(args, k) {
    nm <- names(args)
    if (is.null(nm)) nm <- rep("", length(args))
    which(nm == "")[k]
  }
  # A term option by name or position, as in statnet's term signatures.
  get_arg <- function(args, name, position, default = NULL) {
    nm <- names(args)
    if (!is.null(nm) && name %in% nm) return(args[[name]])
    i <- positional(args, position)
    if (length(i) && !is.na(i)) return(args[[i]])
    default
  }
  decay_of <- function(args, what) {
    fixed <- get_arg(args, "fixed", 2L, TRUE)
    if (!isTRUE(fixed)) {
      stop(sprintf(paste0("Term `%s` is supported only with `fixed = TRUE` (a fixed decay ",
                          "parameter); curved terms are not estimated."), what), call. = FALSE)
    }
    d <- get_arg(args, "decay", 1L, NULL)
    if (is.null(d)) {
      stop(sprintf("Term `%s` needs a `decay` parameter, as in `%s(0.5, fixed = TRUE)`.", what, what),
           call. = FALSE)
    }
    if (!is.numeric(d) || length(d) != 1L || !is.finite(d)) {
      stop(sprintf("The decay of `%s` must be a single finite number.", what), call. = FALSE)
    }
    as.double(d)
  }
  gw_name <- function(base, decay, otp) {
    sprintf("%s.%sfixed.%s", base, if (otp && directed) "OTP." else "", format(decay))
  }

  for (term in terms) {
    what <- if (is.name(term)) {
      as.character(term)
    } else if (is.call(term) && is.name(term[[1]])) {
      as.character(term[[1]])
    } else {
      NULL
    }
    known <- c(setdiff(names(TERM_CODES), "nodematch_diff"), "memory", "delrecip")
    if (is.null(what) || !what %in% known) {
      stop(sprintf("Unknown model term `%s`. See ?ergm_fit for the supported terms.",
                   paste(deparse(term), collapse = "")), call. = FALSE)
    }
    raw_args <- if (is.call(term)) as.list(term)[-1] else list()
    # Attribute names and matrices are evaluated in the formula's
    # environment; term options likewise.
    args <- lapply(raw_args, function(a) eval(a, env))
    nargs <- length(args)
    switch(what,
      edges = add(0L, name = "edges"),
      mutual = {
        need_directed("mutual")
        add(1L, name = "mutual")
      },
      nodecov = , nodeicov = , nodeocov = {
        if (what != "nodecov") need_directed(what)
        attr <- get_arg(args, "attr", 1L)
        k <- add_numeric(attr, what)
        add(TERM_CODES[[what]], attr = k, name = paste0(what, ".", attr))
      },
      nodefactor = , nodeifactor = , nodeofactor = {
        if (what != "nodefactor") need_directed(what)
        attr <- get_arg(args, "attr", 1L)
        cat <- add_categorical(attr, what)
        keep <- select_levels(cat$levels, get_arg(args, "levels", 2L, -1), what)
        if (!length(keep)) {
          stop(sprintf("Term `%s(\"%s\")` keeps no level; see `levels`.", what, attr), call. = FALSE)
        }
        # One term per level keeps names aligned with statistics.
        for (l in keep) {
          add(TERM_CODES[[what]], attr = cat$index, levels = l - 1L,
              name = paste0(what, ".", attr, ".", cat$levels[l]))
        }
      },
      nodematch = {
        attr <- get_arg(args, "attr", 1L)
        diff <- get_arg(args, "diff", 2L, FALSE)
        check_flag(diff, "diff")
        cat <- add_categorical(attr, what)
        if (diff) {
          keep <- select_levels(cat$levels, get_arg(args, "levels", 3L, NULL), what)
          for (l in keep) {
            add(9L, attr = cat$index, levels = l - 1L,
                name = paste0("nodematch.", attr, ".", cat$levels[l]))
          }
        } else {
          add(8L, attr = cat$index, name = paste0("nodematch.", attr))
        }
      },
      absdiff = {
        attr <- get_arg(args, "attr", 1L)
        pow <- get_arg(args, "pow", 2L, 1)
        if (!is.numeric(pow) || length(pow) != 1L || !is.finite(pow) || pow <= 0) {
          stop("The `pow` of `absdiff` must be a positive number.", call. = FALSE)
        }
        k <- add_numeric(attr, what)
        add(10L, param = pow, attr = k,
            name = paste0("absdiff", if (pow != 1) format(pow) else "", ".", attr))
      },
      edgecov = {
        cov <- get_arg(args, "x", 1L)
        attrname <- get_arg(args, "attrname", 2L, NULL)
        label <- if (!is.null(raw_args[["x"]])) paste(deparse(raw_args[["x"]]), collapse = "") else {
          i <- positional(raw_args, 1L)
          if (!is.na(i)) paste(deparse(raw_args[[i]]), collapse = "") else "x"
        }
        if (!is.null(attrname)) label <- attrname
        m <- dyadic_matrix(cov, attrname, nrow(x$nodes), x)
        add(11L, attr = add_dyadic(m), name = paste0("edgecov.", label))
      },
      kstar = , istar = , ostar = {
        if (what == "kstar") need_undirected(what, "`istar()` or `ostar()`") else need_directed(what)
        ks <- get_arg(args, "k", 1L)
        check_counts(ks, what, min = 1)
        for (k in ks) add(TERM_CODES[[what]], param = k, name = paste0(what, k))
      },
      degree = , idegree = , odegree = {
        if (what == "degree") need_undirected(what, "`idegree()` or `odegree()`") else need_directed(what)
        ds <- get_arg(args, "d", 1L)
        check_counts(ds, what, min = 0)
        for (d in ds) add(TERM_CODES[[what]], param = d, name = paste0(what, d))
      },
      triangle = add(15L, name = "triangle"),
      ttriple = {
        need_directed(what)
        add(16L, name = "ttriple")
      },
      ctriple = {
        need_directed(what)
        add(17L, name = "ctriple")
      },
      gwesp = , gwdsp = , gwnsp = {
        d <- decay_of(args, what)
        add(TERM_CODES[[what]], param = d, name = gw_name(what, d, TRUE))
      },
      gwdegree = , gwidegree = , gwodegree = {
        if (what == "gwdegree") need_undirected(what, "`gwidegree()` or `gwodegree()`") else need_directed(what)
        d <- decay_of(args, what)
        base <- c(gwdegree = "gwdeg", gwidegree = "gwideg", gwodegree = "gwodeg")[[what]]
        add(TERM_CODES[[what]], param = d, name = gw_name(base, d, FALSE))
      },
      isolates = add(27L, name = "isolates"),
      memory = , delrecip = {
        if (is.null(lag)) {
          stop(sprintf("Term `%s` is a temporal term; use it in tergm_fit().", what), call. = FALSE)
        }
        prev <- as_binary_matrix(lag)
        if (what == "memory") {
          type <- get_arg(args, "type", 1L, "stability")
          check_choice(type, c("stability", "autoregression", "innovation", "loss"), "type")
          m <- switch(type,
            stability = 2 * prev - 1,
            autoregression = prev,
            innovation = 1 - prev,
            loss = -prev
          )
          add(11L, attr = add_dyadic(m), name = paste0("edgecov.memory.", type))
        } else {
          need_directed("delrecip")
          add(11L, attr = add_dyadic(t(prev)), name = "edgecov.delrecip")
        }
      }
    )
    if (what %in% c("edges", "mutual", "triangle", "ttriple", "ctriple", "isolates") && nargs) {
      stop(sprintf("Term `%s` takes no arguments.", what), call. = FALSE)
    }
  }
  if (!length(out)) stop("The model needs at least one term.", call. = FALSE)
  list(terms = out, names = names, numeric = numeric, categorical = categorical, dyadic = dyadic)
}

node_attribute <- function(x, attr, what) {
  if (!is.character(attr) || length(attr) != 1L || is.na(attr)) {
    stop(sprintf("Term `%s` needs the name of a node attribute (a column of `x$nodes`) as a string.",
                 what), call. = FALSE)
  }
  if (is.null(x$nodes[[attr]])) {
    stop(sprintf("Node attribute `%s` (term `%s`) is not a column of `x$nodes`; it has: %s.",
                 attr, what, paste0("`", names(x$nodes), "`", collapse = ", ")), call. = FALSE)
  }
  x$nodes[[attr]]
}

# Positions (1-based) of the levels a factor term keeps.
select_levels <- function(lv, spec, what) {
  k <- length(lv)
  if (is.null(spec) || isTRUE(spec)) return(seq_len(k))
  if (is.numeric(spec) && length(spec) && !anyNA(spec)) {
    if (all(spec < 0)) return(setdiff(seq_len(k), -spec))
    if (all(spec > 0) && all(spec <= k)) return(as.integer(spec))
  }
  if (is.character(spec) && !anyNA(spec)) {
    pos <- match(spec, lv)
    if (anyNA(pos)) {
      stop(sprintf("Term `%s`: level(s) %s not found; the levels are %s.", what,
                   paste0("\"", spec[is.na(pos)], "\"", collapse = ", "),
                   paste0("\"", lv, "\"", collapse = ", ")), call. = FALSE)
    }
    return(pos)
  }
  stop(sprintf(paste0("Term `%s`: `levels` must be NULL (all levels), positions among the ",
                      "sorted levels (negative to drop), or level names."), what), call. = FALSE)
}

check_counts <- function(v, what, min) {
  if (!is.numeric(v) || !length(v) || anyNA(v) || any(v < min) || any(v != floor(v))) {
    stop(sprintf("Term `%s` needs whole numbers of at least %d.", what, min), call. = FALSE)
  }
}

# A dyadic covariate as an n-by-n matrix aligned with x$nodes.
dyadic_matrix <- function(cov, attrname, n, x) {
  if (inherits(cov, "netframe")) {
    if (nrow(cov$nodes) != n) {
      stop("An `edgecov` netframe must have the same nodes as the model's network.", call. = FALSE)
    }
    ids <- as.character(x$nodes[[x$roles$id]])
    m <- as.matrix(cov, attr = attrname)
    if (!is.null(rownames(m)) && setequal(rownames(m), ids) && !anyDuplicated(ids)) {
      m <- m[ids, ids, drop = FALSE]
    }
    cov <- m
  }
  if (!is.matrix(cov) || nrow(cov) != n || ncol(cov) != n || !(is.numeric(cov) || is.logical(cov))) {
    stop(sprintf("`edgecov` needs an %d-by-%d numeric matrix or a netframe on the same nodes.", n, n),
         call. = FALSE)
  }
  if (anyNA(cov) || any(!is.finite(cov))) {
    stop("`edgecov` matrices must not contain missing or infinite values.", call. = FALSE)
  }
  storage.mode(cov) <- "double"
  unname(cov)
}

as_binary_matrix <- function(x) {
  m <- as.matrix(x)
  m <- unname((m != 0) * 1)
  diag(m) <- 0
  m
}
