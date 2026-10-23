# Stochastic actor-oriented models: panel data, effects, estimation by the
# method of moments, simulation, goodness of fit, and the methods of the
# fit.

# Network panel data for stochastic actor-oriented models (man/saom_data.Rd).
saom_data <- function(networks, covariates = NULL, dyadic = NULL) {
  if (!is.list(networks) || length(networks) < 2L ||
      !all(vapply(networks, inherits, TRUE, "netframe"))) {
    stop("`networks` must be a list of at least two netframes observed at successive moments.",
         call. = FALSE)
  }
  for (k in seq_along(networks)) check_netframe(networks[[k]])
  if (!all(vapply(networks, function(x) x$directed, TRUE))) {
    stop(paste0("Every network must be directed; stochastic actor-oriented models for ",
                "non-directed networks are not implemented."), call. = FALSE)
  }
  first <- networks[[1L]]
  ids <- first$nodes[[first$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  n <- length(ids)
  if (n < 2L) stop("The networks need at least two actors.", call. = FALSE)
  waves <- vector("list", length(networks))
  for (k in seq_along(networks)) {
    x <- networks[[k]]
    xid <- x$nodes[[x$roles$id]]
    pos <- match(xid, ids)
    if (nrow(x$nodes) != n || anyNA(pos)) {
      stop(sprintf(paste0("Network %d has different nodes from the first network; every wave ",
                          "must observe the same actors (use the same `nodes` table)."), k),
           call. = FALSE)
    }
    a <- matrix(FALSE, n, n)
    src <- pos[x$index$src + 1L]
    tgt <- pos[x$index$tgt + 1L]
    keep <- src != tgt
    a[cbind(src[keep], tgt[keep])] <- TRUE
    dimnames(a) <- list(as.character(ids), as.character(ids))
    waves[[k]] <- a
  }
  cov <- saom_actor_covariates(covariates, networks, ids)
  dy <- saom_dyadic_covariates(dyadic, ids, n)
  structure(list(n = n, waves = waves, nodes = first$nodes, id = first$roles$id, ids = ids,
                 covariates = cov, dyadic = dy, cache = new.env(parent = emptyenv())),
            class = "saom_data")
}

saom_actor_covariates <- function(covariates, networks, ids) {
  n <- length(ids)
  waves <- length(networks)
  if (is.null(covariates)) return(list())
  out <- list()
  if (is.character(covariates)) {
    for (name in covariates) {
      cols <- lapply(networks, function(x) {
        if (!name %in% names(x$nodes)) {
          stop(sprintf("Covariate `%s` is not a column of every network's node table.", name),
               call. = FALSE)
        }
        v <- x$nodes[[name]][match(ids, x$nodes[[x$roles$id]])]
        if (!is.numeric(v) && !is.logical(v)) {
          stop(sprintf("Covariate `%s` must be numeric or logical, not %s.", name,
                       describe_class(v)), call. = FALSE)
        }
        as.double(v)
      })
      m <- do.call(cbind, cols)
      # The last wave's values are never used (they would apply to a period
      # after the last observation), so only the others decide whether the
      # covariate changes.
      used <- m[, -ncol(m), drop = FALSE]
      changing <- any(apply(used, 1L, function(r) length(unique(r)) > 1L))
      out[[name]] <- if (changing) used else m[, 1L]
    }
  } else if (is.list(covariates)) {
    if (is.null(names(covariates)) || !all(nzchar(names(covariates)))) {
      stop("`covariates` given as a list must be named.", call. = FALSE)
    }
    for (name in names(covariates)) {
      v <- covariates[[name]]
      if (is.matrix(v) && ncol(v) == 1L) v <- v[, 1L]
      if (is.matrix(v)) {
        if (nrow(v) != n || !ncol(v) %in% c(waves, waves - 1L)) {
          stop(sprintf(paste0("Changing covariate `%s` must have one row per actor (%d) and one ",
                              "column per period (%d)."), name, n, waves - 1L), call. = FALSE)
        }
        # One column per period; a column for the last wave is ignored.
        out[[name]] <- matrix(as.double(v[, seq_len(waves - 1L)]), n, waves - 1L)
      } else {
        if (!is.atomic(v) || length(v) != n) {
          stop(sprintf("Covariate `%s` must have one value per actor (%d).", name, n),
               call. = FALSE)
        }
        out[[name]] <- as.double(v)
      }
    }
  } else {
    stop("`covariates` must be NULL, a character vector of node columns, or a named list.",
         call. = FALSE)
  }
  for (name in names(out)) {
    if (anyNA(out[[name]]) || any(!is.finite(out[[name]]))) {
      stop(sprintf("Covariate `%s` has missing or non-finite values; these are not supported.",
                   name), call. = FALSE)
    }
  }
  out
}

saom_dyadic_covariates <- function(dyadic, ids, n) {
  if (is.null(dyadic)) return(list())
  if (!is.list(dyadic) || is.null(names(dyadic)) || !all(nzchar(names(dyadic)))) {
    stop("`dyadic` must be NULL or a named list of square matrices or netframes.", call. = FALSE)
  }
  out <- list()
  for (name in names(dyadic)) {
    m <- dyadic[[name]]
    if (inherits(m, "netframe")) {
      mid <- m$nodes[[m$roles$id]]
      if (length(mid) != n || anyNA(match(ids, mid))) {
        stop(sprintf("Dyadic covariate `%s` must be on the same actors as the networks.", name),
             call. = FALSE)
      }
      m <- as.matrix(m, attr = m$roles$weight)
      m <- m[as.character(ids), as.character(ids), drop = FALSE]
    }
    if (!is.matrix(m) || nrow(m) != n || ncol(m) != n || !(is.numeric(m) || is.logical(m))) {
      stop(sprintf("Dyadic covariate `%s` must be a numeric %d by %d matrix.", name, n, n),
           call. = FALSE)
    }
    m <- matrix(as.double(m), n, n)
    diag(m) <- 0
    if (any(!is.finite(m))) {
      stop(sprintf("Dyadic covariate `%s` has missing or non-finite values; these are not supported.",
                   name), call. = FALSE)
    }
    out[[name]] <- m
  }
  out
}

# See man/saom_data.Rd.
print.saom_data <- function(x, ...) {
  cat(sprintf("<saom_data> %d waves of a directed network on %d actors\n", length(x$waves), x$n))
  ties <- vapply(x$waves, function(a) as.double(sum(a)), 0)
  cat("ties per wave:", paste(ties, collapse = ", "), "\n")
  d <- vapply(seq_len(length(x$waves) - 1L),
              function(m) as.double(sum(x$waves[[m]] != x$waves[[m + 1L]])), 0)
  cat("tie variables changed per period:", paste(d, collapse = ", "), "\n")
  if (length(x$covariates)) {
    kind <- vapply(x$covariates, function(v) if (is.matrix(v)) "changing" else "constant", "")
    cat("actor covariates:", paste0(names(kind), " (", kind, ")", collapse = ", "), "\n")
  }
  if (length(x$dyadic)) cat("dyadic covariates:", paste(names(x$dyadic), collapse = ", "), "\n")
  invisible(x)
}

# The engine handle for a data object, built on first use.
saom_handle <- function(data) {
  if (!inherits(data, "saom_data")) {
    stop("`data` must be created by saom_data().", call. = FALSE)
  }
  # The handle is an external pointer, which does not survive saveRDS(); a
  # dead handle is rebuilt. The cache is an environment, shared by copies.
  ptr <- data$cache$ptr
  if (!is.null(ptr) && .Call(C_saom_alive, ptr)) return(ptr)
  n <- data$n
  adjacency <- as.double(unlist(lapply(data$waves, function(a) as.double(a)), use.names = FALSE))
  actor <- lapply(data$covariates, function(v) as.double(v))
  dyadic <- lapply(data$dyadic, function(m) as.double(m))
  ptr <- engine_call("saom_data",
    .Call(C_saom_new, n, length(data$waves), adjacency, unname(actor), unname(dyadic)),
    edge_rows = FALSE)
  data$cache$ptr <- ptr
  ptr
}

# ---- Effects -----------------------------------------------------------------

SAOM_EFFECTS <- data.frame(
  effect = c("density", "recip", "transTrip", "transTies", "transRecTrip", "cycle3",
             "gwespFF", "gwespFB", "gwespBF", "gwespBB", "inPop", "inPopSqrt", "outPop",
             "outPopSqrt", "inAct", "inActSqrt", "outAct", "outActSqrt", "egoX", "altX",
             "simX", "sameX", "egoXaltX", "X"),
  code = 0:23,
  needs = c(rep("none", 6L), rep("parameter", 4L), rep("none", 8L), rep("actor", 5L), "dyadic"),
  label = c("outdegree (density)", "reciprocity", "transitive triplets", "transitive ties",
            "transitive recipr. triplets", "3-cycles", "GWESP I -> K -> J", "GWESP I -> K <- J",
            "GWESP I <- K -> J", "GWESP I <- K <- J", "indegree - popularity",
            "indegree - popularity (sqrt)", "outdegree - popularity",
            "outdegree - popularity (sqrt)", "indegree - activity", "indegree - activity (sqrt)",
            "outdegree - activity", "outdegree - activity (sqrt)", "%s ego", "%s alter",
            "%s similarity", "same %s", "%s ego x %s alter", "%s"),
  stringsAsFactors = FALSE
)

# Effects of a stochastic actor-oriented model (man/saom_effects.Rd).
saom_effects <- function(formula, data = NULL) {
  if (is.character(formula)) {
    terms <- lapply(formula, as.name)
  } else if (inherits(formula, "formula")) {
    if (length(formula) != 2L) {
      stop("`formula` must be one-sided, as in `~ density + recip`.", call. = FALSE)
    }
    terms <- split_terms(formula[[2L]])
  } else if (is.data.frame(formula) && all(c("effect", "covariate", "parameter") %in% names(formula))) {
    return(saom_check_effects(formula, data))
  } else {
    stop("`effects` must be a one-sided formula, a character vector, or an effects data frame.",
         call. = FALSE)
  }
  rows <- lapply(terms, parse_effect_term)
  eff <- do.call(rbind, rows)
  saom_check_effects(eff, data)
}

# See man/saom_effects.Rd.
saom_effect_names <- function() SAOM_EFFECTS$effect

split_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("+"))) {
    if (length(expr) == 2L) return(split_terms(expr[[2L]]))
    return(c(split_terms(expr[[2L]]), split_terms(expr[[3L]])))
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("("))) return(split_terms(expr[[2L]]))
  list(expr)
}

parse_effect_term <- function(term) {
  if (is.name(term)) {
    name <- as.character(term)
    args <- list()
  } else if (is.call(term) && is.name(term[[1L]])) {
    name <- as.character(term[[1L]])
    args <- as.list(term)[-1L]
  } else {
    stop(sprintf("Cannot read the effect `%s`.", paste(deparse(term), collapse = "")),
         call. = FALSE)
  }
  row <- SAOM_EFFECTS[SAOM_EFFECTS$effect == name, , drop = FALSE]
  if (!nrow(row)) {
    stop(sprintf("Unknown effect `%s`; see saom_effect_names().", name), call. = FALSE)
  }
  covariate <- NA_character_
  parameter <- NA_real_
  argnames <- names(args)
  if (is.null(argnames)) argnames <- rep("", length(args))
  for (k in seq_along(args)) {
    a <- args[[k]]
    value <- if (is.name(a)) as.character(a) else if (is.character(a) || is.numeric(a)) a else
      stop(sprintf("Cannot read the argument `%s` of effect `%s`.",
                   paste(deparse(a), collapse = ""), name), call. = FALSE)
    if (argnames[k] == "parameter" || (argnames[k] == "" && is.numeric(value))) {
      parameter <- as.double(value)
    } else if (argnames[k] %in% c("", "covariate")) {
      covariate <- as.character(value)
    } else {
      stop(sprintf("Unknown argument `%s` of effect `%s`.", argnames[k], name), call. = FALSE)
    }
  }
  data.frame(effect = name, covariate = covariate, parameter = parameter, stringsAsFactors = FALSE)
}

saom_check_effects <- function(eff, data) {
  eff <- data.frame(effect = as.character(eff$effect), covariate = as.character(eff$covariate),
                    parameter = as.double(eff$parameter), stringsAsFactors = FALSE)
  if (!nrow(eff)) stop("The model needs at least one effect.", call. = FALSE)
  info <- SAOM_EFFECTS[match(eff$effect, SAOM_EFFECTS$effect), , drop = FALSE]
  if (anyNA(info$effect)) {
    stop(sprintf("Unknown effect `%s`; see saom_effect_names().",
                 eff$effect[is.na(info$effect)][1L]), call. = FALSE)
  }
  for (k in seq_len(nrow(eff))) {
    needs <- info$needs[k]
    if (needs %in% c("actor", "dyadic")) {
      if (is.na(eff$covariate[k]) || !nzchar(eff$covariate[k])) {
        stop(sprintf("Effect `%s` needs a covariate, as in `%s(\"name\")`.", eff$effect[k],
                     eff$effect[k]), call. = FALSE)
      }
      if (!is.null(data)) {
        pool <- if (needs == "actor") names(data$covariates) else names(data$dyadic)
        if (!eff$covariate[k] %in% pool) {
          stop(sprintf("Effect `%s` uses the %s covariate `%s`, which the data do not have%s.",
                       eff$effect[k], needs, eff$covariate[k],
                       if (length(pool)) sprintf(" (available: %s)", paste(pool, collapse = ", "))
                       else ""), call. = FALSE)
        }
      }
    } else if (!is.na(eff$covariate[k]) && nzchar(eff$covariate[k])) {
      stop(sprintf("Effect `%s` does not take a covariate.", eff$effect[k]), call. = FALSE)
    }
    if (needs == "parameter") {
      if (is.na(eff$parameter[k])) eff$parameter[k] <- 69
      if (!is.finite(eff$parameter[k]) || eff$parameter[k] < 0) {
        stop(sprintf("Effect `%s` needs a non-negative `parameter` (100 alpha).", eff$effect[k]),
             call. = FALSE)
      }
    } else if (!is.na(eff$parameter[k])) {
      stop(sprintf("Effect `%s` does not take a parameter.", eff$effect[k]), call. = FALSE)
    }
  }
  eff$label <- vapply(seq_len(nrow(eff)), function(k) {
    lab <- info$label[k]
    if (info$needs[k] == "actor" || info$needs[k] == "dyadic") {
      lab <- if (eff$effect[k] == "egoXaltX") sprintf(lab, eff$covariate[k], eff$covariate[k])
      else sprintf(lab, eff$covariate[k])
    } else if (info$needs[k] == "parameter") {
      lab <- sprintf("%s (%s)", lab, format(eff$parameter[k]))
    }
    lab
  }, "")
  eff$name <- vapply(seq_len(nrow(eff)), function(k) {
    if (info$needs[k] %in% c("actor", "dyadic")) paste0(eff$effect[k], ".", eff$covariate[k])
    else if (info$needs[k] == "parameter") paste0(eff$effect[k], ".", format(eff$parameter[k]))
    else eff$effect[k]
  }, "")
  if (anyDuplicated(eff$name)) {
    stop(sprintf("Effect `%s` is included twice.", eff$name[anyDuplicated(eff$name)]),
         call. = FALSE)
  }
  rownames(eff) <- NULL
  eff
}

# The engine's three parallel vectors for an effects table.
saom_spec <- function(eff, data) {
  info <- SAOM_EFFECTS[match(eff$effect, SAOM_EFFECTS$effect), , drop = FALSE]
  cov <- vapply(seq_len(nrow(eff)), function(k) {
    switch(info$needs[k],
      actor = match(eff$covariate[k], names(data$covariates)) - 1L,
      dyadic = match(eff$covariate[k], names(data$dyadic)) - 1L,
      -1L)
  }, 1L)
  par <- ifelse(is.na(eff$parameter), 0, eff$parameter)
  list(codes = as.integer(info$code), covariates = as.integer(cov), parameters = as.double(par))
}

# Target statistics of an effects specification: a matrix with one row
# per period and one column per effect (the statistics of the network at
# the end of the period with the covariates at its start). Used in tests
# and articles; not exported.
saom_targets <- function(data, effects) {
  eff <- saom_effects(effects, data)
  spec <- saom_spec(eff, data)
  periods <- length(data$waves) - 1L
  v <- engine_call("saom_targets",
    .Call(C_saom_targets, saom_handle(data), spec$codes, spec$covariates, spec$parameters),
    edge_rows = FALSE)
  m <- matrix(v, nrow = periods, ncol = nrow(eff), byrow = TRUE)
  colnames(m) <- eff$name
  m
}

# Change statistics of every effect for one actor (1-based) and every
# alter on one wave, with the covariates of one period: a matrix with one
# row per alter and one column per effect. Not exported.
saom_contributions <- function(data, effects, actor, wave = 1L, period = 1L) {
  eff <- saom_effects(effects, data)
  spec <- saom_spec(eff, data)
  v <- engine_call("saom_contributions",
    .Call(C_saom_contributions, saom_handle(data), spec$codes, spec$covariates,
          spec$parameters, data$n, as.integer(wave) - 1L, as.integer(period) - 1L,
          as.integer(actor) - 1L), edge_rows = FALSE)
  m <- matrix(v, nrow = data$n, ncol = nrow(eff), byrow = TRUE)
  dimnames(m) <- list(as.character(data$ids), eff$name)
  m
}

# ---- Estimation ----------------------------------------------------------------

# Settings of the SAOM estimation algorithm (man/saom_control.Rd).
saom_control <- function(conditional = TRUE, nsub = 4, n3 = 1000, n1 = NULL, firstg = 0.2,
                         reduceg = 0.5, diagonalize = 0.2, truncation = 5, theta_bound = 50,
                         seed = NULL, threads = 0L, keep_simulations = TRUE) {
  conditional <- check_flag(conditional, "conditional")
  keep_simulations <- check_flag(keep_simulations, "keep_simulations")
  whole <- function(v, name, min = 0) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < min || v != floor(v)) {
      stop(sprintf("`%s` must be a whole number of at least %d.", name, min), call. = FALSE)
    }
    as.double(v)
  }
  number <- function(v, name) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(sprintf("`%s` must be a positive number.", name), call. = FALSE)
    }
    as.double(v)
  }
  nsub <- whole(nsub, "nsub")
  n3 <- whole(n3, "n3", 2)
  n1 <- if (is.null(n1)) 0 else whole(n1, "n1", 2)
  firstg <- number(firstg, "firstg")
  reduceg <- number(reduceg, "reduceg")
  if (reduceg > 1) stop("`reduceg` must be in (0, 1].", call. = FALSE)
  if (!is.numeric(diagonalize) || length(diagonalize) != 1L || is.na(diagonalize) ||
      diagonalize < 0 || diagonalize > 1) {
    stop("`diagonalize` must be a number in [0, 1].", call. = FALSE)
  }
  truncation <- number(truncation, "truncation")
  theta_bound <- number(theta_bound, "theta_bound")
  threads <- resolve_threads(threads)
  structure(list(conditional = conditional, nsub = nsub, n3 = n3, n1 = n1, firstg = firstg,
                 reduceg = reduceg, diagonalize = as.double(diagonalize), truncation = truncation,
                 theta_bound = theta_bound, seed = seed, threads = threads,
                 keep_simulations = keep_simulations),
            class = "saom_control")
}

# Fit a stochastic actor-oriented model (man/saom_fit.Rd).
saom_fit <- function(data, effects, control = saom_control(), initial = NULL) {
  if (!inherits(data, "saom_data")) stop("`data` must be created by saom_data().", call. = FALSE)
  if (!inherits(control, "saom_control")) {
    stop("`control` must be created by saom_control().", call. = FALSE)
  }
  eff <- saom_effects(effects, data)
  spec <- saom_spec(eff, data)
  periods <- length(data$waves) - 1L
  p <- nrow(eff) + (if (control$conditional) 0L else periods)
  if (!is.null(initial)) {
    if (!is.numeric(initial) || length(initial) != p) {
      stop(sprintf("`initial` must be NULL or a numeric vector of length %d.", p), call. = FALSE)
    }
    initial <- as.double(initial)
    initial[is.na(initial)] <- NaN
  }
  seed <- check_seed(control$seed)
  options <- c(control$conditional, control$n1, control$nsub, control$n3, control$firstg,
               control$reduceg, control$diagonalize, control$truncation, control$theta_bound,
               control$keep_simulations)
  ptr <- saom_handle(data)
  res <- engine_call("saom_fit",
    .Call(C_saom_estimate, ptr, spec$codes, spec$covariates, spec$parameters,
          as.double(options), initial, seed, control$threads), edge_rows = FALSE)
  meta <- res[[18L]]
  n <- data$n
  rate_names <- if (control$conditional) character() else paste0("rate.", seq_len(periods))
  par_names <- c(rate_names, eff$name)
  mat <- function(v, ncol) matrix(v, ncol = ncol, byrow = TRUE)
  vcov <- mat(res[[3L]], p)
  dimnames(vcov) <- list(par_names, par_names)
  derivative <- mat(res[[4L]], p)
  dimnames(derivative) <- list(par_names, par_names)
  stat_cov <- mat(res[[5L]], p)
  dimnames(stat_cov) <- list(par_names, par_names)
  traj <- mat(res[[12L]], p + 1L)
  trajectory <- data.frame(iteration = seq_len(nrow(traj)), subphase = as.integer(traj[, 1L]),
                           traj[, -1L, drop = FALSE])
  names(trajectory) <- c("iteration", "subphase", par_names)
  deviations <- mat(res[[14L]], p)
  colnames(deviations) <- par_names
  times <- mat(res[[15L]], periods)
  colnames(times) <- paste0("period", seq_len(periods))
  aux_len <- as.integer(meta[8L])
  aux <- list(observed = res[[17L]],
              simulated = if (length(res[[16L]])) mat(res[[16L]], aux_len) else NULL, n = n)
  structure(list(
    coefficients = stats::setNames(res[[1L]], par_names),
    se = stats::setNames(res[[2L]], par_names),
    vcov = vcov,
    tconv = stats::setNames(res[[6L]], par_names),
    tconv_max = res[[7L]],
    rates = data.frame(period = seq_len(periods), estimate = res[[10L]], se = res[[11L]]),
    effects = eff,
    targets = stats::setNames(res[[8L]], par_names),
    means = stats::setNames(res[[9L]], par_names),
    derivative = derivative,
    statistic_vcov = stat_cov,
    deviations = deviations,
    phase2_deviations = {
      d <- mat(res[[13L]], p)
      colnames(d) <- par_names
      d
    },
    trajectory = trajectory,
    times = times,
    aux = aux,
    conditional = control$conditional,
    iterations = c(phase1 = meta[4L], phase2 = meta[5L], phase3 = meta[6L]),
    n = n,
    waves = length(data$waves),
    seed = seed,
    control = control,
    call = match.call()
  ), class = "saom_fit")
}

# See man/saom_fit.Rd.
coef.saom_fit <- function(object, ...) object$coefficients

# See man/saom_fit.Rd.
vcov.saom_fit <- function(object, ...) object$vcov

saom_table <- function(x) {
  est <- x$coefficients
  se <- x$se
  tab <- cbind(Estimate = est, `Std. Error` = se, `t ratio` = est / se, `Conv. t` = x$tconv)
  rownames(tab) <- if (x$conditional) x$effects$label else
    c(paste("rate period", seq_len(nrow(x$rates))), x$effects$label)
  tab
}

# See man/saom_fit.Rd.
print.saom_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("\nStochastic actor-oriented model (method of moments, %s estimation)\n",
              if (x$conditional) "conditional" else "unconditional"))
  cat(sprintf("%d waves of a directed network on %d actors; %d phase 3 simulations\n\n",
              x$waves, x$n, as.integer(x$iterations[["phase3"]])))
  if (x$conditional) {
    rates <- cbind(Estimate = x$rates$estimate, `Std. Error` = x$rates$se)
    rownames(rates) <- paste("rate period", x$rates$period)
    print(rates, digits = digits)
    cat("\n")
  }
  print(saom_table(x), digits = digits)
  cat(sprintf("\nOverall maximum convergence ratio: %s\n", format(x$tconv_max, digits = digits)))
  invisible(x)
}

# See man/saom_fit.Rd.
summary.saom_fit <- function(object, ...) {
  converged <- is.finite(object$tconv_max) && object$tconv_max < 0.25 &&
    all(abs(object$tconv) < 0.1)
  structure(c(object, list(converged = converged)), class = c("summary.saom_fit", "saom_fit"))
}

# See man/saom_fit.Rd.
print.summary.saom_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.saom_fit(x, digits = digits)
  cat(if (x$converged) "Convergence is adequate (all |t| < 0.1 and ratio < 0.25).\n" else
    "Convergence is not yet adequate; re-estimate from these estimates (`initial = coef(fit)`).\n")
  cat(sprintf("\nIterations: phase 1 %d, phase 2 %d, phase 3 %d; seed %s\n",
              as.integer(x$iterations[["phase1"]]), as.integer(x$iterations[["phase2"]]),
              as.integer(x$iterations[["phase3"]]), format(x$seed, scientific = FALSE)))
  cat("\nCovariance matrix of the estimates (correlations below the diagonal):\n")
  v <- x$vcov
  show <- v
  sd <- sqrt(diag(v))
  cor <- v / outer(sd, sd)
  show[lower.tri(show)] <- cor[lower.tri(cor)]
  print(show, digits = digits)
  cat("\nObserved and mean simulated statistics:\n")
  print(cbind(observed = x$targets, simulated = x$means), digits = digits)
  invisible(x)
}

# Plot estimation diagnostics of a SAOM fit (man/plot.saom_fit.Rd).
plot.saom_fit <- function(x, which = c("trajectory", "deviations"), ...) {
  which <- match.arg(which)
  if (which == "trajectory") {
    tr <- x$trajectory
    if (!nrow(tr)) stop("The fit has no phase 2 iterations to plot (nsub = 0).", call. = FALSE)
    pars <- names(x$coefficients)
    labels <- rownames(saom_table(x))
    old <- graphics::par(mfrow = grDevices::n2mfrow(length(pars)), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(old))
    bounds <- which(diff(tr$subphase) != 0) + 0.5
    for (k in seq_along(pars)) {
      graphics::plot(tr$iteration, tr[[pars[k]]], type = "l", xlab = "phase 2 iteration",
                     ylab = "value", main = labels[k], ...)
      graphics::abline(v = bounds, lty = 2, col = plot_gray)
      graphics::abline(h = x$coefficients[[k]], col = plot_accent)
    }
  } else {
    sd <- apply(x$deviations, 2L, stats::sd)
    z <- sweep(x$deviations, 2L, ifelse(sd > 0, sd, 1), "/")
    old <- graphics::par(mar = c(4, 10, 2, 1))
    on.exit(graphics::par(old))
    labels <- if (x$conditional) x$effects$label else
      c(paste("rate period", seq_len(nrow(x$rates))), x$effects$label)
    graphics::boxplot(z, horizontal = TRUE, names = labels, las = 1,
                      xlab = "phase 3 deviation (standard deviations)", ...)
    graphics::abline(v = 0, col = plot_accent)
  }
  invisible(x)
}

# ---- Simulation ----------------------------------------------------------------

# Simulate network panels from a stochastic actor-oriented model (man/saom_simulate.Rd).
saom_simulate <- function(data, effects, theta, rates = NULL, n = 1, seed = NULL) {
  if (!inherits(data, "saom_data")) stop("`data` must be created by saom_data().", call. = FALSE)
  eff <- saom_effects(effects, data)
  spec <- saom_spec(eff, data)
  periods <- length(data$waves) - 1L
  if (!is.numeric(theta) || length(theta) != nrow(eff) || anyNA(theta)) {
    stop(sprintf("`theta` must be a numeric vector with one value per effect (%d).", nrow(eff)),
         call. = FALSE)
  }
  conditional <- is.null(rates)
  if (!conditional) {
    if (!is.numeric(rates) || anyNA(rates) || !length(rates) %in% c(1L, periods) || any(rates <= 0)) {
      stop(sprintf("`rates` must be positive, one per period (%d) or a single value.", periods),
           call. = FALSE)
    }
    rates <- as.double(rep_len(rates, periods))
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != floor(n)) {
    stop("`n` must be a positive whole number.", call. = FALSE)
  }
  seed <- check_seed(seed)
  ptr <- saom_handle(data)
  nn <- data$n
  ids <- data$ids
  lapply(seq_len(n), function(s) {
    res <- engine_call("saom_simulate",
      .Call(C_saom_simulate, ptr, spec$codes, spec$covariates, spec$parameters, nn,
            as.double(theta), rates, conditional, seed, as.double(s - 1)), edge_rows = FALSE)
    nets <- lapply(seq_len(periods), function(m) {
      a <- matrix(res[[1L]][(m - 1L) * nn * nn + seq_len(nn * nn)], nn, nn)
      idx <- which(a != 0, arr.ind = TRUE)
      edges <- data.frame(from = ids[idx[, 1L]], to = ids[idx[, 2L]], stringsAsFactors = FALSE)
      netframe(edges, nodes = data$nodes, directed = TRUE, id = data$id)
    })
    stats <- matrix(res[[2L]], ncol = nrow(eff), byrow = TRUE)
    colnames(stats) <- eff$name
    attr(nets, "statistics") <- stats
    attr(nets, "times") <- res[[3L]]
    nets
  })
}

# ---- Goodness of fit -----------------------------------------------------------

# Goodness of fit of a SAOM by auxiliary statistics (man/saom_gof.Rd).
saom_gof <- function(fit, statistic = c("indegree", "outdegree", "triad_census", "geodesic"),
                     levels = NULL, cumulative = TRUE) {
  if (!inherits(fit, "saom_fit")) stop("`fit` must be an saom_fit.", call. = FALSE)
  statistic <- match.arg(statistic)
  cumulative <- check_flag(cumulative, "cumulative")
  if (is.null(fit$aux$simulated)) {
    stop("The fit kept no simulations; estimate with saom_control(keep_simulations = TRUE).",
         call. = FALSE)
  }
  n <- fit$aux$n
  sims <- fit$aux$simulated
  obs <- fit$aux$observed
  pick <- function(cols, names) {
    o <- obs[cols]
    s <- sims[, cols, drop = FALSE]
    names(o) <- names
    colnames(s) <- names
    list(observed = o, simulated = s)
  }
  if (statistic %in% c("indegree", "outdegree")) {
    if (is.null(levels)) levels <- 0:min(8L, n - 1L)
    if (!is.numeric(levels) || anyNA(levels) || any(levels < 0 | levels > n - 1L)) {
      stop(sprintf("`levels` must be degrees between 0 and %d.", n - 1L), call. = FALSE)
    }
    offset <- if (statistic == "indegree") 0L else n
    counts <- pick(offset + seq_len(n), 0:(n - 1L))
    values <- if (cumulative) {
      list(observed = cumsum(counts$observed)[levels + 1L],
           simulated = t(apply(counts$simulated, 1L, cumsum))[, levels + 1L, drop = FALSE])
    } else {
      list(observed = counts$observed[levels + 1L],
           simulated = counts$simulated[, levels + 1L, drop = FALSE])
    }
    labels <- if (cumulative) paste0("<= ", levels) else as.character(levels)
  } else if (statistic == "triad_census") {
    values <- pick(2L * n + 1:16, c("003", "012", "102", "021D", "021U", "021C", "111D",
                                    "111U", "030T", "030C", "201", "120D", "120U", "120C",
                                    "210", "300"))
    labels <- names(values$observed)
  } else {
    if (is.null(levels)) levels <- c(seq_len(min(5L, n - 1L)), Inf)
    if (!is.numeric(levels) || anyNA(levels) || any(levels < 1)) {
      stop("`levels` must be distances of at least 1 (Inf for unreachable pairs).", call. = FALSE)
    }
    counts <- pick(2L * n + 16L + seq_len(n), c(seq_len(n - 1L), Inf))
    idx <- ifelse(is.infinite(levels), n, pmin(levels, n))
    values <- if (cumulative) {
      list(observed = cumsum(counts$observed)[idx],
           simulated = t(apply(counts$simulated, 1L, cumsum))[, idx, drop = FALSE])
    } else {
      list(observed = counts$observed[idx], simulated = counts$simulated[, idx, drop = FALSE])
    }
    labels <- if (cumulative) paste0("<= ", levels) else as.character(levels)
  }
  observed <- stats::setNames(as.double(values$observed), labels)
  simulated <- matrix(values$simulated, nrow = nrow(sims))
  colnames(simulated) <- labels
  centre <- colMeans(simulated)
  cov <- stats::cov(simulated)
  ginv <- generalized_inverse(cov)
  mhd <- function(v) {
    d <- v - centre
    as.numeric(t(d) %*% ginv %*% d)
  }
  distance <- mhd(observed)
  distances <- apply(simulated, 1L, mhd)
  structure(list(statistic = statistic, observed = observed, simulated = simulated,
                 cumulative = cumulative, distance = distance, distances = distances,
                 p = mean(distances >= distance)),
            class = "saom_gof")
}

# Moore-Penrose inverse of a symmetric matrix, by singular value
# decomposition; singular values below tol times the largest are dropped.
generalized_inverse <- function(m, tol = sqrt(.Machine$double.eps)) {
  if (!all(is.finite(m))) return(matrix(0, nrow(m), ncol(m)))
  s <- svd(m)
  keep <- s$d > tol * max(s$d, 0)
  if (!any(keep)) return(matrix(0, nrow(m), ncol(m)))
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# See man/saom_gof.Rd.
print.saom_gof <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  what <- switch(x$statistic, indegree = "in-degree distribution",
                 outdegree = "out-degree distribution", triad_census = "triad census",
                 geodesic = "geodesic distribution")
  cat(sprintf("\nGoodness of fit: %s%s\n", what, if (x$cumulative && x$statistic != "triad_census")
    " (cumulative)" else ""))
  cat(sprintf("Mahalanobis distance %s; Monte Carlo p-value %s (%d simulations)\n\n",
              format(x$distance, digits = digits), format(x$p, digits = digits),
              nrow(x$simulated)))
  q <- apply(x$simulated, 2L, stats::quantile, probs = c(0.05, 0.5, 0.95))
  print(rbind(observed = x$observed, `simulated 5%` = q[1L, ], `simulated median` = q[2L, ],
              `simulated 95%` = q[3L, ]), digits = digits)
  invisible(x)
}

# See man/saom_gof.Rd.
plot.saom_gof <- function(x, violin = TRUE, ...) {
  violin <- check_flag(violin, "violin")
  k <- ncol(x$simulated)
  what <- switch(x$statistic, indegree = "in-degree", outdegree = "out-degree",
                 triad_census = "triad type", geodesic = "distance")
  ylim <- range(c(x$simulated, x$observed), na.rm = TRUE)
  graphics::plot(NA, xlim = c(0.5, k + 0.5), ylim = ylim, xaxt = "n", xlab = what,
                 ylab = "count",
                 main = sprintf("Goodness of fit: p = %s", format(x$p, digits = 3)), ...)
  graphics::axis(1, at = seq_len(k), labels = colnames(x$simulated), las = if (k > 8) 2 else 1)
  for (j in seq_len(k)) {
    v <- x$simulated[, j]
    if (violin && length(unique(v)) > 1L) {
      d <- stats::density(v, from = min(v), to = max(v))
      w <- 0.4 * d$y / max(d$y)
      graphics::polygon(c(j - w, rev(j + w)), c(d$x, rev(d$x)), col = "gray85", border = "gray40")
    } else {
      graphics::boxplot(v, at = j, add = TRUE, axes = FALSE, col = "gray85", outline = FALSE)
    }
  }
  graphics::lines(seq_len(k), x$observed, col = plot_accent, lwd = 2)
  graphics::points(seq_len(k), x$observed, col = plot_accent, pch = 19)
  invisible(x)
}
