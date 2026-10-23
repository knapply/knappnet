#' Network panel data for stochastic actor-oriented models
#'
#' Collects two or more observations of a directed network on the same
#' actors, with actor and dyadic covariates, into the data object that
#' [saom_fit()] estimates stochastic actor-oriented models from (Snijders
#' 2001; Snijders, van de Bunt, and Steglich 2010).
#'
#' @param networks A list of at least two directed [netframe()]s observed at
#'   successive moments (for example the output of [snapshots()]), all on
#'   the same node identifiers. Nodes are matched by identifier and ordered
#'   as in the first network. Each network is read as a binary adjacency
#'   matrix: parallel edges collapse to one tie, self-loops are dropped,
#'   and weights are ignored.
#' @param covariates Actor covariates: `NULL`, a character vector naming
#'   columns of the node tables (a column whose values are the same in every
#'   wave is a constant covariate; one that changes is a changing covariate
#'   whose value during a period is its value at the start of that period),
#'   or a named list of numeric vectors (one value per actor, in the order
#'   of the first network's nodes) and matrices with one row per actor and
#'   one column per period (one fewer than the waves; a column for the last
#'   wave may be given and is ignored, as in RSiena). Changing covariates
#'   need at least three waves. Values must not be missing.
#' @param dyadic Dyadic covariates: `NULL` or a named list of square
#'   numeric matrices (rows and columns in the order of the first
#'   network's nodes; the diagonal is ignored) or netframes on the same
#'   nodes, read as adjacency matrices of their weights (or of ties when
#'   unweighted).
#'
#' @details
#' Covariates are centered as RSiena centers them (Ripley et al. 2026,
#' section 5.3.2): actor covariates by their mean over all actors and all
#' periods, dyadic covariates by the mean of their off-diagonal cells, and
#' the similarity scores of `simX` by their mean over all ordered pairs of
#' distinct actors and periods. `sameX` is not centered. For changing
#' covariates RSiena computes the similarity mean differently (its
#' definition is not documented), so `simX` of a changing covariate can
#' differ slightly from RSiena's.
#'
#' Missing tie or covariate values, structural zeros and ones, and changes
#' of composition (actors joining or leaving) are not supported; see
#' "Not implemented" in [saom_fit()].
#'
#' @return An object of class `saom_data`: a list with `n` (actors),
#'   `waves` (a list of logical adjacency matrices), `nodes` (the node table
#'   of the first network), `covariates`, and `dyadic`.
#' @seealso [saom_fit()], [saom_simulate()], [snapshots()].
#' @references
#' Snijders, T. A. B. (2001). The statistical evaluation of social network
#' dynamics. *Sociological Methodology*, 31(1), 361-395.
#' \doi{10.1111/0081-1750.00099}
#'
#' Snijders, T. A. B., van de Bunt, G. G., and Steglich, C. E. G. (2010).
#' Introduction to stochastic actor-based models for network dynamics.
#' *Social Networks*, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}
#'
#' Ripley, R. M., Snijders, T. A. B., Boda, Z., and Voros, A. (2026).
#' *Manual for RSiena*. University of Oxford, Department of Statistics and
#' Nuffield College; University of Groningen, Department of Sociology.
#' @examples
#' w1 <- netframe(data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "a")),
#'                nodes = data.frame(name = c("a", "b", "c", "d"), smoke = c(1, 0, 1, 0)))
#' w2 <- netframe(data.frame(from = c("a", "b", "c", "d"), to = c("b", "c", "a", "a")),
#'                nodes = w1$nodes)
#' saom_data(list(w1, w2), covariates = "smoke")
#' @export
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

#' @export
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

#' Effects of a stochastic actor-oriented model
#'
#' Turns a model formula into the table of effects [saom_fit()] estimates,
#' and lists the effects available. Effects carry the short names of RSiena
#' (Ripley et al. 2026, chapter 13), and their definitions match RSiena's
#' exactly; see Details.
#'
#' @param formula A one-sided formula whose terms are effect names, for
#'   example `~ density + recip + transTrip + simX("smoke")`. Structural
#'   effects are bare names; covariate effects are calls with the covariate
#'   name as a string or symbol; the GWESP effects accept an internal
#'   parameter, `gwespFF(parameter = 69)` (the default, giving
#'   \eqn{\alpha = 0.69 = \log 2}). A character vector of structural effect
#'   names is also accepted.
#' @param data Optional [saom_data()] whose covariates the formula may use;
#'   covariate names are then checked.
#'
#' @details
#' With \eqn{x_{ij} = 1} for a tie from i to j, \eqn{x_{i+}} the out-degree,
#' \eqn{x_{+i}} the in-degree, \eqn{v_i} a centered actor covariate and
#' \eqn{w_{ij}} a centered dyadic covariate, the evaluation function of
#' actor i is \eqn{f_i(x) = \sum_k \beta_k s_{ik}(x)} with these statistics:
#' \tabular{lll}{
#'   \code{density} \tab out-degree \tab \eqn{x_{i+}} \cr
#'   \code{recip} \tab reciprocity \tab \eqn{\sum_j x_{ij} x_{ji}} \cr
#'   \code{transTrip} \tab transitive triplets \tab \eqn{\sum_{j,h} x_{ij} x_{ih} x_{hj}} \cr
#'   \code{transTies} \tab transitive ties \tab \eqn{\sum_j x_{ij} \max_h x_{ih} x_{hj}} \cr
#'   \code{transRecTrip} \tab transitive reciprocated triplets \tab \eqn{\sum_{j,h} x_{ij} x_{ji} x_{ih} x_{hj}} \cr
#'   \code{cycle3} \tab 3-cycles \tab \eqn{\sum_{j,h} x_{ij} x_{jh} x_{hi}} \cr
#'   \code{gwespFF} etc. \tab geometrically weighted edgewise shared partners \tab
#'     \eqn{\sum_j x_{ij} e^\alpha \{1 - (1 - e^{-\alpha})^{esp_{ij}}\}} \cr
#'   \code{inPop}, \code{inPopSqrt} \tab in-degree popularity \tab \eqn{\sum_j x_{ij} x_{+j}}, \eqn{\sum_j x_{ij} \sqrt{x_{+j}}} \cr
#'   \code{outPop}, \code{outPopSqrt} \tab out-degree popularity \tab \eqn{\sum_j x_{ij} x_{j+}}, \eqn{\sum_j x_{ij} \sqrt{x_{j+}}} \cr
#'   \code{inAct}, \code{inActSqrt} \tab in-degree activity \tab \eqn{x_{i+} x_{+i}}, \eqn{x_{i+} \sqrt{x_{+i}}} \cr
#'   \code{outAct}, \code{outActSqrt} \tab out-degree activity \tab \eqn{x_{i+}^2}, \eqn{x_{i+}^{1.5}} \cr
#'   \code{egoX} \tab covariate ego \tab \eqn{v_i x_{i+}} \cr
#'   \code{altX} \tab covariate alter \tab \eqn{\sum_j x_{ij} v_j} \cr
#'   \code{simX} \tab covariate similarity \tab \eqn{\sum_j x_{ij} (sim_{ij} - \widehat{sim})} \cr
#'   \code{sameX} \tab same covariate \tab \eqn{\sum_j x_{ij} I\{v_i = v_j\}} \cr
#'   \code{egoXaltX} \tab covariate ego x alter \tab \eqn{v_i \sum_j x_{ij} v_j} \cr
#'   \code{X} \tab dyadic covariate \tab \eqn{\sum_j x_{ij} w_{ij}}
#' }
#' For the GWESP effects \eqn{esp_{ij}} counts the actors h with
#' \eqn{i \to h \to j} (`gwespFF`), \eqn{i \to h \leftarrow j} (`gwespFB`),
#' \eqn{i \leftarrow h \to j} (`gwespBF`), or \eqn{i \leftarrow h \leftarrow j}
#' (`gwespBB`), and \eqn{\alpha} is `parameter / 100` (Hunter 2007; Ripley
#' et al. 2026). As in RSiena these are elementary effects: the change
#' statistic of a tie is its own weight, not the change in the sum. The
#' similarity is \eqn{sim_{ij} = 1 - |v_i - v_j| / range} with the range
#' over all actors and waves, centered by its mean over all ordered pairs of
#' distinct actors (and waves).
#'
#' @return `saom_effects()`: a data frame with one row per effect: `effect`
#'   (short name), `covariate`, `parameter`, `label`, and `name` (the
#'   coefficient name). `saom_effect_names()`: the short names of all
#'   available effects.
#' @references
#' Ripley, R. M., Snijders, T. A. B., Boda, Z., and Voros, A. (2026).
#' *Manual for RSiena*. University of Oxford, Department of Statistics and
#' Nuffield College; University of Groningen, Department of Sociology.
#'
#' Hunter, D. R. (2007). Curved exponential family models for social
#' networks. *Social Networks*, 29(2), 216-230.
#' \doi{10.1016/j.socnet.2006.08.005}
#' @examples
#' saom_effects(~ density + recip + transTrip + gwespFF(parameter = 69) + simX("smoke"))
#' saom_effect_names()
#' @export
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

#' @rdname saom_effects
#' @export
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

#' Settings of the SAOM estimation algorithm
#'
#' @param conditional Whether to estimate conditionally on the observed
#'   number of changes per period (the default, as in RSiena): each period's
#'   simulation runs until as many tie variables differ from the start of
#'   the period as were observed to differ, and the rate parameters are
#'   estimated afterwards from the simulated period lengths. `FALSE`
#'   simulates one time unit per period and estimates the rates as
#'   parameters of the Robbins-Monro iteration.
#' @param nsub Number of subphases in phase 2.
#' @param n3 Number of simulations in phase 3, which give the standard
#'   errors, convergence diagnostics, and goodness-of-fit statistics.
#'   1000 is adequate for exploration; use more for publication.
#' @param n1 Number of phase 1 simulations; `NULL` for \eqn{7 + 3p}.
#' @param firstg Initial gain of the Robbins-Monro steps.
#' @param reduceg Factor by which the gain shrinks after each subphase.
#' @param diagonalize Weight in \eqn{[0, 1]} of the diagonal of the
#'   derivative matrix in the step matrix (RSiena's default 0.2; 1 uses
#'   only the diagonal).
#' @param truncation Phase 2 deviations larger than this many standard
#'   deviations (estimated in phase 1) are truncated.
#' @param theta_bound Estimation stops with an error when a parameter
#'   exceeds this in absolute value, a sign of divergence.
#' @param seed Seed for the simulations, a whole number in \[0, 2^53).
#'   `NULL` draws one from R's random number generator, so [set.seed()]
#'   makes results reproducible.
#' @param threads Number of threads for phases 1 and 3; `0` uses all
#'   available cores (see "Threads" in [netframe()]). Results do not depend
#'   on the thread count except in the last bits of floating-point sums.
#' @param keep_simulations Whether to keep the auxiliary statistics of every
#'   phase 3 simulation for [saom_gof()] (\eqn{n_3 (3n + 16)} numbers).
#' @return A list of class `saom_control`.
#' @seealso [saom_fit()]
#' @export
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

#' Fit a stochastic actor-oriented model
#'
#' Estimates a stochastic actor-oriented model (SAOM; Snijders 2001;
#' Snijders, van de Bunt, and Steglich 2010; Snijders 2017) of the changes
#' between the waves of a network panel by the method of moments, using the
#' three-phase Robbins-Monro stochastic approximation algorithm of RSiena.
#'
#' @param data An [saom_data()] object.
#' @param effects The effects of the evaluation function: a one-sided
#'   formula such as `~ density + recip + transTrip + simX("smoke")`, a
#'   character vector of structural effects, or the data frame returned by
#'   [saom_effects()], which documents every effect and its definition.
#' @param control Algorithm settings from [saom_control()].
#' @param initial Initial parameter values: `NULL` for the standard initial
#'   values, or a numeric vector with one value per effect (rates first,
#'   one per period, when `control$conditional` is `FALSE`); `NA` keeps the
#'   standard value for that parameter.
#'
#' @details
#' The model is a continuous-time Markov chain of *ministeps* (Snijders
#' 2001): actor \eqn{i} gets opportunities to change at rate \eqn{\rho_m}
#' during period \eqn{m}, and at each opportunity adds or drops one
#' outgoing tie, or changes nothing, choosing the resulting network
#' \eqn{x'} with probability proportional to \eqn{\exp f_i(x')}, where the
#' evaluation function \eqn{f_i(x) = \sum_k \beta_k s_{ik}(x)} is a weighted
#' sum of the effects in [saom_effects()]. The method of moments chooses
#' \eqn{\beta} so that the statistics \eqn{\sum_m \sum_i s_{ik}} of networks
#' simulated from each wave to the next equal, on average, the observed
#' statistics of the next wave (with changing covariates at their values at
#' the start of each period).
#'
#' Estimation follows RSiena's algorithm (Ripley et al. 2026, chapter 7 and
#' the document *Siena algorithms*): phase 1 estimates the derivative of
#' the expected statistics with respect to the parameters by the
#' score-function method (Schweinberger and Snijders 2007) from
#' \eqn{7 + 3p} simulations; phase 2 runs `nsub` subphases of
#' Robbins-Monro steps \eqn{\theta \leftarrow \theta - a D^{-1} (S - s)},
#' each subphase ending when the deviations have crossed zero and setting
#' \eqn{\theta} to its average, with the gain halved between subphases;
#' phase 3 simulates `n3` times at the final \eqn{\theta} to obtain the
#' convergence t-ratios \eqn{(\bar S_k - s_k) / sd(S_k)} (which should be
#' below 0.1 in absolute value), the overall maximum convergence ratio
#' \eqn{\sqrt{(\bar S - s)' \Sigma^{-1} (\bar S - s)}} (which should be below
#' 0.25), and the standard errors by the delta method,
#' \eqn{cov(\hat\theta) = D^{-1} \Sigma D^{-T}}. Conditional estimation (the
#' default) simulates each period until the observed number of tie changes
#' is reached and estimates the rate of each period as the mean simulated
#' period length in phase 3, reporting the standard deviation of those
#' lengths as its standard error (Snijders 2001, section 4.2; RSiena may
#' compute the latter differently).
#'
#' Standard initial values follow Snijders (2005, equation 11.41): the
#' out-degree parameter starts at \eqn{\log(c^+ d / (c^- (1 - d))) / 2},
#' with \eqn{c^+} ties created, \eqn{c^-} ties dissolved, and \eqn{d} the
#' mean density; other effects start at 0 and rates at twice the observed
#' changes per actor. The simulations use the seeded, counter-based
#' generator of the engine, so a seed reproduces the estimate exactly for
#' any number of threads.
#'
#' @section Not implemented:
#' Behavior dynamics and co-evolution, several dependent networks,
#' non-directed networks, creation and endowment effects, rate effects
#' beyond a constant rate per period, missing data, structural zeros and
#' ones, composition change, maximum likelihood and Bayesian estimation,
#' score-type tests, and time heterogeneity tests. For these, see RSiena.
#'
#' @return An object of class `saom_fit` with `coefficients`, `se`, `vcov`,
#'   `tconv` (convergence t-ratios), `tconv_max`, `rates` (a data frame
#'   with the rate parameter and standard error of each period), `effects`,
#'   `targets`, `means` (phase 3 means of the statistics), `deviations`
#'   (phase 3 deviations, one row per simulation), `trajectory` (phase 2
#'   parameter values by iteration), `times` (phase 3 simulated period
#'   lengths), `aux` (auxiliary statistics for [saom_gof()]), `iterations`,
#'   `seed`, `control`, and `call`. Methods: `print`, `summary`, `coef`,
#'   `vcov`, `plot`.
#' @seealso [saom_data()], [saom_effects()], [saom_control()],
#'   [saom_gof()], [saom_simulate()], [plot.saom_fit()].
#' @references
#' Snijders, T. A. B. (2001). The statistical evaluation of social network
#' dynamics. *Sociological Methodology*, 31(1), 361-395.
#' \doi{10.1111/0081-1750.00099}
#'
#' Snijders, T. A. B. (2005). Models for longitudinal network data. In
#' P. J. Carrington, J. Scott, and S. Wasserman (Eds.), *Models and Methods
#' in Social Network Analysis* (pp. 215-247). Cambridge University Press.
#' \doi{10.1017/CBO9780511811395.011}
#'
#' Snijders, T. A. B., van de Bunt, G. G., and Steglich, C. E. G. (2010).
#' Introduction to stochastic actor-based models for network dynamics.
#' *Social Networks*, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}
#'
#' Snijders, T. A. B. (2017). Stochastic actor-oriented models for network
#' dynamics. *Annual Review of Statistics and Its Application*, 4, 343-363.
#' \doi{10.1146/annurev-statistics-060116-054035}
#'
#' Schweinberger, M., and Snijders, T. A. B. (2007). Markov models for
#' digraph panel data: Monte Carlo-based derivative estimation.
#' *Computational Statistics & Data Analysis*, 51(9), 4465-4483.
#' \doi{10.1016/j.csda.2006.07.014}
#'
#' Robbins, H., and Monro, S. (1951). A stochastic approximation method.
#' *The Annals of Mathematical Statistics*, 22(3), 400-407.
#' \doi{10.1214/aoms/1177729586}
#'
#' Ripley, R. M., Snijders, T. A. B., Boda, Z., and Voros, A. (2026).
#' *Manual for RSiena*. University of Oxford, Department of Statistics and
#' Nuffield College; University of Groningen, Department of Sociology.
#' @examples
#' # A small panel simulated from known parameters, then re-estimated.
#' set.seed(1)
#' start <- netframe(data.frame(from = sample(20, 40, TRUE), to = sample(20, 40, TRUE)),
#'                   nodes = data.frame(id = 1:20, smoke = rep(0:1, 10)))
#' panel <- saom_simulate(saom_data(list(start, start), covariates = "smoke"),
#'                        ~ density + recip + simX("smoke"), theta = c(-2, 1.5, 0.5),
#'                        rates = 4, seed = 1)[[1]]
#' data <- saom_data(c(list(start), panel), covariates = "smoke")
#' fit <- saom_fit(data, ~ density + recip + simX("smoke"),
#'                 control = saom_control(n3 = 200, nsub = 2, seed = 2, threads = 2))
#' fit
#' summary(fit)
#' coef(fit)
#' # Unconditional estimation treats the rates as ordinary parameters.
#' saom_fit(data, ~ density + recip,
#'          control = saom_control(conditional = FALSE, n3 = 200, nsub = 2, seed = 3,
#'                                 threads = 2))
#' @export
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
  seed <- qap_seed(control$seed)
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

#' @export
coef.saom_fit <- function(object, ...) object$coefficients

#' @export
vcov.saom_fit <- function(object, ...) object$vcov

saom_table <- function(x) {
  est <- x$coefficients
  se <- x$se
  tab <- cbind(Estimate = est, `Std. Error` = se, `t ratio` = est / se, `Conv. t` = x$tconv)
  rownames(tab) <- if (x$conditional) x$effects$label else
    c(paste("rate period", seq_len(nrow(x$rates))), x$effects$label)
  tab
}

#' @export
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

#' @export
summary.saom_fit <- function(object, ...) {
  converged <- is.finite(object$tconv_max) && object$tconv_max < 0.25 &&
    all(abs(object$tconv) < 0.1)
  structure(c(object, list(converged = converged)), class = c("summary.saom_fit", "saom_fit"))
}

#' @export
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

#' Plot estimation diagnostics of a SAOM fit
#'
#' `which = "trajectory"` draws the parameter values after every phase 2
#' iteration, with the subphases separated by dashed lines and the final
#' estimate as a horizontal line: a well-behaved estimation settles down
#' within each subphase. `which = "deviations"` draws boxplots of the
#' phase 3 deviations of the simulated statistics from their targets,
#' divided by their standard deviations; they should be centered on 0.
#'
#' @param x An [saom_fit()].
#' @param which `"trajectory"` or `"deviations"`.
#' @param ... Passed to the plotting functions.
#' @return `x`, invisibly.
#' @examples
#' set.seed(2)
#' start <- netframe(data.frame(from = sample(15, 30, TRUE), to = sample(15, 30, TRUE)),
#'                   nodes = data.frame(id = 1:15))
#' panel <- saom_simulate(saom_data(list(start, start)), ~ density + recip,
#'                        theta = c(-1.5, 1), rates = 3, seed = 1)[[1]]
#' fit <- saom_fit(saom_data(c(list(start), panel)), ~ density + recip,
#'                 control = saom_control(n3 = 100, nsub = 2, seed = 2, threads = 2))
#' plot(fit)
#' plot(fit, which = "deviations")
#' @export
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
      graphics::abline(v = bounds, lty = 2, col = "grey60")
      graphics::abline(h = x$coefficients[[k]], col = "firebrick")
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
    graphics::abline(v = 0, col = "firebrick")
  }
  invisible(x)
}

# ---- Simulation ----------------------------------------------------------------

#' Simulate network panels from a stochastic actor-oriented model
#'
#' Simulates the change from each wave of `data` to the next under given
#' parameters, as [saom_fit()] does internally: conditional on the observed
#' number of changes (`rates = NULL`) or for one time unit per period with
#' the given rates.
#'
#' @inheritParams saom_fit
#' @param theta Effect parameters, one per effect, in the order of
#'   `saom_effects(effects)`.
#' @param rates Rate parameters, one per period (or one value recycled), or
#'   `NULL` to simulate until the observed number of changes is reached.
#' @param n Number of panels to simulate.
#' @param seed Seed, a whole number in \[0, 2^53); `NULL` draws one from
#'   R's random number generator.
#' @return A list of `n` simulated panels. Each panel is a list with one
#'   [netframe()] per period, the network at the end of the period, on the
#'   nodes of `data`; the attribute `"statistics"` holds the effect
#'   statistics of each period's network and `"times"` the simulated period
#'   lengths.
#' @seealso [saom_fit()], [saom_gof()]
#' @examples
#' w1 <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "a")),
#'                nodes = data.frame(name = c("a", "b", "c", "d")))
#' sims <- saom_simulate(saom_data(list(w1, w1)), ~ density + recip, theta = c(-1, 1),
#'                       rates = 2, n = 2, seed = 7)
#' sims[[1]][[1]]
#' attr(sims[[1]], "statistics")
#' @export
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
  seed <- qap_seed(seed)
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

#' Goodness of fit of a SAOM by auxiliary statistics
#'
#' Compares a statistic of the observed networks (summed over the waves
#' after the first) with its distribution over the networks simulated in
#' phase 3 of the estimation, as proposed by Lospinoso and Snijders (2019)
#' and implemented in RSiena's `sienaGOF()`. The discrepancy is the
#' Mahalanobis distance of the observed statistic from the mean of the
#' simulated ones, using a generalized inverse of their covariance matrix,
#' and the p-value is the proportion of simulations whose own distance is
#' at least as large.
#'
#' @param fit An [saom_fit()] estimated with `keep_simulations = TRUE`.
#' @param statistic Which auxiliary statistic: the in-degree or out-degree
#'   distribution, the triad census (Holland and Leinhardt 1976), or the
#'   distribution of geodesic distances between ordered pairs of actors.
#' @param levels For the degree distributions, the degrees to report
#'   (default 0 to 8); for geodesics, the distances (default 1 to 5 and
#'   `Inf` for unreachable pairs). Ignored for the triad census.
#' @param cumulative Whether degree and geodesic counts are cumulative
#'   (actors with degree at most k; pairs at distance at most d), as in
#'   RSiena.
#' @return An object of class `saom_gof` with `observed`, `simulated` (one
#'   row per simulation), `p` (the Monte Carlo p-value), `distance`
#'   (the observed Mahalanobis distance), and `distances` (one per
#'   simulation). Methods: `print`, `plot`.
#' @references
#' Lospinoso, J., and Snijders, T. A. B. (2019). Goodness of fit for
#' stochastic actor-oriented models. *Methodological Innovations*, 12(3).
#' \doi{10.1177/2059799119884282}
#'
#' Holland, P. W., and Leinhardt, S. (1976). Local structure in social
#' networks. *Sociological Methodology*, 7, 1-45. \doi{10.2307/270703}
#' @examples
#' set.seed(3)
#' start <- netframe(data.frame(from = sample(15, 30, TRUE), to = sample(15, 30, TRUE)),
#'                   nodes = data.frame(id = 1:15))
#' panel <- saom_simulate(saom_data(list(start, start)), ~ density + recip,
#'                        theta = c(-1.5, 1), rates = 3, seed = 1)[[1]]
#' fit <- saom_fit(saom_data(c(list(start), panel)), ~ density + recip,
#'                 control = saom_control(n3 = 100, nsub = 2, seed = 2, threads = 2))
#' gof <- saom_gof(fit, "indegree")
#' gof
#' plot(gof)
#' plot(saom_gof(fit, "triad_census"))
#' @export
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

#' @export
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

#' @rdname saom_gof
#' @param x An `saom_gof` object.
#' @param violin Whether to draw the simulated distributions as violins
#'   (kernel densities) rather than boxplots.
#' @param ... Passed to the plotting functions.
#' @export
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
      graphics::polygon(c(j - w, rev(j + w)), c(d$x, rev(d$x)), col = "grey85", border = "grey40")
    } else {
      graphics::boxplot(v, at = j, add = TRUE, axes = FALSE, col = "grey85", outline = FALSE)
    }
  }
  graphics::lines(seq_len(k), x$observed, col = "firebrick", lwd = 2)
  graphics::points(seq_len(k), x$observed, col = "firebrick", pch = 19)
  invisible(x)
}
