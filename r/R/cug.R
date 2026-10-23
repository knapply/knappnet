#' Conditional uniform graph tests
#'
#' Compares a network statistic with its distribution over random networks
#' that share some feature of the observed network (Anderson, Butts, and
#' Carley 1999), as sna's `cug.test()` does. The null model is uniform over
#' networks with the same number of nodes (`conditioning = "size"`, that is
#' \eqn{G(n, 1/2)}), the same number of edges (`"edges"`, \eqn{G(n, m)}),
#' or the same dyad census (`"dyad.census"`, the U|MAN distribution of
#' Holland and Leinhardt 1976: the same numbers of mutual, asymmetric, and
#' null dyads).
#'
#' @details
#' `statistic` is any function that takes a netframe (and `...`) and returns
#' a single number, such as [transitivity()], [centralization()], or
#' `function(g) mean(degree(g))`. It is applied to the observed network as
#' given and to `reps` random networks on `1, 2, ..., n` nodes, drawn in
#' parallel with the generators of [random_gnp()], [random_gnm()], and
#' [random_uman()] from a seeded counter-based generator, so the null
#' distribution for a given `seed` does not depend on `threads`; the
#' statistic itself runs in R, one network at a time.
#'
#' The random networks are simple and directed if and only if the observed
#' network is. The observed network's number of edges and dyad census are
#' counted on its simple reading (parallel edges once, self-loops ignored),
#' as sna does on an adjacency matrix with `diag = FALSE`; a network with
#' self-loops or parallel edges gets a warning, since the statistic sees
#' them while the null networks never have them. In an undirected network
#' every tie is a mutual dyad, so conditioning on the dyad census is the
#' same as conditioning on the number of edges and is treated so. Weights
#' and attributes are not carried to the random networks.
#'
#' The p-values are the proportions of random networks whose statistic is
#' at least (`p_greater`) and at most (`p_less`) the observed one, as sna
#' reports them (`Pr(X >= Obs)` and `Pr(X <= Obs)`). A statistic that is
#' `NA` or `NaN` on some random networks is dropped from the distribution
#' for those networks, and the number dropped is reported.
#'
#' @inheritParams betweenness
#' @param statistic A function of a netframe returning one number.
#' @param conditioning What the random networks share with the observed one:
#'   `"size"`, `"edges"`, or `"dyad.census"`.
#' @param reps Number of random networks.
#' @param ... Further arguments passed to `statistic`.
#' @param seed Seed for the random networks, a whole number between 0 and
#'   \eqn{2^{53}}; `NULL` draws one from R's random number generator, so
#'   [set.seed()] makes results reproducible.
#' @param object A `cug_test`.
#'
#' @references
#' Anderson, B. S., Butts, C., & Carley, K. (1999). The interaction of size
#' and density with graph-level indices. *Social Networks*, 21(3), 239-267.
#' \doi{10.1016/S0378-8733(99)00011-8}
#'
#' Holland, P. W., & Leinhardt, S. (1976). Local structure in social
#' networks. *Sociological Methodology*, 7, 1-45. \doi{10.2307/270703}
#'
#' @return An object of class `cug_test`: a list with `observed` (the
#'   statistic on `x`), `null` (its values on the random networks),
#'   `p_greater` and `p_less`, `conditioning`, `reps`, `dropped` (random
#'   networks on which the statistic was undefined), and `statistic` (the
#'   expression given). `summary()` adds quantiles of the null distribution.
#' @export
cug_test <- function(x, statistic, conditioning = c("size", "edges", "dyad.census"),
                     reps = 1000L, ..., seed = NULL, threads = 0L) {
  check_netframe(x)
  if (!is.function(statistic)) {
    stop("`statistic` must be a function of a netframe returning one number.", call. = FALSE)
  }
  conditioning <- match.arg(conditioning)
  reps <- check_count(reps, "reps")
  if (reps < 1L) stop("`reps` must be at least 1.", call. = FALSE)
  seed <- groups_seed(seed)
  threads <- resolve_threads(threads)
  label <- deparse(substitute(statistic))[[1L]]

  n <- nrow(x$nodes)
  s <- edge_summary(x)
  if (s$loops || s$parallel) {
    warning("The network has self-loops or parallel edges; the random networks are simple, ",
            "and the observed number of edges and dyad census ignore them.", call. = FALSE)
  }
  census <- simple_dyads(x)
  observed <- as_statistic(statistic(x, ...), "the observed network")

  # Every tie of an undirected network is a mutual dyad, so its dyad census
  # is its edge count.
  conditioning_used <- if (conditioning == "dyad.census" && !x$directed) "edges" else conditioning
  draw <- switch(conditioning_used,
    size = function(k, s) .Call(C_random_gnp, n, 0.5, x$directed, FALSE, k, s, threads),
    edges = function(k, s) {
      .Call(C_random_gnm, n, census$edges, x$directed, FALSE, k, s, threads)
    },
    dyad.census = function(k, s) {
      .Call(C_random_uman, n, census$mutual, census$asym, k, s, threads)
    }
  )
  nodes <- data.frame(id = seq_len(n))
  # All replicates are drawn in one engine call: replicate r uses random
  # stream (seed, r), so the null distribution for a seed is the same on
  # any number of threads.
  lists <- engine_call("cug_test", draw(reps, seed), edge_rows = FALSE)
  null <- vapply(lists, function(e) {
    g <- netframe(data.frame(from = e[[1L]], to = e[[2L]]), nodes, directed = x$directed)
    as_statistic(statistic(g, ...), "a random network")
  }, numeric(1))
  dropped <- sum(is.na(null))
  null <- null[!is.na(null)]
  structure(list(
    statistic = label, observed = observed, null = null,
    p_greater = if (length(null)) mean(null >= observed) else NaN,
    p_less = if (length(null)) mean(null <= observed) else NaN,
    conditioning = conditioning, reps = reps, dropped = dropped,
    n = n, directed = x$directed
  ), class = "cug_test")
}

# The statistic must be one number; NA and NaN are allowed and dropped.
as_statistic <- function(value, where) {
  if (!is.numeric(value) || length(value) != 1L) {
    stop(sprintf("`statistic` must return a single number; on %s it returned %s of length %d.",
                 where, describe_class(value), length(value)), call. = FALSE)
  }
  as.double(value)
}

# Distinct non-loop ties of the simple reading, and the dyad census of a
# directed network (mutual, asymmetric; null is the rest).
simple_dyads <- function(x) {
  src <- x$index$src
  tgt <- x$index$tgt
  keep <- src != tgt
  src <- src[keep]
  tgt <- tgt[keep]
  n <- nrow(x$nodes)
  if (x$directed) {
    arcs <- unique(src * n + tgt)
    reverse <- (arcs %% n) * n + arcs %/% n
    mutual_arcs <- sum(reverse %in% arcs)
    list(edges = length(arcs), mutual = mutual_arcs %/% 2L, asym = length(arcs) - mutual_arcs)
  } else {
    pairs <- unique(pmin(src, tgt) * n + pmax(src, tgt))
    list(edges = length(pairs), mutual = length(pairs), asym = 0L)
  }
}

#' @rdname cug_test
#' @export
print.cug_test <- function(x, ...) {
  cat("Conditional uniform graph test\n\n")
  cat(sprintf("Statistic: %s\n", x$statistic))
  cat(sprintf("Conditioning: %s (%s, %d nodes, %d random networks%s)\n",
              switch(x$conditioning, size = "size", edges = "size and number of edges",
                     dyad.census = "size and dyad census"),
              if (x$directed) "directed" else "undirected", x$n, x$reps,
              if (x$dropped) sprintf(", %d undefined and dropped", x$dropped) else ""))
  cat(sprintf("Observed value: %s\n", format(x$observed, digits = 6)))
  cat(sprintf("Pr(X >= Obs): %s\n", format(x$p_greater, digits = 4)))
  cat(sprintf("Pr(X <= Obs): %s\n", format(x$p_less, digits = 4)))
  invisible(x)
}

#' @rdname cug_test
#' @export
summary.cug_test <- function(object, ...) {
  q <- if (length(object$null)) {
    stats::quantile(object$null, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), names = TRUE)
  } else {
    stats::setNames(rep(NaN, 7L), c("0%", "2.5%", "25%", "50%", "75%", "97.5%", "100%"))
  }
  structure(c(unclass(object), list(quantiles = q, mean = mean(object$null),
                                    sd = stats::sd(object$null))),
            class = "summary.cug_test")
}

#' @export
print.summary.cug_test <- function(x, ...) {
  print.cug_test(structure(x, class = "cug_test"))
  cat("\nNull distribution:\n")
  cat(sprintf("  mean %s, sd %s\n", format(x$mean, digits = 6), format(x$sd, digits = 6)))
  print(x$quantiles, digits = 6)
  invisible(x)
}
