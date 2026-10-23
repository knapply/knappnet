# Conditional uniform graph tests: the test, its print, summary, and plot
# methods.

# Conditional uniform graph tests (man/net_cug_test.Rd).
net_cug_test <- function(x, statistic, conditioning = c("size", "edges", "dyad.census"),
                     reps = 1000L, ..., seed = NULL, threads = 0L) {
  check_netframe(x)
  if (!is.function(statistic)) {
    stop("`statistic` must be a function of a netframe returning one number.", call. = FALSE)
  }
  conditioning <- match.arg(conditioning)
  reps <- check_count(reps, "reps", min = 1L)
  seed <- check_seed(seed)
  threads <- resolve_threads(threads)
  label <- deparse(substitute(statistic))[[1L]]

  n <- nrow(x$nodes)
  s <- netframe_summary(x)
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
  lists <- engine_call("net_cug_test", draw(reps, seed), edge_rows = FALSE)
  null <- vapply(lists, function(e) {
    g <- net_frame(data.frame(source = e[[1L]], target = e[[2L]]), nodes,
                   directed = x$directed, quiet = TRUE)
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

# See man/net_cug_test.Rd.
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

# See man/net_cug_test.Rd.
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

print.summary.cug_test <- function(x, ...) {
  print.cug_test(structure(x, class = "cug_test"))
  cat("\nNull distribution:\n")
  cat(sprintf("  mean %s, sd %s\n", format(x$mean, digits = 6), format(x$sd, digits = 6)))
  print(x$quantiles, digits = 6)
  invisible(x)
}

# See man/net_cug_test.Rd.
plot.cug_test <- function(x, ..., breaks = "Sturges") {
  pal <- plot_qualitative(8L)
  null <- x$null
  if (!length(null)) stop("The null distribution is empty; nothing to plot.", call. = FALSE)
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- sprintf("CUG test of %s (conditioning on %s)",
                                               x$statistic, x$conditioning)
  if (is.null(dots$xlab)) dots$xlab <- x$statistic
  xlim <- range(c(null, x$observed), finite = TRUE)
  xlim <- xlim + c(-1, 1) * 0.05 * max(diff(xlim), 1e-8)
  h <- graphics::hist(null, breaks = breaks, plot = FALSE)
  do.call(graphics::plot, c(list(h, col = pal[[2L]], border = "white", xlim = xlim), dots))
  graphics::abline(v = x$observed, col = pal[[5L]], lwd = 2.5)
  graphics::legend(if (x$observed > stats::median(null)) "topleft" else "topright",
                   bty = "n", cex = 0.85, lty = c(1, NA, NA), lwd = c(2.5, NA, NA),
                   col = c(pal[[5L]], NA, NA),
                   legend = c(sprintf("observed = %s", format(x$observed, digits = 4)),
                              sprintf("Pr(X >= obs) = %s", format(x$p_greater, digits = 3)),
                              sprintf("Pr(X <= obs) = %s", format(x$p_less, digits = 3))))
  invisible(x)
}
