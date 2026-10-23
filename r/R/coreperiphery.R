#' Core-periphery structure
#'
#' Fits the core-periphery models of Borgatti and Everett (1999), as in
#' UCINET's Network > Core/Periphery > Categorical and Continuous routines.
#' In the discrete model every node is in the core or the periphery; the
#' ideal pattern has every core-core pair tied and every periphery-periphery
#' pair untied, and the core is chosen to maximize the Pearson correlation
#' between the adjacency matrix and the pattern. In the continuous model
#' every node has a coreness \eqn{c_i \ge 0} and the pattern is
#' \eqn{c_i c_j}, so that ties are expected between core nodes, less so
#' between core and peripheral nodes, and least between peripheral nodes;
#' coreness is fitted by least squares (MINRES, UCINET's method) or by
#' maximizing the correlation.
#'
#' @details
#' Pairs of distinct nodes enter the fit: each unordered pair once in
#' undirected networks and both ordered pairs in directed ones. A pair's
#' value is 1 for a tie (parallel edges count once) or, with
#' `weighted = TRUE`, the total weight of the edges between the two nodes in
#' that direction; weights must be non-negative. Self-loops are ignored. A
#' network whose pairs all have the same value (no ties, or every tie) has
#' no core-periphery structure and is an error.
#'
#' **Discrete model.** `links` says what the pattern expects of
#' core-periphery pairs: `"missing"` leaves them out of the correlation
#' (Borgatti and Everett's equation 4, the UCINET default), `"ones"` expects
#' them tied (equation 2), and `"zeros"` expects them untied (equation 3).
#' The search is tabu search (Glover 1989) over single-node moves from
#' `starts` random cores, as in [factions()]; UCINET's routine also uses a
#' combinatorial search whose details are not published. The best core over
#' all starts is returned, and the validation suite checks it against
#' exhaustive search on small networks. The correlation is undefined for
#' cores whose pattern has no variance (for example a core of one node with
#' `links = "missing"` in a two-node network); such cores are never
#' returned, and if no core has a defined correlation the result is `NaN`.
#'
#' **Continuous model, `fit = "minres"`.** Coreness minimizes the residual
#' sum of squares \eqn{\sum_{i \ne j} (a_{ij} - c_i c_j)^2} over
#' \eqn{c \ge 0}: the minimum residual method of factor analysis (Comrey
#' 1962) with one factor and the diagonal left out, which Borgatti and
#' Everett suggest and Boyd, Fitzgerald, Mahutga, and Smith (2010) describe
#' for this model; it is what UCINET's continuous routine computes. Each
#' coreness in turn is set to the value that minimizes the residuals given
#' the others (clamped at 0), and sweeps continue until no coreness changes
#' by more than `tol` relative to the largest, or `max_iter` sweeps have
#' been made. The start is coreness proportional to degree. On data such as
#' a star, whose leaves share no ties, the residual sum has no minimizer (it
#' keeps falling as the hub's coreness grows); the sweeps then stop at
#' `max_iter`, and the unit-norm coreness, which still converges, is
#' reported.
#'
#' **Continuous model, `fit = "correlation"`.** Coreness maximizes the
#' Pearson correlation between the data and the pattern, the criterion
#' Borgatti and Everett state, by the quasi-Newton method of Broyden,
#' Fletcher, Goldfarb, and Shanno with an Armijo line search (Nocedal and
#' Wright 2006, chapters 3 and 6) on the exact gradient, over \eqn{t} with
#' \eqn{c = t^2} so that coreness stays non-negative (without the
#' constraint the product pattern also fits two separate dense groups with
#' coreness of opposite signs). The correlation has several local maxima in
#' general, so the search runs from `starts` points (degree plus one, then
#' random points from the seeded generator) and keeps the best. On many
#' networks the supremum is only approached as all coreness values become
#' equal, where the pattern reduces to the additive \eqn{c_i + c_j} up to
#' scale; the returned coreness is then nearly uniform and the correlation
#' is that of the best additive fit. MINRES does not have this degeneracy,
#' which is why it is the default.
#'
#' Both fits are scaled to unit Euclidean norm (the correlation does not
#' change under scaling, and UCINET reports MINRES coreness so). The
#' validation suite checks MINRES against scipy's bounded quasi-Newton
#' minimizer and the correlation fit against scipy's BFGS from many starts.
#' Starts run in parallel, and the result for a given `seed` does not depend
#' on `threads`.
#'
#' @inheritParams betweenness
#' @param method `"discrete"` (categorical core and periphery) or
#'   `"continuous"` (a coreness score for every node).
#' @param links For the discrete model, the pattern value of core-periphery
#'   pairs: `"missing"`, `"ones"`, or `"zeros"`; see Details.
#' @param fit For the continuous model: `"minres"` (least squares, UCINET)
#'   or `"correlation"`; see Details.
#' @param weighted Whether to use edge weights (total weight between two
#'   nodes) as the pair values instead of binary ties; `NULL` uses them when
#'   the network has a weight column.
#' @param starts Number of starting points (random cores for the discrete
#'   model, points for the continuous correlation fit).
#' @param max_stall Steps without improvement after which a discrete start
#'   ends.
#' @param tenure Tabu tenure for the discrete search; `NULL` chooses
#'   `max(7, n / 4)`.
#' @param max_iter Iteration limit per start (sweeps, for MINRES) for the
#'   continuous model.
#' @param tol Convergence tolerance of MINRES: the largest change in
#'   coreness over a sweep, relative to the largest coreness, below which
#'   sweeps stop.
#' @param seed Seed for the starts, a whole number between 0 and
#'   \eqn{2^{53}}; `NULL` draws one from R's random number generator.
#'
#' @references
#' Borgatti, S. P., & Everett, M. G. (1999). Models of core/periphery
#' structures. *Social Networks*, 21(4), 375-395.
#' \doi{10.1016/S0378-8733(99)00019-2}
#'
#' Boyd, J. P., Fitzgerald, W. J., Mahutga, M. C., & Smith, D. A. (2010).
#' Computing continuous core/periphery structures for social relations data
#' with MINRES/SVD. *Social Networks*, 32(2), 125-137.
#' \doi{10.1016/j.socnet.2009.09.003}
#'
#' Comrey, A. L. (1962). The minimum residual method of factor analysis.
#' *Psychological Reports*, 11(1), 15-18. \doi{10.2466/pr0.1962.11.1.15}
#'
#' Glover, F. (1989). Tabu search, part I. *ORSA Journal on Computing*,
#' 1(3), 190-206. \doi{10.1287/ijoc.1.3.190}
#'
#' Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.).
#' Springer. \doi{10.1007/978-0-387-40065-5}
#'
#' @return A list of class `core_periphery` with `method`; `core`, a logical
#'   vector aligned with the rows of `x$nodes` (`TRUE` for core nodes;
#'   discrete model only); `fit`, `coreness`, a numeric vector aligned with
#'   `x$nodes`, `residual` (the MINRES residual sum of squares, `NaN` for
#'   the correlation fit), and `iterations` (continuous model only);
#'   `correlation`, the correlation between the data and the fitted
#'   pattern; and `network`, the netframe fitted. It has `print` and `plot`
#'   methods.
#' @export
core_periphery <- function(x, method = c("discrete", "continuous"),
                           links = c("missing", "ones", "zeros"),
                           fit = c("minres", "correlation"), weighted = NULL,
                           starts = 20L, max_stall = 200L, tenure = NULL, max_iter = 1000L,
                           tol = 1e-10, seed = NULL, threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  links <- match.arg(links)
  fit <- match.arg(fit)
  weighted <- resolve_additive_weights(x, weighted)
  starts <- check_count(starts, "starts")
  max_stall <- check_count(max_stall, "max_stall")
  tenure <- if (is.null(tenure)) 0L else check_count(tenure, "tenure")
  max_iter <- check_count(max_iter, "max_iter")
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0 || !is.finite(tol)) {
    stop("`tol` must be a single non-negative number.", call. = FALSE)
  }
  seed <- groups_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  if (method == "discrete") {
    code <- switch(links, missing = 0L, ones = 1L, zeros = 2L)
    res <- engine_call("core_periphery",
                       .Call(C_core_periphery, ptr, code, weighted, starts, max_stall, tenure,
                             seed, threads))
    out <- list(method = "discrete", links = links, core = res[[1L]], correlation = res[[2L]],
                network = x)
  } else {
    res <- engine_call("core_periphery",
                       .Call(C_core_periphery_continuous, ptr,
                             if (fit == "minres") 0L else 1L, weighted, max_iter,
                             as.double(tol), starts, seed, threads))
    out <- list(method = "continuous", fit = fit, coreness = res[[1L]],
                correlation = res[[2L]], residual = res[[3L]], iterations = res[[4L]],
                network = x)
  }
  structure(out, class = "core_periphery")
}

#' @export
print.core_periphery <- function(x, ...) {
  if (x$method == "discrete") {
    cat(sprintf("Discrete core-periphery fit (core-periphery pairs %s): %d core, %d periphery\n",
                switch(x$links, missing = "left out", ones = "expected tied",
                       zeros = "expected untied"),
                sum(x$core), sum(!x$core)))
  } else {
    cat(sprintf("Continuous core-periphery fit (%s): coreness of %d nodes (unit norm)\n",
                if (x$fit == "minres") "MINRES" else "maximum correlation",
                length(x$coreness)))
    if (x$fit == "minres") {
      cat(sprintf("Residual sum of squares: %s\n", format(x$residual, digits = 4)))
    }
  }
  cat(sprintf("Correlation with the ideal pattern: %s\n", format(x$correlation, digits = 4)))
  invisible(x)
}

#' Small-world quotient
#'
#' Watts and Strogatz (1998) call a network a small world when its
#' clustering is much higher than that of a random network of the same
#' size and density while its mean distance is about as short. The
#' small-world quotient (Uzzi and Spiro 2005; the "small-world-ness" of
#' Humphries and Gurney 2008) combines the two ratios:
#' \deqn{Q = \frac{C / C_{random}}{L / L_{random}}}{Q = (C / C_random) / (L / L_random)}
#' Here \eqn{C_{random}} and \eqn{L_{random}} are the means over `reps`
#' simulated \eqn{G(n, m)} networks with the same number of nodes and edges
#' (see [random_gnm()]), rather than the approximations
#' \eqn{C_{random} \approx k/n} and \eqn{L_{random} \approx \ln n / \ln k}
#' used in some textbooks, which fail for small or sparse networks.
#'
#' @details
#' Everything is computed on the undirected simple graph underlying the
#' network (directions ignored, parallel edges once, self-loops ignored),
#' without weights, and the random networks are undirected and simple with
#' as many edges as that graph. \eqn{C} is the global transitivity (three
#' times the triangles over the connected triples; `local = FALSE`, igraph's
#' `transitivity(type = "global")`) or the mean local clustering
#' coefficient of Watts and Strogatz with nodes of fewer than two neighbors
#' counted as 0 (`local = TRUE`, NetworkX's `average_clustering()` and
#' igraph's `transitivity(type = "average", isolates = "zero")`). \eqn{L}
#' is the mean distance over pairs of nodes joined by a path (igraph's
#' `mean_distance(unconnected = TRUE)`). Parts that are undefined (no
#' connected triples, no connected pairs) give `NaN`.
#'
#' Replicate \eqn{r} uses random stream `(seed, r)` and replicates run in
#' parallel, so the result for a given `seed` does not depend on `threads`.
#'
#' @inheritParams betweenness
#' @param local Whether to use the mean local clustering coefficient
#'   (`TRUE`) or global transitivity (`FALSE`).
#' @param reps Number of random networks.
#' @param seed Seed for the random networks, a whole number between 0 and
#'   \eqn{2^{53}}; `NULL` draws one from R's random number generator.
#'
#' @references
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#' 'small-world' networks. *Nature*, 393(6684), 440-442.
#' \doi{10.1038/30918}
#'
#' Uzzi, B., & Spiro, J. (2005). Collaboration and creativity: The small
#' world problem. *American Journal of Sociology*, 111(2), 447-504.
#' \doi{10.1086/432782}
#'
#' Humphries, M. D., & Gurney, K. (2008). Network 'small-world-ness': A
#' quantitative method for determining canonical network equivalence. *PLoS
#' ONE*, 3(4), e0002051. \doi{10.1371/journal.pone.0002051}
#'
#' @return A list of class `small_world` with `clustering` and
#'   `mean_distance` of the network, `random_clustering` and
#'   `random_mean_distance` (means over the random networks), `quotient`,
#'   `reps`, and `local`, with `print` and `plot` methods.
#' @export
small_world <- function(x, local = FALSE, reps = 100L, seed = NULL, threads = 0L) {
  check_netframe(x)
  local <- check_flag(local, "local")
  reps <- check_count(reps, "reps")
  if (reps < 1L) stop("`reps` must be at least 1.", call. = FALSE)
  seed <- groups_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  v <- engine_call("small_world", .Call(C_small_world, ptr, local, reps, seed, threads))
  structure(list(clustering = v[[1L]], mean_distance = v[[2L]], random_clustering = v[[3L]],
                 random_mean_distance = v[[4L]], quotient = v[[5L]], reps = reps,
                 local = local, n = nrow(x$nodes)),
            class = "small_world")
}

#' @export
print.small_world <- function(x, ...) {
  cat(sprintf("Small-world quotient: %s\n", format(x$quotient, digits = 4)))
  cat(sprintf("  %s: %s observed, %s in %d random networks (ratio %s)\n",
              if (x$local) "mean local clustering" else "transitivity",
              format(x$clustering, digits = 4), format(x$random_clustering, digits = 4),
              x$reps, format(x$clustering / x$random_clustering, digits = 4)))
  cat(sprintf("  mean distance: %s observed, %s in random networks (ratio %s)\n",
              format(x$mean_distance, digits = 4), format(x$random_mean_distance, digits = 4),
              format(x$mean_distance / x$random_mean_distance, digits = 4)))
  invisible(x)
}

#' Compactness and distance-weighted fragmentation
#'
#' Borgatti (2006) measures how far a network is from falling apart with the
#' distance-weighted fragmentation
#' \deqn{{}^{D}F = 1 - \frac{\sum_{i \ne j} 1 / d_{ij}}{n (n - 1)}}{DF = 1 - sum_{i != j} (1 / d_ij) / (n (n - 1))}
#' where \eqn{1 / d_{ij} = 0} for nodes that cannot reach each other. Its
#' complement, the mean reciprocal distance over ordered pairs of distinct
#' nodes, is the network's compactness (the global efficiency of Latora and
#' Marchiori 2001). Fragmentation is 0 for a complete network and 1 for a
#' network without ties; unlike [fragmentation()] it distinguishes networks
#' that are connected but stretched out. UCINET reports both in
#' Network > Cohesion > Distance.
#'
#' Directed networks use directed distances over ordered pairs (`mode`
#' `"out"`: from \eqn{i} to \eqn{j}; `"in"`: to each node; `"all"`:
#' ignoring direction). Distances are counted in steps, or, with weights, as
#' sums of lengths as in [distances()] (weights read as distances or as
#' strengths converted to lengths \eqn{1/w}; they must be strictly
#' positive). Self-loops never shorten a path and parallel edges count
#' through the shortest. Networks with fewer than two nodes give `NaN`.
#' Compactness matches NetworkX's `global_efficiency()` on undirected
#' unweighted networks.
#'
#' @inheritParams closeness
#'
#' @references
#' Borgatti, S. P. (2006). Identifying sets of key players in a social
#' network. *Computational and Mathematical Organization Theory*, 12(1),
#' 21-34. \doi{10.1007/s10588-006-7084-x}
#'
#' Latora, V., & Marchiori, M. (2001). Efficient behavior of small-world
#' networks. *Physical Review Letters*, 87(19), 198701.
#' \doi{10.1103/PhysRevLett.87.198701}
#'
#' @return A number.
#' @export
compactness <- function(x, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
                        threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("compactness", .Call(C_compactness, ptr, mode_code(mode), w$weighted, w$code,
                                   threads))
}

#' @rdname compactness
#' @export
distance_weighted_fragmentation <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                                            weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  1 - engine_call("distance_weighted_fragmentation",
                  .Call(C_compactness, ptr, mode_code(mode), w$weighted, w$code, threads))
}

#' Average degree and the degree distribution
#'
#' `average_degree()` is the mean of [degree()] over the nodes: for an
#' undirected network \eqn{2m / n}, and for a directed one \eqn{m / n} in
#' each direction (`mode` `"out"` or `"in"`) or \eqn{2m / n} for total
#' degree (`"all"`). UCINET reports it in Network > Cohesion > Density
#' (Average Degree). `degree_distribution()` tabulates the degrees.
#'
#' Parallel edges each count and self-loops count as in [degree()]. Weights
#' are ignored; use [strength()] for weighted degree.
#'
#' @inheritParams degree
#' @param cumulative Whether to add the proportion of nodes with at least
#'   each degree.
#'
#' @return `average_degree()`: a number, `NaN` for a network without nodes.
#'   `degree_distribution()`: a data frame with one row per observed degree
#'   value, holding `degree`, `count` (nodes with that degree), `proportion`
#'   (count over the number of nodes), and, with `cumulative = TRUE`,
#'   `at_least` (proportion of nodes with that degree or more).
#' @export
average_degree <- function(x, mode = c("all", "out", "in"), loops = TRUE) {
  check_netframe(x)
  mode <- match.arg(mode)
  loops <- check_flag(loops, "loops")
  d <- degree(x, mode = mode, loops = loops)
  if (!length(d)) return(NaN)
  mean(d)
}

#' @rdname average_degree
#' @export
degree_distribution <- function(x, mode = c("all", "out", "in"), loops = TRUE,
                                cumulative = FALSE) {
  check_netframe(x)
  mode <- match.arg(mode)
  loops <- check_flag(loops, "loops")
  cumulative <- check_flag(cumulative, "cumulative")
  d <- degree(x, mode = mode, loops = loops)
  values <- sort(unique(d))
  count <- vapply(values, function(v) sum(d == v), integer(1))
  out <- data.frame(degree = values, count = count, proportion = count / length(d))
  if (cumulative) out$at_least <- rev(cumsum(rev(out$proportion)))
  rownames(out) <- NULL
  out
}
