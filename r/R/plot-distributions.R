#' Plot the degree distribution
#'
#' Draws how many nodes have each degree, as a bar chart of counts or as
#' the complementary cumulative distribution \eqn{P(D \ge d)}{P(D >= d)},
#' which on logarithmic axes is the usual display for judging heavy tails
#' (Newman 2003; Clauset, Shalizi, and Newman 2009).
#'
#' @param x A [netframe()].
#' @param mode,loops,convention As in [degree()].
#' @param type `"histogram"` for counts per degree, or `"ccdf"` for the
#'   share of nodes with at least each degree.
#' @param log Whether to use logarithmic axes. Nodes of degree 0 cannot be
#'   shown on a logarithmic degree axis and are left out, with a note in
#'   the subtitle.
#' @param ... Further arguments passed to [graphics::barplot()] (histogram
#'   on linear axes) or [graphics::plot()] (otherwise), such as `main`.
#'
#' @details
#' Degrees are counted by [degree()], so parallel edges each count and a
#' self-loop adds 2 in an undirected network, unless `loops = FALSE`. The
#' histogram has one bar per degree value from 0 to the largest (not binned
#' intervals); with more than 60 distinct values, or with `log = TRUE`, it
#' is drawn as points instead. The CCDF is a step function over the
#' observed degrees; a straight line on logarithmic axes suggests a power
#' law, but a fitted model and a goodness-of-fit test (Clauset et al. 2009)
#' are needed to claim one.
#'
#' @references
#' Newman, M. E. J. (2003). The structure and function of complex networks.
#' *SIAM Review*, 45(2), 167-256. \doi{10.1137/S003614450342480}
#'
#' Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
#' distributions in empirical data. *SIAM Review*, 51(4), 661-703.
#' \doi{10.1137/070710111}
#'
#' @return A data frame with one row per observed degree, invisibly:
#'   `degree`, `count`, and `ccdf` (the share of nodes with at least that
#'   degree).
#' @seealso [degree()], [plot_scores()].
#' @examples
#' g <- example_network("karate")
#' plot_degree_distribution(g)
#' plot_degree_distribution(g, type = "ccdf", log = TRUE)
#' @export
plot_degree_distribution <- function(x, mode = c("all", "out", "in"), loops = TRUE,
                                     type = c("histogram", "ccdf"), log = FALSE,
                                     convention = c("igraph", "networkx", "sna", "ucinet"),
                                     ...) {
  check_netframe(x)
  mode <- match.arg(mode)
  type <- match.arg(type)
  convention <- match.arg(convention)
  log <- check_flag(log, "log")
  d <- degree(x, mode = mode, loops = loops, convention = convention)
  n <- length(d)
  if (n == 0L) stop("The network has no nodes.", call. = FALSE)
  values <- sort(unique(d))
  count <- as.integer(table(factor(d, levels = values)))
  ccdf <- rev(cumsum(rev(count))) / n
  out <- data.frame(degree = values, count = count, ccdf = ccdf)
  dots <- list(...)
  label <- if (x$directed && mode != "all") paste0(mode, "-degree") else "degree"
  zero <- sum(d == 0)
  if (is.null(dots$main)) dots$main <- sprintf("Degree distribution (%d nodes)", n)
  if (is.null(dots$xlab)) dots$xlab <- label
  old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  if (type == "histogram" && !log && max(values) <= 60) {
    all_values <- seq.int(0L, max(values))
    heights <- integer(length(all_values))
    heights[match(values, all_values)] <- count
    if (is.null(dots$ylab)) dots$ylab <- "nodes"
    if (is.null(dots$col)) dots$col <- plot_blue
    if (is.null(dots$border)) dots$border <- NA
    if (is.null(dots$names.arg)) dots$names.arg <- as.character(all_values)
    if (is.null(dots$space)) dots$space <- 0.2
    do.call(graphics::barplot, c(list(height = heights), dots))
  } else {
    keep <- rep(TRUE, length(values))
    if (log) keep <- values > 0
    y <- if (type == "histogram") count else ccdf
    if (is.null(dots$ylab)) {
      dots$ylab <- if (type == "histogram") "nodes" else "share of nodes with at least this degree"
    }
    if (is.null(dots$log)) dots$log <- if (log) "xy" else ""
    if (is.null(dots$pch)) dots$pch <- 16
    if (is.null(dots$col)) dots$col <- plot_blue
    if (is.null(dots$type)) dots$type <- if (type == "ccdf") "s" else "p"
    if (type == "ccdf" && is.null(dots$ylim) && !log) dots$ylim <- c(0, 1)
    if (!any(keep)) stop("Every node has degree 0; nothing to draw on logarithmic axes.",
                         call. = FALSE)
    do.call(graphics::plot, c(list(x = values[keep], y = y[keep]), dots))
    if (type == "ccdf") graphics::points(values[keep], y[keep], pch = 16, col = plot_blue)
    if (log && zero > 0) {
      graphics::mtext(sprintf("%d node%s of degree 0 not shown", zero, if (zero == 1L) "" else "s"),
                      side = 3, line = 0.3, cex = 0.8, col = "gray30")
    }
  }
  invisible(out)
}

#' Plot the distribution of component sizes
#'
#' Draws how many components of each size the network has, with the
#' largest component, which most measures of connectivity are about,
#' marked.
#'
#' @param x A [netframe()].
#' @param mode `"weak"` or `"strong"`, as in [components()].
#' @param log Whether to use logarithmic axes (sizes and counts).
#' @param ... Further arguments passed to [graphics::barplot()] (or
#'   [graphics::plot()] with `log = TRUE`), such as `main`.
#'
#' @return A data frame with one row per distinct size, invisibly:
#'   `size` and `count` (the number of components of that size).
#' @seealso [components()], [fragmentation()].
#' @examples
#' g <- example_network("noordin", relation = "friendship")
#' components(g)$no
#' plot_component_sizes(g)
#' @export
plot_component_sizes <- function(x, mode = c("weak", "strong"), log = FALSE, ...) {
  check_netframe(x)
  mode <- match.arg(mode)
  log <- check_flag(log, "log")
  comp <- components(x, mode)
  if (comp$no == 0L) stop("The network has no nodes.", call. = FALSE)
  sizes <- sort(unique(comp$csize))
  count <- as.integer(table(factor(comp$csize, levels = sizes)))
  out <- data.frame(size = sizes, count = count)
  dots <- list(...)
  if (is.null(dots$main)) {
    dots$main <- sprintf("%d %s component%s", comp$no, mode, if (comp$no == 1L) "" else "s")
  }
  if (is.null(dots$xlab)) dots$xlab <- "component size (nodes)"
  if (is.null(dots$ylab)) dots$ylab <- "components"
  old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  if (log) {
    if (is.null(dots$log)) dots$log <- "xy"
    if (is.null(dots$pch)) dots$pch <- 16
    if (is.null(dots$col)) dots$col <- plot_blue
    do.call(graphics::plot, c(list(x = sizes, y = count), dots))
  } else {
    cols <- rep(plot_blue, length(sizes))
    cols[length(sizes)] <- plot_accent
    if (is.null(dots$col)) dots$col <- cols
    if (is.null(dots$border)) dots$border <- NA
    if (is.null(dots$names.arg)) dots$names.arg <- as.character(sizes)
    do.call(graphics::barplot, c(list(height = count), dots))
  }
  graphics::mtext(sprintf("largest: %d of %d nodes", max(comp$csize), length(comp$membership)),
                  side = 3, line = 0.3, cex = 0.8, col = "gray30")
  invisible(out)
}

#' Plot the triad census
#'
#' Draws the [triad_census()] as a bar chart and, by default, the counts
#' expected under the U|MAN null model of Holland and Leinhardt (1976),
#' which conditions on the dyad census: all networks with the same numbers
#' of mutual, asymmetric, and null dyads are equally likely. Comparing the
#' two shows which triadic structures (transitive triads, cycles, closed
#' triangles) are over- or under-represented beyond what the density and
#' reciprocity of the network already imply.
#'
#' @param x A [netframe()].
#' @param expected Whether to draw the expected counts under U|MAN beside
#'   the observed ones.
#' @param log Whether to use a logarithmic count axis (counts of 0 are
#'   drawn at the bottom).
#' @param threads Number of threads for the census; `0` uses all available
#'   cores.
#' @param ... Further arguments passed to [graphics::barplot()], such as
#'   `main`.
#'
#' @details
#' With \eqn{M}, \eqn{A}, and \eqn{N} mutual, asymmetric, and null dyads
#' among the \eqn{D = n(n-1)/2} dyads, the three dyads of a random triad
#' under U|MAN have \eqn{m} mutual, \eqn{a} asymmetric, and \eqn{k} null
#' dyads with probability
#' \deqn{\binom{M}{m} \binom{A}{a} \binom{N}{k} / \binom{D}{3},}{choose(M, m) choose(A, a) choose(N, k) / choose(D, 3),}
#' and each orientation of the asymmetric dyads is equally likely, so the
#' expected count of a triad type is \eqn{\binom{n}{3}}{choose(n, 3)} times
#' that probability times the share of orientations giving the type
#' (1/4, 1/4, 1/2 for 021D, 021U, 021C; 1/2, 1/2 for 111D, 111U; 3/4, 1/4
#' for 030T, 030C; 1/4, 1/4, 1/2 for 120D, 120U, 120C; 1 otherwise). These
#' are the expected values of Holland and Leinhardt (1976; Wasserman and
#' Faust 1994, ch. 14) and were checked against the average census over
#' every digraph on five nodes. The dyad census reads the network as
#' [triad_census()] does: self-loops ignored, parallel edges counted once,
#' and every edge of an undirected network mutual, so only the four
#' undirected types (003, 102, 201, 300) are drawn for undirected networks.
#'
#' @references
#' Holland, P. W., & Leinhardt, S. (1976). Local structure in social
#' networks. *Sociological Methodology*, 7, 1-45. \doi{10.2307/270703}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return A data frame with one row per triad type drawn, invisibly:
#'   `type`, `observed`, and `expected` (`NA` when `expected = FALSE`).
#' @seealso [triad_census()], [transitivity()], [reciprocity()].
#' @examples
#' g <- netframe(data.frame(from = c(1, 2, 3, 3, 4, 5, 5, 6, 1),
#'                          to   = c(2, 3, 1, 4, 5, 3, 6, 4, 6)))
#' triad_census(g)
#' plot_triad_census(g, threads = 2)
#' plot_triad_census(example_network("karate"), log = TRUE, threads = 2)
#' @export
plot_triad_census <- function(x, expected = TRUE, log = FALSE, threads = 0L, ...) {
  check_netframe(x)
  expected <- check_flag(expected, "expected")
  log <- check_flag(log, "log")
  observed <- triad_census(x, threads = threads)
  exp <- if (expected) expected_triad_census(x) else rep(NA_real_, 16L)
  names(exp) <- names(observed)
  types <- names(observed)
  if (!x$directed) types <- c("003", "102", "201", "300")
  out <- data.frame(type = types, observed = unname(observed[types]),
                    expected = unname(exp[types]), stringsAsFactors = FALSE)
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- "Triad census"
  if (is.null(dots$ylab)) dots$ylab <- "triads"
  if (is.null(dots$xlab)) dots$xlab <- "triad type"
  if (is.null(dots$las)) dots$las <- 2
  if (is.null(dots$border)) dots$border <- NA
  heights <- if (expected) rbind(observed = out$observed, expected = out$expected) else out$observed
  if (is.null(dots$col)) dots$col <- if (expected) c(plot_blue, plot_gray) else plot_blue
  if (is.null(dots$names.arg)) dots$names.arg <- types
  if (expected && is.null(dots$beside)) dots$beside <- TRUE
  if (log) {
    if (is.null(dots$log)) dots$log <- "y"
    floor_ <- 0.5
    heights[!is.na(heights) & heights < floor_] <- floor_
    if (is.null(dots$ylim)) dots$ylim <- c(floor_, max(2, heights, na.rm = TRUE) * 1.5)
  }
  old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  do.call(graphics::barplot, c(list(height = heights), dots))
  if (expected) {
    graphics::legend("topright", bty = "n", cex = 0.8, fill = c(plot_blue, plot_gray),
                     border = NA, legend = c("observed", "expected under U|MAN"))
  }
  invisible(out)
}

# Expected triad census under the U|MAN distribution (Holland and Leinhardt
# 1976): uniform over digraphs with the dyad census of `x`. Returned in the
# order of triad_census().
expected_triad_census <- function(x) {
  n <- nrow(x$nodes)
  if (n < 3L) return(rep(0, 16L))
  src <- as.double(x$index$src)
  tgt <- as.double(x$index$tgt)
  keep <- src != tgt
  src <- src[keep]
  tgt <- tgt[keep]
  lo <- pmin(src, tgt)
  hi <- pmax(src, tgt)
  key <- lo * n + hi
  if (x$directed) {
    forward <- unique(key[src < tgt])
    backward <- unique(key[src > tgt])
    m <- length(intersect(forward, backward))
    a <- length(forward) + length(backward) - 2 * m
  } else {
    m <- length(unique(key))
    a <- 0
  }
  d <- choose(n, 2)
  k <- d - m - a
  triads <- choose(n, 3)
  p <- function(mm, aa, kk) {
    choose(m, mm) * choose(a, aa) * choose(k, kk) / choose(d, 3)
  }
  c(
    triads * p(0, 0, 3),                              # 003
    triads * p(0, 1, 2),                              # 012
    triads * p(1, 0, 2),                              # 102
    triads * p(0, 2, 1) / 4, triads * p(0, 2, 1) / 4, # 021D, 021U
    triads * p(0, 2, 1) / 2,                          # 021C
    triads * p(1, 1, 1) / 2, triads * p(1, 1, 1) / 2, # 111D, 111U
    triads * p(0, 3, 0) * 3 / 4,                      # 030T
    triads * p(0, 3, 0) / 4,                          # 030C
    triads * p(2, 0, 1),                              # 201
    triads * p(1, 2, 0) / 4, triads * p(1, 2, 0) / 4, # 120D, 120U
    triads * p(1, 2, 0) / 2,                          # 120C
    triads * p(2, 1, 0),                              # 210
    triads * p(3, 0, 0)                               # 300
  )
}

#' Ranked dot plot of node scores
#'
#' Draws scores such as centralities as a dot plot sorted from the highest
#' down, labelled with node identifiers: the quickest way to read off who
#' leads on a measure and by how much.
#'
#' @param scores A numeric vector of scores, one per node, as returned by
#'   [degree()], [betweenness()], and the other centrality functions.
#' @param labels Labels for the scores: a [netframe()] whose nodes the
#'   scores are aligned with (its identifiers are used), a character
#'   vector of the same length, or `NULL` for `names(scores)` or, failing
#'   that, positions.
#' @param top How many of the highest scores to show; `Inf` shows all.
#' @param decreasing Whether the highest score is at the top. `FALSE` ranks
#'   from the lowest instead.
#' @param ... Further arguments passed to [graphics::dotchart()], such as
#'   `main` or `xlab`.
#'
#' @details
#' Missing and `NaN` scores (undefined values) are left out and counted in
#' the subtitle. Ties keep the order of `scores`.
#'
#' @return A data frame of the scores drawn, in the order drawn from the
#'   top, invisibly: `label` and `score`.
#' @seealso [plot_centrality()] to draw the scores on the network.
#' @examples
#' g <- example_network("florentine")
#' plot_scores(betweenness(g), g, main = "Betweenness")
#' plot_scores(closeness(g), g, top = 5, xlab = "closeness")
#' @export
plot_scores <- function(scores, labels = NULL, top = 25, decreasing = TRUE, ...) {
  if (!is.numeric(scores)) stop("`scores` must be a numeric vector.", call. = FALSE)
  decreasing <- check_flag(decreasing, "decreasing")
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top < 1) {
    stop("`top` must be a positive number.", call. = FALSE)
  }
  k <- length(scores)
  lab <- if (inherits(labels, "netframe")) {
    check_netframe(labels)
    node_ids(labels)
  } else if (is.null(labels)) {
    if (is.null(names(scores))) as.character(seq_len(k)) else names(scores)
  } else {
    as.character(labels)
  }
  if (length(lab) != k) {
    stop("`labels` must give one label per score.", call. = FALSE)
  }
  ok <- !is.na(scores)
  dropped <- sum(!ok)
  s <- scores[ok]
  lab <- lab[ok]
  if (!length(s)) stop("There are no defined scores to plot.", call. = FALSE)
  o <- order(s, decreasing = decreasing)
  o <- o[seq_len(min(length(o), top))]
  s <- s[o]
  lab <- lab[o]
  dots <- list(...)
  if (is.null(dots$pch)) dots$pch <- 19
  if (is.null(dots$color)) dots$color <- plot_blue
  if (is.null(dots$lcolor)) dots$lcolor <- "gray85"
  if (is.null(dots$xlab)) dots$xlab <- "score"
  if (is.null(dots$cex)) dots$cex <- if (length(s) > 15) 0.75 else 0.9
  left <- label_lines(lab, dots$cex) + 1.5
  old <- graphics::par(mar = c(4.6, max(4.1, left), if (is.null(dots$main)) 2.1 else 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  # dotchart() draws the first element at the bottom, so reverse the order
  # to put the highest score at the top.
  do.call(graphics::dotchart, c(list(x = rev(s), labels = rev(lab)), dots))
  if (dropped > 0L || length(o) < sum(ok)) {
    notes <- character()
    if (length(o) < sum(ok)) notes <- c(notes, sprintf("top %d of %d", length(o), sum(ok)))
    if (dropped > 0L) notes <- c(notes, sprintf("%d undefined left out", dropped))
    graphics::mtext(paste(notes, collapse = "; "), side = 3, line = 0.3, cex = 0.8,
                    col = "gray30")
  }
  invisible(data.frame(label = lab, score = s, stringsAsFactors = FALSE))
}
