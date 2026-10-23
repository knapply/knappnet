#' Draw a network with nodes sized and colored by a score
#'
#' Draws the network with [plot.netframe()], sizing and coloring each node
#' by a score such as a centrality, with a legend named after the score.
#'
#' @param x A [netframe()].
#' @param scores A numeric vector aligned with the rows of `x$nodes`, as
#'   returned by [degree()], [betweenness()], and the other centrality
#'   functions, or the name of a numeric column of `x$nodes`.
#' @param title The name of the score, used for the legend.
#' @param color Whether to color the nodes by the score as well as sizing
#'   them. `FALSE` keeps the default node color.
#' @param ... Further arguments passed to [plot.netframe()], such as
#'   `layout`, `seed`, `node_label`, or `main`.
#'
#' @details
#' Node area grows linearly with the score and the color follows the
#' sequential palette of [plot.netframe()]; undefined scores (`NA`, `NaN`)
#' are drawn at the smallest size in light gray.
#'
#' @return The layout used, invisibly, as [plot.netframe()] returns.
#' @seealso [plot_scores()] for a ranked dot plot of the same scores.
#' @examples
#' g <- example_network("florentine")
#' plot_centrality(g, betweenness(g), title = "betweenness", seed = 1)
#' plot_centrality(g, degree(g), title = "degree", color = FALSE, layout = "circle")
#' @export
plot_centrality <- function(x, scores, title = "score", color = TRUE, ...) {
  check_netframe(x)
  color <- check_flag(color, "color")
  if (is.character(scores) && length(scores) == 1L && scores %in% names(x$nodes)) {
    if (identical(title, "score")) title <- scores
    scores <- x$nodes[[scores]]
  }
  if (!is.numeric(scores) || length(scores) != nrow(x$nodes)) {
    stop("`scores` must be a numeric vector with one value per node, or the name of a numeric column of `x$nodes`.",
         call. = FALSE)
  }
  if (!is.character(title) || length(title) != 1L || is.na(title)) {
    stop("`title` must be a single string.", call. = FALSE)
  }
  w <- with_node_column(x, as.numeric(scores), title)
  invisible(plot(w$x, node_size = w$name, node_color = if (color) w$name else NULL, ...))
}

#' Draw a network colored by a partition
#'
#' Draws the network with [plot.netframe()], coloring each node by its
#' group in a partition: communities from [cluster_louvain()], positions
#' from [concor()] or [regular_equivalence()], components, or any node
#' attribute.
#'
#' @param x A [netframe()].
#' @param membership The group of each node: a vector aligned with the rows
#'   of `x$nodes` (only equality of values matters), or the name of a column
#'   of `x$nodes`.
#' @param title The name of the partition, used for the legend.
#' @param ... Further arguments passed to [plot.netframe()], such as
#'   `layout`, `seed`, `node_size`, or `main`.
#'
#' @details
#' Groups are colored with the qualitative palette of [plot.netframe()] in
#' order of their first appearance in `x$nodes`, so group numbers from
#' [cluster_louvain()] (numbered in order of first appearance) keep their
#' order. More than 15 groups are listed in the legend as a count. When
#' `membership` carries a modularity (as [cluster_louvain()]'s result does)
#' and no `main` is given, the title reports the number of groups and the
#' modularity.
#'
#' @return The layout used, invisibly, as [plot.netframe()] returns.
#' @seealso [cluster_louvain()], [plot_matrix()] for the same partition as
#'   a sorted adjacency matrix, [plot.blockmodel()].
#' @examples
#' g <- example_network("karate")
#' groups <- cluster_louvain(g, seed = 1)
#' plot_partition(g, groups, title = "community", seed = 1)
#' plot_partition(g, "faction", seed = 1, main = "Factions")
#' @export
plot_partition <- function(x, membership, title = "group", ...) {
  check_netframe(x)
  if (is.character(membership) && length(membership) == 1L && membership %in% names(x$nodes) &&
      identical(title, "group")) {
    title <- membership
  }
  if (!is.character(title) || length(title) != 1L || is.na(title)) {
    stop("`title` must be a single string.", call. = FALSE)
  }
  q <- attr(membership, "modularity")
  g <- plot_groups(x, membership, "membership")
  w <- with_node_column(x, g, title)
  dots <- list(...)
  if (is.null(dots$main) && !is.null(q) && is.numeric(q) && length(q) == 1L) {
    dots$main <- sprintf("%d groups, modularity %s", nlevels(g), format(q, digits = 3))
  }
  invisible(do.call(plot, c(list(x = w$x, node_color = w$name), dots)))
}

#' Draw the k-core shells of a network
#'
#' Places nodes on concentric rings by [coreness()], the innermost ring
#' holding the nodes of the highest core, and colors them by coreness, the
#' layered view of core structure that Seidman (1983) describes.
#'
#' @param x A [netframe()].
#' @param mode For directed networks: `"all"`, `"out"`, or `"in"`, as in
#'   [coreness()].
#' @param seed Seed for the layout whose angles order the nodes within each
#'   ring (see Details), as in [layout_network()].
#' @param threads Number of threads for that layout; `0` uses all available
#'   cores.
#' @param ... Further arguments passed to [plot.netframe()], such as
#'   `node_label` or `main`.
#'
#' @details
#' Each distinct coreness value gets a ring, the highest in the center (a
#' single node of the highest core sits at the center itself). Within a
#' ring nodes are evenly spaced, ordered by their angle in a stress layout
#' of the whole network ([layout_network()]), so that nodes close in the
#' network tend to be close on the ring. Rings are drawn as dotted circles
#' and the coreness values are listed in a legend inside the plot.
#'
#' @references
#' Seidman, S. B. (1983). Network structure and minimum degree. *Social
#' Networks*, 5(3), 269-287. \doi{10.1016/0378-8733(83)90028-X}
#'
#' @return The layout used, invisibly: a matrix with columns `x` and `y`.
#' @seealso [coreness()].
#' @examples
#' g <- example_network("karate")
#' plot_coreness(g, seed = 1, threads = 2, main = "Karate club: k-cores")
#' @export
plot_coreness <- function(x, mode = c("all", "out", "in"), seed = NULL, threads = 0L, ...) {
  check_netframe(x)
  mode <- match.arg(mode)
  n <- nrow(x$nodes)
  k <- coreness(x, mode = mode)
  dots <- list(...)
  xy <- if (n > 0L) {
    layout_network(x, seed = seed, threads = threads)
  } else {
    cbind(x = numeric(), y = numeric())
  }
  levels_ <- sort(unique(k), decreasing = TRUE)
  out <- matrix(0, n, 2L, dimnames = list(NULL, c("x", "y")))
  radii <- numeric(length(levels_))
  for (i in seq_along(levels_)) {
    members <- which(k == levels_[i])
    r <- if (i == 1L && length(members) == 1L) 0 else i
    radii[i] <- r
    if (r == 0) next
    ang <- atan2(xy[members, 2] - mean(xy[, 2]), xy[members, 1] - mean(xy[, 1]))
    members <- members[order(ang)]
    theta <- seq(0, 2 * pi, length.out = length(members) + 1L)[-1L] + pi / 2
    out[members, 1] <- r * cos(theta)
    out[members, 2] <- r * sin(theta)
  }
  # The legend is drawn here rather than by plot.netframe() so that the
  # margins are known and the rings can be added in the same coordinates.
  mapped <- map_colors(k)
  args <- c(list(x = x, layout = out, node_color = I(mapped$col), legend = FALSE), dots)
  if (is.null(args$node_size)) args$node_size <- 1
  drawn <- do.call(plot, args)
  if (n > 0L) {
    top <- if (!is.null(dots$main)) 2.2 else 0.4
    old <- graphics::par(mar = c(0.4, 0.4, top, 0.4), xpd = NA)
    on.exit(graphics::par(old), add = TRUE)
    theta <- seq(0, 2 * pi, length.out = 181L)
    for (r in radii[radii > 0]) {
      graphics::lines(r * cos(theta), r * sin(theta), col = "gray80", lty = 3)
    }
    graphics::legend("topleft", legend = mapped$labels, title = "coreness", pch = 21,
                     pt.bg = mapped$key, col = "gray40", pt.cex = 1.6, bty = "n", cex = 0.75,
                     title.adj = 0)
  }
  invisible(drawn)
}

#' Draw an ego network
#'
#' Extracts the ego network of one node with [ego_network()] and draws it
#' with ego enlarged and colored apart from its alters.
#'
#' @inheritParams ego_network
#' @param ... Further arguments passed to [plot.netframe()], such as
#'   `layout`, `seed`, or `main`.
#'
#' @return The ego network (a [netframe()]), invisibly.
#' @seealso [ego_network()], [ego_summary()].
#' @examples
#' g <- example_network("florentine")
#' plot_ego(g, "Medici", seed = 1)
#' plot_ego(g, "Medici", order = 2, seed = 1, main = "Within two steps")
#' @export
plot_ego <- function(x, node, order = 1L, mode = c("all", "out", "in"), ...) {
  check_netframe(x)
  mode <- match.arg(mode)
  e <- ego_network(x, node, order = order, mode = mode, include_ego = TRUE)
  is_ego <- node_ids(e) == as.character(node)
  role <- factor(ifelse(is_ego, "ego", "alter"), levels = c("ego", "alter"))
  w <- with_node_column(e, role, "role")
  dots <- list(...)
  args <- c(list(x = w$x, node_color = w$name), dots)
  if (is.null(args$node_size)) args$node_size <- I(ifelse(is_ego, 2.4, 1.3))
  do.call(plot, args)
  invisible(e)
}

#' Robustness of a network to node removal
#'
#' Removes nodes one at a time in a given order and tracks how the network
#' falls apart: the share of nodes in the largest component, or the
#' fragmentation of Borgatti (2006), after each removal. Curves for several
#' removal strategies (highest degree first, highest betweenness first,
#' random) can be compared with the removal of a chosen set such as the
#' [key_players()], as in the attack-tolerance analyses of Albert, Jeong,
#' and Barabási (2000) and Holme et al. (2002).
#'
#' @param x A [netframe()].
#' @param nodes Node identifiers to remove, in order, for example
#'   `key_players(x, k)$set`. `NULL` draws only the strategies.
#' @param strategy Removal orders to draw, any of `"degree"`,
#'   `"betweenness"`, and `"random"`; `NULL` for none.
#' @param measure `"largest_component"` (its share of the original nodes)
#'   or `"fragmentation"` ([fragmentation()]).
#' @param steps How many removal counts to evaluate along each strategy
#'   curve, spread evenly from none to all but one node (the curve for
#'   `nodes` is evaluated after every removal). The network is rebuilt at
#'   each point, so the cost is `steps` times that of [components()].
#' @param mode `"weak"` or `"strong"` components, as in [components()].
#' @param seed Seed for the random order (see [set.seed()]).
#' @param threads Number of threads for betweenness; `0` uses all available
#'   cores.
#' @param label Legend label for the curve of `nodes`.
#' @param ... Further arguments passed to [graphics::plot()], such as
#'   `main`.
#'
#' @details
#' Strategy orders are fixed at the start from the initial degrees or
#' betweenness scores (Holme et al.'s "initial" attacks); they are not
#' recalculated after each removal. Degree counts every edge end and
#' betweenness is the unnormalized score of [betweenness()] with its
#' default weighting. Ties are broken by node order. Removing a node
#' removes its edges; the remaining network keeps every other node, so the
#' largest-component share is measured against the original number of
#' nodes.
#'
#' @references
#' Albert, R., Jeong, H., & Barabási, A.-L. (2000). Error and attack
#' tolerance of complex networks. *Nature*, 406(6794), 378-382.
#' \doi{10.1038/35019019}
#'
#' Holme, P., Kim, B. J., Yoon, C. N., & Han, S. K. (2002). Attack
#' vulnerability of complex networks. *Physical Review E*, 65(5), 056109.
#' \doi{10.1103/PhysRevE.65.056109}
#'
#' Borgatti, S. P. (2006). Identifying sets of key players in a social
#' network. *Computational and Mathematical Organization Theory*, 12(1),
#' 21-34. \doi{10.1007/s10588-006-7084-x}
#'
#' @return A data frame with one row per evaluated point, invisibly:
#'   `curve` (the strategy or `label`), `removed` (nodes removed),
#'   `fraction` (of the original nodes), and `value` (the measure).
#' @seealso [key_players()], [components()], [fragmentation()].
#' @examples
#' g <- example_network("jemaah_islamiyah")
#' kp <- key_players(g, 3, "fragmentation", seed = 1, threads = 2)
#' plot_robustness(g, nodes = kp$set, label = "key players", seed = 1, threads = 2)
#' plot_robustness(g, strategy = "degree", measure = "fragmentation", threads = 2)
#' @export
plot_robustness <- function(x, nodes = NULL, strategy = c("degree", "betweenness", "random"),
                            measure = c("largest_component", "fragmentation"), steps = 25,
                            mode = c("weak", "strong"), seed = NULL, threads = 0L,
                            label = "selected nodes", ...) {
  check_netframe(x)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  if (!is.null(strategy)) {
    if (!is.character(strategy) || anyNA(strategy) ||
        !all(strategy %in% c("degree", "betweenness", "random"))) {
      stop("`strategy` must be some of \"degree\", \"betweenness\", and \"random\", or NULL.",
           call. = FALSE)
    }
    strategy <- unique(strategy)
  }
  if (!is.numeric(steps) || length(steps) != 1L || is.na(steps) || steps < 2) {
    stop("`steps` must be a number of at least 2.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) || label %in% strategy) {
    stop("`label` must be a single string that is not the name of a strategy.", call. = FALSE)
  }
  n <- nrow(x$nodes)
  if (n < 2L) stop("Robustness needs a network with at least two nodes.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  evaluate <- function(removed) {
    keep <- setdiff(seq_len(n), removed)
    sub <- induced_netframe(x, keep)
    if (measure == "largest_component") {
      comp <- components(sub, mode)
      if (comp$no == 0L) 0 else max(comp$csize) / n
    } else {
      if (length(keep) < 2L) NaN else fragmentation(sub, mode)
    }
  }
  curve_for <- function(removal, counts) {
    vapply(counts, function(k) evaluate(removal[seq_len(k)]), numeric(1))
  }
  counts <- unique(round(seq(0, n - 1, length.out = min(steps, n))))
  curves <- list()
  for (s in strategy) {
    removal <- switch(s,
      degree = order(degree(x), decreasing = TRUE),
      betweenness = order(betweenness(x, threads = threads), decreasing = TRUE),
      random = sample.int(n)
    )
    curves[[s]] <- data.frame(curve = s, removed = counts, fraction = counts / n,
                              value = curve_for(removal, counts), stringsAsFactors = FALSE)
  }
  if (!is.null(nodes)) {
    pos <- set_positions(x, nodes, "nodes", allow_empty = FALSE)
    k <- 0:length(pos)
    curves[[label]] <- data.frame(curve = label, removed = k, fraction = k / n,
                                  value = curve_for(pos, k), stringsAsFactors = FALSE)
  }
  if (!length(curves)) stop("Give `nodes`, a `strategy`, or both.", call. = FALSE)
  out <- do.call(rbind, unname(curves))
  rownames(out) <- NULL
  dots <- list(...)
  ylab <- if (measure == "largest_component") "largest component (share of nodes)" else "fragmentation"
  if (is.null(dots$main)) dots$main <- "Robustness to node removal"
  if (is.null(dots$xlab)) dots$xlab <- "nodes removed (share)"
  if (is.null(dots$ylab)) dots$ylab <- ylab
  if (is.null(dots$xlim)) dots$xlim <- c(0, max(out$fraction))
  if (is.null(dots$ylim)) dots$ylim <- c(0, 1)
  names_ <- names(curves)
  cols <- plot_qualitative(length(names_))
  if (!is.null(nodes)) cols[length(cols)] <- plot_accent
  old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  do.call(graphics::plot, c(list(x = 0, y = 0, type = "n"), dots))
  graphics::grid(col = "gray90", lty = 1)
  for (i in seq_along(names_)) {
    cv <- curves[[i]]
    graphics::lines(cv$fraction, cv$value, col = cols[i], lwd = 2)
    graphics::points(cv$fraction, cv$value, col = cols[i], pch = 16, cex = 0.7)
  }
  graphics::legend("topright", legend = names_, col = cols, lwd = 2, pch = 16, bty = "n",
                   cex = 0.85)
  invisible(out)
}
