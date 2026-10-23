#' Structural holes: constraint, effective size, efficiency, and hierarchy
#'
#' Burt's (1992) measures of how far a node's contacts are redundant with
#' one another. `constraint()` is high when a node's contacts are few or
#' strongly tied to each other; `effective_size()` is the number of
#' contacts, discounted for redundancy; `efficiency()` is effective size
#' divided by the number of contacts; `hierarchy()` measures how far the
#' constraint on a node is concentrated in a single contact.
#'
#' A node's contacts are its in- and out-neighbours other than itself.
#' Proportional tie strength is \eqn{p_{ij} = z_{ij} / \sum_k z_{ik}} with
#' \eqn{z_{ij} = a_{ij} + a_{ji}}, where \eqn{a_{ij}} is the total weight of
#' the edges from i to j (their number in a binary network), so parallel
#' edges add up and a mutual pair of arcs is a stronger tie than a single
#' arc. Self-loops are ignored. Constraint is
#' \eqn{C_i = \sum_j (p_{ij} + \sum_q p_{iq} p_{qj})^2} and effective size
#' is \eqn{\sum_j (1 - \sum_q p_{iq} m_{jq})} with
#' \eqn{m_{jq} = z_{jq} / \max_k z_{jk}}, summing over contacts j and q,
#' \eqn{q \ne j}. In binary undirected networks without parallel edges
#' effective size equals \eqn{n - 2t/n} for n contacts with t ties among
#' them (Borgatti 1997). Hierarchy is the Coleman-Theil disorder index of
#' the constraint each contact contributes, \eqn{c_{ij} = (p_{ij} + \sum_q
#' p_{iq} p_{qj})^2}: with \eqn{N} contacts and \eqn{r_j = c_{ij} / (C_i /
#' N)}, \eqn{H_i = \sum_j r_j \ln r_j / (N \ln N)}, which is 0 when every
#' contact constrains the node equally and approaches 1 as one contact
#' accounts for all of the constraint; this is the definition in UCINET's
#' structural holes routine. All four are `NaN` for nodes without
#' contacts, and hierarchy also for nodes with a single contact, where
#' \eqn{N \ln N = 0}.
#'
#' Weights are read as tie strengths (volume of contact) and used as given,
#' whatever the network's `weight_type`; they must be positive. These
#' conventions match igraph's `constraint()` and NetworkX's `constraint()`
#' and `effective_size()` (which, unlike this package, count self-loops).
#'
#' @inheritParams betweenness
#' @param weighted Whether to use edge weights as tie strengths. `NULL` (the
#'   default) uses them when the network has a weight column.
#'
#' @references
#' Burt, R. S. (1992). *Structural Holes: The Social Structure of
#' Competition*. Harvard University Press. \doi{10.4159/9780674029095}
#'
#' Burt, R. S. (2004). Structural holes and good ideas. *American Journal
#' of Sociology*, 110(2), 349-399. \doi{10.1086/421787}
#'
#' Borgatti, S. P. (1997). Structural holes: Unpacking Burt's redundancy
#' measures. *Connections*, 20(1), 35-38.
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
constraint <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 0L, weighted, threads, "constraint")
}

#' @rdname constraint
#' @export
effective_size <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 1L, weighted, threads, "effective_size")
}

#' @rdname constraint
#' @export
efficiency <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 2L, weighted, threads, "efficiency")
}

#' @rdname constraint
#' @export
hierarchy <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 3L, weighted, threads, "hierarchy")
}

structural_holes <- function(x, measure, weighted, threads, fn) {
  check_netframe(x)
  weighted <- positions_weighted(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_structural_holes, ptr, measure, weighted, threads))
}

#' Gould-Fernandez brokerage roles
#'
#' Counts, for each node, the ordered pairs of other nodes it brokers in
#' each of the five roles of Gould and Fernandez (1989). Node j brokers
#' (i, k) when \eqn{i \to j \to k}{i -> j -> k} and there is no tie
#' \eqn{i \to k}{i -> k}. With the groups of i, j, and k written in order:
#' coordinator (A, A, A), itinerant or consultant (A, B, A), representative
#' (A, A, B), gatekeeper (A, B, B), and liaison (A, B, C).
#'
#' The roles are defined on binary ties: weights are ignored, parallel
#' edges count once, and self-loops are ignored. An undirected edge counts
#' as a tie in both directions, so in an undirected network each unordered
#' pair a node brokers is counted twice. The counts match `raw.nli` from
#' `sna::brokerage()`.
#'
#' @inheritParams betweenness
#' @param groups The group of each node: a vector aligned with the rows of
#'   `x$nodes`, or the name of a column of `x$nodes`. Only equality of
#'   values matters. Missing values are not allowed.
#'
#' @references
#' Gould, R. V., & Fernandez, R. M. (1989). Structures of mediation: A
#' formal approach to brokerage in transaction networks. *Sociological
#' Methodology*, 19, 89-126. \doi{10.2307/270949}
#'
#' @return A data frame with one row per node: `node` (the identifier),
#'   then the counts `coordinator`, `itinerant`, `representative`,
#'   `gatekeeper`, `liaison`, and `total`.
#' @export
brokerage <- function(x, groups, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  n <- nrow(x$nodes)
  if (is.character(groups) && length(groups) == 1L && groups %in% names(x$nodes)) {
    groups <- x$nodes[[groups]]
  }
  if (length(groups) != n) {
    stop(sprintf("`groups` must have one value per node (%d), or name a column of `x$nodes`.",
                 n), call. = FALSE)
  }
  if (anyNA(groups)) {
    stop("`groups` must not contain missing values.", call. = FALSE)
  }
  if (!is.atomic(groups)) {
    stop("`groups` must be an atomic vector (integer, character, factor, ...).",
         call. = FALSE)
  }
  codes <- match(groups, unique(groups)) - 1L
  ptr <- engine_graph(x)
  m <- engine_call("brokerage", .Call(C_brokerage, ptr, codes, threads))
  out <- data.frame(node = x$nodes[[x$roles$id]], stringsAsFactors = FALSE)
  roles <- c("coordinator", "itinerant", "representative", "gatekeeper",
             "liaison", "total")
  for (r in seq_along(roles)) out[[roles[[r]]]] <- m[, r]
  out
}

#' Structural equivalence
#'
#' Compares every pair of nodes by their tie profiles: the ties each node
#' sends (its row of the adjacency matrix) and receives (its column).
#' Structurally equivalent nodes (Lorrain and White 1971) have identical
#' profiles.
#'
#' Nodes i and j are compared over every other node k, so the diagonal and
#' the pair's own ties \eqn{a_{ij}}, \eqn{a_{ji}} are left out (Wasserman
#' and Faust 1994, sec. 9.2):
#' \describe{
#'   \item{`"euclidean"`}{\eqn{\sqrt{\sum_{k \ne i,j} (a_{ik} - a_{jk})^2 +
#'     (a_{ki} - a_{kj})^2}} (Burt 1976).}
#'   \item{`"correlation"`}{the Pearson correlation of the two profiles over
#'     the same \eqn{2(n-2)} entries (Wasserman and Faust 1994). It is 0
#'     when undefined (a profile without variance). The diagonal is 1, or 0
#'     for a node whose profile is constant.}
#'   \item{`"hamming"`}{the sum of absolute differences, which in a binary
#'     network counts the ties that differ.}
#' }
#' These match `sna::sedist()` with its defaults (`diag = FALSE`). In a
#' binary network parallel edges count once; with weights, \eqn{a_{ij}} is
#' the total weight of the edges from i to j, used as given. An undirected
#' edge is a tie in both directions. Self-loops are ignored.
#'
#' @inheritParams betweenness
#' @param method `"euclidean"`, `"correlation"`, or `"hamming"`.
#' @param weighted Whether ties are total edge weights rather than 0 or 1.
#'   `NULL` (the default) uses weights when the network has a weight column.
#'
#' @references
#' Lorrain, F., & White, H. C. (1971). Structural equivalence of individuals
#' in social networks. *Journal of Mathematical Sociology*, 1(1), 49-80.
#' \doi{10.1080/0022250X.1971.9989788}
#'
#' Burt, R. S. (1976). Positions in networks. *Social Forces*, 55(1),
#' 93-122. \doi{10.2307/2577097}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return A symmetric numeric matrix with node identifiers as dimnames.
#' @export
structural_equivalence <- function(x, method = c("euclidean", "correlation", "hamming"),
                                   weighted = NULL, threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  code <- switch(method, euclidean = 0L, correlation = 1L, hamming = 2L)
  weighted <- positions_weighted(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  m <- engine_call("structural_equivalence",
                   .Call(C_structural_equivalence, ptr, code, weighted, threads))
  ids <- as.character(x$nodes[[x$roles$id]])
  dimnames(m) <- list(ids, ids)
  m
}

#' CONCOR blockmodel
#'
#' Partitions nodes into positions with CONCOR (Breiger, Boorman, and Arabie
#' 1975): the matrix of profile correlations (as in
#' `structural_equivalence(x, "correlation")`) is correlated with itself
#' repeatedly until its entries reach +1 or -1, which splits the nodes in
#' two. Each block is split again, starting from its submatrix of the
#' initial correlations, `depth` times, giving at most \eqn{2^{depth}}
#' blocks of nodes with non-constant profiles.
#'
#' Nodes whose profile is constant (isolates, for example) cannot be
#' correlated; they are set aside and given blocks of their own, one per
#' distinct profile. Iteration stops when every entry has absolute value at
#' least `1 - tolerance`, when an iteration changes no entry by more than
#' \eqn{10^{-12}} (a fixed point other than +1 or -1, as in some highly
#' symmetric networks), or after `max_iter` iterations. Nodes are split by
#' the sign of their converged correlation with the block's first node
#' (positive with it, zero or negative apart). Blocks are numbered in
#' depth-first order of the splits, then the blocks of constant profiles.
#'
#' CONCOR can be numerically fragile. A warning is given when a split did
#' not converge, or when some node's correlation with the block's first
#' node came within \eqn{10^{-6}} of zero during the iteration: rounding
#' error may then decide that node's block, so another implementation or
#' platform may give a different partition.
#'
#' @inheritParams structural_equivalence
#' @param depth Number of rounds of splitting (at least 1).
#' @param tolerance Convergence tolerance, between 0 and 1.
#' @param max_iter Maximum number of iterations per split.
#'
#' @references
#' Breiger, R. L., Boorman, S. A., & Arabie, P. (1975). An algorithm for
#' clustering relational data with applications to social network analysis
#' and comparison with multidimensional scaling. *Journal of Mathematical
#' Psychology*, 12(3), 328-383. \doi{10.1016/0022-2496(75)90028-0}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return An integer vector of block numbers (starting at 1) aligned with
#'   the rows of `x$nodes`.
#' @export
concor <- function(x, depth = 1L, weighted = NULL, tolerance = 1e-7, max_iter = 100L,
                   threads = 0L) {
  check_netframe(x)
  if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) || depth < 1 ||
      depth != floor(depth) || depth > .Machine$integer.max) {
    stop("`depth` must be a whole number of at least 1.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !isTRUE(tolerance > 0) ||
      !isTRUE(tolerance < 1)) {
    stop("`tolerance` must be a number between 0 and 1.", call. = FALSE)
  }
  max_iter <- check_count(max_iter, "max_iter")
  weighted <- positions_weighted(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  r <- engine_call("concor", .Call(C_concor, ptr, as.integer(depth), weighted,
                                   as.double(tolerance), max_iter, threads))
  if (r[[2]][[1]] > 0L) {
    warning(sprintf("%d CONCOR split(s) did not converge to +1/-1 within `max_iter` iterations.",
                    r[[2]][[1]]), call. = FALSE)
  }
  if (r[[2]][[2]] > 0L) {
    warning(sprintf(paste0("%d CONCOR split(s) had a correlation within 1e-6 of zero; ",
                           "rounding error may have decided some nodes' blocks."),
                    r[[2]][[2]]), call. = FALSE)
  }
  r[[1]]
}

# Whether a measure that reads weights as given should use them.
positions_weighted <- function(x, weighted) {
  check_flag_or_null(weighted, "weighted")
  weighted <- if (is.null(weighted)) !is.null(x$roles$weight) else isTRUE(weighted)
  if (weighted) require_weights(x)
  weighted
}
