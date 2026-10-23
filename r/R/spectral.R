#' PageRank
#'
#' The stationary distribution of a random walk that, at each step, follows
#' one of the current node's out-ties with probability `damping` and
#' otherwise jumps to a node chosen at random (Brin and Page 1998).
#'
#' @inheritParams eigenvector_centrality
#' @param damping Probability of following a tie rather than jumping,
#'   between 0 and 1.
#' @param personalized Jump probabilities: a numeric vector with one
#'   nonnegative value per node (aligned with `x$nodes`), or the name of
#'   such a column of `x$nodes`. It is rescaled to sum to 1. `NULL` jumps
#'   to every node with equal probability.
#' @param directed Whether to follow the direction of ties in directed
#'   networks. `FALSE` treats every tie as reciprocal.
#' @param weighted Whether to choose ties in proportion to their weights,
#'   used as given whatever the network's `weight_type`. `NULL` uses weights
#'   when the network has a weight column. They must not be negative.
#' @param dangling Where a walker at a node without out-ties (or whose
#'   out-ties all have weight 0) jumps: `"personalized"` by the
#'   personalization vector, as igraph and NetworkX do, or `"uniform"` to
#'   any node with equal probability.
#' @param convention `"igraph"` or `"networkx"`; they differ only in how an
#'   undirected self-loop counts (see Details).
#' @param tol Convergence tolerance: iteration stops when the scores, which
#'   sum to 1, change by at most `tol` in total.
#' @param max_iter Iteration cap. Not converging within it is an error.
#'
#' @details
#' Scores are computed by power iteration from equal scores. The error
#' after stopping is at most `damping / (1 - damping)` times `tol`. Each
#' step runs in parallel on large networks.
#'
#' Parallel edges add up, as their number or their weights. In undirected
#' networks a self-loop counts as two tie ends under the igraph convention
#' (so a walker at the node returns to it with probability 2/degree) and as
#' one under NetworkX's; in directed networks it counts once. The treatment
#' of dangling nodes follows Langville and Meyer (2004). sna and UCINET have
#' no PageRank, so those conventions are not available.
#'
#' Results match igraph's `page_rank()` and NetworkX's `pagerank()` to
#' within their convergence tolerances.
#'
#' @references
#' Brin, S., & Page, L. (1998). The anatomy of a large-scale hypertextual
#' Web search engine. *Computer Networks and ISDN Systems*, 30(1-7),
#' 107-117. \doi{10.1016/S0169-7552(98)00110-X}
#'
#' Langville, A. N., & Meyer, C. D. (2004). Deeper inside PageRank.
#' *Internet Mathematics*, 1(3), 335-380.
#' \doi{10.1080/15427951.2004.10129091}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`, summing to
#'   1.
#' @export
pagerank <- function(x, damping = 0.85, personalized = NULL, directed = TRUE, weighted = NULL,
                     dangling = c("personalized", "uniform"),
                     convention = c("igraph", "networkx"), tol = 1e-12, max_iter = 10000L,
                     threads = 0L) {
  check_netframe(x)
  dangling <- match.arg(dangling)
  convention <- match.arg(convention)
  if (!is.numeric(damping) || length(damping) != 1L || is.na(damping) ||
      damping < 0 || damping > 1) {
    stop("`damping` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(personalized)) {
    personalized <- node_values(x, personalized, "personalized")
  }
  directed <- check_flag(directed, "directed")
  weighted <- positions_weighted(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("pagerank",
    .Call(C_pagerank, ptr, as.double(damping), directed, weighted, personalized,
          if (dangling == "personalized") 0L else 1L,
          if (convention == "igraph") 0L else 1L, as.double(tol), as.double(max_iter),
          threads))
}

#' Katz centrality
#'
#' Counts the walks that end at each node, with a walk of length k weighted
#' by \eqn{\alpha^k} (Katz 1953): the solution of
#' \eqn{x = \alpha A^T x + \beta}{x = alpha t(A) x + beta}.
#'
#' @inheritParams pagerank
#' @param alpha The attenuation factor. It must be below
#'   \eqn{1/\lambda_{max}}, the reciprocal of the largest eigenvalue of the
#'   adjacency matrix, for the walk counts to converge; larger values are an
#'   error.
#' @param beta The exogenous term: one number for every node, a vector with
#'   one value per node, or the name of a column of `x$nodes`.
#' @param mode For directed networks: `"in"` counts walks ending at the node
#'   (NetworkX, igraph), `"out"` walks starting there, and `"all"` ignores
#'   direction.
#' @param normalized Whether to scale the scores to Euclidean length 1.
#'   `NULL` uses the convention's default: `TRUE` for NetworkX, `FALSE` for
#'   igraph.
#' @param convention `"igraph"` (`alpha_centrality()`) or `"networkx"`
#'   (`katz_centrality()`); see Details.
#' @param tol Convergence tolerance of the iteration, relative to the
#'   largest score.
#'
#' @details
#' Expanding the solution gives \eqn{x = \beta (1 + \alpha A^T 1 +
#' \alpha^2 (A^T)^2 1 + \dots)}{x = beta (1 + alpha t(A) 1 + alpha^2
#' t(A)^2 1 + ...)}, a weighted count of walks; this is NetworkX's
#' `katz_centrality()` and igraph's `alpha_centrality()` (Bonacich and Lloyd
#' 2001). Katz's original index leaves out the walks of length 0, which
#' gives \eqn{x/\beta - 1} for constant \eqn{\beta}. The series converges
#' only for \eqn{\alpha < 1/\lambda_{max}}; knappnet checks this with
#' rigorous bounds on \eqn{\lambda_{max}} (Collatz 1942) and stops with an
#' error otherwise, where NetworkX's direct solver returns numbers that are
#' no longer walk counts. Scores are found by iterating
#' \eqn{x \leftarrow \alpha A^T x + \beta}{x <- alpha t(A) x + beta}, in
#' parallel on large networks.
#'
#' Parallel edges add up. Under the igraph convention self-loops are
#' ignored and scores are not normalized; under NetworkX's a self-loop adds
#' its weight once to the diagonal and scores have unit length. Weights are
#' used as given and must not be negative.
#'
#' @references
#' Katz, L. (1953). A new status index derived from sociometric analysis.
#' *Psychometrika*, 18(1), 39-43. \doi{10.1007/BF02289026}
#'
#' Bonacich, P., & Lloyd, P. (2001). Eigenvector-like measures of centrality
#' for asymmetric relations. *Social Networks*, 23(3), 191-201.
#' \doi{10.1016/S0378-8733(01)00038-7}
#'
#' Collatz, L. (1942). Einschließungssatz für die charakteristischen Zahlen
#' von Matrizen. *Mathematische Zeitschrift*, 48(1), 221-226.
#' \doi{10.1007/BF01180013}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
katz_centrality <- function(x, alpha = 0.1, beta = 1, mode = c("in", "out", "all"),
                            weighted = NULL, normalized = NULL,
                            convention = c("igraph", "networkx"), tol = 1e-12,
                            max_iter = 10000L, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  convention <- match.arg(convention)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha < 0) {
    stop("`alpha` must be a single nonnegative number.", call. = FALSE)
  }
  if (is.character(beta) || length(beta) != 1L) beta <- node_values(x, beta, "beta")
  if (!is.numeric(beta) || anyNA(beta)) {
    stop("`beta` must be numeric without missing values.", call. = FALSE)
  }
  check_flag_or_null(normalized, "normalized")
  if (is.null(normalized)) normalized <- convention == "networkx"
  weighted <- positions_weighted(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("katz_centrality",
    .Call(C_katz, ptr, as.double(alpha), as.double(beta), mode_code(mode), weighted,
          if (convention == "igraph") 2L else 1L, normalized, as.double(tol),
          as.double(max_iter), threads))
}

#' Bonacich power centrality
#'
#' Bonacich's (1987) family of measures
#' \eqn{c(\alpha, \beta) = \alpha (I - \beta A)^{-1} A 1}{c(alpha, beta) =
#' alpha (I - beta A)^-1 A 1}: a weighted count of the walks starting at
#' each node, with a walk of length k weighted by \eqn{\beta^{k-1}}.
#'
#' @inheritParams pagerank
#' @param exponent The attenuation \eqn{\beta}. Positive values make a node
#'   powerful when its contacts are well connected; negative values when
#'   its contacts are weak (Bonacich 1987).
#' @param mode For directed networks: `"out"` counts walks along out-ties
#'   (Bonacich, igraph, sna), `"in"` along in-ties, `"all"` ignores
#'   direction.
#' @param loops Whether self-loops enter the diagonal of the adjacency
#'   matrix (the number of loops, or their weights). igraph and sna ignore
#'   the diagonal by default.
#' @param rescale Whether to divide the scores by their sum instead of
#'   scaling them so that their squares sum to the number of nodes.
#' @param convention `"igraph"` counts every edge, so parallel edges add
#'   up; `"sna"` reads the network as an adjacency matrix, so parallel edges
#'   collapse into one tie (keeping the largest weight).
#' @param tol Convergence tolerance of the iterative solver, relative to
#'   the largest score.
#'
#' @details
#' \eqn{\alpha} is chosen so that the squared scores sum to the number of
#' nodes, as in igraph's `power_centrality()` and sna's `bonpow()`.
#' When \eqn{|\beta| \lambda_{max} < 1} (checked with Collatz's bounds) the
#' walk series converges and the scores are found by iterating
#' \eqn{y \leftarrow A 1 + \beta A y}{y <- A 1 + beta A y}. Otherwise, as
#' with the common default \eqn{\beta = 1}, the series diverges, but the
#' linear system usually still has a solution; it is then solved by
#' Gaussian elimination with partial pivoting, which is limited to 5,000
#' nodes. If \eqn{1/\beta} is an eigenvalue of \eqn{A} the system is
#' singular and an error is reported. Without ties every score is `NaN`
#' (\eqn{\alpha} would be 0/0).
#'
#' In directed networks the default counts walks along out-ties, so a node
#' without out-ties scores 0 (Bonacich 1987). Weights are used as given.
#'
#' @references
#' Bonacich, P. (1987). Power and centrality: A family of measures.
#' *American Journal of Sociology*, 92(5), 1170-1182. \doi{10.1086/228631}
#'
#' Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.).
#' Johns Hopkins University Press.
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
power_centrality <- function(x, exponent = 1, mode = c("out", "in", "all"), weighted = NULL,
                             loops = FALSE, rescale = FALSE,
                             convention = c("igraph", "sna"), tol = 1e-12,
                             max_iter = 10000L, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  convention <- match.arg(convention)
  if (!is.numeric(exponent) || length(exponent) != 1L || !is.finite(exponent)) {
    stop("`exponent` must be a single finite number.", call. = FALSE)
  }
  loops <- check_flag(loops, "loops")
  rescale <- check_flag(rescale, "rescale")
  if (convention == "sna" && loops) {
    stop("The sna convention reads the network without self-loops; use `loops = FALSE`.",
         call. = FALSE)
  }
  weighted <- positions_weighted(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- if (convention == "sna") engine_graph_simple(x) else engine_graph(x)
  engine_call("power_centrality",
    .Call(C_power_centrality, ptr, as.double(exponent), mode_code(mode), weighted,
          if (loops) 1L else 2L, rescale, as.double(tol), as.double(max_iter), threads),
    edge_rows = convention == "igraph")
}

#' Hub and authority scores
#'
#' Kleinberg's (1999) hubs and authorities: a good hub points to good
#' authorities, and a good authority is pointed to by good hubs. Authority
#' scores are the leading eigenvector of \eqn{A^T A}{t(A) A} and hub scores
#' that of \eqn{A A^T}{A t(A)}.
#'
#' @inheritParams pagerank
#' @param scale `"max"` (largest score 1, igraph), `"sum"` (scores sum to
#'   1, NetworkX), or `"unit"` (Euclidean length 1). `NULL` uses the
#'   convention's choice.
#' @param convention `"igraph"` or `"networkx"`: the scaling, and how an
#'   undirected self-loop counts (twice or once on the diagonal).
#' @param tol Convergence tolerance (relative residual of the eigenvector).
#' @param max_iter Cap on matrix-vector products. Not converging within it
#'   is an error.
#'
#' @details
#' The authority vector is computed with the Lanczos method on
#' \eqn{A^T A}{t(A) A}, which is never formed, and the hub vector is
#' \eqn{A a}{A a}. Parallel edges add up; weights are used as given and
#' must not be negative.
#'
#' The scores are unique only when the largest singular value of \eqn{A} is
#' simple. When it is not (for example in bipartite or disconnected
#' networks, or a directed cycle), a warning is given and the result is the
#' limit of Kleinberg's iteration \eqn{a \leftarrow A^T h}{a <- t(A) h},
#' \eqn{h \leftarrow A a}{h <- A a} from equal hub scores, which is also what
#' NetworkX's `hits()` returns; igraph's answer then depends on its
#' starting vector. Without ties every score is `NaN`.
#'
#' @references
#' Kleinberg, J. M. (1999). Authoritative sources in a hyperlinked
#' environment. *Journal of the ACM*, 46(5), 604-632.
#' \doi{10.1145/324133.324140}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
hub_score <- function(x, weighted = NULL, scale = NULL, convention = c("igraph", "networkx"),
                      tol = 1e-12, max_iter = 10000L, threads = 0L) {
  hits_common(x, "hub_score", weighted, scale, match.arg(convention), tol, max_iter,
              threads)
}

#' @rdname hub_score
#' @export
authority_score <- function(x, weighted = NULL, scale = NULL,
                            convention = c("igraph", "networkx"), tol = 1e-12,
                            max_iter = 10000L, threads = 0L) {
  hits_common(x, "authority_score", weighted, scale, match.arg(convention), tol, max_iter,
              threads)
}

hits_common <- function(x, fn, weighted, scale, convention, tol, max_iter, threads) {
  check_netframe(x)
  if (is.null(scale)) scale <- if (convention == "igraph") "max" else "sum"
  scale <- check_choice(scale, c("max", "sum", "unit"), "scale")
  weighted <- positions_weighted(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call(fn,
    .Call(C_hits, ptr, weighted, if (convention == "igraph") 0L else 1L,
          switch(scale, max = 0L, sum = 1L, unit = 2L), as.double(tol),
          as.double(max_iter), threads))
  if (!res[[4]] && nrow(x$nodes) > 1L && !all(is.nan(res[[1]]))) {
    warning("The largest singular value of the adjacency matrix is not simple, so hub ",
            "and authority scores are not unique; the result is the limit of ",
            "Kleinberg's iteration from equal scores.", call. = FALSE)
  }
  if (fn == "hub_score") res[[1]] else res[[2]]
}

check_iteration <- function(tol, max_iter) {
  if (!is.numeric(tol) || length(tol) != 1L || !isTRUE(tol > 0)) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !isTRUE(max_iter >= 1) ||
      max_iter != floor(max_iter)) {
    stop("`max_iter` must be a positive whole number.", call. = FALSE)
  }
}

# Numeric node values: a vector aligned with x$nodes, or the name of a
# column of x$nodes.
node_values <- function(x, values, what) {
  if (is.character(values) && length(values) == 1L) {
    if (!values %in% names(x$nodes)) {
      stop(sprintf("`%s` names no column of `x$nodes`.", what), call. = FALSE)
    }
    values <- x$nodes[[values]]
  }
  if (!is.numeric(values) || length(values) != nrow(x$nodes) || anyNA(values)) {
    stop(sprintf(paste0("`%s` must be a numeric vector with one value per node (%d) and ",
                        "no missing values, or the name of such a column of `x$nodes`."),
                 what, nrow(x$nodes)), call. = FALSE)
  }
  as.double(values)
}
