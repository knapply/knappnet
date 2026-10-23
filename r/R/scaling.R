#' Classical multidimensional scaling
#'
#' Places nodes in `k` dimensions so that their Euclidean distances
#' approximate a matrix of dissimilarities, by classical (metric) scaling
#' (Torgerson 1952; Gower 1966): the squared dissimilarities are
#' double-centred, and the coordinates are the leading eigenvectors of the
#' result scaled by the square roots of their eigenvalues. This is the
#' method of `stats::cmdscale()`, whose results it reproduces up to the sign
#' of each dimension, and of UCINET's metric MDS, used in Everton (2012,
#' Appendix 3) to draw networks by their geodesic or structural-equivalence
#' distances (see [mds_distances()]). `mds_stress()` reports how well any
#' configuration fits a dissimilarity matrix.
#'
#' @param d A square numeric matrix of dissimilarities (symmetric,
#'   non-negative, zero on the diagonal), or a `"dist"` object.
#' @param k The number of dimensions, between 1 and one less than the
#'   number of points.
#' @param points A numeric matrix with one row per point and one column per
#'   dimension.
#'
#' @details
#' Any dissimilarity can be scaled, Euclidean or not; for non-Euclidean
#' input some eigenvalues are negative, and only dimensions with a
#' positive eigenvalue get coordinates, so fewer than `k` columns can be
#' returned (an eigenvalue counts as positive when it exceeds \eqn{10^{-12}}
#' times the largest; `cmdscale()` keeps any positive value, so it can
#' return extra dimensions of rounding error). The sign of each dimension
#' is arbitrary in the method; here it is fixed so that the first
#' coordinate of each dimension that is not negligible is positive, which
#' makes results reproducible across platforms. Dimensions with equal
#' eigenvalues are not unique: any rotation within them fits equally well.
#'
#' Goodness of fit is `cmdscale()`'s pair of ratios, the sum of the first
#' `k` eigenvalues over the sum of the absolute values of all eigenvalues
#' and over the sum of the positive ones (equal for Euclidean input).
#' Stress is Kruskal's (1964) stress formula 1,
#' \deqn{\sqrt{\sum_{i<j} (d_{ij} - \delta_{ij})^2 / \sum_{i<j} d_{ij}^2},}{sqrt(sum (d_ij - delta_ij)^2 / sum d_ij^2),}
#' with \eqn{d} the distances between the points and \eqn{\delta} the
#' dissimilarities themselves (metric stress, so no monotone regression is
#' involved; it is 0 when the configuration reproduces the dissimilarities
#' exactly). Non-metric scaling (Kruskal 1964) is not provided; see
#' `MASS::isoMDS()`.
#'
#' The eigenproblem is solved in the compiled engine by Householder
#' reduction to tridiagonal form and the implicit QL algorithm (Golub and
#' Van Loan 2013), so results do not depend on R's LAPACK build.
#'
#' @references
#' Torgerson, W. S. (1952). Multidimensional scaling: I. Theory and method.
#' *Psychometrika*, 17(4), 401-419. \doi{10.1007/BF02288916}
#'
#' Gower, J. C. (1966). Some distance properties of latent root and vector
#' methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.
#' \doi{10.1093/biomet/53.3-4.325}
#'
#' Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness
#' of fit to a nonmetric hypothesis. *Psychometrika*, 29(1), 1-27.
#' \doi{10.1007/BF02289565}
#'
#' Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.).
#' Johns Hopkins University Press.
#'
#' Everton, S. F. (2012). *Disrupting Dark Networks*. Cambridge University
#' Press. \doi{10.1017/CBO9781139136877}
#'
#' @return `classical_mds()` returns a list with `points` (a matrix with one
#'   row per point, rownames from `d`, and up to `k` columns `Dim1`,
#'   `Dim2`, ...), `eigenvalues` (all of them, in decreasing order), `gof`
#'   (the two goodness-of-fit ratios), and `stress`. `mds_stress()` returns
#'   a number.
#' @seealso [mds_distances()], [layout_network()] for drawing networks.
#' @examples
#' g <- example_network("florentine")
#' fit <- classical_mds(mds_distances(g), k = 2)
#' fit$stress
#' head(fit$points)
#' plot(fit$points, type = "n")
#' text(fit$points, labels = rownames(fit$points), cex = 0.7)
#' @export
classical_mds <- function(d, k = 2) {
  d <- as_dissimilarity(d)
  n <- nrow(d)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || k != floor(k)) {
    stop("`k` must be a positive whole number.", call. = FALSE)
  }
  if (n < 2L) {
    stop("Classical scaling needs at least two points.", call. = FALSE)
  }
  if (k > n - 1L) {
    stop(sprintf("`k` must be at most %d, one less than the number of points.", n - 1L),
         call. = FALSE)
  }
  k <- as.integer(k)
  r <- engine_call("classical_mds", .Call(C_classical_mds, d, n, k), edge_rows = FALSE)
  dims <- r[[2L]]
  points <- r[[1L]][, seq_len(dims), drop = FALSE]
  dimnames(points) <- list(rownames(d), if (dims) paste0("Dim", seq_len(dims)))
  stress <- if (dims) {
    engine_call("classical_mds", .Call(C_mds_stress, d, n, points, dims), edge_rows = FALSE)
  } else {
    NaN
  }
  list(points = points, eigenvalues = r[[3L]], gof = r[[4L]], stress = stress)
}

#' @rdname classical_mds
#' @export
mds_stress <- function(d, points) {
  d <- as_dissimilarity(d)
  n <- nrow(d)
  if (!is.numeric(points) || !is.matrix(points)) {
    stop("`points` must be a numeric matrix with one row per point.", call. = FALSE)
  }
  if (nrow(points) != n) {
    stop(sprintf("`points` has %d rows but `d` has %d points.", nrow(points), n),
         call. = FALSE)
  }
  if (anyNA(points) || any(!is.finite(points))) {
    stop("`points` must be finite.", call. = FALSE)
  }
  storage.mode(points) <- "double"
  engine_call("mds_stress", .Call(C_mds_stress, d, n, points, ncol(points)), edge_rows = FALSE)
}

# A dissimilarity matrix as a square double matrix. Structural checks
# (symmetry, zero diagonal, finiteness) are made by the engine.
as_dissimilarity <- function(d) {
  if (inherits(d, "dist")) d <- as.matrix(d)
  if (!is.numeric(d) || !is.matrix(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a square numeric matrix of dissimilarities, or a \"dist\" object.",
         call. = FALSE)
  }
  if (anyNA(d)) {
    stop("`d` must not contain missing values.", call. = FALSE)
  }
  storage.mode(d) <- "double"
  d
}

#' Distances between nodes for scaling
#'
#' Dissimilarity matrices of the kinds Everton (2012, Appendix 3) scales
#' with [classical_mds()]: geodesic distances, or structural-equivalence
#' distances between nodes' tie profiles.
#'
#' @inheritParams betweenness
#' @param method `"geodesic"` for shortest-path lengths with direction
#'   ignored (see [distances()]); `"euclidean"`, `"hamming"`, or
#'   `"correlation"` for [structural_equivalence()] distances, with
#'   correlations \eqn{r} converted to the dissimilarity \eqn{1 - r}.
#' @param unreachable For geodesic distances, what to use for pairs with no
#'   path between them: `"max_plus_one"` (the default) uses the longest
#'   finite distance plus 1, as UCINET does for closeness; `"diameter"`
#'   uses the longest finite distance; or a number.
#'
#' @details
#' Geodesic distances ignore direction so that the matrix is symmetric, as
#' classical scaling requires. Weights are handled as in [distances()]:
#' distances, or strengths converted to lengths 1/weight. Structural
#' equivalence distances use ties as given (total weights when the
#' network is weighted, or 0/1 otherwise).
#'
#' @references
#' Everton, S. F. (2012). *Disrupting Dark Networks*. Cambridge University
#' Press. \doi{10.1017/CBO9781139136877}
#'
#' @return A symmetric numeric matrix with node identifiers as dimnames.
#' @seealso [classical_mds()], [distances()], [structural_equivalence()].
#' @export
mds_distances <- function(x, method = c("geodesic", "euclidean", "hamming", "correlation"),
                          unreachable = "max_plus_one", weighted = NULL, weights_as = NULL,
                          threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  if (method == "geodesic") {
    d <- distances(x, mode = "all", weighted = weighted, weights_as = weights_as,
                   threads = threads)
    if (any(is.infinite(d))) {
      finite <- d[is.finite(d)]
      longest <- if (length(finite)) max(finite) else 0
      fill <- if (is.numeric(unreachable)) {
        if (length(unreachable) != 1L || is.na(unreachable) || unreachable < 0) {
          stop("`unreachable` must be \"max_plus_one\", \"diameter\", or a non-negative number.",
               call. = FALSE)
        }
        unreachable
      } else {
        switch(check_choice(unreachable, c("max_plus_one", "diameter"), "unreachable"),
               max_plus_one = longest + 1, diameter = longest)
      }
      d[is.infinite(d)] <- fill
    }
    return(d)
  }
  m <- structural_equivalence(x, method = method, weighted = weighted, threads = threads)
  if (method == "correlation") {
    m <- 1 - m
    diag(m) <- 0
  }
  m
}
