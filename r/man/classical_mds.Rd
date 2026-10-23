\name{classical_mds}
\alias{classical_mds}
\alias{mds_stress}
\title{Classical multidimensional scaling}
\usage{
classical_mds(d, k = 2)

mds_stress(d, points)
}
\arguments{
\item{d}{A square numeric matrix of dissimilarities (symmetric,
non-negative, zero on the diagonal), or a \code{"dist"} object.}
\item{k}{The number of dimensions, between 1 and one less than the number
of points.}
\item{points}{A numeric matrix with one row per point and one column per
dimension.}
}
\value{\code{classical_mds()} returns a list with \code{points} (a matrix
with one row per point, rownames from \code{d}, and up to \code{k} columns
\code{Dim1}, \code{Dim2}, ...), \code{eigenvalues} (all of them, in
decreasing order), \code{gof} (the two goodness-of-fit ratios), and
\code{stress}. \code{mds_stress()} returns a number.}
\description{
Places points in \code{k} dimensions so that their Euclidean distances
approximate a matrix of dissimilarities, by classical (metric) scaling
(Torgerson 1952; Gower 1966): the squared dissimilarities are
double-centred, and the coordinates are the leading eigenvectors of the
result scaled by the square roots of their eigenvalues. This is the method
of \code{stats::cmdscale()}, whose results it reproduces up to the sign of
each dimension, and of UCINET's metric MDS, used in Everton (2012,
Appendix 3) to draw networks by their geodesic or structural-equivalence
distances (see \code{\link{mds_distances}}). \code{mds_stress()} reports
how well any configuration fits a dissimilarity matrix.
}
\details{
Any dissimilarity can be scaled, Euclidean or not; for non-Euclidean input
some eigenvalues are negative, and only dimensions with a positive
eigenvalue get coordinates, so fewer than \code{k} columns can be returned
(an eigenvalue counts as positive when it exceeds \eqn{10^{-12}} times the
largest; \code{cmdscale()} keeps any positive value, so it can return extra
dimensions of rounding error). The sign of each dimension is arbitrary in
the method; here it is fixed so that the first coordinate of each
dimension that is not negligible is positive, which makes results
reproducible across platforms. Dimensions with equal eigenvalues are not
unique: any rotation within them fits equally well.

Goodness of fit is \code{cmdscale()}'s pair of ratios, the sum of the first
\code{k} eigenvalues over the sum of the absolute values of all
eigenvalues and over the sum of the positive ones (equal for Euclidean
input). Stress is Kruskal's (1964) stress formula 1,
\deqn{\sqrt{\sum_{i<j} (d_{ij} - \delta_{ij})^2 / \sum_{i<j} d_{ij}^2},}{sqrt(sum (d_ij - delta_ij)^2 / sum d_ij^2),}
with \eqn{d} the distances between the points and \eqn{\delta} the
dissimilarities themselves (metric stress, so no monotone regression is
involved; it is 0 when the configuration reproduces the dissimilarities
exactly). Non-metric scaling (Kruskal 1964) is not provided; see
\code{MASS::isoMDS()}.

The eigenproblem is solved in the compiled engine by Householder reduction
to tridiagonal form and the implicit QL algorithm (Golub and Van Loan
2013), so results do not depend on R's LAPACK build. The validation suite
compares coordinates, eigenvalues, and fit with
\code{cmdscale(d, k, eig = TRUE)}.
}
\references{
Torgerson, W. S. (1952). Multidimensional scaling: I. Theory and method.
\emph{Psychometrika}, 17(4), 401-419. \doi{10.1007/BF02288916}

Gower, J. C. (1966). Some distance properties of latent root and vector
methods used in multivariate analysis. \emph{Biometrika}, 53(3-4),
325-338. \doi{10.1093/biomet/53.3-4.325}

Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness of
fit to a nonmetric hypothesis. \emph{Psychometrika}, 29(1), 1-27.
\doi{10.1007/BF02289565}

Golub, G. H., & Van Loan, C. F. (2013). \emph{Matrix Computations} (4th
ed.). Johns Hopkins University Press.

Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge
University Press. \doi{10.1017/CBO9781139136877}
}
\seealso{\code{\link{mds_distances}}, \code{\link{layout_network}} for
drawing networks.}
\examples{
g <- example_network("florentine")
fit <- classical_mds(mds_distances(g), k = 2)
fit$stress
head(fit$points)
plot(fit$points, type = "n")
text(fit$points, labels = rownames(fit$points), cex = 0.7)

# The same configuration as stats::cmdscale(), up to the sign of each axis.
cm <- stats::cmdscale(mds_distances(g), k = 2)
all.equal(abs(unname(cm)), abs(unname(fit$points)))
}
