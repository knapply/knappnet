\name{plot.classical_mds}
\alias{plot.classical_mds}
\title{Plot a multidimensional scaling}
\usage{
\method{plot}{classical_mds}(x, which = c("configuration", "shepard", "scree"), d = NULL,
     dims = c(1, 2), labels = NULL, ...)
}
\arguments{
\item{x}{A \code{\link{classical_mds}} result.}
\item{which}{\code{"configuration"}, \code{"shepard"}, or \code{"scree"}.}
\item{d}{For the Shepard diagram, the dissimilarity matrix (or
\code{"dist"} object) that was scaled. The result does not keep it.}
\item{dims}{The two dimensions to draw in the configuration.}
\item{labels}{Whether to label the points with the rownames of
\code{x$points}; \code{NULL} labels at most 60 points.}
\item{...}{Further arguments passed to \code{\link[graphics]{plot}}
(configuration and Shepard diagram) or \code{\link[graphics]{barplot}}
(scree plot), such as \code{main}.}
}
\value{\code{x}, invisibly.}
\description{
Draws the configuration of a \code{\link{classical_mds}} with its point
labels, a Shepard diagram of fitted distances against the input
dissimilarities, or a scree plot of the eigenvalues, the three displays
used to judge a scaling solution (Kruskal 1964; Borg and Groenen 2005).
}
\details{
The configuration is drawn with equal scales on both axes, so distances in
the plot are the fitted distances. The Shepard diagram plots the Euclidean
distance between every pair of points in the configuration (over all its
dimensions) against the dissimilarity scaled, with the line of perfect
fit; the stress \code{x$stress} summarizes the scatter. The scree plot
shows every eigenvalue in decreasing order, negative ones (the sign of
non-Euclidean input) in a contrasting color, and marks the dimensions
kept.
}
\references{
Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness of
fit to a nonmetric hypothesis. \emph{Psychometrika}, 29(1), 1-27.
\doi{10.1007/BF02289565}

Borg, I., & Groenen, P. J. F. (2005). \emph{Modern Multidimensional
Scaling: Theory and Applications} (2nd ed.). Springer.
\doi{10.1007/0-387-28981-X}
}
\seealso{\code{\link{classical_mds}}, \code{\link{mds_distances}}.}
\examples{
g <- example_network("florentine")
d <- mds_distances(g)
fit <- classical_mds(d, k = 2)
plot(fit, main = "Florentine families by geodesic distance")
plot(fit, which = "shepard", d = d)
plot(fit, which = "scree")
}
