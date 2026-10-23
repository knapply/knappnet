\name{mds_distances}
\alias{mds_distances}
\title{Distances between nodes for scaling}
\description{
Dissimilarity matrices of the kinds Everton (2012, Appendix 3) scales with
\code{\link{classical_mds}}: geodesic distances, or structural-equivalence
distances between nodes' tie profiles.
}
\usage{
mds_distances(x, method = c("geodesic", "euclidean", "hamming", "correlation"),
              unreachable = "max_plus_one", weighted = NULL, weights_as = NULL,
              threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{method}{\code{"geodesic"} for shortest-path lengths with direction
ignored (see \code{\link{distances}}); \code{"euclidean"},
\code{"hamming"}, or \code{"correlation"} for
\code{\link{structural_equivalence}} distances, with correlations \eqn{r}
converted to the dissimilarity \eqn{1 - r}.}
\item{unreachable}{For geodesic distances, what to use for pairs with no
path between them: \code{"max_plus_one"} (the default) uses the longest
finite distance plus 1, as UCINET does for closeness; \code{"diameter"}
uses the longest finite distance; or a number.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}. (geodesic distances only).}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
Geodesic distances ignore direction so that the matrix is symmetric, as
classical scaling requires. Weights are handled as in
\code{\link{distances}}: distances, or strengths converted to lengths
1/weight. Structural equivalence distances use ties as given (total
weights when the network is weighted, or 0/1 otherwise).
}
\value{A symmetric numeric matrix with node identifiers as dimnames.}
\references{
Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge
University Press. \doi{10.1017/CBO9781139136877}
}
\seealso{\code{\link{classical_mds}}, \code{\link{distances}},
\code{\link{structural_equivalence}}.}
\examples{
g <- example_network("florentine")
mds_distances(g)[1:4, 1:4]
mds_distances(g, "euclidean")[1:4, 1:4]
classical_mds(mds_distances(g, "correlation"))$stress
}
