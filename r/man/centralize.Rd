\name{centralize}
\alias{centralize}
\title{Centralization of any scores}
\description{
Freeman's (1978) centralization index for a vector of node scores:
\eqn{\sum_i (c_{\max} - c_i)}{sum_i (c_max - c_i)}, the sum of differences
between the largest score and every score, divided by
\code{theoretical_max}, the largest value that sum can take for a network
with the same number of nodes under the measure's definition.
\code{\link{centralization}} computes the index with the theoretical maxima
of degree, betweenness, closeness, and eigenvector centrality; this
function serves any other measure whose maximum is known, or reports the
unnormalized sum when it is not.
}
\usage{
centralize(scores, theoretical_max = NULL)
}
\arguments{
\item{scores}{A numeric vector of node scores.}
\item{theoretical_max}{The largest possible value of the sum of differences
for \code{length(scores)} nodes, or \code{NULL} to return the sum itself.
It depends on the measure, its normalization, and the network's direction;
see the table in \code{\link{centralization}}.}
}
\value{A number: \code{NaN} if \code{scores} is empty or contains \code{NaN}.}
\references{
Freeman, L. C. (1978). Centrality in social networks: Conceptual
clarification. \emph{Social Networks}, 1(3), 215-239.
\doi{10.1016/0378-8733(78)90021-7}
}
\seealso{\code{\link{centralization}}.}
\examples{
g <- netframe(data.frame(from = c(1, 1, 1, 1), to = c(2, 3, 4, 5)), directed = FALSE)
d <- degree(g)
centralize(d, theoretical_max = (5 - 1) * (5 - 2))  # a star: 1
centralization(g, "degree")                        # the same
centralize(d)                                      # the unnormalized sum
}
