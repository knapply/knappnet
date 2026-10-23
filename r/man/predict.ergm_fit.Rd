\name{predict.ergm_fit}
\alias{predict.ergm_fit}
\title{Tie probabilities from an exponential random graph model}
\description{
For every free dyad, the probability of a tie given the rest of the
network as observed,
\deqn{P(Y_{ij} = 1 \mid y^c_{ij}) = \mathrm{logit}^{-1}(\theta' \delta_{ij}),}{P(Y_ij = 1 | rest) = logistic(theta' delta_ij),}
where \eqn{\delta_{ij}} is the dyad's change statistic. These are the
probabilities that enter the pseudolikelihood (Strauss and Ikeda 1990).
For a dyad-independent model they are the marginal tie probabilities; for a
dependent model they are conditional on the observed state of every other
dyad and are not marginal probabilities.
}
\usage{
\method{predict}{ergm_fit}(object, type = c("response", "link"), coef = NULL, ...)
}
\arguments{
\item{object}{A fitted model from \code{\link{ergm_fit}}.}
\item{type}{\code{"response"} for probabilities, \code{"link"} for their
logits.}
\item{coef}{Coefficients to use instead of the fitted ones.}
\item{...}{Not used.}
}
\value{
A data frame with one row per free dyad of the observed network:
\code{from} and \code{to} (node identifiers; undirected dyads appear once,
with \code{from} first in \code{x$nodes}), \code{tie} (whether the dyad is
a tie in the observed network), and \code{p} (or \code{link}).
}
\references{
Strauss, D., & Ikeda, M. (1990). Pseudolikelihood estimation for social
networks. \emph{Journal of the American Statistical Association}, 85(409),
204-212. \doi{10.1080/01621459.1990.10475327}
}
\seealso{\code{\link{ergm_fit}}}
\examples{
flo <- example_network("florentine")
fit <- ergm_fit(flo ~ edges + triangle, estimate = "MPLE")
p <- predict(fit)
head(p[order(-p$p), ])
}
