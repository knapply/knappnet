\name{average_degree}
\alias{average_degree}
\alias{degree_distribution}
\title{Average degree and the degree distribution}
\usage{
average_degree(x, mode = c("all", "out", "in"), loops = TRUE)

degree_distribution(x, mode = c("all", "out", "in"), loops = TRUE, cumulative = FALSE)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks: \code{"all"} (in-degree plus
out-degree), \code{"out"}, or \code{"in"}. Ignored for undirected
networks.}
\item{loops}{Whether self-loops count toward degree, as in
\code{\link{degree}}.}
\item{cumulative}{Whether to add the proportion of nodes with at least each
degree.}
}
\value{
\code{average_degree()}: a number, \code{NaN} for a network without nodes.
\code{degree_distribution()}: a data frame with one row per observed degree
value, holding \code{degree}, \code{count} (nodes with that degree),
\code{proportion} (count over the number of nodes), and, with
\code{cumulative = TRUE}, \code{at_least} (the proportion of nodes with that
degree or more).
}
\description{
\code{average_degree()} is the mean of \code{\link{degree}} over the nodes:
for an undirected network \eqn{2m / n}, and for a directed one \eqn{m / n}
in each direction (\code{mode} \code{"out"} or \code{"in"}) or \eqn{2m / n}
for total degree (\code{"all"}). UCINET reports it in Network > Cohesion >
Density (Average Degree). \code{degree_distribution()} tabulates the
degrees.
}
\details{
Parallel edges each count and self-loops count as in \code{\link{degree}}
(twice in undirected networks when \code{loops = TRUE}). Weights are
ignored; use \code{\link{strength}} for weighted degree. The cumulative
column is the empirical complementary distribution function
\eqn{P(K \ge k)}, the form in which heavy-tailed degree distributions are
usually drawn (Clauset, Shalizi, and Newman 2009); fitting a power law to
it is not provided here.
}
\references{
Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
distributions in empirical data. \emph{SIAM Review}, 51(4), 661-703.
\doi{10.1137/070710111}
}
\seealso{\code{\link{degree}}, \code{\link{graph_density}}}
\examples{
g <- example_network("karate")
average_degree(g)
dd <- degree_distribution(g, cumulative = TRUE)
dd
plot(dd$degree, dd$at_least, log = "xy", xlab = "degree k", ylab = "P(K >= k)")
}
