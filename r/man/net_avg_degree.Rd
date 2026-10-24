\name{net_avg_degree}
\alias{net_avg_degree}
\alias{net_degree_dist}
\title{Average degree and the degree distribution}
\description{
\code{net_avg_degree()} is the mean of \code{\link{net_degree}} over the nodes:
for an undirected network \eqn{2m / n}, and for a directed one \eqn{m / n}
in each direction (\code{mode} \code{"out"} or \code{"in"}) or \eqn{2m / n}
for total degree (\code{"all"}). UCINET reports it in Network > Cohesion >
Density (Average Degree). \code{net_degree_dist()} tabulates the
degrees.
}
\usage{
net_avg_degree(x, mode = c("all", "out", "in"), loops = TRUE)

net_degree_dist(x, mode = c("all", "out", "in"), loops = TRUE, cumulative = FALSE)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{mode}{For directed networks: \code{"all"} (in-degree plus
out-degree), \code{"out"}, or \code{"in"}. Ignored for undirected
networks.}
\item{loops}{Whether self-loops count toward degree, as in
\code{\link{net_degree}}.}
\item{cumulative}{Whether to add the proportion of nodes with at least each
degree.}
}
\details{
Parallel edges each count and self-loops count as in \code{\link{net_degree}}
(twice in undirected networks when \code{loops = TRUE}). Weights are
ignored; use \code{\link{net_strength}} for weighted degree. The cumulative
column is the empirical complementary distribution function
\eqn{P(K \ge k)}, the form in which heavy-tailed degree distributions are
usually drawn (Clauset, Shalizi, and Newman 2009); fitting a power law to
it is not provided here.
}
\value{
\code{net_avg_degree()}: a number, \code{NaN} for a network without nodes.
\code{net_degree_dist()}: a data frame with one row per observed degree
value, holding \code{degree}, \code{count} (nodes with that degree),
\code{proportion} (count over the number of nodes), and, with
\code{cumulative = TRUE}, \code{at_least} (the proportion of nodes with that
degree or more).
}
\references{
Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
distributions in empirical data. \emph{SIAM Review}, 51(4), 661-703.
\doi{10.1137/070710111}
}
\seealso{\code{\link{net_degree}}, \code{\link{net_density}}}
\examples{
g <- net_example("karate")
net_avg_degree(g)
dd <- net_degree_dist(g, cumulative = TRUE)
dd
plot(dd$degree, dd$at_least, log = "xy", xlab = "degree k", ylab = "P(K >= k)")
}
