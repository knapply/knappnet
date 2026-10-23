\name{aggregate_network}
\alias{aggregate_network}
\title{Aggregate a temporal network into a weighted network}
\description{
Collapses the contacts between each pair of nodes (ordered pairs in
directed networks, unordered otherwise) into one edge carrying the number
of contacts, their total duration, and the first and last times. With
\code{from} and \code{to}, only contacts active in \eqn{[from, to)} count
and durations are clipped to the window, so an instantaneous contact adds
to the count but not to the duration.
}
\usage{
aggregate_network(x, weight = c("count", "duration"), from = NULL,
  to = NULL)
}
\arguments{
\item{x}{A temporal network (\code{\link{temporal_network}}).}
\item{weight}{Which column becomes the weight: \code{"count"} or
\code{"duration"}. Either is a tie strength
(\code{weight_type = "strength"}), so path-based measures use 1/weight as
the length.}
\item{from, to}{Optional window; both or neither.}
}
\details{
Self-loops are kept as edges from a node to itself. Edges appear in the
order of each pair's first contact; in undirected networks the endpoints
are given in the order of that first contact. The original weights, if
any, are not used. All nodes are kept, with their columns. The result is
the time-aggregated (static) network of Holme and Saramäki (2012), with
contact counts or total contact durations as weights.
}
\value{
A static \code{\link{netframe}} whose edges have columns for the
endpoints, \code{count}, \code{duration} (in the unit of the times),
\code{first} (earliest start), and \code{last} (latest end).
}
\references{
Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
}
\seealso{
\code{\link{network_window}}, \code{\link{snapshots}},
\code{\link{rolling}}.
}
\examples{
calls <- data.frame(from = c("a", "b", "a", "c"), to = c("b", "a", "c", "a"),
                    time = c(0, 10, 20, 25), end = c(5, 12, 20, 40))
g <- temporal_network(calls, end = "end", directed = FALSE)
aggregate_network(g)$edges
aggregate_network(g, weight = "duration", from = 0, to = 30)$edges
# Weighted static measures on the aggregate.
degree(aggregate_network(g), weighted = TRUE)
}
