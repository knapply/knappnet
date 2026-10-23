\name{rolling}
\alias{rolling}
\title{Apply a measure over a series of windows}
\description{
Computes a node-level (or network-level) measure on every snapshot of a
temporal network (\code{\link{snapshots}}) and returns the results as one
tidy data frame, ready for plotting or modelling change over time.
}
\usage{
rolling(x, fun, width = NULL, step = width, ..., breaks = NULL,
  from = NULL, to = NULL, nodes = c("all", "present", "active"))
}
\arguments{
\item{x}{A temporal network (\code{\link{temporal_network}}).}
\item{fun}{A function, or the name of one, taking a netframe as its first
argument.}
\item{width}{Window length: a number in the unit of the times (seconds for
date-times, days for dates) or a \code{difftime}.}
\item{step}{Distance between window starts, as \code{width}.}
\item{...}{Further arguments to \code{fun}.}
\item{breaks}{Alternatively to \code{width}, increasing times: windows are
\eqn{[breaks_k, breaks_{k+1})}.}
\item{from, to}{With \code{width}: the first window starts at \code{from}
(default: the earliest start time), and windows start up to and including
\code{to} (default: the latest end time).}
\item{nodes}{Which nodes each snapshot keeps: \code{"all"} (the default),
\code{"present"}, or \code{"active"}; see \code{\link{network_at}}.}
}
\details{
Snapshots keep every node by default, so node-level results line up
across windows. \code{fun} receives each snapshot (a temporal network
holding the contacts active in the window) and must return either one
value per node of the snapshot or a single value; any knappnet measure
works, for example \code{\link{degree}}, \code{\link{betweenness}}, or
\code{\link{graph_density}}. Static measures see each contact as an edge;
aggregate inside \code{fun} (as in the examples) to count each pair of
nodes once. Windows are formed as in \code{\link{snapshots}}; the rows
are ordered by window and then by node.
}
\value{
A data frame with columns \code{start} and \code{end} (of each window),
\code{node} (the node identifier, or \code{NA} for single values), and
\code{value}.
}
\seealso{
\code{\link{snapshots}}, \code{\link{aggregate_network}}.
}
\examples{
logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3"),
                     host = c("h1", "h1", "h2", "h2", "h3", "h1"),
                     time = c(0, 1, 5, 6, 12, 13))
g <- temporal_network(logins)
rolling(g, degree, width = 5)
rolling(g, function(s) graph_density(aggregate_network(s)), width = 10, step = 5)
# Arguments pass through to the measure.
rolling(g, degree, width = 10, step = 5, mode = "in")
}
