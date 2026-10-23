\name{temporal_closeness}
\alias{temporal_closeness}
\title{Temporal closeness}
\description{
Pan and Saramäki's (2011) temporal closeness centrality
\deqn{C_i = \frac{1}{n - 1} \sum_{j \ne i} \frac{1}{\tau_{ij}}}{C_i = (1 / (n - 1)) sum_{j != i} 1 / tau_ij}
where \eqn{\tau_{ij}} is the average temporal distance from \eqn{i} to
\eqn{j} over all starting times in the window (their Eq. 1; see
\code{type = "average"} in \code{\link{temporal_distances}}). Nodes that
cannot reach \eqn{j} at all contribute 0, as in harmonic centrality. High
values mean that information leaving the node at a random time reaches
the others quickly.
}
\usage{
temporal_closeness(x, start = NULL, end = NULL, min_wait = 0,
  max_wait = Inf, strict = FALSE, traversal = 0, threads = 0L)
}
\arguments{
\item{x}{A temporal network (\code{\link{temporal_network}}).}
\item{start, end}{The window over which starting times are averaged;
\code{NULL} uses the earliest start time and the latest arrival of any
contact.}
\item{min_wait, max_wait}{Waiting times allowed at each intermediate
node; see \code{\link{temporal_distances}}.}
\item{strict, traversal}{Path rules, as in \code{\link{temporal_distances}}.
Pan and Saramäki's own rule is \code{strict = TRUE},
\code{traversal = "duration"}: each contact starts strictly after the
previous one ends. The defaults allow simultaneous contacts to chain, as
in Holme and Saramäki (2012).}
\item{threads}{Number of threads, a whole number; \code{0} uses all
available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
The average temporal distance is computed exactly from each pair's
fastest-path profile (the walks that leave latest for each arrival), with
one pass over the contacts per node, in parallel. Values are in
reciprocal units of time (per second for date-times, per day for dates).
Self-loops, parallel contacts, direction, and waiting rules are treated
as in \code{\link{temporal_distances}}; weights and node presence
intervals are not used. Results are checked against a brute-force
integration over starting times.
}
\value{
A numeric vector aligned with the rows of \code{x$nodes}; \code{NaN} for
every node when there are fewer than two nodes or the window has no
length.
}
\references{
Pan, R. K., & Saramäki, J. (2011). Path lengths, correlations, and
centrality in temporal networks. \emph{Physical Review E}, 84(1), 016105.
\doi{10.1103/PhysRevE.84.016105}

Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
}
\seealso{
\code{\link{temporal_distances}}, \code{\link{closeness}} and
\code{\link{harmonic_centrality}} for static networks.
}
\examples{
msgs <- data.frame(from = c("a", "b", "b", "a"), to = c("b", "c", "c", "c"),
                   time = c(1, 1, 3, 6))
g <- temporal_network(msgs)
# Simultaneous contacts chain (Holme and Saramaki 2012).
temporal_closeness(g, start = 0, end = 10, threads = 2)
# Strictly increasing times (Pan and Saramaki 2011).
temporal_closeness(g, start = 0, end = 10, strict = TRUE, threads = 2)
# The average temporal distances behind the first result.
temporal_distances(g, "average", start = 0, end = 10, threads = 2)
}
