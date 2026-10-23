\name{temporal_reach}
\alias{temporal_reach}
\alias{reachability_ratio}
\title{Temporal reachability}
\description{
\code{temporal_reach()} counts, for every node, the other nodes it can
reach by time-respecting walks (\code{mode = "out"}: the size of its
\emph{set of influence}) or that can reach it (\code{mode = "in"}: its
\emph{source set}), after Holme (2005). \code{reachability_ratio()} is
the share of ordered pairs of distinct nodes joined by a time-respecting
walk (Holme 2005), which is 1 in a network where everyone can eventually
reach everyone in time order.
}
\usage{
temporal_reach(x, mode = c("out", "in"), start = NULL, end = NULL,
  min_wait = 0, max_wait = Inf, strict = FALSE, traversal = 0,
  threads = 0L)

reachability_ratio(x, start = NULL, end = NULL, min_wait = 0,
  max_wait = Inf, strict = FALSE, traversal = 0, threads = 0L)
}
\arguments{
\item{x}{A temporal network (\code{\link{temporal_network}}).}
\item{mode}{\code{"out"} or \code{"in"}.}
\item{start, end}{The window: walks leave at or after \code{start} and
arrive at or before \code{end}; \code{NULL} uses the earliest start time
and the latest arrival of any contact.}
\item{min_wait, max_wait}{Waiting times allowed at each intermediate
node; see \code{\link{temporal_distances}}.}
\item{strict}{If \code{TRUE}, waiting times must be strictly greater
than \code{min_wait}, so with \code{min_wait = 0} contacts at the same
time never chain.}
\item{traversal}{Time to traverse a contact: a number (or
\code{difftime}), \code{0} by default, or \code{"duration"}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all
available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
Walks follow the rules described in \code{\link{temporal_distances}},
including the window and waiting times: by default times along a walk are
non-decreasing (Kempe, Kleinberg, and Kumar 2002; Holme and Saramäki
2012), and with \code{strict = TRUE} they strictly increase (Pan and
Saramäki 2011). Self-loops never count; in undirected networks
\code{mode} makes no difference. A node never counts itself, so
\code{temporal_reach()} ranges from 0 to \eqn{n - 1} and
\code{reachability_ratio()} is their sum over \eqn{n(n - 1)}. Each node
costs one pass over the contacts (64 nodes per pass without
\code{max_wait}), in parallel. Results are checked against a brute-force
search over all walks.
}
\value{
\code{temporal_reach()} returns a numeric vector aligned with the rows of
\code{x$nodes}; \code{reachability_ratio()} a number (\code{NaN} with
fewer than two nodes).
}
\references{
Holme, P. (2005). Network reachability of real-world contact sequences.
\emph{Physical Review E}, 71(4), 046119.
\doi{10.1103/PhysRevE.71.046119}

Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
}
\seealso{
\code{\link{temporal_distances}}, \code{\link{reachability}} for static
networks.
}
\examples{
msgs <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"), time = c(1, 2, 1))
g <- temporal_network(msgs)
temporal_reach(g, threads = 2)
temporal_reach(g, mode = "in", threads = 2)
reachability_ratio(g, threads = 2)

# Simultaneous contacts: a -> b and b -> c both at time 1.
same <- temporal_network(data.frame(from = c("a", "b"), to = c("b", "c"), time = c(1, 1)))
# Non-decreasing times (Kempe et al. 2002; Holme and Saramaki 2012): a reaches c.
temporal_reach(same, strict = FALSE, threads = 2)
# Strictly increasing times (Pan and Saramaki 2011): a reaches only b.
temporal_reach(same, strict = TRUE, threads = 2)
}
