\name{net_temporal_distances}
\alias{net_temporal_distances}
\title{Time-respecting paths}
\description{
Optimal time-respecting walks between nodes of a temporal network:
earliest arrival, latest departure, fastest, and shortest (fewest
contacts), or Pan and Saramäki's (2011) average temporal distance.
}
\usage{
net_temporal_distances(x, type = c("earliest_arrival", "latest_departure",
  "fastest", "shortest", "average"), from = NULL, to = NULL,
  start = NULL, end = NULL, min_wait = 0, max_wait = Inf, strict = FALSE,
  traversal = 0, threads = 0L)
}
\arguments{
\item{x}{A temporal network (\code{\link{net_temporal}}).}
\item{type}{\code{"earliest_arrival"}, \code{"latest_departure"},
\code{"fastest"}, \code{"shortest"}, or \code{"average"}; see Details.}
\item{from, to}{Identifiers of the origins (rows) and destinations
(columns); \code{NULL} uses every node.}
\item{start, end}{The window: walks leave at or after \code{start} and
arrive at or before \code{end}, given as times of the network's type.
\code{NULL} uses the earliest start time and the latest arrival of any
contact.}
\item{min_wait, max_wait}{Waiting times allowed at each intermediate
node, from arriving to leaving again (a number in the unit of the times,
or a \code{difftime}). \code{max_wait = Inf} sets no limit.}
\item{strict}{If \code{TRUE}, waiting times must be strictly greater than
\code{min_wait}; with \code{min_wait = 0}, consecutive contacts must have
strictly increasing times, so contacts at the same time never chain.}
\item{traversal}{Time to traverse a contact: a number (or
\code{difftime}), \code{0} by default, or \code{"duration"} for each
contact's own duration (its end minus its start), as for a flight or a
phone call that must finish before the information moves on.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
}
\details{
\strong{Walks.} Following Wu et al. (2014), each contact from \eqn{u} to
\eqn{v} starting at \eqn{t} is a \emph{temporal arc} that leaves \eqn{u}
at \eqn{t} and arrives at \eqn{v} at \eqn{t + \lambda}, where
\eqn{\lambda} is the traversal time. A time-respecting walk is a sequence
of arcs, each leaving the node where the previous one arrived, such that
the waiting time \eqn{w} between arriving and leaving satisfies
\eqn{w \ge}{w >=} \code{min_wait} (or \eqn{w >} \code{min_wait} with
\code{strict = TRUE}) and \eqn{w \le}{w <=} \code{max_wait}. The default,
non-decreasing times (\code{min_wait = 0}, \code{strict = FALSE}), is the
definition of Kempe, Kleinberg, and Kumar (2002) and of Holme and
Saramäki (2012), under which contacts at the same time can chain
(\eqn{a \to b}{a -> b} and \eqn{b \to c}{b -> c} at the same time let
\eqn{a} reach \eqn{c}). With \code{strict = TRUE} times strictly
increase, as in Pan and Saramäki (2011), who also require each contact to
start after the previous one ends (\code{strict = TRUE},
\code{traversal = "duration"}). The walk may leave its origin at any time
in the window; waiting rules apply only at intermediate nodes. Undirected
contacts can be traversed either way; self-loops never help; parallel
contacts are separate contacts. Interval contacts are traversed from
their start time; to let walks join a long contact part-way, split it
into shorter contacts. Node presence intervals and weights are not used.

\strong{Values} for walks from \eqn{i} to \eqn{j} (row \eqn{i}, column
\eqn{j}):
\describe{
\item{\code{"earliest_arrival"}}{the earliest arrival time at \eqn{j}
(Wu et al. 2014); \code{start} (the window start) on the diagonal.}
\item{\code{"latest_departure"}}{the latest time to leave \eqn{i} and
still reach \eqn{j} by \code{end}; \code{end} on the diagonal,
\code{-Inf} without a walk.}
\item{\code{"fastest"}}{the smallest elapsed time from leaving \eqn{i} to
arriving at \eqn{j}.}
\item{\code{"shortest"}}{the fewest contacts.}
\item{\code{"average"}}{Pan and Saramäki's (2011, Eq. 1) average
temporal distance: the time from a starting time \eqn{t} to the earliest
arrival at \eqn{j}, averaged over all \eqn{t} in the window, where a start
after the last walk waits for the first walk of the window repeated after
its end (their periodic boundary condition). It needs a window of
positive length (otherwise \code{NaN}).}
}
Unreachable pairs get \code{Inf} (\code{-Inf} for latest departure).
Times are numbers on the scale of the data: as given for numeric times,
seconds since 1970-01-01 UTC for date-times, and days since then for
dates. Lengths of time are in the unit of the times. Without
\code{max_wait}, an optimal walk never needs to repeat a node, so the
values are those of time-respecting paths; with \code{max_wait}, waiting
can only be cut short by moving on, and the optimum is taken over walks.

\strong{Algorithms.} Contacts are sorted once by time into an edge stream
and each origin costs one pass over it (Wu et al. 2014): latest departure
is earliest arrival in reversed time, fastest walks keep the latest
departure that reaches each arrival, and fewest contacts the smallest
count. Waiting limits use a heap of pending arrivals and a monotone queue
of usable ones. Earliest arrival without \code{max_wait} runs 64 origins
per pass. Origins run in parallel. Latest departure costs one pass per
destination, and \code{"average"} one pass per node of the network.
Results are checked against a brute-force search over all walks.
}
\value{
A numeric matrix with one row per origin and one column per destination,
named by node identifiers.
}
\references{
Kempe, D., Kleinberg, J., & Kumar, A. (2002). Connectivity and inference
problems for temporal networks. \emph{Journal of Computer and System
Sciences}, 64(4), 820-842. \doi{10.1006/jcss.2002.1829}

Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}

Pan, R. K., & Saramäki, J. (2011). Path lengths, correlations, and
centrality in temporal networks. \emph{Physical Review E}, 84(1), 016105.
\doi{10.1103/PhysRevE.84.016105}

Wu, H., Cheng, J., Huang, S., Ke, Y., Lu, Y., & Xu, Y. (2014). Path
problems in temporal graphs. \emph{Proceedings of the VLDB Endowment},
7(9), 721-732. \doi{10.14778/2732939.2732945}
}
\seealso{
\code{\link{net_temporal_reach}} for counts of reachable nodes,
\code{\link{net_temporal_closeness}} for a centrality built on average
temporal distances, \code{\link{net_distances}} for static paths.
}
\examples{
# b passes a message from a to c only if it hears it first.
msgs <- data.frame(source = c("a", "b", "b", "a"), target = c("b", "c", "c", "c"),
                   time = c(1, 1, 3, 6))
g <- net_temporal(msgs)
net_temporal_distances(g, threads = 2)
# Non-decreasing times (Kempe et al. 2002; Holme and Saramaki 2012): the
# contacts a -> b and b -> c at time 1 chain, so a reaches c at time 1.
net_temporal_distances(g, from = "a", strict = FALSE, threads = 2)
# Strictly increasing times (Pan and Saramaki 2011): a reaches c at 3.
net_temporal_distances(g, from = "a", strict = TRUE, threads = 2)
net_temporal_distances(g, "latest_departure", threads = 2)
net_temporal_distances(g, "fastest", threads = 2)
net_temporal_distances(g, "shortest", threads = 2)
# At most one time unit of waiting at b rules out a -> b (1), b -> c (3),
# so a reaches c only by the direct contact at 6.
net_temporal_distances(g, from = "a", strict = TRUE, max_wait = 1, threads = 2)

# Contacts that last: b hears from a during [0, 2) and calls c at 1 and 3.
calls <- data.frame(source = c("a", "b", "b"), target = c("b", "c", "c"),
                    time = c(0, 1, 3), end = c(2, 2, 4))
h <- net_temporal(calls, end = "end")
# Traversal takes no time: the message is at b from time 0.
net_temporal_distances(h, from = "a", threads = 2)
# Each contact must finish first (Pan and Saramaki 2011): c hears at 4.
net_temporal_distances(h, from = "a", strict = TRUE, traversal = "duration",
                   threads = 2)
}
