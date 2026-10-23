\name{network_at}
\alias{network_at}
\alias{network_window}
\title{A temporal network at a time or in a window}
\usage{
network_at(x, time, nodes = c("all", "present", "active"))

network_window(x, from, to, nodes = c("all", "present", "active"))
}
\arguments{
\item{x}{A temporal network (\code{\link{temporal_network}}).}
\item{time}{A time, of the type of the network's times.}
\item{from, to}{Start (included) and end (excluded) of the window.}
\item{nodes}{Which nodes to keep: \code{"all"} (the default) keeps every
node, so results of measures stay aligned across times; \code{"present"}
keeps the nodes present according to the network's node presence
intervals (every node if it has none); \code{"active"} keeps only nodes
with a contact in the result. Endpoints of the contacts kept are always
kept.}
}
\value{A temporal network (\code{\link{netframe}}).}
\description{
\code{network_at()} keeps the contacts active at one time and
\code{network_window()} those active at some time in the window
\eqn{[from, to)}, following the half-open rule of
\code{\link{temporal_network}}. Both return a temporal network with all
columns and roles, contacts in their original order, and contacts
unchanged (not clipped to the window).
}
\details{
A contact from \eqn{s} to \eqn{e > s} is active at \eqn{t} when
\eqn{s \le t < e}{s <= t < e} and in \eqn{[from, to)} when \eqn{s < to}
and \eqn{e > from}; an instantaneous contact (\eqn{e = s}) is active only
at \eqn{t = s} and in a window when \eqn{from \le s < to}{from <= s < to}.
Node presence intervals follow the same rule. Both functions run in
compiled code; for a series of windows use \code{\link{snapshots}}, which
sorts the contacts once for the whole series.
}
\seealso{
\code{\link{snapshots}} for a series of windows,
\code{\link{aggregate_network}} to collapse a window into weighted edges.
}
\examples{
calls <- data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d"),
                    time = c(1, 2, 2, 5), end = c(3, 2, 4, 6))
g <- temporal_network(calls, end = "end", directed = FALSE)
network_at(g, 2)$edges
# The call from a to b ends at 3, so it is no longer active then.
network_at(g, 3)$edges
network_window(g, 3, 5, nodes = "active")$edges
network_window(g, 3, 5, nodes = "active")$nodes
}
