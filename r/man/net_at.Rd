\name{net_at}
\alias{net_at}
\alias{net_window}
\title{A temporal network at a time or in a window}
\description{
\code{net_at()} keeps the contacts active at one time and
\code{net_window()} those active at some time in the window
\eqn{[from, to)}, following the half-open rule of
\code{\link{net_temporal}}. Both return a temporal network with all
columns and roles, contacts in their original order, and contacts
unchanged (not clipped to the window).
}
\usage{
net_at(x, time, nodes = c("all", "present", "active"))

net_window(x, from, to, nodes = c("all", "present", "active"))
}
\arguments{
\item{x}{A temporal network (\code{\link{net_temporal}}).}
\item{time}{A time, of the type of the network's times.}
\item{from, to}{Start (included) and end (excluded) of the window.}
\item{nodes}{Which nodes to keep: \code{"all"} (the default) keeps every
node, so results of measures stay aligned across times; \code{"present"}
keeps the nodes present according to the network's node presence
intervals (every node if it has none); \code{"active"} keeps only nodes
with a contact in the result. Endpoints of the contacts kept are always
kept.}
}
\details{
A contact from \eqn{s} to \eqn{e > s} is active at \eqn{t} when
\eqn{s \le t < e}{s <= t < e} and in \eqn{[from, to)} when \eqn{s < to}
and \eqn{e > from}; an instantaneous contact (\eqn{e = s}) is active only
at \eqn{t = s} and in a window when \eqn{from \le s < to}{from <= s < to}.
Node presence intervals follow the same rule. Both functions run in
compiled code; for a series of windows use \code{\link{net_snapshots}}, which
sorts the contacts once for the whole series.
}
\value{A temporal network (\code{\link{net_frame}}).}
\seealso{
\code{\link{net_snapshots}} for a series of windows,
\code{\link{net_aggregate}} to collapse a window into weighted edges.
}
\examples{
calls <- data.frame(source = c("a", "a", "b", "c"), target = c("b", "c", "c", "d"),
                    time = c(1, 2, 2, 5), end = c(3, 2, 4, 6))
g <- net_temporal(calls, end = "end", directed = FALSE)
net_at(g, 2)$edges
# The call from a to b ends at 3, so it is no longer active then.
net_at(g, 3)$edges
net_window(g, 3, 5, nodes = "active")$edges
net_window(g, 3, 5, nodes = "active")$nodes
}
