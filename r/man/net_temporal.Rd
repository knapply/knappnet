\name{net_temporal}
\alias{net_temporal}
\alias{net_is_temporal}
\alias{temporal_net_frame}
\alias{print.temporal_net_frame}
\alias{print.spatiotemporal_net_frame}
\title{Temporal networks}
\description{
\code{net_temporal()} creates a network whose edges are timed
contacts: each edge has a start time and, optionally, an end time or a
duration, or belongs to a panel wave, and nodes may have presence spells.
The result is a \code{\link[=net_frame]{net_frame}} of class
\code{temporal_net_frame} with extra column roles, so every static measure
still works on it (treating each contact as an edge), and the temporal
functions listed under \sQuote{See also} use the times.
\code{net_is_temporal()} tells whether a net_frame is one.
}
\usage{
net_temporal(edges, nodes = NULL, start = "time", end = NULL, duration = NULL,
             wave = NULL, node_start = NULL, node_end = NULL, directed = TRUE,
             source = "source", target = "target", id = "id", weight = NULL,
             weight_type = NULL, quiet = FALSE, graph = list())

net_is_temporal(x)

\method{print}{temporal_net_frame}(x, ..., n = 6L, width = getOption("width"))

\method{print}{spatiotemporal_net_frame}(x, ..., n = 6L, width = getOption("width"))
}
\arguments{
\item{edges}{A data frame of contacts, one row per contact, or a
\code{\link{net_frame}} whose edge table holds the times; its nodes,
endpoints, weights, and other roles are kept.}
\item{nodes}{A data frame of nodes, or \code{NULL}; see
\code{\link{net_frame}}. Must be \code{NULL} when \code{edges} is a
net_frame.}
\item{start}{The name of the column of \code{edges} holding each contact's
start time. Times may be numbers (in any unit), \code{Date}, or
\code{POSIXct} date-times. The default is the column named \code{"time"};
a net_frame that already has times keeps them.}
\item{end}{Optional name of the column of end times, of the same type as
\code{start}.}
\item{duration}{Optional name of a column of durations, an alternative to
\code{end}: numbers in the unit of the times (seconds for date-times, days
for dates) or a \code{difftime}. Without \code{end} or \code{duration},
every contact is instantaneous.}
\item{wave}{For panel data, the name of an integer column of \code{edges}
giving the wave each edge was observed in, instead of \code{start}: wave
\eqn{k} is the spell \eqn{[k, k + 1)}, so every temporal function applies
unchanged. See \code{\link{net_panel}} to build one from a list of
networks.}
\item{node_start, node_end}{Optional names of columns of \code{nodes}
holding when each node is present, of the same type as \code{start} (wave
numbers with \code{wave}). A missing value means present from the beginning
(start) or until the end (end).}
\item{directed}{Whether contacts are directed. Ignored when \code{edges}
is a net_frame, which keeps its own direction.}
\item{source, target, id, weight, weight_type, quiet, graph}{As in
\code{\link{net_frame}}.}
\item{x}{Any object (a temporal net_frame for \code{print()}).}
\item{..., n, width}{As in \code{\link[=net_frame]{print.net_frame}}.}
}
\details{
\strong{Times.} Times are stored as given in \code{x$edges} and
\code{x$nodes}; the column roles are recorded in \code{x$roles}
(\code{start}, \code{end} or \code{duration}, \code{node_start},
\code{node_end}, and the time type). Internally dates and date-times are
converted to seconds since 1970-01-01 UTC, so a date is its midnight in
UTC. Arguments that are times (such as \code{from} and \code{to} of
\code{\link{net_window}}) take the type of the data: numbers for
numeric times, and dates or date-times for dates and date-times.
Arguments that are lengths of time (window widths, waiting times) take a
\code{difftime}, or a number in the unit of the data: seconds for
date-times, days for dates. Results that are lengths of time use the same
unit.

\strong{When a contact is active.} Intervals are half-open, as the spells
of networkDynamic (Butts et al. 2026): a contact from \eqn{s} to
\eqn{e > s} is active at times \eqn{s \le t < e}{s <= t < e}, and an
instantaneous contact (\eqn{e = s}) only at \eqn{t = s}. A contact is
active in a window \eqn{[a, b)} when it is active at some time in it:
\eqn{s < b} and \eqn{e > a} for an interval, \eqn{a \le s < b}{a <= s < b}
for an instantaneous contact. So back-to-back contacts \eqn{[1, 2)} and
\eqn{[2, 3)} never overlap, consecutive windows of a tumbling series never
share an instantaneous contact, and an interval belongs to every window it
overlaps. Node presence intervals follow the same rule.

\strong{Contacts and edges.} Self-loops and repeated contacts between the
same nodes are allowed and kept; each contact is one edge of the net_frame.
Static measures applied directly to a temporal network see every contact
as an edge (a multigraph); use \code{\link{net_aggregate}} for one
weighted edge per pair of nodes, or \code{\link{net_window}} and
\code{\link{net_snapshots}} for the contacts of a period.

\strong{Classes.} The result has class \code{c("temporal_net_frame",
"net_frame")}, or \code{c("spatiotemporal_net_frame", "spatial_net_frame",
"temporal_net_frame", "net_frame")} when the nodes also carry coordinates or
geometry. The temporal functions check the class and refuse other
networks. Transformations that keep the time columns keep the class
(\code{\link{net_subnetwork}}, \code{\link{net_window}}, ...); those
that drop them, such as \code{\link{net_aggregate}}, return a static
network. The printed card says how time is recorded (contacts, edge
spells, or waves), the time range, and whether nodes have presence spells,
and adds a line on the time class.

\strong{Log data.} Tables of events such as authentication logs or network
flows (one row per event, with source, destination, and time) are contact
sequences: pass them, or a net_frame built from them, straight to
\code{net_temporal()}. Millions of contacts are fine; slicing and
paths run in compiled code after a single sort by time.
}
\value{
\code{net_temporal()} returns a net_frame of class \code{temporal_net_frame}
with time roles; \code{net_is_temporal()} returns \code{TRUE} or
\code{FALSE}.
}
\section{Missing values and invalid input}{
Start times must not be missing or infinite and, with \code{end} or
\code{duration}, neither may end times or durations; an edge must not end
before it starts, nor be active outside the presence spell of either
endpoint. The engine reports every such problem at once, with the rows of
\code{edges} or \code{nodes}. Start and end columns must have the same
type. Everything else is checked as in \code{\link{net_frame}}.
}
\references{
Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}

Butts, C. T., Leslie-Cook, A., Krivitsky, P. N., & Bender-deMoll, S.
(2026). networkDynamic: Dynamic extensions for network objects. R package
version 0.12.0. \doi{10.32614/CRAN.package.networkDynamic}
}
\seealso{
\code{\link{net_at}}, \code{\link{net_window}},
\code{\link{net_snapshots}}, \code{\link{net_aggregate}},
\code{\link{net_temporal_distances}}, \code{\link{net_temporal_reach}},
\code{\link{net_temporal_closeness}}, \code{\link{net_rolling}},
\code{\link{net_burstiness}}, and \code{\link{net_panel}} for panel waves.
}
\examples{
# Instantaneous contacts with numeric times.
msgs <- data.frame(source = c("a", "b", "b", "a"), target = c("b", "c", "c", "c"),
                   time = c(1, 1, 3, 6))
g <- net_temporal(msgs)
g
net_is_temporal(g)
# Every contact is an edge for static measures.
net_degree(g)

# Interval contacts with date-times, and node presence intervals.
calls <- data.frame(
  caller = c("a", "b", "a"), callee = c("b", "c", "c"),
  began = as.POSIXct(c("2024-05-01 10:00", "2024-05-01 10:20", "2024-05-01 11:00"),
                     tz = "UTC"),
  ended = as.POSIXct(c("2024-05-01 10:05", "2024-05-01 10:30", "2024-05-01 11:30"),
                     tz = "UTC"))
people <- data.frame(id = c("a", "b", "c"),
                     joined = as.POSIXct(c("2024-05-01 09:00", "2024-05-01 10:10", NA),
                                         tz = "UTC"))
h <- net_temporal(calls, people, start = "began", end = "ended", node_start = "joined",
                  directed = FALSE, source = "caller", target = "callee")
h
class(h)
}
