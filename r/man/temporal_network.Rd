\name{temporal_network}
\alias{temporal_network}
\alias{is_temporal}
\title{Temporal networks}
\usage{
temporal_network(edges, nodes = NULL, start = "time", end = NULL,
  duration = NULL, directed = TRUE, from = 1L, to = 2L, id = 1L,
  node_start = NULL, node_end = NULL, weight = NULL,
  weight_type = c("distance", "strength"), graph = list())

is_temporal(x)
}
\arguments{
\item{edges}{A data frame of contacts, one row per contact, or a
\code{\link{netframe}} whose edge table holds the times; its nodes,
endpoints, weights, and other roles are kept.}
\item{nodes}{A data frame of nodes, or \code{NULL}; see
\code{\link{netframe}}. Must be \code{NULL} when \code{edges} is a
netframe.}
\item{start}{Column of \code{edges} holding each contact's start time, by
name or position. Times may be numbers (in any unit), \code{Date}, or
\code{POSIXct} date-times. By default the column named \code{"time"}, or,
if there is none, the column named \code{"start"}; a netframe that already
has times keeps them.}
\item{end}{Optional column of end times, of the same type as
\code{start}.}
\item{duration}{Optional column of durations, an alternative to
\code{end}: numbers in the unit of the times (seconds for date-times, days
for dates) or a \code{difftime}. Without \code{end} or \code{duration},
every contact is instantaneous.}
\item{directed}{Whether contacts are directed. Ignored when \code{edges}
is a netframe, which keeps its own direction.}
\item{from, to, id, weight, weight_type, graph}{As in
\code{\link{netframe}}.}
\item{node_start, node_end}{Optional columns of \code{nodes} holding when
each node is present, of the same type as \code{start}. A missing value
means present from the beginning (start) or until the end (end).}
\item{x}{Any object.}
}
\value{
\code{temporal_network()} returns a \code{\link{netframe}} with time
roles; \code{is_temporal()} returns \code{TRUE} or \code{FALSE}.
}
\description{
\code{temporal_network()} creates a network whose edges are timed
contacts: each edge has a start time and, optionally, an end time or a
duration, and nodes may have presence intervals. The result is an ordinary
\code{\link{netframe}} with extra column roles, so every static measure
still works on it (treating each contact as an edge), and the temporal
functions listed under \sQuote{See also} use the times.
\code{is_temporal()} tells whether a netframe has times.
}
\details{
\strong{Times.} Times are stored as given in \code{x$edges} and
\code{x$nodes}; the column roles are recorded in \code{x$roles}
(\code{start}, \code{end} or \code{duration}, \code{node_start},
\code{node_end}, and the time type). Internally dates and date-times are
converted to seconds since 1970-01-01 UTC, so a date is its midnight in
UTC. Arguments that are times (such as \code{from} and \code{to} of
\code{\link{network_window}}) take the type of the data: numbers for
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
same nodes are allowed and kept; each contact is one edge of the netframe.
Static measures applied directly to a temporal network see every contact
as an edge (a multigraph); use \code{\link{aggregate_network}} for one
weighted edge per pair of nodes, or \code{\link{network_window}} and
\code{\link{snapshots}} for the contacts of a period.

\strong{Log data.} Tables of events such as authentication logs or network
flows (one row per event, with source, destination, and time) are contact
sequences: pass them, or a netframe built from them, straight to
\code{temporal_network()}. Millions of contacts are fine; slicing and
paths run in compiled code after a single sort by time.
}
\section{Missing values and invalid input}{
Start times must not be missing and, with \code{end} or \code{duration},
neither may end times or durations; ends must not precede starts. Errors
name the rows of \code{edges}. Start and end columns must have the same
type. Everything else is checked as in \code{\link{netframe}}.
}
\references{
Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}

Butts, C. T., Leslie-Cook, A., Krivitsky, P. N., & Bender-deMoll, S.
(2026). networkDynamic: Dynamic extensions for network objects. R package
version 0.12.0. \doi{10.32614/CRAN.package.networkDynamic}
}
\seealso{
\code{\link{network_at}}, \code{\link{network_window}},
\code{\link{snapshots}}, \code{\link{aggregate_network}},
\code{\link{temporal_distances}}, \code{\link{temporal_reach}},
\code{\link{temporal_closeness}}, \code{\link{rolling}},
\code{\link{burstiness}}.
}
\examples{
# Instantaneous contacts with numeric times.
msgs <- data.frame(from = c("a", "b", "b", "a"), to = c("b", "c", "c", "c"),
                   time = c(1, 1, 3, 6))
g <- temporal_network(msgs)
g
is_temporal(g)
# Every contact is an edge for static measures.
degree(g)

# Interval contacts with date-times, and node presence intervals.
calls <- data.frame(
  caller = c("a", "b", "a"), callee = c("b", "c", "c"),
  began = as.POSIXct(c("2024-05-01 10:00", "2024-05-01 10:20", "2024-05-01 11:00"),
                     tz = "UTC"),
  ended = as.POSIXct(c("2024-05-01 10:05", "2024-05-01 10:30", "2024-05-01 11:30"),
                     tz = "UTC"))
people <- data.frame(name = c("a", "b", "c"),
                     joined = as.POSIXct(c("2024-05-01 09:00", "2024-05-01 10:10", NA),
                                         tz = "UTC"))
h <- temporal_network(calls, people, start = "began", end = "ended",
                      node_start = "joined", directed = FALSE)
h$roles[c("start", "end", "node_start", "time_class")]
}
