\name{net_snapshots}
\alias{net_snapshots}
\title{Snapshots of a temporal network}
\description{
Splits a temporal network into a series of windows, sliding (\code{step}
smaller than \code{width}), tumbling (\code{step} equal to \code{width},
the default), or with gaps (\code{step} larger), or between given
\code{breaks}. Each window \eqn{[a, b)} holds the contacts active in it,
as in \code{\link{net_window}}.
}
\usage{
net_snapshots(x, width = NULL, step = width, breaks = NULL, from = NULL,
  to = NULL, nodes = c("all", "present", "active"))
}
\arguments{
\item{x}{A temporal network (\code{\link{net_temporal}}).}
\item{width}{Window length: a number in the unit of the times (seconds for
date-times, days for dates) or a \code{difftime}.}
\item{step}{Distance between window starts, as \code{width}.}
\item{breaks}{Alternatively to \code{width}, increasing times: windows are
\eqn{[breaks_k, breaks_{k+1})}.}
\item{from, to}{With \code{width}: the first window starts at \code{from}
(default: the earliest start time), and windows start up to and including
\code{to} (default: the latest end time), so with \code{step <= width}
every contact falls in some window.}
\item{nodes}{Which nodes each snapshot keeps: \code{"all"} (the default),
\code{"present"}, or \code{"active"}; see \code{\link{net_at}}.}
}
\details{
Contacts are sorted by time once and the whole series is computed in one
pass, with a heap of active contacts, in \eqn{O(m \log m)}{O(m log m)}
time plus the size of the output, so series over millions of contacts are
fast. An interval contact appears in every window it overlaps; an
instantaneous contact at \eqn{t} is in the windows with
\eqn{a \le t < b}{a <= t < b}. Snapshots keep every column and role of
\code{x}, with contacts in their original order and not clipped to the
window.
}
\value{
A list of temporal networks, one per window, with attribute
\code{"windows"}, a data frame of window \code{start} and \code{end}
times.
}
\seealso{
\code{\link{net_rolling}} to apply a measure to every snapshot,
\code{\link{net_window}} for one window.
}
\examples{
logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2"),
                     host = c("h1", "h1", "h2", "h2", "h3"),
                     time = c(0, 1, 5, 6, 12))
g <- net_temporal(logins)
s <- net_snapshots(g, width = 5)
attr(s, "windows")
vapply(s, function(w) nrow(w$edges), numeric(1))
# Sliding windows of width 6 every 3 time units.
vapply(net_snapshots(g, width = 6, step = 3), function(w) nrow(w$edges), numeric(1))
# Windows between given times.
attr(net_snapshots(g, breaks = c(0, 6, 13)), "windows")
}
