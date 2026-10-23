\name{net_plot_timeline}
\alias{net_plot_timeline}
\title{Plot the timeline of contacts}
\description{
Draws when each pair of nodes (or each node) is in contact: one row per
pair, with a tick for every instantaneous contact and a bar for every
interval contact, the raw-data view of a contact sequence (Holme and
Saramäki 2012, Fig. 1).
}
\usage{
net_plot_timeline(x, by = c("pair", "node"), max_rows = 50, ...)
}
\arguments{
\item{x}{A temporal network (\code{\link{net_temporal}}).}
\item{by}{\code{"pair"} for one row per pair of nodes (ordered pairs in a
directed network) or \code{"node"} for one row per node, on which every
contact of the node appears.}
\item{max_rows}{How many rows to draw at most: the pairs or nodes with the
most contacts are kept and the rest are counted in the subtitle.}
\item{...}{Further arguments passed to \code{\link[graphics]{plot}}, such
as \code{main} or \code{xlab}.}
}
\details{
Rows are sorted by number of contacts, the most active at the top. The
time axis is labelled in the type of the network's times (numbers, dates,
or date-times). Interval contacts are drawn from their start to their end;
an instantaneous contact is a tick at its time.
}
\value{A data frame of the rows drawn, from the top, invisibly: \code{row}
(the pair or node label) and \code{contacts}.}
\references{
Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
}
\seealso{\code{\link{net_temporal}}, \code{\link{net_plot_rolling}},
\code{\link{net_burstiness}}.}
\examples{
calls <- data.frame(source = c("a", "a", "b", "c", "a", "b"),
                    target   = c("b", "c", "c", "d", "b", "c"),
                    time = c(1, 2, 2, 5, 7, 9), end = c(3, 2, 4, 6, 7, 10))
g <- net_temporal(calls, end = "end", directed = FALSE)
net_plot_timeline(g)
net_plot_timeline(g, by = "node", main = "Contacts per node")
}
