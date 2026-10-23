\name{net_node_spells}
\alias{net_node_spells}
\alias{net_nodes_at}
\title{Node presence spells and node attributes over time}
\description{
In a temporal network a node may be present only during some spells (it
joins, leaves, and may return), and its attributes may change at known
times. One spell per node is given by the \code{node_start} and
\code{node_end} columns of \code{\link{net_frame}}; several spells per
node by its long \code{node_spells} table; attribute changes by its long
\code{node_changes} table. \code{net_node_spells()} returns every node's
presence spells, and \code{net_nodes_at()} returns the node table with
the attribute values in force at a given time.
}
\usage{
net_node_spells(x)

net_nodes_at(x, time)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{time}{One time, of the class of the network's times (a number,
\code{Date}, or \code{POSIXct}).}
}
\details{
A spell is the half-open interval \eqn{[start, end)}; a missing start or
end leaves it open on that side, and a node with no spell is present
throughout. Two spells of one node may not overlap (\eqn{[0, 5)} and
\eqn{[5, 9)} do not) unless \code{overlapping_spells = TRUE} is declared.
Every edge spell or contact must lie within a presence spell of each
endpoint; violating rows are listed in the error. This is the vertex
activity of networkDynamic (Butts et al.) and the composition change of
RSiena (Ripley et al.).

A row of \code{node_changes} sets the attribute columns it carries for its
node from its time on (inclusive); the node table holds the values before
the first change. \code{net_nodes_at()} applies every change at or before
\code{time} in time order; a column that appears only in
\code{node_changes} is added, missing before its first change.
}
\value{
\code{net_node_spells()}: a data frame with the identifier column and
\code{start} and \code{end} (of the class of the network's times;
\code{-Inf} and \code{Inf} for open sides), one row per spell, ordered by
node and start. \code{net_nodes_at()}: a data frame like \code{x$nodes}.
}
\references{
Butts, C. T., Leslie-Cook, A., Krivitsky, P. N., & Bender-deMoll, S.
(2026). \emph{networkDynamic: Dynamic Extensions for Network Objects}. R
package version 0.12.0. \doi{10.32614/CRAN.package.networkDynamic}

Ripley, R. M., Snijders, T. A. B., Boda, Z., & Voros, A. (2026).
\emph{Manual for RSiena}. University of Oxford.
\url{https://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf}
}
\seealso{
\code{\link{net_temporal}}, \code{\link{net_frame}}, and
\code{\link{net_plot_snapshots}}, which greys out absent nodes.
}
\examples{
nodes <- data.frame(id = c("ann", "bob", "cy"), status = "guest")
spells <- data.frame(id = c("ann", "ann", "bob"), start = c(0, 8, 0), end = c(5, 20, 20))
changes <- data.frame(id = c("ann", "bob"), time = c(8, 12), status = c("member", "staff"))
edges <- data.frame(source = c("ann", "ann", "bob"), target = c("bob", "bob", "cy"),
                    from = c(1, 12, 3), to = c(4, 15, 3))
g <- net_frame(edges, nodes, directed = FALSE, start = "from", end = "to",
               node_spells = spells, node_changes = changes)
g
net_node_spells(g)
net_nodes_at(g, 10)
net_nodes_at(g, 13)

# An edge while a node is absent is an error that names the row.
edges$to[1] <- 10
try(net_frame(edges, nodes, directed = FALSE, start = "from", end = "to",
              node_spells = spells))

# Absent nodes are drawn grey in each window.
net_plot_snapshots(g, breaks = c(0, 5, 10, 20), seed = 1)
}
