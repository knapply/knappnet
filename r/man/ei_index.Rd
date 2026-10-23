\name{ei_index}
\alias{ei_index}
\title{E-I index}
\usage{
ei_index(x, groups, scope = c("global", "group", "node"),
         direction = c("all", "out", "in"),
         convention = c("krackhardt", "homophily"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{groups}{Group of each node: the name of a column of \code{x$nodes}, or
a vector with one value per node. Missing values are not allowed.}
\item{scope}{\code{"global"} for the whole network, \code{"group"} for each
group, or \code{"node"} for each node.}
\item{direction}{For groups and nodes of directed networks: count ties sent
(\code{"out"}), received (\code{"in"}), or both (\code{"all"}).}
\item{convention}{How group-level ties are counted: \code{"krackhardt"} (see
Details) or \code{"homophily"}, which reproduces the R package homophily.}
}
\value{For \code{"global"}, a number. For \code{"group"}, a vector named by
group in order of first appearance. For \code{"node"}, a vector aligned with
the rows of \code{x$nodes}. External and internal tie counts are attached as
attributes \code{"external"} and \code{"internal"}.}
\description{
Krackhardt and Stern's (1988) E-I index compares the ties that cross group
boundaries (external, E) with those inside groups (internal, I):
\deqn{\frac{E - I}{E + I}}{(E - I) / (E + I)}
It ranges from -1 (every tie inside a group) to 1 (every tie between groups).
}
\details{
Network: E is the number of ties between nodes of different groups and I the
number within groups. Group (Krackhardt and Stern 1988): I is the number of
ties among the group's members and E the number of ties between a member and a
non-member (with \code{direction = "out"} or \code{"in"}, only those sent or
received by members). Node: the node's own ties to members of its group (I)
and to others (E); in directed networks ties sent and received both count with
\code{direction = "all"}, so a mutual dyad counts twice.

Ties are binary: weights are ignored and each parallel edge counts as a
separate tie (a warning is given when there are parallel edges, as homophily
does). A self-loop is an internal tie. An index is \code{NaN} where there are
no ties.

\code{convention = "homophily"} matches homophily's \code{ei_index()}. Network
and node indices agree with the default except that homophily counts an
undirected self-loop twice at its node. For groups, homophily takes I from the
diagonal of its folded mixing matrix (see \code{\link{mixing_matrix}}) and E
from the group's row of the upper triangle, so a tie between two groups counts
as external only for the group that appears first in \code{x$nodes}. UCINET
also reports the index for groups and nodes; its documentation could not be
consulted, so agreement with UCINET is not claimed for those levels.
}
\references{
Krackhardt, D., & Stern, R. N. (1988). Informal networks and organizational
crises: An experimental simulation. \emph{Social Psychology Quarterly}, 51(2),
123-140. \doi{10.2307/2786835}

Knapp, B. G. (2019). homophily: Measuring network homophily. R package,
\url{https://github.com/knapply/homophily}.
}
\seealso{
\code{\link{mixing_matrix}}, \code{\link{homophily}}
}
\examples{
staff <- data.frame(name = c("ann", "bob", "cat", "dan", "eve"),
                    dept = c("sales", "sales", "eng", "eng", "hr"))
g <- netframe(data.frame(from = c("ann", "ann", "bob", "cat", "dan"),
                         to = c("bob", "cat", "dan", "dan", "eve")),
              nodes = staff, directed = FALSE)
ei_index(g, "dept")
ei_index(g, "dept", scope = "node")

# Group level: the two conventions count ties between groups differently.
ei_index(g, "dept", scope = "group")                           # Krackhardt and Stern
ei_index(g, "dept", scope = "group", convention = "homophily") # homophily
}
