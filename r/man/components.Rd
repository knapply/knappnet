\name{components}
\alias{components}
\alias{fragmentation}
\title{Connected components and fragmentation}
\usage{
components(x, mode = c("weak", "strong"))

fragmentation(x, mode = c("weak", "strong"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{\code{"weak"} or \code{"strong"}.}
}
\value{
\code{components()} returns a list with \code{membership} (the component of
each node, aligned with the rows of \code{x$nodes}), \code{csize} (the size
of each component), and \code{no} (the number of components), as igraph
does. \code{fragmentation()} returns a number, \code{NaN} for networks with
fewer than two nodes.
}
\description{
\code{components()} finds the weak or strong components of a network.
\code{fragmentation()} is the proportion of ordered pairs of distinct nodes
that lie in different components (Borgatti 2006):
\deqn{F = 1 - \frac{\sum_k s_k (s_k - 1)}{n (n - 1)}}{F = 1 - sum_k s_k (s_k - 1) / (n (n - 1))}
where \eqn{s_k} are the component sizes.
}
\details{
In a weak component every node can reach every other when the direction of
edges is ignored; in a strong component every node reaches every other along
directed paths. Undirected networks have one kind of component, whichever
\code{mode} is chosen. Isolated nodes are components of size 1. Weights,
parallel edges, and self-loops do not matter.

Strong components are found with Tarjan's (1972) algorithm, written without
recursion, in time linear in the size of the network; weak components with
union-find over the edge list (Tarjan 1975), in near-linear time. Components are numbered in order of their
first node in \code{x$nodes}. Membership matches igraph's
\code{components()} and NetworkX's component functions up to the numbering.

Fragmentation of weak components equals UCINET's normalized heterogeneity of
component sizes and its "Fragmentation" for undirected data.
}
\references{
Tarjan, R. (1972). Depth-first search and linear graph algorithms.
\emph{SIAM Journal on Computing}, 1(2), 146-160. \doi{10.1137/0201010}

Tarjan, R. E. (1975). Efficiency of a good but not linear set union
algorithm. \emph{Journal of the ACM}, 22(2), 215-225.

Borgatti, S. P. (2006). Identifying sets of key players in a social network.
\emph{Computational and Mathematical Organization Theory}, 12(1), 21-34.
\doi{10.1007/s10588-006-7084-x}
}
\examples{
g <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)), directed = FALSE)
components(g)
fragmentation(g)
}
