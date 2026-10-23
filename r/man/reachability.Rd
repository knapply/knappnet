\name{reachability}
\alias{reachability}
\title{Reachability}
\description{
The number of other nodes each node can reach along a path.
}
\usage{
reachability(x, mode = c("out", "in", "all"), threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks: \code{"out"} counts the nodes a node
reaches, \code{"in"} the nodes that reach it, and \code{"all"} the nodes
joined to it when directions are ignored (the size of its weak component,
less one). Ignored for undirected networks.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
Weights, parallel edges, and self-loops do not matter. Undirected networks
and \code{mode = "all"} take time linear in the size of the network;
otherwise one breadth-first search runs from every node, in parallel. The
full reachability matrix is
\code{is.finite(distances(x, mode = mode, weighted = FALSE))}.

Results match igraph's \code{ego_size()} with an unlimited order (less one)
and NetworkX's \code{descendants()}, \code{ancestors()}, and
\code{node_connected_component()} on every network in the validation suite.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\examples{
g <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)))
reachability(g)
reachability(g, "in")
reachability(g, "all")
}
