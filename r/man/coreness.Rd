\name{coreness}
\alias{coreness}
\title{Coreness (k-cores)}
\usage{
coreness(x, mode = c("all", "out", "in"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks, \code{"all"} (in-degree plus
out-degree), \code{"out"}, or \code{"in"}. Ignored for undirected
networks.}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
The coreness (core number) of a node is the largest \eqn{k} such that the
node belongs to the \eqn{k}-core: the largest subgraph in which every node
has degree at least \eqn{k} (Seidman 1983). It is computed with the
\eqn{O(m)} algorithm of Batagelj and Zaversnik (2003).
}
\details{
In a directed network, \code{mode = "out"} builds cores on out-degree,
\code{"in"} on in-degree, and \code{"all"} on their sum, so a mutual pair of
arcs counts twice.

Parallel edges each count. A self-loop adds 2 to its node's degree in an
undirected network, and 1 to out-degree and 1 to in-degree in a directed
one; a loop stays inside every subgraph containing its node. Weights are
ignored: cores are defined on binary ties. Isolated nodes have coreness 0.

These conventions match \code{igraph::coreness()}. NetworkX's
\code{core_number()} rejects self-loops and multigraphs; on simple
networks without loops all three agree.
}
\references{
Seidman, S. B. (1983). Network structure and minimum degree.
\emph{Social Networks}, 5(3), 269-287. \doi{10.1016/0378-8733(83)90028-X}

Batagelj, V., & Zaversnik, M. (2003). An O(m) algorithm for cores
decomposition of networks. arXiv preprint cs/0310049.
\doi{10.48550/arXiv.cs/0310049}
}
\examples{
g <- netframe(data.frame(from = c(1, 2, 3, 3), to = c(2, 3, 1, 4)),
              directed = FALSE)
coreness(g)
}
