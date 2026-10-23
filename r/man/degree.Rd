\name{degree}
\alias{degree}
\alias{strength}
\title{Degree and strength}
\description{
\code{degree()} counts edge endpoints at each node; it is the same for
weighted and unweighted networks. \code{strength()} sums the weights of those
edges instead (Barrat et al. 2004). \code{degree(x, weighted = TRUE)} is the
same as \code{strength(x)}.
}
\usage{
degree(x, mode = c("all", "out", "in"), loops = TRUE, weighted = FALSE,
       convention = c("igraph", "networkx", "sna", "ucinet"))

strength(x, mode = c("all", "out", "in"), loops = TRUE)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks, \code{"all"}, \code{"out"}, or
\code{"in"}. Ignored for undirected networks.}
\item{loops}{Whether self-loops count.}
\item{weighted}{Whether to sum edge weights.}
\item{convention}{\code{"igraph"} (the default) and \code{"networkx"} count
every edge; \code{"sna"} and \code{"ucinet"} count distinct ties between
distinct nodes.}
}
\details{
Parallel edges each count. In undirected networks a self-loop adds 2 to its
node; in directed networks it adds 1 to out-degree and 1 to in-degree. These
match igraph and NetworkX. Set \code{loops = FALSE} to ignore self-loops,
which matches the sna default (\code{diag = FALSE}).

\code{convention = "sna"} or \code{"ucinet"} reads the network as an
adjacency matrix, as those tools do: parallel edges collapse into one tie
(keeping the largest weight) and self-loops are ignored, whatever
\code{loops} says. See \code{\link{conventions}}.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004). The
architecture of complex weighted networks. \emph{Proceedings of the National
Academy of Sciences}, 101(11), 3747-3752. \doi{10.1073/pnas.0400087101}
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 2), to = c(2, 3, 3)))
degree(g, mode = "out")

# A self-loop on node 1 and a repeated edge 1 -> 2: the conventions differ
# in whether loops and parallel edges count.
h <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(1, 2, 2, 3)))
degree(h)                          # igraph: loops count twice, parallel edges each count
degree(h, convention = "networkx") # NetworkX
degree(h, convention = "sna")      # sna: adjacency matrix, loops ignored
degree(h, convention = "ucinet")   # UCINET: adjacency matrix, diagonal ignored
}
