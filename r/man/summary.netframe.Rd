\name{summary.netframe}
\alias{summary.netframe}
\alias{print.summary.netframe}
\title{Summarize a network}
\usage{
\method{summary}{netframe}(object, paths = NULL, weighted = NULL, threads = 0L, ...)

\method{print}{summary.netframe}(x, digits = 3L, ...)
}
\arguments{
\item{object}{A \code{\link{netframe}}.}
\item{paths}{Whether to compute the mean distance and diameter, which need a
shortest-path search from every node of the largest component. \code{NULL}
(the default) computes them when that component has at most 5,000 nodes and
otherwise skips them, saying so.}
\item{weighted}{Whether the mean distance and diameter use edge weights
(which must then be positive). \code{NULL} uses them when the network has a
weight column; see \code{\link{distances}}.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
\item{...}{Not used.}
\item{x}{A \code{summary.netframe} object.}
\item{digits}{Significant digits for proportions and means.}
}
\value{
\code{summary()} returns an object of class \code{summary.netframe}: a list
holding the values described in Details, with a \code{print()} method.
\code{print()} returns its argument invisibly.
}
\description{
The standard descriptive table of a network: its size, ties, components,
cohesion, distances, and degree distribution.
}
\details{
The summary reports:
\describe{
  \item{Nodes, edges}{Rows of \code{x$nodes} and \code{x$edges}.}
  \item{Self-loops, multiple edges}{Edges from a node to itself, and edges
    repeating an earlier edge between the same (ordered, if directed) pair
    of nodes, as igraph's \code{sum(which_multiple())}.}
  \item{Isolates}{Nodes without an edge to another node.}
  \item{Density}{\code{\link{graph_density}()}.}
  \item{Components}{The number of weak components
    (\code{\link{components}()}) and the size of the largest.}
  \item{Reciprocity}{\code{\link{reciprocity}()}, for directed networks.}
  \item{Transitivity}{\code{\link{transitivity}()}.}
  \item{Mean distance, diameter}{\code{\link{mean_distance}()} and
    \code{\link{diameter}()} with \code{unconnected = TRUE}, following edge
    directions. When the network has more than one component they are
    computed on the largest weak component (the first, if several are
    equally large), and the printout says so.}
  \item{Degree}{The mean, smallest, and largest \code{\link{degree}()} (all
    edges; a self-loop counts twice in an undirected network).}
  \item{Degree centralization}{\code{\link{centralization}()} of degree with
    its defaults.}
}
Everything except the mean distance and diameter takes time roughly linear
in the size of the network, so networks with hundreds of thousands of nodes
are summarized quickly.
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 2, 3, 5), to = c(2, 3, 3, 4, 6)),
              directed = FALSE)
summary(g)
s <- summary(g)
s$transitivity
}
