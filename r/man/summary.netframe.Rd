\name{summary.netframe}
\alias{summary.netframe}
\alias{print.summary.netframe}
\title{Summarize a network}
\description{
The standard descriptive table of a network: its size, ties, components,
cohesion, distances, and degree distribution.
}
\usage{
\method{summary}{netframe}(object, paths = NULL, weighted = NULL, threads = 0L, ...)

\method{print}{summary.netframe}(x, digits = 3L, ...)
}
\arguments{
\item{object}{A \code{\link{net_frame}}.}
\item{paths}{Whether to compute the mean distance and diameter, which need a
shortest-path search from every node of the largest component. \code{NULL}
(the default) computes them when that component has at most 5,000 nodes and
otherwise skips them, saying so.}
\item{weighted}{Whether the mean distance and diameter use edge weights (which
must then be positive). \code{NULL} (the default) uses them when the network
has a weight column; see \code{\link{net_distances}}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
\item{...}{Not used.}
\item{x}{A \code{summary.netframe} object.}
\item{digits}{Significant digits for proportions and means.}
}
\details{
The summary reports:
\describe{
  \item{Nodes, edges}{Rows of \code{x$nodes} and \code{x$edges}.}
  \item{Self-loops, multiple edges}{Edges from a node to itself, and edges
    repeating an earlier edge between the same (ordered, if directed) pair
    of nodes, as igraph's \code{sum(which_multiple())}.}
  \item{Isolates}{Nodes without an edge to another node.}
  \item{Density}{\code{\link{net_density}()}.}
  \item{Components}{The number of weak components
    (\code{\link{net_components}()}) and the size of the largest.}
  \item{Reciprocity}{\code{\link{net_reciprocity}()}, for directed networks.}
  \item{Transitivity}{\code{\link{net_transitivity}()}.}
  \item{Mean distance, diameter}{\code{\link{net_mean_distance}()} and
    \code{\link{net_diameter}()} with \code{unconnected = TRUE}, following edge
    directions. When the network has more than one component they are
    computed on the largest weak component (the first, if several are
    equally large), and the printout says so.}
  \item{Degree}{The mean, smallest, and largest \code{\link{net_degree}()} (all
    edges; a self-loop counts twice in an undirected network).}
  \item{Degree centralization}{\code{\link{net_centralization}()} of degree with
    its defaults.}
}
Everything except the mean distance and diameter takes time roughly linear
in the size of the network, so networks with hundreds of thousands of nodes
are summarized quickly.
}
\value{
\code{summary()} returns an object of class \code{summary.netframe}: a list
holding the values described in Details, with a \code{print()} method.
\code{print()} returns its argument invisibly.
}
\examples{
g <- net_frame(data.frame(source = c(1, 1, 2, 3, 5), target = c(2, 3, 3, 4, 6)),
              directed = FALSE)
summary(g)
s <- summary(g)
s$transitivity
}
