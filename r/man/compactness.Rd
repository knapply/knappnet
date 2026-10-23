\name{compactness}
\alias{compactness}
\alias{distance_weighted_fragmentation}
\title{Compactness and distance-weighted fragmentation}
\usage{
compactness(x, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
            threads = 0L)

distance_weighted_fragmentation(x, mode = c("out", "in", "all"), weighted = NULL,
                                weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks: \code{"out"} uses paths from each node,
\code{"in"} paths to it, and \code{"all"} ignores direction.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{weights_as}{How to read weights: \code{"distance"} (larger keeps nodes
farther apart) or \code{"strength"} (larger brings them closer; lengths
\eqn{1/w}). \code{NULL} uses the network's \code{weight_type}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A number.}
\description{
Borgatti (2006) measures how far a network is from falling apart with the
distance-weighted fragmentation
\deqn{{}^{D}F = 1 - \frac{\sum_{i \ne j} 1 / d_{ij}}{n (n - 1)}}{DF = 1 - sum_{i != j} (1 / d_ij) / (n (n - 1))}
where \eqn{1 / d_{ij} = 0} for nodes that cannot reach each other. Its
complement, the mean reciprocal distance over ordered pairs of distinct
nodes, is the network's compactness (the global efficiency of Latora and
Marchiori 2001). Fragmentation is 0 for a complete network and 1 for a
network without ties; unlike \code{\link{fragmentation}} it distinguishes
networks that are connected but stretched out. UCINET reports both in
Network > Cohesion > Distance.
}
\details{
Directed networks use directed distances over ordered pairs (\code{mode}
\code{"out"}: from \eqn{i} to \eqn{j}; \code{"in"}: to each node;
\code{"all"}: ignoring direction). Distances are counted in steps, or, with
weights, as sums of lengths as in \code{\link{distances}} (weights read as
distances or as strengths converted to lengths \eqn{1/w}; they must be
strictly positive). Self-loops never shorten a path and parallel edges
count through the shortest. Networks with fewer than two nodes give
\code{NaN}. Compactness matches NetworkX's \code{global_efficiency()} on
undirected unweighted networks and the mean of \code{1 / igraph::distances()}
off the diagonal otherwise.
}
\references{
Borgatti, S. P. (2006). Identifying sets of key players in a social
network. \emph{Computational and Mathematical Organization Theory}, 12(1),
21-34. \doi{10.1007/s10588-006-7084-x}

Latora, V., & Marchiori, M. (2001). Efficient behavior of small-world
networks. \emph{Physical Review Letters}, 87(19), 198701.
\doi{10.1103/PhysRevLett.87.198701}
}
\seealso{\code{\link{fragmentation}}, \code{\link{distances}}}
\examples{
# A path 1 - 2 - 3 and an isolate 4.
g <- netframe(data.frame(from = c(1, 2), to = c(2, 3)), nodes = data.frame(id = 1:4),
              directed = FALSE)
compactness(g)                       # (4 * 1 + 2 * 1/2) / 12
distance_weighted_fragmentation(g)
fragmentation(g)                     # plain fragmentation only sees components

# Direction matters: a directed path reaches fewer ordered pairs.
d <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
compactness(d)                       # out
compactness(d, mode = "all")         # ignoring direction
}
