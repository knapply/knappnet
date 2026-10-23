\name{closeness}
\alias{closeness}
\alias{harmonic_centrality}
\title{Closeness and harmonic centrality}
\description{
Closeness is the reciprocal of the sum of a node's shortest-path lengths to
the other nodes (Sabidussi 1966; Freeman 1978). Harmonic centrality sums the
reciprocals of those lengths instead, so unreachable nodes simply contribute
0 (Marchiori and Latora 2000; Rochat 2009; Boldi and Vigna 2014).
}
\usage{
closeness(x, mode = NULL, weighted = NULL, normalized = NULL, threads = 0L,
          weights_as = NULL, convention = c("igraph", "networkx", "sna", "ucinet"))

harmonic_centrality(x, mode = NULL, weighted = NULL, normalized = NULL,
                    threads = 0L, weights_as = NULL,
                    convention = c("igraph", "networkx", "sna", "ucinet"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks: \code{"out"} uses paths from the node,
\code{"in"} paths to it, and \code{"all"} ignores direction. \code{NULL}
uses the convention's default: \code{"in"} for NetworkX and \code{"out"}
otherwise.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{normalized}{Whether to normalize (see Details). \code{NULL} uses the
convention's default: unnormalized under igraph's convention and for
NetworkX's harmonic centrality, normalized otherwise. NetworkX's closeness
exists only in normalized form.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
\item{weights_as}{Override the network's \code{weight_type}:
\code{"distance"} or \code{"strength"}.}
\item{convention}{Whose definition to follow; see Details and
\code{\link{conventions}}.}
}
\details{
Let \eqn{n} be the number of nodes, \eqn{r} the number of other nodes a node
reaches, and \eqn{S} the sum of the lengths to them. Tools differ in how they
handle nodes that cannot be reached:
\describe{
  \item{\code{"igraph"}}{Only reachable nodes count: \eqn{1/S}, or \eqn{r/S}
    normalized. A node that reaches no other gets \code{NaN}.}
  \item{\code{"networkx"}}{The Wasserman and Faust (1994, p. 201) index
    \eqn{(r/S)(r/(n-1))}, which is always normalized; 0 if \eqn{r = 0}.}
  \item{\code{"sna"}}{\eqn{(n-1)/S} (or \eqn{1/S} unnormalized), and 0 for
    any node that does not reach every other node.}
  \item{\code{"ucinet"}}{\eqn{(n-1)/S'}, where \eqn{S'} counts every missing
    length as one more than the longest shortest path in the network,
    UCINET's default ("max observed distance plus 1").}
}
All four agree on strongly connected networks. With the sna and UCINET
conventions the network is read as an adjacency matrix (parallel edges
collapse, self-loops are ignored) and weights are not used, as in those
tools.

Harmonic centrality is \eqn{\sum 1/d}{sum 1/d}, divided by \eqn{n-1} when normalized;
this is igraph's \code{harmonic_centrality()}, NetworkX's
\code{harmonic_centrality()} (never normalized, paths to the node), sna's
\code{closeness(cmode = "suminvdir")}, and UCINET's reciprocal closeness
(both normalized). The convention only changes the defaults of \code{mode}
and \code{normalized} for harmonic centrality.

Self-loops never lie on a shortest path, and parallel edges do not change
path lengths (with weights the shortest one counts). Weights are handled as
in \code{\link{betweenness}}: distances, or strengths converted to lengths
1/weight, and they must be positive. When every weight is the same, lengths
are that weight times the number of steps. Searches from each node run in
parallel.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Sabidussi, G. (1966). The centrality index of a graph. \emph{Psychometrika},
31(4), 581-603. \doi{10.1007/BF02289527}

Freeman, L. C. (1978). Centrality in social networks: Conceptual
clarification. \emph{Social Networks}, 1(3), 215-239.
\doi{10.1016/0378-8733(78)90021-7}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods and
Applications}. Cambridge University Press. \doi{10.1017/CBO9780511815478}

Marchiori, M., & Latora, V. (2000). Harmony in the small-world.
\emph{Physica A}, 285(3-4), 539-546. \doi{10.1016/S0378-4371(00)00311-3}

Rochat, Y. (2009). Closeness centrality extended to unconnected graphs: The
harmonic centrality index. 6th Conference on Applications of Social Network
Analysis (ASNA 2009), Zurich.

Boldi, P., & Vigna, S. (2014). Axioms for centrality. \emph{Internet
Mathematics}, 10(3-4), 222-262. \doi{10.1080/15427951.2013.865686}
}
\examples{
# A path 1-2-3 and a separate pair 4-5: each tool treats unreachable
# nodes differently.
g <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)), directed = FALSE)
closeness(g)                          # igraph: reachable nodes only
closeness(g, convention = "networkx") # NetworkX: Wasserman-Faust scaling
closeness(g, convention = "sna")      # sna: 0 unless every node is reached
closeness(g, convention = "ucinet")   # UCINET: missing distances = longest + 1

harmonic_centrality(g)                          # igraph: not normalized
harmonic_centrality(g, convention = "networkx") # NetworkX: not normalized
harmonic_centrality(g, convention = "sna")      # sna: normalized
harmonic_centrality(g, convention = "ucinet")   # UCINET: normalized
}
