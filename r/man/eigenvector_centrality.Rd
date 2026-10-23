\name{eigenvector_centrality}
\alias{eigenvector_centrality}
\title{Eigenvector centrality}
\usage{
eigenvector_centrality(x, mode = NULL, weighted = NULL, scale = NULL,
                       convention = c("igraph", "networkx", "sna", "ucinet"),
                       tol = 1e-12, max_iter = 10000L, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks: \code{"in"} sums the scores of nodes with
ties \emph{to} the node (igraph with \code{directed = TRUE}, NetworkX),
\code{"out"} the nodes it has ties to (sna), and \code{"all"} ignores
direction (igraph's R default, \code{directed = FALSE}, and UCINET, which
symmetrizes). \code{NULL} uses the convention's choice.}
\item{weighted}{Whether to use edge weights, as given (larger means a
stronger tie), whatever the network's \code{weight_type}. \code{NULL} uses
them when the network has a weight column. They must not be negative.}
\item{scale}{\code{"max"} (largest score 1, igraph) or \code{"unit"}
(Euclidean norm 1, NetworkX, sna, and UCINET). \code{NULL} uses the
convention's choice.}
\item{convention}{Whose conventions to follow; see Details and
\code{\link{conventions}}.}
\item{tol}{Convergence tolerance: iteration stops when no score (scaled so
the largest is 1) changes by more than \code{tol}.}
\item{max_iter}{Iteration cap. Not converging within it is an error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Each node's score is proportional to the sum of its neighbors' scores: the
principal eigenvector of the adjacency matrix (Bonacich 1972, 1987).
}
\details{
Scores are computed by power iteration on \eqn{A + I} (von Mises and
Pollaczek-Geiringer 1929), starting from equal scores. The shift by the
identity changes no eigenvector but guarantees convergence on bipartite and
other periodic networks, where plain iteration on \eqn{A} cycles. Each
iteration costs time linear in the size of the network and runs in parallel
on large networks.

Parallel edges add up. In undirected networks a self-loop adds 2 to its
diagonal entry (two edge endpoints) under the igraph convention and 1 under
NetworkX's; directed loops add 1. The sna and UCINET conventions read the
network as an adjacency matrix: parallel edges collapse (keeping the largest
weight) and self-loops are ignored.

The eigenvector is unique (up to scale) and positive only when the network
is connected (strongly connected, for directed networks with \code{mode}
\code{"in"} or \code{"out"}); otherwise a warning is given. For undirected
networks (and \code{mode = "all"}) each component is then computed
separately: components whose largest eigenvalue is below the network's score
0, and components that share it each contribute their unit eigenvector
\eqn{u} times \eqn{\sum u}{sum(u)}, the limit of power iteration from equal
scores. For directed networks the result is the limit of the iteration; if
components sharing the largest eigenvalue are joined by a directed path, the
iteration does not converge and an error is reported. Convergence is also
slow when the two largest eigenvalues are close. NetworkX refuses
disconnected networks, and so does the \code{"networkx"} convention. When more than
one component has ties, UCINET reports a different measure (lambda-squared
scores computed within components), which is not implemented here, so the
\code{"ucinet"} convention stops.

Networks without cycles have no positive eigenvalue. Following igraph, nodes
whose score feeds no other node's score (nodes without outgoing ties for
\code{mode = "in"}, without incoming ties for \code{"out"}, and every node in
a network without edges) then score 1 before scaling and all others 0.

Results match igraph's \code{eigen_centrality()}, NetworkX's
\code{eigenvector_centrality_numpy()}, and sna's \code{evcent()} to within
their convergence tolerances on connected networks.
}
\references{
Bonacich, P. (1972). Factoring and weighting approaches to status scores and
clique identification. \emph{Journal of Mathematical Sociology}, 2(1),
113-120. \doi{10.1080/0022250X.1972.9989806}

Bonacich, P. (1987). Power and centrality: A family of measures.
\emph{American Journal of Sociology}, 92(5), 1170-1182. \doi{10.1086/228631}

von Mises, R., & Pollaczek-Geiringer, H. (1929). Praktische Verfahren der
Gleichungsauflösung. \emph{Zeitschrift für Angewandte Mathematik und
Mechanik}, 9(2), 152-164. \doi{10.1002/zamm.19290090206}
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 3)), directed = FALSE)
eigenvector_centrality(g)                          # igraph: largest score 1
eigenvector_centrality(g, convention = "networkx") # NetworkX: unit length
eigenvector_centrality(g, convention = "sna")      # sna: unit length
eigenvector_centrality(g, convention = "ucinet")   # UCINET: unit length

# In directed networks the tools follow different ties: igraph and
# NetworkX use incoming ties, sna outgoing ties, and UCINET ignores
# direction.
d <- netframe(data.frame(from = c(1, 2, 3, 1, 4), to = c(2, 3, 1, 4, 1)))
eigenvector_centrality(d)                          # igraph
eigenvector_centrality(d, convention = "networkx") # NetworkX
eigenvector_centrality(d, convention = "sna")      # sna
eigenvector_centrality(d, convention = "ucinet")   # UCINET
}
