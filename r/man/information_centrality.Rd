\name{information_centrality}
\alias{information_centrality}
\title{Information centrality}
\description{
Stephenson and Zelen's (1989) centrality from the information carried by
every path between a node and the others, not only the shortest. With
\eqn{A} the symmetric adjacency matrix, \eqn{B = L + J} (the Laplacian
plus a matrix of ones) and \eqn{C = B^{-1}},
\deqn{I_i = \frac{1}{c_{ii} + (T - 2R)/n},}{I_i = 1 / (c_ii + (T - 2R)/n),}
where \eqn{T} is the trace of \eqn{C} and \eqn{R} any of its row sums.
}
\usage{
information_centrality(x, weighted = NULL, convention = c("sna", "networkx"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{weighted}{Whether to use edge weights as tie strengths (conductances:
a stronger tie carries more information), whatever the network's
\code{weight_type}. \code{NULL} uses them when the network has a weight
column. They must be positive.}
\item{convention}{\code{"sna"} returns \eqn{I_i} as below, the value of
\code{sna::infocent()} and UCINET; \code{"networkx"} returns \eqn{I_i / n},
the value of NetworkX's \code{information_centrality()}, which computes
current-flow closeness (Brandes and Fleischer 2005 show the two agree up
to that factor).}
}
\details{
The measure is defined for undirected ties: a directed network is
symmetrized with \eqn{a_{ij} = \max(a_{ij}, a_{ji})}{a_ij = max(a_ij,
a_ji)}, sna's default for \code{infocent()}. Self-loops are ignored.
Without weights parallel edges count once (the adjacency matrix has a 1);
with weights \eqn{a_{ij}}{a_ij} is the total weight of the edges between
\eqn{i} and \eqn{j}.

Isolated nodes score 0 and are left out of the matrix, as in sna, so
\eqn{n} is the number of nodes with ties. The formula needs those nodes to
be connected (otherwise \eqn{B} is singular); when they are not, every
node with ties gets \code{NaN}. NetworkX refuses such networks altogether,
so under the \code{"networkx"} convention they are an error, and \eqn{n}
is the number of nodes in the network, with no isolates allowed.

The inverse is computed by Cholesky factorization in the compiled engine.
The validation suite compares the result with \code{sna::infocent()} and
NetworkX's \code{information_centrality()}.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Stephenson, K., & Zelen, M. (1989). Rethinking centrality: Methods and
examples. \emph{Social Networks}, 11(1), 1-37.
\doi{10.1016/0378-8733(89)90016-6}

Brandes, U., & Fleischer, D. (2005). Centrality measures based on current
flow. In \emph{STACS 2005} (Lecture Notes in Computer Science 3404, pp.
533-544). Springer. \doi{10.1007/978-3-540-31856-9_44}

Butts, C. T. (2008). Social network analysis with sna. \emph{Journal of
Statistical Software}, 24(6), 1-51. \doi{10.18637/jss.v024.i06}
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 3)),
              directed = FALSE)
information_centrality(g)                           # sna::infocent(), UCINET
information_centrality(g, convention = "networkx")  # nx.information_centrality()
}
