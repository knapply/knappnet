\name{lambda_sets}
\alias{lambda_sets}
\alias{print.lambda_sets}
\alias{plot.lambda_sets}
\title{Lambda sets and edge connectivity}
\description{
The edge connectivity \eqn{\lambda(u, v)} of two nodes is the smallest
number of edges whose removal disconnects them, which equals the number of
edge-disjoint paths between them (Menger's theorem; Ford and Fulkerson
1956). A lambda set is a set of nodes any two of which are more
edge-connected to each other than either is to any node outside (Borgatti,
Everett, and Shirey 1990). Because
\eqn{\lambda(a, c) \ge \min(\lambda(a, b), \lambda(b, c))}, the relation
"\eqn{\lambda(a, b) \ge k}" is an equivalence for every \eqn{k \ge 1}, and
its classes with at least two members are exactly the lambda sets;
increasing \eqn{k} gives the nested hierarchy that UCINET's Network >
Subgroups > Lambda Sets routine reports.
}
\usage{
lambda_sets(x)

\method{print}{lambda_sets}(x, ...)

\method{plot}{lambda_sets}(x, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}} (for \code{print} and \code{plot}, a
\code{lambda_sets} object).}
\item{...}{For \code{plot}, arguments passed to
\code{\link[stats]{plot.hclust}}, such as \code{main}; ignored by
\code{print}.}
}
\details{
All pairwise connectivities come from \eqn{n - 1} maximum-flow
computations with Gusfield's (1990) method (a Gomory-Hu tree; Gomory and
Hu 1961), each flow found by breadth-first augmenting paths (Edmonds and
Karp 1972). The undirected simple graph underlying the network is used:
directions, parallel edges, self-loops, and weights are ignored. Nodes in
different components have connectivity 0. The connectivity matrix matches
igraph's \code{edge_connectivity(g, source, target)} and NetworkX's
\code{edge_connectivity(G, s, t)} on the simple undirected graph, and the
validation suite checks it against edge-disjoint paths counted directly.
}
\value{
A list of class \code{lambda_sets} with \code{connectivity}, the
\eqn{n \times n} matrix of pairwise edge connectivities (zero diagonal) with
node identifiers as dimnames, and \code{sets}, a list of the lambda sets
from the lowest to the highest connectivity: each element is a list with
\code{lambda} (the connectivity \eqn{k}) and \code{members} (node
identifiers). A set that is the same at several levels is listed once, with
\code{lambda} the highest connectivity its members share. \code{plot} draws
the nested sets as a dendrogram (single linkage on the connectivities, whose
heights are labelled with the connectivity) and returns \code{x}
invisibly.
}
\references{
Borgatti, S. P., Everett, M. G., & Shirey, P. R. (1990). LS sets, lambda
sets and other cohesive subsets. \emph{Social Networks}, 12(4), 337-357.
\doi{10.1016/0378-8733(90)90014-Z}

Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
\emph{Canadian Journal of Mathematics}, 8, 399-404.
\doi{10.4153/CJM-1956-045-5}

Gomory, R. E., & Hu, T. C. (1961). Multi-terminal network flows.
\emph{Journal of the Society for Industrial and Applied Mathematics}, 9(4),
551-570. \doi{10.1137/0109047}

Gusfield, D. (1990). Very simple methods for all pairs network flow
analysis. \emph{SIAM Journal on Computing}, 19(1), 143-155.
\doi{10.1137/0219009}

Edmonds, J., & Karp, R. M. (1972). Theoretical improvements in algorithmic
efficiency for network flow problems. \emph{Journal of the ACM}, 19(2),
248-264. \doi{10.1145/321694.321699}
}
\seealso{\code{\link{n_cliques}}, \code{\link{bicomponents}},
\code{\link{bridges}}}
\examples{
# Two triangles joined by a single bridge.
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5), to = c(2, 3, 3, 4, 5, 6, 6)),
              directed = FALSE)
ls <- lambda_sets(g)
ls$connectivity
# Every node is 1-connected to every other; the triangles are 2-connected.
ls
plot(ls)
}
