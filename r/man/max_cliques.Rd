\name{max_cliques}
\alias{max_cliques}
\alias{clique_number}
\title{Maximal cliques and the clique number}
\description{
A clique is a set of nodes that are all pairwise adjacent; it is maximal
when no other node is adjacent to every member. \code{max_cliques()} lists
the maximal cliques with the Bron-Kerbosch algorithm (Bron and Kerbosch
1973) using the pivot rule of Tomita, Tanaka, and Takahashi (2006) over a
degeneracy ordering (Eppstein, Loffler, and Strash 2010).
\code{clique_number()} is the size of the largest clique.
}
\usage{
max_cliques(x, min_size = 1L, threads = 0L)

clique_number(x, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{min_size}{Smallest clique size to return.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
Cliques are found in the undirected simple graph underlying the network:
edge directions are ignored (an arc in either direction makes two nodes
adjacent), parallel edges count once, and self-loops are ignored. Weights
are ignored. An isolated node is a maximal clique of size 1, and a network
without nodes has clique number 0.

\code{min_size} drops maximal cliques with fewer members; it does not list
smaller subsets of larger cliques.

These conventions match \code{igraph::max_cliques()} and
\code{igraph::clique_num()}, and NetworkX's \code{find_cliques()} applied
to the undirected simple graph.
}
\value{
\code{max_cliques()}: a list of vectors of node identifiers (values of the
identifier column of \code{x$nodes}), one per clique, each in node order;
the cliques are in lexicographic order of node positions.
\code{clique_number()}: a single integer, the size of the largest clique.
}
\references{
Bron, C., & Kerbosch, J. (1973). Algorithm 457: Finding all cliques of an
undirected graph. \emph{Communications of the ACM}, 16(9), 575-577.
\doi{10.1145/362342.362367}

Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time
complexity for generating all maximal cliques and computational
experiments. \emph{Theoretical Computer Science}, 363(1), 28-42.
\doi{10.1016/j.tcs.2006.06.015}

Eppstein, D., Loffler, M., & Strash, D. (2010). Listing all maximal cliques
in sparse graphs in near-optimal time. In \emph{Algorithms and Computation
(ISAAC 2010)}, Lecture Notes in Computer Science 6506, 403-414.
\doi{10.1007/978-3-642-17517-6_36}

Matula, D. W., & Beck, L. L. (1983). Smallest-last ordering and clustering
and graph coloring algorithms. \emph{Journal of the ACM}, 30(3), 417-427.
\doi{10.1145/2402.322385}
}
\examples{
g <- netframe(data.frame(from = c("a", "b", "a", "c"), to = c("b", "c", "c", "d")),
              directed = FALSE)
max_cliques(g)
max_cliques(g, min_size = 3)
clique_number(g)
}
