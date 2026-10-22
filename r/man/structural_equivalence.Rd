\name{structural_equivalence}
\alias{structural_equivalence}
\title{Structural equivalence}
\usage{
structural_equivalence(x, method = c("euclidean", "correlation", "hamming"),
                       weighted = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{method}{\code{"euclidean"}, \code{"correlation"}, or
\code{"hamming"}.}
\item{weighted}{Whether ties are total edge weights rather than 0 or 1.
\code{NULL} (the default) uses weights when the network has a weight
column.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A symmetric numeric matrix with node identifiers as dimnames.}
\description{
Compares every pair of nodes by their tie profiles: the ties each node sends
(its row of the adjacency matrix) and receives (its column).
Structurally equivalent nodes (Lorrain and White 1971) have identical
profiles.
}
\details{
Nodes i and j are compared over every other node k, so the diagonal and the
pair's own ties \eqn{a_{ij}}{a_ij} and \eqn{a_{ji}}{a_ji} are left out
(Wasserman and Faust 1994, sec. 9.2):
\describe{
  \item{\code{"euclidean"}}{\eqn{\sqrt{\sum_{k \ne i,j} (a_{ik} - a_{jk})^2 +
    (a_{ki} - a_{kj})^2}}{sqrt(sum_{k != i,j} (a_ik - a_jk)^2 + (a_ki -
    a_kj)^2)} (Burt 1976).}
  \item{\code{"correlation"}}{The Pearson correlation of the two profiles
    over the same \eqn{2(n-2)} entries (Wasserman and Faust 1994). It is 0
    when undefined (a profile without variance over those entries). The
    diagonal is 1, or 0 for a node whose profile is constant.}
  \item{\code{"hamming"}}{The sum of absolute differences over the same
    entries, which in a binary network counts the ties that differ.}
}
These match \code{sna::sedist()} with its defaults (\code{diag = FALSE}).

In a binary network parallel edges count once. With weights,
\eqn{a_{ij}}{a_ij} is the total weight of the edges from i to j, used as
given whatever the network's \code{weight_type}. An undirected edge is a tie
in both directions, so each tie appears twice in a profile, as with a
symmetric adjacency matrix in sna. Self-loops are ignored.
}
\references{
Lorrain, F., & White, H. C. (1971). Structural equivalence of individuals in
social networks. \emph{Journal of Mathematical Sociology}, 1(1), 49-80.
\doi{10.1080/0022250X.1971.9989788}

Burt, R. S. (1976). Positions in networks. \emph{Social Forces}, 55(1),
93-122. \doi{10.2307/2577097}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4)))
structural_equivalence(g)
structural_equivalence(g, "correlation")
}
