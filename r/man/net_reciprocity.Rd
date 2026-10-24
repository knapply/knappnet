\name{net_reciprocity}
\alias{net_reciprocity}
\title{Reciprocity}
\description{
How far ties in a directed network are returned. Reciprocity is a property
of dyads, so the network is read as an adjacency matrix: parallel arcs count
once and self-loops are ignored.
}
\usage{
net_reciprocity(x, mode = c("default", "ratio", "dyadic", "correlation"))
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{mode}{Which index; see Details.}
}
\details{
With \eqn{M} mutual, \eqn{A} asymmetric, and \eqn{N} null dyads among the
\eqn{n(n-1)/2} pairs of nodes:
\describe{
  \item{\code{"default"}}{The share of arcs that are reciprocated,
    \eqn{2M / (2M + A)} (Newman, Forrest, and Balthrop 2002). This is
    igraph's \code{net_reciprocity()}, sna's
    \code{grecip(measure = "edgewise")}, and NetworkX's
    \code{overall_reciprocity()}.}
  \item{\code{"ratio"}}{Mutual dyads over connected dyads,
    \eqn{M / (M + A)}: igraph's \code{mode = "ratio"} and sna's
    \code{"dyadic.nonnull"}.}
  \item{\code{"dyadic"}}{Symmetric dyads (mutual or null) over all dyads,
    \eqn{(M + N) / (M + A + N)}: sna's \code{"dyadic"} (Wasserman and Faust
    1994).}
  \item{\code{"correlation"}}{Garlaschelli and Loffredo's (2004)
    \eqn{\rho = (r - p) / (1 - p)}{rho = (r - p) / (1 - p)}, where \eqn{r}
    is the \code{"default"} index and \eqn{p} the density: the correlation
    between \eqn{x_{ij}}{x_ij} and \eqn{x_{ji}}{x_ji} over ordered pairs of
    distinct nodes, which is sna's \code{"correlation"}. It is 0 when ties
    are returned as often as chance would have it, positive when more often,
    and negative when less.}
}
A ratio whose denominator is 0 (a network without arcs, or a complete one
for \code{"correlation"}) is \code{NaN}; sna reports 1 for the correlation
in those cases. An undirected network is read as a
symmetric directed one, so every tie is reciprocated. Weights are ignored.

igraph's \code{net_reciprocity()} also ignores self-loops by default
(\code{ignore.loops = TRUE}); the validation suite compares it on the
simplified network, where parallel arcs are merged as here.
}
\value{A number.}
\references{
Newman, M. E. J., Forrest, S., & Balthrop, J. (2002). Email networks and the
spread of computer viruses. \emph{Physical Review E}, 66(3), 035101.
\doi{10.1103/PhysRevE.66.035101}

Garlaschelli, D., & Loffredo, M. I. (2004). Patterns of link reciprocity in
directed networks. \emph{Physical Review Letters}, 93(26), 268701.
\doi{10.1103/PhysRevLett.93.268701}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods and
Applications}. Cambridge University Press. \doi{10.1017/CBO9780511815478}
}
\examples{
g <- net_frame(data.frame(source = c("a", "b", "b", "c"), target = c("b", "a", "c", "d")))
net_reciprocity(g)
net_reciprocity(g, "ratio")
net_reciprocity(g, "correlation")
}
