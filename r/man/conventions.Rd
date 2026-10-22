\name{conventions}
\alias{conventions}
\title{Conventions of igraph, NetworkX, sna, and UCINET}
\description{
The same measure name means different numbers in different programs.
Functions whose results depend on such choices take a \code{convention}
argument, \code{"igraph"} (the default), \code{"networkx"}, \code{"sna"}, or
\code{"ucinet"}, which reproduces that tool. This page records the
differences; the engine records the same table in its source
(\code{conventions.rs}).
}
\details{
\eqn{n} is the number of nodes, \eqn{r_i} the number of other nodes node
\eqn{i} reaches, and \eqn{S_i} the sum of shortest-path lengths to them.
\tabular{lllll}{
  \tab igraph \tab NetworkX \tab sna \tab UCINET \cr
  Ties \tab every edge \tab every edge \tab adjacency matrix: parallel edges collapse, diagonal ignored \tab as sna \cr
  Parallel edges in shortest paths \tab distinct paths \tab one path (the shortest edge) \tab one path \tab one path \cr
  Normalized degree \tab \eqn{d/(n-1)} \tab \eqn{d/(n-1)} \tab none \tab \eqn{d/(n-1)} \cr
  Normalized betweenness \tab pairs not involving the node \tab same \tab none \tab same, as a percentage \cr
  Closeness \tab \eqn{1/S_i} or \eqn{r_i/S_i}; \code{NaN} if \eqn{r_i = 0} \tab \eqn{(r_i/S_i)(r_i/(n-1))}, paths to the node \tab \eqn{(n-1)/S_i}; 0 if any node is unreachable \tab \eqn{(n-1)/S_i}, unreachable lengths set to the longest plus 1 \cr
  Harmonic closeness \tab optional division by \eqn{n-1} \tab not normalized, paths to the node \tab divided by \eqn{n-1} \tab divided by \eqn{n-1} \cr
  Eigenvector scale \tab maximum 1 \tab Euclidean norm 1 \tab norm 1 \tab norm 1 \cr
  Eigenvector direction \tab in-ties (\code{directed = TRUE}) \tab in-ties \tab out-ties \tab symmetrized \cr
  Eigenvector, disconnected \tab computed \tab refused \tab computed \tab lambda-squared scores \cr
  Centralization maxima \tab see \code{\link{centralization}} \tab none \tab see \code{\link{centralization}} \tab see \code{\link{centralization}}
}
UCINET conventions are those of version 6.849, taken from its help pages and
from the port of its routines by its authors (Borgatti, Agneessens, and
Everett 2026), which states each formula from UCINET's own source and checks
it against UCINET's output. UCINET itself is not run in this package's
validation suite until golden files produced with it are supplied.
}
\seealso{
\code{\link{degree}}, \code{\link{betweenness}}, \code{\link{closeness}},
\code{\link{harmonic_centrality}}, \code{\link{eigenvector_centrality}},
\code{\link{centralization}}, \code{\link{graph_density}}
}
\references{
Analytic Technologies. UCINET 6 for Windows help: Closeness
(\url{http://www.analytictech.com/ucinet/help/oeky47.htm}), Betweenness
(\url{http://www.analytictech.com/ucinet/help/pvn0n9.htm}), Eigenvector
(\url{http://www.analytictech.com/ucinet/help/24d9zjj.htm}), and Multiple
Measures (\url{http://www.analytictech.com/ucinet/help/hs4214.htm}).

Borgatti, S. P., Agneessens, F., & Everett, M. G. (2026). xucinet:
UCINET-style social network analysis in R, version 2.0.0.9000.
\url{https://github.com/StephenBorgatti/xucinet}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods and
Applications}. Cambridge University Press. \doi{10.1017/CBO9780511815478}
}
\keyword{documentation}
