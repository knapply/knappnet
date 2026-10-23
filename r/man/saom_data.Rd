\name{saom_data}
\alias{saom_data}
\alias{print.saom_data}
\title{Network panel data for stochastic actor-oriented models}
\usage{
saom_data(networks, covariates = NULL, dyadic = NULL)
}
\arguments{
\item{networks}{A list of at least two directed \code{\link{netframe}}s
observed at successive moments (for example the output of
\code{\link{snapshots}}), all on the same node identifiers. Nodes are
matched by identifier and ordered as in the first network. Each network is
read as a binary adjacency matrix: parallel edges collapse to one tie,
self-loops are dropped, and weights are ignored.}
\item{covariates}{Actor covariates: \code{NULL}, a character vector naming
columns of the node tables (a column whose values are the same in every
wave is a constant covariate; one that changes is a changing covariate
whose value during a period is its value at the start of that period), or
a named list of numeric vectors (one value per actor, in the order of the
first network's nodes) and matrices with one row per actor and one column
per wave. Values must not be missing.}
\item{dyadic}{Dyadic covariates: \code{NULL} or a named list of square
numeric matrices (rows and columns in the order of the first network's
nodes; the diagonal is ignored) or netframes on the same nodes, read as
adjacency matrices of their weights (or of ties when unweighted).}
}
\value{
An object of class \code{"saom_data"}: a list with \code{n} (the number of
actors), \code{waves} (a list of logical adjacency matrices), \code{nodes}
(the node table of the first network), \code{covariates}, and
\code{dyadic}.
}
\description{
Collects two or more observations of a directed network on the same
actors, with actor and dyadic covariates, into the data object that
\code{\link{saom_fit}} estimates stochastic actor-oriented models from
(Snijders 2001; Snijders, van de Bunt, and Steglich 2010).
}
\details{
Covariates are centered as RSiena centers them (Ripley et al. 2026, section
5.3.2): actor covariates by their mean over all actors and all waves
supplied, dyadic covariates by the mean of their off-diagonal cells, and the
similarity scores of \code{simX} by their mean over all ordered pairs of
distinct actors and waves. \code{sameX} is not centered.

Missing tie or covariate values, structural zeros and ones, and changes of
composition (actors joining or leaving) are not supported; see
\sQuote{Not implemented} in \code{\link{saom_fit}}.
}
\references{
Snijders, T. A. B. (2001). The statistical evaluation of social network
dynamics. \emph{Sociological Methodology}, 31(1), 361-395.
\doi{10.1111/0081-1750.00099}

Snijders, T. A. B., van de Bunt, G. G., and Steglich, C. E. G. (2010).
Introduction to stochastic actor-based models for network dynamics.
\emph{Social Networks}, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}

Ripley, R. M., Snijders, T. A. B., Boda, Z., and Voros, A. (2026).
\emph{Manual for RSiena}. University of Oxford, Department of Statistics and
Nuffield College; University of Groningen, Department of Sociology.
\url{https://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf}
}
\seealso{
\code{\link{saom_fit}}, \code{\link{saom_simulate}}, \code{\link{snapshots}}.
}
\examples{
w1 <- netframe(data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "a")),
               nodes = data.frame(name = c("a", "b", "c", "d"), smoke = c(1, 0, 1, 0)))
w2 <- netframe(data.frame(from = c("a", "b", "c", "d"), to = c("b", "c", "a", "a")),
               nodes = w1$nodes)
saom_data(list(w1, w2), covariates = "smoke")
}
