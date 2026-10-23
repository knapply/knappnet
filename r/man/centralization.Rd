\name{centralization}
\alias{centralization}
\title{Freeman centralization}
\usage{
centralization(x, measure = c("degree", "betweenness", "closeness"), mode = NULL,
               loops = TRUE, normalized = TRUE,
               convention = c("igraph", "networkx", "sna", "ucinet"),
               threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{measure}{\code{"degree"}, \code{"betweenness"}, or
\code{"closeness"}.}
\item{mode}{For directed networks: which degree, or which direction of paths
for closeness. Ignored for betweenness and undirected networks. \code{NULL}
uses igraph's defaults (\code{"all"} for degree, \code{"out"} for
closeness), or out-degree under the UCINET convention, which reports only
in- and out-degree centralization.}
\item{loops}{Whether self-loops count toward degree (igraph convention
only).}
\item{normalized}{Whether to divide by the theoretical maximum.}
\item{convention}{Whose theoretical maxima and closeness definition to use:
\code{"igraph"}, \code{"sna"}, or \code{"ucinet"}. NetworkX has no
centralization index, so \code{"networkx"} is an error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A number.}
\description{
How far the distribution of a centrality is from the most unequal one
possible: the sum of differences between the largest score and every score,
divided by the largest value that sum can take in a network of the same
size, which a star attains (Freeman 1978).
}
\details{
Degree uses edge counts and betweenness is unnormalized shortest-path
betweenness; weights are not used. Closeness centralization uses normalized
closeness as each tool defines it (see \code{\link{closeness}}); under the
igraph convention it is \code{NaN} if some node reaches no other, and under
the UCINET convention unless every node reaches every other, as in those
tools.

The theoretical maxima are:
\tabular{llll}{
  Measure \tab igraph \tab sna \tab UCINET \cr
  degree, undirected \tab \eqn{(n-1)(n-2)}, or \eqn{(n-1)n} with loops \tab \eqn{(n-1)(n-2)} \tab \eqn{(n-1)(n-2)} \cr
  in- or out-degree \tab \eqn{(n-1)^2}, or \eqn{(n-1)n} with loops \tab \eqn{(n-1)^2} \tab \eqn{(n-1)^2} \cr
  total degree \tab \eqn{2(n-1)(n-2)}, or \eqn{2(n-1)^2} with loops \tab \eqn{2(n-1)(n-2)} \tab (not reported) \cr
  betweenness \tab \eqn{(n-1)^2(n-2)}, halved if undirected \tab same \tab same \cr
  closeness, undirected or \code{"all"} \tab \eqn{(n-1)(n-2)/(2n-3)} \tab same \tab same \cr
  closeness, in or out \tab \eqn{(n-1)^2/n} \tab same \tab \eqn{(n-1)(n-2)/(2n-3)}
}
With the sna and UCINET conventions the network is read as an adjacency
matrix: parallel edges collapse and self-loops are ignored. The igraph
convention matches \code{centr_degree()}, \code{centr_betw()}, and
\code{centr_clo()}; the sna convention matches \code{sna::centralization()}
with \code{degree}, \code{betweenness}, and \code{closeness}. UCINET reports
centralization as a percentage in some routines; this function always
returns a proportion. UCINET versions before 6.849 used \eqn{(n-1)(n-2)} for
directed degree centralization.
}
\references{
Freeman, L. C. (1978). Centrality in social networks: Conceptual
clarification. \emph{Social Networks}, 1(3), 215-239.
\doi{10.1016/0378-8733(78)90021-7}
}
\examples{
star <- netframe(data.frame(from = 1, to = 2:6), directed = FALSE)
centralization(star)
centralization(star, "betweenness")

# The theoretical maxima differ between tools. NetworkX has no
# centralization index, so its convention is an error.
centralization(star)                          # igraph
centralization(star, convention = "sna")      # sna
centralization(star, convention = "ucinet")   # UCINET
centralization(star, "closeness")                        # igraph
centralization(star, "closeness", convention = "sna")    # sna
centralization(star, "closeness", convention = "ucinet") # UCINET
}
