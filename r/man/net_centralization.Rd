\name{net_centralization}
\alias{net_centralization}
\title{Freeman centralization}
\description{
How far the distribution of a centrality is from the most unequal one
possible: the sum of differences between the largest score and every score,
divided by the largest value that sum can take in a network of the same
size, which a star attains (Freeman 1978).
}
\usage{
net_centralization(x, measure = c("degree", "betweenness", "closeness", "eigenvector"),
               mode = NULL,
               loops = TRUE, normalized = TRUE,
               convention = c("igraph", "networkx", "sna", "ucinet"),
               threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{measure}{\code{"degree"}, \code{"betweenness"}, \code{"closeness"},
or \code{"eigenvector"}.}
\item{mode}{For directed networks: which degree, or which direction of paths
for closeness. Ignored for betweenness and undirected networks, and always
\code{"all"} for eigenvector centrality. \code{NULL} uses igraph's defaults
(\code{"all"} for degree and eigenvector centrality, \code{"out"} for
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
\code{\link{net_frame}}).}
}
\details{
Degree uses edge counts and betweenness is unnormalized shortest-path
betweenness; weights are not used. Closeness centralization uses normalized
closeness as each tool defines it (see \code{\link{net_closeness}}); under the
igraph convention it is \code{NaN} if some node reaches no other, and under
the UCINET convention unless every node reaches every other, as in those
tools.

Eigenvector centralization follows igraph's \code{centr_eigen()}: scores
scaled so the largest is 1 (see \code{\link{net_eigen}}),
direction ignored, and the theoretical maximum \eqn{n - 2}, attained by a
network of one edge and \eqn{n - 2} isolated nodes (scores 1, 1, and 0),
which igraph's documentation gives as the most centralized network for this
measure; a star then scores less than 1 (0.75 for ten nodes). It is
\code{NaN} for fewer than three nodes. When more than one component has
ties the eigenvector depends on more than the network (see
\code{\link{net_eigen}}; a warning is given), and so does its
centralization. sna and UCINET compute eigenvector scores with unit norm
and use maxima of their own that could not be verified from their
documentation, so only the igraph convention is provided for this measure.

The theoretical maxima are:
\tabular{llll}{
  Measure \tab igraph \tab sna \tab UCINET \cr
  degree, undirected \tab \eqn{(n-1)(n-2)}, or \eqn{(n-1)n} with loops \tab \eqn{(n-1)(n-2)} \tab \eqn{(n-1)(n-2)} \cr
  in- or out-degree \tab \eqn{(n-1)^2}, or \eqn{(n-1)n} with loops \tab \eqn{(n-1)^2} \tab \eqn{(n-1)^2} \cr
  total degree \tab \eqn{2(n-1)(n-2)}, or \eqn{2(n-1)^2} with loops \tab \eqn{2(n-1)(n-2)} \tab (not reported) \cr
  betweenness \tab \eqn{(n-1)^2(n-2)}, halved if undirected \tab same \tab same \cr
  closeness, undirected or \code{"all"} \tab \eqn{(n-1)(n-2)/(2n-3)} \tab same \tab same \cr
  closeness, in or out \tab \eqn{(n-1)^2/n} \tab same \tab \eqn{(n-1)(n-2)/(2n-3)} \cr
  eigenvector (largest score 1) \tab \eqn{n-2} \tab (not provided) \tab (not provided)
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
\value{A number.}
\references{
Freeman, L. C. (1978). Centrality in social networks: Conceptual
clarification. \emph{Social Networks}, 1(3), 215-239.
\doi{10.1016/0378-8733(78)90021-7}
}
\seealso{\code{\link{net_centralize}} for the index of any vector of scores.}
\examples{
star <- net_frame(data.frame(source = 1, target = 2:6), directed = FALSE)
net_centralization(star)
net_centralization(star, "betweenness")

# The theoretical maxima differ between tools. NetworkX has no
# centralization index, so its convention is an error.
net_centralization(star)                          # igraph
net_centralization(star, convention = "sna")      # sna
net_centralization(star, convention = "ucinet")   # UCINET
net_centralization(star, "closeness")                        # igraph
net_centralization(star, "closeness", convention = "sna")    # sna
net_centralization(star, "closeness", convention = "ucinet") # UCINET

# Eigenvector centralization (igraph's centr_eigen() only): a star of
# ten nodes scores 0.75, and one edge among ten nodes scores 1.
net_centralization(net_frame(data.frame(source = 1, target = 2:10), directed = FALSE), "eigenvector")
net_centralization(net_frame(data.frame(source = 1, target = 2), nodes = data.frame(id = 1:10),
                        directed = FALSE), "eigenvector")

# Centralization of any scores, given the theoretical maximum:
net_centralize(net_degree(star), theoretical_max = 5 * 4)
}
