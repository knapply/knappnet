\name{bridging_centrality}
\alias{bridging_centrality}
\alias{bridging_coefficient}
\title{Bridging centrality}
\usage{
bridging_centrality(x, threads = 0L)

bridging_coefficient(x)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Hwang et al.'s (2008) bridging centrality finds nodes that lie between
densely connected regions: it is betweenness times the bridging
coefficient
\deqn{BC(v) = \frac{1 / d(v)}{\sum_{u \in N(v)} 1 / d(u)},}{BC(v) = (1 / d(v)) / sum_u 1 / d(u),}
where \eqn{d} is the number of distinct neighbours and the sum runs over
the neighbours of \eqn{v}. The coefficient is large for nodes of low
degree whose neighbours have high degree.
}
\details{
The bridging coefficient ignores direction, weights, self-loops, and
parallel edges, and is \code{NaN} for isolated nodes.
\code{bridging_centrality()} multiplies it by unnormalized shortest-path
betweenness computed on binary ties (\code{\link{betweenness}} with
\code{weighted = FALSE}), which respects direction in directed networks,
as Hwang et al. define it on the undirected networks they study. Neither
igraph, NetworkX, nor sna implements the measure.
}
\references{
Hwang, W., Kim, T., Ramanathan, M., & Zhang, A. (2008). Bridging
centrality: Graph mining from element level to group level. In
\emph{Proceedings of the 14th ACM SIGKDD International Conference on
Knowledge Discovery and Data Mining} (pp. 336-344). ACM.
\doi{10.1145/1401890.1401934}
}
\seealso{\code{\link{betweenness}}, \code{\link{constraint}}.}
\examples{
# Two triangles joined by a path through node 4: node 4 bridges them.
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 5, 5, 6),
                         to   = c(2, 3, 3, 4, 5, 6, 7, 7)), directed = FALSE)
bridging_coefficient(g)
bridging_centrality(g)
}
