\name{net_compare_partitions}
\alias{net_compare_partitions}
\title{Compare two partitions}
\description{
Distance or similarity between two partitions of the same nodes, for
example the output of two community detection methods or a method and
known groups.
}
\usage{
net_compare_partitions(a, b, method = c("vi", "nmi", "split_join", "rand",
                                        "adjusted_rand", "ami"))
}
\arguments{
\item{a, b}{Two vectors of group labels of the same length, one per node
(any atomic type or factors; only equality of labels matters). Missing
values are an error.}
\item{method}{The measure; see \sQuote{Measures}.}
}
\details{
\strong{Measures.} With \eqn{n} nodes, \eqn{n_{ij}} the number of nodes in
group \eqn{i} of \code{a} and group \eqn{j} of \code{b}, and entropies and
mutual information \eqn{I} in natural units:
\describe{
  \item{\code{"vi"}}{Variation of information (Meila 2007),
    \eqn{H(a) + H(b) - 2 I(a, b)}, a distance (0 for identical
    partitions). igraph's \code{compare(method = "vi")}.}
  \item{\code{"nmi"}}{Normalized mutual information (Danon et al. 2005),
    \eqn{2 I(a, b) / (H(a) + H(b))}, a similarity in \eqn{[0, 1]}; 1 when
    both partitions have zero entropy. igraph's \code{method = "nmi"}.}
  \item{\code{"split_join"}}{The split-join distance of van Dongen (2000):
    the sum of the projection distances \eqn{n - \sum_i \max_j n_{ij}} of
    \code{a} from \code{b} and of \code{b} from \code{a}, in numbers of
    nodes. The sum of the two values of igraph's
    \code{split_join_distance()}, and igraph's
    \code{compare(method = "split.join")}.}
  \item{\code{"rand"}}{The Rand index (Rand 1971): the fraction of the
    \eqn{n (n - 1) / 2} node pairs that both partitions put together or
    both put apart. igraph's \code{method = "rand"}.}
  \item{\code{"adjusted_rand"}}{The Rand index corrected for chance under
    the hypergeometric model (Hubert and Arabie 1985); 1 when its
    denominator vanishes, which happens exactly when both partitions are a
    single group or both are all singletons. igraph's
    \code{method = "adjusted.rand"}.}
  \item{\code{"ami"}}{Adjusted mutual information (Vinh, Epps, and Bailey
    2010), \eqn{(I - E[I]) / (\bar{H} - E[I])} with \eqn{E[I]} the expected
    mutual information under the hypergeometric model and \eqn{\bar{H}} the
    arithmetic mean of the two entropies (the normalization that matches
    \code{"nmi"}). 1 in the same two cases as the adjusted Rand index.
    igraph has no counterpart.}
}

Without nodes every measure is \code{NaN}. The Rand indices are \code{NaN}
for a single node, where igraph raises an error.
}
\value{A single number.}
\references{
Meila, M. (2007). Comparing clusterings: an information based distance.
\emph{Journal of Multivariate Analysis}, 98(5), 873--895.
\doi{10.1016/j.jmva.2006.11.013}

Danon, L., Diaz-Guilera, A., Duch, J., & Arenas, A. (2005). Comparing
community structure identification. \emph{Journal of Statistical
Mechanics: Theory and Experiment}, 2005(09), P09008.
\doi{10.1088/1742-5468/2005/09/P09008}

van Dongen, S. (2000). Performance criteria for graph clustering and
Markov cluster experiments. Technical Report INS-R0012, Centrum voor
Wiskunde en Informatica, Amsterdam.

Rand, W. M. (1971). Objective criteria for the evaluation of clustering
methods. \emph{Journal of the American Statistical Association}, 66(336),
846--850. \doi{10.1080/01621459.1971.10482356}

Hubert, L., & Arabie, P. (1985). Comparing partitions. \emph{Journal of
Classification}, 2(1), 193--218. \doi{10.1007/BF01908075}

Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic measures
for clusterings comparison: Variants, properties, normalization and
correction for chance. \emph{Journal of Machine Learning Research}, 11,
2837--2854.
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_leiden}},
\code{\link{net_label_prop}}}
\examples{
g <- net_example("karate")
found <- net_louvain(g, seed = 1)
truth <- g$nodes$faction
net_compare_partitions(found, truth, "vi")              # igraph compare(method = "vi")
net_compare_partitions(found, truth, "nmi")             # igraph method = "nmi"
net_compare_partitions(found, truth, "split_join")      # igraph method = "split.join"
net_compare_partitions(found, truth, "rand")            # igraph method = "rand"
net_compare_partitions(found, truth, "adjusted_rand")   # igraph method = "adjusted.rand"
net_compare_partitions(found, truth, "ami")             # no igraph counterpart

# Labels are arbitrary: relabeled partitions are identical.
net_compare_partitions(found, 10 * found, "vi")
}
