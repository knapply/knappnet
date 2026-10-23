\name{cluster_hierarchical}
\alias{cluster_hierarchical}
\alias{cut_clusters}
\title{Hierarchical clustering of nodes}
\description{
Agglomerative hierarchical clustering (Johnson 1967) of the nodes from any
distance or similarity matrix, such as \code{\link{distances}} or
\code{\link{structural_equivalence}}, as in UCINET's Tools > Cluster
Analysis > Hierarchical routine. \code{cluster_hierarchical()} wraps
\code{\link[stats]{hclust}}; \code{cut_clusters()} cuts the tree into groups
numbered in order of first appearance.
}
\usage{
cluster_hierarchical(d, method = c("average", "single", "complete", "ward"),
                     similarity = FALSE)

cut_clusters(tree, k = NULL, height = NULL)
}
\arguments{
\item{d}{A square numeric matrix of distances (or similarities) between
nodes.}
\item{method}{Linkage: \code{"average"}, \code{"single"},
\code{"complete"}, or \code{"ward"}.}
\item{similarity}{Whether \code{d} holds similarities (larger means closer)
rather than distances.}
\item{tree}{The result of \code{cluster_hierarchical()} (or any
\code{hclust}).}
\item{k}{Number of groups to cut into, or \code{NULL} to cut at
\code{height}.}
\item{height}{Height at which to cut, when \code{k} is \code{NULL}.}
}
\details{
A similarity matrix (\code{similarity = TRUE}, for example correlations
between rows of the adjacency matrix) is turned into distances by
subtracting it from its largest value. Unreachable pairs (\code{Inf} in a
distance matrix) are given one more than the largest finite distance, so
they join last; missing values are an error. The matrix must be square and
symmetric; its dimnames label the leaves.

The \code{"average"}, \code{"single"}, and \code{"complete"} methods are
UPGMA, nearest neighbor, and farthest neighbor (Sokal and Michener 1958;
Johnson 1967). \code{"ward"} is Ward's (1963) method on squared distances,
\code{hclust()}'s \code{"ward.D2"}. Ties between equal distances are broken
as \code{hclust()} breaks them.
}
\value{
\code{cluster_hierarchical()}: an object of class \code{hclust}, with the
row names of \code{d} as labels. \code{cut_clusters()}: an integer vector of
group labels in the order of the rows of \code{d}, numbered 1, 2, ... in
order of first appearance.
}
\references{
Johnson, S. C. (1967). Hierarchical clustering schemes.
\emph{Psychometrika}, 32(3), 241-254. \doi{10.1007/BF02289588}

Sokal, R. R., & Michener, C. D. (1958). A statistical method for evaluating
systematic relationships. \emph{University of Kansas Science Bulletin}, 38,
1409-1438.

Ward, J. H. (1963). Hierarchical grouping to optimize an objective
function. \emph{Journal of the American Statistical Association}, 58(301),
236-244. \doi{10.1080/01621459.1963.10500845}
}
\seealso{\code{\link{distances}}, \code{\link{structural_equivalence}},
\code{\link{concor}}}
\examples{
g <- example_network("florentine")
# Cluster families by geodesic distance.
tree <- cluster_hierarchical(distances(g), method = "average")
plot(tree)
cut_clusters(tree, k = 3)

# Or by structural equivalence (a distance), or correlation (a similarity).
tree2 <- cluster_hierarchical(structural_equivalence(g), method = "complete")
cut_clusters(tree2, height = 2)
}
