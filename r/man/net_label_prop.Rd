\name{net_label_prop}
\alias{net_label_prop}
\title{Label propagation communities}
\description{
Communities by label propagation (Raghavan, Albert, and Kumara 2007):
every node starts with its own label (or a label you give it) and nodes
repeatedly adopt the label that most of their neighbors carry, until every
node carries a most frequent label. Four variants of the update and
stopping rules are offered, each reproducing a routine of a reference
tool.
}
\usage{
net_label_prop(x, variant = c("dominance", "retention", "fast", "semi_synchronous"),
               mode = c("all", "out", "in"), weighted = NULL, initial = NULL,
               fixed = NULL, seed = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}; the \code{"semi_synchronous"} variant
needs an undirected one.}
\item{variant}{Update and stopping rule; see \sQuote{Variants}.}
\item{mode}{For directed networks, the direction labels travel:
\code{"out"} propagates labels along the edges (a node reads its
in-neighbors), \code{"in"} against them (it reads its out-neighbors), and
\code{"all"} ignores directions, as igraph advises for directed networks.
Undirected networks ignore it.}
\item{weighted}{Whether to use edge weights, as given. \code{NULL} (the
default) uses them when the network has a weight column. They must not be
negative.}
\item{initial}{Starting labels: a vector with one value per node (any
atomic type) or the name of a column of \code{x$nodes}; \code{NA} marks a
node without a label. \code{NULL} gives every node its own label.}
\item{fixed}{A logical vector with one value per node (or the name of such
a column of \code{x$nodes}) marking nodes whose label never changes; they
must be labeled in \code{initial}.}
\item{seed}{Seed for the random node orders and tie breaking, a whole
number between 0 and \eqn{2^{53}}; \code{NULL} draws one from R's random
number generator, so \code{\link{set.seed}()} makes results reproducible.
The \code{"semi_synchronous"} variant uses no random numbers.}
}
\details{
\strong{Variants.} A node carries a \emph{dominant} label when no label has
a larger total weight among the neighbors it reads.
\describe{
  \item{\code{"dominance"}}{igraph's default (\code{cluster_label_prop()}).
    Each sweep visits the nodes in a new random order and every node
    draws uniformly among its dominant labels; the method stops after a
    sweep at whose end every node carries a dominant label, the stopping
    rule of Raghavan et al.}
  \item{\code{"retention"}}{NetworkX's \code{asyn_lpa_communities()}. A
    node keeps its label when it is dominant and otherwise draws among the
    dominant labels; the method stops after a sweep that changes no
    label.}
  \item{\code{"fast"}}{Fast label propagation (Traag and Subelj 2023), as
    igraph's \code{cluster_label_prop(lpa_variant = "fast")} and
    NetworkX's \code{fast_label_propagation_communities()}. Nodes wait in
    a queue and a node whose label changes puts the neighbors that read it
    back in the queue; the method stops when the queue is empty.}
  \item{\code{"semi_synchronous"}}{Cordasco and Gargano (2010), as
    NetworkX's \code{label_propagation_communities()}. Nodes are colored
    greedily and each round updates the color classes in turn; a node takes
    its unique dominant label and among several keeps its own if dominant
    and otherwise takes the largest. It is deterministic, ignores
    \code{seed}, and converges on every undirected network.}
}

Self-loops are ignored (a node never reads its own label; NetworkX counts
a loop as one neighbor), and parallel edges add their weights (NetworkX
counts them once). Weighted ties count with their weight. Randomized
variants stop after at most 100,000 sweeps even if the stopping rule is not
met, which can happen when ties make labels oscillate. Nodes still
unlabeled at the end (no labeled node reaches them) each form a group of
their own; isolated nodes keep their own group. Because the methods are
randomized, partitions differ from those of the reference tools even for
the same variant; the validation suite compares the mean modularity over
many seeds and recomputes the stopping rule from its definition.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance.}
\references{
Raghavan, U. N., Albert, R., & Kumara, S. (2007). Near linear time
algorithm to detect community structures in large-scale networks.
\emph{Physical Review E}, 76(3), 036106. \doi{10.1103/PhysRevE.76.036106}

Cordasco, G., & Gargano, L. (2010). Community detection via
semi-synchronous label propagation algorithms. In \emph{2010 IEEE
International Workshop on Business Applications of Social Network Analysis
(BASNA)}, 1--8. \doi{10.1109/BASNA.2010.5730298}

Traag, V. A., & Subelj, L. (2023). Large network community detection by
fast label propagation. \emph{Scientific Reports}, 13, 2701.
\doi{10.1038/s41598-023-29610-z}
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
net_label_prop(g, "dominance", seed = 1)       # igraph cluster_label_prop()
net_label_prop(g, "retention", seed = 1)       # NetworkX asyn_lpa_communities()
net_label_prop(g, "fast", seed = 1)            # igraph lpa_variant = "fast";
                                               # NetworkX fast_label_propagation_communities()
net_label_prop(g, "semi_synchronous")          # NetworkX label_propagation_communities()

# Directed networks: labels along the edges, against them, or ignoring them.
d <- net_example("krackhardt", relation = "advice")
net_label_prop(d, mode = "out", seed = 1)      # igraph mode = "out"
net_label_prop(d, mode = "in", seed = 1)       # igraph mode = "in"
net_label_prop(d, mode = "all", seed = 1)      # igraph mode = "all"

# Semi-supervised: two members keep the label of their faction.
start <- rep(NA, 34)
start[c(1, 34)] <- c("Mr. Hi", "Officer")
groups <- net_label_prop(g, initial = start, fixed = !is.na(start), seed = 1)
net_plot_partition(g, groups, seed = 1)
}
