\name{transitivity}
\alias{transitivity}
\alias{clustering}
\alias{average_clustering}
\title{Transitivity and clustering coefficients}
\description{
\code{transitivity()} is the share of connected triples that are closed:
three times the number of triangles divided by the number of paths of length
two (Newman 2003). \code{clustering()} gives each node's local clustering
coefficient, the share of pairs of its neighbors that are themselves
adjacent (Watts and Strogatz 1998), and \code{average_clustering()} their
mean.
}
\usage{
transitivity(x, convention = c("igraph", "networkx"), threads = 0L)

clustering(x, weighted = FALSE, method = NULL,
           convention = c("igraph", "networkx"), threads = 0L)

average_clustering(x, weighted = FALSE, method = NULL,
                   convention = c("igraph", "networkx"), threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{convention}{How a coefficient that is 0 / 0 is reported, and which
weighted coefficient is the default: \code{"igraph"} (the default) or
\code{"networkx"}; see Details and \code{\link{conventions}}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
\item{weighted}{Whether to compute a weighted local coefficient (see Details).
The binary coefficient is the default, whatever the network's weights, as in
igraph and NetworkX.}
\item{method}{The weighted coefficient: \code{"barrat"} (Barrat et al. 2004,
igraph) or \code{"onnela"} (Onnela et al. 2005, NetworkX). \code{NULL} uses
the convention's choice.}
}
\details{
All three use the simple undirected reading of the network: directions are
ignored, parallel edges count once, and self-loops are ignored, as in
igraph's \code{transitivity()}. A node's degree \eqn{k} is its number of
distinct neighbors other than itself. Triangles are listed once each after
ordering nodes by degree (Chiba and Nishizeki 1985), in parallel.

The coefficient of a node with fewer than two neighbors, and the
transitivity of a network without connected triples, are 0 / 0. Under the
igraph convention they are \code{NaN}, and \code{average_clustering()} leaves
those nodes out, as igraph's
\code{transitivity(type = "average", isolates = "NaN")} does (\code{NaN} if
no node has two neighbors). Under the NetworkX convention they are 0, and the
average includes them as 0, as NetworkX's \code{clustering()},
\code{transitivity()}, and \code{average_clustering()} do (and igraph with
\code{isolates = "zero"}). The two averages can differ greatly in networks
with many leaves (Kaiser 2008).

Weighted coefficients use, for each pair of adjacent nodes, the total weight
of the edges joining them in either direction, read as tie strengths
whatever the network's \code{weight_type}; weights must be positive. With
\eqn{s_i} the sum of node \eqn{i}'s weights, Barrat's coefficient is
\deqn{\frac{1}{s_i (k_i - 1)} \sum_{j,h} \frac{w_{ij} + w_{ih}}{2} a_{ij}
a_{ih} a_{jh}}{1 / (s_i (k_i - 1)) sum_(j,h) (w_ij + w_ih) / 2 a_ij a_ih a_jh}
(igraph's \code{type = "barrat"}), and Onnela's is
\deqn{\frac{1}{k_i (k_i - 1)} \sum_{j,h} (\hat w_{ij} \hat w_{ih} \hat
w_{jh})^{1/3}}{1 / (k_i (k_i - 1)) sum_(j,h) (w_ij w_ih w_jh)^(1/3)}
with every weight divided by the largest one (NetworkX's weighted
\code{clustering()}); both sums run over ordered pairs of neighbors. On
networks without parallel edges these are the published definitions; with
equal weights both reduce to the binary coefficient.

Results match igraph (on the simple undirected reading) and NetworkX on every
network in the validation suite.
}
\value{
\code{transitivity()} and \code{average_clustering()} return a number;
\code{clustering()} returns a numeric vector aligned with the rows of
\code{x$nodes}.
}
\references{
Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of 'small-world'
networks. \emph{Nature}, 393(6684), 440-442. \doi{10.1038/30918}

Newman, M. E. J. (2003). The structure and function of complex networks.
\emph{SIAM Review}, 45(2), 167-256. \doi{10.1137/S003614450342480}

Barrat, A., Barthélemy, M., Pastor-Satorras, R., & Vespignani, A. (2004).
The architecture of complex weighted networks. \emph{Proceedings of the
National Academy of Sciences}, 101(11), 3747-3752.
\doi{10.1073/pnas.0400087101}

Onnela, J.-P., Saramäki, J., Kertész, J., & Kaski, K. (2005). Intensity and
coherence of motifs in weighted complex networks. \emph{Physical Review E},
71(6), 065103. \doi{10.1103/PhysRevE.71.065103}

Kaiser, M. (2008). Mean clustering coefficients: The role of isolated nodes
and leafs on clustering measures for small-world networks. \emph{New Journal
of Physics}, 10(8), 083042. \doi{10.1088/1367-2630/10/8/083042}

Chiba, N., & Nishizeki, T. (1985). Arboricity and subgraph listing
algorithms. \emph{SIAM Journal on Computing}, 14(1), 210-223.
\doi{10.1137/0214017}
}
\examples{
# A triangle with a pendant node and an isolate.
g <- netframe(data.frame(from = c(1, 2, 3, 3), to = c(2, 3, 1, 4)),
              nodes = data.frame(id = 1:5), directed = FALSE)
transitivity(g)
clustering(g)
clustering(g, convention = "networkx")
average_clustering(g)
average_clustering(g, convention = "networkx")
}
