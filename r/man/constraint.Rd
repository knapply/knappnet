\name{constraint}
\alias{constraint}
\alias{effective_size}
\alias{efficiency}
\title{Structural holes: constraint, effective size, and efficiency}
\usage{
constraint(x, weighted = NULL, threads = 0L)

effective_size(x, weighted = NULL, threads = 0L)

efficiency(x, weighted = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{weighted}{Whether to use edge weights as tie strengths. \code{NULL}
(the default) uses them when the network has a weight column.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Burt's (1992) measures of how far a node's contacts are redundant with one
another. \code{constraint()} is high when a node's contacts are few or
strongly tied to each other; \code{effective_size()} is the number of
contacts, discounted for redundancy; \code{efficiency()} is effective size
divided by the number of contacts.
}
\details{
A node's contacts are its in- and out-neighbours other than itself.
Proportional tie strength is \eqn{p_{ij} = z_{ij} / \sum_k z_{ik}}{p_ij =
z_ij / sum_k z_ik} with \eqn{z_{ij} = a_{ij} + a_{ji}}{z_ij = a_ij + a_ji},
where \eqn{a_{ij}}{a_ij} is the total weight of the edges from i to j (their
number in a binary network), so parallel edges add up and a mutual pair of
arcs is a stronger tie than a single arc. An undirected network gives the
same results as the directed network with every edge reciprocated.
Self-loops are ignored.

Constraint is
\deqn{C_i = \sum_j \left(p_{ij} + \sum_{q \ne j} p_{iq} p_{qj}\right)^2}{C_i
= sum_j (p_ij + sum_{q != j} p_iq p_qj)^2}
and effective size is
\deqn{E_i = \sum_j \left(1 - \sum_{q \ne j} p_{iq} m_{jq}\right), \quad
m_{jq} = z_{jq} / \max_k z_{jk},}{E_i = sum_j (1 - sum_{q != j} p_iq m_jq),
m_jq = z_jq / max_k z_jk,}
where j and q range over the contacts of i (Burt 1992, 2004). In binary
undirected networks without parallel edges, effective size equals
\eqn{n - 2t/n} for a node with n contacts and t ties among them (Borgatti
1997). All three measures are \code{NaN} for nodes without contacts.

Weights are read as tie strengths (volume of contact) and used as given,
whatever the network's \code{weight_type}; they must be positive. These
conventions match igraph's \code{constraint()} and NetworkX's
\code{constraint()} and \code{effective_size()}. NetworkX counts
self-loops in a node's tie volume; this package, like igraph, does not.
}
\references{
Burt, R. S. (1992). \emph{Structural Holes: The Social Structure of
Competition}. Harvard University Press. \doi{10.4159/9780674029095}

Burt, R. S. (2004). Structural holes and good ideas. \emph{American Journal
of Sociology}, 110(2), 349-399. \doi{10.1086/421787}

Borgatti, S. P. (1997). Structural holes: Unpacking Burt's redundancy
measures. \emph{Connections}, 20(1), 35-38.
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 3)),
              directed = FALSE)
constraint(g)
effective_size(g)
efficiency(g)
}
