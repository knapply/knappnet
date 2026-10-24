\name{net_spinglass}
\alias{net_spinglass}
\title{Spinglass communities}
\description{
Finds communities as the ground state of a Potts model, approached by
simulated annealing (Reichardt and Bornholdt 2006), with the extension to
negative ties of Traag and Bruggeman (2009). Its arguments follow
igraph's \code{cluster_spinglass()}.
}
\usage{
net_spinglass(x, spins = 25L, parupdate = FALSE, start_temp = 1, stop_temp = 0.01,
              cool_fact = 0.99, update_rule = c("config", "simple"), gamma = 1,
              implementation = c("orig", "neg"), gamma_minus = 1, weighted = NULL,
              seed = NULL)
}
\arguments{
\item{x}{A connected \code{\link{net_frame}}; directions are ignored.}
\item{spins}{Largest number of groups (spin states), a positive whole
number, capped at the number of nodes.}
\item{parupdate}{Whether to update all nodes at once in each sweep, from
the previous configuration, instead of one at a time in random order. Not
available with \code{implementation = "neg"}.}
\item{start_temp, stop_temp}{Starting and stopping temperatures, in units of
edge weight; \code{stop_temp} must not exceed \code{start_temp}.}
\item{cool_fact}{Cooling factor, strictly between 0 and 1: the temperature
is multiplied by it after every sweep.}
\item{update_rule}{Null model: \code{"config"} (the configuration model,
for which, when \code{gamma = 1}, minus the energy is
modularity times the total weight, up to a constant) or \code{"simple"} (a random network of the same
density).}
\item{gamma}{Resolution parameter \eqn{\gamma} of the positive ties
(non-negative).}
\item{implementation}{\code{"orig"} (Reichardt and Bornholdt, positive
weights only) or \code{"neg"} (Traag and Bruggeman, which also accepts
negative weights).}
\item{gamma_minus}{Resolution parameter of the negative ties, used only
with \code{implementation = "neg"} (non-negative).}
\item{weighted}{Whether to use edge weights, as given. \code{NULL} (the
default) uses them when the network has a weight column. With
\code{"orig"} they must not be negative.}
\item{seed}{Seed for the starting spins and the sweeps, a whole number
between 0 and \eqn{2^{53}}; \code{NULL} draws one from R's random number
generator, so \code{\link{set.seed}()} makes results reproducible.}
}
\details{
With spin \eqn{s_i} for each node, the method minimizes
\deqn{H = - \sum_{i < j} (A_{ij} - \gamma p_{ij}) \delta(s_i, s_j),}{
H = - sum_{i<j} (A_ij - gamma p_ij) delta(s_i, s_j),}
with \eqn{p_{ij} = k_i k_j / 2m} for the configuration model and
\eqn{p_{ij} = 2m / (n (n - 1))} for the simple one (\eqn{k_i} the weighted
degree, \eqn{m} the total weight). The signed form adds, for negative
weights, a second term with their magnitudes, their own null model and
degrees, and \code{gamma_minus}. Updates are heat-bath (Gibbs) updates: a
node moves to spin \eqn{s} with probability proportional to
\eqn{\exp(-\Delta H(s) / T)}; \eqn{T} starts at \code{start_temp} and is
multiplied by \code{cool_fact} after each sweep until it falls to
\code{stop_temp}. A last zero-temperature pass moves nodes to their best
spin until no move lowers \eqn{H}, so the result is a local minimum.

The network is taken as its underlying simple undirected network: parallel
edges and reciprocal arcs add their weights and self-loops are ignored.
Disconnected networks are an error, as in igraph; analyze each component
on its own. igraph does not document its annealing schedule, so partitions
are not expected to match igraph's, only in quality (the validation suite
compares the mean modularity over many seeds); the same \code{seed} gives
the same result on every platform.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with the
attributes
\describe{
  \item{\code{"modularity"}}{Modularity of the partition on the underlying
    undirected network; \code{NaN} when a weight is negative or the network
    has no edges.}
  \item{\code{"energy"}}{The energy \eqn{H} of the final configuration.}
}
}
\references{
Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community
detection. \emph{Physical Review E}, 74(1), 016110.
\doi{10.1103/PhysRevE.74.016110}

Traag, V. A., & Bruggeman, J. (2009). Community detection in networks with
positive and negative links. \emph{Physical Review E}, 80(3), 036115.
\doi{10.1103/PhysRevE.80.036115}

Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by
simulated annealing. \emph{Science}, 220(4598), 671--680.
\doi{10.1126/science.220.4598.671}
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_modularity}},
\code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
# Null model (igraph update.rule) and implementation, one call each.
net_spinglass(g, update_rule = "config", seed = 1)   # igraph default: "config", "orig"
net_spinglass(g, update_rule = "simple", seed = 1)   # igraph update.rule = "simple"
groups <- net_spinglass(g, implementation = "neg", seed = 1)   # igraph implementation = "neg"
net_plot_partition(g, groups, seed = 1)

# Fewer spin states allow fewer groups.
max(net_spinglass(g, spins = 3, seed = 1))
}
