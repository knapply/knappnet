\name{random_gnp}
\alias{random_gnp}
\alias{random_gnm}
\alias{random_configuration}
\alias{random_uman}
\alias{random_small_world}
\alias{random_preferential}
\alias{random_sbm}
\title{Random networks}
\usage{
random_gnp(n, p, directed = FALSE, loops = FALSE, reps = 1L, seed = NULL,
           threads = 0L)

random_gnm(n, m, directed = FALSE, loops = FALSE, reps = 1L, seed = NULL,
           threads = 0L)

random_configuration(degrees, in_degrees = NULL, simple = FALSE, max_tries = 1000L,
                     reps = 1L, seed = NULL, threads = 0L)

random_uman(n, mutual, asym, reps = 1L, seed = NULL, threads = 0L)

random_small_world(n, k, p, reps = 1L, seed = NULL, threads = 0L)

random_preferential(n, m, reps = 1L, seed = NULL, threads = 0L)

random_sbm(sizes, probs, directed = FALSE, loops = FALSE, reps = 1L, seed = NULL,
           threads = 0L)
}
\arguments{
\item{n}{Number of nodes.}
\item{p}{Tie probability (\code{random_gnp()}) or rewiring probability
(\code{random_small_world()}).}
\item{m}{Number of edges (\code{random_gnm()}) or edges added by each new
node (\code{random_preferential()}, at least 1 and less than \code{n}).}
\item{directed}{Whether to draw a directed network.}
\item{loops}{Whether self-loops are possible.}
\item{degrees}{Degree of every node (undirected), or out-degree of every
node when \code{in_degrees} is given (directed).}
\item{in_degrees}{In-degree of every node, for a directed network; the sums
of \code{degrees} and \code{in_degrees} must agree.}
\item{simple}{Whether to repeat the matching until the network is simple.}
\item{max_tries}{Attempts at a simple matching before giving up.}
\item{mutual, asym}{Numbers of mutual and asymmetric dyads; the remaining
dyads are null.}
\item{k}{Neighbors per node in the ring lattice (even, less than
\code{n}).}
\item{sizes}{Block sizes (nodes are numbered block by block).}
\item{probs}{A square matrix of tie probabilities between blocks, symmetric
for undirected networks.}
\item{reps}{Number of networks to draw.}
\item{seed}{Seed, a whole number between 0 and \eqn{2^{53}}; \code{NULL}
draws one from R's random number generator, so \code{\link{set.seed}()}
makes results reproducible.}
\item{threads}{Number of threads for drawing replicates in parallel, a whole
number; \code{0} uses all available cores, or 2 during \verb{R CMD check}
(see \sQuote{Threads} in \code{\link{netframe}}).}
}
\value{
A netframe with nodes \code{1, 2, ..., n} (integer identifiers in column
\code{id}) when \code{reps = 1}, otherwise a list of \code{reps} netframes.
}
\description{
Seeded, reproducible random networks from the classic models:
\describe{
  \item{\code{random_gnp()}}{\eqn{G(n, p)}: each possible pair of nodes is
    tied independently with probability \code{p} (Gilbert 1959).}
  \item{\code{random_gnm()}}{\eqn{G(n, m)}: a network drawn uniformly from
    those with exactly \code{m} edges (Erdős and Rényi 1959).}
  \item{\code{random_configuration()}}{The configuration model: the edge
    ends ("stubs") of a degree sequence are matched uniformly at random
    (Bollobás 1980; Molloy and Reed 1995; Newman, Strogatz, and Watts
    2001).}
  \item{\code{random_uman()}}{U|MAN: a directed network drawn uniformly from
    those with given numbers of mutual, asymmetric, and null dyads (Holland
    and Leinhardt 1976), the null model of \code{\link{cug_test}} with
    \code{conditioning = "dyad.census"}.}
  \item{\code{random_small_world()}}{The Watts-Strogatz (1998) model: a ring
    lattice joining each node to its \code{k} nearest neighbors, each edge
    rewired with probability \code{p}.}
  \item{\code{random_preferential()}}{The Barabási-Albert (1999) model: from
    a star on \code{m + 1} nodes, each new node joins \code{m} existing
    nodes drawn with probability proportional to their degree.}
  \item{\code{random_sbm()}}{A stochastic block model (Holland, Laskey, and
    Leinhardt 1983): nodes in blocks of the given sizes, a node of block
    \eqn{r} and a node of block \eqn{s} tied independently with probability
    \code{probs[r, s]}.}
}
}
\details{
Unless stated otherwise the networks are simple: no self-loops, no parallel
edges. \code{random_gnp()} and \code{random_gnm()} can include self-loops on
request (\code{n} more possible edges, as igraph's \code{sample_gnp(loops =
TRUE)}). \code{random_configuration()} returns the multigraph by default
(loops and parallel edges kept, so every node has exactly its degree, a
loop counting 2 in an undirected network); with \code{simple = TRUE} the
whole matching is repeated until it has no loop or parallel edge, which
gives every simple network with the degree sequence the same probability,
and stops with an error after \code{max_tries} attempts for sequences whose
simple matchings are rare. Erasing loops and parallel edges instead would
change the degrees and is not offered.

In the Watts-Strogatz model each lattice edge is considered once, lap by lap
(\code{j = 1, ..., k/2} around the ring), and with probability \code{p} its
far end is moved to a node drawn uniformly from those not already adjacent,
as in NetworkX's \code{watts_strogatz_graph()}; the network stays simple
with \code{n k / 2} edges. The preferential attachment model starts from a
star on \code{m + 1} nodes with node 1 at the center, as NetworkX's
\code{barabasi_albert_graph()}, and has \code{m (n - m)} edges. \eqn{G(n,
p)} and block models are sampled by geometric skipping in time
proportional to the output, and \eqn{G(n, m)} by Floyd's sampling without
replacement (Batagelj and Brandes 2005; Bentley and Floyd 1987).

Replicate \eqn{r} of a batch uses random stream \code{(seed, r)} of a
counter-based generator, so a seed reproduces every replicate on any
platform and thread count, and the first replicate of a batch of
\code{reps} is the network \code{reps = 1} gives. The exact constraints
(edge counts, degrees, dyad censuses) are checked by the validation suite
(\code{validation/cug_check.py}).
}
\references{
Gilbert, E. N. (1959). Random graphs. \emph{Annals of Mathematical
Statistics}, 30(4), 1141-1144. \doi{10.1214/aoms/1177706098}

Erdős, P., & Rényi, A. (1959). On random graphs I. \emph{Publicationes
Mathematicae Debrecen}, 6, 290-297.

Bollobás, B. (1980). A probabilistic proof of an asymptotic formula for the
number of labelled regular graphs. \emph{European Journal of
Combinatorics}, 1(4), 311-316. \doi{10.1016/S0195-6698(80)80030-8}

Molloy, M., & Reed, B. (1995). A critical point for random graphs with a
given degree sequence. \emph{Random Structures & Algorithms}, 6(2-3),
161-180. \doi{10.1002/rsa.3240060204}

Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2001). Random graphs
with arbitrary degree distributions and their applications. \emph{Physical
Review E}, 64(2), 026118. \doi{10.1103/PhysRevE.64.026118}

Holland, P. W., & Leinhardt, S. (1976). Local structure in social networks.
\emph{Sociological Methodology}, 7, 1-45. \doi{10.2307/270703}

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
'small-world' networks. \emph{Nature}, 393(6684), 440-442.
\doi{10.1038/30918}

Barabási, A.-L., & Albert, R. (1999). Emergence of scaling in random
networks. \emph{Science}, 286(5439), 509-512.
\doi{10.1126/science.286.5439.509}

Holland, P. W., Laskey, K. B., & Leinhardt, S. (1983). Stochastic
blockmodels: First steps. \emph{Social Networks}, 5(2), 109-137.
\doi{10.1016/0378-8733(83)90021-7}

Batagelj, V., & Brandes, U. (2005). Efficient generation of large random
networks. \emph{Physical Review E}, 71(3), 036113.
\doi{10.1103/PhysRevE.71.036113}

Bentley, J., & Floyd, B. (1987). Programming pearls: A sample of
brilliance. \emph{Communications of the ACM}, 30(9), 754-757.
\doi{10.1145/30401.315746}
}
\seealso{\code{\link{cug_test}}, \code{\link{small_world}}}
\examples{
g <- random_gnp(20, 0.2, seed = 1, threads = 2)
g
nrow(random_gnm(20, 30, seed = 1, threads = 2)$edges)

# A configuration-model network with the karate club's degrees.
karate <- example_network("karate")
cm <- random_configuration(degree(karate), simple = TRUE, seed = 1, threads = 2)
all(degree(cm) == degree(karate))

# U|MAN with 5 mutual and 10 asymmetric dyads on 8 nodes.
reciprocity(random_uman(8, 5, 10, seed = 1, threads = 2))

sw <- random_small_world(30, 4, 0.1, seed = 1, threads = 2)
ba <- random_preferential(30, 2, seed = 1, threads = 2)
blocks <- random_sbm(c(10, 10), matrix(c(0.5, 0.05, 0.05, 0.5), 2), seed = 1, threads = 2)

# Several replicates at once; replicate 1 is the single draw.
reps <- random_gnm(10, 15, reps = 3, seed = 7, threads = 2)
identical(reps[[1]]$edges, random_gnm(10, 15, seed = 7, threads = 2)$edges)
}
