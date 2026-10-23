#' Random networks
#'
#' Seeded, reproducible random networks from the classic models, returned
#' as netframes with nodes `1, 2, ..., n`:
#' \describe{
#'   \item{`random_gnp()`}{\eqn{G(n, p)}: each possible pair of nodes is
#'     tied independently with probability `p` (Gilbert 1959).}
#'   \item{`random_gnm()`}{\eqn{G(n, m)}: a network drawn uniformly from
#'     those with exactly `m` edges (Erdős and Rényi 1959).}
#'   \item{`random_configuration()`}{The configuration model: the edge ends
#'     ("stubs") of a degree sequence are matched uniformly at random
#'     (Bollobás 1980; Molloy and Reed 1995; Newman, Strogatz, and Watts
#'     2001).}
#'   \item{`random_uman()`}{U|MAN: a directed network drawn uniformly from
#'     those with given numbers of mutual, asymmetric, and null dyads
#'     (Holland and Leinhardt 1976), the null model of [cug_test()] with
#'     `conditioning = "dyad.census"`.}
#'   \item{`random_small_world()`}{The Watts-Strogatz (1998) model: a ring
#'     lattice joining each node to its `k` nearest neighbors, each edge
#'     rewired with probability `p`.}
#'   \item{`random_preferential()`}{The Barabási-Albert (1999) model: from a
#'     star on `m + 1` nodes, each new node joins `m` existing nodes drawn
#'     with probability proportional to their degree.}
#'   \item{`random_sbm()`}{A stochastic block model (Holland, Laskey, and
#'     Leinhardt 1983): nodes in blocks of the given sizes, a node of block
#'     \eqn{r} and a node of block \eqn{s} tied independently with
#'     probability `probs[r, s]`.}
#' }
#'
#' @details
#' Unless stated otherwise the networks are simple: no self-loops, no
#' parallel edges. `random_gnp()` and `random_gnm()` can include self-loops
#' on request (`n` more possible edges, as igraph's `sample_gnp(loops =
#' TRUE)`). `random_configuration()` returns the multigraph by default
#' (loops and parallel edges kept, so every node has exactly its degree, a
#' loop counting 2 in an undirected network); with `simple = TRUE` the whole
#' matching is repeated until it has no loop or parallel edge, which gives
#' every simple network with the degree sequence the same probability, and
#' stops with an error after `max_tries` attempts for sequences whose simple
#' matchings are rare. Erasing loops and parallel edges instead would change
#' the degrees and is not offered.
#'
#' In the Watts-Strogatz model each lattice edge is considered once, in lap
#' order (`j = 1, ..., k/2` around the ring), and with probability `p` its
#' far end is moved to a node drawn uniformly from those not already
#' adjacent, as in NetworkX's `watts_strogatz_graph()`; the network stays
#' simple with `n k / 2` edges. The preferential attachment model starts
#' from a star on `m + 1` nodes with node 1 at the center, as NetworkX's
#' `barabasi_albert_graph()`, and has `m (n - m)` edges. \eqn{G(n, p)} and
#' block models are sampled by geometric skipping in time proportional to
#' the output, and \eqn{G(n, m)} by Floyd's sampling without replacement
#' (Batagelj and Brandes 2005; Bentley and Floyd 1987).
#'
#' Replicate \eqn{r} of a batch uses random stream `(seed, r)` of a
#' counter-based generator, so a seed reproduces every replicate on any
#' platform and thread count, and the first replicate of a batch of `reps`
#' is the network `reps = 1` gives.
#'
#' @param n Number of nodes.
#' @param p Tie probability (`random_gnp()`) or rewiring probability
#'   (`random_small_world()`).
#' @param m Number of edges (`random_gnm()`) or edges added by each new node
#'   (`random_preferential()`, at least 1 and less than `n`).
#' @param directed Whether to draw a directed network.
#' @param loops Whether self-loops are possible.
#' @param degrees Degree of every node (undirected), or out-degree of every
#'   node when `in_degrees` is given (directed).
#' @param in_degrees In-degree of every node, for a directed network; the
#'   sums of `degrees` and `in_degrees` must agree.
#' @param simple Whether to repeat the matching until the network is simple.
#' @param max_tries Attempts at a simple matching before giving up.
#' @param mutual,asym Numbers of mutual and asymmetric dyads; the remaining
#'   dyads are null.
#' @param k Neighbors per node in the ring lattice (even, less than `n`).
#' @param sizes Block sizes (nodes are numbered block by block).
#' @param probs A square matrix of tie probabilities between blocks,
#'   symmetric for undirected networks.
#' @param reps Number of networks to draw.
#' @param seed Seed, a whole number between 0 and \eqn{2^{53}}; `NULL`
#'   draws one from R's random number generator, so [set.seed()] makes
#'   results reproducible.
#' @param threads Worker threads for drawing replicates in parallel; see
#'   [netframe()].
#'
#' @references
#' Gilbert, E. N. (1959). Random graphs. *Annals of Mathematical
#' Statistics*, 30(4), 1141-1144. \doi{10.1214/aoms/1177706098}
#'
#' Erdős, P., & Rényi, A. (1959). On random graphs I. *Publicationes
#' Mathematicae Debrecen*, 6, 290-297.
#'
#' Bollobás, B. (1980). A probabilistic proof of an asymptotic formula for
#' the number of labelled regular graphs. *European Journal of
#' Combinatorics*, 1(4), 311-316. \doi{10.1016/S0195-6698(80)80030-8}
#'
#' Molloy, M., & Reed, B. (1995). A critical point for random graphs with a
#' given degree sequence. *Random Structures & Algorithms*, 6(2-3), 161-180.
#' \doi{10.1002/rsa.3240060204}
#'
#' Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2001). Random graphs
#' with arbitrary degree distributions and their applications. *Physical
#' Review E*, 64(2), 026118. \doi{10.1103/PhysRevE.64.026118}
#'
#' Holland, P. W., & Leinhardt, S. (1976). Local structure in social
#' networks. *Sociological Methodology*, 7, 1-45. \doi{10.2307/270703}
#'
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#' 'small-world' networks. *Nature*, 393(6684), 440-442.
#' \doi{10.1038/30918}
#'
#' Barabási, A.-L., & Albert, R. (1999). Emergence of scaling in random
#' networks. *Science*, 286(5439), 509-512.
#' \doi{10.1126/science.286.5439.509}
#'
#' Holland, P. W., Laskey, K. B., & Leinhardt, S. (1983). Stochastic
#' blockmodels: First steps. *Social Networks*, 5(2), 109-137.
#' \doi{10.1016/0378-8733(83)90021-7}
#'
#' Batagelj, V., & Brandes, U. (2005). Efficient generation of large random
#' networks. *Physical Review E*, 71(3), 036113.
#' \doi{10.1103/PhysRevE.71.036113}
#'
#' Bentley, J., & Floyd, B. (1987). Programming pearls: A sample of
#' brilliance. *Communications of the ACM*, 30(9), 754-757.
#' \doi{10.1145/30401.315746}
#'
#' @return A netframe with nodes `1, 2, ..., n` (integer identifiers in
#'   column `id`) when `reps = 1`, otherwise a list of `reps` netframes.
#' @export
random_gnp <- function(n, p, directed = FALSE, loops = FALSE, reps = 1L, seed = NULL,
                       threads = 0L) {
  n <- check_count(n, "n")
  p <- check_probability(p, "p")
  directed <- check_flag(directed, "directed")
  loops <- check_flag(loops, "loops")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_gnp",
                       .Call(C_random_gnp, n, p, directed, loops, args$reps, args$seed,
                             args$threads), edge_rows = FALSE)
  random_result(lists, n, directed)
}

#' @rdname random_gnp
#' @export
random_gnm <- function(n, m, directed = FALSE, loops = FALSE, reps = 1L, seed = NULL,
                       threads = 0L) {
  n <- check_count(n, "n")
  m <- check_count(m, "m")
  directed <- check_flag(directed, "directed")
  loops <- check_flag(loops, "loops")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_gnm",
                       .Call(C_random_gnm, n, m, directed, loops, args$reps, args$seed,
                             args$threads), edge_rows = FALSE)
  random_result(lists, n, directed)
}

#' @rdname random_gnp
#' @export
random_configuration <- function(degrees, in_degrees = NULL, simple = FALSE, max_tries = 1000L,
                                 reps = 1L, seed = NULL, threads = 0L) {
  degrees <- check_degrees(degrees, "degrees")
  if (!is.null(in_degrees)) {
    in_degrees <- check_degrees(in_degrees, "in_degrees")
    if (length(in_degrees) != length(degrees)) {
      stop("`degrees` and `in_degrees` must have the same length.", call. = FALSE)
    }
  }
  simple <- check_flag(simple, "simple")
  max_tries <- check_count(max_tries, "max_tries")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_configuration",
                       .Call(C_random_configuration, degrees, in_degrees, simple, max_tries,
                             args$reps, args$seed, args$threads), edge_rows = FALSE)
  random_result(lists, length(degrees), !is.null(in_degrees))
}

#' @rdname random_gnp
#' @export
random_uman <- function(n, mutual, asym, reps = 1L, seed = NULL, threads = 0L) {
  n <- check_count(n, "n")
  mutual <- check_count(mutual, "mutual")
  asym <- check_count(asym, "asym")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_uman",
                       .Call(C_random_uman, n, mutual, asym, args$reps, args$seed,
                             args$threads), edge_rows = FALSE)
  random_result(lists, n, TRUE)
}

#' @rdname random_gnp
#' @export
random_small_world <- function(n, k, p, reps = 1L, seed = NULL, threads = 0L) {
  n <- check_count(n, "n")
  k <- check_count(k, "k")
  p <- check_probability(p, "p")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_small_world",
                       .Call(C_random_small_world, n, k, p, args$reps, args$seed,
                             args$threads), edge_rows = FALSE)
  random_result(lists, n, FALSE)
}

#' @rdname random_gnp
#' @export
random_preferential <- function(n, m, reps = 1L, seed = NULL, threads = 0L) {
  n <- check_count(n, "n")
  m <- check_count(m, "m")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_preferential",
                       .Call(C_random_preferential, n, m, args$reps, args$seed, args$threads),
                       edge_rows = FALSE)
  random_result(lists, n, FALSE)
}

#' @rdname random_gnp
#' @export
random_sbm <- function(sizes, probs, directed = FALSE, loops = FALSE, reps = 1L, seed = NULL,
                       threads = 0L) {
  sizes <- check_degrees(sizes, "sizes")
  if (!is.matrix(probs) || !is.numeric(probs) || nrow(probs) != ncol(probs) ||
      nrow(probs) != length(sizes)) {
    stop("`probs` must be a square numeric matrix with one row per block.", call. = FALSE)
  }
  if (anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop("`probs` must hold probabilities between 0 and 1.", call. = FALSE)
  }
  directed <- check_flag(directed, "directed")
  loops <- check_flag(loops, "loops")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_sbm",
                       .Call(C_random_sbm, sizes, matrix(as.double(probs), nrow(probs)),
                             directed, loops, args$reps, args$seed, args$threads),
                       edge_rows = FALSE)
  random_result(lists, sum(sizes), directed)
}

random_common <- function(reps, seed, threads) {
  reps <- check_count(reps, "reps")
  if (reps < 1L) stop("`reps` must be at least 1.", call. = FALSE)
  list(reps = reps, seed = groups_seed(seed), threads = resolve_threads(threads))
}

check_probability <- function(p, what) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
    stop(sprintf("`%s` must be a single probability between 0 and 1.", what), call. = FALSE)
  }
  as.double(p)
}

check_degrees <- function(d, what) {
  if (!is.numeric(d) || anyNA(d) || any(d < 0) || any(d != floor(d)) ||
      any(d > .Machine$integer.max)) {
    stop(sprintf("`%s` must be a vector of non-negative whole numbers.", what),
         call. = FALSE)
  }
  as.integer(d)
}

# Turns the engine's list of list(from, to) (1-based node positions) into
# netframes on nodes 1..n.
random_result <- function(lists, n, directed) {
  nodes <- data.frame(id = seq_len(n))
  nets <- lapply(lists, function(e) {
    netframe(data.frame(from = e[[1L]], to = e[[2L]]), nodes, directed = directed)
  })
  if (length(nets) == 1L) nets[[1L]] else nets
}
