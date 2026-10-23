# Random networks (man/random_gnp.Rd).
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

# See man/random_gnp.Rd.
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

# See man/random_gnp.Rd.
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

# See man/random_gnp.Rd.
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

# See man/random_gnp.Rd.
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

# See man/random_gnp.Rd.
random_preferential <- function(n, m, reps = 1L, seed = NULL, threads = 0L) {
  n <- check_count(n, "n")
  m <- check_count(m, "m")
  args <- random_common(reps, seed, threads)
  lists <- engine_call("random_preferential",
                       .Call(C_random_preferential, n, m, args$reps, args$seed, args$threads),
                       edge_rows = FALSE)
  random_result(lists, n, FALSE)
}

# See man/random_gnp.Rd.
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
