#' Bicomponents (blocks)
#'
#' A bicomponent is a maximal set of nodes that stays connected when any one
#' of its nodes is removed; two bicomponents share at most one node, a cut
#' point (Wasserman and Faust 1994, section 4.2). Found with the depth-first
#' search of Hopcroft and Tarjan (1973) in time linear in the size of the
#' network.
#'
#' The undirected simple graph underlying the network is used: an edge in
#' either direction joins two nodes, parallel edges count once, self-loops
#' and weights are ignored. A bridge and its two endpoints form a
#' bicomponent of two nodes, so every edge lies in exactly one bicomponent;
#' isolated nodes belong to none. This matches igraph's
#' `biconnected_components()` and NetworkX's `biconnected_components()`,
#' and UCINET's Network > Regions > Bi-Component routine. See
#' [cut_points()] for the cut points themselves.
#'
#' @inheritParams graph_density
#'
#' @references
#' Hopcroft, J., & Tarjan, R. (1973). Algorithm 447: Efficient algorithms
#' for graph manipulation. *Communications of the ACM*, 16(6), 372-378.
#' \doi{10.1145/362248.362272}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return A list of vectors of node identifiers (values of
#'   `x$nodes[[id]]`), each bicomponent's members in node order, in
#'   lexicographic order of member positions.
#' @export
bicomponents <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  idx <- engine_call("bicomponents", .Call(C_bicomponents, ptr))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

#' n-cliques, n-clans, and k-plexes
#'
#' Relaxations of the clique (see [max_cliques()]) that admit groups whose
#' members are not all directly tied:
#' \describe{
#'   \item{n-clique}{A maximal set of nodes every two of which are joined by
#'     a path of length at most `n` in the whole network (Luce 1950; Alba
#'     1973). The paths may run through nodes outside the set. 1-cliques are
#'     cliques.}
#'   \item{n-clan}{An n-clique whose own subgraph has diameter at most `n`,
#'     so that the short paths stay inside the group (Mokken 1979).}
#'   \item{k-plex}{A maximal set of nodes each of which is adjacent to all
#'     but at most `k - 1` of the other members (Seidman and Foster 1978).
#'     1-plexes are cliques.}
#' }
#'
#' @details
#' All three use the undirected simple graph underlying the network:
#' directions, parallel edges, self-loops, and weights are ignored, as in
#' UCINET's Network > Subgroups routines for symmetric data. Members are
#' listed in node order and groups in lexicographic order of member
#' positions.
#'
#' n-cliques are the maximal cliques of the graph joining two nodes when
#' their distance is at most `n`, which is built with one breadth-first
#' search per node (in parallel) and searched with the Bron-Kerbosch
#' algorithm; n-clans are the n-cliques that pass the diameter test. With
#' `min_size` at most 1 an isolated node is an n-clique and n-clan of its
#' own. `min_size` only drops smaller groups; it never reports non-maximal
#' subsets.
#'
#' k-plexes are listed with a Bron-Kerbosch set enumeration (every subset of
#' a k-plex is a k-plex, which makes the scheme exact), after dropping the
#' nodes with fewer than `min_size - k` neighbors, which cannot belong to a
#' k-plex of `min_size` nodes (Seidman and Foster 1978). Any set of at most
#' `k` nodes is a k-plex, so small values of `min_size` report many trivial
#' sets; `min_size` below `2k - 1` also admits disconnected k-plexes. The
#' default, `2k - 1`, is the smallest size at which every k-plex is
#' connected with diameter at most 2. The number of maximal k-plexes can
#' grow exponentially with the network, so the search stops with an error
#' after `max_steps` steps rather than running indefinitely; raise the
#' limit, `k`, or `min_size` as needed. Whether the limit is reached does
#' not depend on `threads`.
#'
#' @inheritParams betweenness
#' @param n Largest path length between members: 1 for cliques, 2 for
#'   2-cliques, and so on.
#' @param k Each member may have at most `k - 1` non-neighbors among the
#'   other members (`k >= 1`).
#' @param min_size Smallest group to report. For k-plexes `NULL` means
#'   `2k - 1`.
#' @param max_steps Limit on the number of search steps for k-plexes.
#'
#' @references
#' Luce, R. D. (1950). Connectivity and generalized cliques in sociometric
#' group structure. *Psychometrika*, 15(2), 169-190.
#' \doi{10.1007/BF02289199}
#'
#' Alba, R. D. (1973). A graph-theoretic definition of a sociometric
#' clique. *Journal of Mathematical Sociology*, 3(1), 113-126.
#' \doi{10.1080/0022250X.1973.9989826}
#'
#' Mokken, R. J. (1979). Cliques, clubs and clans. *Quality and Quantity*,
#' 13(2), 161-173. \doi{10.1007/BF00139635}
#'
#' Seidman, S. B., & Foster, B. L. (1978). A graph-theoretic generalization
#' of the clique concept. *Journal of Mathematical Sociology*, 6(1),
#' 139-154. \doi{10.1080/0022250X.1978.9989883}
#'
#' Bron, C., & Kerbosch, J. (1973). Algorithm 457: Finding all cliques of an
#' undirected graph. *Communications of the ACM*, 16(9), 575-577.
#' \doi{10.1145/362342.362367}
#'
#' @return A list of vectors of node identifiers (values of
#'   `x$nodes[[id]]`), one per group.
#' @export
n_cliques <- function(x, n = 2L, min_size = 1L, threads = 0L) {
  check_netframe(x)
  n <- check_count(n, "n")
  if (n < 1L) stop("`n` must be at least 1.", call. = FALSE)
  min_size <- check_count(min_size, "min_size")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("n_cliques", .Call(C_n_cliques, ptr, n, min_size, FALSE, threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

#' @rdname n_cliques
#' @export
n_clans <- function(x, n = 2L, min_size = 1L, threads = 0L) {
  check_netframe(x)
  n <- check_count(n, "n")
  if (n < 1L) stop("`n` must be at least 1.", call. = FALSE)
  min_size <- check_count(min_size, "min_size")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("n_clans", .Call(C_n_cliques, ptr, n, min_size, TRUE, threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

#' @rdname n_cliques
#' @export
k_plexes <- function(x, k = 2L, min_size = NULL, max_steps = 1e7, threads = 0L) {
  check_netframe(x)
  k <- check_count(k, "k")
  if (k < 1L) stop("`k` must be at least 1.", call. = FALSE)
  min_size <- if (is.null(min_size)) 2L * k - 1L else check_count(min_size, "min_size")
  if (!is.numeric(max_steps) || length(max_steps) != 1L || is.na(max_steps) || max_steps < 1) {
    stop("`max_steps` must be a number of at least 1.", call. = FALSE)
  }
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("k_plexes",
                     .Call(C_k_plexes, ptr, k, min_size, as.double(max_steps), threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

#' Lambda sets and edge connectivity
#'
#' The edge connectivity \eqn{\lambda(u, v)} of two nodes is the smallest
#' number of edges whose removal disconnects them, which equals the number
#' of edge-disjoint paths between them (Menger's theorem; Ford and Fulkerson
#' 1956). A lambda set is a set of nodes any two of which are more
#' edge-connected to each other than either is to any node outside
#' (Borgatti, Everett, and Shirey 1990). Because
#' \eqn{\lambda(a, c) \ge \min(\lambda(a, b), \lambda(b, c))}, the
#' relation "\eqn{\lambda(a, b) \ge k}" is an equivalence for every
#' \eqn{k \ge 1}, and its classes with at least two members are exactly the
#' lambda sets; increasing \eqn{k} gives the nested hierarchy that UCINET's
#' Network > Subgroups > Lambda Sets routine reports.
#'
#' All pairwise connectivities come from \eqn{n - 1} maximum-flow
#' computations with Gusfield's (1990) method (a Gomory-Hu tree; Gomory and
#' Hu 1961), each flow found by breadth-first augmenting paths (Edmonds and
#' Karp 1972). The undirected simple graph underlying the network is used:
#' directions, parallel edges, self-loops, and weights are ignored. Nodes
#' in different components have connectivity 0. The connectivity matrix
#' matches igraph's `edge_connectivity(g, u, v)` and NetworkX's
#' `edge_connectivity(G, u, v)` on the simple undirected graph.
#'
#' @inheritParams graph_density
#'
#' @references
#' Borgatti, S. P., Everett, M. G., & Shirey, P. R. (1990). LS sets, lambda
#' sets and other cohesive subsets. *Social Networks*, 12(4), 337-357.
#' \doi{10.1016/0378-8733(90)90014-Z}
#'
#' Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
#' *Canadian Journal of Mathematics*, 8, 399-404.
#' \doi{10.4153/CJM-1956-045-5}
#'
#' Gomory, R. E., & Hu, T. C. (1961). Multi-terminal network flows.
#' *Journal of the Society for Industrial and Applied Mathematics*, 9(4),
#' 551-570. \doi{10.1137/0109047}
#'
#' Gusfield, D. (1990). Very simple methods for all pairs network flow
#' analysis. *SIAM Journal on Computing*, 19(1), 143-155.
#' \doi{10.1137/0219009}
#'
#' Edmonds, J., & Karp, R. M. (1972). Theoretical improvements in
#' algorithmic efficiency for network flow problems. *Journal of the ACM*,
#' 19(2), 248-264. \doi{10.1145/321694.321699}
#'
#' @return A list of class `lambda_sets` with `connectivity`, the
#'   \eqn{n \times n} matrix of pairwise edge connectivities (zero diagonal)
#'   with node identifiers as dimnames, and `sets`, a list of the lambda
#'   sets from the lowest to the highest connectivity: each element is a
#'   list with `lambda` (the connectivity \eqn{k}) and `members` (node
#'   identifiers). A set that is the same at several levels is listed once,
#'   with `lambda` the highest connectivity its members share. It has
#'   `print` and `plot` (dendrogram) methods.
#' @export
lambda_sets <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  m <- engine_call("lambda_sets", .Call(C_lambda_connectivity, ptr))
  ids <- x$nodes[[x$roles$id]]
  dimnames(m) <- list(ids, ids)
  n <- nrow(m)
  sets <- list()
  seen <- character()
  top <- if (n) max(m) else 0
  for (k in seq_len(top)) {
    done <- logical(n)
    for (v in seq_len(n)) {
      if (done[v]) next
      members <- which(m[v, ] >= k | seq_len(n) == v)
      done[members] <- TRUE
      if (length(members) < 2L) next
      key <- paste(members, collapse = " ")
      if (key %in% seen) {
        # The same set at a higher level: record its true connectivity.
        sets[[match(key, seen)]]$lambda <- k
        next
      }
      seen <- c(seen, key)
      sets[[length(sets) + 1L]] <- list(lambda = k, members = ids[members])
    }
  }
  structure(list(connectivity = m, sets = sets), class = "lambda_sets")
}

#' @export
print.lambda_sets <- function(x, ...) {
  cat(sprintf("Lambda sets of %d nodes: %d sets\n", nrow(x$connectivity), length(x$sets)))
  for (s in x$sets) {
    members <- format(s$members)
    shown <- if (length(members) > 12L) {
      paste0(paste(members[1:12], collapse = ", "), ", ... (", length(members), " nodes)")
    } else {
      paste(members, collapse = ", ")
    }
    cat(sprintf("  lambda = %d: %s\n", s$lambda, shown))
  }
  invisible(x)
}

#' Factions
#'
#' Partitions the nodes into `k` groups that best match the ideal of
#' cohesive, mutually separate groups, as UCINET's Network > Subgroups >
#' Factions routine does (Borgatti, Everett, and Freeman 2002; Hanneman and
#' Riddle 2005, chapter 11). The ideal pattern has every pair within a group
#' tied and no pair across groups tied. With `fit = "errors"` the partition
#' minimizes the number of errors, absent ties within groups plus present
#' ties between groups (the Hamming distance to the ideal, the criterion of
#' a binary blockmodel with complete diagonal blocks; Doreian, Batagelj, and
#' Ferligoj 2005); with `fit = "correlation"` it maximizes the Pearson
#' correlation between the adjacency matrix and the ideal pattern (the
#' criterion Borgatti and Everett 1999 use for core-periphery models).
#'
#' @details
#' The adjacency matrix is binary: parallel edges count once, self-loops are
#' ignored, and weights are not used. Undirected networks compare each
#' unordered pair once; directed networks compare both ordered pairs.
#' Every group is non-empty (`1 <= k <= n`). The correlation is undefined
#' when the network has no ties or every possible tie, when `k = 1`, or
#' when `k = n`; `fit = "correlation"` then stops with an error, and the
#' `correlation` attribute of an errors fit is `NaN`.
#'
#' Finding the best partition is a combinatorial optimization problem.
#' UCINET's documentation does not describe its search in enough detail to
#' reproduce it, so this function uses tabu search (Glover 1989) over moves
#' of single nodes between groups from `starts` random partitions: at each
#' step the best move that is not tabu is made (a tabu move is allowed when
#' it improves on the best partition found so far), a node that leaves a
#' group may not return for `tenure` steps, and a start ends after
#' `max_stall` steps without improvement. The best partition over all
#' starts is returned, the earliest start winning ties. Like any heuristic
#' it can miss the global optimum; the validation suite checks the fit
#' against exhaustive search on small networks. Starts run in parallel, and
#' the result for a given `seed` does not depend on `threads`.
#'
#' @inheritParams betweenness
#' @param k Number of groups.
#' @param fit `"errors"` (UCINET's default, Hamming criterion) or
#'   `"correlation"`.
#' @param starts Number of random starting partitions.
#' @param max_stall Steps without improvement after which a start ends.
#' @param tenure Tabu tenure (steps a node may not return to a group it
#'   left); `NULL` chooses `max(7, n / 4)`.
#' @param seed Seed for the starting partitions, a whole number between 0
#'   and \eqn{2^{53}}; `NULL` draws one from R's random number generator,
#'   so [set.seed()] makes results reproducible.
#'
#' @references
#' Borgatti, S. P., Everett, M. G., & Freeman, L. C. (2002). *UCINET for
#' Windows: Software for Social Network Analysis*. Analytic Technologies.
#'
#' Hanneman, R. A., & Riddle, M. (2005). *Introduction to Social Network
#' Methods*. University of California, Riverside.
#'
#' Doreian, P., Batagelj, V., & Ferligoj, A. (2005). *Generalized
#' Blockmodeling*. Cambridge University Press.
#' \doi{10.1017/CBO9780511584176}
#'
#' Borgatti, S. P., & Everett, M. G. (1999). Models of core/periphery
#' structures. *Social Networks*, 21(4), 375-395.
#' \doi{10.1016/S0378-8733(99)00019-2}
#'
#' Glover, F. (1989). Tabu search, part I. *ORSA Journal on Computing*,
#' 1(3), 190-206. \doi{10.1287/ijoc.1.3.190}
#'
#' @return An integer vector of group labels aligned with the rows of
#'   `x$nodes`, numbered 1, 2, ... in order of first appearance, with the
#'   number of errors in attribute `"errors"` and the correlation with the
#'   ideal pattern in attribute `"correlation"`.
#' @export
factions <- function(x, k = 2L, fit = c("errors", "correlation"), starts = 20L,
                     max_stall = 200L, tenure = NULL, seed = NULL, threads = 0L) {
  check_netframe(x)
  k <- check_count(k, "k")
  fit <- match.arg(fit)
  starts <- check_count(starts, "starts")
  max_stall <- check_count(max_stall, "max_stall")
  tenure <- if (is.null(tenure)) 0L else check_count(tenure, "tenure")
  seed <- groups_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call("factions",
                     .Call(C_factions, ptr, k, if (fit == "errors") 0L else 1L, starts,
                           max_stall, tenure, seed, threads))
  structure(res[[1L]], errors = res[[2L]], correlation = res[[3L]])
}

# A seed is a whole number in [0, 2^53]; NULL draws one from R's generator
# (two draws give 52 random bits, so set.seed() controls the result).
groups_seed <- function(seed) {
  if (is.null(seed)) {
    return(floor(stats::runif(1L) * 2^26) * 2^26 + floor(stats::runif(1L) * 2^26))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0 ||
      seed > 2^53 || seed != floor(seed)) {
    stop("`seed` must be a whole number between 0 and 2^53.", call. = FALSE)
  }
  as.double(seed)
}

#' Girvan-Newman, leading eigenvector, and greedy modularity communities
#'
#' Three classic partitioning methods for undirected networks, each
#' returning the partition with its modularity (see [modularity()]):
#' \describe{
#'   \item{`cluster_girvan_newman()`}{Repeatedly removes the edge with the
#'     highest edge betweenness, recomputing betweenness after every removal,
#'     so that the network falls apart into components (Girvan and Newman
#'     2002); of the partitions met along the way, the one with the highest
#'     modularity is returned (Newman and Girvan 2004). Edge betweenness is
#'     accumulated with Brandes' (2001) algorithm, recomputed only in the
#'     component that lost an edge.}
#'   \item{`cluster_leading_eigen()`}{Divides the network in two by the signs
#'     of the leading eigenvector of the modularity matrix, and each group
#'     again by its generalized modularity matrix, until no division raises
#'     modularity (Newman 2006).}
#'   \item{`cluster_fast_greedy()`}{Starts with every node in its own group
#'     and repeatedly merges the two connected groups whose merger raises
#'     modularity the most, returning the best partition along the way
#'     (Clauset, Newman, and Moore 2004; Newman 2004).}
#' }
#'
#' @details
#' All three accept only undirected networks, as the igraph functions they
#' correspond to (`cluster_edge_betweenness()`, `cluster_leading_eigen()`,
#' and `cluster_fast_greedy()`) do; a directed network is an error rather
#' than being silently symmetrized. Groups are numbered in order of their
#' first node.
#'
#' Girvan-Newman reads the network as a simple graph (parallel edges count
#' once, self-loops are ignored) and counts hops, without weights; the
#' reported modularity is that of the simple graph, as NetworkX's
#' `girvan_newman()` with `modularity()` on a `Graph` gives. The edge
#' removed at each step is the one with the largest betweenness; values
#' within a relative \eqn{10^{-9}} of the largest count as tied and the tie
#' goes to the edge whose endpoints come first in node order. igraph
#' removes the tied edge with the lowest edge id and NetworkX the first in
#' its iteration order, so on networks with tied betweenness the three can
#' split differently. The chosen partition is the first of highest
#' modularity in removal order. Cost is \eqn{O(m^2 n)} in the worst case;
#' searches from each source run in parallel.
#'
#' The leading eigenvector method and greedy merging treat weights (when
#' used) as in [modularity()]: they are added as given and must be
#' non-negative, parallel edges add up, and a self-loop of weight \eqn{w}
#' adds \eqn{2w} to its node's diagonal entry and degree. The leading
#' eigenvector of each group's matrix comes from a dense symmetric
#' eigendecomposition (Householder tridiagonalization and implicit QL;
#' Martin, Reinsch, and Wilkinson 1968; Bowdler et al. 1968), which costs
#' \eqn{O(n^3)} time and \eqn{O(n^2)} memory for the first division, so
#' networks above 3000 nodes are refused. A group is left undivided when the
#' leading eigenvalue is not positive (larger than \eqn{10^{-10}} times a
#' bound on the eigenvalues) or when the division does not raise
#' modularity. The eigenvector is scaled so that its largest component is
#' positive; nodes with positive components form one side and the others
#' the other. Results match igraph's `cluster_leading_eigen()` except where
#' an eigenvalue is degenerate (the network has a symmetry) or a component
#' is numerically zero, where the division is not unique.
#'
#' Greedy merging compares the change in modularity exactly when weights
#' are whole numbers and `resolution = 1`, and breaks ties toward the
#' groups with the lowest-numbered members. igraph and NetworkX break ties
#' by their own internal order, so on networks with tied mergers the
#' partitions can differ while having similar modularity; the validation
#' suite verifies that any such difference is a pure tie-break. Mergers
#' continue past the first that lowers modularity until no two groups are
#' connected, as in igraph, and the first partition of highest modularity
#' is returned. `resolution` is the parameter of Reichardt and Bornholdt
#' (2006), as in NetworkX's `greedy_modularity_communities()`.
#'
#' A network without edges puts every node in its own group with modularity
#' `NaN`.
#'
#' @inheritParams modularity
#' @param threads Worker threads for the betweenness computations; see
#'   [netframe()].
#'
#' @references
#' Girvan, M., & Newman, M. E. J. (2002). Community structure in social and
#' biological networks. *Proceedings of the National Academy of Sciences*,
#' 99(12), 7821-7826. \doi{10.1073/pnas.122653799}
#'
#' Newman, M. E. J., & Girvan, M. (2004). Finding and evaluating community
#' structure in networks. *Physical Review E*, 69(2), 026113.
#' \doi{10.1103/PhysRevE.69.026113}
#'
#' Brandes, U. (2001). A faster algorithm for betweenness centrality.
#' *Journal of Mathematical Sociology*, 25(2), 163-177.
#' \doi{10.1080/0022250X.2001.9990249}
#'
#' Newman, M. E. J. (2006). Modularity and community structure in networks.
#' *Proceedings of the National Academy of Sciences*, 103(23), 8577-8582.
#' \doi{10.1073/pnas.0601602103}
#'
#' Clauset, A., Newman, M. E. J., & Moore, C. (2004). Finding community
#' structure in very large networks. *Physical Review E*, 70(6), 066111.
#' \doi{10.1103/PhysRevE.70.066111}
#'
#' Newman, M. E. J. (2004). Fast algorithm for detecting community
#' structure in networks. *Physical Review E*, 69(6), 066133.
#' \doi{10.1103/PhysRevE.69.066133}
#'
#' Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community
#' detection. *Physical Review E*, 74(1), 016110.
#' \doi{10.1103/PhysRevE.74.016110}
#'
#' Martin, R. S., Reinsch, C., & Wilkinson, J. H. (1968). Householder's
#' tridiagonalization of a symmetric matrix. *Numerische Mathematik*,
#' 11(3), 181-195. \doi{10.1007/BF02161841}
#'
#' Bowdler, H., Martin, R. S., Reinsch, C., & Wilkinson, J. H. (1968). The
#' QR and QL algorithms for symmetric matrices. *Numerische Mathematik*,
#' 11(4), 293-306. \doi{10.1007/BF02166681}
#'
#' @return An integer vector of group labels aligned with the rows of
#'   `x$nodes`, numbered 1, 2, ... in order of first appearance, with the
#'   partition's modularity in attribute `"modularity"`.
#'   `cluster_girvan_newman()` also records the edges removed, in order, as
#'   a two-column data frame of node identifiers in attribute `"removed"`,
#'   the modularity before any removal and after each removal in attribute
#'   `"path"`, and the divisive hierarchy as an `hclust` object in attribute
#'   `"tree"` (heights count edge removals from the end, so `plot()` of it
#'   draws the dendrogram). `cluster_leading_eigen()` records the number of
#'   divisions in attribute `"splits"`.
#' @export
cluster_girvan_newman <- function(x, threads = 0L) {
  check_netframe(x)
  require_undirected(x, "Girvan-Newman clustering")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_girvan_newman", .Call(C_girvan_newman, ptr, threads))
  ids <- x$nodes[[x$roles$id]]
  removed <- data.frame(from = ids[res[[3L]]], to = ids[res[[4L]]], stringsAsFactors = FALSE)
  structure(res[[1L]], modularity = res[[2L]], removed = removed, path = res[[5L]],
            tree = removal_tree(length(ids), res[[3L]], res[[4L]], as.character(ids)))
}

# The divisive hierarchy of Girvan-Newman as an hclust object: read in
# reverse, the removals are merges. The edge removed last joins two
# singletons at height 1, the edge removed first makes the last merge at
# height m, and components that were never connected are joined above
# that, at height m + 1.
removal_tree <- function(n, from, to, labels) {
  m <- length(from)
  parent <- seq_len(n)
  find <- function(v) {
    while (parent[v] != v) v <- parent[v]
    v
  }
  # hclust's convention: singletons are negative, earlier merges positive.
  cluster <- -seq_len(n)
  merge <- matrix(0L, 0L, 2L)
  height <- numeric()
  join <- function(a, b, h) {
    ra <- find(a)
    rb <- find(b)
    if (ra == rb) return(FALSE)
    merge <<- rbind(merge, c(cluster[ra], cluster[rb]))
    height <<- c(height, h)
    parent[rb] <<- ra
    cluster[ra] <<- nrow(merge)
    TRUE
  }
  for (i in rev(seq_len(m))) join(from[i], to[i], m - i + 1)
  roots <- unique(vapply(seq_len(n), find, integer(1)))
  if (length(roots) > 1L) {
    for (r in roots[-1L]) join(roots[[1L]], r, m + 1)
  }
  if (n < 2L) {
    return(structure(list(merge = merge, height = height, order = seq_len(n), labels = labels,
                          method = "girvan-newman"), class = "hclust"))
  }
  tree <- structure(list(merge = merge, height = height, order = integer(), labels = labels,
                         method = "girvan-newman"), class = "hclust")
  # Leaf order from the merge structure, as stats computes it.
  tree$order <- leaf_order(merge)
  tree
}

leaf_order <- function(merge) {
  # Depth-first over the merge tree with an explicit stack (no recursion,
  # so large trees cannot exhaust R's expression nesting).
  out <- integer()
  stack <- nrow(merge)
  while (length(stack)) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (node < 0) {
      out <- c(out, -node)
    } else {
      stack <- c(stack, merge[node, 2L], merge[node, 1L])
    }
  }
  out
}

#' @rdname cluster_girvan_newman
#' @export
cluster_leading_eigen <- function(x, weighted = NULL) {
  check_netframe(x)
  require_undirected(x, "The leading eigenvector method")
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_leading_eigen", .Call(C_leading_eigen, ptr, weighted))
  structure(res[[1L]], modularity = res[[2L]], splits = res[[3L]])
}

#' @rdname cluster_girvan_newman
#' @export
cluster_fast_greedy <- function(x, resolution = 1, weighted = NULL) {
  check_netframe(x)
  require_undirected(x, "Greedy modularity optimization")
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_fast_greedy", .Call(C_fast_greedy, ptr, resolution, weighted))
  structure(res[[1L]], modularity = res[[2L]])
}

require_undirected <- function(x, what) {
  if (x$directed) {
    stop(what, " needs an undirected network; rebuild it with ",
         "netframe(..., directed = FALSE) or symmetrize() it.", call. = FALSE)
  }
}

#' Hierarchical clustering of nodes
#'
#' Agglomerative (Johnson 1967) hierarchical clustering of the nodes from
#' any distance or similarity matrix, such as [distances()] or
#' [structural_equivalence()], as in UCINET's Tools > Cluster Analysis >
#' Hierarchical routine. `cluster_hierarchical()` wraps [stats::hclust()];
#' `cut_clusters()` cuts the tree into groups numbered in order of first
#' appearance.
#'
#' @details
#' A similarity matrix (`similarity = TRUE`, for example correlations
#' between rows of the adjacency matrix) is turned into distances by
#' subtracting it from its largest finite value. Unreachable pairs
#' (`Inf` in a distance matrix) are given one more than the largest finite
#' distance, so they join last; missing values are an error. The matrix
#' must be square and symmetric; its dimnames label the leaves.
#'
#' The `"average"`, `"single"`, and `"complete"` methods are UPGMA, nearest
#' neighbor, and farthest neighbor (Sokal and Michener 1958; Johnson 1967).
#' `"ward"` is Ward's (1963) method on squared distances, [stats::hclust()]'s
#' `"ward.D2"`.
#'
#' @param d A square numeric matrix of distances (or similarities) between
#'   nodes.
#' @param method Linkage: `"average"`, `"single"`, `"complete"`, or
#'   `"ward"`.
#' @param similarity Whether `d` holds similarities (larger means closer)
#'   rather than distances.
#' @param tree The result of `cluster_hierarchical()` (or any `hclust`).
#' @param k Number of groups to cut into, or `NULL` to cut at `height`.
#' @param height Height at which to cut, when `k` is `NULL`.
#'
#' @references
#' Johnson, S. C. (1967). Hierarchical clustering schemes. *Psychometrika*,
#' 32(3), 241-254. \doi{10.1007/BF02289588}
#'
#' Sokal, R. R., & Michener, C. D. (1958). A statistical method for
#' evaluating systematic relationships. *University of Kansas Science
#' Bulletin*, 38, 1409-1438.
#'
#' Ward, J. H. (1963). Hierarchical grouping to optimize an objective
#' function. *Journal of the American Statistical Association*, 58(301),
#' 236-244. \doi{10.1080/01621459.1963.10500845}
#'
#' @return `cluster_hierarchical()`: an object of class `hclust`, with the
#'   row names of `d` as labels. `cut_clusters()`: an integer vector of group
#'   labels in the order of the rows of `d`, numbered 1, 2, ... in order of
#'   first appearance.
#' @export
cluster_hierarchical <- function(d, method = c("average", "single", "complete", "ward"),
                                 similarity = FALSE) {
  method <- match.arg(method)
  similarity <- check_flag(similarity, "similarity")
  if (!is.matrix(d) || !is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a square numeric matrix.", call. = FALSE)
  }
  if (anyNA(d)) stop("`d` must not contain missing values.", call. = FALSE)
  if (!isTRUE(all.equal(d, t(d), check.attributes = FALSE))) {
    stop("`d` must be symmetric.", call. = FALSE)
  }
  finite <- d[is.finite(d)]
  if (similarity) {
    if (any(!is.finite(d))) stop("A similarity matrix must be finite.", call. = FALSE)
    d <- max(finite) - d
  } else if (any(is.infinite(d))) {
    d[is.infinite(d)] <- if (length(finite)) max(finite) + 1 else 1
  }
  labels <- rownames(d)
  if (is.null(labels)) labels <- as.character(seq_len(nrow(d)))
  dist <- stats::as.dist(d)
  attr(dist, "Labels") <- labels
  stats::hclust(dist, method = if (method == "ward") "ward.D2" else method)
}

#' @rdname cluster_hierarchical
#' @export
cut_clusters <- function(tree, k = NULL, height = NULL) {
  if (!inherits(tree, "hclust")) {
    stop("`tree` must be the result of cluster_hierarchical() (an `hclust` object).",
         call. = FALSE)
  }
  if (is.null(k) == is.null(height)) {
    stop("Give exactly one of `k` and `height`.", call. = FALSE)
  }
  if (!is.null(k)) {
    k <- check_count(k, "k")
    if (k < 1L || k > length(tree$order)) {
      stop("`k` must be between 1 and the number of nodes.", call. = FALSE)
    }
    groups <- stats::cutree(tree, k = k)
  } else {
    if (!is.numeric(height) || length(height) != 1L || is.na(height)) {
      stop("`height` must be a single number.", call. = FALSE)
    }
    groups <- stats::cutree(tree, h = height)
  }
  out <- match(groups, unique(groups))
  names(out) <- NULL
  out
}
