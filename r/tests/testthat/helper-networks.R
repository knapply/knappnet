# Small networks shared by several test files. testthat sources helper-*.R
# before the tests, so these are available everywhere; a test file defines
# any network only it uses.

# Krackhardt's kite (Krackhardt 1990) as an undirected edge list on nodes
# 1 to 10, and as a netframe.
kite_edges <- function() {
  data.frame(
    source = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    target   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
}

kite <- function() net_frame(kite_edges(), nodes = data.frame(id = 1:10), directed = FALSE)

# A path 1 - 2 - 3, an edge 4 - 5, and an isolate 6.
pieces <- function() {
  net_frame(data.frame(source = c(1, 2, 4), target = c(2, 3, 5)),
           nodes = data.frame(id = 1:6), directed = FALSE)
}

# A random simple directed network on n nodes with about m edges.
random_digraph <- function(n, m, seed) {
  set.seed(seed)
  e <- data.frame(source = sample(n, m, TRUE), target = sample(n, m, TRUE))
  e <- e[e$source != e$target, ]
  e[!duplicated(e), ]
}
