\name{factions}
\alias{factions}
\title{Factions}
\usage{
factions(x, k = 2L, fit = c("errors", "correlation"), starts = 20L,
         max_stall = 200L, tenure = NULL, seed = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{k}{Number of groups.}
\item{fit}{\code{"errors"} (UCINET's default, the Hamming criterion) or
\code{"correlation"}; see Description.}
\item{starts}{Number of random starting partitions.}
\item{max_stall}{Steps without improvement after which a start ends.}
\item{tenure}{Tabu tenure (steps a node may not return to a group it
left); \code{NULL} chooses \code{max(7, n / 4)}.}
\item{seed}{Seed for the starting partitions, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{
An integer vector of group labels aligned with the rows of \code{x$nodes},
numbered 1, 2, ... in order of first appearance, with the number of errors
in attribute \code{"errors"} and the correlation with the ideal pattern in
attribute \code{"correlation"}.
}
\description{
Partitions the nodes into \code{k} groups that best match the ideal of
cohesive, mutually separate groups, as UCINET's Network > Subgroups >
Factions routine does (Borgatti, Everett, and Freeman 2002; Hanneman and
Riddle 2005, chapter 11). The ideal pattern has every pair within a group
tied and no pair across groups tied. With \code{fit = "errors"} the
partition minimizes the number of errors, absent ties within groups plus
present ties between groups (the Hamming distance to the ideal, the
criterion of a binary blockmodel with complete diagonal blocks; Doreian,
Batagelj, and Ferligoj 2005); with \code{fit = "correlation"} it maximizes
the Pearson correlation between the adjacency matrix and the ideal pattern
(the criterion Borgatti and Everett 1999 use for core-periphery models).
}
\details{
The adjacency matrix is binary: parallel edges count once, self-loops are
ignored, and weights are not used. Undirected networks compare each
unordered pair once; directed networks compare both ordered pairs. Every
group is non-empty (\code{1 <= k <= n}). The correlation is undefined when
the network has no ties or every possible tie, when \code{k = 1}, or when
\code{k = n}; \code{fit = "correlation"} then stops with an error, and the
\code{correlation} attribute of an errors fit is \code{NaN}.

Finding the best partition is a combinatorial optimization problem.
UCINET's documentation does not describe its search in enough detail to
reproduce it, so this function uses tabu search (Glover 1989) over moves of
single nodes between groups from \code{starts} random partitions: at each
step the best move that is not tabu is made (a tabu move is allowed when it
improves on the best partition found so far), a node that leaves a group may
not return for \code{tenure} steps, and a start ends after
\code{max_stall} steps without improvement. The best partition over all
starts is returned, the earliest start winning ties. Like any heuristic it
can miss the global optimum; the validation suite checks the fit against
exhaustive search on small networks. Starts run in parallel, and the result
for a given \code{seed} does not depend on \code{threads}.

The two criteria often agree but not always: the error count treats a
missing tie within a group and a present tie between groups alike, while
the correlation weighs them by how rare ties are and by the sizes of the
groups, so on a sparse network it tolerates missing ties within a large
group more readily.
}
\references{
Borgatti, S. P., Everett, M. G., & Freeman, L. C. (2002). \emph{UCINET for
Windows: Software for Social Network Analysis}. Analytic Technologies.

Hanneman, R. A., & Riddle, M. (2005). \emph{Introduction to Social Network
Methods}. University of California, Riverside.

Doreian, P., Batagelj, V., & Ferligoj, A. (2005). \emph{Generalized
Blockmodeling}. Cambridge University Press. \doi{10.1017/CBO9780511584176}

Borgatti, S. P., & Everett, M. G. (1999). Models of core/periphery
structures. \emph{Social Networks}, 21(4), 375-395.
\doi{10.1016/S0378-8733(99)00019-2}

Glover, F. (1989). Tabu search, part I. \emph{ORSA Journal on Computing},
1(3), 190-206. \doi{10.1287/ijoc.1.3.190}
}
\seealso{\code{\link{core_periphery}}, \code{\link{cluster_girvan_newman}},
\code{\link{modularity}}}
\examples{
# Two triangles joined by one tie, plus a pendant on the second.
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5, 6), to = c(2, 3, 3, 4, 5, 6, 6, 7)),
              directed = FALSE)
f <- factions(g, k = 2, seed = 1, threads = 2)
f
attr(f, "errors")

# The two criteria can prefer different partitions. Here a square
# 1-5-4-6-1 has a pendant 3 on node 4, and node 2 is an isolate: counting
# errors pairs the isolate with the pendant (UCINET's default), while the
# correlation criterion leaves the isolate alone.
h <- netframe(data.frame(from = c(1, 1, 3, 4, 4), to = c(5, 6, 4, 5, 6)),
              nodes = data.frame(id = 1:6), directed = FALSE)
factions(h, k = 2, fit = "errors", seed = 1, threads = 2)       # UCINET Factions (Hamming)
factions(h, k = 2, fit = "correlation", seed = 1, threads = 2)  # correlation with the ideal
}
