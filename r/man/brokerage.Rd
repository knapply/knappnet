\name{brokerage}
\alias{brokerage}
\title{Gould-Fernandez brokerage roles}
\description{
Counts, for each node, the ordered pairs of other nodes it brokers in each
of the five roles of Gould and Fernandez (1989).
}
\usage{
brokerage(x, groups, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{groups}{The group of each node: a vector aligned with the rows of
\code{x$nodes}, or the name of a column of \code{x$nodes}. Only equality of
values matters. Missing values are not allowed.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
Node j brokers the ordered pair (i, k) when \eqn{i \to j \to k}{i -> j -> k}
and there is no tie \eqn{i \to k}{i -> k}, with i, j, and k distinct. With
the groups of i, j, and k written in order, the roles are:
\describe{
  \item{coordinator}{(A, A, A): all three in the same group.}
  \item{itinerant}{(A, B, A): the broker mediates between two members of
    another group; also called consultant.}
  \item{representative}{(A, A, B): the broker passes ties from its own group
    to another.}
  \item{gatekeeper}{(A, B, B): the broker passes ties from another group
    into its own.}
  \item{liaison}{(A, B, C): three different groups.}
}

The roles are defined on binary ties: weights are ignored, parallel edges
count once, and self-loops are ignored. An undirected edge counts as a tie
in both directions, so in an undirected network each unordered pair a node
brokers is counted twice, once in each direction. The counts match
\code{raw.nli} from \code{sna::brokerage()} (columns \code{w_I}, \code{w_O},
\code{b_IO}, \code{b_OI}, \code{b_O}, \code{t}).
}
\value{A data frame with one row per node: \code{node} (the identifier),
then the counts \code{coordinator}, \code{itinerant}, \code{representative},
\code{gatekeeper}, \code{liaison}, and \code{total}.}
\references{
Gould, R. V., & Fernandez, R. M. (1989). Structures of mediation: A formal
approach to brokerage in transaction networks. \emph{Sociological
Methodology}, 19, 89-126. \doi{10.2307/270949}
}
\examples{
g <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "d")),
              nodes = data.frame(id = c("a", "b", "c", "d"),
                                 team = c("x", "x", "y", "y")))
brokerage(g, "team")
}
