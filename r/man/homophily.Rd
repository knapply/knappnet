\name{homophily}
\alias{homophily}
\title{Within-group tie shares}
\description{
For each group, the share of its members' ties that stay inside the group,
against the share expected if ties ignored groups: the baseline homophily of
McPherson, Smith-Lovin, and Cook (2001).
}
\usage{
homophily(x, groups, direction = c("all", "out", "in"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{groups}{Group of each node: the name of a column of \code{x$nodes}, or
a vector with one value per node. Missing values are not allowed.}
\item{direction}{For directed networks: ties sent (\code{"out"}), received
(\code{"in"}), or both (\code{"all"}) by members.}
}
\details{
Ties are counted from the members' side: a tie between two members counts
once for each, so \code{internal} is twice the number of internal ties in an
undirected network (and with \code{direction = "all"} in a directed one).
\code{observed} is internal / (internal + external), \code{NaN} for a group
without ties. \code{expected} is \eqn{(n_g - 1)/(n - 1)}, the share of a
member's possible contacts who belong to its group, which is the expected
share when every tie is equally likely; \code{observed} above \code{expected}
indicates ties concentrated within the group. Self-loops are ignored,
parallel edges each count, and weights are not used.
}
\value{A data frame with one row per group (in order of first appearance):
\code{group}, \code{size}, \code{internal}, \code{external},
\code{observed}, and \code{expected}.}
\references{
McPherson, M., Smith-Lovin, L., & Cook, J. M. (2001). Birds of a feather:
Homophily in social networks. \emph{Annual Review of Sociology}, 27, 415-444.
\doi{10.1146/annurev.soc.27.1.415}
}
\seealso{
\code{\link{ei_index}}, \code{\link{assortativity_nominal}}
}
\examples{
staff <- data.frame(name = c("ann", "bob", "cat", "dan"),
                    dept = c("sales", "sales", "eng", "eng"))
g <- netframe(data.frame(from = c("ann", "ann", "bob", "cat"),
                         to = c("bob", "cat", "dan", "dan")),
              nodes = staff, directed = FALSE)
homophily(g, "dept")
}
