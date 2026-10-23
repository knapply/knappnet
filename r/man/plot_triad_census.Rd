\name{plot_triad_census}
\alias{plot_triad_census}
\title{Plot the triad census}
\description{
Draws the \code{\link{triad_census}} as a bar chart and, by default, the
counts expected under the U|MAN null model of Holland and Leinhardt
(1976), which conditions on the dyad census: all networks with the same
numbers of mutual, asymmetric, and null dyads are equally likely.
Comparing the two shows which triadic structures (transitive triads,
cycles, closed triangles) are over- or under-represented beyond what the
density and reciprocity of the network already imply.
}
\usage{
plot_triad_census(x, expected = TRUE, log = FALSE, threads = 0L, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{expected}{Whether to draw the expected counts under U|MAN beside the
observed ones.}
\item{log}{Whether to use a logarithmic count axis (counts of 0 are drawn
at the bottom).}
\item{threads}{Number of threads for the census, a whole number; \code{0} uses
all available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
\item{...}{Further arguments passed to \code{\link[graphics]{barplot}},
such as \code{main}.}
}
\details{
With \eqn{M}, \eqn{A}, and \eqn{N} mutual, asymmetric, and null dyads
among the \eqn{D = n(n-1)/2} dyads, the three dyads of a random triad
under U|MAN have \eqn{m} mutual, \eqn{a} asymmetric, and \eqn{k} null
dyads with probability
\deqn{{M \choose m} {A \choose a} {N \choose k} / {D \choose 3},}{choose(M, m) choose(A, a) choose(N, k) / choose(D, 3),}
and each orientation of the asymmetric dyads is equally likely, so the
expected count of a triad type is \eqn{{n \choose 3}}{choose(n, 3)} times
that probability times the share of orientations giving the type (1/4,
1/4, 1/2 for 021D, 021U, 021C; 1/2, 1/2 for 111D, 111U; 3/4, 1/4 for 030T,
030C; 1/4, 1/4, 1/2 for 120D, 120U, 120C; 1 otherwise). These are the
expected values of Holland and Leinhardt (1976; Wasserman and Faust 1994,
ch. 14) and were checked against the average census over every digraph on
five nodes. The dyad census reads the network as \code{\link{triad_census}}
does: self-loops ignored, parallel edges counted once, and every edge of
an undirected network mutual, so only the four undirected types (003, 102,
201, 300) are drawn for undirected networks.
}
\value{A data frame with one row per triad type drawn, invisibly:
\code{type}, \code{observed}, and \code{expected} (\code{NA} when
\code{expected = FALSE}).}
\references{
Holland, P. W., & Leinhardt, S. (1976). Local structure in social
networks. \emph{Sociological Methodology}, 7, 1-45. \doi{10.2307/270703}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\seealso{\code{\link{triad_census}}, \code{\link{transitivity}},
\code{\link{reciprocity}}.}
\examples{
g <- netframe(data.frame(from = c(1, 2, 3, 3, 4, 5, 5, 6, 1),
                         to   = c(2, 3, 1, 4, 5, 3, 6, 4, 6)))
triad_census(g)
plot_triad_census(g, threads = 2)
plot_triad_census(example_network("karate"), log = TRUE, threads = 2)
}
